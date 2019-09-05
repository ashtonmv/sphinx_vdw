// ---------------------------------------------------------------------------
//
//      The ab-initio based multiscale library
//
//                  S / P H I / n X
//
//      Copyright:  Max-Planck-Institute for Iron Research
//                  40237 Duesseldorf, Germany
//
//      Contact:    https://sxlib.mpie.de
//      Authors:    see sphinx/AUTHORS
//      License:    see sphinx/LICENSE
//
// ---------------------------------------------------------------------------

#include <SxSpinConstraint.h>
#include <SxFermi.h>
#include <SxPAWSet.h>
#include <SxPAWRho.h>
#include <SxRho.h>
#include<SxHamiltonian.h>
#include <SxPAWHamiltonian.h>
#include <SxSpectrum.h>
#include <SxPartialWaveBasis.h>
#include <SxAtomicStructure.h>
#include <SxProjector.h>
#include <SxRBasis.h>
#include <SxGkBasis.h>
#include <SxGBasis.h>
#include <SxMesh3D.h>
#include <SxSymMatrix.h>
#include <SxSimpleParser.h>
#include <SxFileParser.h>

namespace Timer {
   enum SpinConstraintTimers { OptimizeNu, NuIntro , PAWInt, OmegaInit,
                               OmegaInit1, OmegaInit2,
                               SpinsDiag, NuDiagH, OmegaRot, NuRest};
}

SX_REGISTER_TIMERS (Timer::SpinConstraintTimers)
{
   using namespace Timer;
   regTimer (OptimizeNu, "optimize nu");
   regTimer (NuIntro   , "opt. nu startup");
   regTimer (PAWInt    , "spin constr. PAW int");
   regTimer (OmegaInit , "setup omega_nn");
   regTimer (OmegaInit1 , "setup omega_nn 1");
   regTimer (OmegaInit2 , "setup omega_nn 2");
   regTimer (SpinsDiag , "spinsDiag total");
   regTimer (NuDiagH   , "spinsDiag H diag");
   regTimer (OmegaRot  , "spinsDiag rot. omega");
   regTimer (NuRest    , "spinsDiag rest");
}

SxSpinConstraint::SxSpinConstraint ( const SxPtr<SxPAWPot> &pawPotPtrIn)
   : wvPtr (NULL)
{
   pawPotPtr = pawPotPtrIn;

}

void SxSpinConstraint::read (const SxSymbolTable *topLvl,
                             const SxAtomicStructure &structure)
{
   // --- create variables for spin constraint
   targetSpin.resize (structure.getNAtoms());
   targetSpin.set (0.);
   constrainedAtom.resize (structure.getNAtoms ());
   constrainedAtom.set (false);
   // --- read the target spins from input.sx
   SYMBOLPARSE(topLvl)  {
      FOREACH_SYMBOLGROUP("spinConstraint")  {
         SxString constraintFile = SYMBOLGET("file") || "";
         if (constraintFile.getSize () > 0)  {
            SxFileParser fp(constraintFile);
            fp.topic ("spin constraints"); // get nicer errors
            SX_LOOP(ia)  {
               fp.skipWhite ();
               if (fp.reads("X") || fp.reads("x"))  {
                  constrainedAtom(ia) = false;
                  fp.nextLine ();
               } else {
                  constrainedAtom(ia) = true;
                  fp >> targetSpin(ia);
               }
            }
         } else {
            const SxArray<SxString> &labels = structure.getLabels ();

            SxString labelConstraint = SYMBOLGET("label");
            double   spinConst       = SYMBOLGET("constraint");

            for (int iAtom = 0; iAtom < structure.getNAtoms(); iAtom++) {
               if ( labelConstraint == labels(iAtom) ) {
                  targetSpin(iAtom) = spinConst;
                  constrainedAtom(iAtom) = true;
               }
            }
         }
      }
   }
   cout << SX_SEPARATOR;
   cout << "| Spin constraints" << endl;
   cout << SX_SEPARATOR;
   SX_LOOP(ia)  {
      sxprintf("| Atom %3d ",int(ia+1));
      if (constrainedAtom(ia))
         cout << "spin=" << targetSpin(ia);
      else
        cout << "not constrained.";
      cout << endl;
   }
   cout << SX_SEPARATOR;
   checkSym (structure);
}

SxArray<SxMatrix<Complex16> >
SxSpinConstraint::getOmegaNN (int iSpin, int ik) const
{
   SX_CLOCK (Timer::OmegaInit);
   SX_CHECK (pawPotPtr);
   SX_CHECK (wvPtr);
   const SxDiracVec<TPrecCoeffG> psi = (*wvPtr)(iSpin, ik);
   int nAtom = (int)targetSpin.getSize ();
   int nStates = (int)psi.nCols ();
   const SxPAWPot &pawPot = *pawPotPtr;
   const SxArray<SxMatrix<Double> > &omegaIJ = pawPot.omegaPAW;
   const SxAtomicStructure &structure = wvPtr->getGkBasisPtr ()->getTau ();

   SxArray<SxMatrix<Complex16> > result(nAtom);

   ssize_t ng = psi.getBasis<SxPAWBasis> ().gBasis->ng;
   ssize_t ipBasis = ng;
   SxVector<TPrecCoeffG> projAtom; // keep memory allocated
   SxVector<TPrecCoeffG> projOmega;

   // --- loop over atoms
   for (int iAtom = 0; iAtom < nAtom; iAtom++) {
      int is = structure.getISpecies(iAtom);
      int npt = pawPot.getNProjType (is);
      int npl = pawPot.getNProj(is);
      if (!constrainedAtom(iAtom)) {
         ipBasis += npl;
         continue;
      }
      // --- get compact copy of projectors for this atom
      projAtom. reformat(npl, nStates);
      for (ssize_t iState = 0; iState < nStates; ++iState)
         for (ssize_t ipl = 0; ipl < npl; ++ipl)
            projAtom(ipl, iState) = psi(ipBasis + ipl, iState);

      //SX_START_TIMER(Timer::OmegaInit1);
      // --- apply PAW volume operator on projAtom (and transpose)
      projOmega.reformat (nStates, npl);
      projOmega.set (0.);
      for (int ipt = 0, ipl0 = 0; ipt < npt; ipt++) {
         int li = pawPot.lPhi(is)(ipt);
         int nm = 2 * li + 1;
         for (int jpt = 0, jpl0 = 0; jpt < npt; jpt++) {
            int lj = pawPot.lPhi(is)(jpt);
            if (li == lj) {
               double integral = omegaIJ(is)(ipt, jpt);
               if (iSpin == 1) integral = -integral;
               for (int im = 0; im < nm; ++im)  {
                  for (int iState = 0; iState < nStates; ++iState)
                     projOmega(iState, ipl0 + im) += integral * projAtom(jpl0 + im, iState).conj ();
               }
            }
            jpl0 += 2 * lj + 1;
         }
         ipl0 += nm;
      }
      //SX_STOP_TIMER(Timer::OmegaInit1);
      //SX_CLOCK (Timer::OmegaInit2);
      // get matrix elements
      result(iAtom) = projOmega ^ projAtom;
      // next atom
      ipBasis += npl;
   }
   return result;
} 

void SxSpinConstraint::setNuA(double value)
{
   SX_CHECK (targetSpin.getSize () > 0);
   nuA.resize (targetSpin.getSize ());
   nuA.set (value);
}
void SxSpinConstraint::setConstraints(const SxVector<Double> &targetSpinIn,
                                      const SxAtomicStructure &str)
{
   SX_CHECK(targetSpinIn.getSize () == str.getNAtoms (),
            targetSpinIn.getSize (), str.getNAtoms ());
   targetSpin.resize (str.getNAtoms ());
   SX_LOOP(ia)
      targetSpin(ia) = constrainedAtom(ia) ? targetSpinIn(ia) : 0.;
   checkSym (str);
}

void SxSpinConstraint::checkSym (const SxAtomicStructure &str) const
{
   if (!str.cell.symGroupPtr) return;
   const SxSymGroup &syms =  *str.cell.symGroupPtr;
   if (syms.getNSymmorphic () == 1) return;

   SxGrid grid(str, 10);
   bool fatal = false;
   SxArray2<bool> fatalIJ;
   for (int iSym = 0; iSym < syms.getNSymmorphic (); ++iSym)  {
      SxAtomInfo::ConstPtr map = str.match (grid,
                                            syms.getSymmorphic (iSym) ^ str);
      if (!map)  {
         cout << "Symmetry inconsistency!" << endl;
         cout << "Please check the symmetries of your structure." << endl;
         cout << "If you are not close to the numerical noise limit ("
              << str.epsEqual << ")," << endl
              << "contact the developers with your structure file" << endl;
         SX_EXIT;
      }
      SX_LOOP(ia)  {
         int ja = (int)map->parentMap(ia);
         if (fabs(targetSpin(ia) - targetSpin(ja)) > 1e-3
             || (constrainedAtom(ia)? 0 : 1) != (constrainedAtom(ja)? 0 : 1))  {
            if (!fatal)  {
               cout << endl << SX_SEPARATOR;
               cout << "ERROR: Symmetry inconsistent spin constraints:" << endl;
               cout << SX_SEPARATOR;
               int nAtoms = str.getNAtoms ();
               fatalIJ.reformat (nAtoms, nAtoms);
               fatalIJ.set (false);
            }
            if (!fatalIJ(ia,ja))  {
               cout << "Atom " << (ia + 1) << " @ " << str(ia) << " and" << endl
                    << "atom " << (ja + 1) << " @ " << str(ja) << " are symmetry equivalent"
                    << endl;
               cout << "Target spin " << (ia + 1) << ": " << targetSpin(ia);
               if (!constrainedAtom(ia)) cout << " (not constrained)";
               cout << endl;
               cout << "Target spin " << (ja + 1) << ": " << targetSpin(ja);
               if (!constrainedAtom(ja)) cout << " (not constrained)";
               cout << endl;
               fatalIJ(ia, ja) = fatalIJ(ja, ia) = true;
            }
            fatal = true;
         }
      }
   }
   if (fatal)  {
      cout << "=================" << endl;
      cout << "--- WHAT TO DO?" << endl;
      cout << "=================" << endl;
      cout << "* Fix the target spins to a symmetry-compatible set"
           << endl;
      cout << "* use labels to reduce the symmetry" << endl;
      cout << "* specify symmetry in the input file" << endl;
      cout << "  To switch off all symmetries, insert" << endl
           << "     symmetry { operator { S = [[1,0,0],[0,1,0],[0,0,1]]; } }"
           << endl << "  into the structure {} group" << endl;
      SX_QUIT;
   }
   cout << "Spin constraints: symmetry check OK" << endl;
}

double SxSpinConstraint::calcKappaOpt (const SxVector<Double> &MSpinIn, const SxVector<Double> &MSpinPlusIn, const double kappa)
{
   if ((MSpinIn - targetSpin).normSqr () / double(MSpinIn.getSize ()) > 1e-2)  {
      cout << "MSpinIN   and MSpinPlusIN   = " << MSpinIn << "  " << MSpinPlusIn << endl;
      cout << "targetSpin  =  " << targetSpin << endl;
   }
   double sumK = 0, sumK2 = 0; 
   sumK = dot( (MSpinIn - targetSpin) , (MSpinPlusIn - MSpinIn) );
   sumK2 = (MSpinPlusIn - MSpinIn).normSqr();
   return - (sumK * kappa)/(sumK2);
}

SxVector<Double>
SxSpinConstraint::computeNu (const SxPtr<SxPAWSet> &wavesPtr, SxFermi &fermi,
                             const Real8 ekt, double epsilon)
{
   SX_CLOCK(Timer::OptimizeNu);
   SX_START_TIMER (Timer::NuIntro);
   SxPAWSet &waves = *wavesPtr;
   double kappa = 0.01, kappaOpt = 0, meanError;
   int nSpin = fermi.getNSpin ();
   int nk = fermi.getNk ();
   int nStates = fermi.getNStates ();
   Eps eps0(nStates,nSpin,nk);
   ssize_t nAtom = nuA.getSize();

   wvPtr = &waves;
   for (int iSpin = 0; iSpin < nSpin; iSpin++) {
         for (int ik = 0; ik < nk; ik++) {
            SX_MPI_LEVEL("waves-k");
            if (SxLoopMPI::myWork(ik)) {
               eps0(iSpin,ik) = fermi.eps(iSpin,ik).getCopy();
            }//LoopMPI
         }//ik
   }//iSpin
   SxVector<Double> MSpin(nAtom), MSpinPlus(nAtom);
   // --- iteration loop for spin constraint
   int nStep = 100;
   SxVector<Double> searchOld, deltaSpinOld;
   double meanOld = 0.;
   SX_STOP_TIMER(Timer::NuIntro);
   for (int iStep = 0; iStep < nStep; iStep++) {
      //cout << "nuA=" << nuA << endl;
      MSpin = getSpinsDiag(eps0, fermi, ekt);
      //cout << "MSpin=" << MSpin << endl;
      SxVector<Double> deltaSpin = MSpin - targetSpin;
      SX_LOOP(ia) if (!constrainedAtom(ia)) deltaSpin(ia) = 0.;
      SxVector<Double> searchNu = deltaSpin.getCopy ();

      meanError = deltaSpin.normSqr() * 1./double(nAtom);
      if ( meanError < epsilon ) break;

      if (iStep > 0) {
         // Fletcher-Reeves
         double gamma = meanError / meanOld;
         // Polak-Ribiere
         //double gamma = (meanError - dot(deltaSpinOld,deltaSpin)/double(nAtom)) / meanOld;
         // Hestenes-Stiefel
         //double gamma = dot(deltaSpinOld-deltaSpin,deltaSpin)/dot(deltaSpinOld - deltaSpin,searchOld);
         cout << "gamma=" << gamma << " R_M^2=" << meanError << endl;
         searchNu += gamma * searchOld;
      } else {
         cout << "gamma=X R_M^2=" << meanError << endl;
      }

      nuA += kappa * searchNu;
      MSpinPlus = getSpinsDiag(eps0, fermi, ekt);

      // --- calculate optimal kappa with test Spin^(i+1)
      kappaOpt = calcKappaOpt(MSpin, MSpinPlus, kappa);
      cout << "kappa opt=" << kappaOpt << endl;
      // --- update nuA for next iteration
      if (kappaOpt * searchNu.maxval () > 3. * ekt) {
         // note: factor 3 * ekt is not tested
         kappaOpt = 3. * ekt / searchNu.maxval ();
         cout << "Restricting step size" << endl;
      }
      nuA += (kappaOpt - kappa) * searchNu;
      //cout << "delta Nu=" << nuA << endl;
      // check for abort criterion
      searchOld = searchNu;
      deltaSpinOld = deltaSpin.getCopy ();
      meanOld = meanError;

      // --- adjust trial step
      double g = 1.5 * fabs(kappaOpt) / kappa;
      if (g > 2.) g = 2.;
      if (g < 0.5) g = 0.5;
      kappa *= pow(g,0.7);
      cout << "kappa=" << kappa << endl;
   } // iStep

   // rotate waves
   getSpinsDiag (eps0, fermi, ekt, wavesPtr);
   wvPtr = NULL;
   return nuA;
}

SxVector<Double> SxSpinConstraint::getSpinsDiag ( const Eps &eps0In, SxFermi &fermiIn, const double &ektIn , SxPtr<SxPAWSet> wavesPtr ) 
{
   SX_CHECK (wvPtr);
   SX_CLOCK (Timer::SpinsDiag);
   int nSpin = fermiIn.getNSpin ();
   int nk = fermiIn.getNk ();
   int nStates = fermiIn.getNStates ();
   ssize_t nAtom = nuA.getSize();
   SxSymMatrix<Complex16>::Eigensystem eig;
   SxArray<Eps> omega(nAtom);
   SX_MPI_LEVEL("waves-k");
   // --- selective resize of omega
   SX_LOOP(iAtom)  {
      // only constrained atoms
      if (constrainedAtom(iAtom))  {
         omega(iAtom).bundle.resize (nk);
         for (int ik = 0; ik < nk; ik++) {
            // only ik's for this MPI task
            if (SxLoopMPI::myWork(ik))  {
               omega(iAtom).bundle(ik).resize (nSpin);
               SX_LOOP(iSpin) omega(iAtom).bundle(ik)(iSpin).resize (nStates);
            }
         }
      }
   }
   for( int iSpin = 0; iSpin < nSpin; iSpin++) {
      for (int ik = 0; ik < nk; ik++) {
         if (!SxLoopMPI::myWork(ik)) continue;
         // --- set up new Hamiltonian matrix for current deltaNu
         SxSymMatrix<Complex16> Hnew(nStates);
         Hnew.set (0.);
         // diagonal elements
         SX_LOOP(n) Hnew(n,n) += eps0In(n,iSpin,ik);
         // spin constraint contribution
         SxArray<SxMatrix<Complex16> > omegaNN = getOmegaNN(iSpin,ik);
         SX_LOOP(iAtom)  {
            if (!constrainedAtom(iAtom)) continue;
            for (ssize_t i = 0; i < nStates; ++i)
               for (ssize_t j = 0; j <=i ; ++j)
                  Hnew(j,i) +=  nuA(iAtom) * omegaNN(iAtom)(j,i);
         }

         // --- diagonalize
         { 
            SX_CLOCK (Timer::NuDiagH);
            eig = Hnew.eigensystem();
            SX_LOOP(in) fermiIn.eps(in,iSpin,ik) = eig.vals(in);
         }

         // --- rotate omegann according to new eigenfunctions
         SX_START_TIMER(Timer::OmegaRot);
         for( int iAtom = 0; iAtom < nAtom; iAtom++) {
            if (!constrainedAtom(iAtom)) continue;
            // get the new state-specific omegas as diagonal elements
            // of rotated omegann matrix
            // rotate from right
            SxMatrix<Complex16> omegannEig = omegaNN(iAtom) ^ eig.vecs;
            // compute only diagonals from rotated from left
            SX_LOOP(in) 
               omega(iAtom)(in,iSpin,ik) = dot(eig.vecs.colRef(in),
                                               omegannEig.colRef (in));
         }
         SX_STOP_TIMER(Timer::OmegaRot);

         if (wavesPtr) (*wavesPtr)(iSpin, ik).rotate (toVector(eig.vecs));
      }//ik
   }//iSpin
   SX_CLOCK(Timer::NuRest);
   // tell every ik after MPI for fermiDistribution
   fermiIn.eps.synMPI();
   // --- get the new fermi occupations for current deltaNu
   fermiIn.fermiDistribution(ektIn);

   SxVector<Double> MSpin(nAtom);
   for (int iAtom = 0; iAtom < nAtom; iAtom++) {
      MSpin(iAtom) = 0;
      if (!constrainedAtom(iAtom)) continue;
      for (int iSpin = 0; iSpin < nSpin; iSpin++) {
         for (int ik = 0; ik < nk; ik++) {
            if (SxLoopMPI::myWork(ik)) {
               for (int in = 0; in < fermiIn.getNStates(ik); ++in)  {
                   MSpin(iAtom) += omega(iAtom)(in,iSpin,ik) * fermiIn.focc(in,iSpin,ik) * fermiIn.kpPtr->weights(ik);
               }//in
            }//LoopMPI
         }//ik
      }//iSpin
   }//iAtom
   SX_MPI_SOURCE ("waves-k", TaskGroupMaster);
   SX_MPI_TARGET (TopLevel, TaskGroupAll);
   SxLoopMPI::sum(MSpin);
   //SX_LOOP(iAtom) cout << "MSpin(" << iAtom << ") = " << setprecision(8) << MSpin(iAtom) << endl;

   return MSpin;
}
