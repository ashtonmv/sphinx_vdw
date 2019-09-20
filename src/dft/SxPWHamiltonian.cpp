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

#include <SxPWHamiltonian.h>
#include <SxRadBasis.h>
#include <SxAtomicOrbitals.h>
#include <SxQuantumNumbers.h>
#include <SxProjector.h>
#include <SxConstants.h>
#include <SxDipoleCorrZ.h>
#include <SxNeighbors.h>
#include <SxSpeciesRef.h>
#include <SxYlm.h>

SxPWHamiltonian::SxPWHamiltonian ()
   : SxHamiltonian (),
     wavesPtr(NULL)
{
   xcPtr = SxPtr<SxXC>::create ();
   xcMeshDensity = 1;
   ekt = 0.;
   eKin = eHartreeElec = eHartreeGauss = eHartree = eLocPs = eNl
        = eEwald = eSelf = eTotal = eDoubleCounting = 0.;
   contrib = CALC_DEFAULT;
   nlBlockSize = -1;
   dipoleCorrection = false;
   zField = 0.;
   zAlignLast = 0.;
   printHartree=true;
}

SxPWHamiltonian::SxPWHamiltonian (const SxPWSet &wavesIn,
                                  const RhoR &rhoRIn,
                                  const SxPseudoPot &psPotIn,
                                  const SxAtomicStructure &str,
                                  bool  rsProjNew)
   : SxHamiltonian (),
     rho (rhoRIn),
     wavesPtr (&wavesIn),
     psPot (psPotIn)
{
   xcPtr = SxPtr<SxXC>::create(psPotIn, &rho.rBasisPtr->getGBasis (), 
                               int(rhoRIn.getSize()));
   xcMeshDensity = 1;
   ekt = 0.;
   eKin = eHartreeElec = eHartreeGauss = eHartree = eLocPs = eNl
        = eEwald = eSelf = eTotal = eDoubleCounting = 0.;
   ekt         = 0.;

   nlBlockSize = -1;

   contrib     = CALC_DEFAULT;
   gBasisPtr   = &(rho.rBasisPtr->getGBasis ());

   structure = str;
   
   eTotal = eKin = eHartree = eExt = eLocPs = eNl = eEwald = eSelf = 0.;

   const SxGBasis &G   = *gBasisPtr;
   const SxGkBasis &Gk =  wavesIn.getGkBasis();

   printHartree = true;
   dipoleCorrection = false;
   zField = 0.;
   zAlignLast = 0.;
   calcForces = false;
   
   computePhiGauss (G);
   computeESelf ();
   computePhiLocPseudo (G);
   computePhiNonLocal (Gk, rsProjNew);
   
   // --- resize vEffR
   int nSpin = (int)rhoRIn.getSize();
   xcPtr->nSpin = nSpin;
   size_t nElem = rhoRIn(0).getSize ();
   vEffR.resize (nSpin);
   for (int iSpin=0; iSpin < nSpin; iSpin++)
      vEffR(iSpin).resize (nElem);

}

SxPWHamiltonian::SxPWHamiltonian (const SxGBasis  &G,
                                  const SxGkBasis &Gk,
                                  const SxPtr<SxPseudoPot> &potPtr,
                                  const SxString  &rhoFile,
                                  const SxSymbolTable *table,
                                  const SxPtr<SxPWSet> &wavesPtrIn)
   : wavesPtr(wavesPtrIn ? wavesPtrIn.getPtr () : NULL),
     psPot (*potPtr),
     gBasisPtr (&G),
     xcPtr (SxPtr<SxXC>::create ())
{
   xcMeshDensity = 1;
   ekt = 0.;
   eKin = eHartreeElec = eHartreeGauss = eHartree = eLocPs = eNl
        = eEwald = eSelf = eTotal = eDoubleCounting = 0.;
   contrib = CALC_DEFAULT;
   nlBlockSize = -1;
   dipoleCorrection = false;
   zField = 0.;
   zAlignLast = 0.;
   printHartree=true;
   SX_NO_MPI;
   // read XC
   read(table);
   // set structure
   SX_CHECK (G.structPtr);
   structure = *G.structPtr;
   const SxRBasis &R = G.getRBasis ();
   // read Rho
   rho.rBasisPtr = &R;
   rho.readRho (rhoFile);

   // --- setup necessary parts of Hamiltonian
   int nSpin = SxHamiltonian::getNSpin (table);
   vEffR.resize (nSpin);
   for (int iSpin = 0; iSpin < nSpin; ++iSpin)
      vEffR(iSpin).resize (R.fft3d.mesh.getSize ());

   computePhiLocPseudo(G);
   computePhiNonLocal (Gk);
   computePhiGauss (G);
   computeRhoGauss (G);

   xcPtr->nSpin = nSpin;
   xcPtr->init();

   if (psPot.nlcc) xcPtr->computePhiCore (psPot, &G);
   if (contrib & CALC_X) xcPtr->enableExchange ();
   if (contrib & CALC_C) xcPtr->enableCorrelation ();

   Contrib fullContrib = contrib;
   contrib = CALC_EFF;
   calcForces = false;
   compute (SxFermi(), true);
   contrib = fullContrib;
}


SxPWHamiltonian::~SxPWHamiltonian ()
{

}

void SxPWHamiltonian::setXCMeshDensity(int in)
{
   xcMeshDensity = in;
}

void SxPWHamiltonian::computeRho (const Focc &focc, const SxPsiSet &psiSet)
{
   const SxPW *waves = dynamic_cast<const SxPW *>(&psiSet);
   SX_CHECK (waves);
   rho.computeRho(focc, *waves);
   rho.normalizeRho ();
}

void SxPWHamiltonian::compute (const SxFermi &fermi,
                               bool tauChanged, bool /*rhoChanged*/)
{
   SX_CHECK (rho.rBasisPtr);
   SX_CHECK (gBasisPtr);
   const SxRBasis    &R = *rho.rBasisPtr;
   const SxGBasis    &G = *gBasisPtr;

   if (tauChanged)  {

      // (0) --- update force symmetrizer

      // (1) --- artificial Gaussian screening charge
      computeRhoGauss (*gBasisPtr);

      // (2) --- local pseudopotential
      SxMeshG vLocG (G.ng, 0.);
      vLocG.setBasis (&G);
      for (int iSpecies=0; iSpecies < structure.getNSpecies (); iSpecies++)
         vLocG += G.structureFactors(iSpecies)
                * phiLocPseudo(iSpecies);
      vLocR = R.symmetrize (vLocG.to(R));

      // (3) --- exchange correlation
      xcPtr->computeXC ();

      // (3) --- Ewald contribution
      if (contrib & CALC_SCR) computeScrPotential ();
   }
   
   update (fermi);
}

void SxPWHamiltonian::update (const SxFermi &fermi)
{
   SX_CLOCK (Timer::hUpdate);

   if (contrib != CALC_DEFAULT && contrib)  printContribMask ();

   int nSpin = (int)rho.rhoR.getSize ();
   SxDiracVec<TPrecForcesG> fArg;

   PrecEnergy eDoubleCountingOld = eDoubleCounting;
   eTotal = eKin = eHartree = eLocPs = eNl = eDoubleCounting = 0.;
   fermiPtr = &fermi;
   if (calcForces)  {
      fNl      = structure.getNewStr ();
      fHartree = structure.getNewStr ();
      fLocPs   = structure.getNewStr ();
      fXC      = structure.getNewStr ();
      fNl.set      (Coord(0.,0.,0.));
      fHartree.set (Coord(0.,0.,0.));
      fLocPs.set   (Coord(0.,0.,0.));
      fXC.set      (Coord(0.,0.,0.));
      if (fScr.getSize () == 0)  {
         fScr = structure.getNewStr ();
         fScr.set (Coord(0., 0., 0.));
      }
   }


   // --- compute total charge density in Fourier space
   SxMeshR rhoRTotal (rho(0).getBasisPtr());
   rhoRTotal.set(0.);
   for (int iSpin=0; iSpin < nSpin; iSpin++) rhoRTotal += rho(iSpin);

   SX_CHECK (rho.rBasisPtr);
   SX_CHECK (gBasisPtr);
   const SxRBasis &R    = *rho.rBasisPtr;
   const SxGBasis &G    = *gBasisPtr;
   SxMeshG rhoG = rhoRTotal.to(G);

   SX_MPI_TARGET (TopLevel, TaskGroupAll);
   SX_MPI_SOURCE(CurrentLevel, TaskGroupMaster);
   // --- k-dependent contributions: kinetic & non-local pp
   if ((contrib & CALC_KIN) || (contrib & CALC_NL)) {
      SX_CHECK (wavesPtr);
      SxLaplacian L;
      const SxPWSet &waves =  getWavesRef ();
      int nk           =  waves.getNk();
      int nStates       =  waves.getNStates ();
      const SxGkBasis &Gk  =  waves.getGkBasis();
      
      SX_MPI_LEVEL("waves-k");
      for (int ik=0; ik < nk; ik++)  {
         if (!SxLoopMPI::myWork(ik)) continue;
         for (int iSpin=0; iSpin < nSpin; iSpin++)  {
            // --- kinetic contribution
            if (contrib & CALC_KIN) {
               SX_CLOCK (Timer::Kin);
               for (int i=0; i < nStates; i++)  {
                  eKin += 0.5 * Gk.weights(ik) * fermi.focc(i,iSpin,ik)
                        * (waves(i,iSpin,ik) | L | waves(i,iSpin,ik) );
               }
            }
            // --- non-local potentials
            if (contrib & CALC_NL)
               computeVNl (iSpin, ik, fermi.focc);
         }
      }
      // MPI summations
      eKin = SxLoopMPI::sum (eKin);
      eNl = SxLoopMPI::sum (eNl);
      SxLoopMPI::sum (fNl);
   }

   // --- Hartree contribution
   if (contrib & CALC_HARTREE)  {
      SX_CLOCK (Timer::Hartree);
      int ng = G.ng - 1;
      SxMeshG rhoHartreeG, vHartreeG(G); vHartreeG(0) = 0.;
      
      rhoHartreeG            = rhoG  - rhoGaussG;
      // sxprintf ("rhoG    : %.12f\n", rhoG(0).re * sqrt(structure.cell.volume));
      // sxprintf ("rhoGauss: %.12f\n", rhoGaussG(0).re * sqrt(structure.cell.volume));
      double charge = rhoHartreeG(0) * sqrt(structure.cell.volume);
      if (fabs(charge) > 1e-6)  {
         cout << "Charge = " << charge << endl;
      } else {
         charge = 0.;
      }
      if (fabs(charge) < 1e-10 || (!dipoleCorrection))
         rhoHartreeG(0) = 0.; // enforce charge neutrality
      vHartreeG(SxIdx(1,ng)) = FOUR_PI * rhoHartreeG(SxIdx(1,ng))
                               / G.g2(SxIdx(1,ng));
      rhoHartreeG.setBasis (gBasisPtr); // TODO: UGLY!!!
      vHartreeR    = vHartreeG.to (R);

      if (calcForces)  {
         SX_LOOP2(is,ia)  {
            fArg = vHartreeG.conj()   // vHartreeG(G=0) = 0!!!
                 * phiGaussG(is)
                 * G.getPhaseFactors((int)is,(int)ia);
            SX_LOOP(d)
               fHartree.ref(is,ia)(d) -= (PrecForcesR)((I*G.gVec.colRef(d)) 
                                                        * fArg).sum();
         }
      }

      double eField = 0.;
      if (dipoleCorrection)  {
         SxDipoleCorrZ dipCorr (rhoRTotal, rhoHartreeG.to(R), charge);
         dipCorr.extraField = zField;
         SxMeshR vDipole(R);
         vDipole.set (0.);
         dipCorr.correctPotential (&vDipole);
         SxRho(vDipole * HA2EV).writeRho ("vDipole.sxb");
         vHartreeR += vDipole;
         if (calcForces) dipCorr.correctForces (&fHartree, structure, psPot);
         PsiG vDipoleG = vDipole.to (G);
         // 1/2 of field energy missing in eHartree due to 1/2 prefactor
         if (fabs(charge) < 1e-6)  {
            eField           = 0.5 * zField * dipCorr.dipole;
            eDoubleCounting += 0.5 * zField * dipCorr.dipole;
         }

         eDoubleCounting += 0.5 * dot(rhoHartreeG, vDipoleG).re;
         eDoubleCounting -= dot(rhoG, vDipoleG).re;

         double vAvgLeft = 0.;
         for (int i = 0; i < R.fft3d.mesh(0); i++)  {
            for (int j = 0; j < R.fft3d.mesh(1); j++)  {
               vAvgLeft += vHartreeR(R.fft3d.mesh.getMeshIdx(i,j,dipCorr.rhoMinZ, SxMesh3D::Positive));
            }
         }
         vAvgLeft /= R.fft3d.mesh(0) * R.fft3d.mesh(1);
         // align at zAlign
         if (dipCorr.dipole * charge > 0.)  {
            vAvgLeft += dipCorr.dipole / charge * dipCorr.extraField;
         } else if (fabs(charge) > 0.)  {
            vAvgLeft += (dipCorr.dipole / charge + R.cell(2,2)) * dipCorr.extraField;
         }

         cout << "vAlign=" << vAvgLeft << endl;
         eField += 0.5 * charge * vAvgLeft;
         eDoubleCounting += 0.5 * charge * vAvgLeft;
         // compute alignment energy for z = 0
         if (fabs(charge) > 0.)  {
            double zAlignNow = dipCorr.getZalign ();
            double c = fabs(R.cell(2,2));
            while (zAlignNow - zAlignLast > 0.8 * c) zAlignNow -= c;
            while (zAlignNow - zAlignLast < -0.8 * c) zAlignNow += c;
            if (fabs(zAlignNow - zAlignLast) > 0.3 * c && zAlignNow != 0.)  {
               cout << "WARNING: big shift in zAlign from " << zAlignLast
                    << " to " << zAlignNow << endl;
               cout << "         The total energy may jump!" << endl;
            }
            zAlignLast = zAlignNow;
            cout << "zAlign=" << zAlignNow << endl;
            double A = structure.cell.volume / c;
            double QL = dipCorr.extraField /(-FOUR_PI) * A;
            double QR = -(QL+charge);
            double UL = -zAlignNow * dipCorr.extraField;
            double fieldFromCharge = - FOUR_PI * charge / A;
            double UR = -zAlignNow * (dipCorr.extraField + fieldFromCharge);
            cout << "QL=" << QL << " QR=" << QR 
                 << " UL=" << (UL * HA2EV) << " eV UR=" << (UR * HA2EV) 
                 << " eV" << endl;
            double eAlign0 = 0.5 * charge * (UR - UL);
            eField += eAlign0;
            eDoubleCounting += eAlign0;
            vAvgLeft += UR;
         }

         // align on zAlign
         vHartreeR -= vAvgLeft;
         eField -= 0.5 * charge * vAvgLeft;
         eDoubleCounting -= charge * vAvgLeft;
         eDoubleCounting += rhoG(0).re * sqrt(structure.cell.volume) * vAvgLeft;
      }

      eHartree     = 0.5 * tr(rhoHartreeG.to(R) * vHartreeR)
                   + eField;
      eHartreeElec = TWO_PI
                   * (rhoG(SxIdx(1,ng)).absSqr() / G.g2(SxIdx(1,ng))).sum();
      eDoubleCounting -= eHartreeElec;
      sxprintf ("# eDC = %20.16f\n", eDoubleCounting);
      eDoubleCounting += TWO_PI * (rhoGaussG(SxIdx(1,ng)).absSqr () / G.g2(SxIdx(1,ng))).sum ();


   }

   // --- local pseudopotential
   if (contrib & CALC_LOC)  {
      SX_CLOCK (Timer::locPs);
      eLocPs = 0.;
      SX_LOOP(iSpin)
         eLocPs += tr(vLocR * rho(iSpin));
      if (calcForces)  {
         SX_LOOP2(is,ia)  {
            fArg = rhoG.conj() * phiLocPseudo(is) 
                 * G.getPhaseFactors((int)is,(int)ia);
            SX_LOOP(d)
               fLocPs.ref(is,ia)(d)
                  += (PrecForcesR)( (I * G.gVec.colRef(d)) * fArg).sum();
         }
      }
   }

   // --- external potential
   SxAtomicStructure fExt;
   if (contrib & CALC_EXT)  {
      eExt = 0.;
      if (!(   xcPtr->xcFunctional == SxXC::EXX
            || xcPtr->xcFunctional == SxXC::EXX_LDA
            || xcPtr->xcFunctional == SxXC::EXX_PBE
            || xcPtr->xcFunctional == SxXC::EXX_PBE_WB
            || xcPtr->xcFunctional == SxXC::Slater
            || xcPtr->xcFunctional == SxXC::KLI))  {
         SX_LOOP(iSpin)
            eExt += tr(vExtR(iSpin) * rho(iSpin));
         if (vExtActsOnNuclei)  {
            // minus-sign: nuclear charge in e-units (cf. rhoHartreeG)
            double eExtNuc = -tr(vExtR(0) * rhoGaussG.to(R).real ());
            eExt += eExtNuc;
            eDoubleCounting += eExtNuc;
         }
      }  else  {
         eExt = xcPtr->eXc;
      }
      if (calcForces && vExtActsOnNuclei)  {
         if (fExt.getSize () == 0) fExt = structure.getNewStr ();
         fExt.set (0., 0., 0.);
         PsiG vExtG = vExtR(0).to(G);
         SX_LOOP2(is,ia)  {
            fArg = vExtG.conj() * phiGaussG(is) 
                 * G.getPhaseFactors((int)is,(int)ia);
            SX_LOOP(d)
               //                    (-iG * fArg ).re
               fExt.ref(is,ia)(d) += (G.gVec.colRef(d)*fArg).sum().im;
         }
      }
   }

   // --- exchange-correlation
   if (contrib & (CALC_X | CALC_C))  {
      SX_CLOCK (Timer::XC);
	   if (xcMeshDensity >= 2) {
			cout << "Averaging XC-Potential ..." << endl;
         
         const SxPWSet &waves =  getWavesRef();
         double eXc1 = 0.;
			int iG;
			double sp;
			SxComplex16 phase;
			SxVector<Complex16> phaseVec(G.ng);
			SxVector<Double> Rvec(3);
			SxVector<Double> Gvec(3);
			SxList<SxMeshR> vXc1; vXc1.resize(nSpin);
			SxList<SxMeshG> vXc1G; vXc1G.resize(nSpin);
			SxList<SxMeshG> vXcG; vXcG.resize(nSpin);
			SxList<SxMeshG> myrhoG; myrhoG.resize(nSpin);
			SxList<SxMeshG> rhoGorig; rhoGorig.resize(nSpin);
			
			double length 
				= (structure.cell (0,0) / R.fft3d.mesh(0) /(double)xcMeshDensity);
			
         for (int iSpin = 0; iSpin < nSpin; iSpin++) {
				vXc1(iSpin).resize (xcPtr -> vXc(iSpin).getSize ()); 
            vXc1(iSpin).setBasis(&R);vXc1(iSpin).set (0.);
            vXc1G(iSpin) = SxMeshG(G.ng); vXc1G(iSpin).setBasis(&G);
            myrhoG(iSpin) = SxMeshG(G.ng); myrhoG(iSpin).setBasis(&G);
            rhoGorig(iSpin) = SxMeshG(G.ng); rhoGorig(iSpin).setBasis(&G);
         }
		
			
			for (int iSpin = 0; iSpin < nSpin; iSpin++) {
				rhoGorig(iSpin) = rho.rhoR(iSpin).to(G);
			}

			
			for (int x = 0; x < xcMeshDensity; x++) {
				for (int y = 0; y < xcMeshDensity; y++) {
					for (int z = 0; z < xcMeshDensity; z++) {
					
						Rvec(0) = length*(double)x;
						Rvec(1) = length*(double)y;
						Rvec(2) = length*(double)z;
						
						for (iG = 0; iG < G.ng; iG++) { 
								Gvec(0) =  -G.gVec.row(iG)(0);
								Gvec(1) =  -G.gVec.row(iG)(1);
								Gvec(2) =  -G.gVec.row(iG)(2);
                        sp = (Gvec*Rvec). sum ();
								phase.re = cos(sp); phase.im = sin(sp);
                        phaseVec(iG) = phase;
						}

						
						for (int iSpin = 0; iSpin < nSpin; iSpin++) {
                     for (int i = 0; i < myrhoG(iSpin).getSize (); i++) {
                        myrhoG(iSpin)(i) = phaseVec(i)*rhoGorig(iSpin)(i);
                     }
                     rho.rhoR(iSpin) = myrhoG(iSpin).to(R);
						}
                 
                  xcPtr->updateXC (rho.rhoR, wavesPtr, fermiPtr);
				
						eXc1 += xcPtr->eXc/(double)(pow((double)xcMeshDensity, 3.));
						
                  for (iG = 0; iG < G.ng; iG++) { 
								Gvec(0) =  G.gVec.row(iG)(0);
								Gvec(1) =  G.gVec.row(iG)(1);
								Gvec(2) =  G.gVec.row(iG)(2);
								sp = Gvec(0)*Rvec(0) + Gvec(1)*Rvec(1) + Gvec(2)*Rvec(2);
								phase.re = cos(sp);
								phase.im = sin(sp);
								phaseVec(iG) = phase;
						}
						
				
						for (int iSpin = 0; iSpin < nSpin; iSpin++) {
							vXc1G(iSpin) = xcPtr -> vXc(iSpin).to(G);
                     for (int i = 0; i < vXc1G(iSpin).getSize (); i++) {
                        vXc1G(iSpin)(i) = vXc1G(iSpin)(i)*phaseVec(i);
                     }
							xcPtr -> vXc(iSpin) = vXc1G(iSpin).to(R);
						}
                  
						for (int iSpin = 0; iSpin < nSpin; iSpin++) {
                     for (int i = 0; i < vXc1(iSpin).getSize (); i++) {
                        vXc1(iSpin)(i) = vXc1(iSpin)(i) + (xcPtr -> vXc(iSpin))(i)/pow((double)xcMeshDensity, 3.);
                     }
						}
                 
					}
				}	
				
			}	

			xcPtr -> eXc = eXc1;

		
			for (int i = 0; i < nSpin; i++) {
            for (int j = 0; j < vXc1(i).getSize(); j++) {
				xcPtr -> vXc(i)(j) = vXc1(i)(j);
            }
			}

         for (int iSpin = 0; iSpin < waves.getNSpin (); iSpin++) {
            rho.rhoR(iSpin) = rhoGorig(iSpin).to(R);
         }
		cout << "...ready" << endl;
		}       
      else {
         xcPtr->updateXC (rho.rhoR, wavesPtr, fermiPtr);
      }
      if (calcForces && xcPtr->nlcc)  {
         SxMeshR vXcTotal = (nSpin == 1) ? xcPtr->vXc(0)
                                         : xcPtr->vXc(0) + xcPtr->vXc(1);
         PsiG vXcG = vXcTotal.to (G);

         SX_LOOP(is)  {
            if (xcPtr->phiCore(is).getSize () > 0)  {
               SX_LOOP(ia)  {
                  fArg = xcPtr->phiCore(is)
                       * vXcG.conj()
                       * G.getPhaseFactors((int)is,(int)ia);
                  SX_LOOP(d)
                     fXC.ref(is,ia)(d)
                        += (PrecForcesR)( (I * G.gVec.colRef(d)) * fArg).sum();
               }
            }
         }
      }
      eDoubleCounting += xcPtr->eXc;
      SX_LOOP(iSpin)
         eDoubleCounting -= dot(rho.rhoR(iSpin),xcPtr->vXc(iSpin)) * R.dOmega;
   }

   // --- assemble total energy
   eTotal = (contrib & CALC_KIN     ? eKin       : 0.)
          + (contrib & CALC_HARTREE ? eHartree   : 0.)
          + (contrib & CALC_XC      ? xcPtr->eXc : 0.)
          + (contrib & CALC_LOC     ? eLocPs     : 0.)
          + (contrib & CALC_EXT     ? eExt       : 0.)
          + (contrib & CALC_NL      ? eNl        : 0.)
          +  eEwald - eSelf;

   if (applyVDWCorrection) {
      vdwCorrection.compute();
      eVDW = vdwCorrection.totalEnergy;
      eTotal = eTotal + eVDW;
   }
   eDoubleCounting += eEwald - eSelf;

   if (wavesPtr && fermi.getNk () > 0) {
      PrecEnergy eHarrisFoulkes = 0.;
      const SxGkBasis &Gk  =  wavesPtr->getGkBasis();
      SX_LOOP3(ik,iSpin,i) 
         eHarrisFoulkes += fermi.focc(i,iSpin,ik) * fermi.eps(i,iSpin,ik) 
                         * Gk.weights(ik);
      eHarrisFoulkes += eDoubleCountingOld;
      sxprintf ("# eTotal = %20.16f\n", eTotal);
      sxprintf ("# eHarrisFoulkes = %20.16f\n", eHarrisFoulkes);
   }
   
   // --- assemble effective potential
   for (int iSpin=0; iSpin < nSpin; iSpin++)  {
      vEffR(iSpin) = 0.;
      vEffR(iSpin).setBasis (rho(0).getBasisPtr()); // TODO: ugly
      if (contrib & CALC_HARTREE)  vEffR(iSpin) += vHartreeR;
      if (contrib & CALC_LOC)      vEffR(iSpin) += vLocR;
      if (contrib & CALC_XC)       vEffR(iSpin) += xcPtr->vXc(iSpin);
      if (contrib & CALC_EXT)      vEffR(iSpin) += vExtR(iSpin);
      vEffR(iSpin) = R.symmetrize (vEffR(iSpin));
   }

   // --- assemble total forces
   if (calcForces)  {
      fTotal = fHartree + fLocPs + fNl + fXC + fScr;
      if (contrib & CALC_EXT && vExtActsOnNuclei) fTotal += fExt;
      cout << fTotal << endl;
   }
}



PsiG SxPWHamiltonian::operator* (const PsiG &psiG) const
{
   SX_CLOCK (Timer::hPsi);

   // --- TODO: ugly!!!
   int i       =  psiG.handle->auxData.i;

   if (i != -1)   return H_PsiG  (psiG); // apply H to only single state 
   else           return H_PsiGI (psiG); // apply H to all states at once
}

PsiG SxPWHamiltonian::H_PsiG (const PsiG &psiG) const
{
   SX_CHECK (psiG.getBasisPtr());
   SX_CHECK (dynamic_cast<const SxGBasis *>(psiG.getBasisPtr()));
   SX_CHECK (rho.rBasisPtr);

   const SxGBasis &Gk = *dynamic_cast<const SxGBasis *>(psiG.getBasisPtr());
   const SxRBasis &R  = *rho.rBasisPtr;

   // --- TODO: ugly!!!
   int i       =  psiG.handle->auxData.i;
   int iSpin   =  psiG.handle->auxData.iSpin;
   int ik      =  psiG.handle->auxData.ik;

   PsiG dPsiG (psiG.getSize(), 0.);
   dPsiG.setBasis (&Gk);

   // --- TODO: ugly!!!
   dPsiG.handle->auxData.i     = i;
   dPsiG.handle->auxData.iSpin = iSpin;
   dPsiG.handle->auxData.ik    = ik;

   if (contrib & CALC_KIN)   {
      SX_CLOCK (Timer::tPsi);
      dPsiG += (0.5 * Gk.g2) * psiG;
   }
   if (   (contrib & CALC_HARTREE) || (contrib & CALC_XC)
       || (contrib & CALC_LOC)     || (contrib & CALC_RHO)
       || (contrib & CALC_EXT) )
   {
      SX_CLOCK (Timer::vEffPsi);
      dPsiG += (Gk | (vEffR(iSpin) * (R | psiG)) );
   }

   if (contrib & CALC_NL && phiOrbNl.getSize () > 0 && !rsProj2)    {
      SX_CLOCK (Timer::vNlPsi);
      PsiG pPsi = *nlProj | psiG;
      SX_LOOP(io) pPsi(io) /= eKB(phiOrbNl(io));
      dPsiG += Gk | pPsi;
   }
   
   if (contrib & CALC_NL)    {
      // --- EES-G real-space method
      SxVector<TPrecCoeffG> phiPsi;
      int is;
      // loop over reference orbitals
      for (int iPhi = 0; iPhi < phiNl(ik).getSize (); ++iPhi)  {
         is = phiNl(ik)(iPhi).handle->auxData.is;
         SX_CHECK (is >= 0 && is < structure.getNSpecies (),
                   is, structure.getNSpecies ());
         if (psPot.realSpace(is))  {
            SxAtomicStructure atoms = SxSpeciesRef(is) | structure; 
            phiPsi = rsProjector->project (phiNl(ik)(iPhi), psiG, atoms);
            dPsiG += rsProjector->gradient (phiPsi / eKB(iPhi), phiNl(ik)(iPhi),
                                            atoms);
         }
      }
      // --- new real-space projectors
      if (rsProj2)  {
         SxDiracVec<TPrecCoeffG> pPsi;
         pPsi = rsProj2->project (psiG);
         for (int iPhi = 0; iPhi < pPsi.getSize (); ++iPhi)
            pPsi(iPhi) /= eKB(phiOrbNl(iPhi));
         dPsiG += rsProj2->gradient (pPsi, Gk);
      }
   }
   
   return dPsiG;
}


PsiG SxPWHamiltonian::H_PsiGI (const PsiG &psiGI) const
{
   SX_CHECK (dynamic_cast<const SxGBasis *>(psiGI.getBasisPtr()));
   SX_CHECK (rho.rBasisPtr);

   int i=0, nStates = (int)psiGI.nCols();
   const SxGBasis &Gk = *dynamic_cast<const SxGBasis *>(psiGI.getBasisPtr());
   const SxRBasis &R  = *rho.rBasisPtr;

   // --- TODO: ugly!!!
   int iSpin   =  psiGI.handle->auxData.iSpin;
   int ik      =  psiGI.handle->auxData.ik;

   SxDiracMat<TPrecCoeffG> dPsiGI (psiGI.nRows(), psiGI.nCols());
   dPsiGI.handle->auxData = psiGI.handle->auxData;
   dPsiGI.set (0.);
   PsiG dPsiG, psiG;
#ifdef USE_FFT2d1d
   bool directVeff = true;
#else
   bool directVeff = false;
#endif
   if (directVeff && 
      (contrib & (CALC_HARTREE | CALC_XC | CALC_LOC | CALC_RHO | CALC_EXT)))
   {
      SX_CLOCK (Timer::vEffPsi);
      dPsiGI += Gk.convolute (psiGI, vEffR(iSpin));
   }

   for (i=0; i < nStates; i++)  {
      dPsiG = dPsiGI.colRef(i);
      psiG  = psiGI.colRef(i);

      // --- TODO: ugly!!!
      //psiG.handle->auxData.i     = i;
      //psiG.handle->auxData.iSpin = iSpin;
      //psiG.handle->auxData.ik    = ik;

      if (contrib & CALC_KIN)   {
         SX_CLOCK (Timer::tPsi);
         dPsiG += (0.5 * Gk.g2) * psiG;
      }
      if ((   (contrib & CALC_HARTREE) || (contrib & CALC_XC)
           || (contrib & CALC_LOC)     || (contrib & CALC_RHO)
           || (contrib & CALC_EXT) ) && (!directVeff))
      {
         SX_CLOCK (Timer::vEffPsi);
         dPsiG += (Gk | (vEffR(iSpin) * (R | psiG)) );
      }
      if (contrib & CALC_NL)    {
         // --- EES-G real-space method
         SxVector<TPrecCoeffG> phiPsi;
         int is;
         // loop over reference orbitals
         for (int iPhi = 0; iPhi < phiNl(ik).getSize (); ++iPhi)  {
            is = phiNl(ik)(iPhi).handle->auxData.is;
            SX_CHECK (is >= 0 && is < structure.getNSpecies (),
                      is, structure.getNSpecies ());
            if (psPot.realSpace(is))  {
               SxAtomicStructure atoms = SxSpeciesRef(is) | structure; 
               phiPsi = rsProjector->project (phiNl(ik)(iPhi), psiG, atoms);
               dPsiG += rsProjector->gradient (phiPsi / eKB(iPhi),
                                               phiNl(ik)(iPhi), atoms);
            }
         }
      }
   }

   if ((contrib & CALC_NL) && phiOrbNl.getSize () > 0 && !rsProj2)  {
      SX_CLOCK (Timer::vNlPsi);
      PsiG pPsi = *nlProj | psiGI;
      SX_LOOP2(iState, io) pPsi(io,iState) /= eKB(phiOrbNl(io));
      dPsiGI += Gk | pPsi;
   }

   // --- new real-space projectors
   if (rsProj2)  {
      int nPhi = (int)phiOrbNl.getSize ();
      SxDiracVec<TPrecCoeffG> pPsi;
      pPsi = rsProj2->project (psiGI);
      for (i=0; i < nStates; i++)  {
         for (int iPhi = 0; iPhi < nPhi; ++iPhi)
            pPsi(iPhi,i) /= eKB(phiOrbNl(iPhi));
         //dPsiGI.colRef(i) += rsProj2->gradient (pPsi.colRef(i), Gk);
      }
      dPsiGI += rsProj2->gradient (pPsi, Gk);
   }

   return dPsiGI;
}


SxDiracVec<TPrecCoeffG::TReal> 
SxPWHamiltonian::preconditioner (const PsiG &dPsi,
                                 Preconditioner type) const
{
   const SxGBasis *Gk = dynamic_cast<const SxGBasis *>(dPsi.getBasisPtr());
   SX_CHECK (Gk);
   
   SxDiracVec<TPrecCoeffG::TReal> x, x2, x3, x4, n, K;

   if (type == Payne)  {
      // M. P. Teter, M. C. Payne, and D. C. Allan.  Solution of Schrodinger’s
      // equation for large systems.  Phys. Rev. B, 40, 12255 (1989).
      int nCol = (int)dPsi.nCols ();
      SX_CHECK (nCol > 0, nCol);
      double kin = 0.;
      SxLaplacian L;
      // average kinetic energy over all states (usually 1)
      for (int i = 0; i < nCol; ++i)
         kin += (dPsi.colRef (i) | L | dPsi.colRef(i) );
      kin /= nCol;

      //kin *= 0.5;

      /*
      x  = Gk->g2 / kin;
      x2 = x.sqr();
      x3 = x.cub();
      x4 = x2.sqr();
      n  = 27. + 18.*x + 12.*x2 + 8.*x3;
      K  = n / (n + 16.*x4);
      K /= kin; // make K approach 1/(G^2 - kin)
      */
      double kinInv = 1. / kin;
      K.resize (Gk->ng);
      K.setBasis (Gk);
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for (ssize_t ig = 0; ig < Gk->ng ; ++ig)   {
         double x_ = Gk->g2(ig) * kinInv;
         double x2_ = x_ * x_;
         double x3_ = x_ * x2_;
         double x4_ = x2_ * x2_;
         double n_ = 27. + 18. * x_ + 12. * x2_ + 8. * x3_;
         K(ig) = n_ / (kin * (n_ + 16. * x4_));
      }
      return K;
   }  else  {
      const SxPWSet &waves =  getWavesRef();
      int nStates =  waves.getNStates();
      x = 0.5 * Gk->g2 / (0.5 * eKin / nStates);
      n = 1.+x*(1.+x*(1.+x*(1.+x*(1.+x*(1.+x*(1.+x*(1.+x)))))));
      return n / (1. + x*n);
   }
}




// ==========================================================================
//    H A R T R E E  contribution
// ==========================================================================
void SxPWHamiltonian::computePhiGauss (const SxGBasis &G)
{
   SX_CLOCK (Timer::phiGauss);

   const int s = SxQuantumNumbers::s;
   Real8 zv, rGauss;

   SxRadBasis          r (psPot.rad, psPot.logDr);  // |r>
   SxDiracVec<Double>  func, rad;

   phiGaussG.resize (structure.getNSpecies ());
   for (int iSpecies = 0; iSpecies < structure.getNSpecies (); iSpecies++)  {
      // --- numerical transform
      zv     = psPot.valenceCharge(iSpecies);
      rGauss = psPot.rGauss(iSpecies);
      rad    = SxDiracVec<Double> (psPot.rad(iSpecies));

      func   = exp ( -(rad / rGauss).sqr() );
      func  *= zv / (sqrt(PI*PI*PI) * rGauss*rGauss*rGauss );
      func  /= SxYlm::getYlmNormFactor(s,0);

      // TODO: --- ugly!!!
      func.handle->auxData.is    = iSpecies;
      func.handle->auxData.l     = s;
      func.handle->auxData.m     = 0;
      func.setBasis (&r);
      
      //phiGaussG(iSpecies) = SUM (r, ( G | r ) * ( r | func ));
      // --- analytic transform
      double volume = structure.cell.volume;
      phiGaussG(iSpecies) = zv/sqrt(volume) * exp(-0.25 * rGauss*rGauss * G.g2);
   }
}


void SxPWHamiltonian::computeRhoGauss (const SxGBasis &G)
{
   rhoGaussG = SxMeshG (G.ng, 0.);
   rhoGaussG.setBasis (&G);
   for (int iSpecies=0; iSpecies < structure.getNSpecies (); iSpecies++)  {
      rhoGaussG += G.structureFactors(iSpecies)
                 * phiGaussG(iSpecies);
   }
   const SxRBasis &R = G.getRBasis ();
   rhoGaussG = G | (R.symmetrize (R | rhoGaussG));

   int ng = G.ng - 1;
   eHartreeGauss = TWO_PI
                 * ( rhoGaussG(SxIdx(1,ng)).absSqr() / G.g2(SxIdx(1,ng)) ).sum();
}


void SxPWHamiltonian::computeESelf ()
{
   eSelf = 0.;
   for (int iSpecies=0; iSpecies < structure.getNSpecies (); iSpecies++)  {
      eSelf += (psPot.valenceCharge(iSpecies) * psPot.valenceCharge(iSpecies))
             *  structure.getNAtoms(iSpecies) 
             /  psPot.rGauss(iSpecies);
   }

   eSelf /= sqrt (TWO_PI);
}


// ==========================================================================
//   L O C A L   P S E U D O P O T E N T I A L   contribution
// ==========================================================================
void SxPWHamiltonian::computePhiLocPseudo (const SxGBasis &G)
{
   SX_CLOCK (Timer::phiLocPs);
   double zv, rGauss, omega = structure.cell.volume;
   SxRadBasis          r (psPot.rad, psPot.logDr);  // |r>
   SxDiracVec<Double>  func, rad, pot;
   
   int l;

   const int s = SxQuantumNumbers::s;

   phiLocPseudo.resize(structure.getNSpecies ());
   double avgV = 0;;
   cout <<SX_SEPARATOR;
   for (int iSpecies=0; iSpecies < structure.getNSpecies (); iSpecies++)  {
      l      = psPot.lLoc(iSpecies);
      zv     = psPot.valenceCharge(iSpecies);
      rGauss = psPot.rGauss(iSpecies);

      rad    = SxDiracVec<Double> (psPot.rad(iSpecies));
      pot    = SxDiracVec<Double> (psPot.pseudoPot(iSpecies)(l));

      func   = pot + zv / rad * derf (rad / rGauss);
      func  /= SxYlm::getYlmNormFactor(s,0);  // single atoms are spherical (s-like)

      // TODO: --- ugly!!!
      func.handle->auxData.is    = iSpecies;
      func.handle->auxData.l     = s;
      func.handle->auxData.m     = 0;
      func.setBasis (&r);
      
      phiLocPseudo(iSpecies) = SUM (r, ( G | r ) * ( r | func ) );

      cout << "| " << psPot.chemName(iSpecies) << ": "
           << "Delta V Omega/atom = " 
           << (phiLocPseudo(iSpecies)(0) * sqrt(omega) * HA2EV)
           << " eV bohr^3" << endl;
      avgV += phiLocPseudo(iSpecies)(0) * structure.getNAtoms(iSpecies);
   }
   sxprintf ("| Average electrostatic potential [eV]: %.12f\n",
             avgV / sqrt(omega) * HA2EV);
   cout << SX_SEPARATOR;
}


// ==========================================================================
//   N  O N  -  L O C A L   P S E U D O P O T E N T I A L  contribution
// ==========================================================================
void SxPWHamiltonian::computePhiNonLocal (const SxGkBasis &Gk, bool rsProjNew)
{
   SX_CLOCK (Timer::phiNonLoc);

   int ik, iSpecies, iAtom, l, m, phiIdx;
   double diag, lDr;
   SxPtr<SxRadBasis> rPtr = SxPtr<SxRadBasis>::create(psPot.rad, psPot.logDr); // |r>
   SxRadBasis &r = *rPtr;
   SxAtomicOrbitals          orbitals (psPot.pseudoPsi, rPtr);           // |mu>
   SxDiracVec<Double>        psi_pot, rad, psi, pot;                  // :r
   SxList<int>               phiOrbList;                 //:{is,n,l,m,ik}
   SxList<SxQuantumNumbers>  psiOrbList;               //:{is,ia,n,l,m,ik}
   double unit = printHartree ? 1.  :  HA2EV;
   SxString u  = printHartree ? "H" : "eV";
   bool rsProjEESG = false;

   // --- create lookup for {is,ia,n,l,m,ik} |-> {is,n,l,m,ik}
   int nSpecies = structure.getNSpecies ();
   int phiIdxStart;
   for (iSpecies=0, phiIdx=0; iSpecies < nSpecies; iSpecies++) {
      phiIdxStart = phiIdx;
      for (iAtom=0;iAtom<structure.getNAtoms(iSpecies);iAtom++) {
         phiIdx = phiIdxStart;
         for (l = SxQuantumNumbers::s; l <= psPot.lMax(iSpecies); l++)  {
            if (l != psPot.lLoc(iSpecies))  {
               for (m=-l; m <= l; m++, phiIdx++)  {
                  if (!psPot.realSpace(iSpecies))  {
                     phiOrbList << phiIdx;
                     psiOrbList << SxQuantumNumbers (iSpecies, iAtom, -1, l,
                                                     m, -1, -1);
                  } else {
                     rsProjEESG = true;
                  }
               }
            }
         }
      }
   }
   if (rsProjNew) rsProjEESG = false;

   phiOrbNl = phiOrbList;
   psiOrbNl = psiOrbList;

   // --- memory for phiNl and eKB
   eKB.resize   (phiIdx);
   phiNl.resize (Gk.nk);
   if (rsProjEESG) {
      SX_MPI_LEVEL("waves-k");
      for (ik=0; ik < Gk.nk; ik++)  {
         if (SxLoopMPI::myWork (ik))
            phiNl(ik).resize (phiIdx);
      }
   }

   cout << SX_SEPARATOR;
   cout << "| Kleinman-Bylander energies:\n";
   cout << SX_SEPARATOR;

   SxArray<SxArray<SxDiracVec<TReal8> > > projR(nSpecies);

   // --- compute eKB and phiNl
   phiIdx = 0;
   for (iSpecies=0; iSpecies < psPot.getNSpecies (); iSpecies++)  {
      cout << "| " << psPot.chemName(iSpecies) << ": ";
      if (!psPot.realSpace(iSpecies))  {
         int nProj = 0;
         for (l = SxQuantumNumbers::s; l <= psPot.lMax(iSpecies); l++)
            if ( l != psPot.lLoc(iSpecies) ) nProj++;
         projR(iSpecies).resize (nProj);
      }
      int iProj = 0;
      for (l = SxQuantumNumbers::s; l <= psPot.lMax(iSpecies); l++)  {
         if ( l != psPot.lLoc(iSpecies) )  {

            // --- retrieve pseudopotential information
            rad = psPot.rad(iSpecies);
            psi = psPot.pseudoPsi(iSpecies)(l);
            pot = psPot.pseudoPot(iSpecies)(l);
            lDr = psPot.logDr(iSpecies);

            // --- compute diagonal elements (Kleinman-Bylander energy)
            diag = (psi.sqr() * rad.cub() * pot).integrate(lDr);
            cout << "  " << SxQuantumNumbers::lToChar(l) << "=" 
                         << unit * diag << " " << u;
            
            // TODO: ugly!!!
            psi_pot = psi * pot;
            psi_pot.setBasis (&r);
            psi_pot.handle->auxData.is = iSpecies;
            psi_pot.handle->auxData.l  = l;
            if (!psPot.realSpace(iSpecies))
               projR(iSpecies)(iProj++) = psi_pot;
            for (m=-l; m <= l; m++)  {
               eKB(phiIdx) = diag;

               if (rsProjEESG)  {
                  SX_MPI_LEVEL("waves-k");
                  // TODO: ugly!!!
                  psi_pot.handle->auxData.m  = m;
                  for (ik=0; ik < Gk.nk; ik++)  {
                     if (SxLoopMPI::myWork(ik))
                        phiNl(ik)(phiIdx) = Gk(ik) | psi_pot;
                  }
               }
               phiIdx++;
            }
         }  else  {
            cout << "  " << SxQuantumNumbers::lToChar(l) << "=local";
         }
      }
      cout << endl;
   }
   cout << SX_SEPARATOR;
   nlProj = SxPtr<SxAOBasis>::create (Gk, *rPtr, projR);
}



void SxPWHamiltonian::computeVNl (int iSpin, int ik, const Focc &focc)
{
   SX_CHECK (wavesPtr);
   SX_CLOCK (Timer::nonLoc);

   const SxPWSet &waves =  getWavesRef();
   bool allStatesAvailable = dynamic_cast<const SxPW *>(&waves); 
   const SxGkBasis &Gk  =  waves.getGkBasis();
   int i, nStates       =  waves.getNStates();
   SxDiracVec<TPrecCoeffG> allPotPsiFac;
   double nlFac;

   int nNlOrb  = (int)psiOrbNl.getSize();
   if (nNlOrb == 0 && !psPot.useRealSpace ()) return;
   if (! calcForces)  {
      
      if (nNlOrb > 0 && !rsProj2)  {
         // --- get projections < mu | Vps | Psi > 
         if (allStatesAvailable)  {
            const SxPW &wavesRef = (const SxPW &)waves;
            allPotPsiFac = *nlProj | wavesRef(iSpin, ik);
         } else {
            // NOTE: inefficient, should do blocking
            allPotPsiFac.reformat (nNlOrb, nStates);
            for (i = 0; i < nStates; ++i)
               allPotPsiFac.colRef (i) = *nlProj | waves(i,iSpin,ik);
         }
         // --- get nonlocal energy
         // 
         //               |< mu | Vps | Psi_n >|^2 
         // E  =  Sum     -------------------------  * focc(n)
         //                  < mu | Vps | mu >
         for (int iNlOrb=0; iNlOrb < nNlOrb; iNlOrb++)  {  //:is,ia,n,l,m

            nlFac = Gk.weights(ik) / eKB(phiOrbNl(iNlOrb));
            eNl += nlFac * (allPotPsiFac.row(iNlOrb).absSqr ()
                            * focc(iSpin, ik)                 ).sum ();
         }
      }

      if (psPot.useRealSpace ()) {
         // --- EES-G real-space method
         SxVector<TPrecCoeffG> phiPsi;
         int is;
         // loop over reference orbitals
         for (int iPhi = 0; iPhi < phiNl(ik).getSize (); ++iPhi)  {
            is = phiNl(ik)(iPhi).handle->auxData.is;
            SX_CHECK (is >= 0 && is < structure.getNSpecies (),
                      is, structure.getNSpecies ());
            if (psPot.realSpace(is))  {
               SxAtomicStructure atoms = SxSpeciesRef(is) | structure;
               for (i = 0; i < nStates; ++i)  {
                  // skip unoccupied states
                  if (fabs(focc(i, iSpin, ik)) < 1e-12) continue;
                  // --- get projections < mu | Vps | Psi >
                  phiPsi = rsProjector->project (phiNl(ik)(iPhi),
                                                 waves(i,iSpin,ik),
                                                 atoms);
                  // --- get nonlocal energy (see above)
                  eNl += Gk.weights(ik) * focc(i, iSpin, ik)
                         * phiPsi.absSqr ().sum () / eKB (iPhi);
               }
            }
         }
      }
      // --- new real-space projectors
      if (rsProj2)  {
         SxDiracVec<TPrecCoeffG> pPsi;
         if (allStatesAvailable)  {
            pPsi = rsProj2->project (waves(iSpin, ik));
            int nPhi = (int)phiOrbNl.getSize ();
            for (i = 0; i < nStates; ++i)  {
               // skip unoccupied states
               if (fabs(focc(i, iSpin, ik)) < 1e-12) continue;
               for (int iPhi = 0; iPhi < nPhi; ++iPhi)  {
                  // --- get nonlocal energy (see above)
                  eNl += Gk.weights(ik) * focc(i, iSpin, ik)
                         * pPsi(iPhi, i).absSqr () / eKB (phiOrbNl(iPhi));
               }
            }
         } else {
            for (i = 0; i < nStates; ++i)  {
               // skip unoccupied states
               if (fabs(focc(i, iSpin, ik)) < 1e-12) continue;
               pPsi = rsProj2->project (waves(i, iSpin, ik));
               for (int iPhi = 0; iPhi < pPsi.getSize (); ++iPhi)  {
                  // --- get nonlocal energy (see above)
                  eNl += Gk.weights(ik) * focc(i, iSpin, ik)
                         * pPsi(iPhi).absSqr () / eKB (phiOrbNl(iPhi));
               }
            }
         }
      }
   } else {
      if (nNlOrb > 0 && !rsProj2)  {
         // --- get projections <phi|psi> = < mu | Vps | Psi >
         //     and             <phi|G|psi>
         //     as a ((1+3) x nOrb, n) matrix

         allPotPsiFac.reformat (4 * nNlOrb, nStates);
         if (allStatesAvailable)  {
            const SxPW &wavesRef = (const SxPW &)waves;
            PsiG pPsi = *nlProj | wavesRef(iSpin, ik);
            SxArray<SxDiracMat<TPrecCoeffG> > pPsiGrad 
               = nlProj->gradProject (wavesRef(iSpin, ik));
            for (int iState = 0; iState < nStates; ++iState)  {
               for (int iNlOrb = 0; iNlOrb < nNlOrb; ++iNlOrb)  {
                  allPotPsiFac(4 * iNlOrb, iState) = pPsi(iNlOrb, iState);
                  for (int iDir = 0; iDir < 3; iDir++)
                     allPotPsiFac(4 * iNlOrb + 1 + iDir, iState) 
                        = I * pPsiGrad(iDir)(iNlOrb, iState);
               }
            }
         } else {
            for (int iState = 0; iState < nStates; ++iState)  {
               // --- inefficient
               PsiG psi = waves(iState, iSpin, ik);
               PsiG pPsi = *nlProj | psi;
               SxArray<SxDiracMat<TPrecCoeffG> > pPsiGrad 
                  = nlProj->gradProject (psi);
               for (int iNlOrb = 0; iNlOrb < nNlOrb; ++iNlOrb)  {
                  allPotPsiFac(4 * iNlOrb, iState) = pPsi(iNlOrb);
                  for (int iDir = 0; iDir < 3; iDir++)
                     allPotPsiFac(4 * iNlOrb + 1 + iDir, iState) 
                        = I * pPsiGrad(iDir)(iNlOrb);
               }
            }
         }
         
         // --- sum over projectors to get nonlocal energy and forces
         SxQuantumNumbers orb; 
         int d;
         SxDiracVec<TPrecCoeffG> potPsiFac;
         for (int iNlOrb=0; iNlOrb < nNlOrb; iNlOrb++)  {  //:is,ia,n,l,m

            // get <mu | Vps | Psi_n> for current mu and all n
            potPsiFac = allPotPsiFac.row(4 * iNlOrb);
            
            //               |< mu | Vps | Psi_n >|^2 
            // E  =  Sum     -------------------------  * focc(n)
            //                  < mu | Vps | mu >
            nlFac = Gk.weights(ik) / eKB(phiOrbNl(iNlOrb));
            eNl += nlFac * (  potPsiFac.absSqr () 
                            * focc(iSpin, ik)     ).sum ();
         
            // --- forces (at atom that phi belongs to)
            //               <mu|Vps|iG|Psi_n><Psi_n|Vps|mu>
            // F  = -2Re Sum -------------------------------  * focc(n)
            //                      < mu | Vps | mu >
            orb = psiOrbNl(iNlOrb);
            potPsiFac *= focc(iSpin, ik);
            for (d=0; d < 3; d++)  {  // {:xyz}
               fNl.ref(orb.iSpecies,orb.iAtom)(d) -= 2.
                  * nlFac * dot(allPotPsiFac.row(4 * iNlOrb + d + 1),
                                potPsiFac                       ).im;
            }
         }
      }
      if (psPot.useRealSpace ()) {
         // --- EES-G real-space method
         SxVector<TPrecCoeffG> phiPsi;
         int is;
         // loop over reference orbitals
         for (int iPhi = 0; iPhi < phiNl(ik).getSize (); ++iPhi)  {
            is = phiNl(ik)(iPhi).handle->auxData.is;
            SX_CHECK (is >= 0 && is < structure.getNSpecies (),
                      is, structure.getNSpecies ());
            if (psPot.realSpace(is))  {
               SxAtomicStructure atoms = SxSpeciesRef(is) | structure;
               double preFac;
               for (i = 0; i < nStates; ++i)  {
                  // skip unoccupied states
                  if (fabs(focc(i, iSpin, ik)) < 1e-12) continue;
                  // --- get projections < mu | Vps | Psi >
                  phiPsi = rsProjector->projectGrad (phiNl(ik)(iPhi),
                                                     waves(i,iSpin,ik),
                                                     atoms);
                  // --- get nonlocal energy (see above)
                  preFac = Gk.weights(ik) / eKB (iPhi)
                         * focc(i, iSpin, ik);
                  eNl += preFac * phiPsi.row(0).absSqr ().sum ();
                  int ia, iDir;
                  for (ia = 0; ia < structure.getNAtoms(is); ++ia)  {
                     for (iDir = 0; iDir < 3; ++iDir)  {
                        fNl.ref(is,ia)(iDir) -= 2. * preFac
                                              * (phiPsi(iDir+1,ia).conj ()
                                                 * phiPsi(0,ia)).re;
                                       
                     }
                  }
               }
            }
         }
      }
      if (rsProj2)  {
         // forces not yet implemented
         SX_EXIT;
      }
   }

}



// ==========================================================================
//   I O N  -  I O N  contribution (realspace part of the Ewald summation)
// ==========================================================================
void SxPWHamiltonian::computeScrPotential ()
{
   const SxAtomicStructure &str = structure;

   double rGaussMax = 0.;
   for (int is = 0; is < psPot.getNSpecies (); ++is)
      if (rGaussMax < psPot.rGauss(is)) rGaussMax = psPot.rGauss(is);
   const Real8 rEwald = 20. * rGaussMax; // more than enough

   // --- for neighbor search
   SxGrid grid (str, 10);
   SxNeighbors neighbors;
   int neighborParams = SxNeighbors::StoreRel + SxNeighbors::StoreDistSqr;

   eEwald = 0.;
   if (calcForces)  {
     if (fScr.getSize () == 0) fScr = structure.getNewStr ();
     fScr.set (Coord(0., 0., 0.));
   }

   // --- loop over all atoms
   int is, ia, js;
   SxIdx jsRange;
   Real8 zvI, zvJ, rGaussI, rGaussJ, rGaussNorm;
   SxVector<Double> dist, arg, dEsr;
   SxAtomicStructure dForce;
   for (is=0; is < str.getNSpecies (); is++)  {  // --- species i
      zvI     = psPot.valenceCharge(is);
      rGaussI = psPot.rGauss(is);
      for (ia=0; ia < str.getNAtoms(is); ia++)  {  // --- atom i

         // compute neighbors within rEwald
         neighbors.compute (grid, str, str.getAtom(is,ia), 
                            rEwald, neighborParams);

         for (js=0; js < str.getNSpecies (); js++)  { // --- species j
            
            // are there neighbors of that species?
            if (neighbors.relPositions.getNAtoms(js) == 0) continue;
            
            zvJ      = psPot.valenceCharge(js);
            rGaussJ  = psPot.rGauss(js);


            jsRange = neighbors.relPositions.getRange(js); // atoms j
            dist  = sqrt(neighbors.distSqr(jsRange)); // for all atoms j
            
            // --- check for crashing atoms
            {
               int minIdx;
               if (dist.minval (&minIdx) < 0.30)  {
                  int ja = neighbors.idx(minIdx) - str.atomInfo->offset(js);
                  sxprintf ("Ewald: Atoms are closer than 0.3 Bohr!\n");
                  sxprintf ("is=%d ia=%d <--> is=%d ia=%d\n", 
                          is,ia,js,ja);
                  sxprintf ("dist = %g\n", dist(minIdx));
                  SX_EXIT;
               }
            }

            // --- energy
            rGaussNorm = 1. / sqrt(rGaussI * rGaussI  +  rGaussJ * rGaussJ);
            arg     = dist * rGaussNorm;  // for all atoms j
            dEsr    = (0.5 * zvI * zvJ) * derfc (arg) / dist; // for all atoms j
            eEwald += dEsr.sum ();

            if (calcForces)  {
               // --- forces
               SxSpeciesRef jsFilter (js, SxSpeciesRef::MapToParentOfFiltered);

               double fac = zvI * zvJ * rGaussNorm / SQRT_PI;

               // forces for atoms j
               dForce = (   dEsr + fac * exp(-arg.sqr())  )
                        / neighbors.distSqr(jsRange)
                      * (jsFilter | neighbors.relPositions);

               fScr            += dForce;        // forces for atoms j
               fScr.ref(is,ia) -= dForce.sum (); // opposite force for atom i
            }
         }
      }
   }

}





// ==========================================================================
//   U T I L I T Y    Functions
// ==========================================================================
PrecEnergy SxPWHamiltonian::getEnergy (const SxPsiSet &psiSet,
                                       const SxFermi &fermi)
{
   wavesPtr = dynamic_cast<const SxPWSet*>(&psiSet);
   SX_CHECK (wavesPtr);
   update (fermi);
   return eTotal;
}

void SxPWHamiltonian::writeRho (const SxString &filename) const
{
   rho.writeRho (filename);
}

SxDensity& SxPWHamiltonian::getRho ()
{
   return rho;
}

// TODO: to be called update ()
void SxPWHamiltonian::set (const SxPWSet &psiSet, const SxFermi &fermi)
{
   //TODO --- Bad pseudopotential tight binding initialization
   //         fermi is used for waves and lcao initialization.
   //         Therefore, fermi may have more spin channels than init-waves, but not less  
   SX_CHECK(fermi.getNStates() == psiSet.getNStates(), fermi.getNStates(), psiSet.getNStates());
   SX_CHECK(fermi.getNSpin() >= psiSet.getNSpin(), fermi.getNSpin(), psiSet.getNSpin());
   SX_CHECK(fermi.getNk() == psiSet.getNk(), fermi.getNk(), psiSet.getNk());
   
   wavesPtr = &psiSet;
   
   compute (fermi, true, true);
}

SxDiracSymMat<TPrecCoeffG>
SxPWHamiltonian::getMatrix (const SxGkBasis &gkBasis,
                            int iSpin, int ik, int nG_)
{
   const SxGBasis &G = *gBasisPtr;
   const SxGBasis *GkPtr = &( gkBasis(ik) );
   SX_FFT_REAL     fftFactorRtoG = G.fft3d(0).scaleFor; // fft factor R --> G
   int             iG, iGp, nG = (nG_ == -1) ? GkPtr->ng : nG_;
   int             idxG;

   SX_CHECK (nG <= GkPtr->ng, nG, GkPtr->ng);

   // TODO: Should be just one entity.
   SxDiracSymMat<TPrecCoeffG> H(nG), K, Veff(nG), Vnl(nG);

   // --- kinetic part
   K = K.identity (0.5 * GkPtr->g2, nG);

   // --- effective potential
   SxMeshG VeffG = ( G | vEffR (iSpin) ) * fftFactorRtoG;

   for (iG = 0; iG < nG; iG++)  {
      for (iGp = iG; iGp < nG; iGp++)  {
         idxG = G.getIdxGDiff (iG, GkPtr, iGp, GkPtr);
         Veff(iG,iGp) = (idxG == -1) ? (PrecCoeffG)0. : VeffG(idxG);
      }  // :iGp
   }  // :iG

   // --- non-local pseudopotential
   //                            _____
   //                            \      < G'+k | Vps | mu >  < mu | Vps | G+k >
   //     < G'+k | Vnl | G+k > =  >     ---------------------------------------
   //                            /                 < mu | Vps | mu >
   //                           /__mu__
   //                            _____  
   //                            \      VnlG(G'+k)   VnlG^* (G+k)
   //                          =  >    ---------------------------
   //                            /                Ekb
   //                           /__mu__
   Vnl.set(0.);

   SX_LOOP(iNlOrb)  {
      ssize_t phiIdx = phiOrbNl(iNlOrb);
      double diag   = eKB(phiIdx);       // < mu | Vps | mu >

      PsiG phiNlGk = nlProj->getAOinG (ik, (int)iNlOrb);

      for (iG = 0; iG < nG; iG++)  {
         for (iGp = iG; iGp < nG; iGp++)  {
            Vnl(iG,iGp) += phiNlGk(iG) * phiNlGk(iGp).conj() / diag;
         }  // :iGp
      }  // iG
   }  // :iNlOrb

   // TODO: Use only H!
   H = K + Veff + Vnl;
   
   H.handle->parameters |= (IS_TRIANGULAR | UPPER_RIGHT);
   H.reshape (nG,nG);

   SX_CHECK (H.handle->parameters & UPPER_RIGHT);
   SX_CHECK (H.handle->parameters & IS_TRIANGULAR);

   return H;
}

const SxSymbolTable *
SxPWHamiltonian::getHamiltonianGroup (const SxSymbolTable *table)
{
   return (table->containsGroup ("PWHamiltonian")) ? table->getGroup("PWHamiltonian") : NULL;
}


void SxPWHamiltonian::read (const SxSymbolTable *table)
{
   const SxSymbolTable *hamiltonian = NULL;
   try  {
      hamiltonian = table->getGroup("PWHamiltonian");
//    if (hamiltonian->contains("nEmptyStates"))
//       nEmptyStates = hamiltonian->get("nEmptyStates")->toInt();
//    if (hamiltonian->contains("nExcessElectrons"))
//       nExcessElectrons = hamiltonian->get("nExcessElectrons")->toReal();
//    if (hamiltonian->contains("ekt"))
//       ekt = hamiltonian->get("ekt")->toReal();

      if (hamiltonian->containsGroup("vdwCorrection")) {
         bool applyVDWCorrection = true;
         SxVDW vdwCorrection(structure, table);
      }

      ekt = hamiltonian->get("ekt")->toReal() / HA2EV;

      xcPtr->read (hamiltonian, structure.cell);
      if (xcPtr->nSpin == 0)
         xcPtr->nSpin = SxHamiltonian::getNSpin (table);

      if (hamiltonian->contains ("hContrib"))   {
         contrib = (Contrib)hamiltonian->get("hContrib")->toInt();
         if ( !(contrib & CALC_X) )  xcPtr->disableExchange ();
         if ( !(contrib & CALC_C) )  xcPtr->disableCorrelation ();
      }
      

      if (hamiltonian->containsGroup("vExt"))  {
         SxSymbolTable *vExt = hamiltonian->getGroup("vExt");
         SxString vExtFile = vExt->get("file")->toString();
         vExtActsOnNuclei = true;
         if (vExt->contains("actOnNuclei"))  {
            vExtActsOnNuclei = vExt->get("actOnNuclei")->toAttribute ();
         }
         if (vExtActsOnNuclei && xcPtr->nSpin == 2)  {
            cout << "No classical external potential for spin-polarized case!" << endl;
            cout << "Spin-polarized potentials cannot act on nuclei" << endl;
            SX_QUIT;
         }

         try  {
            SxBinIO io (vExtFile, SxBinIO::BINARY_READ_ONLY);
            SX_CHECK (rho.rBasisPtr);
            const SxRBasis &R = *rho.rBasisPtr;
            int nSpinFile = io.getDimension ("nMeshes");
            int nSpin = xcPtr->nSpin;
            if (nSpinFile != nSpin)  {
               cout << "Spin inconsistency when reading external potential '"
                    << vExtFile << "'.\n";
               cout << "There should be " << nSpin << "spin channel";
               if (nSpin > 1) cout << 's';
               cout << ", but found " << nSpinFile << " in file." << endl;
               SX_QUIT;
            }
            SxMatrix3<Double> cellFile;
            SxVector3<Int>    dim;
            vExtR.resize (nSpin);
            vExtR = io.readMesh (&cellFile, &dim);
            if ((cellFile - R.cell).absSqr ().sum () > 1e-7)  {
               cout << "cell inconsistency in external potential." << endl;
               cout << "Expected " << R.cell << endl;
               cout << "Found " << cellFile << endl;
               SX_EXIT;
            }
            if (!(dim == R.fft3d.mesh))  {
               cout << "Mesh inconsistency in external potential." << endl;
               cout << "Expected " << R.fft3d.mesh << ", found " << dim << endl;
               SX_EXIT;
            }
            for (int iSpin = 0; iSpin < nSpin; ++iSpin)  
               vExtR(iSpin).setBasis (&R);
         } catch (SxException e) {
            e.print ();
            SX_EXIT;
         }

         if ( !(contrib & CALC_EXT) )  {
            printContribMask ();
            sxprintf ("Error:   if an external potential is provided \n"
                    "         hContrib must contain CALC_V_EXT.\n");
            sxprintf ("Example: hContrib = CALC_DEFAULT + CALC_V_EXT;\n");
            SX_QUIT;
         }
      }

      if (hamiltonian->contains("nlBlockSize"))  {
         nlBlockSize = hamiltonian->get("nlBlockSize")->toInt ();
      }
      if (hamiltonian->contains("dipoleCorrection"))  {
         dipoleCorrection = hamiltonian->get("dipoleCorrection")
                            ->toAttribute ();
      }
      zField = (hamiltonian->contains ("zField"))
             ? hamiltonian->get("zField")->toReal () / HA2EV
             : 0.;

   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   
   if (psPot.useRealSpace ())  {
      SxVector3<Int> mesh, minMesh;
      int splineOrder;
      try  {
         SxSymbolTable *nlEES = hamiltonian->getGroup ("nlEES");
         double eCut = SxGBasis::getECut (table);
         minMesh = SxGBasis::getMeshSize (eCut, structure.cell, 0.5);
         if (nlEES->contains ("mesh"))  {
            mesh = SxVector3<Int>(nlEES->get("mesh")->toIntList ());
            if (!SxGBasis::isCommensurableMesh(mesh, structure.cell, true))  {
               cout << "ERROR: The provided nlEES mesh is not commensurable!\n";
               SxGBasis::getCommensurableMesh (structure.cell, mesh);
               SX_QUIT;
            }
         } else {
            double meshAccuracy = nlEES->get("meshAccuracy")->toReal ();
            mesh = SxGBasis::getCommensurableMesh (eCut, structure.cell,
                                                   meshAccuracy);
         }
         splineOrder = nlEES->get("splineOrder")->toInt ();
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
      for (int dir = 0; dir < 3; ++dir)  {
         if (mesh(dir) < minMesh(dir))  {
            cout << "ERROR: the nlEES mesh is too small." << endl;
            cout << "Minimum size: " << minMesh << endl;
            cout << "Actual size:  " << mesh << endl;
            SX_QUIT;
         }
         if (splineOrder >= mesh(dir))  {
            cout << "ERROR: the nlEES splineOrder (" << splineOrder
                 << ") must be smaller than the mesh (" << mesh << ").\n";
            SX_QUIT;
         }
      }
      rsProjector = SxPtr<SxEESGProj>::create (structure.cell, mesh,
                                               splineOrder);
   }

   // --- new real-space projectors
   if (hamiltonian->containsGroup ("rsProj"))  {
      int nSpecies = structure.getNSpecies ();
      // --- setup projectors
      SxRadBasis rad (psPot.rad, psPot.logDr);
      SxArray<SxArray<SxDiracVec<Double> > > proj(nSpecies);
      for (int is = 0; is < nSpecies; ++is)  {
         int npt = (int)psPot.pseudoPsi(is).getSize ();
         proj(is).resize (npt-1);
         for (int ipt = 0, ip = 0; ipt < npt; ipt++)  {
            if (psPot.lLoc(is) == ipt) continue;
            proj(is)(ip) = psPot.pseudoPsi(is)(ipt) 
                         * psPot.pseudoPot(is)(ipt);
            proj(is)(ip).handle->auxData.is = is;
            proj(is)(ip).handle->auxData.l = ipt;
            proj(is)(ip).setBasis (&rad);
            ip++;
         }
      }

      // create real-space projectors
      rsProj2 = SxPtr<SxRSProj>::create (hamiltonian->getGroup ("rsProj"),
                                         structure,
                                         proj);
   }

}

bool SxPWHamiltonian::rereadTable (const SxSymbolTable *table)
{
   SX_CHECK (table);
   SX_CHECK (xcPtr);
   bool changed;
   // xc mesh
   changed = xcPtr->rereadTable (table);

   // --- dipole correction
   if (table->contains ("dipoleCorrection"))  {
      changed = true;
      bool dipole = table->get ("dipoleCorrection")->toAttribute ();
      cout << "| Dipole correction is now switched "
           << (dipole ? "on" : "off") << "." << endl;

      // is dipole correction being switched?
      if (dipole != dipoleCorrection)  {
         cout << "| Recomputing potentials ..." << endl;
         SX_CHECK (rho.rBasisPtr);
         const SxRBasis &R    = *rho.rBasisPtr;
         SxMeshR rhoRTotal (R);
         rhoRTotal.set(0.);
         for (int iSpin=0; iSpin < rho.rhoR.getSize (); iSpin++)
            rhoRTotal += rho(iSpin);
         SxMeshR rhoHartree = rhoRTotal - rhoGaussG.to (R);
         SX_CHECK (R.cell.volume > 0., R.cell.volume);
         double charge = rhoHartree.sum () * R.cell.volume
                         / (double)R.getNElements ();
         rhoHartree -= charge / R.cell.volume;
         if (fabs(charge) < 1e-5) charge = 0.;
         SxDipoleCorrZ dipCorr (rhoRTotal, rhoHartree, charge);
         SxMeshR vDip(R);
         vDip.set (0.);
         dipCorr.correctPotential (&vDip);
         for (int iSpin = 0; iSpin < vEffR.getSize (); ++iSpin)  {
            if (dipole)
               vEffR(iSpin) += vDip;
            else
               vEffR(iSpin) -= vDip;
         }
      }

      dipoleCorrection = dipole;
   }

   return changed;
}

void SxPWHamiltonian::backToDefault (const SxSymbolTable *table)
{
   SX_CHECK (xcPtr);
   // xc basis
   xcPtr->rereadTable (NULL);

   // --- dipole correction
   bool dipole = withDipoleCorrection (table->topLevel ());
   if (dipole != dipoleCorrection)  {
      cout << "| Dipole correction is now switched "
           << (dipole ? "on" : "off") << "." << endl;
      cout << "| Recomputing potentials ..." << endl;
      SX_CHECK (rho.rBasisPtr);
      const SxRBasis &R    = *rho.rBasisPtr;
      SxMeshR rhoRTotal (R);
      rhoRTotal.set(0.);
      for (int iSpin=0; iSpin < rho.rhoR.getSize (); iSpin++)
         rhoRTotal += rho(iSpin);
      SxMeshR rhoHartree = rhoRTotal - rhoGaussG.to (R);
      SX_CHECK (R.cell.volume > 0., R.cell.volume);
      double charge = rhoHartree.sum () * R.cell.volume
                      / (double)R.getNElements ();
      rhoHartree -= charge / R.cell.volume;
      if (fabs(charge) < 1e-5) charge = 0.;
      SxDipoleCorrZ dipCorr (rhoRTotal, rhoHartree, charge);
      SxMeshR vDip(R);
      vDip.set (0.);
      dipCorr.correctPotential (&vDip);
      for (int iSpin = 0; iSpin < vEffR.getSize (); ++iSpin)  {
         if (dipole)
            vEffR(iSpin) += vDip;
         else
            vEffR(iSpin) -= vDip;
      }
   }
   dipoleCorrection = dipole;
}


void SxPWHamiltonian::validateNStates (int nStates)
{
   int nStatesMax = wavesPtr->getGkBasis().nGkMin;

   if (nStates > nStatesMax)  {
      sxprintf ("Error: input parameter(s) \"nEmptyStates\" or "
              "\"nEmptyStatesChi\" in (PW)Hamil-\n"
              "       tonian group chosen too large. Please reduce it/them or "
              "try with a\n"
              "       larger cutoff \"eCut\". The number of states you want to "
              "compute must\n"
              "       not exceed %d.\n", nStatesMax);
      SX_QUIT;
   }
}

void SxPWHamiltonian::printEnergies () const
{
   double UNIT = printHartree ? 1.  :  HA2EV;
   SxString u  = printHartree ? "H" : "eV";

   cout << endl;
   cout << SX_SEPARATOR;
   sxprintf ("| eKin                = %14.8f %s\n", UNIT * eKin, u.ascii());
   sxprintf ("| Self energy         = %14.8f %s\n", UNIT * eSelf, u.ascii());
   sxprintf ("| e Hartree           = %14.8f %s\n", UNIT * eHartree, u.ascii());
   sxprintf ("| Electrostatic       = %14.8f %s\n", UNIT *(eHartree- eSelf), 
                                                  u.ascii());
   sxprintf ("| e local pseudo      = %14.8f %s\n", UNIT * eLocPs, u.ascii());
   sxprintf ("| e non-local         = %14.8f %s\n", UNIT * eNl, u.ascii());
   sxprintf ("| xc energy           = %14.8f %s\n", UNIT * xcPtr->eXc, u.ascii());
   //sxprintf ("| xc potential energy = %14.7f %s\n", UNIT * xcPtr->eVxc,
   //                                               u.ascii());
   sxprintf ("| scr. energy         = %14.8f %s\n", UNIT * eEwald, u.ascii());
   if (contrib & CALC_EXT)
      sxprintf ("| ext. energy         = %14.8f %s\n", UNIT * eExt, u.ascii());

   sxprintf ("| TOTAL ENERGY        = %14.8f %s\n", UNIT * eTotal, u.ascii());

   cout << SX_SEPARATOR;
//   fermiPtr->printOccupation ();
// if (calcForces)  {
//    printForces ();
//    cout << SX_SEPARATOR;
// }
   cout.flush ();
}

void SxPWHamiltonian::writeVElStat (const SxString &file) const
{
   if (!((contrib & CALC_HARTREE) && (contrib & CALC_LOC))) return;
   try  {

      const SxRBasis *R = dynamic_cast<const SxRBasis *>(rho.rBasisPtr);
      SX_CHECK(R);
      
      SxBinIO io (file, SxBinIO::BINARY_WRITE_ONLY);

      const SxVector3<Int> mesh = R->getMesh ();
      io.writeMesh ((vLocR + vHartreeR)*HA2EV, R->cell, mesh);  
      io.setMode (SxBinIO::WRITE_DATA);
      io.writeMesh ((vLocR + vHartreeR)*HA2EV, R->cell, mesh);  

      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

}

void SxPWHamiltonian::printContribMask () const
{
   cout << endl;
   cout << SX_SEPARATOR;
   cout << "| Hamiltonian consists of the following contributions:\n";
   if (contrib & CALC_KIN)      cout << "|    kinetic part" << endl;
   if (contrib & CALC_HARTREE)  {
                                cout << "|    Hartree part" << endl;
      if (dipoleCorrection)     cout << "|    dipole correction" << endl;
   }
   if (contrib & CALC_LOC)      cout << "|    LocPs part" << endl;
   if (contrib & CALC_X)        cout << "|    exchange part" << endl;
   if (contrib & CALC_C)        cout << "|    correlation part" << endl;
   if (contrib & CALC_NL)       cout << "|    non-loc part" << endl;
   if (contrib & CALC_EXT) {
                                cout << "|    " 
                                     << (vExtActsOnNuclei ? "classical"
                                                          : "xc-like")
                                     << " external potential" << endl;
   }
   if (contrib & CALC_RHO)      cout << "|    rho will be updated." << endl;
   cout << SX_SEPARATOR;
   cout.flush ();

}

SxMatrix<TPrecCoeffG> 
SxPWHamiltonian::getNlMatrix (const SxMuPW &mu, int ik,
                              const SxAtomicStructure &str) const
{
   SX_CHECK (rsProjector);
   SX_CHECK (ik >= 0 && ik < mu.getNk (), ik, mu.getNk ());
   int nOrb = mu.getNStates ();
   SxMatrix<TPrecCoeffG> res(nOrb, nOrb);
   const SxArray<PsiG> &muRef = mu.muPhi(ik);
   res.set (0.);

   int nMu = (int)muRef.getSize ();
   SxArray<SxVector<TPrecCoeffG> > pPhiMu(nMu);
   SxArray<SxConstPtr<SxAtomInfo> > info(nMu);
   SX_CHECK (wavesPtr);
   Coord kVec = wavesPtr->getGkBasis ().getK (ik); // TODO get k from G basis
   SxAtomicStructure atPhi, atMu;
   int ia, ka;
   for (int iPhi = 0; iPhi < phiNl(ik).getSize (); ++iPhi)  {
      atPhi = SxSpeciesRef(phiNl(ik)(iPhi).handle->auxData.is) | str;

      // --- set up <phi @ ka|mu @ ia>
      for (int iMu = 0; iMu < nMu; ++iMu)  {
         atMu = SxSpeciesRef(muRef(iMu).handle->auxData.is) | str;
         
         // --- set up relative distances
         SxAtomicStructure distPhiMu;
         for (ia = 0; ia < atMu.getNAtoms (); ++ia)  {
            distPhiMu.newSpecies ();
            for (ka = 0; ka < atPhi.getNAtoms (); ++ka)  {
               distPhiMu.addAtom (atPhi(ka) - atMu(ia));
            }
         }
         distPhiMu.endCreation ();
         SX_CHECK (distPhiMu.getNSpecies () == atMu.getNAtoms (),
                   distPhiMu.getNSpecies (), atMu.getNAtoms ());
         // store info
         info(iMu) = distPhiMu.atomInfo;
         
         // get projection
         pPhiMu(iMu) = rsProjector->project (phiNl(ik)(iPhi), muRef(iMu),
                                             distPhiMu);

         // --- add mu phase exp(-ikr_mu) for consistency with mu(i,iSpin,ik)
         {
            SxVector<TPrecCoeffG>::Iterator pIt = pPhiMu(iMu).begin ();
            SxComplex16 phase;
            for (ia = 0; ia < atMu.getNAtoms (); ++ia)  {
               phase = SxEESGProj::getPhase (-(kVec ^ atMu.constRef (ia)));
               for (ka = 0; ka < atPhi.getNAtoms (); ++ka) *pIt++ *= phase;
            }
         }
         // phi phase exp(ikr_phi) cancels out
      }

      // --- add phi contribution to Hamiltonian
      for (int iOrb = 0; iOrb < nOrb; ++iOrb)  {
         int iMu = mu.getRefIdx (iOrb);
         SX_CHECK (iMu != -1);
         int iAtom = mu.psiOrb(iOrb).iAtom;
         for (int jOrb = 0; jOrb < nOrb; ++jOrb)  {
            int jMu = mu.getRefIdx(jOrb);
            SX_CHECK (jMu != -1);
            int jAtom = mu.psiOrb(jOrb).iAtom;
            for (ka = 0; ka < atPhi.getNAtoms (); ++ka)  {
               int pIdxI = info(iMu)->offset(iAtom) + ka;
               int pIdxJ = info(jMu)->offset(jAtom) + ka;

               // sum_ka <mu @ ia|phi @ ka> E^{-1}(phi) <phi @ ka| mu @ ja>
               res(iOrb, jOrb) += pPhiMu(iMu)(pIdxI).conj ()
                                * pPhiMu(jMu)(pIdxJ) / eKB(iPhi);
            }
         }
      }
   }
   VALIDATE_VECTOR(res);
   return res;
}



