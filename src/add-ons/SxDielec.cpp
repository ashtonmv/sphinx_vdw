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

#include <SxConfig.h>
#include <SxDielec.h>
#include <SxConstants.h>
#include <SxMatrix3.h>
#include <SxRBasis.h>
#include <SxRho.h>

/*
Coord kgrid;
int nff;
*/
#ifndef SX_STANDALONE
#include <SxTimer.h>
#include <SxProjector.h>

SxMatrix3<TPrecCoeffG>
SxDielec::compute (const SxPW    &waves, 
                   const SxFermi &fermi,
                   const SxPerturbK &kp,
                   double imagFreq,
                   const SxCell &cell,
                   bool   computeWing)
{
   SxArray<double> frequencies(1);
   frequencies(0) = imagFreq;
   compute (waves,fermi,kp,frequencies,cell,computeWing);
   return head(0);
}

/*
double ffIntegrate (double e0, double beta, Coord ep)
{
   if (beta == 0.) return 0.;
   double res = 0.;
   double e, f;
   ep /= kgrid;
   const int n = nff;
   for (int x = -n; x <= n; ++x)  {
      for (int y = -n; y <= n; ++y)  {
         for (int z = -n; z <= n; ++z)  {
            e = e0 + (ep(0) * x + ep(1) * y + ep(2) * z)/double(2 * n + 1);
            f = 1. / ( 1. + exp (beta * e));
            res += f * (1. - f);
         }
      }
   }
   return res / (double(2 * n + 1) * double(2 * n + 1) * double(2 * n + 1));
}
*/

void
SxDielec::compute (const SxPW    &waves, 
                   const SxFermi &fermi,
                   const SxPerturbK &kp,
                   SxArray<double> imagFreq,
                   const SxCell &cell,
                   bool   computeWing)
{
   enum {TOTALTIME,CPVTIME, FFTG2R, SUM_COND, CV_PROD_R, SYMWING, nTimers};
   SxTimer timer (nTimers);
   timer.setName (TOTALTIME, "total time");
   timer.setName (CPVTIME, "<c|p|v>");
   timer.setName (FFTG2R, "G->r FFT");
   timer.setName (SUM_COND, "[c(r) g(c,v)]^*");
   timer.setName (CV_PROD_R, "[c(r)g(c,v)]^* v(r)");
   timer.setName (SYMWING, "symmetrize wings");

   timer.start (TOTALTIME);
   int nk = waves.getNk ();
   int iSpin = 0;
   int nOmega = int(imagFreq.getSize ());

   int nSym = cell.symGroupPtr->getNSymmorphic ();
   SX_CHECK (nSym > 0, nSym);

   
   // --- determine number of valence and conduction states
   int nStates    = waves.getNStates (),
       nVal       = fermi.getNValenceBands (), 
       nCond      = fermi.getNConductionBands (),
       iCondStart = nStates - nCond;
   
   cout << "nVal = " << nVal << "; nCond = " << nCond << endl;
   
   if (nCond < 3 * nVal)  {
      cout << SX_SEPARATOR;
      cout << "| WARNING: the number of conduction bands is very low!" << endl;
      cout << SX_SEPARATOR;
   }
   SX_CHECK (nVal > 0, nVal);

   // matrix elements (nCond x nVal x 3) of perturbation operator k_i*p
   // SxDiracMat<TPrecCoeffG> kpElements(nVal*nCond, 3);
   SxDiracMat<TPrecCoeffG> kpElements;
   SxDiracMat<TPrecCoeffG> aux(nVal * nCond, 3);
   // frequency factors
   SxDiracMat<Double> freqHead(nCond, nVal);
   SxDiracMat<Double> freqWing;
   SxDiracMat<Double>::Iterator freqHeadIt, freqWingIt;

   // bands in R-space
   SxDiracMat<TPrecCoeffG> valBandsR;
   PsiG wCond;
   // contribution to wing from current k-point
   SxDiracMat<TPrecCoeffG> contrib;
   int meshSize = 0, iw, iDir;

   if (computeWing)  {
      const SxRBasis &R = dynamic_cast<const SxGBasis *> 
                          (waves(0,0).getBasisPtr ())->getRBasis ();
      meshSize = R.fft3d.meshSize;
      freqWing.reformat (nCond, nVal);

      wings.resize (nOmega);
      for (iw = 0; iw < nOmega; iw++)  {
         wings(iw).reformat(meshSize,3);
         wings(iw).set(0.);
         wings(iw).setBasis (&R);
      }

      contrib.reformat(meshSize,3);
      contrib.setBasis (&R);

      valBandsR.reformat (meshSize, nVal);
   }

   SxDiracMat<Double> result(3,3);
   head.resize (nOmega);
   for (iw = 0; iw < nOmega; iw++) head(iw).set(0.);

   int r,c;
   double epsV, foccV, kWeight, wSquared;
   double dEps, dFocc, dEpsW2;
   SxDiracVec<TPrecEps>::Iterator epsIt, foccIt;
   const SxGBasis *gBasis = NULL;
   bool zeroFreq;

   // --- loop over k-points
   for (int ik = 0; ik < nk; ik++)  {

      (cout << "kp-Mat: ik = " << (ik+1) << endl).flush ();

      // --- get <psi_v|p|psi_c>
      {
         timer.start (CPVTIME);
         PsiG wavesRef, wVal; // WARNING: dangerous SxIdx-type references
         wavesRef = waves(iSpin, ik);
         int ng = int(wavesRef.nRows ());
         wavesRef.reshape (ng * nStates);
         wVal = wavesRef(SxIdx(0, ng * nVal - 1));
         wVal.reshape (ng, nVal);
         wCond = PsiG (); // free reference
         wCond = wavesRef(SxIdx(ng * iCondStart, ng * nStates - 1));
         wCond.reshape (ng, nCond);
         wavesRef.reshape (ng, nStates);
         kpElements = kp.getMatrixElements (wCond, wVal);
         timer.stop (CPVTIME);
         // get G+k basis (needed for FFTs below)
         gBasis = dynamic_cast<const SxGBasis *>(wavesRef.getBasisPtr ());
         SX_CHECK (gBasis);
      }

      if (computeWing)  {
         // --- FT bands
         const SxRBasis &R = dynamic_cast<const SxGBasis *> 
                             (waves(iSpin,ik).getBasisPtr ())->getRBasis ();
         (cout << "Fourier transform bands..." << endl).flush ();
         timer.start(FFTG2R);
         for (int iv = 0; iv < nVal; iv++)
            valBandsR.colRef(iv) <<= (R | waves(iv,iSpin,ik)).conj ();
         timer.stop(FFTG2R);
      }
         
      for (iw = 0; iw < nOmega; iw++)  {
         result.set (0.);

         // --- get frequency factor on imaginary frequency axis

         wSquared = imagFreq(iw) * imagFreq(iw);
         zeroFreq = fabs(imagFreq(iw)) < 1e-8;

         freqHeadIt = freqHead.begin ();
         if (computeWing)
            freqWingIt = freqWing.begin ();
         for (int iv = 0; iv < nVal; iv++)  {
            epsV  = fermi.eps(iv,iSpin,ik);
            foccV = fermi.focc(iv,iSpin,ik);
            epsIt  = fermi.eps(iSpin,ik).begin ();  epsIt  += iCondStart;
            foccIt = fermi.focc(iSpin,ik).begin (); foccIt += iCondStart;
            for (int ic = 0; ic < nCond; ++ic, ++epsIt, ++foccIt) {
               if (ic + iCondStart < iv)  {
                  // avoid double counting
                  *freqHeadIt++ = 0.;
                  if (computeWing)
                     *freqWingIt++ = 0.;
                  continue;
               }
               if (ic + iCondStart == iv && !zeroFreq)  {
                  // intraband terms
                  // (0.5)^2 from 1/nSpin for the focc
                  // 0.5 to compensate double counting factor for interband
                  //     terms included in the prefactor below

                  *freqHeadIt++ = beta * 0.125 * foccV * (2.-foccV) / wSquared;
                  /*
                  SX_CHECK (fabs(kpElements.colRef(0)(ic + iv*nCond).im) < 1e-10);
                  SX_CHECK (fabs(kpElements.colRef(1)(ic + iv*nCond).im) < 1e-10);
                  SX_CHECK (fabs(kpElements.colRef(2)(ic + iv*nCond).im) < 1e-10);
                  *freqHeadIt++ = 0.5 * beta 
                                * ffIntegrate (epsV - fermi.eFermi, beta,
  cell.getReciprocalCell ().transpose () 
  ^ Coord (kpElements.colRef(0)(ic + iv*nCond),
           kpElements.colRef(1)(ic + iv*nCond),
           kpElements.colRef(2)(ic + iv*nCond)))
                                  / wSquared;
                  */
                  if (computeWing)
                     *freqWingIt++ = 0.;
                  continue;
               }
               // eps_c - eps_v
               dEps = *epsIt - epsV;
               if (fabs(dEps) < 1e-8 /* && zeroFreq */)  {
                  // this would be true metallic screening for omega=0
                  // for others, result is 0 from foccs anyway
                  // avoid divergence
                  *freqHeadIt++ = 0.;
                  if (computeWing)
                     *freqWingIt++ = 0.;
                  continue;
               }
               dFocc = -0.5 * (*foccIt - foccV); // 0.5 = 1/nSpin
               
               dEpsW2 = dEps * dEps + wSquared;
               
               // head frequency factor
               *freqHeadIt++ = dFocc/(dEps * dEpsW2);
               
               // wing frequency factor
               if (computeWing)
                  *freqWingIt++ = dFocc / dEpsW2;
            }
         }

         // --- calculate head

         // multiply with frequency factor for each direction
         aux.copy(kpElements);
         for (iDir = 0; iDir < 3; iDir++)
            aux.colRef(iDir) *= freqHead;
      
         // sum over all valence/conduction band pairs
         kWeight = waves.getGkBasis ().weights(ik);
         result += kWeight * (aux.adjoint () ^ kpElements).real ();

         // transcribe to SxMatrix3
         // TODO: implement SxDiracMat::toMatrix3
         for (r = 0; r < 3; r++)
            for (c = 0; c < 3; c++)
               head(iw)(r,c) += result(r,c);

         // --- calculate wing
      
         if (computeWing)  {
            const SxRBasis &R = dynamic_cast<const SxGBasis *> 
                                (waves(iSpin,ik).getBasisPtr ())->getRBasis ();
            // multiply with frequency factor for each direction
            aux.copy(kpElements);
            for (iDir = 0; iDir < 3; iDir++)
               aux.colRef(iDir) *= freqWing;

            aux.reshape (nCond, nVal * 3);
         
            SxDiracVec<TPrecCoeffG> gCond, wingCondSum, wingCondSumR;
            gCond.reformat(nCond, 3);

            contrib.set(0.);

            cout << "sum over valence bands";
            for (int iv = 0; iv < nVal; iv++)  {
               (cout << '.').flush ();

               // pick <all c|p(3)|current v>
               gCond.colRef(0) <<= aux.colRef(iv);
               gCond.colRef(1) <<= aux.colRef(iv + nVal);
               gCond.colRef(2) <<= aux.colRef(iv + 2*nVal);

               // sum over conduction bands
               timer.start (SUM_COND);
               
               // [sum c] <G|c><c|p|v> * f(eps_c - eps_v)
               wingCondSum = (wCond ^ gCond);
               wingCondSum.setBasis(gBasis);
               
               timer.stop (SUM_COND);
               for (iDir = 0; iDir < 3; iDir++)  {

                  timer.start (FFTG2R);
                  // [sum G] <r|G><G|c><c|p|v> * f(eps_c - eps_v)
                  //        =<r|c><c|p|v> f(eps_c - eps_v) 
                  wingCondSumR = wingCondSum.colRef(iDir).to(R);
                  timer.stop (FFTG2R);

                  timer.start (CV_PROD_R);
                  // [sum v] <r|c><c|p|v><v|r> * f(eps_c - eps_v)
                  contrib.colRef(iDir) += wingCondSumR * valBandsR.colRef(iv);
                  timer.stop (CV_PROD_R);
               }
            
            }
            (cout << endl).flush ();
            // multiply with weight and add to result
            // make use of k / -k symmetry:
            //    <r|nk> = <r|n(-k)>^*
            // <vk|p|ck> = -<v(-k)|p|c(-k)>^*
            wings(iw).plus_assign_ax (kWeight, contrib.imag ());
         
            // restore original aux shape
            aux.reshape (nCond * nVal,3);

         } // if computeWing

      } // iw
   } // ik

   // --- symmetrize and "add" prefactors
   int i, iSym;

   SxArray<SxArray<ssize_t> > idxRotated;
   if (computeWing)  {
      // --- set up S(ir) -> irRot 
      SX_CHECK (dynamic_cast<const SxGBasis *>(waves(0,0).getBasisPtr ()));
      SxFFT3d &fft = dynamic_cast<const SxGBasis *>(waves(0,0).getBasisPtr ())
                     ->getRBasis ().fft3d;
      SxArray<SxMatrix3<Int> > symOpRel;
      SxVector3<Int> rVec;
      symOpRel.resize (nSym);
      idxRotated.resize (nSym);
      SxCell meshCell( cell.col(0) / fft.mesh(0),
                       cell.col(1) / fft.mesh(1),
                       cell.col(2) / fft.mesh(2));
      for (iSym = 0; iSym < nSym; iSym++)  {
         symOpRel(iSym) = meshCell.carToRel(cell.symGroupPtr
                                            ->getSymmorphic(iSym));
         idxRotated(iSym).resize (meshSize);
      }
      for (int ir = 0; ir < meshSize; ir++)  {
         rVec = fft.mesh.getMeshVec(ir, SxMesh3D::Origin);
         for (iSym = 0; iSym < nSym; iSym++)  {
            // rotate r
            idxRotated(iSym)(ir) = fft.mesh.getMeshIdx(symOpRel(iSym) ^ rVec,
                                                       SxMesh3D::Unknown);
         }
      }

   }

   SxMatrix3<Double> symResult;
   double prefactor     = 16. * PI / cell.volume;
   double prefactorWing = - 16. * PI / sqrt(cell.volume);

   for (iw = 0; iw < nOmega; iw++)  {

      // --- symmetrize head
      symResult.set (0.);
      for (iSym = 0; iSym < nSym; iSym++) {
         const SymMat &S = cell.symGroupPtr->getSymmorphic(iSym);
         symResult += S ^ head(iw) ^ S.transpose ();
      }

      // multiply with prefactor / nSym
      symResult *= prefactor / double(nSym); 
   
      // add diagonal 1
      for (i = 0; i < 3; i++) symResult(i,i) += 1.;

      head(iw) = symResult;

      // --- symmetrize wings
      if (computeWing)  {
      
         timer.start(SYMWING);
         int ir;
         contrib.set (0.);
         SxDiracMat<Double> wingRotated, wingSym;
         wingSym.reformat (meshSize, 3);
         wingSym.set (0.);
         for (iSym = 0; iSym < nSym; iSym++)  {
            // rotate <w(q)>
            SxDiracMat<Double> STrans (cell.symGroupPtr->getSymmorphic(iSym)
                                       .transpose ());
            wingRotated = wings(iw) ^ STrans;
            for (i = 0; i < 3; i++)
               for (ir = 0; ir < meshSize; ir++)
                  wingSym(idxRotated(iSym)(ir),i) += wingRotated(ir,i);
         } // iSym
         wings(iw) = wingSym * (prefactorWing / double(nSym));
         timer.stop(SYMWING);

      } // if computeWing
      
   } // iw
   timer.stop(TOTALTIME);
   timer.print (TOTALTIME);
   
}

#else  /* SX_STANDALONE */

#include <SxCLI.h>
#include <SxBinIO.h>
#include <SxProjector.h>

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();
   // --- init S/PHI/nX timers

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.authors = "C. Freysoldt";
   SxCLI::CliArg *opt;
   cli.preUsageMessage =
      "This add-on calculates the dielectric tensor via k*p perturbation "
      "theory.";
   SxString wavesFile = 
      cli.option ("-w|--waves","file","waves file with many unoccupied states")
      .toString ("waves.sxb");

   cli.newGroup ("nonlocal pseudopotential");
   bool nonlocalContributions =
      cli.option("--nl|--nonlocal","include non-local contributions").toBool ();
   SxString inFile =
      cli.option ("-i|--input","file","S/PHI/nX input file (for nl potentials)")
      .toString ("input.sx");
   int blockSize
      = cli.option ("--blocksize","set blocksize. Smaller blocksizes may reduce "
                    "the memory consumption, larger blocksizes may speed up the "
                    " <phi|psi> part").toInt (64);

   cli.setGroup (SxCLI::generalGroup);
   
   int nStates = cli.option("-n|--nStates", "number", 
                            "number of states to use. It must not exceed "
                            "the number of states in the waves file.")
                 .toInt (0);
   cli.last ().defaultValue = "default: use all states";

   SxString epsFile 
      = cli.option ("-e|--eps", "eps file", 
                    "use energies in <eps file> rather than those waves file")
        .toString ("");

   int fermiGroup = cli.newGroup ("Fermi distribution");
   
   double ekt = cli.option("--ekt","energy [eV]",
                           "recalculate Fermi distribution and include "
                           "intraband terms (for metals at omega>0)")
                .toDouble () / HA2EV;

   /*
   kgrid = SxVector3<Int>(cli.option("--kgrid", "vect", "k grid").toIntList3 (",x"));

   nff = cli.option ("--nff", "number", "refinement mesh").toInt (5);
   */
   
   bool nointra = cli.option ("--nointra", "no intraband").toBool ();
                           
   cli.setGroup (SxCLI::generalGroup);

   bool wings = cli.option ("--wings", "calculate wings, too").toBool ();
   bool smallMesh 
      = cli.option ("--half", "use smaller FFT mesh for wings").toBool ();

   int sfMode = cli.newGroup ("single frequency");
   int nOmega = 1;
   double frequency
      = cli.option ("-f|--freq", "energy", 
                    "imaginary frequency/energy (Hartree)").toDouble (0.,0.);

   SxString wingFile
      = cli.option ("-s|--savepolar", "file",
                    "save wings of polarizability as mesh file").toString ("");
   if (!cli.error && !wings && wingFile.getSize () > 0)  {
      cout << "Wings can't be saved if they are not calculated." << endl;
      cout << "Specify --wings" << endl;
      cli.setError ();
   }

   cli.newGroup ("multiple frequency");
   cli.excludeGroup(sfMode);
   bool readFreq 
      = cli.option ("--readfreq", "read multiple frequencies from input")
        .required ().toBool ();

   opt = &cli.option("-o|--output", "file", "output file name");
   SxString outputFile = opt->toString ("");
   opt->defaultValue = "(default: screen/stdout)";

   cli.setGroup (SxCLI::generalGroup);

   opt = &cli.option ("-p|--precision","int",
                      "number of decimals to print on screen");
   if (opt->exists ())
      cout.precision (opt->toInt (false, 0));

   cli.version ("1.4");
   cli.finalize ();

   // --- read input

   SxGBasis G;
   SxRBasis R;
   SxPW waves(wavesFile, SxPW::ReadOnDemand);
   SxPtr<SxGkBasis> gkBasisPtr; 
   SxFermi fermi;
   SxAtomicStructure structure;
   SxPerturbK kp;
   double ecut = 0.;

   if (nStates > 0)
      waves.changeNStates (nStates);
   try {
      SxBinIO io (wavesFile, SxBinIO::BINARY_READ_ONLY);
      if (wings)  {
         gkBasisPtr = SxPtr<SxGkBasis>::create (io, wings);
         waves.setGkBasisPtr (gkBasisPtr);
      } else  {
         gkBasisPtr = SxPtr<SxGkBasis>::create (io);
         waves.setGkBasisPtr (gkBasisPtr);
      }
      if (nStates > 0)  {
         int nSpin = io.getDimension("nSpin");
         fermi = SxFermi (0., nStates, nSpin, *gkBasisPtr);
         fermi.read (io, true);
      } else {
         fermi.read (io);
         fermi.kpPtr = &*gkBasisPtr;
      }
      
      // Read cell including symmetries
      structure.read (io);
      
      if (wings)  {
         SxVector3<Int> mesh;
         io.read ("meshDim", &mesh);
         io.read("eCut", &ecut);

         R.set (mesh,structure.cell);
         G.set (mesh, structure.cell, ecut);
         G.registerRBasis (R);
         R.registerGBasis (G);
         SxGkBasis &gkBasis = *gkBasisPtr;
         for (int ik = 0; ik < gkBasis.getNk (); ++ik) 
            gkBasis(ik).registerRBasis (R);
      }

   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   
   SxGkBasis &gkBasis = *gkBasisPtr;

   if (waves.getNSpin () == 2)  {
      cout << "No spin polarized version yet." << endl;
      SX_EXIT;
   }
   
   if (epsFile.getSize () > 0)  {
      int i, kpInFile, statesInFile, ik;
      int nSpin = 1, iSpin = 0;
      SxFermi::peekSpectrumFile (epsFile, &kpInFile, &statesInFile);
      SxFermi epsFromFile (0., statesInFile, nSpin, gkBasis);
      if (kpInFile != gkBasis.nk)  {
         cout << "'" << epsFile << "' has " << kpInFile;
         cout << " kpoints, but " << gkBasis.nk << " in " << wavesFile << endl;
         SX_QUIT;
      }
      epsFromFile.readSpectrum (epsFile, &structure.cell, &gkBasis);
      for (ik = 0; ik < gkBasis.nk; ik++)  {
         if (statesInFile < fermi.getNStates ())  {
            i = statesInFile - 1;
            double shift = epsFromFile.eps(i,iSpin,ik) 
                           - fermi.eps(i, iSpin, ik);
            cout << "ik = " << (ik + 1); 
            cout << ": shifting states " << (statesInFile + 1) << " to ";
            cout << fermi.getNStates () << " by " << (shift * HA2EV) << " eV\n";
            // shift higher DFT states like highest QP state known
            for (i = statesInFile; i < fermi.getNStates (); i++)
               fermi.eps(i, iSpin, ik) += shift;
         }
         // copy eps
         for (i = 0; i < statesInFile && i < fermi.getNStates (); i++)
            fermi.eps(i, iSpin, ik) = epsFromFile.eps (i,iSpin,ik);
      }
   }

   SxDielec dielec;

   // --- Recalculate Fermi distribution
   if (cli.groupAvailable (fermiGroup) )  {
      // get current number of electrons
      double nEl = 0.;
      for (int ik = 0; ik < fermi.getNk (); ik++)
         nEl += fermi.focc(0/* iSpin */,ik).sum () * gkBasis.weights(ik);
      cout << "nElectrons = " << nEl << endl;
      fermi.nElectrons = nEl;
      // set kpPtr
      fermi.kpPtr = &gkBasis;
      // recalculate
      fermi.fermiDistribution (ekt);
      if (ekt > 1e-8 && !nointra) dielec.beta = 1./ekt;
   }
   fermi.updateValCon ();
   
   if (nonlocalContributions)  {
      SxParser parser;
      SxParser::Table table = parser.read (inFile);
      SxPseudoPot psPot(&*table);
      gkBasis.changeTau(structure);
      kp.set (psPot, gkBasis, structure);
      kp.blockSize = blockSize;
   }

   // change to smaller FFT mesh?
   if (wings && smallMesh)  {
      // FFT mesh dimensions with ecut / 2
      SxVector3<Int> mesh 
         = SxGBasis::getCommensurableMesh (0.5 * ecut, structure.cell);

      SxFFT3d newFFT(SxFFT3d::Both, mesh, structure.cell.volume);
      
      for (int ik = 0; ik < gkBasis.getNk (); ik++)
         gkBasis(ik).replaceMesh (newFFT);
      G.replaceMesh(newFFT);
      R.set (mesh, structure.cell);
   }


   // --- output

   SxArray<double> frequencies;
   if (!readFreq)  {
      SxMatrix3<TPrecCoeffG> polarMat;
      polarMat = dielec.compute (waves, fermi, kp, frequency, structure.cell, 
                                 wings);

      // --- save wings of polarizability
      if (wingFile.getSize () > 0 && wings)  {
         RhoR wingR(3);
         for (int iDir = 0; iDir < 3; iDir++)  {
            wingR(iDir) = dielec.wings(0).colRef(iDir);
            wingR(iDir).setBasis (&R);
         }
         SxRho (wingR).writeRho(wingFile);
      }

      cout << "The dielectric tensor without local field effects is: " << endl;
      for (int a = 0; a < 3; ++a)  {
         for (int b = 0; b < 3; ++b)
            if (fabs(polarMat(a,b).re) < 5e-7) polarMat(a,b).re = 0.;
         sxprintf ("%9.6f %9.6f %9.6f\n", polarMat(a, 0).re,
                                          polarMat(a, 1).re,
                                          polarMat(a, 2).re);
      }
      cout << "The average is " << (polarMat.trace ().re / 3.) << endl;
   } else {
      cout << "Number of frequencies :";
      cin >> nOmega; 
      cout << endl;
      if (nOmega < 1)  {
         cout << "Invalid number of frequency points." << endl;
         SX_QUIT;
      }
      frequencies.resize (nOmega);
      for (int i = 0; i < nOmega; i++)  {
         cout << "freq. no. " << (i+1) << ": ";
         cin >> frequencies(i);
         cout << endl;
      }
      dielec.compute (waves, fermi, kp, frequencies, structure.cell, wings);

      if (outputFile.getSize () == 0)  {
         cout << "The dielectric tensor without local field effects is: ";
         cout << endl;
         for (int i = 0; i < nOmega; i++)  {
            cout << frequencies(i) << ": ";
            for (int a = 0; a < 3; ++a)  {
               sxprintf ("%9.6f %9.6f %9.6f\n", dielec.head(i)(a, 0),
                                                dielec.head(i)(a, 1),
                                                dielec.head(i)(a, 2));
            }
            cout << "The average is ";
            cout << (dielec.head(i).trace () / 3.) << endl;
         }
      }
         
   }

   // timer for k.p
   kp.printTimer ();
      
   FILE *outFile = NULL;
   if (outputFile.getSize () > 0)  {
      outFile = fopen(outputFile.ascii (), "w");
      if (!outFile)  {
         cout << "Can't open file '" << outputFile << "'." << endl;
         SX_EXIT;
      }
      fprintf(outFile,"%6u%6u\n",nOmega,G.ng);
   }

   if (wings)  {
      SxArray<SxDiracVec<TPrecCoeffG> > wingsG(3);
      // TODO: implement SxVector<Complex> /= SxVector<Double>
      // SxDiracVec<Double> g2;
      SX_CHECK (G.g2.getSize () == G.ng, G.g2.getSize (), G.ng);
      SX_CHECK (G.n123(0)(0) == 0, G.n123(0)(0));
      SxVector3<Int> vec;
      SxDiracVec<TPrecCoeffG> wingInR;
      int ng = G.ng, iDir, ig, iw;
      SxDiracVec<TPrecG> absGInv(ng);
      SxIdx nonZero(1,ng-1);
      absGInv(0) = 0.;
      absGInv(nonZero) = sqrt(1. / G.g2(nonZero) );

      for (iw = 0; iw < nOmega; iw++)  {
         cout << "iw = " << (iw+1) << endl;
         if (outFile) {
            fprintf(outFile,"%20.10f\n",frequencies(iw));
            fprintf(outFile,"%20.10f",dielec.head(iw)(0,0));
            fprintf(outFile,"%20.10f",dielec.head(iw)(0,1));
            fprintf(outFile,"%20.10f",dielec.head(iw)(0,2));
            fprintf(outFile,"%20.10f",dielec.head(iw)(1,1));
            fprintf(outFile,"%20.10f",dielec.head(iw)(1,2));
            fprintf(outFile,"%20.10f\n",dielec.head(iw)(2,2));
         }
         for (iDir = 0; iDir < 3; iDir++)  {
            // dielec.wings contains only imaginary part
            wingInR = I * dielec.wings(iw).colRef(iDir);
            wingInR.setBasis (&R);
            wingsG(iDir) = (G | wingInR );
            wingsG(iDir) *= absGInv;
         }

         for (ig = 0; ig < ng; ig++)  {
            vec = G.fft3d(0).mesh.getMeshVec (G.n123(0)(ig), SxMesh3D::Origin);
            if (outFile)
               fprintf(outFile,"%5i%5i%5i\n",vec(0),vec(1),vec(2));
            else
               cout << vec << ": ";
            for (iDir = 0; iDir < 3; iDir++)
               if (outFile)  {
                  fprintf(outFile, "%20.10e%20.10e", wingsG(iDir)(ig).re,
                          wingsG(iDir)(ig).im);
               } else {
                  cout << ' ' << wingsG(iDir)(ig);
               }
            if (outFile) 
               fprintf(outFile,"\n");
            else
               cout << endl;
         }
         cout << endl;

      } // iw

   } // computeWing
   if (outFile) fclose(outFile);

}

#endif /* SX_STANDALONE */
