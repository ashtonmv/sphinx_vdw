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

#include <SxCLI.h>
#include <SxPW.h>
#include <SxBinIO.h>
#include <SxFermi.h>
#include <SxPseudoPot.h>
#include <SxPerturbK.h>
#include <SxParser.h>


int main (int argc, char **argv)
{

   initSPHInXMath ();

   SxCLI cli (argc,argv);
   cli.authors = "C. Freysoldt";
   cli.preUsageMessage = 
      "This add-on interpolates bandstructure via perturbation theory.\n"
      "The first-order derivative in H is exact, the 2-nd order ignores"
      " the contribution of the non-local pseudopotential. The interpolated "
      "band structure may therefore have parabolic noise.";
   SxString wavesFile = cli.option ("-w|--waves","file","waves file")
                        .toString("waves.sxb");

   int blockSize
      = cli.option ("--blocksize","set blocksize. Smaller blocksizes may reduce "
                    "the memory consumption, larger blocksizes may speed up the "
                    " <phi|psi> part").toInt (64);
   bool linear 
      = cli.option ("--linear","use only linear expansion. For flat bands you "
                    "can see a quite drastic quadratic error. This might be "
                    "to a missing cancellation for the quadratic terms, since "
                    "only the kinetic part is taken into account.").toBool ();

   SxString inputFile = cli.option ("-i|--input","file","S/PHI/nX input file")
                        .toString ("input.sx");

   cli.finalize ();

   SxPW waves(wavesFile, SxPW::ReadOnDemand);
   SxFermi fermi;

   SxAtomicStructure structure;

   try {
      SxBinIO io (wavesFile, SxBinIO::BINARY_READ_ONLY);
      fermi.read (io);
      structure.read(io);
      waves.setGkBasisPtr (SxPtr<SxGkBasis>::create (io));
      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   SxGkBasis &gkBasis = waves.getGkBasis();
   fermi.kpPtr = &gkBasis;

   int nStates = waves.getNStates ();
   int nk = waves.getNk (), nSpin = waves.getNSpin ();
   for (int ik = 0; ik < nk; ++ik)
      gkBasis(ik).memMode = SxGBasis::SaveMemory;

   SxDiracMat<TPrecCoeffG> npn(nk,3);

   // --- initialize k dot p
   SxPerturbK kp;
   SxKPoints kNew;
   {
      SxParser parser;
      SxParser::Table table = parser.read(inputFile);
      SxPseudoPot psPot(&*table);
      kNew = SxKPoints(structure.cell, &*table);
      gkBasis.changeTau(structure);
      kp.set (psPot, gkBasis, structure);
      kp.blockSize = blockSize;
   }

   int nkNew = kNew.getNk ();
   SxCell recCell = structure.cell.getReciprocalCell ();
   SxAtomicStructure dk(recCell, nkNew), kOld(recCell, nk), dkSym;
   for (int ik = 0; ik < nk; ++ik)
      kOld.setAtom(ik, gkBasis.getK (ik));

   // --- find out which old k-points are used for the new k-points
   Coord k,deltaK;
   double dk2;
   SxArray<SxList<int> > whichNewK(nk);
   int nSym = recCell.symGroupPtr->getNSymmorphic (), idMin, whichOldK=-1;
   for (int ik = 0; ik < nkNew; ++ik)  {
      k = kNew.getK (ik);
      dk2 = 1e10; // huge
      // --- find nearest old k-point
      for (int iSym = 0; iSym < nSym; ++iSym)  {
         dkSym = (recCell.symGroupPtr->getSymmorphic(iSym) ^ k) - kOld;
         dkSym.map (recCell, SxCell::WignerSeitz);
         if (dk2 > dkSym.absSqr ().minval (&idMin))  {
            deltaK = dkSym.getAtom(idMin);
            dk.setAtom(ik,deltaK);
            dk2 = deltaK.normSqr ();
            whichOldK = idMin;
         }
      }
      whichNewK(whichOldK).append (ik);
      cout << "new k " << (ik+1) 
           << " will be computed from old k " << (whichOldK+1) 
           << " (|dk|=" << sqrt(dk2) << ")." << endl;
   }

   
   // --- compute band energies
   SxFermi fermiNew(0., nStates, nSpin, kNew);

   SxDiracMat<TPrecCoeffG> kpMatElem;
   SxDiracMat<TPrecCoeffG> ham(nStates, nStates);
   SxList<int>::Iterator itK, itKEnd;
   for (int ik = 0; ik < nk; ik++)  {
      (cout << "old ik = " << (ik+1) << endl).flush ();
      for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
         if (whichNewK(ik).getSize () == 0) continue;
         PsiGI psi = waves(iSpin,ik);
         kpMatElem = kp.getMatrixElements (psi, psi);
         itKEnd = whichNewK(ik).end ();

         for (itK = whichNewK(ik).begin (); itK != itKEnd; ++itK)  {
            cout << "computing new k " << (*itK + 1) << endl;
            ham.set (0.);
            // --- diagonal
            for (int i = 0; i < nStates; ++i)   {
               ham(i,i) = fermi.eps(i,iSpin,ik);
               // 2nd order kinetic energy correction
               if (!linear) ham(i,i) += + 0.5 * dk(*itK).normSqr ();
            }
            // --- k dot p correction
            for (int iDir = 0; iDir < 3; ++iDir)
               ham.plus_assign_ax (dk(*itK)(iDir), kpMatElem.colRef(iDir));
            // --- diagonalize ham and save eigenvalues
            fermiNew.eps(iSpin,*itK) = ham.eigenvalues ();
         }

      }
      
   }

   kp.printTimer ();
   fermiNew.writeSpectrum("eps_kdotp","dat");
   
}
