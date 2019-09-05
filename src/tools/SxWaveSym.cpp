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
#include <SxBinIO.h>
#include <SxPW.h>
#include <SxFermi.h>
#include <SxPAWOverlap.h>
#include <SxPWOverlap.h>
#include <SxRotation.h>

// author: C. Freysoldt
/** 
   Applying a symmetrie operation S in r-space corresponds 
   to applying the transpose of the transformation matrix \f$\mathbf S\f$
   in r-space to the G-vectors:

   \f[ S f(R) := f(\mathbf{S} R) \f]
   \f[ f(R) = \sum c_G \cdot e^{G \cdot R} \f]
   \f[ S f(R) = \sum c_G \cdot e^{G \cdot (\mathbf{S} R)} \f]
   Using the definition of an adjoint operator we get immediately
   \f[ G \cdot (\mathbf{S} R) = (\mathbf{S}^\dagger G) \cdot R \f]
   Because we are working with real matrices, the adjoint matrix is the
   transpose of the matrix.

   Using a non-orthogonal basis aMat for R doesn't change the formulae
   if a corresponding basis bMat is used, where
   \f[ (aMat)^T = \lambda \cdot (bMat)^{-1} ~ \lambda \epsilon R \f]
   i.e. if R = aMat R', G = bMat G' we can replace R by R', G by G' and
   \f$\mathbf S\f$ by the symmetry matrix \f$\mathbf S'\f$ for the 
   relative coordinates:
   \f[ \mathbf S' = aMat \mathbf S aMat^{-1} \f] 
   the formulae above are still fine. This is easily shown by writing out
   all the necessary basis transformations.

   The symmetry operations in G space are applied in the form of mappings:
   If we have a list of G vectors (our basis), we have
   \f[ \mathbf S G_i = G_j \f]
   and all we have to store is the mapping i->j.

   If we have a G+k basis, the whole stuff becomes
   \f[ \mathbf S (G + k) = \mathbf S G + \mathbf S k = G' + \Delta G + k \f]
   where \f$\Delta G\f$ is the shift of k after the symmetry transformation.
   Only symmetry operations where \f$\Delta G\f$ is a lattice vector are
   allowed.
   

*/

int main (int argc, char** argv)  {

   // --- command line parsing
   SxCLI cli(argc, argv);

   SxString wavesFile = cli.option ("-w", "file", "S/PHI/nX waves file")
                    .toString ("waves.sxb");
   int ik = cli.option ("-k", "number", "k-point number").toInt (1,1) - 1;

   int pawGroup = cli.newGroup ("PAW");
   SxString inputFile = cli.option ("-i", "file", "SPHInX input file (for PAW)")
                     .toString ("input.sx");
   bool nonSymmorphic
      = cli.option ("-t|--nonsymmorphic", "include nonsymmorphic symmetries")
        .toBool ();

   cli.finalize ();

   initSPHInXMath ();

   // --- read waves
   SxAtomicStructure structure;
   SxPtr<SxGkBasis> gkBasisPtr;
   SxFermi fermi;
   SxMesh3D mesh;
   SxSpeciesData specData;

   try  {
      SxBinIO io(wavesFile, SxBinIO::BINARY_READ_ONLY);
      structure.read (io);
      gkBasisPtr = SxPtr<SxGkBasis>::create (io, false);
      fermi.read (io);
      fermi.kpPtr = &*gkBasisPtr;
      specData.read (io);
      io.read ("meshDim", &mesh);
      io.close ();
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
   structure.print (specData);
   SxPW waves(wavesFile, SxPW::ReadOnDemand);
   waves.setGkBasisPtr (gkBasisPtr);
   SxGkBasis &gkBasis = *gkBasisPtr;

   SxPtr<SxOverlapBase> SPtr;
   if (cli.groupAvailable (pawGroup))  {
      SxParser parser;
      SxConstPtr<SxSymbolTable> table = parser.read (inputFile);
      if (table->containsGroup("pawPot"))   {
         SxPtr<SxPAWPot> pawPotPtr =SxPtr<SxPAWPot>::create (&*table);
         SxRadBasis radBasis(pawPotPtr->rad, pawPotPtr->logDr);
         gkBasis.changeTau (structure);
         SxPtr<SxPartialWaveBasis> pBasis 
            = SxPtr<SxPartialWaveBasis>::create (pawPotPtr, structure);
         pBasis->createProjBasis (gkBasis);
         SPtr = SxPtr<SxPAWOverlap>::create (pBasis, pawPotPtr);
      } else {
         cout << "WARNING: no pawPot group found in " << inputFile << endl;
         cout << "Will use standard plane-wave norm" << endl;
      }
   }
   // --- use plane-wave overlap
   if (!SPtr) SPtr = SxPtr<SxPWOverlap>::create ();
   // --- setup symmetry mapping in G space
   SxCell recCell = structure.cell.getReciprocalCell ();
   int iSym, ig;

   if (nonSymmorphic) structure.updateSymGroup ();
   structure.cell.symGroupPtr->print ();

   // --- look for symmetry operations compatible with kVec
   // i.e.  SymOp ^ k = G' + k, with lattice vector G'

   int nOp = structure.cell.symGroupPtr->getSize ();
   int ng = gkBasis(ik).ng;
   SxVector3<TPrecG> kVec, kTransformed;
   // get kVec in relative coordinates
   kVec = recCell.carToRel (gkBasis.getK(ik));
   SxVector3<Int> gShift; // G' in SymOp ^ k = G' + k
   SxList<SxMatrix3<Int> > symOpList; // transformed into array later
   SxList<int> symIdList;
   SxList<SxVector3<Int> > gShiftList; // to array later

   SymMat S;
   for (iSym = 0; iSym < nOp; iSym++)  {
      S = recCell.carToRel (structure.cell.symGroupPtr->getRot(iSym));
      kTransformed = S ^ kVec;
      kTransformed -= kVec; // get G'
      SX_LOOP(iDir) gShift(iDir) = (int)round(kTransformed(iDir));
      //gShift = SxVector3<Int> (kTransformed);
      cout << gShift << "==" << kTransformed << "? ";
      if ((gShift - kTransformed).normSqr () < 1e-12)  {
         symOpList  << S;
         symIdList  << (iSym + 1);
         gShiftList << gShift;
         cout << "yes";
      } else {
         cout << "no";
      }
      cout << endl;
   }

   // --- arrays of symmetry operations and shift vectors
   nOp = int(symOpList.getSize ());
   SX_CHECK (gShiftList.getSize () == nOp, gShiftList.getSize (), nOp);
   SxArray<SxMatrix3<Int> > symOp (symOpList);
   SxArray<SxVector3<Int> > gShifts (gShiftList);
      
   // --- get symmetry operation mappings
   SxArray<SxArray<ssize_t> > symInIdx (nOp); // mapping
   SxArray<PsiG> transPhase (nOp); // translation phase for non-symmorphic ops

   SxVector<Int> idxMesh (mesh.getSize ());
   SxVector<TPrecFFTIdx> &n123 = gkBasis(ik).n123(0);
   idxMesh.set (-1); // DEBUG
   for (ig = 0; ig < ng; ig++) idxMesh(n123(ig)) = ig;
   
   SxVector3<Int> gVecRel, result; // G before and after transformation
   for (iSym = 0; iSym < nOp; iSym++)  {
      SxArray<ssize_t> &symIdx = symInIdx(iSym);
      SxVector3<Double> trans = (*structure.cell.symGroupPtr)(iSym).trans;
      symIdx.resize (ng);
      transPhase(iSym).resize (ng);
      for (ig = 0; ig < ng; ig++)  {
         // --- get x, y, z from n123 (from ig)
         gVecRel = mesh.getMeshVec (n123(ig), SxMesh3D::Origin);
         // symmetry transformation
         result = (symOp(iSym) ^ gVecRel) + gShifts(iSym);
         // --- convert x,y,z -> n123 (getFFTIdx) and n123 -> jg (idxMesh)
         symIdx(ig) = idxMesh(mesh.getMeshIdx (result, SxMesh3D::Origin));
         SX_CHECK (symIdx(ig) >= 0, iSym, ig);
         transPhase(iSym)(ig) = exp (-I * gkBasis(ik).getG(ig) ^ trans);
      }
   }
   idxMesh.resize (0);

   // --- calculate symmetry characters of wavefunctions
   PsiG psi, psiRot;
   SxComplex16 ovl = 0;
   int from = 0, iState; 
   double epsOld = fermi.eps(0, 0, ik), eps;// deal with degenerated states

   // Overlap operator
   SxOverlap O(SPtr);
   for (iSym = 0; iSym < nOp; iSym++)  {
      cout << "iSym = " << symIdList(iSym) << ": " << endl;
      cout << SxRotation(structure.cell.symGroupPtr->getRot(symIdList(iSym)-1)).getName () << " ";
      Coord t = structure.cell.carToRel((*structure.cell.symGroupPtr)(symIdList(iSym)-1).trans);
      SX_LOOP(iDir) if (fabs(t(iDir)) < 1e-10) t(iDir) = 0.;
      cout << " t=" << t << endl;
      for (iState = 0; iState < waves.getNStates (); iState++)  {
         eps = fermi.eps (iState, 0, ik);
         if (fabs(epsOld - eps) > 1e-6 / HA2EV)  {
            // no degeneration
            ovl = 0;
            from = iState;
            epsOld = eps;
         }
         psi = waves(iState, 0, ik);
         psiRot = psi.getSorted (symInIdx(iSym));
         psiRot *= transPhase(iSym);
         ovl += (psi | O | psiRot);
         // suppress output if next state (if any) has same energy
         if ((iState + 1) != waves.getNStates ()
             && fabs(fermi.eps (iState + 1, 0, ik) - eps) <= 1e-6 / HA2EV) 
            continue;
         // --- output
         cout << (from + 1);
         if (from < iState) cout << "-" << (iState + 1);
         cout << ": ";
         if  (ovl.im > 1e-8) 
            cout << ovl; 
         else 
            cout << ((fabs(ovl.re) < 1e-12) ? 0. : ovl.re); 
         cout << endl;
      }
   }
}
