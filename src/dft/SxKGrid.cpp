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

#include <SxKGrid.h>

SxKGrid::SxKGrid ()
{
   /* empty */
}


SxKGrid::SxKGrid (const SxKPoints &in, const SxCell &realSpaceCell)
   : SxKPoints (in)
{
   setup (realSpaceCell);
}

SxKGrid::SxKGrid (const SxAtomicStructure &str,
                  const SxSymbolTable *table)
   : SxKPoints (str.cell, table)
{
   setup (str.cell);
}

SxKGrid::SxKGrid (const SxBinIO &io, const SxCell &realSpaceCell)
{
   SxKPoints::read (io);
   setup (realSpaceCell);
}


void SxKGrid::setup (const SxCell &cell)
{
   SX_CHECK (cell.volume > 0.);
   SX_CHECK (cell.symGroupPtr);
   SxCell rCell (cell.getReciprocalCell ());
   const SxArray<SymMat> S = cell.symGroupPtr->getSymmorphic ();
   int iSym, nSym = int(S.getSize ());
   SX_CHECK (nSym > 0);
   Coord k, k0 = getK(0), diff;

   // --- get Monkhorst-Pack grid points by rotating symmetrized k-points
   SxList<Coord> lattice;
   for (int ik = 0; ik < nk; ++ik)  {
      k = getK(ik) - k0;
      for (iSym = 0; iSym < nSym; ++iSym)  {
         diff = S(iSym) ^ k;
         rCell.map (&diff, SxCell::Origin);
         lattice << diff;
      }
   }
   // and add the reciprocal cell vectors (needed for folding = 1)
   lattice << rCell.col(0) << rCell.col(1) << rCell.col(2);

   // --- get Monkhorst-Pack lattice from points
   int nBasis = subCell.setFromVectors(lattice);
   SX_CHECK (nBasis == 3);
   if (nBasis != 3)  {
      cout << "SxKGrid: failed to process k-points..." << endl;
      SX_EXIT;
   }
    
   // --- get the folding: maximum component in relative coordinates
   SxVector3<Int> rel;
   int i,j;
   foldingMP = 0;
   for (i = 0; i < 3; ++i)  {
      rel = subCell.carToRel (rCell.col(i));
      if ((subCell.relToCar(rel) - rCell.col(i)).normSqr () > 1e-6)  {
         cout << "Can't identify k<->R fft mesh." << endl;
         cout << "Guessed k-mesh cell: " << subCell << endl;
         cout << "reciprocal cell basis vector " << (i+1) << " isn't integer "
                 "but" << endl << subCell.carToRel (rCell.col(i)) << endl;
         SX_QUIT;
      }
      for (j = 0; j < 3; ++j)
         if (foldingMP(j) < abs(rel(j)) )
             foldingMP(j) = abs(rel(j));
   }

   // --- setup FFT
   int meshSize = foldingMP.product ();
   fftRk.dir = SxFFT::None;
   fftRk.setMesh (foldingMP, subCell.volume);
   
   // --- determine symmetrized k-point index for each FFT mesh point
   Coord sK;
   symK.resize (meshSize, 0);
   // debug
   SxVector<Int> count(nk); count = 0;
   
   int ik;
   for (i = 0; i < meshSize; ++i)  {
      k = subCell.relToCar(fftRk.mesh.getMeshVec(i, SxMesh3D::Positive)) + k0;
      
      for (iSym = 0; iSym < nSym; ++iSym)  {
         sK = S(iSym) ^ k;
         for (ik = 0; ik < nk; ++ik)  {
            diff = sK - getK(ik);
            rCell.map (&diff, SxCell::Origin);
            if (diff.normSqr () < 1e-7)  {
               symK(i)=(ik+1);
               break;
            }
            diff = sK + getK(ik);
            rCell.map (&diff, SxCell::Origin);
            if (diff.normSqr () < 1e-7)  {
               symK(i)=-(ik+1);
               break;
            }
         }
         if (ik < nk) break;
      }
      // debug
      count (abs(symK(i))-1)++;
   }

   // debug
   cout << (count / weights) << endl;
}

const SxFFT3d& SxKGrid::setupFFT (const SxFFT::Directions dir)
{
   SX_CHECK (fftRk.mesh.product () > 0);
   SX_CHECK (subCell.volume > 0.);
   fftRk.dir = dir;
   fftRk.setMesh (fftRk.mesh, subCell.volume);
   return fftRk;
}

void SxKGrid::setMesh (SxDiracVec<Complex16> *mesh, 
                       const SxDiracVec<Complex16> &values)
{
   SX_CHECK (mesh);
   SX_CHECK (values.getSize () == nk);
   int n = (int)symK.getSize ();
   SX_CHECK (n > 0); // has it been setup??
   SX_CHECK (mesh->getSize () == symK.getSize ());
   int ik;
   SxDiracVec<Complex16>::Iterator it = mesh->begin ();
   for (int i = 0; i < n; ++i, ++it)  {
      if ((ik = symK(i)) > 0)
         *it = values(ik-1);
      else
         *it = values(-ik-1);
   }
}




