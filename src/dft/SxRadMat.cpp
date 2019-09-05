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

#include <SxRadMat.h>
#include <SxGrid.h>
#include <SxPAWPot.h>
#include <SxPartialWaveBasis.h>

/*
References:
 1. P. E. Bloechl, Phys. Rev. B 50, 17953 (1994).
 2. C. Freysoldt, PAW implementation notes
*/

void SxRadMat::resize (const SxConstPtr<SxAtomInfo> &info,
                       int nSpin)
{
   SX_CHECK (info);
   atomInfo = info;
   int nSite = atomInfo->nAtoms.sum ();
   rho.reformat (nSpin, nSite);
}

void SxRadMat::resize (const SxConstPtr<SxAtomInfo> &info,
                       const SxPAWPot &pawPot,
                       int nSpin)
{
   SX_CHECK (pawPot.getNSpecies () == info->nSpecies,
             pawPot.getNSpecies (), info->nSpecies);
   int nSpecies = pawPot.getNSpecies ();
   resize (info, nSpin);
   for (int iSpin = 0; iSpin < getNSpin (); ++iSpin)  {
      for (int is = 0, iSite = 0; is < nSpecies; ++is)  {
         int nAtoms = info->nAtoms(is);
         int nOrb = pawPot.getNProj (is);
         for (int ia = 0; ia < nAtoms; ++ia, iSite++)  {
            rho(iSpin, iSite).reformat (nOrb, nOrb);
         }
      }
   }
}

void SxRadMat::resize (const SxRadMat &in)
{
   atomInfo = in.atomInfo;
   rho.reformat (in.getNSpin (), in.getNSite ());
   for (int iSpin = 0; iSpin < getNSpin (); ++iSpin)  {
      for (int iSite = 0; iSite < getNSite (); ++iSite)  {
         ssize_t nOrb = in.rho(iSpin, iSite).nRows ();
         SX_CHECK (in.rho(iSpin, iSite).nCols () == nOrb,
                   in.rho(iSpin, iSite).nCols (), nOrb);
         rho(iSpin, iSite).reformat (nOrb, nOrb);
      }
   }
}


void SxRadMat::set (double x)
{
   for (int iSpin = 0; iSpin < getNSpin (); ++iSpin)
      for (int iSite = 0; iSite < getNSite (); ++iSite)
         rho(iSpin, iSite).set (x);
}

void SxRadMat::symmetrize (const SxAtomicStructure &str,
                           const SxPAWPot &pot,
                           const SxYlmRotGroup &ylmRot)
{
   SX_CHECK (str.cell.symGroupPtr);
   const SxSymGroup &syms = *str.cell.symGroupPtr;
   int nSym = syms.getNSymmorphic ();
   SxVector<Int> symAtomIdx(nSym);
   SxGrid grid (str, 10);

   for (int is = 0; is < str.getNSpecies (); ++is)  {
      for (int ia = 0; ia < str.getNAtoms (is); ++ia)  {
         // --- find equivalent atoms
         for (int iSym = 0; iSym < nSym; ++iSym)  {
            symAtomIdx(iSym) = str.find (syms.getSymmorphic(iSym) 
                  ^ str.getAtom(is,ia), grid)
                  -str.atomInfo->offset(is);
         }
         // check if symmetrization has been performed already
         if (symAtomIdx.minval () < ia) continue;

         for (int iSpin = 0; iSpin < getNSpin (); ++iSpin)  {
            // --- loop of partial wave types
            for (int ipt = 0; ipt < pot.getNProjType (is); ++ipt)  {
               int l = pot.lPhi(is)(ipt), nm = 2 * l + 1;
               int offset = pot.offset(is)(ipt);
               for (int ipt2 = 0; ipt2 < pot.getNProjType (is); ++ipt2)  {
                  int l2 = pot.lPhi(is)(ipt2), nm2 = 2 * l2 + 1;
                  int offset2 = pot.offset(is)(ipt2);

                  // --- construct symmetrized Dij (type1,type2)
                  SxMatrix<Double> sym (nm, nm2), rotRho(nm,nm2);
                  sym.set (0.);
                  for (int iSym = 0; iSym < nSym; ++iSym)  {
                     // --- collect rotRho
                     SxMatrix<Double> &rotD 
                        = (*this)(iSpin, is, symAtomIdx(iSym));
                     const SxMatrix<Double> &Dl  = ylmRot(iSym)(l);
                     const SxMatrix<Double> &Dl2 = ylmRot(iSym)(l2);
                     for (int m = 0; m < nm; ++m)
                        for (int m2 = 0; m2 < nm2; ++m2)
                           rotRho(m,m2) = rotD(offset + m, offset2 + m2);
                     // --- symmetrize 
                     sym += Dl.transpose () ^ rotRho ^ Dl2;
                  }
                  sym /= double(nSym);
                  /*
                  if (l == 1 && l2 == 1)  {
                     cout << "is=" << is << "; ia=" << ia << endl;
                     cout << "l=1 matrix" << endl;
                     cout << sym << endl;
                     cout << sym.eigenvalues () << endl;
                     cout << sym.eigenvectors () << endl;
                  }
                  */

                  // --- distribute symmetrized matrix
                  for (int iSym = 0; iSym < nSym; ++iSym)  {
                     SxMatrix<Double> &rotD 
                        = (*this)(iSpin, is, symAtomIdx(iSym));
                     const SxMatrix<Double> &Dl  = ylmRot(iSym)(l);
                     const SxMatrix<Double> &Dl2 = ylmRot(iSym)(l2);
                     // rotate with symmetry
                     rotRho = Dl ^ sym ^ Dl2.transpose ();
                     // distribute
                     for (int m = 0; m < nm; ++m)
                        for (int m2 = 0; m2 < nm2; ++m2)
                           rotD(offset + m, offset2 + m2) = rotRho(m,m2);
                  }
               }
            }
         }
      }
   }
}

SxDiracVec<Complex16>
SxRadMat::operator^ (const SxDiracVec<Complex16> &p) const
{
   // TODO: check loop order with respect to performance
   //       iState loop should probably outmost
   // or try alternative
   // loop is,ia 
   //   SxIdx range(offset, offset+npl-1)   
   //   setup complex Hij matrix
   //   res.rowRef (range) <<= Hij ^ p.rowRef (range)
   // => need rowRef as submatrix class
   SX_CHECK (p.handle);
   const SxPartialWaveBasis &pBasis = p.getBasis<SxPartialWaveBasis> ();
   SxConstPtr<SxPAWPot> potPtr = pBasis.getPotPtr ();
   int iSpin = p.handle->auxData.iSpin;
   SX_CHECK(iSpin >= 0 && iSpin < getNSpin (), iSpin, getNSpin ());

   SxDiracMat<Complex16> res(p.nRows (), p.nCols ());
   res.set (0.);
   res.setBasis (&pBasis);
   res.handle->auxData = p.handle->auxData;

   int offset = 0;
   for (int is = 0; is < getNSpecies (); ++is)  {
      int nProjLocal = potPtr->getNProj (is);

      // --- set up i^l 
      SxVector<Complex16> il(nProjLocal); // i^(-l)  (:ip)
      for (int ipt = 0, ipl = 0; ipt < potPtr->lPhi(is).getSize (); ++ipt)  {
         int l = potPtr->lPhi(is)(ipt);
         SxComplex16 c = (l & 1) ? I : SxComplex16(1.);
         if (l & 2) c=-c;
         for (int m = -l; m <=l; ++m, ++ipl)  {
            il(ipl) = c;
         }
      }

      // --- loop over atoms
      for (int ia = 0; ia < getNAtoms(is); ++ia, offset+=nProjLocal)  {

         const SxMatrix<Double>& Hloc = (*this)(iSpin,is,ia);
         // --- loops over projectors
         for (int ipl = 0; ipl < nProjLocal; ++ipl)  {
            for (int jpl = 0; jpl < nProjLocal; ++jpl)  {

               // potential energy * i^l projector phases (Ref. 1)
               SxComplex16 Hij = Hloc(ipl,jpl) 
                               * (il(ipl) * il(jpl).conj ());
               
               // compute sum_j H(R)ij <pj|psi>
               for (int iState = 0; iState < p.nCols (); ++iState)  {
                  // + projector i^l phase-correction 
                  res(offset+ipl, iState) += Hij * p(offset+jpl,iState);
               }
            }
         }
      }
   }
   return res;
}

double tr (const SxRadMat &D, const SxRadMat &A)
{
   SX_CHECK (D.getNSpin () == A.getNSpin (), D.getNSpin (), A.getNSpin ());
   SX_CHECK (D.getNSite () == A.getNSite (), D.getNSite (), A.getNSite ());
   double res = 0.;
   for (int iSpin = 0; iSpin < D.getNSpin (); ++iSpin)
      for (int iSite = 0; iSite < D.getNSite (); ++iSite)
         res += dot(D.rho(iSpin, iSite), A.rho(iSpin,iSite));
   SX_CHECK_NUM (res);
   return res;
}

void SxRadMat::readRho (const SxBinIO &io)
{
   SxPtr<SxAtomInfo> info;
   try {
      int nSpecies = io.getDimension ("nSpecies");
      int nSpin    = io.getDimension ("nMeshes");
      info = info.create (nSpecies);
      SxVector<Int> npl(nSpecies);

      io.read ("nAtoms", &info->nAtoms, nSpecies);
      info->setupOffset ();
      io.read ("npl", &npl, nSpecies);

      resize (info, nSpin);
      for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
         int offset = 0;
         for (int is = 0; is < nSpecies; ++is)  {
            int N = npl(is);
            for (int ia = 0; ia < info->nAtoms(is); ++ia)  {
               SxMatrix<Double> &Dij = (*this)(iSpin, is, ia);
               Dij.reformat (N, N);
               io.readRow ("Dij", &Dij, iSpin, offset);
               offset += N*N;
            }
         }
      }
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

}

void SxRadMat::writeRho (SxBinIO &io) const
{
   int nSpecies = getNSpecies ();
   SxVector<Int> npl(nSpecies);
   int nCoeff = 0;
   // --- get dimensions
   for (int is = 0; is < getNSpecies (); ++is) {
      npl(is) = int( (*this)(0,is,0).nCols () );
      nCoeff += getNAtoms(is) * npl(is) * npl(is);
   }

   try  {
      io.addDimension ("nSpecies", nSpecies);
      if (io.ncMode == SxBinIO::WRITE_DATA || !io.contains ("nAtoms"))
         io.write ("nAtoms", atomInfo->nAtoms, "nSpecies");
      io.write ("npl", npl, "nSpecies");
      io.addDimension ("nDij", nCoeff);
      if (io.ncMode == SxBinIO::WRITE_HEADER)  {
         io.addDoubleVar ("Dij", "nMeshes", "nDij");
      } else {
         for (int iSpin = 0; iSpin < getNSpin (); ++iSpin)  {
            int offset = 0;
            for (int is = 0; is < nSpecies; ++is)  {
               for (int ia = 0; ia < atomInfo->nAtoms(is); ++ia)  {
                  io.writeRow ("Dij", (*this)(iSpin,is,ia), iSpin, offset);
                  offset += npl(is) * npl(is);
               }
            }
         }
      }
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

}
