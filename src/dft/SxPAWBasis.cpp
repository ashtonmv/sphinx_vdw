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

#include <SxPAWBasis.h>
#include <SxRBasis.h>
#include <SxAOBasis.h>

void SxPAWBasis::print () const
{
   gBasis->print ();
   pBasis->print ();
}

SxDiracVec<TGBasisType>
SxPAWBasis::toPWBasis (const SxGBasis *gIn,
                       const SxDiracVec<Complex16> &psi) const
{
   SX_CHECK (gBasis);
   SX_CHECK (gIn == gBasis.getPtr ());
   int nPsi = (int)psi.nCols ();
   if (nPsi == 0) nPsi = 1;
   SX_CHECK (psi.getSize () > 0);
   SxDiracVec<TGBasisType> res;
   res.reformat (gBasis->getNElements (), nPsi);
   res.handle->auxData = psi.handle->auxData;
   res.setBasis (gIn);
   int N = (int)getNElements ();
   // --- copy G-vector part into new vector
   for (int n = 0; n < nPsi; ++n)
      res.colRef (n) <<= psi( SxIdx(N * n, N * n + gBasis->ng - 1) );
   return res;
}

SxDiracVec<TGBasisType>
SxPAWBasis::toPartials (const SxPartialWaveBasis *pIn,
                        const SxDiracVec<Complex16> &psi) const
{
   SX_CHECK (gBasis);
   SX_CHECK (pBasis);
   SX_CHECK (pIn == pBasis.getPtr ());
   int nPsi = (int)psi.nCols ();
   if (nPsi == 0) nPsi = 1;
   SX_CHECK (psi.getSize () > 0);
   SxDiracVec<TGBasisType> res;
   res.reformat (pBasis->getNElements (), nPsi);
   res.handle->auxData = psi.handle->auxData;
   res.setBasis (pIn);
   int N = (int)getNElements ();
   // --- copy partial wave part into new vector
   for (int n = 0; n < nPsi; ++n)
      res.colRef (n) <<= psi( SxIdx(N * n + gBasis->ng, N * (n+1) - 1) );
   return res;
}

SxDiracVec<TGBasisType>
SxPAWBasis::toRBasis (const SxRBasis *rIn,
                      const SxDiracVec<Complex16> &psi) const
{
   ssize_t nPsi = psi.nCols ();
   if (nPsi == 0) nPsi = 1;
   SX_CHECK (psi.getSize () > 0);
   SxDiracVec<TGBasisType> res;
   res.reformat (rIn->getNElements (), nPsi);
   res.handle->auxData = psi.handle->auxData;
   res.setBasis (rIn);
   int N = (int)getNElements ();
   // get G coefficients -> 
   for (int n = 0; n < nPsi; ++n)  {
      SxIdx gPart(N * n, N * n + gBasis->ng - 1);
      SxDiracVec<TGBasisType> gCoeff = psi(gPart);
      res.colRef (n) <<= gBasis->toRealSpace (rIn, gCoeff);
   }
   return res;
}


SxComplex16 SxPAWBasis::scalarProduct (const SxDiracVec<Complex16> &x,
                                       const SxDiracVec<Complex16> &y) const
{
   SX_CHECK (x.getSize () == getNElements (),
             x.getSize (), getNElements ());
   SxIdx gVecs (0, gBasis->ng - 1);
   if (&x == &y) return x(gVecs).normSqr ();
   return dot (x(gVecs), y(gVecs));
}

SxDiracVec<Complex16> SxPAWBasis::toAO (const SxAOBasis *basis,
                                        const SxDiracVec<Complex16> &in) const
{
   int ik = in.handle->auxData.ik;
#ifndef NDEBUG
   int nk = (int)basis->refOrbitals.getSize ();
   SX_CHECK (ik >= 0 && ik < nk, ik, nk);
#endif
   if (basis->SPtr) return basis->fromPWBasis(basis->SPtr->apply(in), ik);
   return basis->fromPWBasis(in, ik);
}

SxDiracVec<Complex16> SxPAWBasis::toAny (const SxBasis *basis,
                                       const SxDiracVec<Complex16> &in) const
{
   // --- identity ?
   if (dynamic_cast<const SxPAWBasis*> (basis))  {
      SX_CHECK (basis == this);
      return in;
   }
   // --- G basis?
   if (const SxGBasis* gk = dynamic_cast<const SxGBasis*> (basis))
      return toPWBasis (gk, in);

   if (const SxAOBasis* ao = dynamic_cast<const SxAOBasis*> (basis))
      return toAO (ao, in);

   SX_EXIT; return in;
}
                       

