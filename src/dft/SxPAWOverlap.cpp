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

#include <SxPAWOverlap.h>
#include <SxProjector.h>
#include <SxPAWBasis.h>
#include <SxPAWHamiltonian.h>

SxComplex16 SxPAWOverlap::dot (const PsiG &x, const PsiG &y) const
{
   SX_CLOCK (Timer::PAWdot);
   SX_CHECK (pBasis);
   SX_START_TIMER (Timer::PAWdotP);
   SxDiracVec<Complex16> px = *pBasis | x,
                         py = *pBasis | y;
   SX_STOP_TIMER (Timer::PAWdotP);
   return (x | y) + dotPartial (px, py);
}

SxComplex16 SxPAWOverlap::dotPartial (const SxDiracVec<Complex16> &px,
                                      const SxDiracVec<Complex16> &py) const
{
   SX_CHECK (pawPot);
   SxComplex16 res = 0.;
   // --- partial wave correction (ref. 1, Eq. (30))
   int nProjGlobal = (int)px.getSize ();
   for (int ipg = 0; ipg < nProjGlobal; ++ipg)  {
      const SxAOBasis::OrbitalIndex &id = pBasis->projectors->orbitalMap(ipg);
      const SxArray<int>   &lPhi = pawPot->lPhi(id.is);
      const SxArray<int> &offset = pawPot->offset(id.is);
      int ipt = pBasis->projectors->refOrbMap(id.is)(id.io).n;
      const SxDiracMat<Double> &dover = pawPot->deltaS(id.is);
      int nProjLocal = int(lPhi.getSize ());
      // --- loop over projectors with same (l,m) on this atom
      for (int jpt = 0; jpt < nProjLocal; ++jpt)  {
         if (lPhi(jpt) == lPhi(ipt))  {
            // same l
            int offJ = offset(jpt) - offset(ipt);
            res += dover(ipt, jpt) * px(ipg).conj () * py(ipg + offJ);
         }
      }
   }
   return res;
}

double SxPAWOverlap::normSqr (const PsiG &x) const
{
   SX_CHECK (pBasis);
   SX_CLOCK (Timer::PAWdot);
   SX_START_TIMER (Timer::PAWdotP);
   SxDiracVec<Complex16> px = *pBasis | x;
   SX_STOP_TIMER (Timer::PAWdotP);
   return (x | x) + dotPartial (px, px);
}

/// Normalize states x
void SxPAWOverlap::normalize (SxDiracVec<Complex16> *xPtr) const
{
   SX_CHECK (xPtr);
   int n = (int)xPtr->nCols ();
   if (n == 0) n = 1;
   SX_CLOCK (Timer::PAWdot);
   SX_START_TIMER (Timer::PAWdotP);
   SxDiracVec<Complex16> px = *pBasis | *xPtr;
   SX_STOP_TIMER (Timer::PAWdotP);

   for (int i = 0; i < n; ++i)  {
      SxDiracVec<Complex16> Xi = xPtr->colRef(i),
                            pi = px.colRef(i);
      double nrm2 = (Xi | Xi) + dotPartial(pi, pi);
      SX_CHECK (nrm2 > 0., nrm2);
      Xi *= 1. / sqrt(nrm2);
   }
}
                         
SxDiracMat<Complex16> 
SxPAWOverlap::getMatrix (const SxDiracVec<Complex16> &x,
                         const SxDiracVec<Complex16> &y) const
{
   SX_CHECK (x.getBasisPtr () == y.getBasisPtr ());
   SX_CHECK (pBasis);
   // --- get G basis
   const SxGBasis *gPtr = dynamic_cast<const SxGBasis *> (x.getBasisPtr ());
   if (!gPtr) gPtr = x.getBasis<SxPAWBasis> ().gBasis.getPtr ();
   SX_CHECK (gPtr);

   if (&x == &y)  {
      // --- self-overlap
      SxDiracVec<TPrecCoeffG> px = *pBasis | x;
      return x.overlap (x, gPtr->ng) + getDeltaS (px, px);
   }

   if (x.nRows () == gPtr->ng)  {
      // --- G basis
      if (x.nCols () == y.nCols () && x.nCols () < 16)  {
         // --- put x and y into one matrix and apply S
         int ng = gPtr->ng, nx = (int)x.nCols (), ny = (int)y.nCols ();
         SxIdx xPart(0, nx * ng - 1),
               yPart(ny * ng, (nx+ny)*ng - 1);
         SxDiracMat<TPrecCoeffG> xy(ng , nx + ny);
         xy.setBasis (*gPtr);
         xy(xPart) <<= x;
         xy(yPart) <<= y;
         SxDiracVec<TPrecCoeffG> pxy = *pBasis | xy;
         return x.overlap (y) 
                + getDeltaS (pxy(xPart).reshape (ng, nx),
                             pxy(yPart).reshape (ng, ny));
      } else if (x.nCols () > y.nCols ())  {
         // apply S to y
         return x.overlap (apply(y));
      } else {
         // apply S to x
         return apply (x).overlap (y);
      }
   }
   // PAW basis
   return x.overlap (y, gPtr->ng) + getDeltaS (*pBasis | x, *pBasis | y);
}

SxDiracMat<TPrecCoeffG> 
SxPAWOverlap::getDeltaS (const SxDiracVec<TPrecCoeffG> &px, 
                         const SxDiracVec<TPrecCoeffG> &py) const
{
   SxDiracMat<Complex16> Sni(px.nCols (), px.nRows ());
   Sni.set (0.);

   int offset = 0;
   int nOrb = int(pBasis->getNElements ());
   for (int is = 0; is < pawPot->getNSpecies (); ++is)  {
      int nProjLocal = pawPot->getNProj (is);

      // get 1-center corrections
      const SxDiracMat<Double> &dS = pawPot->deltaS(is);

      // --- loop over atoms
      for (int ia = 0; 
           offset < nOrb && pBasis->projectors->orbitalMap(offset).is == is;
           ++ia, offset+=nProjLocal)
      {

         // --- loops over projectors
         for (int ipl = 0; ipl < nProjLocal; ++ipl)  {
            const SxAOBasis::AoIndex &idI 
               = pBasis->projectors->refOrbMap(is)(ipl);
            for (int jpl = 0; jpl < nProjLocal; ++jpl)  {
               const SxAOBasis::AoIndex &idJ 
                  = pBasis->projectors->refOrbMap(is)(jpl);
               if (idI.l != idJ.l || idI.m != idJ.m) continue;

               // --- Ref. 1, Eq. 30, |pi>(...)<pj term
               // compute sum_j S(R)ij <pj|psi'>
               // TODO:
               // pxAdj = px.adjoint ()
               // Sni.colRef (offset + ipl)
               // .plus_assign_ax (dS(..), pxAdj.colRef(offset + jpl)
               for (int iState = 0; iState < px.nCols (); ++iState)  {
                  Sni(iState, offset + ipl) += px(offset+jpl,iState).conj () 
                                               * dS(idI.n, idJ.n);
               }
            }

         }
      }
   }
   // compute sum_j [sum_i <psi|pi> S(R)ij] <pj|psi'>]
   return Sni ^ py;
}

void SxPAWOverlap::orthonormalize (PsiGI *psiPtr, 
                                   SxOrthoMethod how) const
{
   SX_CLOCK (Timer::ortho);
   SX_CHECK (psiPtr);
   PsiGI &psi = *psiPtr;
   int nStates = (int)psi.nCols ();

   // --- get G basis
   const SxGBasis *gPtr = dynamic_cast<const SxGBasis *> (psi.getBasisPtr ());
   if (!gPtr) gPtr = psi.getBasis<SxPAWBasis> ().gBasis.getPtr ();
   SX_CHECK (gPtr);

   // --- compute current overlap matrix S
   SxDiracMat<TPrecCoeffG> S = getMatrix (psi, psi);

   if (how == GramSchmidt)  {
      // --- Gram/Schmidt orthogonalization
      SxDiracMat<TPrecCoeffG> U = S.choleskyDecomposition ().adjoint ();
      psi.rotate (U.inverse ());
   }  else  {
      // --- diagonalization, according to LOEWDIN orthogonalization scheme
      SxDiracMat<TPrecCoeffG> U, Ieps;      // S being the overlap matrix,
                                            // U is S^-1/2, and
      SxDiracSymMat<TPrecCoeffG> symS;
      SxDiracSymMat<TPrecCoeffG>::Eigensystem eig;

      symS = SxDiracSymMat<TPrecCoeffG> (nStates);
      for (int iState = 0; iState < nStates; iState++)  {
         for (int jState = iState; jState < nStates; jState++)  {
            symS(iState,jState)  = S(iState,jState);
         }
      }

      eig   = symS.eigensystem();
      Ieps  = Ieps.identity(1. / sqrt(eig.vals) );
      U     = eig.vecs ^ Ieps ^ eig.vecs.adjoint();

      psi.rotate (U);

   }
}

void SxPAWOverlap::setOrthogonal (PsiG *X, const PsiG &orthoStates) const
{
   SX_CHECK (X);
   if (orthoStates.getSize () == 0) return;

   SX_CHECK (X->nRows () == orthoStates.nRows (),
             X->nRows (), orthoStates.nRows ());
   SX_CHECK (X->getBasisPtr () == orthoStates.getBasisPtr ());

   SX_CLOCK (Timer::ortho);

   SxDiracVec<TPrecCoeffG> &psi = *X;
   int nPsi = (int)psi.nCols ();
   if (nPsi == 0) nPsi = 1;
   int nOrtho = (int)orthoStates.nCols ();

   // for low n, the direct (convential) BLAS1 approach is faster
   // because there is a certain overhead, and BLAS2 routines may
   // not be fastest for BLAS1-like tasks
   // the cross-over point should be tested
   // this should be the zdotc / gemv crossover for nPsi == 1
   // this should be the gemv / gemm crossover for nPsi > 1
   bool blas1 = (nOrtho <= 150 && nPsi == 1);
   const SxPAWBasis *pawBasis 
      = dynamic_cast<const SxPAWBasis*> (psi.getBasisPtr ());
   if (pawBasis)  {
      if (blas1) {
         for (int j=0; j < nOrtho; j++)  {
            // psi -= |j><j|S|psi>
            PrecCoeffG scp = dot(orthoStates.colRef(j), psi);
            psi.plus_assign_ax (-scp, orthoStates.colRef(j));
         }
      } else {
         const SxGBasis &gk = *pawBasis->gBasis;
         const SxPartialWaveBasis &p = *pawBasis->pBasis;
         SxDiracVec<TPrecCoeffG> scps;
         //scps = (gk | orthoStates).adjoint () ^ (gk | psi);
         scps = orthoStates.overlap (psi, gk.ng);
         scps += getDeltaS (p | orthoStates, p | psi);
         psi -= (orthoStates ^ scps);
      }
   } else {
      // --- |G+k> basis only: do minimal amount of projections on the fly
      if (nPsi < nOrtho)  {
         // apply S operator to psi
         SxDiracVec<TPrecCoeffG> Spsi = apply (psi);
         if (blas1)  {
            for (int j=0; j < nOrtho; j++)  {
               SxDiracVec<TPrecCoeffG> psiJ;
               // psi -= |j><j|S|psi>
               psiJ = orthoStates.colRef(j);
               psi.plus_assign_ax (-::dot(psiJ, Spsi), psiJ);
            }
         } else {
            SxDiracVec<TPrecCoeffG> scps;
            scps = (Spsi.adjoint () ^ orthoStates).adjoint ();
            psi -= (orthoStates ^ scps);
         }
      } else {
         // apply S operator to orthoStates
         SxDiracVec<TPrecCoeffG> Sortho = apply (orthoStates);
         if (blas1) {
            for (int j=0; j < nOrtho; j++)  {
               // psi -= |j><j|S|psi>
               PrecCoeffG scp = ::dot(Sortho.colRef(j), psi);
               psi.plus_assign_ax (-scp, orthoStates.colRef(j));
            }
         } else {
             SxDiracVec<TPrecCoeffG> scps;
             scps = Sortho.adjoint () ^ psi;
             psi -= (orthoStates ^ scps);
         }
      }
   }
}

PsiG SxPAWOverlap::apply (const PsiG &psi) const
{
   SX_CHECK (pawPot);
   SX_CHECK (psi.handle);
   int ik = psi.handle->auxData.ik;
   SxDiracVec<Complex16> p = *pBasis | psi;
   
   SxDiracMat<Complex16> Sin(p.nRows (), p.nCols ());
   Sin.set (0.);
   int nOrb = int(pBasis->getNElements ());
   int offset = 0;
   for (int is = 0; is < pawPot->getNSpecies (); ++is)  {
      int nProjLocal = pawPot->getNProj (is);

      // get 1-center corrections
      const SxDiracMat<Double> &dS = pawPot->deltaS(is);

      // --- loop over atoms
      for (int ia = 0;
           offset < nOrb && pBasis->projectors->orbitalMap(offset).is == is;
           ++ia, offset+=nProjLocal)
      {

         // --- loops over projectors
         for (int ipl = 0; ipl < nProjLocal; ++ipl)  {
            const SxAOBasis::AoIndex &idI 
               = pBasis->projectors->refOrbMap(is)(ipl);
            for (int jpl = 0; jpl < nProjLocal; ++jpl)  {
               const SxAOBasis::AoIndex &idJ 
                  = pBasis->projectors->refOrbMap(is)(jpl);
               if (idI.l != idJ.l || idI.m != idJ.m) continue;

               // --- Ref. 1, Eq. 30, |pi>(...)<pj term
               // compute sum_j S(R)ij <pj|psi>
               for (int iState = 0; iState < p.nCols (); ++iState)  {
                  Sin(offset+ipl, iState) += dS(idI.n, idJ.n) 
                                           * p(offset+jpl,iState);
               }
            }
         }
      }
   }
   Sin.setBasis (pBasis.getPtr ());
   Sin.handle->auxData.ik = ik;
   const SxGBasis *gBasis
      = dynamic_cast<const SxGBasis *>(psi.getBasisPtr ());
   if (!gBasis)
      gBasis = psi.getBasis<SxPAWBasis> ().gBasis.getPtr ();
   SX_CHECK (gBasis);
   SxDiracMat<TPrecCoeffG> res = *gBasis | Sin;
   res += (*gBasis | psi);
   res.handle->auxData = psi.handle->auxData;
   res.setBasis (gBasis);
   return res;
}



