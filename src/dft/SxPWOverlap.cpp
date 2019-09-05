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

#include <SxPWOverlap.h>
#include <SxTimer.h>
#include <SxPW.h>

double SxPWOverlap::normSqr (const SxDiracVec<Complex16> &psi) const
{
   return psi.normSqr ();
}

SxComplex16 SxPWOverlap::dot (const SxDiracVec<Complex16> &x,
                              const SxDiracVec<Complex16> &y) const
{
   return ::dot (x, y);
}

/// Apply overlap operator
SxDiracVec<Complex16> 
SxPWOverlap::apply (const SxDiracVec<Complex16> &psi) const
{
   // nothing to be done
   return psi;
}

/// Set states x orthogonal to states y
void SxPWOverlap::setOrthogonal (PsiG *X, const PsiG &orthoStates) const
{
   // by C. Freysoldt
   // This routine is heavily commented on technical details
   // This is because I'm implementing it now, but the XPress library
   // will need changes here, but maybe not by me
   SX_CLOCK (Timer::ortho);

   SX_CHECK (X);
   if (orthoStates.getSize () == 0) return;

   SX_CHECK (X->nRows () == orthoStates.nRows (),
             X->nRows (), orthoStates.nRows ());

   SxDiracVec<TPrecCoeffG> &psi = *X;
   int nPsi = (int)psi.nCols ();
   if (nPsi == 0) nPsi = 1;
   int nOrtho = (int)orthoStates.nCols ();

   int j;
   SxDiracVec<TPrecCoeffG> psiJ;
   PrecCoeffG scp;
   if (nOrtho > 150 || nPsi > 1)  {
      // for low n, the direct (convential) BLAS1 approach is faster
      // because there is a certain overhead, and BLAS2 routines may
      // not be fastet for BLAS1-like tasks
      // the cross-over point should be tested
      // this should be the zdotc / gemv(transpose='c') crossover for nPsi == 1
      // this should be the gemv / gemm crossover for nPsi > 1

      // --- get <j|psi> for all 0 <= j < nOrtho
      SxDiracVec<TPrecCoeffG> scps;
      if (nPsi > nOrtho) 
         // --- variant 1 (proper for XPress)
         scps = orthoStates.adjoint () ^ psi;
      else 
         // --- variant 2
         scps = (psi.adjoint () ^ orthoStates).adjoint ();
      
      // --- subtract |j><j|psi>
      // BLAS2 (nPsi == 1) or BLAS3 version
      psi -= (orthoStates ^ scps);
      
   } else {
      // --- conventional BLAS1 approach
      for (j=0; j < nOrtho; j++)  {
         psiJ = orthoStates.colRef(j);
//       psi -= psiJ * scp;
         psi.plus_assign_ax (-::dot(psiJ, psi), psiJ); // psi -= |j><j|psi>
      }
   }
}

/// Orthonormalize states x
void SxPWOverlap::orthonormalize (PsiGI *psiPtr, 
                                  SxOrthoMethod how) const
{
   SX_CHECK (psiPtr);
   PsiGI &psi = *psiPtr;
   int nStates = (int)psi.nCols ();

   if (how == GramSchmidt)  {
      // --- Gram/Schmidt orthogonalization
      SX_CLOCK (Timer::ortho);
      SxDiracMat<TPrecCoeffG> S = psi.overlap (psi);
      SxDiracMat<TPrecCoeffG> U = S.choleskyDecomposition ().adjoint ();
      psi.rotate (U.inverse ());
   }  else  {
   // --- diagonalization, according to LOEWDIN orthogonalization scheme
      SX_CLOCK (Timer::ortho);
      PsiGI psiGI;
      int i, j;

      SxDiracMat<TPrecCoeffG> U, Ieps;      // S being the overlap matrix,
                                               // U is U^-1/2, and
      SxDiracSymMat<TPrecCoeffG> S;
      SxDiracSymMat<TPrecCoeffG>::Eigensystem eig;

      S = SxDiracSymMat<TPrecCoeffG> (nStates);
      {
         SX_CLOCK (Timer::loewdinS);
         SxDiracMat<TPrecCoeffG> fullS = psi.overlap (psi);
         for (j=0; j < nStates; j++)  {
            for (i=0; i <= j; i++)  {
               S(i,j)  = fullS(i,j);
            }
         }
      }

      eig   = S.eigensystem();
      Ieps  = Ieps.identity(1. / sqrt(eig.vals) );
      U     = eig.vecs ^ Ieps ^ eig.vecs.adjoint();

      psi.rotate (U);

   }
}
