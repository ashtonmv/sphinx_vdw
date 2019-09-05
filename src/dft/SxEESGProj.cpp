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

#include <SxEESGProj.h>
#include <SxGBasis.h>
#include <SxProjector.h>

// --- ref1: H.-S. Lee, M. E. Tuckerman, G. J. Martyna
//           Efficient Evaluation of Nonlocal Pseudopotentials via
//           Euler Exponential spline interpolation
//           ChemPhysChem 2005, 6, 1827-1835

SxEESGProj::SxEESGProj (const SxCell &cell, const SxVector3<Int> &mesh,
                        int order)
: splineOrder (order),
  cachedG(NULL)
{
   SX_CHECK  (order > 1, order);
   SX_CHECK (order < mesh(0), order, mesh(0));
   SX_CHECK (order < mesh(1), order, mesh(1));
   SX_CHECK (order < mesh(2), order, mesh(2));
   sBasisPtr = SxPtr<SxRBasis>::create (mesh(0), mesh(1), mesh(2), cell);
   SxFFT1d fftx(SxFFT::Forward, mesh(0));
   SxFFT1d ffty(SxFFT::Forward, mesh(1));
   SxFFT1d fftz(SxFFT::Forward, mesh(2));
   dpx = getDp (splineOrder, fftx);
   dpy = getDp (splineOrder, ffty);
   dpz = getDp (splineOrder, fftz);
}

void SxEESGProj::cardBSplineAll (double xx, int p, SxVector<Double> &res)
{
   SX_CLOCK (CardBSpline);
   SX_CHECK (res.getSize () >= p, res.getSize (), p);
   if (p == 1) {
      // needed for derivative of p=2
      res(0) = 1.;
      return;
   }
   SX_CHECK (p >= 2, p);
   SX_CHECK (xx >= 0. && xx <= 1., xx);
   // Start with p = 2
   res(0) = xx;
   res(1) = 1. - xx;

   // iterate upward
   for (int pp = 2; pp < p; ++pp)  {
      res(pp) = 0.;
      double inv = 1. / pp;
      double invX = inv * xx + 1.;
      double invPp = 1. + inv;
      for (int i = pp; i > 0; --i, invX -= inv)  {
         // compute cardBSpline(xx+i,pp+1)
         // ref 1, Eq. 12
         res(i) = invX         * res(i)    // cardBSpline(xx+i, pp)
                + (invPp-invX) * res(i-1); // cardBSpline(xx+i-1, pp)
      }
      res(0) *= invX;
   }
}

double SxEESGProj::cardBSpline (double x, int p, SxVector<Double> &work)
{
   SX_CHECK (work.getSize () >= p, work.getSize (), p);
   if (x <= 0. || x >= p) return 0.;
   cardBSplineAll (x - floor(x), p, work);
   SX_CHECK_NUM (work((int)floor(x)));
   return work((int)floor(x));
}

SxVector<Complex16> SxEESGProj::getDp (int p, SxFFT1d &fft)
{
   SX_CLOCK (DpSetup);
   // This function evaluates Eq. 11 of ref. 1
   SX_CHECK (p > 1, p);
   int s, n = fft.meshSize;
   SX_CHECK (n > p, n, p);
   // Equation (11): summand
   SxVector<Double> mp(p);
   cardBSplineAll (0., p, mp);
   for (s = 0; s < p; ++s) fft.setElement(s, mp(s));
   for (     ; s < n; ++s) fft.setElement(s, 0.);
   // sum in [] of Eq. 11
   fft.fftForward ();
   SxVector<Complex16> res(n);
   // Equation (11)
   double phaseFac = TWO_PI * p / double(n);
   for (int g = 0; g < n; ++g)  {
      res(g) = getPhase (g * phaseFac) /
         reinterpret_cast<const SxComplex16 &>(fft.outArray[g]);
   }
   if (n % 2 == 0) res(n/2) = 0.; // |G| < n/2
   return res;
}

SxVector<Complex16> 
SxEESGProj::projectGrad (const SxDiracVec<Complex16> &refPhi,
                         const SxDiracVec<Complex16> &psi,
                         const SxAtomicStructure     &centers) const
{
   SX_CLOCK (Projection);
   const SxGBasis *gPtr = dynamic_cast<const SxGBasis*>(psi.getBasisPtr ());
   SX_CHECK (gPtr);
   const SxRBasis &S = *sBasisPtr;
   // make sure that S is registered at G
   gPtr->registerBasis (S);

   SxDiracVec<Complex16> pG(*gPtr), pS;
   {
      // --- set up pG(ig) = psi(ig).conj () * phi(ig) * prod_dir d_p(dir,ig)
      // Equation 15b of ref. 1
      SX_CLOCK (PhiPsiDp);
      SxDiracVec<Complex16>::Iterator phiIt = refPhi.begin (),
                                      psiIt = psi.begin (),
                                      pGIt  = pG.begin ();
      if (cachedG != gPtr) updateG (gPtr);
      ssize_t iPacked = 0;
      for (int ig = 0; ig < gPtr->ng; ++ig, iPacked+=3)  {
         // Equation 15b
         *pGIt++ = (*phiIt++).conj () * (*psiIt++)
                 * dpx(packedGrel(iPacked    ))
                 * dpy(packedGrel(iPacked + 1))
                 * dpz(packedGrel(iPacked + 2));
      }
   }

   // Equation 15a
   {
      SX_CLOCK (GtoS);
      pS = ( S | pG );
   }

   SxMatrix<Complex16> res(4, centers.getNAtoms ());
   SxVector<Double> mpx(splineOrder), mpy(splineOrder), mpz(splineOrder),
                    work(splineOrder);
   SxVector<Double> mpxf(splineOrder), mpyf(splineOrder), mpzf(splineOrder);
   for (int ia = 0; ia < centers.getNAtoms (); ++ia)  {
      Coord pos = S.cell.carToRel(centers.getAtom(ia)) * S.fft3d.mesh;
      // --- set up phase factors
      // "phase" factors for Eq. 14
      {
         SX_CLOCK (ProjMp);
         int i;
         SxVector3<Double> posOffset = pos - ceil(pos) + 1.;
         cardBSplineAll(posOffset(0), splineOrder, work);
         for (i = 0; i < splineOrder; ++i) mpx(i) = work(splineOrder-1-i);
         cardBSplineAll(posOffset(1), splineOrder, work);
         for (i = 0; i < splineOrder; ++i) mpy(i) = work(splineOrder-1-i);
         cardBSplineAll(posOffset(2), splineOrder, work);
         for (i = 0; i < splineOrder; ++i) mpz(i) = work(splineOrder-1-i);
         // --- mp(u) gradients with respect to u
         // Eq. (13) for u = splineOrder-1-i
         work(splineOrder-1) = 0.;
         cardBSplineAll(posOffset(0), splineOrder - 1, work);
         for (i = 0; i < splineOrder-1; ++i) 
            mpxf(i) = work(splineOrder-1-i) - work(splineOrder-2-i);
         mpxf(splineOrder-1) = work(0);
         
         cardBSplineAll(posOffset(1), splineOrder - 1, work);
         for (i = 0; i < splineOrder-1; ++i) 
            mpyf(i) = work(splineOrder-1-i) - work(splineOrder-2-i);
         mpyf(splineOrder-1) = work(0);
         
         cardBSplineAll(posOffset(2), splineOrder - 1, work);
         for (i = 0; i < splineOrder-1; ++i) 
            mpzf(i) = work(splineOrder-1-i) - work(splineOrder-2-i);
         mpzf(splineOrder-1) = work(0);
      }


      // --- find boundaries
      SxVector3<Int> from = (SxVector3<Int>(ceil(pos)) + (-splineOrder))
                             % S.fft3d.mesh;
      for (int d = 0; d < 3; ++d) if (from(d) < 0) from(d) += S.fft3d.mesh(d);

      // --- perform sum of Eq. 14
      {
         SX_CLOCK (ProjSumS);
         SxComplex16 sum = 0., fx = 0., fy = 0., fz = 0., sumz, sumzf;
         SxDiracVec<Complex16>::Iterator pIt;
         SxVector<Double>::Iterator mpIt, mpfIt;
         int zUp = min(S.fft3d.mesh(2)-from(2), splineOrder);
         int z;
         for (int x = 0; x < splineOrder; ++x)  {
            for (int y = 0; y < splineOrder; ++y)  {
               sumz = 0.;
               sumzf = 0.;
               pIt = pS.begin ();
               pIt += S.fft3d.mesh.getMeshIdx(from(0)+x,from(1)+y,from(2),
                                              SxMesh3D::Unknown);
               mpIt = mpz.begin ();
               mpfIt = mpzf.begin ();
               for (z = 0; z < zUp; ++z, ++pIt)  {
                  sumz  += *pIt * *mpIt++;
                  sumzf += *pIt * *mpfIt++;
               }
               pIt += -S.fft3d.mesh(2);
               for (; z < splineOrder; ++z, ++pIt)  {
                  sumz  += *pIt * *mpIt++;
                  sumzf += *pIt * *mpfIt++;
               }
               sum += mpx (x) * mpy (y) * sumz;
               fx  += mpxf(x) * mpy (y) * sumz;
               fy  += mpx (x) * mpyf(y) * sumz;
               fz  += mpx (x) * mpy (y) * sumzf;
            }
         }
         res(0,ia) = sqrt(S.cell.volume) * sum;
         {
            // backtransform force element from relative to absolute
            // coordinates
            // note: cartesian vec x, car->rel transformation mat A
            // d E = x^T dE/dx = (A x)^T (dE/ d(A x))
            // => dE / dx = A^T dE/d(Ax)
            // NOTE:  A = cell.inv is taken from SxCell!
            SxVector3<Complex16> fAbs 
               = (S.cell.inverse ().transpose ()
                  ^ (SxVector3<Complex16>(fx, fy, fz) * S.fft3d.mesh));
            // prefactor
            fAbs *= sqrt(S.cell.volume);
            res(1,ia) = fAbs(0);
            res(2,ia) = fAbs(1);
            res(3,ia) = fAbs(2);
         }

      }
   }
   VALIDATE_VECTOR(res);
   return res;

}

SxVector<Complex16>
SxEESGProj::project (const SxDiracVec<Complex16> &refPhi,
                     const SxDiracVec<Complex16> &psi,
                     const SxAtomicStructure     &centers) const
{
   SX_CLOCK (Projection);
   const SxGBasis *gPtr = dynamic_cast<const SxGBasis*>(psi.getBasisPtr ());
   SX_CHECK (gPtr);
   const SxRBasis &S = *sBasisPtr;
   // make sure that S is registered at G
   gPtr->registerBasis (S);

   SxDiracVec<Complex16> pG(*gPtr), pS;
   {
      // --- set up pG(ig) = psi(ig).conj () * phi(ig) * prod_dir d_p(dir,ig)
      // Equation 15b of ref. 1
      SX_CLOCK (PhiPsiDp);
      SxDiracVec<Complex16>::Iterator phiIt = refPhi.begin (),
                                      psiIt = psi.begin (),
                                      pGIt  = pG.begin ();
      if (cachedG != gPtr) updateG (gPtr);
      ssize_t iPacked = 0;
      for (int ig = 0; ig < gPtr->ng; ++ig, iPacked+=3)  {
         // Equation 15b
         *pGIt++ = (*phiIt++).conj () * (*psiIt++)
                 * dpx(packedGrel(iPacked    ))
                 * dpy(packedGrel(iPacked + 1))
                 * dpz(packedGrel(iPacked + 2));
      }
   }

   // Equation 15a
   {
      SX_CLOCK (GtoS);
      pS = ( S | pG );
   }

   SxVector<Complex16> res(centers.getNAtoms ());
   SxVector<Double> mpx(splineOrder), mpy(splineOrder), mpz(splineOrder),
                    work(splineOrder);
   for (int ia = 0; ia < centers.getNAtoms (); ++ia)  {
      Coord pos = S.cell.carToRel(centers.getAtom(ia)) * S.fft3d.mesh;
      // --- set up phase factors
      // "phase" factors for Eq. 14
      {
         SX_CLOCK (ProjMp);
         int i;
         SxVector3<Double> posOffset = pos - ceil(pos) + 1.;
         cardBSplineAll(posOffset(0), splineOrder, work);
         for (i = 0; i < splineOrder; ++i) mpx(i) = work(splineOrder-1-i);
         cardBSplineAll(posOffset(1), splineOrder, work);
         for (i = 0; i < splineOrder; ++i) mpy(i) = work(splineOrder-1-i);
         cardBSplineAll(posOffset(2), splineOrder, work);
         for (i = 0; i < splineOrder; ++i) mpz(i) = work(splineOrder-1-i);
      }


      // --- find boundaries
      SxVector3<Int> from = (SxVector3<Int>(ceil(pos)) + (-splineOrder))
                             % S.fft3d.mesh;
      for (int d = 0; d < 3; ++d) if (from(d) < 0) from(d) += S.fft3d.mesh(d);

      // --- perform sum of Eq. 14
      {
         SX_CLOCK (ProjSumS);
         SxComplex16 sum = 0., sumz;
         SxDiracVec<Complex16>::Iterator pIt;
         SxVector<Double>::Iterator mpIt;
         int zUp = min(S.fft3d.mesh(2)-from(2), splineOrder);
         int z;
         for (int x = 0; x < splineOrder; ++x)  {
            for (int y = 0; y < splineOrder; ++y)  {
               sumz = 0.;
               pIt = pS.begin ();
               pIt += S.fft3d.mesh.getMeshIdx(from(0)+x,from(1)+y,from(2),
                                              SxMesh3D::Unknown);
               mpIt = mpz.begin ();
               for (z = 0; z < zUp; ++z)
                  sumz += *pIt++ * *mpIt++;
               pIt += -S.fft3d.mesh(2);
               for (; z < splineOrder; ++z)
                  sumz += *pIt++ * *mpIt++;
               sum += mpx(x) * mpy(y) * sumz;
            }
         }
         res(ia) = sqrt(S.cell.volume) * sum;
      }
   }
   VALIDATE_VECTOR(res);
   return res;

}

SxDiracVec<Complex16>
SxEESGProj::gradient (const SxVector<Complex16>   &proj,
                      const SxDiracVec<Complex16> &refPhi,
                      const SxAtomicStructure     &centers) const
{
   SX_CLOCK (Gradient);
   const SxGBasis *gPtr = dynamic_cast<const SxGBasis*>(refPhi.getBasisPtr ());
   SX_CHECK (gPtr);
   const SxRBasis &S = *sBasisPtr;

   SxDiracVec<Complex16> pS(S), res;
   SxVector<Double> mpx(splineOrder), mpy(splineOrder), mpz(splineOrder),
                    work(splineOrder);

   pS.set (0.);
   for (int ia = 0; ia < centers.getNAtoms (); ++ia)  {
      Coord pos = S.cell.carToRel(centers.getAtom(ia)) * S.fft3d.mesh;
      // --- set up phase factors
      // "phase" factors for Eq. 17b
      {
         SX_CLOCK (ProjMp);
         int i;
         SxVector3<Double> posOffset = pos - ceil(pos) + 1.;
         cardBSplineAll(posOffset(0), splineOrder, work);
         for (i = 0; i < splineOrder; ++i) mpx(i) = work(splineOrder-1-i);
         cardBSplineAll(posOffset(1), splineOrder, work);
         for (i = 0; i < splineOrder; ++i) mpy(i) = work(splineOrder-1-i);
         cardBSplineAll(posOffset(2), splineOrder, work);
         for (i = 0; i < splineOrder; ++i) mpz(i) = work(splineOrder-1-i);
      }


      // --- find boundaries
      SxVector3<Int> from = (SxVector3<Int>(ceil(pos)) + (-splineOrder))
                             % S.fft3d.mesh;
      for (int d = 0; d < 3; ++d) if (from(d) < 0) from(d) += S.fft3d.mesh(d);

      // --- perform sum of Eq. 17b
      {
         SX_CLOCK (GradientQ);
         SxComplex16 factor, factorX;
         SxDiracVec<Complex16>::Iterator pIt;
         SxVector<Double>::Iterator mpIt;
         int zUp = min(S.fft3d.mesh(2)-from(2), splineOrder);
         int z;
         for (int x = 0; x < splineOrder; ++x)  {
            factorX = S.fft3d.meshSize / sqrt(S.cell.volume)
                      * mpx(x) * proj(ia);
            for (int y = 0; y < splineOrder; ++y)  {
               factor =  factorX * mpy(y);
               pIt = pS.begin ();
               pIt += S.fft3d.mesh.getMeshIdx(from(0)+x,from(1)+y,from(2),
                                              SxMesh3D::Unknown);
               mpIt = mpz.begin ();
               for (z = 0; z < zUp; ++z)
                  *pIt++ += *mpIt++ * factor;
               pIt += -S.fft3d.mesh(2);
               for (; z < splineOrder; ++z)
                  *pIt++ += *mpIt++ * factor;
            }
         }
      }
   }

   // Equation 17a
   {
      SX_CLOCK (StoG);
      res = ( *gPtr | pS );
   }


   {
      // Equation 16  of ref. 1
      SX_CLOCK (GradientMul);
      if (cachedG != gPtr) updateG (gPtr);
      SxDiracVec<Complex16>::Iterator phiIt = refPhi.begin (),
                                      resIt  = res.begin ();
      ssize_t iPacked = 0;
      for (int ig = 0; ig < gPtr->ng; ++ig, iPacked+=3)  {
         *resIt++ *= (  dpx(packedGrel(iPacked    ))
                      * dpy(packedGrel(iPacked + 1))
                      * dpz(packedGrel(iPacked + 2))).conj ()
                     * (*phiIt++);
      }
   }

   VALIDATE_VECTOR(res);
   return res;

}

void SxEESGProj::updateG (const SxGBasis *gPtr) const
{
   SX_CHECK (gPtr);
   int sBasisIdx = gPtr->getBasisId (sBasisPtr.getPtr ());
   // --- maps ig vectors to 1D phases
   packedGrel.resize (3*gPtr->ng);
   SxVector3<Int> Grel;
   SxVector<TPrecFFTIdx>::Iterator n123It = gPtr->n123(sBasisIdx).begin ();
   // x,y,z are just names. We work in relative coordinates here.
   for (int ig = 0; ig < gPtr->ng; ++ig, ++n123It)  {
      Grel = sBasisPtr->fft3d.mesh.getMeshVec(*n123It, SxMesh3D::Positive);
      packedGrel(3 * ig    ) = Grel(0);
      packedGrel(3 * ig + 1) = Grel(1);
      packedGrel(3 * ig + 2) = Grel(2);
   }
   cachedG = gPtr;
}
