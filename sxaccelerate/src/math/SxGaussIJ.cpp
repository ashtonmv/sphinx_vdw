// ---------------------------------------------------------------------------
//
//           The general purpose cross platform C/C++ framework
//
//                       S x A c c e l e r a t e
//
//           Home:       https://www.sxlib.de
//           License:    Apache 2
//           Authors:    see src/AUTHORS
//
// ---------------------------------------------------------------------------

#include <SxGaussIJ.h>
#include <SxMathLib.h>

SxGaussIJ::SxGaussIJ ()
   : dist2(-1.), workOrder(-1)
{
   // empty
}

void SxGaussIJ::computeDoubleFactorial (int n)
{
   int nOld = (int)doubleFactorial.getSize () - 1;
   if (n > nOld) doubleFactorial.resize (n + 1, true);
   // --- precompute (2l+1)!! = 1 * 3 * 5 * ... (2l+1)
   if (nOld == -1) doubleFactorial(++nOld) = 1.0;
   for(int i = nOld + 1; i < doubleFactorial.getSize (); i++)
      doubleFactorial(i) = double(2*i + 1) * doubleFactorial(i-1);
}

void SxGaussIJ::computeK (double a2, int lmax)
{
   SX_CHECK (a2 > 1e-20, a2);
   // precompute (2l+1)!!
   computeDoubleFactorial (lmax);

   // --- set up K(n,l,dist2/rc2)
   double a = sqrt(a2);

   if (K.nCols () != lmax + 1 || K.nRows () != lmax + 1)
      K = SxMatrix<Double> (lmax+1, lmax+1);
#ifndef NDEBUG
   {
      double nan = sqrt(-1.);
      for (ssize_t i = 0; i < K.getSize (); ++i) K(i) = nan;
   }
#endif

   // --- set up for n == l and n == l+2
   {
      double prefac = SQRT_PI * a * exp (-a2);
      K(0,0) = 0.5 * PI * derf(a);
      for (int l = 0; l < lmax; ++l)  {
         K(l+1,l+1) = 0.5 * PI * doubleFactorial(l) * gammp(1.5 + l, a2);
         if (l + 2 <= lmax)
            K(l+2,l) = prefac * pow (2. * a2, double(l+1));
      }
   }

   // --- recursion for n = l + d
   for (int d = 4; d <= lmax; d+=2)  {
      for (int n = d; n <= lmax; ++n)  {
         int l = n-d;
         K(n,l) = (2 * l + 3) * K(n-1,l+1) - K(n,l+2);
      }
   }
}


void SxGaussIJ::set (const SxVector3<Double> &r, double rc2, int lSum)
{
   dist2 = r.normSqr ();
   computeDoubleFactorial (lSum);
   if (dist2 < 1e-12)  {
      if (K.nCols () != lSum + 1 || K.nRows () != lSum + 1)
         K = SxMatrix<Double> (lSum+1, lSum+1);
#ifndef NDEBUG
      {
         double nan = sqrt(-1.);
         for (ssize_t i = 0; i < K.getSize (); ++i) K(i) = nan;
      }
#endif
      // --- special case R=R'
      double dfac = 1.;
      // rc2 is presently rc1^2 + rc2^2 = 2 rc^2. Switch to rc2;
      rc2 *= 0.5;
      double prefac = 4. * sqrt(TWO_PI);
      for (int l = 0; l <= lSum; ++l)  {
         dfac *= (2 * l + 1);
         K(l,l) = prefac / (dfac * (2 * l + 1) * pow(rc2, l+0.5));
      }
      return;
   }

   // --- set up K(n,l,dist2/rc2)
   computeK (dist2 / rc2, lSum);

   // --- include r^(-n-1) into K
   for (int l = 0; l <= lSum; ++l)
      for (int n = l; n <= lSum; n+=2)
         K(n,l) *= pow (dist2, -0.5 * (n+1));

   // precompute Ylm(r)
   setupYlm (r, lSum); 
}

void SxGaussIJ::setDelta (const SxVector3<Double> &r,
                          double rc2A, double rc2B, int lSum)
{
   dist2 = r.normSqr ();
   computeDoubleFactorial (lSum);
   if (dist2 < 1e-12)  {
      if (K.nCols () != lSum + 1 || K.nRows () != lSum + 1)
         K = SxMatrix<Double> (lSum+1, lSum+1);
#ifndef NDEBUG
      {
         double nan = sqrt(-1.);
         for (ssize_t i = 0; i < K.getSize (); ++i) K(i) = nan;
      }
#endif
      // --- special case R=R'
      double dfac = 1.;
      // rc2 is presently rc1^2 + rc2^2 = 2 rc^2. Switch to rc2;
      rc2A *= 0.5;
      rc2B *= 0.5;
      double prefac = 4. * sqrt(TWO_PI);
      for (int l = 0; l <= lSum; ++l)  {
         dfac *= (2 * l + 1);
         K(l,l) = prefac / (dfac * (2 * l + 1)) 
                * (pow(rc2A, -l-0.5) - pow(rc2B, -l-0.5));
      }
      return;
   }

   /*
   if (K.getSize () == 0) K.reformat (lSum + 1, lSum + 1);
   SxMatrix<Double> K2 = K;
   // --- compute KB in workspace
   if (workspace.getSize () < sqr(lSum + 1)) 
      workspace = SxVector<Double> (sqr(lSum+1));
   K = SxMatrix<Double> ();
   K = workspace(SxIdx(0, sqr(lSum+1)-1));
   K.reshape (lSum+1, lSum+1);
   computeK (dist2 / rc2B, lSum);
   // --- compute KA in K
   K = SxMatrix<Double> ();
   K = K2;
   computeK (dist2 / rc2A, lSum);
   // subtract KB (in workspace) and scale with 1/r^{n+1}
   for (int l = 0; l <= lSum; ++l)  {
      for (int n = l; n <= lSum; n+=2)  {
         K(n,l) -= workspace(n + l * (lSum+1));
         K(n,l) *= pow (dist2, -0.5 * (n+1));
      }
   }
   */

   // --- set up K(n,l,dist2/rc2)
   if (K.nCols () != lSum + 1 || K.nRows () != lSum + 1)
      K = SxMatrix<Double> (lSum+1, lSum+1);
#ifndef NDEBUG
   {
      double nan = sqrt(-1.);
      for (ssize_t i = 0; i < K.getSize (); ++i) K(i) = nan;
   }
#endif

   // --- set up for n == l and n == l+2
   {  // + K(rc2A)
      double a2 = dist2/rc2A, a = sqrt(a2);
      double prefac = SQRT_PI * a * exp (-a2);
      K(0,0) = 0.5 * PI * derf(a);
      double pow2al = 2. * a2;
      for (int l = 0; l < lSum; ++l)  {
         K(l+1,l+1) = 0.5 * PI * doubleFactorial(l) * gammp(1.5 + l, a2);
         if (l + 2 <= lSum)  {
            //K(l+2,l) = prefac * pow (2. * a2, double(l+1));
            K(l+2,l) = prefac * pow2al;
         }
         pow2al *= 2. * a2;
      }
   }
   {  // - K(rc2B)
      double a2 = dist2/rc2B, a = sqrt(a2);
      double prefac = SQRT_PI * a * exp (-a2);
      K(0,0) -= 0.5 * PI * derf(a);
      double pow2al = 2. * a2;
      for (int l = 0; l < lSum; ++l)  {
         K(l+1,l+1) -= 0.5 * PI * doubleFactorial(l) * gammp(1.5 + l, a2);
         if (l + 2 <= lSum)  {
            // K(l+2,l) -= prefac * pow (2. * a2, double(l+1));
            K(l+2,l) -= prefac * pow2al;
         }
         pow2al *= 2. * a2;
      }
   }

   // --- recursion for n = l + d
   for (int d = 4; d <= lSum; d+=2)  {
      for (int n = d; n <= lSum; ++n)  {
         int l = n-d;
         K(n,l) = (2 * l + 3) * K(n-1,l+1) - K(n,l+2);
      }
   }

   // scale with 1/r^{n+1}
   double distInv = 1./sqrt(dist2), rl1 = distInv; 
   for (int l = 0; l <= lSum; ++l)  {
      //double rn1 = pow (dist2, -0.5 * (l+1));
      double rn1 = rl1;
      rl1 *= distInv;
      for (int n = l; n <= lSum; n+=2)  {
         //K(n,l) *= pow (dist2, -0.5 * (n+1));
         K(n,l) *= rn1;
         rn1 /= dist2;
      }
   }

   // precompute Ylm(r)
   setupYlm (r, lSum);
}

void SxGaussIJ::setupYlm (const SxVector3<Double> &r, int lMax)
{
   // --- precompute Ylm(r)
   if (ylmR.getSize () < sqr(lMax+1))
      ylmR = SxVector<Double> (sqr(lMax+1));
   SxVector<Double> ylmRef = ylmR(SxIdx(0,sqr(lMax+1)-1));
   SxYlm::getYlmArray (lMax, r(0), r(1), r(2), &ylmRef);
   for (int l = 0; l <= lMax; ++l)
      for (int m = -l; m <= l; ++m)
         ylmR(SxYlm::combineLm(l,m)) *= SxYlm::getYlmNormFactor(l,m);
}

SxMatrix<Double>
SxGaussIJ::compute (int lmax1, int lmax2, 
                    const SxYlm::SxClebschTable &clebschGordan,
                    int order,
                    SxMatrix<Double> *resPtr)
{
   SX_CHECK (K.nRows () >= lmax1 + lmax2 + order,
             K.nRows (), lmax1 + lmax2 + order);
   SX_CHECK (clebschGordan.getDim(0) >= sqr(lmax1+order+lmax2+1),
             clebschGordan.getDim(0), lmax1+lmax2+order);
   SX_CHECK (clebschGordan.getDim(1) >= sqr(lmax1+1),
             clebschGordan.getDim(1), lmax1+1);
   SX_CHECK (clebschGordan.getDim(2) >= sqr(lmax2+1),
             clebschGordan.getDim(2), lmax2);
   SX_CHECK ((dist2 >= 1e-12) || lmax1 == lmax2, lmax1, lmax2);
   lmax1 += order;
   workOrder = order;
   int nlm1 = sqr(lmax1 + 1);
   workNlm1 = nlm1;
   workNlm2 = sqr(lmax2 + 1);
   int nlm2 = workNlm2 * (order + 1);
   if (workspace.getSize () < nlm1 * nlm2)
      workspace.resize (nlm1 * nlm2);
   SxMatrix<Double> res;
   if (resPtr)  {
      res = *resPtr;
   } else {
      res = SxMatrix<Double> (workspace(SxIdx(0,nlm1 * nlm2 -1)));
      res.reshape (nlm1, nlm2);
   }
   SX_CHECK (res.nRows () == nlm1, res.nRows (), nlm1);
   SX_CHECK (res.nCols () == nlm2, res.nCols (), nlm2);

   if (dist2 < 1e-12)  {
      res.set (0.);
      if (order == 0)  {
         for (int l = 0, lm = 0; l <= lmax1; ++l)
            for (int m = -l; m <= l; ++m, ++lm)
               res(lm, lm) = K(l,l);
      }
      return res;
   }
   // --- now set up res
   for (int io = 0; io <= order; io++)  {
      for (int l2 = 0; l2 <= lmax2; ++l2)  {
         for (int m2 = -l2; m2 <= l2; ++m2)  {
            int lm2 = SxYlm::combineLm(l2,m2);
            for (int l1 = 0; l1 <= lmax1; ++l1)  {
               if (l1 + l2 + 2 * io > lmax1 + lmax2) continue;
               double prefac = 32. * PI;
               prefac /= doubleFactorial(l1) * doubleFactorial(l2);
               for (int m1 = -l1; m1 <= l1; ++m1)  {
                  int lm1 = SxYlm::combineLm(l1,m1);

                  double sumL = 0.;
                  for (int l = abs(l1-l2); l <= l1+l2; l+=2)  {
                     int lm = SxYlm::combineLm (l,-l);
                     double sumM = 0.;
                     for (int m = -l; m <= l; ++m, ++lm)
                        sumM += clebschGordan(lm, lm1, lm2) * ylmR(lm);
                     sumL += K(l1+l2+ 2 * io,l) * sumM;
                  }
                  SX_CHECK_NUM (sumL);
                  res(lm1, lm2 + io*workNlm2) = prefac * sumL;
               }
            }
         }
      }
   }
   if (!resPtr)  {
      res = SxMatrix<Double> ();
      res = SxMatrix<Double> (workspace(SxIdx(0,workNlm1 * workNlm2 -1)));
      res.reshape (workNlm1, workNlm2);
   }
   return res;
}

SxVector3<Double>
SxGaussIJ::getForce (int lm1, int lm2, 
                     const SxYlm::SxClebschTable &clebschGordan)
{
   SX_CHECK (lm1 < 1000000 && lm1 >= 0, lm1);
   SX_CHECK (lm1 < workNlm1, lm1, workNlm1);
   SX_CHECK (workOrder >= 1, workOrder);
   // --- get l1 and m1
   int l1, m1;
   for (l1 = 0; l1 < 1000; ++l1)  {
      m1 = lm1 - SxYlm::combineLm (l1,0);
      if (m1 <= l1) break;
   }
   SX_CHECK (clebschGordan.getDim(0) >= sqr(l1+1),
             clebschGordan.getDim(0), l1+1);
   SX_CHECK  (clebschGordan.getDim(1) >= 4, clebschGordan.getDim(1));
   SX_CHECK (clebschGordan.getDim(2) >= sqr(l1+2),
             clebschGordan.getDim(2), l1);
   SxVector3<Double> force(0., 0., 0.);
   for (int LM = 1; LM < 4; ++LM)  {
      int xyz = LM % 3; // 2 -> y ; 3 -> z ; 4 -> x
      for (int l3 = abs(l1 - 1); l3 <= l1 + 1; l3 += 2)  {
         int io = (l1 - l3 + 1) / 2;
         for (int m3 = -l3; m3 <= l3; ++m3)  {
            int lm3 = SxYlm::combineLm(l3,m3);
            force(xyz) += clebschGordan(lm1, LM, lm3) * doubleFactorial(l3)
                          * workspace (lm3 + (io * workNlm2 + lm2) * workNlm1);
         }
      }
   }
   force *= sqrt(FOUR_PI / 3.) / doubleFactorial(l1);
   return force;
}

SxMatrix3<Double>
SxGaussIJ::getHesse (int lm1, int lm2, 
                     const SxYlm::SxClebschTable &clebschGordan)
{
   SX_CHECK (lm1 < 1000000 && lm1 >= 0, lm1);
   SX_CHECK (lm1 < workNlm1, lm1, workNlm1);
   SX_CHECK (workOrder >= 2, workOrder);
   // --- get l1 and m1
   int l1, m1;
   for (l1 = 0; l1 < 1000; ++l1)  {
      m1 = lm1 - SxYlm::combineLm (l1,0);
      if (m1 <= l1) break;
   }
   SX_CHECK (clebschGordan.getDim(0) >= sqr(l1+1),
             clebschGordan.getDim(0), l1+1);
   SX_CHECK  (clebschGordan.getDim(1) >= 9, clebschGordan.getDim(1));
   SX_CHECK (clebschGordan.getDim(2) >= sqr(l1+3),
             clebschGordan.getDim(2), l1);
   SxMatrix3<Double> H;
   H.set (0.);
   const double sqrt3 = sqrt(3.), sqrt5 = sqrt(5.);
   for (int LM = 4; LM < 9; ++LM)  {
      double res = 0.;
      for (int l3 = abs(l1 - 2); l3 <= l1 + 2; l3 += 2)  {
         int io = (l1 - l3 + 2) / 2;
         for (int m3 = -l3; m3 <= l3; ++m3)  {
            int lm3 = SxYlm::combineLm(l3,m3);
            res += clebschGordan(lm1, LM, lm3) * doubleFactorial(l3)
                 * workspace (lm3 + (io*workNlm2 + lm2) * workNlm1);
         }
      }
      SX_CHECK_NUM(res);
      res /= sqrt5; // prefactor 1/3 done below
      switch (LM) {
         case 4: H(1,0) = H(0,1) = sqrt3 * res; break;
         case 5: H(1,2) = H(2,1) = sqrt3 * res; break;
         case 6: H(0,0) -= res;
                 H(1,1) -= res;
                 H(2,2) += 2. * res;
                 break;
         case 7: H(0,2) = H(2,0) = sqrt3 * res; break;
         case 8: H(0,0) += sqrt3 * res;
                 H(1,1) -= sqrt3 * res;
      }
   }
   double H0 = -workspace (lm1 + (1*workNlm2 + lm2) * workNlm1);
   H0 *= doubleFactorial(l1);
   H0 /= sqrt(FOUR_PI); // = clebschGordan (lm1, 0, lm3=lm1)
   for (int i = 0; i < 3; ++i) H(i,i) += H0;
   H *= sqrt(FOUR_PI) / (3. * doubleFactorial(l1));
   return H;
}
