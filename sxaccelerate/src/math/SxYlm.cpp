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

#include <SxYlm.h>

inline int sxmax(int a, int b, int c)
{
   int ab = (a > b) ? a : b;
   return (ab > c) ? ab : c;
}

inline int sxmin(int a, int b, int c)
{
   int ab = (a < b) ? a : b;
   return (ab < c) ? ab : c;
}

// ---------------------------------------------------------------------------
inline double SxYlm::jsbTaylor (int l, double x)
{
   double result = 0.0;
   //int n = 0;
   if (l==0)   {
      double twoN = 0., x2 = x * x, x4 = x2 * x2, denom = 6.;
      double px = 1./6.;
      // px is x^{4n} / (2n + 3)!
      // denom = (2n + 2)(2n + 3)
      // twoN = 2n
      // result becomes
      // sum_n (-1)^n x^2n / (2n+1)!
      while (px > 1e-14)   {
         result += px * (denom - x2);
         twoN += 4.;
         denom = (twoN + 2.) * (twoN+3.);
         px *= x4 / (twoN * (twoN+1.) * denom);
      }
      return result;
   }
   double delta = 1.0;
   double twoN = 0.;
   double x2 = x*x;
   if (l==1)   {
      double px = x;
      while (fabs(delta) > 1e-14)   {
         delta = px / ((twoN+3.) * (twoN+1.));
         result += delta;
         twoN += 2.;
         //for (int i = 1; i <= 2*n;i++) px *= x / double(i);
         px *= -x2 / ((twoN-1.) * (twoN));
      }
   } else if (l==2) {
      double px = x2;
      while (fabs(delta) > 1e-14)   {
         delta = px / ((twoN+5.) * (twoN+3.) * (twoN+1.));
         result += delta;
         twoN += 2.;
         //for (int i = 1; i <= 2*n;i++) px *= x / double(i);
         px *= -x2 / ((twoN-1.) * (twoN));
      }
   } else {
      SX_EXIT;
   }
   return result;
}
double SxYlm::jsb(int l, double x)  
{
   SX_CHECK (l >= 0, l);
   SX_CHECK (x > -1e-12, x);

   if (l <= 2 && x < 1)  {
      return jsbTaylor(l,x);
   }
   
   if (x < 1e-10) {
      double xl = x, dfl = 1., f = 1.;
      for (int i = 1; i < l; ++i, f+= 2.)  {
         xl *= x; // x^l in the end
         dfl *= f; // (2l-1)!! in the end
      }
      return xl / dfl;
   }

   // --- special cases
   if (l == 0) return sin(x) / x;
   
   double invx = 1./x;
   double s,c;
   sincos(x, &s, &c);

   if (l == 1) return (s*invx - c)*invx;
   if (l == 2) {
      double tix = 3. * invx;
      return ((tix*invx - 1.) * s - tix * c) * invx;
   }

   if (x > 2 * l - 1.) {
      // upward recursion
      double jn2 = s * invx;       // l = 0, j_(n-2)
      double jn1 = (jn2 - c) * invx; // l = 1, j_(n-1)
      double jn = 0., twoNminusOne = 3.;
      for (int n = 2; n <= l; n++)  {
         jn = twoNminusOne * invx * jn1 - jn2;
         jn2 = jn1;
         jn1 = jn;
         twoNminusOne += 2.;
      }
      return jn;
   } else {
      // downward recursion
      double jn2 = 0., jn1 = 1., jn = 0.;
      int n = 2 * l - 1 + int(15. * sqrt(double(l)));
      double twoNplusThree = 2 * n + 3;
      double res = 0.;
      for (; n >= 0; --n)  {
         jn = twoNplusThree *invx * jn1 - jn2;
         if (n == l) res = jn;
         if (n > l)  {
            // renormalize jn1 to 1 to avoid overflow
            jn2 = 1.; // jn1 / jn1
            jn1 = jn / jn1;
         } else {
            // rescale by x (compensated by pow below)
            jn2 = jn1 * x;
            jn1 = jn * x;
         }
         twoNplusThree -= 2.;
      }
      return res * (jsb(0,x) / jn) * pow (x, double(l));
   }
}

SxArray2<double> SxYlm::getJsbDerivCoeffs(int L, int maxOrder)
{
   SxArray2<double> FkL(L + maxOrder + 1, maxOrder + 1);
   FkL.set (0.);
   FkL(L,0) = 1.;
   for (int m = 0; m < maxOrder; ++m)  {
      // k = 0
      double invM = 1./double(m+1);
      FkL(1,m+1) -= FkL(0,m) * invM;
      for (int k = max(1,L-m); k <= L + m; ++k)  {
         FkL(k-1,m+1) += double(k)/double(2*k+1) * FkL(k,m) * invM;
         FkL(k+1,m+1) -= double(k+1)/double(2*k+1) * FkL(k,m) * invM;
      }
   }
   return FkL;
}

void SxYlm::getPlmArray(int lmax, int m, double x, SxVector<Double> *out)
{
   SX_CHECK (out);
   SX_CHECK (lmax >= 0, lmax);
   SX_CHECK (m >= 0 && m <= lmax, m, lmax);
   SX_CHECK (out->getSize () == lmax + 1, out->getSize (), lmax); 
   SX_CHECK (fabs(x) <= 1., x);

   // This is an extremely low-level function, so we do some
   // low-level stuff here.
   // get pointer to m-th element
   double *P = & ((*out)(m));
   
   
   // --- Calculate start value 
   // P(l=m,m) = (2m - 1)!! * sqrt(1-x^2) ^ m 
   {
      double Plm = 1.;
      if (m)  {
         // (1-x)*(1+x) is more precise than 1 - x*x
         double factor = sqrt((1. - x) * (1. + x));
         double dfac = 2. * factor;
         for (int i = 1; i < m; i++)  {
            Plm *= factor;
            factor += dfac;
         }
         Plm *= factor; // i = m
      }
      if (m & 1)
         *P = -Plm;
      else
         *P = Plm;
   }

   // --- now recursion on l
   if (m < lmax)  {
      double h2 = 2 * m + 1; // l + m - 1
      // P(l=m+1,m) = x * (2m + 1) * P(l=m,m)
      P[1] = P[0] * x * h2;
      if (m + 1 < lmax)  {
         double h3 = 2.; // l - m
         double h1 = x * (h2 + h3); // x * (2l - 1)
         for (int l = m+2; l < lmax; ++l, ++P)  {
            // P(l,m) = (x*(2l-1)*P(l-1,m) - (l+m-1)*P(l-2,m))/(l-m)
            P[2] = (h1 * P[1] - h2 * P[0]) / h3;
            h1 += 2. * x;
            h2 += 1.;
            h3 += 1.;
         }
         // l = lmax here so we don't need to update h1, h2 and h3
         P[2] = (h1 * P[1] - h2 * P[0]) / h3;
      } // lmax > m+1
   } // lmax > m
}

inline void getPlmArray(int lmax, double cos_theta, SxVector<Double> *out)
{
   SX_CHECK (out);
   SX_CHECK (lmax >= 0, lmax);
   SX_CHECK (out->getSize () == sqr(lmax + 1), out->getSize (), lmax); 
   SX_CHECK (fabs(cos_theta) <= 1., cos_theta);

   SxVector<Double> P = *out;
   
   // --- Calculate start value 
   // P(l=m,m) = (-1)^m (2m - 1)!! * sqrt(1-x^2) ^ m 
   {
      double Plm = 1.;
      P(0) = Plm;
      if (lmax)  {
         // (1-x)*(1+x) is more precise than 1 - x*x
         double factor = sqrt((1. - cos_theta) * (1. + cos_theta));
         double dfac = 2. * factor;
         for (int m = 1; m < lmax; m++)  {
            Plm *= -factor;
            P(m*(m+2)) = Plm;
            factor += dfac;
         }
         Plm *= -factor; // i = m
         P(lmax*(lmax+2)) = Plm;
      }
   }

   // --- now recursion on l
   for (int m = 0; m < lmax; ++m)  {
      double P0, P1, P2;
      P0 = P(m*(m+2));
      double h2 = 2 * m + 1; // l + m - 1
      // P(l=m+1,m) = x * (2m + 1) * P(l=m,m)
      P((m+4)*m+2) = P1 = P0 * cos_theta * h2;
      if (m + 1 < lmax)  {
         double h3 = 2.; // l - m
         double h1 = cos_theta * (h2 + h3); // x * (2l - 1)
         for (int l = m+2; /* done inside */; ++l)  {
            // P(l,m) = (x*(2l-1)*P(l-1,m) - (l+m-1)*P(l-2,m))/(l-m)
            P2 = (h1 * P1 - h2 * P0) / h3;
            P(l*(l+1)+m)=P2;
            if (l == lmax) break;
            P0 = P1;
            P1 = P2;
            h1 += 2. * cos_theta;
            h2 += 1.;
            h3 += 1.;
         }
      } // lmax > m+1
   } // lmax > m
}

double SxYlm::getYlmNormFactor (int l, int m)
{
   SX_CHECK (l >= 0, l);
   m = abs(m);
   SX_CHECK (m <= l, m, l);
   // The normalization factor is
   // N = sqrt((2l + 1)/4pi * (l-m)! / (l+m)!)
   return sqrt(double(2*l + 1) / FOUR_PI * facFrac(l - m, l + m) );
}

void SxYlm::getYlmArray(int lmax, double x, double y, double z,
                 SxVector<Double> *resultPtr)
{
   SX_CHECK (resultPtr);
   SX_CHECK (lmax >= 0, lmax);
   SxVector<Double> &result = *resultPtr;
#ifndef NDEBUG
   int nl = lmax + 1;
#endif
   SX_CHECK (result.getSize () == nl * nl, result.getSize (), nl);
   if (lmax == 1)  {
      result(0) = 1.;
      double ri = 1./sqrt(x*x+y*y+z*z);
      result(1) = SQRT2 * ri * y;
      result(2) = ri * z;
      result(3) = SQRT2 * ri * x;
      return;
   }

   double cos_theta;
   // for m<>0 extra factor of sqrt(2) for REAL Ylm (compared to COMPLEX)
   SxComplex16 eiphi, eimphi(SQRT2,0.); 
   {
      double s2 = x*x + y*y;
      SX_CHECK (s2 + z*z > 1e-12); // 
      cos_theta = z / sqrt(s2 + z*z);
      if (s2 > 1e-12)  {
         double s = 1./sqrt(s2);
         eiphi.re = x*s;
         eiphi.im = y*s;
      } else {
         eiphi.re = 1.;
         eiphi.im = 0.;
      }

   }
   ::getPlmArray(lmax, cos_theta, &result);

   for (int m = 1; m <= lmax; m++)  {
      eimphi *= eiphi;
      // getPlmArray(lmax, m, cos_theta, &P);
      for (int l = m; l <= lmax; l++)
      {
         // double Plm = P(l);
         double Plm = result(l*(l+1)+m);
         if (m & 1) Plm = -Plm; // (-1)^m, so we get +iy, +iz, +ix
         result(l*(l+1) + m) = Plm * eimphi.re;
         result(l*(l+1) - m) = Plm * eimphi.im;
      }
   }
   VALIDATE_VECTOR(result);
}

double SxYlm::faculty (int n)
{
   SX_CHECK (n >= 0, n);
   if (n < 2) return 1.;
   double res = 2.;
   double fac = 3.;
   for (int i = 3; i <= n; ++i, fac += 1.)
      res *= fac;
   return res;
}

double SxYlm::facFrac (int n, int m)
{
   SX_CHECK (n >= 0, n);
   SX_CHECK (m >= 0, m);
   if (n == m) return 1.;
   if (n > m)  {
      double res = m + 1;
      double fac = res + 1.;
      for (int i = m + 2; i <= n; ++i, fac += 1.)
         res *= fac;
      return res;
   } else {
      double res = n + 1;
      double fac = res + 1.;
      for (int i = n + 2; i <= m; ++i, fac += 1.)
         res *= fac;
      return (1. / res);
   }
}

double SxYlm::wigner3j (int j1, int j2, int j3, int m1, int m2, int m3)
{
   SX_CHECK (j1 >= 0, j1);
   SX_CHECK (j2 >= 0, j2);
   SX_CHECK (j3 >= 0, j3);
   SX_CHECK (abs(m1) <= j1, m1, j1);
   SX_CHECK (abs(m2) <= j2, m2, j2);
   SX_CHECK (abs(m3) <= j3, m3, j3);
   if (m1 + m2 + m3 != 0)
      return 0.;
   int l1 = - j1 + j2 + j3;
   if (l1 < 0) return 0.;
   int l2 = + j1 - j2 + j3;
   if (l2 < 0) return 0.;
   int l3 = + j1 + j2 - j3;
   if (l3 < 0) return 0.;

   int d1 = j1 - m1;
   int d2 = j2 + m2; // like this!
   
   int fromK = sxmax(j2 - j3 - m1, j1 - j3 + m2, 0);
   int toK = sxmin(l3, d1, d2);

   double sum = 0.;
   bool sign = (j1 + j2 + m3 + fromK) & 1; // true = -, false = +

   {
      double summand;
      for (int k = fromK; k <= toK; k++)  {
         summand = facFrac(l1, l1 - d2 + k) * facFrac(l2, l2 - d1 + k);
         summand *= facFrac(l3, l3 - k);
         summand /= faculty(k) * faculty(d1 - k) * faculty(d2 - k);
         if (sign)
            sum -= summand;
         else
            sum += summand;
         sign = !sign;
      }
   }

   double frac;
   frac = facFrac(j1+m1,l1) * facFrac(j2+m2,l2) * facFrac(j3+m3,l3);
   frac *= facFrac (j3-m3, 1+j1+j2+j3); 
   frac *= faculty (j1-m1) * faculty(j2-m2);

   return sum * sqrt(frac);
}

SxYlm::SxClebschTable SxYlm::getClebschGordan (int l1max, int l2max, int l3max,
                               enum SxYlm::YlmType type)
{
   SX_CHECK (l1max >= 0, l1max);
   SX_CHECK (l2max >= 0, l2max);
   SX_CHECK (l3max >= 0, l3max);
   SX_CHECK(type == RealYlm || type == ComplexYlm);
   int n1 = (l1max + 1) * (l1max + 1);
   int n2 = (l2max + 1) * (l2max + 1);
   int n3 = (l3max + 1) * (l3max + 1);

   SxNArray<SxComplex16,3> res(n1, n2, n3);
   SxClebschTable dres(n1, n2, n3);
   int i1,i2,i3;

   /*
      The Clebsch-Gordan coefficients of the complex spherical harmonics are
      related to WIgner 3j symbols by
                                           ( l1 l2 l3 )   ( l1 l2  l3 )
      <l1 l2 m1 m2| l3 m3> = C(l1,l2,l3) * ( 0  0  0  ) * ( m1 m2 -m3 )
      C(l1,l2,l3) = (-1)^m3  sqrt((2l1+1)(2l2+1)(2l3+1)/4pi)
      */
   int l1,l2,l3,m1,m2,m3;
   double C;
   for (l1 = 0; l1 <= l1max; l1++)  {
      for (m1 = -l1; m1 <= l1; m1++)  {
         i1 = combineLm(l1,m1);
         for (l2 = 0; l2 <= l2max; l2++)  {
            for (m2 = -l2; m2 <= l2; m2++)  {
               i2 = combineLm(l2,m2);
               for (l3 = 0; l3 <= l3max; l3++)  {
                  C = sqrt(double((2 * l1 + 1) * (2 * l2 + 1) * (2 * l3 + 1)) 
                           / FOUR_PI) 
                      * wigner3j(l1,l2,l3,0,0,0);
                  for (m3 = -l3; m3 <= l3; m3++)  {
                     i3 = combineLm(l3,m3);
                     res(i1, i2, i3) = C * wigner3j(l1,l2,l3,m1,m2,-m3);
                     // -1 ^ m3
                     if (m3 & 1) 
                        res(i1, i2, i3) = -res(i1, i2, i3);
                     // i^l for even l
                     if ( abs(l1 + l2 - l3) % 4 == 2 )
                        res(i1, i2, i3) = -res(i1, i2, i3);
                  } // m3
               } // l3
            } // m2
         } // l2
      } // m1
   }  // l1

   if (type == RealYlm)  {
      // Now we have to change to real spherical harmonics
      
      SxComplex16 pos,neg;
      double sign;
      int ip,in;
      const double norm = 1. / SQRT2;
      for (l1 = 0; l1 <= l1max; l1++)  {
         for (m1 = 1; m1 <= l1; m1++)  {
            ip = combineLm(l1,m1);
            in = combineLm(l1,-m1);
            sign = (m1 & 1) ? -1. : 1.;
            for (i2 = 0; i2 < n2; i2++)  {
               for (i3 = 0; i3 < n3; i3++)  {
                  pos = norm * (res(ip, i2, i3) + sign * res(in, i2, i3));
                  neg = norm * (res(ip, i2, i3) - sign * res(in, i2, i3));
                  res(ip, i2, i3) = pos;
                  res(in, i2, i3) = -I * neg;
               }
            }
         }
      }
      for (l2 = 0; l2 <= l2max; l2++)  {
         for (m2 = 1; m2 <= l2; m2++)  {
            ip = combineLm(l2,m2);
            in = combineLm(l2,-m2);
            sign = (m2 & 1) ? -1. : 1.;
            for (i1 = 0; i1 < n1; i1++)  {
               for (i3 = 0; i3 < n3; i3++)  {
                  pos = norm * (res(i1, ip, i3) + sign * res(i1, in, i3));
                  neg = norm * (res(i1, ip, i3) - sign * res(i1, in, i3));
                  res(i1, ip, i3) = pos;
                  res(i1, in, i3) = -I * neg;
               }
            }
         }
      }
      for (l3 = 0; l3 <= l3max; l3++)  {
         for (m3 = 1; m3 <= l3; m3++)  {
            ip = combineLm(l3,m3);
            in = combineLm(l3,-m3);
            sign = (m3 & 1) ? -1. : 1.;
            for (i1 = 0; i1 < n1; i1++)  {
               for (i2 = 0; i2 < n2; i2++)  {
                  pos = norm * (res(i1, i2, ip) + sign * res(i1, i2, in));
                  neg = norm * (res(i1, i2, ip) - sign * res(i1, i2, in));
                  res(i1, i2, ip) = pos;
                  res(i1, i2, in) = I * neg;
               }
            }
         }
      }
   }

   for (i1 = 0; i1 < n1; i1++)
      for (i2 = 0; i2 < n2; i2++)
         for (i3 = 0; i3 < n3; i3++)  {
            SX_CHECK (fabs(res(i1, i2, i3).im) < 1e-12, i1, i2, i3);
            dres(i1, i2, i3) = res(i1, i2, i3).re;
         }
   
   return dres;
}

SxArray<SxMatrix<Double> > 
SxYlm::computeYlmRotMatrices (int lmax, const SxMatrix3<Double> &S,
      const SxClebschTable &complexClebschGordan)
{
   SX_CHECK (lmax >= 0, lmax);
   // check that S is unitary
   SX_CHECK (fabs(fabs(S.determinant ()) - 1.) < 1e-6, fabs(fabs(S.determinant ()) - 1.));
   SX_CHECK (complexClebschGordan.getDim (0) >= lmax *lmax,
             complexClebschGordan.getDim (0), lmax);
   SX_CHECK  (complexClebschGordan.getDim (1) >= 4,
              complexClebschGordan.getDim (1));
   SX_CHECK (complexClebschGordan.getDim (2) >= (lmax+1)*(lmax + 1),
             complexClebschGordan.getDim (2), lmax);
   SxArray<SxMatrix<Double> > res(lmax + 1);
   // l = 0
   res(0).reformat (1,1);
   res(0)(0) = 1.;
   if (lmax == 0) return res;
   // l == 1
                              //  py      pz      px
   res(1) = SxMatrix3<Double> ( S(1,1), S(1,2), S(1,0),  // py
                                S(2,1), S(2,2), S(2,0),  // pz
                                S(0,1), S(0,2), S(0,0)); // px
   int M,N,m,n, lM,lN, dm, dn, ldm, ldn;

   SxMatrix<Complex16> D1 = res(1), Dold, D;
   // change from real to complex
   double norm = sqrt(.5);
   SxComplex16 pos, neg;
   for (m = 0; m <= 2; ++m)  {
      pos = D1(m,2); neg = D1(m,0);
      D1(m,0) = norm * ( pos + I * neg);
      D1(m,2) = norm * (-pos + I * neg);
   }
   for (m = 0; m <= 2; ++m)  {
      pos = D1(2,m); neg = D1(0,m);
      D1(0,m) = norm * ( pos - I * neg);
      D1(2,m) = norm * (-pos - I * neg);
   }
   //cout << D1 << endl;

   Dold = D1;
   // recursion
   // D[L]_{MN} = sum_{mn} <LN|(L-1)(N-n) 1n> D[L-1]_{M-m,N-n} D[1]_{mn} 
   //                      <(L-1)(M-m) 1m|LM> * normCG
   // normCG := 4 PI / (2(L-1)+1) / (2*1+1) / Wigner3j (L-1,1,L,0,0,0)^2
   for (int l = 2; l <= lmax; ++l)  {
      double normCG = FOUR_PI 
                    / ((2*l - 1) * (2*1+1) * sqr(wigner3j (l-1,1,l,0,0,0)));
      D.reformat (2*l+1,2*l+1);
      D.set (0.);
      //cout << "l:" << l << endl;
      for (M = -l; M <= l; ++M)  {
         lM = combineLm(l,M);
         for (N = -l; N <= l; ++N)  {
            lN = combineLm (l, N);
            for (m = -1; m <= 1; ++m)  {
               dm = M-m;
               if (abs(dm) > l - 1) continue;
               ldm = combineLm(l-1,dm);
               for (n = -1; n <= 1; ++n)  {
                  dn = (N-n);
                  if (abs(dn) > l - 1) continue;
                  ldn = combineLm(l-1,dn);
                  D(l+M,l+N) += Dold(l-1+dm,l-1+dn) * D1(1+m,1+n)
                              * complexClebschGordan(ldm, combineLm(1,m), lM)
                              * complexClebschGordan(ldn, combineLm(1,n), lN)
                              * normCG;
               }
            }
         }
      }
      /*
      cout << "complex:\n";
      for (M = -l; M <= l; ++M)  {
        for (N = -l; N <= l; ++N)  {
          cout << D(l+M,l+N) << " ";
        }
        cout << endl;
      }
      */
      Dold = SxMatrix<Complex16> (); // bug in copy () for ill-sized vectors??
      Dold.copy (D);
      for (M = -l; M <= l; ++M)  {
         for (N = 1; N <=l; ++N)  {
            double sign = (N & 1) ? -1. : 1.;
            pos = norm * (sign * D(l+M,l+N) + D(l+M,l-N));
            neg = norm * (sign * D(l+M,l+N) - D(l+M,l-N));
            D(l+M,l+N) = pos;
            D(l+M,l-N) = I * neg;
         }
      }
      for (M = 1; M <=l; ++M)  {
         for (N = -l; N <= l; ++N)  {
            double sign = (M & 1) ? -1. : 1.;
            pos = norm * (sign * D(l+M,l+N) + D(l-M,l+N));
            neg = norm * (sign * D(l+M,l+N) - D(l-M,l+N));
            D(l+M,l+N) = pos;
            D(l-M,l+N) = -I * neg;
         }
      }
      res(l) = D.real ();
      /*
      cout << "real:\n";
      for (M = -l; M <= l; ++M)  {
        for (N = -l; N <= l; ++N)  {
          cout << D(l+M,l+N) << " ";
        }
        cout << endl;
      }
      */
      SX_CHECK(D.imag ().normSqr () < 1e-10, D.imag ().normSqr ());
   }
   return res;
}

SxYlmRotGroup
SxYlm::computeYlmRotMatrices (const SxArray<SxMatrix3<Double> > &syms, int lmax)
{
   ssize_t nSym = syms.getSize ();
   SxYlmRotGroup DS(nSym);
   int iSym, jSym, kSym;
   // --- set up rotation matrices
   // what if lmax = zero like in case for hydrogen ???
   SxClebschTable clebschGordan;
   if (!lmax) 
      clebschGordan = getClebschGordan (lmax, 1, lmax,
                                                  SxYlm::ComplexYlm);
   else
      clebschGordan = getClebschGordan (lmax-1, 1, lmax,
                                                  SxYlm::ComplexYlm);
   for (iSym = 0; iSym < nSym; ++iSym)
      DS(iSym) = computeYlmRotMatrices (lmax, syms(iSym), clebschGordan);
   // Validate group representation for L-rotation matrices 
   SxArray<double> dMax(lmax+1);
   dMax.set (0.);
   double diff;
   for (iSym = 0; iSym < nSym; ++iSym)  {
      for (jSym = iSym; jSym < nSym; ++jSym)  {
         SxMatrix3<Double> S = syms(iSym) ^ syms(jSym);
         for (kSym = 0; kSym < nSym; ++kSym)
            if ((S - syms(kSym)).absSqr ().sum () < 1e-12) break;
         if (kSym == nSym)  {
            cout << "iSym=" << iSym << "; jSym = " << jSym << endl;
            cout << "Symmetry group validation error" << endl;
            cout << "S=" << S << endl;
            cout << "S is product of two symmetries, but not in group." << endl;
            SX_EXIT;
         }
         for (int l = 1; l <= lmax; ++l)  {
            diff = ((DS(iSym)(l) ^ DS(jSym)(l)) - DS(kSym)(l)).absSqr ().sum ();
            if (diff > dMax(l)) dMax(l)=diff;
            if (diff > 1e-16)  {
               cout << "Error for rotation group in Ylm-space" << endl;
               cout << "Validation error of D-group for l=" << l << endl;
               cout << "iSym=" << iSym << "; jSym = " << jSym 
                    << "=> " << kSym << endl;
               cout << "l=" << l << ": diff=" << diff << endl;
               SX_EXIT;
            }
         }
      }
   }
   /*
   for (int l = 0; l <= lmax; ++l)
      cout << "Max error for l = " << l << ": " << dMax(l) << endl;
   */
   return DS;
}


