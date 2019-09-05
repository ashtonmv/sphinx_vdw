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

// References
// (1) M. Ernzerhof, J. P. Perdew, Generalized gradient approximation to the
//     angle- and system-averaged exchange hole, J. Chem. Phys. 109,3313 (1998).
// (2) HSE implementation notes
// (3) T. M. Henderson et al., J. Chem. Phys. 131, 044108 (2009)
#include <SxPBExHole.h>
#include <SxMathLib.h>
#include <SxCLI.h>
#include <SxConstants.h>

// Ref. 1, page 3315, text below Eq. (19)
const double SxPBExHole::A =  1.0161144,
             SxPBExHole::B = -0.37170836,
             SxPBExHole::C = -0.077215461,
             SxPBExHole::D =  0.57786348,
             SxPBExHole::E = -0.051955731;
// Ref. 1, page 3320, text below Eq. (A5)
const double SxPBExHole::a1 = 0.00979681,
             SxPBExHole::a2 = 0.0410834,
             SxPBExHole::a3 = 0.187440,
             SxPBExHole::a4 = 0.00120824,
             SxPBExHole::a5 = 0.0347188;

void SxPBExHole::computeHs (double s, double s2, double s4)  {
   // Eq. A5
   //    H = (a1 + a2 * s2) * s2
   //      / (1. + s4 * (a3 + s * a4 + s2 * a5));
   double denom = 1. + s4 * (a3 + s * a4 + s2 * a5);
   Hs = (2. * a1 + 4. * a2 * s2) * s / denom
      - (a1 + a2 * s2) * s2 / (denom * denom) 
        * s2 * s * (4. * a3 + 5. * a4 * s + 6. * a5 * s2);
}

void SxPBExHole::computeFGH (double s, bool deriv)
{
#ifndef NDEBUG
   if (!deriv)  {
      double nan = sqrt(-1.);
      Fs = Gs = Hs = nan;
   }
#endif
   SX_CHECK (s >= 0., s);
   // get s2 = s^2, s4 = s^4
   double s2 = s * s, s4 = s2 * s2;
   computeH (s, s2, s4);
   if (deriv) computeHs (s, s2, s4);
   // Ref 1, Eq. (25)
   F = (H * ( 16. *A*A + 36. * (B - A * D)) - 4. / 3. ) / (36. * C);
   //double F = 6.475 * H + 0.4797;
   if (deriv)
      Fs = Hs * ( 16. *A*A + 36. * (B - A * D)) / (36. * C);

   // ref (2), Eq. (15)
   double Ds = D + H * s2,
          denom = 16. * Ds * Ds * Ds * sqrt(Ds);
          //denom = 16. * pow(Ds, 7./2.);
          // Ref 1, Eq. A2
   double denoms = 0., Dss = 0.;
   if (deriv) {
      Dss = 2. * H * s + Hs * s2;
      denoms = 56. * Ds * Ds * sqrt(Ds) * Dss;
   }
   double Ps, Es, Pss = 0., Ess = 0.;
   double a = SQRT_PI * (Ps = (15. * E
                               + Ds * (6. * C * (1. + F * s2) 
                                         + Ds * (4. * B + 8. * A * Ds))))
              / denom
            - 0.75 * PI * sqrt(A) * (Es = exp (9. * H * s2 / (4. * A))
              * (derfc (1.5 * s * sqrt (H / A)))),
          // Ref 1, Eq. A3
          b = 15. * sqrt(PI) * s2 / denom;
   double limit0 = 3e-3;
   if (deriv && s > limit0)  {
      Pss = Dss * (6. * C * (1. + F * s2) 
                     + Ds * (8. * B + 24. * A * Ds))
          + Ds * (6. * C * (2. * F * s + Fs * s2)); 
      Ess = Es * 9. * Dss / (4. *A)
          - 1.5 * (Hs * s + 2. * H) / (SQRT_PI * sqrt(H*A));
   }
   //cout << "##" << s << '\t' << Es << '\t' << Ess << endl; 
   // --- Ref. 1, Eq A1  
   if (s > 1.)  {
      // avoid numerical noise: replace -0.75*PI by lim_{s->0} a(s) 
      double a0 = SQRT_PI * (15. * E + D * (6. * C + D * (4. * B + 8. * A * D)))
                  / (16. * pow(D, 3.5)) - 0.75 * PI * sqrt(A);
      //G = -(0.75 * PI + a) / (b * E);
      G = (a0 - a) / (b * E);
      if (deriv)  {
         // a = SQRT_PI * Ps / denom - 0.75 * PI * sqrt(A) * Es;
         double as = SQRT_PI * (Pss / denom - Ps * denoms / sqr(denom))
                   - 0.75 * PI * sqrt(A) * Ess;

         // Ref 1, Eq. A3
         double bs = 15. * sqrt(PI) * (2. * s / denom - s2 * denoms / sqr(denom));
         Gs = (a - a0) / (sqr (b) * E) * bs
            - as / (b * E) ;
      }
   } else if (s > limit0)  {
      // a(s) should go to -0.75 * PI + O(s^2)
      // avoid numerical noise: replace -0.75*PI by lim_{s->0} a(s) 
      double P0 = 15. * E + D * (6. * C + D * (4. * B + 8. * A * D));
      double denom0 = 16. * pow(D, 3.5);
      G = (Ps - denom / denom0 * P0) / (15. * s2)
        - 0.05 * sqrt(PI * A) * (Es - 1.) * denom / s2;
      G /= -E;
      if (deriv)  {
         Gs = (Pss - denoms / denom0 * P0) / (15. * s2)
            - 2. * (Ps - denom / denom0 * P0) / (15. * s2 * s)
            - 0.05 * sqrt(PI * A) * (Ess * denom + (Es-1.) * denoms) / s2 
            + 0.1 * sqrt(PI * A) * (Es - 1.) * denom / (s2 * s);
         Gs /= -E;
      }
   } else {
      // numerical fit  (range 0.004 to 0.01) to expression above
      G = 0.50589566 + 1.36993437 * s2;
      if (deriv)  {
         Gs = 1.36993437 * 2. * s;
      }
   }
}

double SxPBExHole::exchangeHolePBE(double s, double y)
{
   computeFGH (s);
   double y2 = y * y, s2 = s*s;
   // Ref. 1, Eq. (24)
   /*
   double J = (- A / (y2 * (1. + 4./9. * A * y2))
               +(A/y2 + B + C * (1. + F * s2) * y2
                 + E * (1. + G * s2) * y2 * y2) * exp(- D * y2))
            * exp (- s2 * H * y2);
   */
   double eDy  = exp (- D * y2),
          eDy1 = (y > 1e-2) ? (eDy - 1.) / y2
                             : D * (-1. + 0.5 * D * y2 * (1. - D/3. * y2));
   double J = (A/(1. + 4./9. * A * y2) * (eDy1 + 4./9. * A * eDy)
               + (B  + C * (1. + F * s2) * y2 + E * (1. + G * s2) * y2 * y2)
                 * exp(- D * y2))
              * exp (- s2 * H * y2);

   return J;
}

void SxPBExHole::test (int argc, char **argv)
{
   SxCLI cli(argc, argv);
   double s0 = cli.option ("-s0", "number", "start value for s").toDouble (0.);
   double smax = cli.option ("-smax", "number", "max value for s").toDouble (3.);
   double ds = cli.option ("-ds", "number", "delta s").toDouble (0.5);
   int yGroup = cli.newGroup ("Test J(s,y)");
   double y0 = cli.option ("-y0", "number", "start value for y").toDouble (0.);
   double ymax = cli.option ("-ymax", "number", "max value for y").toDouble (10.);
   double dy = cli.option ("-dy", "number", "delta y").toDouble ();
   int nuGroup = cli.newGroup ("Test Fx(s,nu)");
   double nu = cli.option ("-nu", "number", "nu").toDouble ();

   cli.finalize ();

   cout.precision (16);
   if (cli.groupAvailable (nuGroup))  {
      cout << "# nu = " << nu << endl;
      cout << "# s Fx(s,nu)  (d Fx/ds)  (d Fx/d nu)" << endl;
      for (double s = s0; s <= smax; s += ds)  {
         cout << s << '\t';
         screenedFx (s, nu);
         cout << Fx << '\t' << Fxs << '\t' << FxNu;
         if (cli.groupAvailable (yGroup))  {
            double sum = 0.;
            for (double y = y0; y < ymax; y += dy)  {
               sum += exchangeHolePBE (s, y) * y * derfc (nu * y);
            }
            sum *= -8./9. * dy;
            cout << '\t' << sum;
         }
         cout << endl;
      }
   } else if (cli.groupAvailable (yGroup))  {
      for (double s = s0; s <= smax; s += ds)  {
         cout << "# s = " << s << endl;
         cout << "# y J(s,y)" << endl;
         for (double y = y0; y < ymax; y += dy)  {
            cout << y << '\t' << exchangeHolePBE (s, y) << endl;
         }
         cout << "&\n";
      }
   } else {
      cout << "# s F G H" << endl;
      for (double s = s0; s < smax; s += ds)  {
         computeFGH (s, true);
         cout << s << '\t' << F << '\t' << G << '\t' << H;// << endl;
         cout << '\t' << Fs << '\t' << Gs << '\t' << Hs << endl;
      }
   }
}

double SxPBExHole::screenedFx (double s, double nu)
{
   SX_CHECK (s>= 0., s);
   if (s < 1e-8)  {
      screenedFx (1e-8, nu);
      Fx  += 0.5 * (Fxs * 1e8) * sqr(s - 1e-8);
      Fxs += (Fxs * 1e8) * (s - 1e-8);
      // Fxnu is accurate to O(s^2)
      return Fx;
   }
   computeFGH (s, true);
   SX_CHECK (nu >= 0., nu);
   double s2 = s * s;
   // ref (2), Eq. (15)
   double Ds = D + H * s2, Ds2 = Ds * Ds,
          Cs = C * (1. + F * s2),
          Es = E * (1. + G * s2),
   // ref (2), Eq. (16)
          eta2 = 1. / (Ds + nu * nu), eta = sqrt(eta2),
          eta2Ds = eta2 * Ds;


   double FxB, FxC, FxE;
   // --- ref (2), Eq. (14)
   FxB = (1. - nu * eta) / (2. * Ds);
   FxC = (1. - nu * eta * (1. + 0.5 * eta2Ds)) / (2. * Ds2);
   FxE = (1. - nu * eta * (1. + eta2Ds * (0.5  + 0.375 * eta2Ds)))
         / ( Ds * Ds2);

   Fx = FxB * B + FxC * Cs + FxE * Es;

   // Ref (3), Eq. A5 (\overline A = 4/9 A, see text below Eq. 4c)
   double alpha = 2.25 *s*s * H / A,
          beta  = 2.25 * nu * nu / A, sBeta = sqrt(beta),
          gamma = 2.25 * D / A,
          sabc = sqrt(alpha + beta + gamma);

   // Ref (3), Eq. A7
   //  I(s,nu) without -8/9 prefactor
   double K = hendersonK (alpha, beta);
   Fx += -A * (log (sBeta + sabc)  - K);
   SX_CHECK_NUM (Fx);

   // s-derivatives of Ds, eta, alpha
   double Dss    = 2. * H * s + Hs * s2,
          etaS   = -0.5 * Dss * eta2 * eta,
          alphaS = 2.25 * s * (2. * H + s * Hs) / A;

   // --- s-derivatives
   Fxs = 0.;
   Fxs  = (-FxB * Dss - 0.5 * nu * etaS) / Ds * B;
   Fxs += (-FxC * 2. * Dss 
           - 0.5 * nu/Ds * (etaS + 0.5 * eta2 * (3. * etaS * Ds + eta * Dss)) )
          / Ds * Cs
        + FxC * C * (F * 2. * s + Fs * s2);
   Fxs += (-FxE * 3. * Dss
           - (          nu *                       etaS 
              + 0.5   * nu * eta2 *          (3. * etaS * Ds + eta * Dss)
              + 0.375 * nu * eta2 * eta2Ds * (5. * etaS * Ds + 2. * eta * Dss)
             ) / Ds2) / Ds
          * Es
        + FxE * E * (G * 2. * s + Gs * s2);

   // Derivative of I(s,nu), see Ref (3), Eq. A10 for K-derivative
   Fxs += -A * (0.5/((sBeta + sabc)*sabc)
                    - K + log (sBeta + sqrt(alpha + beta)))
          * alphaS;
   SX_CHECK_NUM (Fxs);

   // --- nu-derivatives
   double etaNu = - nu * eta * eta2;

   FxNu = 0.;
   FxNu  = - B  * (eta + nu * etaNu) / (2. * Ds);
   FxNu += - Cs / (2. * Ds2) 
         * (                               eta +      etaNu * nu
            + 0.5   * eta2 *        Ds  * (eta + 3. * etaNu * nu));
   FxNu += - Es / (Ds2 * Ds) 
         * (                               eta +      etaNu * nu
            + 0.5   * eta2 *        Ds  * (eta + 3. * etaNu * nu)
            + 0.375 * eta2 * eta2 * Ds2 * (eta + 5. * etaNu * nu));

   // Derivative of I(s,nu), see Ref (2), Eq. 29 for K-derivative
   //Fxnu += -A * (0.5/(sBeta*sabc) - SQRT_PI/(2.*sBeta)
   //          * exp(alpha+beta)*erfc (sqrt(alpha + beta)));
   FxNu += (-A * (0.5/(sBeta*sabc) 
            - SQRT_PI/(2.*sBeta) * erfcexp (sqrt(alpha + beta)))) 
         * 4.5 * nu / A;
   SX_CHECK_NUM (FxNu);

   // --- prefactor
   double pre = -8./9.;
   Fx   *= pre;
   Fxs  *= pre;
   FxNu *= pre;
   return Fx;
}

double SxPBExHole::hendersonK (double alpha, double beta)
{
   const int N = 20;
   SX_CHECK (alpha > 0., alpha);
   SX_CHECK (beta > 0., beta);
   // compute
   // K(a,b) = \int_0^\infty dz exp(-z) log (sqrt(b) + sqrt(a + b + z))
   double sqrtBeta = sqrt(beta), ab = alpha + beta;

   if (gaussLegendre.getSize () != N)  {
      gaussLegendre.setupLegendre (N);
      gaussLaguerre.setupLaguerre (N);
   }
   SX_CHECK (gaussLegendre.getSize () == N, gaussLegendre.getSize ());
   SX_CHECK (gaussLaguerre.getSize () == N, gaussLaguerre.getSize ());

   double K = 0.;
   // --- integrate from 4 to infinity
   for (int i = 0; i < gaussLaguerre.getSize (); ++i)  {
      double z = gaussLaguerre.nodes(i) + 4.;
      K += gaussLaguerre.weights(i) * log (sqrtBeta + sqrt(ab + z));
   }
   K *= exp(-4.);

   // --- integrate from 5/8 to 4
   double dK = 0.;
   for (int i = 0; i < gaussLegendre.getSize (); ++i)  {
      double z = gaussLegendre.nodes(i) * 1.6875 + 2.3125;
      dK += gaussLegendre.weights(i) * log (sqrtBeta + sqrt(ab + z))
            * exp (-z);
   }
   K += 1.6875 * dK;

   // --- integrate from 0 to 5/8
   dK = 0.;
   if (ab > 0.125)  {
      // Gauss-Legendre
      for (int i = 0; i < gaussLegendre.getSize (); ++i)  {
         double z = gaussLegendre.nodes(i) * 0.3125 + 0.3125;
         dK += gaussLegendre.weights(i) * log (sqrtBeta + sqrt(ab + z))
               * exp (-z);
      }
      K += 0.3125 * dK;
   } else {
      // integration by parts
      K -=  exp(-0.625) * log(sqrtBeta + sqrt(ab + 0.625))
           -              log(sqrtBeta + sqrt(ab));
      // => leaves e^{-z}( 1/(a+z) - sqrt(beta)/((a+z)sqrt(ab+z)) )

      // [e^{-z} - P4(z)] (1 -  sqrt(beta)/sqrt(ab+z)) / (a+z)
      for (int i = 0; i < gaussLegendre.getSize (); ++i)  {
         double z = gaussLegendre.nodes(i) * 0.3125 + 0.3125;
         dK += gaussLegendre.weights(i) 
         * (1. - sqrtBeta / sqrt(ab + z))/(alpha + z)
               * (exp (-z) - P4(z));
      }
      K += 0.3125 * 0.5 * dK;

      // --- integral P4(z)
      // (-alpha)^n 
      K += 0.5 * P4 (-alpha) * log ((alpha+0.625)/alpha);
      
      // --- (-1)^n/n! z^n/(a+z)
      // Ref (2), Eq. (19)
      double an = alpha, u = 1. + 0.625 / alpha, facn = 1.;
      for (int n = 1; n <= 4; ++n, an *= alpha)  {
         dK = 0.;
         facn *= n;
         double uk = u, fac = 1./facn;
         for (int k = 1; k <= n; ++k, uk *= u)  {
            fac *= -(n+1-k) / double(k);
            dK += (uk - 1.) * fac / k ;
         }
         K += 0.5 * an * dK;
      }

      // --- (-1)^n/n! z^n/(a+z) * sqrt(beta/(a+b+z))
      K -= 0.5 * (  integralP4 (alpha, beta, 0.625)
                  - integralP4 (alpha, beta, 0.));
   }
   SX_CHECK_NUM (K);
   return K;
}


double SxPBExHole::integralP4 (double alpha, double beta, double z)
{
   double x = alpha + beta;
   double p2, p3, p4;
   // Ref 2, Eq. (23)
   p2 = -alpha      - 2./3.  * x                 + z/3.;
   // Ref 2, Eq. (24)
   p3 = -alpha * p2 - 2./15. * x * (-4. * x + 2. * z) + z*z/5.;
   // Ref 2, Eq. (25)
   p4 = -alpha * p3 - 2./35. * x * (8. * x*x - 4. * x * z + 3.*z*z) + z*z*z/7.;

   double ll = log ( (alpha + z) / sqr(sqrt(x+z) + sqrt(beta)));

   // Ref 2, Eq. (20) with z^n prefactors
   return 2. * sqrt(beta * (x+z)) * (-1. + p2 / 2. - p3 / 6. + p4 / 24.)
          + P4(-alpha) * ll; 
}

