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

#include <SxXCFunctional.h>
#include <SxPrecision.h>
#include <SxString.h>
#include <SxFileIO.h>

//------------------------------------------------------------------------------
// Bibliography:
//
// ref1 - PRB 23, 5048 (1981)
// ref2 - PRB 54, 16533 (1996)
// ref3 - K. Burke's original code, in src/api/k.burkes.f
// ref4 - PRB 33, 8800 (1986)
//        J.P. Perdew, Y. Wang, Accurate and simple density functional for the
//        electronic exchange energy: Generalized gradient approximation
// ref5 - PRB 45, 13244 (1992)
//        J.P. Perdew, Y. Wang, Accurate and simple analytic expression of the
//        electron-gas correlation energy
// ref6 - PRL 77, 3865 (1996)
// ref7 - "Electronic Structure of Solids", 11-20 (1991), 
//        ed. P.Ziesche, H.Eschrig (Akademie Verlag, Berlin)
// ref8 - http://www.sfhingx.de/documentation/lectures/pbe_white_bird.{ps,pdf}
// ref9 - PRB 59, 10031 (1999)
// ref10 - B. Grabowski, diploma thesis
// ref11 - Henderson et al., J. Chem. Phys. 131, 044108 (2009)
//------------------------------------------------------------------------------

// invalidate omegaHSE; must be set before using HSE
double SxXCFunctional::omegaHSE = sqrt(-1.);

// default value for exact-exchange mixing in hybrid functionals
double SxXCFunctional::alphaHybrid = 0.25;

void SxXCFunctional::resize (int nSpinIn)
{
   nSpin = nSpinIn;
   vXc.resize (nSpin);
#ifndef NDEBUG
   double nan = sqrt(-1.);
   eXc = epsXc = nan;
   vXc.set (nan);
#endif
   if (mode & ComputeKernel)  {
      fXc.resize (nSpin);
      for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
         fXc(iSpin).resize (nSpin);
#        ifndef NDEBUG
         fXc(iSpin).set (nan);
#        endif
      }
   }
}

//==============================================================================
// LDA - Perdew/Zunger '81
//==============================================================================
void SxXCFunctional::computeLDA (double rho)
{
   SX_CHECK (nSpin == 1, nSpin);

   const double A  = 0.0311,  B  = -0.048;  // ref1, (C7)
   // --- ref1, (table xii)
   const double g  = -0.1423;
   const double b1 =  1.0529;
   const double b2 =  0.3334;
   const double C  =  0.0020;
   const double D  = -0.0116;

   double rS = getRs(rho);
   if (mode & ComputeC)  {
      // --- correlation part
      double rrs, lrs, bpr;
      if ( rS > 1.0 )  {
         rrs  = sqrt (rS);
         bpr  = 1. + b1 * rrs + b2 * rS;
         epsXc = g / bpr;                 // ref1, (C3, unpol)
         if (mode & (ComputePotential | ComputeKernel))  {

            vXc(0) = (1. + 7./6.*b1 * rrs + 4./3. * b2 * rS)
                   *  epsXc / bpr;             // ref1, (C4, unpol)
            SX_CHECK_NUM(vXc(0));
         }
         if (mode & ComputeKernel)  {
            /*
            fXc(0)(0) = ((1. + 7./6. * b1 * rrs + 4./3. * b2 * rS)
                          * (b1 / 3. * rrs + 2./3. * b2 * rS) / bpr 
                          - (7. / 36. * b1 * rrs + 4. / 9. * b2 * rS))
                      * epsXc / (bpr * rho);
            */
            if (rho > 0.)  {
               fXc(0)(0) = (vXc(0) * (b1 * rrs / 3. + 2. * b2 * rS / 3.) 
                            + epsXc * (-7./36.*b1*rrs -4./9.*b2*rS))
                         / (bpr * rho);
               SX_CHECK_NUM(fXc(0)(0));
            } else {
               fXc(0)(0) = 0.;
            }
         }
      } else {
        lrs = log(rS);
        // ref1, (C5, unpol)
        epsXc = A * lrs + B + C * rS * lrs + D * rS;
        if (mode & ComputePotential)  {
           vXc(0) = A * lrs + (B - 1./3.*A) + 2./3.*C * rS * lrs
                  +  1./3. * (2. * D - C) * rS;  // ref1, (C6, unpol)
           SX_CHECK_NUM(vXc(0));
        }
        if (mode & ComputeKernel)  {
           fXc(0)(0) = - A / (3. * rho) 
                     - 2./9. * C * rS * (lrs + 1.) / rho
                     - (2. * D - C) / 9. * rS / rho;
           SX_CHECK_NUM(fXc(0)(0));
        }
      }
   } else {
      epsXc = 0.;
      vXc(0) = 0.;
      if (mode & ComputeKernel) fXc(0)(0) = 0.;
   }

   if (mode & ComputeX)  {
      // --- exchange part
      PrecRhoR ax = -3. / FOUR_PI * pow(9./4.* PI, 1./3.);
      double epsX = ax / rS;
      if (mode & ComputePotential)  {
         vXc(0) += 4./3. * epsX;
         SX_CHECK_NUM(vXc(0));
      }
      if (mode & ComputeKernel && rho > 0.)  {
         fXc(0)(0) += 4./9. * epsX / rho;
         SX_CHECK_NUM(fXc(0)(0));
      }
      epsXc += epsX;
   }
   
   eXc = epsXc * rho;
   SX_CHECK_NUM(eXc);
}

void SxXCFunctional::test (double rhoMin, double rhoMax, double factor)
{
   SX_CHECK (factor > 1., factor);
   for (double r = rhoMin; r <= rhoMax; r *= factor)  {
      computeLDA (r);
      sxprintf ("%.16f\t%.16f", r, getRs(r));
      if (mode & ComputeEnergy)
         sxprintf("\t%.16f",eXc);
      if (mode & ComputePotential)
         sxprintf("\t%.16f",vXc(0));
      if (mode & ComputeKernel)
         sxprintf("\t%.16f",fXc(0)(0));
      cout << endl;
   }
}

//==============================================================================
// LSDA - Perdew/Zunger '81
//==============================================================================
void SxXCFunctional::computeLSDA (double rhoUp, double rhoDown)
{
   SX_CHECK (nSpin == 2, nSpin);
   SX_CHECK (!(mode & ComputeKernel)); // not available, yet
   static const double EPS = 1e-7;
   static const double A_u  = 0.0311,  B_u  = -0.048;  // ref1, (C7)
   static const double A_p  = 0.01555, B_p  = -0.0269; // ref1, (C11)
   // --- ref1, (table xii)
   static const double g_u  = -0.1423, g_p  = -0.0843;
   static const double b1_u =  1.0529, b1_p =  1.3981;
   static const double b2_u =  0.3334, b2_p =  0.2611;
   static const double C_u  =  0.0020, C_p  =  0.0007;
   static const double D_u  = -0.0116, D_p  = -0.0048;

   static const double xn  = pow (2., 4./3.) - 2.;  // denom. of ref1, (C15)
   static const double axUnpol = -3./FOUR_PI * pow (9./4.*PI, 1./3.);
   static const double axPol   = axUnpol    * pow (2., 1./3.);
   
   // calculate r_s
   double rho = rhoUp + rhoDown;
   double rS = getRs (rho);

   // calculate spin polarization
   double zeta    = getZeta (rhoUp, rhoDown),
          absZeta = fabs(zeta);

   // --- compute f and f'=df/d(zeta)
   double f, fPrime;
   if (absZeta < EPS)  {
     f = fPrime = 0.;
   } else if (absZeta < (1. - EPS) ) {
     f      = ( pow (1.+ zeta, 4./3.) + pow (1. - zeta, 4./3.) - 2. )
             / xn;                                          // ref1, (C15)
     fPrime = ( pow (1.+ zeta, 1./3.) - pow (1. - zeta, 1./3.) ) * 4./3. 
             / xn;                                          // ref5, (A4)
   } else  {
     f      = 1.;
     fPrime = pow (2., 1./3.) * 4./3. / xn;
   }

   double eXcUnpol, eXcPol, vXcUnpol = 0., vXcPol = 0.;

   if (mode & ComputeC)  {
      // --- correlation part
      double lrs, rrs, bprPol, bprUnpol;
      if ( rS > 1. )  {
         rrs       = sqrt (rS);
         bprUnpol  = 1. + b1_u * rrs + b2_u * rS;
         eXcUnpol  = g_u / bprUnpol;          // ref1, (C3, unpol)
         bprPol    = 1. + b1_p * rrs + b2_p * rS;
         eXcPol    = g_p / bprPol;            // ref1, (C3, pol)
         if (mode & ComputePotential)  {
            vXcUnpol  = (1. + 7./6.*b1_u * rrs + 4./3.* b2_u * rS)
                      *  eXcUnpol / bprUnpol;   // ref1, (C4, unpol)
            vXcPol    = (1. + 7./6.*b1_p * rrs + 4./3.* b2_p * rS)
                      *  eXcPol   / bprPol;     // ref1, (C4, unpol)
         }
      } else  {
         lrs       = log  (rS);
         // ref1, (C5, unpol)
         eXcUnpol  = A_u * lrs + B_u + C_u * rS * lrs + D_u * rS; 
         // ref1, (C5, pol)
         eXcPol    = A_p * lrs + B_p + C_p * rS * lrs + D_p * rS; 
         if (mode & ComputePotential)  {
            vXcUnpol  = A_u * lrs + (B_u - 1./3.*A_u) + 2./3.*C_u * rS * lrs
                     +  1./3. * (2. * D_u - C_u) * rS;  // ref1, (C6, unpol)
            vXcPol    = A_p * lrs + (B_p - 1./3.*A_p) + 2./3.*C_p * rS * lrs
                     +  1./3. * (2. * D_p - C_p) * rS;  // ref1, (C6, pol)
         }
      }
   } else {
      eXcUnpol = eXcPol = 0.;
      if (mode & ComputePotential) vXcUnpol = vXcPol = 0.;
   }

   if (mode & ComputeX)  {
      // --- exchange contribution
      eXcUnpol += axUnpol / rS;
      eXcPol   += axPol   / rS;
      if (mode & ComputePotential)  {
         vXcUnpol += 4./3. * axUnpol / rS;
         vXcPol   += 4./3. * axPol   / rS;
      }
   }

   epsXc = eXcUnpol + f * (eXcPol - eXcUnpol); // ref1, (C12)
   eXc = rho * epsXc; // ref1, (C12)

   if (mode & ComputePotential)  {
      double vAvg, dVxc;
      vAvg = vXcUnpol                             // ref1, (C13)
             + f    *          (vXcPol - vXcUnpol) 
             - zeta * fPrime * (eXcPol - eXcUnpol);
      dVxc = fPrime * (eXcPol - eXcUnpol);
      vXc(0) = vAvg + dVxc; // ref1, (C13)
      vXc(1) = vAvg - dVxc; // ref1, (C13)
   }
}


//==============================================================================
// PBE - Perdew, Burke, Enzerhof
//==============================================================================
void SxXCFunctional::computePBE (double rhoUp, double gradRhoNormUp,
                                 double rhoDown, double gradRhoNormDown,
                                 double gradRhoNormTl)
{
   SX_CHECK (!(mode & ComputeKernel)); // not available, yet

   if (eXc_grad.getSize () != nSpin && (mode & ComputePotential))
      eXc_grad.resize (nSpin);
   if (mode & ComputeC)  {
      // --- correlation part
      int saveMode = mode;
      mode &= ~ComputeX; // switch off exchange
      // --- compute LSDA
      if (nSpin == 1)  {
         double rhoHalf = 0.5 * rhoUp;
         computeLSDA_PW (rhoHalf, rhoHalf);
      } else {
         computeLSDA_PW (rhoUp, rhoDown);
      }
      // restore mode
      mode = saveMode;

      // gradient corrections
      static const double ETA = 1e-14;
      double rho  = (nSpin == 1) ? rhoUp : (rhoUp + rhoDown);
      double rs = rsTot; // computed in computeLSDA_PW
      double zeta, phi;
      if  (nSpin == 1) {
         zeta = 0.;
         phi = 1.;
      } else {
         zeta =  getZeta (rhoUp,rhoDown);
         phi = 0.5 * (cbrt(sqr(1. + zeta)) + cbrt(sqr(1. - zeta))); // ref2,(20)
      }
      //double ks   = pow ((192.* fabs(rho)/PI), (1./6.)); //ref 10, (A11b)
      double ks, t;
      if (rho <= 1e-12)  {
         ks = 1e-16;
         t = 0.;
      } else {
         ks = sqrt(cbrt((192.* fabs(rho)/PI))); //ref 10, (A11b)
         t = gradRhoNormTl / (2. * ks * phi * rho);     //ref 10, (A10c)
      }
      double H, eC = epsXc;
      static const double gamma = 0.031090690869654895035;   // ref6,(5f)
      static const double beta  = 0.06672455060314922;       // ref6,(4f)
      double A, A2, phi3, phi4, t2, t4, t6;
      phi3 = phi * phi * phi;
      phi4 = phi * phi3;
      t2   = t * t; 
      t4   = t2 * t2;
      t6   = t4 * t2;
      double expEcUnif = exp (-eC / (gamma*phi3));
      if (eC < 0.)
         A  = beta/gamma / (expEcUnif - 1.);           // ref6,(A8)
      else
         A = 0.;
      A2 = A * A;
      H  = gamma * phi3 * log (1. + beta/gamma * t2 * (1. + A*t2)
                                    / (1. + A*t2 + A2*t4) );      // ref6,(7)

      // gradient correction to correlation energy density
      epsXc += H;
      // recompute eXc
      eXc = epsXc * rho;

      // --- gradient correction to vXc

      // denominator for Eq. A.51,57, ref. 10
      double denom = (1. + A*t2 + A2*t4)
                   * (gamma + beta*t2 + A*gamma*t2 + A*beta*t4 + A2*gamma*t4);
      double denomI = 1. / denom;
      
      double A_phi, A_ecUnif, H_A, H_rs, H_t;
      double A_part = A2 * expEcUnif / (beta * phi4);
      // dA/d phi     Eq. A54 ref. 10
      //A_phi      = -3. * A2 * eC * expEcUnif / (beta * phi4);
      A_phi      = -3. * eC * A_part;
      // dA/d ecUnif  Eq. A52 ref. 10
      //A_ecUnif   = A2 * expEcUnif / (beta * phi3);
      A_ecUnif   = A_part * phi;
      // dH/dA        Eq. A51 ref. 10
      //H_A        = -beta * gamma * phi3 * A * t6 * (2.  +  A*t2) / denom;
      H_A        = -beta * gamma * phi3 * A * t6 * (2. +  A*t2) * denomI;
      // dH/drs       Eq. A50 ref. 10
      H_rs       = H_A * A_ecUnif * ec_rs;
      // dH/d t  Eq. A57 ref. 10
      //H_t        = 2. * beta * gamma * phi3 * t * (1.  +  2.*A*t2) / denom;
      H_t        = 2. * beta * gamma * phi3 * t * (1.  +  2.*A*t2) * denomI;

      for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
         vXc(iSpin) += H - rs * (1./3.) * H_rs - (7. / 6.) * t * H_t;
      }
      if (nSpin == 2)  {
         // d phi/d zeta
         double phi_z  = (  pow( sqr(1.+zeta) + ETA, -1./6.)
                         -  pow( sqr(1.-zeta) + ETA, -1./6.)
                          ) * (1./3.);
         // dH/d zeta    Eq. A53 ref. 10
         double H_z = H_A * (A_phi * phi_z + A_ecUnif * ec_z)
                    + 3. * H/phi * phi_z; // A56 ref. 10
         double zetaDeriv = H_z - H_t * phi_z * t/phi;
         vXc(0) += (1. - zeta) * zetaDeriv;
         vXc(1) -= (1. + zeta) * zetaDeriv;
      }

      // gradient with respect to gradRhoNormTl
      eXc_gradRhoTl =  H_t/(2.*phi*ks);  // ref9 (D.45) 
      SX_CHECK_NUM (eXc_gradRhoTl);

   } else {
      epsXc = 0.;
      eXc = 0.;
      vXc.set(0.);
      eXc_gradRhoTl = 0.;
      if (nSpin == 1) rsTot = getRs (rhoUp);
   }

   if (mode & ComputeX)  {
      // --- exchange part
      for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
         double rho = (iSpin == 0) ? rhoUp : rhoDown;
         if (rho <= 1e-50)  {
            eXc_grad(iSpin) = 0.;
            continue;
         }
         double rhoX, rS;
         if (nSpin == 1)  {
            rhoX = rho;
            rS = rsTot;
         } else  {
            rhoX = 2. * rho; // exchange formula for up+down
            rS = getRs(rhoX);
         }
         double kf = getKf(rS);
         double s = ((iSpin == 0) ? gradRhoNormUp : gradRhoNormDown)
                  / (2. * kf * rho); // ref2, (14)

         //static const double ax  = -0.738558766382022; // -3/4 (3/pi)^{1/3}
         static const double axrs = -0.75 * cbrt(2.25 / (PI * PI));
         static const double mu    = 0.2195149727645171; // ref1;ref2,(15f)
         static const double kappa = 0.8040;             // ref2,(15f)
         double exUnif, epsX, fX;
         //exUnif = ax * pow (fabs(rhoX), 1./3.);        // ref2, (12)
         exUnif = axrs/rS;                               // ref2, (12)
         double denomInv = kappa / (kappa + mu * s*s);
         fX = 1. + kappa - kappa * denomInv; // ref2, (15)
         epsX  = exUnif * fX;                            // ref2, (11)
         epsXc += epsX;
         eXc += epsX * rho;
         if (mode & ComputePotential)  {
            double dfds = 2. * mu * s * denomInv * denomInv;
            // d eXc / d rho
            vXc(iSpin) += (4. / 3.) * epsX // f * d (exUnif * rho) / d rho
                       // exUnif * rho * df/ds * ds/drho
                       -  exUnif * dfds * s * (4./3.); // exUnif * rho * df/drho
            SX_CHECK_NUM(vXc(iSpin));
            // d eXc / d gradNorm
            eXc_grad(iSpin) = exUnif * dfds / (2. * kf);
            SX_CHECK_NUM(eXc_grad(iSpin));
         }
      }
   } else {
      if (mode & ComputePotential)
         eXc_grad.set (0.);
   }
}

void SxXCFunctional::addScreenedExPBE (double rho, double gradRho,
                                       int iSpin, double scaling)
{
   SX_CHECK (!(mode & ComputeKernel)); // not available, yet
   SX_CHECK (iSpin >= 0 && iSpin < nSpin, iSpin, nSpin);
   if (! (mode & ComputeX)) return;

   if (eXc_grad.getSize () != nSpin && (mode & ComputePotential))  {
      eXc_grad.resize (nSpin);
      eXc_grad(iSpin) = 0.;
   }

   if (!xHolePtr) xHolePtr = SxPtr<SxPBExHole>::create ();

   // --- exchange part
   if (rho <= 1e-18) return; 
   double rhoX = (nSpin == 2) ? 2. * rho : rho; // exchange formula for up+down
   double rS = getRs(rhoX);
   double kf = getKf(rS);
   double s = gradRho / (2. * kf * rho); // ref11, Eq. 4a
   static const double sMax = 8.572844;

   double sBar;
   if (s < 50.)
      sBar = s - (1. - exp(-s)) * log(1. + exp(s-sMax)); // ref11, Eq. 9
   else
      sBar = sMax;
   double nu = omegaHSE / kf; // ref 11, below Eq. (7)
   xHolePtr->screenedFx (sBar, nu);

   double epsX = scaling * (-0.75 / PI) * kf * xHolePtr->Fx; // ref 11, Eq. 5
   epsXc += epsX;
   eXc   += epsX * rho; // ref 11, Eq. 5

   if (mode & ComputePotential)  {
      double ds    = -4. / 3. * s,  // rho *     ds/(d rho)
             dNu   = -1. / 3. * nu; // rho * (d nu)/(d rho)
      double dsBar;
      if (s < 50.)  {
         double es = exp(-s), esm = exp(s-sMax);
         dsBar = 1. - es * log (1. + esm) - (1.-es) * esm/(1. + esm);
      } else {
         dsBar = 0.;
      }
      // d eXc / d rho
      vXc(iSpin) += 4. / 3. * epsX
                 +  scaling * (-0.75 / PI) * kf
                    * (xHolePtr->Fxs * dsBar * ds + xHolePtr->FxNu * dNu);
      SX_CHECK_NUM(vXc(iSpin));
      // d eXc / d gradNorm
      eXc_grad(iSpin) += 0.5 * scaling * (-0.75 / PI) * xHolePtr->Fxs * dsBar;
      SX_CHECK_NUM(eXc_grad(iSpin));
   }
}

void SxXCFunctional::computeHSE (double rhoUp, double gradRhoNormUp,
                                 double rhoDown, double gradRhoNormDown,
                                 double gradRhoNormTl)
{
   computePBE (rhoUp, gradRhoNormUp, rhoDown, gradRhoNormDown, gradRhoNormTl);
   addScreenedExPBE (rhoUp, gradRhoNormUp, 0, -alphaHybrid);
   if (nSpin == 2)
      addScreenedExPBE (rhoDown, gradRhoNormDown, 1, -alphaHybrid);
}

void SxXCFunctional::computePBEHybrid (double rhoUp, double gradRhoNormUp,
                                       double rhoDown, double gradRhoNormDown,
                                       double gradNorm, double alphaExact)
{
   SX_CHECK (!(mode & ComputeKernel)); // not available, yet

   int origMode = mode;
   // switch off correlation
   mode &= ~ComputeC;
   computePBE (rhoUp, gradRhoNormUp, rhoDown, gradRhoNormDown, gradNorm);

   // --- keep exchange values
   double epsX = epsXc, eX = eXc;
   SX_CHECK (nSpin == 1 || nSpin == 2, nSpin);
   double vX[2], eX_grad[2];
   for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
      vX[iSpin] = vXc(iSpin);
      eX_grad[iSpin] = eXc_grad(iSpin);
   }

   // switch off exchange
   mode = origMode & (~ComputeX);
   computePBE (rhoUp, gradRhoNormUp, rhoDown, gradRhoNormDown, gradNorm);

   // --- add scaled exchange values
   double scale = 1. - alphaExact;
   epsXc += scale * epsX;
   eXc   += scale * eX;
   for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
      vXc(iSpin)      += scale * vX[iSpin];
      eXc_grad(iSpin) += scale * eX_grad[iSpin];
   }
   
   mode = origMode;
}

//==============================================================================
// LSDA - Perdew/Wang '91
//==============================================================================
void SxXCFunctional::computeLSDA_PW (double rhoUp, double rhoDown)
{
   SX_CHECK (!(mode & ComputeKernel)); // not available, yet

   // check that rhoUp = rhoDown (=0.5 rhoTotal) for unpolarized case.
   // rhoUp and rhoDown must be exactly the same, thus no fabs(rhoUp-rhoDown)
   SX_CHECK (nSpin == 2 || rhoUp == rhoDown, nSpin, rhoUp, rhoDown); 

   static const double xn      = pow (2., 4./3.) - 2.;  // denom. of ref1, (C15)
   static const double xnI     = 1. / xn;
   static const double axUnpol = -3./FOUR_PI * pow (9./4.* PI, 1./3.);
   static const double axPol   = axUnpol    * pow (2., 1./3.); 
   static const double EPS = 1e-7;
   
   // calculate r_s
   double rho = rhoUp + rhoDown;
   double rS = rsTot = getRs (rho);

   // calculate spin polarization
   double zeta    = getZeta (rhoUp, rhoDown),
          absZeta = fabs(zeta);

   // --- compute f and f'=df/d(zeta)
   double f, fz;
   if (absZeta < EPS)  {
     f = fz = 0.;
   } else if (absZeta < (1. - EPS) ) {
     /*
     f  = ( pow (1.+ zeta, 4./3.) + pow (1. - zeta, 4./3.) - 2. )
          / xn;                                          // ref1, (C15); ref5, (9)
     fz = ( pow (1.+ zeta, 1./3.) - pow (1. - zeta, 1./3.) ) * 4./3. 
          / xn;                                          // ref5, (A4)
     */
     double zp1 = 1. + zeta;
     double zm1 = 1. - zeta;
     double zp1_3 = cbrt(zp1);
     double zm1_3 = cbrt(zm1);
     f = (zp1 * zp1_3 + zm1 * zm1_3 - 2.) * xnI;
     fz = (zp1_3 - zm1_3) * (4./3.) * xnI;
   } else  {
     f      = 1.;
     fz = cbrt(2.) * (4./3.) * xnI;
   }

   if (mode & ComputeC)  {
      // --- correlation part
      // --- ref5, (tab 1, col5-unpolar,  col6-polarized,    col7-alpha_c(rs))
      static const double AU  = 0.0310907, AP  = 0.01554535,  Aa  = 0.0168869;
      static const double a1U = 0.21370,   a1P = 0.20548,     a1a = 0.11125;
      static const double b1U = 7.5957,    b1P = 14.1189,     b1a = 10.357;
      static const double b2U = 3.5876,    b2P = 6.1977,      b2a = 3.6231;
      static const double b3U = 1.6382,    b3P = 3.3662,      b3a = 0.88026;
      static const double b4U = 0.49294,   b4P = 0.62517,     b4a = 0.49671;
      
      double eUnpol, ePol, alf, alf_rs, eu_rs, ep_rs;
      // double ec_rs, ec_z; // => class members for PBE functional
      double zeta_3 = zeta*zeta*zeta;
      double zeta_4 = zeta_3 * zeta;

      static const double fzz = 1.7099209341613656176;  // df'/dz, ref5,(A4)|z->0
      
      // --- get e_c_{pol,unpol} and spin stiffness
      alf    = G_PW91 (Aa, a1a, b1a, b2a, b3a, b4a, rS, &alf_rs); 
      alf    = -alf; alf_rs = -alf_rs; // flip sign
      eUnpol = G_PW91 (AU, a1U, b1U, b2U, b3U, b4U, rS, &eu_rs); 
      if (nSpin == 1)
         ePol = ep_rs = 0.;
      else
         ePol = G_PW91 (AP, a1P, b1P, b2P, b3P, b4P, rS, &ep_rs);

      epsXc = eUnpol + alf*f/fzz*(1.-zeta_4)
            + (ePol - eUnpol)*f*zeta_4;                         // ref5,(8)

      if (mode & ComputePotential)  {
         // --- LSD potential
         ec_rs = eu_rs*(1. - f*zeta_4) + ep_rs*f*zeta_4
               + alf_rs * f/fzz * (1. - zeta_4);                // ref5,(A2)
         ec_z  = 4. * zeta_3*f  * ( ePol - eUnpol - alf/fzz)
               +             fz * (zeta_4*ePol - zeta_4*eUnpol
               +                    (1. - zeta_4) * alf/fzz);   // ref5,(A3)

         vXc(0) = epsXc - rS*ec_rs/3. - zeta*ec_z + ec_z;         // ref5,(A1)
         if (nSpin != 1)
            vXc(1) = epsXc - rS*ec_rs/3. - zeta*ec_z - ec_z;    // ref5,(A1)
      }
   } else {
      epsXc = 0.;
      if (mode & ComputePotential) vXc.set (0.);
   }
   
   
   if (mode & ComputeX)  {
      // --- exchange contribution (from Perdew-Zunger '81)
      //SX_EXIT;
      double eXcUnpol, eXcPol, vXcUnpol = 0., vXcPol = 0.;
      eXcUnpol = axUnpol / rS;
      eXcPol   = axPol   / rS;
      if (mode & ComputePotential)  {
         vXcUnpol = 4./3. * axUnpol / rS;
         vXcPol   = 4./3. * axPol   / rS;
      }

      epsXc += eXcUnpol + f * (eXcPol - eXcUnpol); // ref1, (C12)

      if (mode & ComputePotential)  {
         if (nSpin == 1) {
            vXc(0) += vXcUnpol; 
         } else {
            double vAvg, dVxc;
            vAvg = vXcUnpol                             // ref1, (C13)
                + f    *      (vXcPol - vXcUnpol) 
                - zeta * fz * (eXcPol - eXcUnpol);
            dVxc = fz * (eXcPol - eXcUnpol);
            vXc(0) += vAvg + dVxc; // ref1, (C13)
            vXc(1) += vAvg - dVxc; // ref1, (C13)
         }
      }
   }
   
   eXc = rho * epsXc; // ref1, (C12)
}

double SxXCFunctional::G_PW91 (double A, double a1, double b1,  double b2, 
                              double b3, double b4, double rs, double *grs)
{
   double sqrt_rs = sqrt(rs);
   double q0, q1, q1Prime, q2;
   q0      = -2. * A * ( 1. + a1 * rs);                     // ref5,(A6)
   q1      =  2. * A * ( b1*sqrt_rs + b2*rs 
                       + b3*rs*sqrt_rs + b4*rs*rs);         // ref5,(A7), p=1
   q1Prime = A * ( b1 / sqrt_rs + 2.*b2 + 3.*b3*sqrt_rs
                 + 4.*b4*rs);                               // ref5,(A8), p=1
   *grs    = -2. * A * a1 * (q2=log(1.+1./q1)) 
           - q0*q1Prime / (q1*q1 + q1);                     // ref5,(A5)
   return    q0 * q2;
}


void SxXCFunctional::testPW91 ()
{
   // ---- compare result of G with fhi98md's gcor2 subroutine
   double rs, g, eu_rs;
   SxFileIO fp;
   try {
      fp.open ("pw91.dat", "a");
      for (rs=0.1; rs < 5.0; rs+=0.1)  {
         g = G_PW91 (0.0310907,0.21370,7.5957,3.5876,1.6382, 0.49294,rs,&eu_rs);
         fp.write(SxString(rs) + "\t" + g + "\t" + eu_rs + "\n");
      }
      fp.close ();
   } catch (const SxException &e)  {
      e.print ();
      SX_EXIT;
   }
}

