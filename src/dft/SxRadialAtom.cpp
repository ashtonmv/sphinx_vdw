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

#include <SxRadialAtom.h>
#include <SxRadBasis.h>
#include <SxXCFunctional.h>
#include <SxDifferential.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace Timer {
   enum RadialAtomTimer { RadXC, XcSetup, Rho, GetGrad, Grad, Xc, XcG1,
                          XcG2, XcFunc };
}

SX_REGISTER_TIMERS(Timer::RadialAtomTimer)
{
   using namespace Timer;
   regTimer (RadXC, "rad. xc");
   regTimer (XcSetup, "agrid setup");
   regTimer (Rho, "rho/agrid");
   regTimer (GetGrad, "GetGrad");
   regTimer (Grad, "Grad");
   regTimer (Xc, "Xc");
   regTimer (XcG1, "XcG1");
   regTimer (XcG2, "XcG2");
   regTimer (XcFunc, "XcFunc");
}

SxRadialAtom::SxRadialAtom (SxXC::XCFunctional xcIn,
                            SxSphereGrid::GridType gridTypeIn)
   : eXc(-1.),
     gridType (gridTypeIn),
     xcFunctional(xcIn)
{
   // empty
}

SxDiracVec<Double> 
SxRadialAtom::getHartreePotential (const SxDiracVec<Double> &rho)
{
   const SxRadBasis *radPtr 
      = dynamic_cast<const SxRadBasis *>(rho.getBasisPtr ());
   SX_CHECK (radPtr);
   int is = rho.handle->auxData.is;
   int l = rho.handle->auxData.l;
   const SxDiracVec<Double> &radFunc = radPtr->radFunc (is);

   double prefac = FOUR_PI/(2*l+1) * radPtr->logDr(is);

   int ir, nr = (int)rho.getSize ();
   SxDiracVec<Double> res(nr);
   res.handle->auxData = rho.handle->auxData;

   // forward integration (r' < r)
   // 4pi / ((2l+1)r^(l+1)) \int_0^r (r')^3 d(ln r') (r')^(l) rho(r')
   double sum = 0., r, rl1;
   int nend = (nr % 2 == 1) ? nr -1 : nr - 4;

   res(0) = 0.;
   sum    = pow(radFunc(0), l + 3.) * rho(0) / 3.;
   // --- Simpson integration 1/3 + 4/3 + 1/3
   for (ir = 1; ir < nend; ++ir)  {
      r = radFunc(ir);
      rl1 = pow(r, double(l+1));
      res(ir) = prefac * sum / rl1;
      sum += (ir & 1 ? 4. : 2.) / 3. * rl1 * r * r * rho(ir);
   }
   if (nr % 2 == 0)  {
      ir = nr - 4; // 1/3 + 3/8 = 17/24 
      r = radFunc(ir);
      rl1 = pow(r, double(l+1));
      res(ir) = prefac * sum / rl1;
      sum +=  17. / 24. * rl1 * r * r * rho(ir);

      ir = nr - 3; // 9/8
      r = radFunc(ir);
      rl1 = pow(r, double(l+1));
      res(ir) = prefac * sum / rl1;
      sum +=  9. / 8. * rl1 * r * r * rho(ir);

      ir = nr - 2; // 9/8
      r = radFunc(ir);
      rl1 = pow(r, double(l+1));
      res(ir) = prefac * sum / rl1;
      sum +=  9. / 8. * rl1 * r * r * rho(ir);

      ir = nr - 1; // 3/8
      r = radFunc(ir);
      rl1 = pow(r, double(l+1));
      res(ir) = prefac * sum / rl1;
      sum +=  3. / 8. * rl1 * r * r * rho(ir);
   } else {
      ir = nr - 1; // 1/3
      r = radFunc(ir);
      rl1 = pow(r, double(l+1));
      res(ir) = prefac * sum / rl1;
      sum +=  1. / 3. * rl1 * r * r * rho(ir);
   }
   VALIDATE_VECTOR(res);

   // backward integration (r' >= r)
   // 4pi r^l / ((2l+1)) \int_r^\infty (r')^3 d(ln r') 1/(r')^(l+1) rho(r')
   double rl;

   sum = 0.;
   if (nr % 2 == 0)  {
      ir = nr - 1; // 3/8
      r = radFunc(ir);
      rl = pow(r, double(l));
      sum += 3. / 8. * r * r / rl * rho(ir);
      res(ir) += prefac * sum * rl;

      ir = nr - 2; // 9/8
      r = radFunc(ir);
      rl = pow(r, double(l));
      sum += 9. / 8. * r * r / rl * rho(ir);
      res(ir) += prefac * sum * rl;

      ir = nr - 3; // 9/8
      r = radFunc(ir);
      rl = pow(r, double(l));
      sum += 9. / 8. * r * r / rl * rho(ir);
      res(ir) += prefac * sum * rl;

      ir = nr - 4; // 1/3 + 3/8 = 17/24 
      r = radFunc(ir);
      rl = pow(r, double(l));
      sum += 17. / 24. * r * r / rl * rho(ir);
      res(ir) += prefac * sum * rl;
   } else {
      ir = nr - 1; // 1/3
      r = radFunc(ir);
      rl = pow(r, double(l));
      sum += 1. / 3. * r * r / rl * rho(ir);
      res(ir) += prefac * sum * rl;
   }

   // --- Simpson integration 1/3 + 4/3 + 1/3
   for (--ir; ir > 0; --ir)  {
      // --- odd
      r = radFunc(ir);
      rl = pow(r, double(l));
      sum += (ir & 1 ? 4. : 2.) / 3. * r * r / rl * rho(ir);
      res(ir) += prefac * sum * rl;
   }
   // ir = 0 (1/3)
   r = radFunc(ir);
   rl = pow(r, double(l));
   sum += 1. / 3. * r * r / rl * rho(ir);
   res(ir) += prefac * sum * rl;

   VALIDATE_VECTOR(res);
   return res;
}

SxDiracVec<Double> 
SxRadialAtom::getHartreeSmooth(const SxDiracVec<Double> &rho)
{
   cout << "smooth" << endl;
   const SxRadBasis *radPtr 
      = dynamic_cast<const SxRadBasis *>(rho.getBasisPtr ());
   SX_CHECK (radPtr);
   int is = rho.handle->auxData.is;
   int l = rho.handle->auxData.l;
   const SxDiracVec<Double> &radFunc = radPtr->radFunc (is);

   double prefac = FOUR_PI/(2*l+1) * radPtr->logDr(is);

   int ir, nr = (int)rho.getSize ();
   SxDiracVec<Double> res(nr);
   res.handle->auxData = rho.handle->auxData;

   // forward integration (r' < r)
   // 4pi / ((2l+1)r^(l+1)) \int_0^r (r')^3 d(ln r') (r')^(l) rho(r')
   double sum = 0., r, rl1;
   int nend = (nr % 2 == 1) ? nr -1 : nr - 4;

   res(0) = 0.;
   sum    = pow(radFunc(0), l + 3.) * rho(0) / 3.;
   // --- Simpson integration 1/3 + 4/3 + 1/3
   double dVOld = 3. * sum;
   for (ir = 1; ir < nend; ++ir)  {
      r = radFunc(ir);
      rl1 = pow(r, double(l+1));
      double dV = rl1 * r * r * rho(ir);
      res(ir) = prefac * (sum + (ir & 1 ? (dVOld/6. + dV/2.) : (dV/3.))) / rl1;
      sum += (ir & 1 ? 4. : 2.) / 3. * dV;
      dVOld = dV;
   }
   if (nr % 2 == 0)  {
      ir = nr - 4; // 1/3 + 3/8 = 17/24 
      r = radFunc(ir);
      rl1 = pow(r, double(l+1));
      res(ir) = prefac * sum / rl1;
      sum +=  17. / 24. * rl1 * r * r * rho(ir);

      ir = nr - 3; // 9/8
      r = radFunc(ir);
      rl1 = pow(r, double(l+1));
      res(ir) = prefac * sum / rl1;
      sum +=  9. / 8. * rl1 * r * r * rho(ir);

      ir = nr - 2; // 9/8
      r = radFunc(ir);
      rl1 = pow(r, double(l+1));
      res(ir) = prefac * sum / rl1;
      sum +=  9. / 8. * rl1 * r * r * rho(ir);

      ir = nr - 1; // 3/8
      r = radFunc(ir);
      rl1 = pow(r, double(l+1));
      res(ir) = prefac * sum / rl1;
      sum +=  3. / 8. * rl1 * r * r * rho(ir);
   } else {
      ir = nr - 1; // 1/3
      r = radFunc(ir);
      rl1 = pow(r, double(l+1));
      res(ir) = prefac * sum / rl1;
      sum +=  1. / 3. * rl1 * r * r * rho(ir);
   }
   VALIDATE_VECTOR(res);

   // backward integration (r' >= r)
   // 4pi r^l / ((2l+1)) \int_r^\infty (r')^3 d(ln r') 1/(r')^(l+1) rho(r')
   double rl;

   sum = 0.;
   if (nr % 2 == 0)  {
      ir = nr - 1; // 3/8
      r = radFunc(ir);
      rl = pow(r, double(l));
      sum += 3. / 8. * r * r / rl * rho(ir);
      res(ir) += prefac * sum * rl;

      ir = nr - 2; // 9/8
      r = radFunc(ir);
      rl = pow(r, double(l));
      sum += 9. / 8. * r * r / rl * rho(ir);
      res(ir) += prefac * sum * rl;

      ir = nr - 3; // 9/8
      r = radFunc(ir);
      rl = pow(r, double(l));
      sum += 9. / 8. * r * r / rl * rho(ir);
      res(ir) += prefac * sum * rl;

      ir = nr - 4; // 1/3 + 3/8 = 17/24 
      r = radFunc(ir);
      rl = pow(r, double(l));
      sum += 17. / 24. * r * r / rl * rho(ir);
      res(ir) += prefac * sum * rl;
   } else {
      ir = nr - 1; // 1/3
      r = radFunc(ir);
      rl = pow(r, double(l));
      sum += 1. / 3. * r * r / rl * rho(ir);
      res(ir) += prefac * sum * rl;
   }

   // --- Simpson integration 1/3 + 4/3 + 1/3
   dVOld = 0.;
   for (--ir; ir > 0; --ir)  {
      // --- odd
      r = radFunc(ir);
      rl = pow(r, double(l));
      double dV = r * r / rl * rho(ir);
      sum += (ir & 1 ? 4. : 2.) / 3. * dV;
      res(ir) += prefac * (sum - (ir & 1 ? (5./6.*dV - dVOld/6.) : (dV/3.))) * rl;
      //res(ir) += prefac * sum * rl;
      dVOld = dV;
   }
   // ir = 0 (1/3)
   r = radFunc(ir);
   rl = pow(r, double(l));
   sum += 1. / 3. * r * r / rl * rho(ir);
   res(ir) += prefac * sum * rl;

   VALIDATE_VECTOR(res);
   return res;
}

double
SxRadialAtom::computeScreenedHartree (const SxDiracVec<Double> &rho1,
                                      const SxDiracVec<Double> &rho2,
                                      double omega, double dg)
{
   SX_CHECK (dg > 0.);
   SX_CHECK (rho1.getBasisPtr () == rho2.getBasisPtr ());
   const SxRadBasis &rad = rho1.getBasis<SxRadBasis> ();
   SX_CHECK (rho1.handle);
   SX_CHECK (rho2.handle);
   SX_CHECK (rho2.handle->auxData.is == rho1.handle->auxData.is,
             rho2.handle->auxData.is, rho1.handle->auxData.is);
   SX_CHECK (rho2.handle->auxData.l == rho1.handle->auxData.l,
             rho2.handle->auxData.l, rho1.handle->auxData.l);
   SX_CHECK (rho2.handle->auxData.m == rho1.handle->auxData.m,
             rho2.handle->auxData.m, rho1.handle->auxData.m);
   
   int is = rho1.handle->auxData.is,
       l  = rho1.handle->auxData.l;
   SX_CHECK (is >= 0 && is < rad.getNSpecies (), is, rad.getNSpecies ());
   SX_CHECK (l >= 0, l);
   double sum = 0.;
   double logDr = rad.logDr(is);
   SxDiracVec<Double> r3 = rad.radFunc(is).cub (), jsbR3;
   if (omega > 1e-10)  {
      // --- screening part: reciprocal space integration
      double gMax = 13. * omega;
      int ng = int(gMax / dg);
      SX_CHECK (ng > 10 && ng < 1000000, ng);
      for (int ig = 0; ig < ng; ++ig)  {
         double g = ig * dg;
         jsbR3  = SxRadBasis::jsb(l, g * rad.radFunc(is));
         jsbR3 *= r3;

         // HSE implementation notes, Eq. 10
         sum += weightSimpson(ig, ng) 
              * (rho1 * jsbR3).integrate (logDr) // Eq. 9, no prefactor 
              * (rho2 * jsbR3).integrate (logDr) // Eq. 9, no prefactor
              * exp(-g*g/(4. * omega * omega));
      }
      sum *= 8. * dg; // prefactor Eq. 9/10
   }
   double bareHartree = (rho1 * getHartreePotential(rho2) * r3)
                        .integrate (logDr);
   return bareHartree - sum;
}

SxArray<SxRadialMesh>
SxRadialAtom::computeXC (const SxArray<SxRadialMesh> &rho)
{
   SX_CLOCK (Timer::RadXC);
   SX_START_TIMER (Timer::XcSetup);
   const SxRadBasis &radBasis = rho(0).getBasis ();
   const int is = rho(0).getIs ();
   const SxDiracVec<Double> &rad = radBasis.radFunc(is);
   const int nSpin = int(rho.getSize ());
   const int nr = rho(0).getNRad ();
   const int nLm = (int)rho(0).meshData.nCols ();
   const int lmax = (int) lround(sqrt(double(nLm))) - 1;
   SX_CHECK (nLm == sqr(lmax+1), nLm, lmax);
   const double dex = radBasis.logDr(is);

   eXc = 0.;
   SxArray<SxRadialMesh> vXc(nSpin), vXcGradR;
   for (int iSpin = 0; iSpin < nSpin; ++iSpin) {
      vXc(iSpin).resize (rho(0));
      vXc(iSpin).set (0.);
   }

   SxXCFunctional xc(nSpin, SxXCFunctional::ComputeExcVxc);
   SxDiracVec<Double> energy(radBasis(is));

   SxSphereGrid aGrid (gridType, lmax);
   int nOm = (int)aGrid.getSize ();

   bool gga = (   xcFunctional == SxXC::PBE
               || xcFunctional == SxXC::PBE0
               || xcFunctional == SxXC::HSE06);
   bool ldaPW = xcFunctional == SxXC::LDA_PW;
   SxMatrix<Double> YlmXyz, Psilm2;
   SxArray<SxVector3<Double> > Psilm3;
   if (gga)  {
      aGrid.setupYlmPsilm (lmax);
      vXcGradR  .resize (nSpin);
      for (int iSpin = 0; iSpin < nSpin; ++iSpin) {
         vXcGradR(iSpin).resize (rho(0));
         vXcGradR(iSpin).set (0.);
      }

      // --- set up vecY_lm and Psi_lm in various orderings
      YlmXyz = SxMatrix<Double> (nOm, 3 * nLm);
      for (int ilm = 0; ilm < nLm; ++ilm)
         for (int idir = 0; idir < 3; ++idir) 
            YlmXyz.colRef (idir + 3*ilm) <<= aGrid.Ylm.colRef(ilm)
                                             * aGrid.xyz.colRef(idir);
      Psilm2.copy (aGrid.Psilm); 
      Psilm2.reshape (nOm * 3, nLm);
      YlmXyz.reshape (nOm * 3, nLm);
      Psilm3.resize (nLm * nOm);
      for (int ilm = 0; ilm < nLm; ++ilm)
         for (int idir = 0; idir < 3; ++idir) 
            for (int iOm = 0; iOm < nOm; ++iOm)
               Psilm3(ilm + nLm * iOm)(idir) = Psilm2(iOm + nOm * idir, ilm); 

   }
   int numGradOrder = 4;
   SxDifferential radialGrad (numGradOrder);
   SxMatrix<Double> vXcGradRadAll;
   if (gga)  {
      vXcGradRadAll.reformat (nLm, nr * nSpin);
      vXcGradRadAll.set (0.);
   }
   SX_STOP_TIMER (Timer::XcSetup);

#  ifdef USE_OPENMP
#  pragma omp parallel firstprivate(xc)
#  endif
   {
      SX_ALLOC_CACHE;
      SxArray<SxMatrix<Double> > rhoAng(nSpin), vXcAng(nSpin);
      for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
         vXcAng(iSpin).resize (nOm);
      }
      SxVector<Double> vXcLm;

      SxVector<Double> dRhoDr(nLm);
      SxVector<Double> rhoLm(nLm);

      SxMatrix<Double> gradAng, vXcGradRad, vXcGradAng;
      if (gga)  {
         gradAng   .reformat (nOm, nSpin * 3);
         vXcGradRad.reformat (nOm, nSpin);
         vXcGradAng.reformat (nLm, nSpin);
      }
#     ifdef USE_OPENMP
#     pragma omp for schedule(dynamic,10)
#     endif
      for (int ir = 0; ir < nr; ++ir)  {

         // --- set up density on angular grid
         if (gga) {
            gradAng.set (0.);
            gradAng.reshape (nOm * 3, nSpin);
         }
         for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
            //SX_START_TIMER (Timer::Rho);
            for (int ilm = 0; ilm < nLm; ++ilm)
               rhoLm(ilm) = rho(iSpin).meshData(ir,ilm);
            rhoAng(iSpin) = aGrid.Ylm ^ rhoLm;
            //SX_STOP_TIMER (Timer::Rho);
            if (gga)  {
               //SX_CLOCK (Timer::Grad);
               // --- set up gradient
               double r = rad(ir) ;
               for (int ilm = 0; ilm < nLm; ++ilm)  {
                  dRhoDr(ilm)
                     = radialGrad.apply (rho(iSpin).meshData.colRef(ilm), ir)
                     / (dex * r);
               }
               // radial gradient
               gradAng.colRef(iSpin) += YlmXyz ^ dRhoDr;
               // angular gradient
               gradAng.colRef(iSpin) += Psilm2 ^ (rhoLm / r);
               vXcGradAng.set (0.);
            }
            VALIDATE_VECTOR (rhoAng(iSpin));
         }
         if (gga)  {
            gradAng.reshape (nOm, 3 * nSpin);
         }

         // compute xc energy/potential on spherical grid
         energy(ir) = 0.;
         //SX_START_TIMER (Timer::Xc);
         if (nSpin == 1)  {
            double eXc_ = 0.;
            for (int iOm = 0; iOm < nOm; ++iOm)  {
               if (!gga)  {
                  if (!ldaPW)   {
                     // --- LDA Perdew/Zunger '81
                     xc.computeLDA (rhoAng(0)(iOm));
                  } else   {
                     // --- LSDA Perdew/Wang '91
                     double rhoHalf = 0.5 * rhoAng(0)(iOm);
                     xc.computeLSDA_PW(rhoHalf,rhoHalf);
                  }
               } else {
                  // --- GGA functionals
                  SxVector3<Double> grad(gradAng(iOm, 0), gradAng(iOm,1),
                                         gradAng(iOm, 2));
                  double gNorm = grad.norm ();
                  //SX_START_TIMER (Timer::XcFunc);
                  if (xcFunctional == SxXC::PBE)
                     xc.computePBE (rhoAng(0)(iOm), gNorm, -1., -1., gNorm);
                  else if (xcFunctional == SxXC::PBE0)
                     xc.computePBEHybrid (rhoAng(0)(iOm), gNorm, -1., -1.,
                                          gNorm, SxXCFunctional::alphaHybrid);
                  else if (xcFunctional == SxXC::HSE06)
                     xc.computeHSE (rhoAng(0)(iOm), gNorm, -1., -1., gNorm);
                  else { SX_EXIT; }

                  //SX_STOP_TIMER (Timer::XcFunc);
                  if (gNorm > 1e-16)  {
                     double vXcGradNorm = xc.eXc_gradRhoTl + xc.eXc_grad(0);
                     // vXc from radial gradient
                     vXcGradRad(iOm, 0) = vXcGradNorm
                                        * (aGrid.getXyz(iOm) ^ grad) / gNorm;
                     // vXc from angular gradient (without factor 4pi/r)
                     Coord g = (aGrid.weights(iOm) * vXcGradNorm / gNorm)
                             * grad;
                     int offset = nLm * iOm;
                     for (int ilm = 0; ilm < nLm; ++ilm)
                        vXcGradAng(ilm) += (g ^ Psilm3(ilm + offset));
                  } else {
                     vXcGradRad(iOm) = 0.;
                  }
               }

               SX_CHECK_NUM (xc.vXc(0));
               vXcAng(0)(iOm) =  xc.vXc(0) * aGrid.weights(iOm);
               eXc_ += aGrid.weights(iOm) * xc.eXc;
            }
            energy(ir) += FOUR_PI * eXc_;
         } else {
            SxVector3<Double> gradSpin[2];
            double gSpinNorm[2];
            double eXc_ = 0.;
            for (int iOm = 0; iOm < nOm; ++iOm)  {
               if (!gga)  {
                  if (!ldaPW)   {
                     // --- LSDA Perdew/Zunger '81
                     xc.computeLSDA (rhoAng(0)(iOm),rhoAng(1)(iOm));
                  } else   {
                     // --- LSDA Perdew/Wang '91
                     xc.computeLSDA_PW (rhoAng(0)(iOm),rhoAng(1)(iOm));
                  }
               } else {
                  // --- GGA

                  // --- get gradients
                  //SX_START_TIMER (Timer::GetGrad);
                  for (int iSpin = 0; iSpin < nSpin; ++iSpin)
                     for (int idir = 0; idir < 3; ++idir)
                        gradSpin[iSpin](idir) = gradAng(iOm, idir + 3 * iSpin);
                  SxVector3<Double> grad = gradSpin[0] + gradSpin[1];
                  double gNorm   = grad.norm ();
                  gSpinNorm[0] = gradSpin[0].norm ();
                  gSpinNorm[1] = gradSpin[1].norm ();
                  //SX_STOP_TIMER (Timer::GetGrad);


                  // compute functional
                  //SX_START_TIMER (Timer::XcFunc);
                  if (xcFunctional == SxXC::PBE)
                     xc.computePBE (rhoAng(0)(iOm), gSpinNorm[0],
                                    rhoAng(1)(iOm), gSpinNorm[1], gNorm);
                  else if (xcFunctional == SxXC::PBE0)
                     xc.computePBEHybrid (rhoAng(0)(iOm), gSpinNorm[0],
                                          rhoAng(1)(iOm), gSpinNorm[1], gNorm,
                                          SxXCFunctional::alphaHybrid);
                  else if (xcFunctional == SxXC::HSE06)
                     xc.computeHSE (rhoAng(0)(iOm), gSpinNorm[0],
                                    rhoAng(1)(iOm), gSpinNorm[1], gNorm);
                  else { SX_EXIT; }
                  //SX_STOP_TIMER (Timer::XcFunc);

                  //SX_START_TIMER (Timer::XcG1);
                  // ---  vXc from gradients
                  Coord gTot, xyz = aGrid.getXyz(iOm);
                  double vXcGradRadTot;
                  if (gNorm > 1e-16)  {
                     // ---  vXc from total density gradient
                     // vXc from radial gradient (total density)
                     vXcGradRadTot = xc.eXc_gradRhoTl * (xyz ^ grad) / gNorm;
                     // for vXc from angular gradient...
                     gTot = (aGrid.weights(iOm) * xc.eXc_gradRhoTl / gNorm)
                          * grad;
                  } else {
                     vXcGradRadTot = 0.;
                     gTot = 0.;
                  }

                  for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
                     if (gSpinNorm[iSpin] > 1e-16)  {
                        // ---  vXc from spin density gradient
                        // vXc from radial gradient
                        vXcGradRad(iOm,iSpin) = vXcGradRadTot
                                         + xc.eXc_grad(iSpin) / gSpinNorm[iSpin]
                                           * (xyz ^ gradSpin[iSpin]);
                        // vXc from angular gradient (without factor 4pi/r)
                        double facSpin = aGrid.weights(iOm)
                                       * xc.eXc_grad(iSpin) / gSpinNorm[iSpin];
                        Coord g = gTot + facSpin * gradSpin[iSpin];
                        int offset = nLm * iOm;
                        for (int ilm = 0; ilm < nLm; ++ilm)
                           vXcGradAng(ilm, iSpin) += (g ^ Psilm3(ilm + offset));
                     } else if (gNorm > 1e-16) {
                        vXcGradRad(iOm,iSpin) = vXcGradRadTot;
                        // vXc from angular gradient (without factor 4pi/r)
                        int offset = nLm * iOm;
                        for (int ilm = 0; ilm < nLm; ++ilm)
                           vXcGradAng(ilm,iSpin) += (gTot^Psilm3(ilm + offset));
                     }
                  }
                  //SX_STOP_TIMER (Timer::XcG1);
               }
               SX_CHECK_NUM (xc.vXc(0));
               SX_CHECK_NUM (xc.vXc(1));
               vXcAng(0)(iOm) = xc.vXc(0) * aGrid.weights(iOm);
               vXcAng(1)(iOm) = xc.vXc(1) * aGrid.weights(iOm);
               eXc_ += aGrid.weights(iOm) * xc.eXc;
            }
            energy(ir) += FOUR_PI * eXc_;
         }
         //SX_STOP_TIMER (Timer::Xc);
         SX_CHECK_NUM(energy(ir));

         for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
            VALIDATE_VECTOR (vXcAng(iSpin));
            // exchange-correlation potential
            vXcLm = vXcAng(iSpin) ^ aGrid.Ylm;

            // GGA: vXc contribution from angular gradient
            if (gga) vXcLm.plus_assign_ax(1./rad(ir),vXcGradAng.colRef(iSpin));
            for (int lm = 0; lm < nLm; ++lm)  {
               vXc(iSpin).meshData(ir,lm) = FOUR_PI * vXcLm(lm);
            }
            if (gga)  {
               // GGA: vXc contribution from radial gradient
               // get d Exc / d rho'(ilm, ir)  (integration factor 4pi is below)
               vXcLm = (vXcGradRad.colRef(iSpin) * aGrid.weights) ^ aGrid.Ylm;
               vXcGradRadAll.colRef(ir + nr * iSpin) <<= vXcLm;

            }
         }
      }
      // --- now apply contribution of the radial gradient to vXc,
      //     cannot be done within openmp-parallel ir loop because each
      //     ir contributes to its neighbours
      if (gga)  {
#       ifdef USE_OPENMP
#          pragma omp for nowait collapse(2) schedule(dynamic,1)
#        endif
         for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
            for (int lm = 0; lm < nLm; ++lm)  {
               for (int ir = 0; ir < nr; ++ir)  {
                  // now multiply with numerical gradient kernel
                  // d rho'(ir) / d rho(jr) * weight(ir)
                  // weight(ir) = rad(ir)^2 * Simpson-weight(ir,nr)
                  double gradKernel = vXcGradRadAll(lm, ir + nr * iSpin);
                  radialGrad.addVariation (sqr(rad(ir)) * weightSimpson (ir,nr)
                                           * gradKernel * FOUR_PI,
                                           vXcGradR(iSpin).meshData,
                                           ir, lm);
               }
               // add derivative wrt radial grad (without integration weights)
               // new loop, because previous addVariation is non-local in ir
               for (int ir = 0; ir < nr; ++ir)
                  vXc(iSpin).meshData(ir,lm) += vXcGradR(iSpin).meshData(ir,lm)
                     / (dex * rad(ir)*rad(ir)*rad(ir) * weightSimpson (ir,nr));
            }
         }
      }
   }


   energy *= rad.cub ();

   eXc += energy.integrate (radBasis.logDr(is))
              /* + integrate0 (energy) */;
   //cout << "eXcR0=" << energy(0) / 3. << endl;

   return vXc;
}

SxDiracVec<Double>
SxRadialAtom::computeXC (const SxDiracVec<Double> &rad,
                             double logDr,
                             const SxDiracVec<Double> &rho,
                             SxXC::XCFunctional xcFunctional,
                             double *eXcPtr,
                             SxXCFunctional::WhatToCompute xcMode,
                             bool simpsonWeights)
{
   int nr = (int)rad.getSize ();
   int nSpin = 1;
   SX_CHECK (rho.getSize () == nr, rho.getSize (), nr);

   SxDiracVec<Double> vXc(nr), vXcGradR;
   vXc.set (0.);
   if (eXcPtr) *eXcPtr = 0.;

   SxXCFunctional xc(nSpin, xcMode | SxXCFunctional::ComputeEandV);

   bool gga = (   xcFunctional == SxXC::PBE
               || xcFunctional == SxXC::PBE0
               || xcFunctional == SxXC::HSE06);
   bool ldaPW = xcFunctional == SxXC::LDA_PW;
   if (gga)  {
      vXcGradR.resize (nr);
      vXcGradR.set (0.);
   }
   int numGradOrder = 4;
   SxDifferential radialGrad (numGradOrder);

   double dRhoDr, r;

   for (int ir = 0; ir < nr; ++ir)  {
      r = rad(ir);
      if (gga)  {
         dRhoDr = radialGrad.apply (rho, ir) / (logDr * r);
      }

      // compute xc potential on radial grid
      if (!gga)  {
         // --- LDA
         if (!ldaPW)   {
            // --- LDA Perdew/Zunger '81
            xc.computeLDA (rho(ir));
         } else   {
            // --- LSDA Perdew/Wang '91
            double rhoHalf = 0.5 * rho(ir);
            xc.computeLSDA_PW(rhoHalf,rhoHalf);
         }
      } else {
         double gNorm = fabs(dRhoDr);
         if (xcFunctional == SxXC::PBE)
            xc.computePBE (rho(ir), gNorm, -1., -1., gNorm);
         else if (xcFunctional == SxXC::PBE0)
            xc.computePBEHybrid (rho(ir), gNorm, -1., -1., gNorm,
                                 SxXCFunctional::alphaHybrid);
         else if (xcFunctional == SxXC::HSE06)
            xc.computeHSE (rho(ir), gNorm, -1., -1., gNorm);
         else { SX_EXIT; }
         if (gNorm > 1e-16)  {
            // vXc from radial gradient
            double vXcGradRad = xc.eXc_gradRhoTl + xc.eXc_grad(0);
            if (dRhoDr < 0.) vXcGradRad = -vXcGradRad;
            SX_CHECK_NUM (vXcGradRad);
            double weight = r * r;
            if (simpsonWeights) weight *= weightSimpson (ir,nr);
            radialGrad.addVariation (weight * vXcGradRad, vXcGradR, ir);
         }
      }
      vXc(ir) = xc.vXc(0);
      if (eXcPtr)
         *eXcPtr += xc.eXc * rad(ir) * rad(ir) * rad(ir) * weightSimpson(ir,nr);
   }
   if (eXcPtr) *eXcPtr *= FOUR_PI * logDr;

   if (gga)  {
      // add derivative wrt radial grad (without integration weights)
      for (int ir = 0; ir < nr; ++ir)  {
         double weight = logDr * rad(ir)*rad(ir)*rad(ir);
         if (simpsonWeights) weight *= weightSimpson (ir,nr);
         vXc(ir) += vXcGradR(ir) / weight;
      }
   }

   return vXc;
}

SxArray<SxRadialMesh>
SxRadialAtom::computeXC2 (const SxArray<SxRadialMesh> &rho,
                              double *eXcPtr)
{
   SX_CLOCK (Timer::RadXC);
   SX_CHECK (eXcPtr);
   SX_CHECK (rho.getSize () == 1, rho.getSize ());
   const SxRadBasis &radBasis = rho(0).getBasis ();
   int is = rho(0).getIs ();
   const SxDiracVec<Double> &rad = radBasis.radFunc(is);
   int nSpin = int(rho.getSize ());
   int nr = rho(0).getNRad ();
   int nLm = (int)rho(0).meshData.nCols ();

   SxArray<SxRadialMesh> vXc(nSpin);
   for (int iSpin = 0; iSpin < nSpin; ++iSpin) {
      vXc(iSpin).resize (rho(0));
      vXc(iSpin).set (0.);
   }

   SxXCFunctional xc(nSpin,  SxXCFunctional::ComputeExcVxc 
                           | SxXCFunctional::ComputeKernel);
   SxDiracVec<Double> energy(radBasis(is));
   for (int ir = 0; ir < nr; ++ir)  {
      xc.computeLDA (rho(0)(ir,0,0) * SQRT_1_4PI); // * Y_00

      // exchange-correlation energy, l=m=0
      energy(ir) = FOUR_PI * xc.eXc; // 4pi = angular factor
      for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
         // exchange-correlation potential, l=m=0
         vXc(iSpin)(ir,0,0) = xc.vXc(iSpin) * sqrt(FOUR_PI);

         for (int lm = 1; lm < nLm; ++lm)  {
            for (int jSpin = 0; jSpin < nSpin; ++jSpin)
               // --- exchange-correlation potential, l!=0
               /// Ref. 1, Eq. (29), derivative wrt rho_lm
               vXc(iSpin).meshData(ir,lm) += xc.fXc(iSpin)(jSpin) 
                                             * rho(jSpin).meshData(ir,lm);
            {
               SX_EXIT;
               // missing contribution to vXc(l=0)
               // from (d fXc(iSpin)(jSpin) / d rho) rho(iSpin)rho(jSpin)
            }

            // --- exchange-correlation energy, l!=0
            /// Ref. 1, Eq. (29)
            energy(ir) += 0.5 * rho(iSpin).meshData(ir,lm) 
                              * vXc(iSpin).meshData(ir,lm);
         }
      }
   }

   energy *= rad.cub ();
   *eXcPtr += energy.integrate (radBasis.logDr(is))
              /* + integrate0 (energy) */;
   
   return vXc;
}

double SxRadialAtom::integrate0(const SxDiracVec<Double> &x)
{
   const SxRadBasis &rad = x.getBasis<SxRadBasis> ();
   const SxDiracVec<Double> &r = rad.radFunc(x.handle->auxData.is);
   double r0 = r(0);
   double m = (x(0) - x(1)) / (r0 - r(1));
   return r0*r0*r0/3. * (x(0) + 0.25 * m * r0);
}

SxDiracVec<Double> 
SxRadialAtom::laplace (const SxDiracVec<Double> &f,
                       const SxDiracVec<Double> &r,
                       double logDr, int l)
{
   SX_CHECK (l >= 0, l);
   int order = 4;
   SxDifferential grad (order);
   SxDiracVec<Double> res;
   res = grad.apply (r * grad.apply (f)) / r;
   res /= sqr(logDr);
   res.plus_assign_ax (-l *(l+1), f);
   res /= r.sqr ();
   return res;
}

