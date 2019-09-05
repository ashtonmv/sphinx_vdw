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

#include <SxPreconditioner.h>
#include <SxGBasis.h>
#include <SxRho.h>
#include <SxPAWRho.h>
#include <SxRhoFilter.h>
#include <SxProjector.h>
#include <SxHamiltonian.h>

SxPreconditioner::SxPreconditioner ()
: scaling(1.),
  spinScaling (-1.),
  type(SxPreconditioner::Identity),
  kerkerDamping(-1.),
  dielecConstant(-1.)
{
   /* empty */
}

void SxPreconditioner::readTable(const SxSymbolTable *table)
{
   SX_CHECK(table);
   try  {
      if (table->getName() != "preconditioner")  {
         if (table->containsGroup("preconditioner"))  {
            table = table->getGroup("preconditioner");
         } else if (table->contains("kerkerDamping"))  {
            // old style Kerker
            cout << "Warning: deprecate input style for preconditioner." 
                 << endl;
            type = Kerker;
            kerkerDamping = table->get("kerkerDamping")->toReal ();
            scaling = table->get("kerkerScaling")->toReal ();
            return;
         }
      }
      int typeIn = table->get("type")->toInt ();
      if (typeIn < 0 || typeIn >= nTypes)  {
         cout << "Illegal preconditioner type " << typeIn << endl;
         SX_QUIT;
      }
      type = Type(typeIn);

      if (table->contains("scaling"))
         scaling = table->get("scaling")->toReal ();
      else
         scaling = (type == Kerker) ? 0.8 : 1.0;

      if (table->contains("spinScaling"))
         spinScaling = table->get("spinScaling")->toReal ();
      else
         spinScaling = scaling;

      if (type == Kerker)  {
         kerkerDamping = table->contains("kerkerDamping") 
                       ? table->get("kerkerDamping")->toReal ()
                       : 0.6;
      }
      if (type == Freysoldt)  {
         kerkerDamping = table->contains("kerkerDamping") 
                       ? table->get("kerkerDamping")->toReal ()
                       : 0.6;
      }
      if (type == CSRB)  {
         dielecConstant = table->get("dielecConstant")->toReal ();
      }

      /*
      if (type == Elliptic)  {
         bool useDipoleCorrection = false;
         if (table->parent->contains ("dipoleCorrection"))  {
            useDipoleCorrection = table->parent->get ("dipoleCorrection")
                                  ->toAttribute ();
         } else {
            const SxSymbolTable *ham
               = SxHamiltonian::getHamiltonianGroup(table->topLevel ());
            if (ham->contains ("dipoleCorrection"))
               useDipoleCorrection = ham->get("dipoleCorrection")->toAttribute ();
         }
         if (useDipoleCorrection)  {
            dipoleCorr = SxPtr<SxDipoleCorrZ>::create ();
         }
      }
      */
         
   } catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }

}

void SxPreconditioner::useKerker (double scalingIn, double dampingIn)
{
   type          = Kerker;
   scaling       = scalingIn;
   kerkerDamping = dampingIn;
   if (spinScaling == -1.) spinScaling = scalingIn;
}

SxMeshR SxPreconditioner::operator* (const SxMeshR &meshR) const
{
   if (type == Identity) return scaling * meshR;

   if (type == Elliptic) return scaling * solveElliptic (meshR);

   SX_CHECK (dynamic_cast<const SxRBasis*>(meshR.getBasisPtr ()));
   const SxRBasis &R = *dynamic_cast<const SxRBasis*>(meshR.getBasisPtr ());
   const SxGBasis &G = R.getGBasis ();

   return scaling * (meshR + (K * meshR.to(G)).to(R) );
}

PsiG SxPreconditioner::applyElliptic (const PsiG &meshG) const
{
   const SxRBasis &R = ellipticA.getBasis<SxRBasis> ();
   const SxGBasis &G = R.getGBasis ();
   PsiG V(G);
   V(0) = 0.;
   for (int ig = 1; ig < G.ng; ++ig)
      V(ig) = FOUR_PI * meshG(ig) / G.g2(ig);

   SxMeshR VR = R | V;
   if (dipoleCorr)  {
      // inactive - proved not helpful
      SX_EXIT;
      SxMeshR meshR = R | meshG;
      dipoleCorr->update (ellipticB, meshR);
      dipoleCorr->correctPotential (&VR);
      V = G | VR; // backtrafo
   }
   // --- now metallic polarization
   SxDiracVec<Double> inducedRho = VR * ellipticB;
   // --- now dielectric screening
   /*
   NOT TESTED!!!
   SxDiracVec<Double> ellipticA1 = (ellipticA - 1.)/FOUR_PI;
   if (ellipticA1.normSqr () > 1e-12 * ellipticA.getSize ())  {
      PsiG P(G);
      P.set (0.);
      for (int idir = 0; idir < 3; idir++)  {
         const PsiG &gVec = G.gVec.colRef(idir);
         P += (G | ((R | gVec * V) * ellipticA1)) * gVec;
      }
      inducedRho += R | P;
   }
   */
   // metallic neutralization (as if by potential shift)
   inducedRho -= (inducedRho.sum ()/ellipticB.sum ()) * ellipticB;

   return meshG + (G | inducedRho);
}

SxMeshR SxPreconditioner::solveElliptic (const SxMeshR &meshR) const
{
   PsiG Xold;
   const SxRBasis &R = *dynamic_cast<const SxRBasis*>(meshR.getBasisPtr ());
   const SxGBasis &G = R.getGBasis ();
   PsiG meshG = G | (meshR - (meshR.sum ()/ellipticB.sum ()) * ellipticB);
   PsiG solutionG = meshG.getCopy ();

   double norm2 = meshG.normSqr ();
   PsiG epsResultG = applyElliptic(solutionG);
   double dRes2Old = -1.;

   int reInit = 10;
   for (int it = 0; it < 500; ++it)  {
      PsiG dRes = meshG - epsResultG;
      double dRes2 = dRes.normSqr ();
      if (dRes2 < 1e-10 * norm2) {
         cout << "Elliptic solver: converge in it=" << (it+1)
              << " to rel. norm=" << sqrt(dRes2/norm2) << endl;
         break;
      }
      PsiG X = 1. * dRes;
      double gamma = 0.;
      if (it > 0)  {
         gamma = dRes2 / dRes2Old;
         X += gamma * Xold;
      }
      Xold = X;
      dRes2Old = dRes2;
      // line search
      PsiG epsX = applyElliptic (X);
      double step = dRes.normSqr () / dot (epsX, X);
      solutionG.plus_assign_ax (step, X);
      epsResultG = applyElliptic(solutionG);
      // update
      // epsResultG.plus_assign_ax (step, epsX);
      if (fabs(step) < 1e-4 || gamma > 0.95) reInit--;
      if (reInit == 0)  {
         cout << "Elliptic solver CG reinit at it=" << (it+1)
              << " rel. norm=" << sqrt(dRes2/norm2)
              << endl;
         Xold.set (0.);
         reInit = 10;
      }
   }
   return meshR + (R | (solutionG - epsResultG));
   //return meshR + (R | (solutionG - meshG));
}

void SxPreconditioner::setRho (const SxDensity &rho)
{
   if (rho.checkType<SxRho> ())  {
      const SxRho &pwRho = rho.getRef<SxRho> ();
      int nSpin = pwRho.getNSpin ();
      SX_CHECK (nSpin > 0 && nSpin <= 2, nSpin);
      if (nSpin == 1)
         setRho (pwRho.rhoR(0));
      else
         setRho (pwRho.rhoR(0) + pwRho.rhoR(1));
   } else if (rho.checkType<SxPAWRho> ())  {
      // for PAW: setRho with plane-wave density
      setRho (rho.getRef<SxPAWRho> ().pwRho);
   } else {
      SX_EXIT;
   }
}

void SxPreconditioner::setRho (const SxMeshR &rho)
{
   // --- do we need to set up K?
   if (K.getSize () > 0 || type == Identity) return;

   // --- get G basis
   const SxRBasis *rBasisPtr 
      = dynamic_cast<const SxRBasis *>(rho.getBasisPtr ());
   SX_CHECK (rBasisPtr);

   if (type == Elliptic)  {
      if (ellipticA.getSize () == 0)  {
         const SxGBasis &G = rBasisPtr->getGBasis ();
         double beta = 1.;
         SxDiracVec<Double> gauss = exp((-0.5 * beta*beta) * G.g2);
         SxMeshR rho2 = (*rBasisPtr) | (gauss * (G | rho));
         ellipticA.resize (rho.getSize ());
         ellipticB.resize (rho.getSize ());
         ellipticA.setBasis (rBasisPtr);
         ellipticB.setBasis (rBasisPtr);
         SX_LOOP(ir)  {
            if (rho2(ir) > 1e-3)  {
               ellipticA(ir) = 1.;
               ellipticB(ir) = 1./PI * cbrt(3./PI * rho2(ir));
               //ellipticB(ir) = 1.;
            } else {
               ellipticA(ir) = 1.; // should be dielectric constant
               ellipticB(ir) = 0.;
            }
         }
         // --- now broaden
         ellipticA = *rBasisPtr | (gauss * ( G | ellipticA) );
         ellipticB = *rBasisPtr | (gauss * ( G | ellipticB) );
         //SxRho(ellipticA).writeRho ("ellipticA.sxb");
         //SxRho(ellipticB).writeRho ("ellipticB.sxb");
         //SxRho(rho).writeRho ("rhoE.sxb");
      }
      return;
   }

   // --- get TF screening length
   double q2tf;
   if (type == Kerker) {
      // from input
      q2tf = kerkerDamping * kerkerDamping;
   } else  {
      // from average density
      double rhoAvg = rho.sum () / double(rho.getSize ());
      q2tf = 4./PI * cbrt (3. * PI * PI * rhoAvg);
   }
   
   setupK (rBasisPtr->getGBasis (), q2tf);
}

void SxPreconditioner::setupK(const SxGBasis &G, double q2tf)
{
   if (type == Identity) return;  // K not needed

   // --- setup specific
   if (type == Kerker || type == Freysoldt)  {
      // Kerker (scaling is done elsewhere)
      K = G.g2 / (G.g2 + q2tf); // ref2, (82)
   } else if (type == Lindhard)  {
      // Lindhard
      K = SxDiracVec<TReal8>(G);
      for (int ig = 0; ig < G.ng; ++ig)
        K(ig) = invLindhard(G.g2(ig), q2tf);
   } else if (type == CSRB)  {
      // CSRB
      K = SxDiracVec<TReal8>(G);
      for (int ig = 0; ig < G.ng; ++ig)
        K(ig) = invCSRB(G.g2(ig), q2tf, dielecConstant);
   }
   
   K += -1.; // we apply the preconditioner as difference

   if (G.g2(0) < 1e-10) K(0) = 0.; // ensure that K*R[m] does not change norm
}


inline double 
SxPreconditioner::invLindhard(double g2, double q2tf) 
{
   double gByPi = sqrt(g2) / PI;
   double y = 0.5 *  q2tf / gByPi;
   double ym = y - 1.;
   
   if (fabs(ym) < 1e-6) return 0.;
   double yp = y + 1.;
   double q2tfX2f = q2tf + gByPi * ym * yp * log(fabs(yp/ym));
   return g2 / (g2 + 0.5 * q2tfX2f);
}

double 
SxPreconditioner::invCSRB(double g2, double q2tf, double eps0)
{
   // --- ref 1 model
   const double alpha = 1.563; // empirical fit parameter, see ref 4.
   /*
   const double kFermi = PI/4. * q2tf;
   const double plasFreq2 = kFermi * kFermi * q2tf / 3.;
   double invDEps = 1./(eps0 -1.)
               + alpha * g2/q2tf
               + g2*g2/(4. * plasFreq2);
   */
   double g2Byq2 = g2/q2tf;
   double invDEps = 1./(eps0 -1.) + g2Byq2*(alpha + 12. * g2Byq2/(PI_2 * q2tf));

   // return 1. / ( 1. + 1. / invDEps);
   return invDEps/(invDEps+1.);

}

void SxPreconditioner::print () const
{
   cout << "|      scaling factor:         " << scaling << endl;
   if (fabs (spinScaling - scaling) > 1e-10)
      cout << "|      spin scaling factor:    " << spinScaling << endl;
   if (type == Identity) return;
   
   cout << "|      preconditioner type:    ";
   if (type == Lindhard)
      cout << "Lindhard" << endl;
   else if (type == CSRB)  {
      cout << "Cappellini/Del Sole/Reining/Bechstedt" << endl;
      cout << "|      dielectric constant:    " << dielecConstant << endl;
   } else  if (type == Kerker) {
      cout << "Kerker (Thomas-Fermi)" << endl;
      cout << "|      Kerker damping:         " << kerkerDamping << endl;
   } else if (type == Elliptic) {
      cout << "Elliptic" << endl;
   } else {
      cout << int(type) << endl;
   }
}

SxDensity SxPreconditioner::operator* (const SxDensity &in) const
{
   if (in.checkType<SxRho> ())  {
      const SxRho &rho = in.getRef<SxRho> ();
      SxPtr<SxRho> res = SxPtr<SxRho>::create ();
      int nSpin = rho.getNSpin ();
      res->rhoR.resize (nSpin);
      res->rBasisPtr = rho.rBasisPtr;
      SxMeshR rhoTot = (nSpin == 1) ? rho.rhoR(0)
                                    : rho.rhoR(0) + rho.rhoR(1);
      rhoTot = (*this) * rhoTot;
      if (nSpin == 2)  {
         SX_CHECK(spinScaling >= 0., spinScaling);
         SxMeshR rhoSpin = spinScaling * (rho.rhoR(0) - rho.rhoR(1));
         res->rhoR(0) = 0.5 * (rhoTot + rhoSpin);
         res->rhoR(1) = 0.5 * (rhoTot - rhoSpin);
      } else {
         res->rhoR(0) = rhoTot;
      }
      return SxDensity(res);
   } else if (in.checkType<SxPAWRho> ())  {
      // --- startup
      const SxPAWRho &rho = in.getRef<SxPAWRho> ();
      SxPtr<SxPAWRho> res = SxPtr<SxPAWRho>::create ();
      int nSpin = rho.getNSpin ();
      res->pwRho.rhoR.resize (nSpin);
      SX_CHECK (rho.pwRho.rBasisPtr);
      res->pwRho.rBasisPtr = rho.pwRho.rBasisPtr;
      res->potPtr = rho.potPtr;

      // --- precondition plane-wave pseudo-density
      SxMeshR rhoTot = (nSpin == 1) ? rho.pwRho.rhoR(0).getCopy ()
                                    : rho.pwRho.rhoR(0) + rho.pwRho.rhoR(1);
      // TODO: add compensation charges (from Dij) to rhoTot
      const SxGBasis &G = rho.pwRho.rBasisPtr->getGBasis ();
      SxRhoFilter filter (G.g2(G.ng-1));
      SxRho rhoComp = (filter | in).getRef<SxRho> ();
      if (nSpin == 2) rhoComp(0) += rhoComp(1);

      //rhoTot = (*this) * rhoTot;
      // This is what we want to get
      // PRECOND does the Hartree preconditioning
      // scaling is the global scaling factor
      // res = scaling * [PRECOND(Filter(rhoPS + rhoHat))
      //                  + rhoPS - Filter(rhoPS + rhoHat)]
      // = scaling * (rhoPS - Filter(rhoPS + rhoHat))
      // + scaling * PRECOND(Filter(rhoPS + rhoHat))
      rhoTot -= rhoComp(0);
      rhoTot *= scaling;
      // now the scaling * PRECOND(...) part
      rhoTot += (*this) * rhoComp(0);
      if (nSpin == 2)  {
         SX_CHECK(spinScaling >= 0., spinScaling);
         SxMeshR rhoSpin = spinScaling * (rho.pwRho.rhoR(0) - rho.pwRho.rhoR(1));
         res->pwRho.rhoR(0) = 0.5 * (rhoTot + rhoSpin);
         res->pwRho.rhoR(1) = 0.5 * (rhoTot - rhoSpin);
      } else {
         res->pwRho.rhoR(0) = rhoTot;
      }

      bool noScaling = fabs(scaling     -1.) < 1e-10,
           noSpinMix = fabs(spinScaling -1.) < 1e-10 || nSpin == 1;
      if (noScaling && noSpinMix)  {
         // "copy" 1-center density matrix (reference counting!)
         res->Dij = rho.Dij;
         res->blockAO = rho.blockAO;
      } else {
         res->Dij.resize (rho.Dij.atomInfo, nSpin);
         SxMatrix<Double> DijAvg, DijSpin;
         SX_LOOP2(is,ia)  {
            if (nSpin == 1)  { 
               res->Dij(0,is,ia) = scaling * rho.Dij(0,is,ia);
            } else {
               DijAvg  = 0.5 * scaling     * (  rho.Dij(0,is,ia) 
                                              + rho.Dij(1,is,ia)),
               DijSpin = 0.5 * spinScaling * (  rho.Dij(0,is,ia) 
                                              - rho.Dij(1,is,ia));
               res->Dij(0,is,ia) = DijAvg + DijSpin;
               res->Dij(1,is,ia) = DijAvg - DijSpin;
            }
         }
         ssize_t nExtra = rho.blockAO.getSize ();
         res->blockAO.resize (nExtra);
         for (int iAO = 0; iAO < nExtra; ++iAO)  {
            res->blockAO(iAO) = SxPtr<SxBlockDensityMatrix>::create ();
            SxBlockDensityMatrix &resBlock = *res->blockAO(iAO),
                                 &rhoBlock = *rho.blockAO(iAO);
            resBlock.resize (nSpin, rhoBlock.getNSite ());
            SX_LOOP(iSite)  {
               if (nSpin  == 1)  {
                  resBlock(0, iSite) = scaling * rhoBlock(0, iSite);
               } else {
                  SxMatrix<Double> &resUp   = resBlock(0, iSite),
                                   &resDown = resBlock(1, iSite),
                                   &rhoUp   = rhoBlock(0, iSite),
                                   &rhoDown = rhoBlock(1, iSite);
                  resUp = 0.5 * (scaling + spinScaling) * rhoUp;
                  resUp.plus_assign_ax(0.5 * (scaling - spinScaling), rhoDown);
                  resDown = 0.5 * (scaling + spinScaling) * rhoDown;
                  resDown.plus_assign_ax(0.5 * (scaling - spinScaling), rhoUp);
               }
            }
         }
      }
      // TODO: subtract mixed compensation charges (from Dij) from rhoTot
      return SxDensity(res);
   } else {
      SX_EXIT;
   }
}
