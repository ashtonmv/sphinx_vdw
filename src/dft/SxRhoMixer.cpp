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

//Ref1: PRB (54), 11169-11186 (1996)
//Ref2: Comp. Mat. Sci. (6), 15-50 (1996)
//Ref3: M. Palummo, R. Del Sole, L. Reining, F. Bechstedt, G. Cappellini,
//      Solid State Communications 95, 393-398 (1995).
//Ref4: G. Cappellini, R. Del Sole, L. Reining, F. Bechstedt,
//      Phys. Rev. B 47, 9892 (1993)

#include <SxRhoMixer.h>
#include <SxSymMatrix.h>
#include <SxRho.h>
#include <SxPAWRho.h>

SxRhoMixer::SxRhoMixer (MixerType typeIn, double rhoMixingIn, int maxStepsIn)
   : linearVacMixing(false),
     adaptiveScaling(false),
     residueProfile(false),
     maxSteps (maxStepsIn),
     type (typeIn),
     renormModus (SxRhoMixer::RenormOn),
     spinMixing (-1.),
     normR (-1.),
     normRPol (-1.),
     UP (SxXC::UP), DN (SxXC::DOWN)
{
   rhoMixing = rhoMixingIn;
   maxVacDensity   = 3e-5; // not well tested
   vacSmoothFactor = 1.5;  // not well tested
}

SxRhoMixer::SxRhoMixer (const SxSymbolTable *table, bool spin)
   : linearVacMixing(false),
     adaptiveScaling(false),
     residueProfile(false),
     renormModus (SxRhoMixer::RenormOn),
     spinMixing (-1.),
     normR (-1.),
     normRPol (-1.),
     UP (SxXC::UP), DN (SxXC::DOWN)
{
   maxVacDensity   = 3e-5; // not well tested
   vacSmoothFactor = 1.5;  // not well tested
   if (table)  {
      readTable (table);
   } else {
      type = Pulay;
      preconditioner.useKerker (1., 1.);
      rhoMixing = 1.;
   }
   
   if (!spin) spinMixing = -1.;
}

void SxRhoMixer::readTable (const SxSymbolTable *cmd)
{
   SX_CHECK(cmd);
   try  {
      int mixingMethod  = cmd->contains("mixingMethod")
                        ? cmd->get("mixingMethod")->toInt()
                        : 2;
      switch (mixingMethod)  {
         case 0 : type = Linear; break;
         case 1 : type = Linear;
                  preconditioner.readTable (cmd);
                  cout << "Warning: mixingMethod PRECOND_LINEAR is deprecate. "
                          "Use mixingMethod=LINEAR and preconditioner group "
                          "instead." << endl;
                  break;
         case 2 : type = Pulay;
                  maxSteps  = cmd->contains ("nPulaySteps")
                            ? cmd->get("nPulaySteps")->toInt()
                            : 7;
                  break;
         case 3 : type = Pulay;
                  maxSteps  = cmd->contains ("nPulaySteps")
                            ? cmd->get("nPulaySteps")->toInt()
                            : 7;
                  preconditioner.readTable (cmd);
                  cout << "Warning: mixingMethod PRECOND_PULAY is deprecate. "
                          "Use mixingMethod=PULAY  and preconditioner group "
                          "instead." << endl;
                  break;
         default: sxprintf ("Unknown mixing method: %d\n", mixingMethod);
                  SX_QUIT;
      }
      if (cmd->containsGroup ("preconditioner"))
         preconditioner.readTable (cmd->getGroup("preconditioner"));
      else
         preconditioner.useKerker (1., 1.);
      
      // --- set rhoMixing
      if (cmd->contains("rhoMixing"))
         rhoMixing = cmd->get("rhoMixing")->toReal ();
      else if (preconditioner.getType () != SxPreconditioner::Identity)
         rhoMixing = 1.;
      else if (preconditioner.scaling < 1. - 1e-3)
         rhoMixing = 1.;
      else
         rhoMixing = 0.8;
      
      spinMixing = (cmd->contains("spinMixing"))
                 ?  cmd->get("spinMixing")->toReal()
                 :  -1.;

      if (cmd->contains ("linearVacMixing"))
         linearVacMixing = cmd->get("linearVacMixing")->toAttribute ();

      if (cmd->contains("residueProfile"))
         residueProfile = cmd->get("residueProfile")->toAttribute ();
      if (cmd->contains("adaptiveScaling"))
         adaptiveScaling = cmd->get("adaptiveScaling")->toAttribute ();

      if (cmd->containsGroup("filter"))
         filterPtr = SxPtr<SxRhoFilter>::create (cmd->getGroup ("filter"));
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}


SxRhoMixer::~SxRhoMixer ()
{
   // empty
}


void SxRhoMixer::setType (MixerType t)
{
   type = t;
}

void SxRhoMixer::setNormModus (NormHandling modus)
{
   renormModus = modus;
}


void SxRhoMixer::setMaxSteps (int m)
{
   maxSteps = m;
}

int SxRhoMixer::getMaxSteps ()
{
   return maxSteps;
}


void SxRhoMixer::print () const
{
   cout << "|      mixing scheme:          ";
   if (preconditioner.getType () != SxPreconditioner::Identity) 
      cout << "precond. ";
   switch (type)  {
      case Pulay  : cout << "Pulay"; break;
      case Linear : cout << "linear"; break;
      default : cout << "unknown";
   }
   cout << " mixer" << endl;
   
   if (type == Pulay)  {
      cout << "|      mixing steps:           " << maxSteps << endl;
   }
   if (linearVacMixing)  {
      cout << "|      linear vacuum mixing:   on" << endl;
      cout << "|         max. vacuum density: " << maxVacDensity << endl;
      cout << "|         smoothening factor:  "
           << pow(vacSmoothFactor,-log(1e-2)) << endl;
   }
   cout << "|      rho Mixing:             " << rhoMixing << endl;
   if (spinMixing >= 0.)  {
      cout << "|      spin Mixing:            " << spinMixing << endl;
   }
   preconditioner.print ();
}

void SxRhoMixer::addRhoIn (const SxVector<Double> &in)
{
   SxRho rho;
   rho.rhoR.resize (1);
   rho.rhoR(0).resize (in.getSize());
   for (int i=0; i < in.getSize(); ++i)  rho.rhoR(0)(i) = in(i);
   addRhoIn (rho);
}

void SxRhoMixer::addRhoIn (const SxDensity &in)
{
   rhoIn.append (in.getCopy ());
}



void SxRhoMixer::addRhoOut (const SxVector<Double> &out)
{
   SxRho rho;
   rho.rhoR.resize (1);
   rho.rhoR(0).resize (out.getSize());
   for (int i=0; i < out.getSize(); ++i)  rho.rhoR(0)(i) = out(i);
   addRhoOut (rho);
}


void SxRhoMixer::addRhoOut (const SxDensity &out)
{
   SX_CLOCK(Timer::RhoMixing);
   rhoOut << out.getCopy ();
   if (renormModus == RenormOn)
      rhoOut.last ().renormalize ();

   SX_CHECK (rhoIn.getSize() == rhoOut.getSize(),
             rhoIn.getSize(),   rhoOut.getSize());
}


void SxRhoMixer::addResidue (const SxDensity &resIn)
{
   SX_CHECK (rhoIn.getSize() > 0, rhoIn.getSize());
   
   // --- update of rhoOut
   SxDensity out = rhoIn.last ().getCopy ();
   out += resIn;
   rhoOut << out;

   SX_CHECK (rhoIn.getSize() == rhoOut.getSize(),
             rhoIn.getSize(), rhoOut.getSize());
}


SxDensity SxRhoMixer::getRhoIn (int i) const
{
   SX_CHECK (i >= 0 && i < rhoIn.getSize(), i, rhoIn.getSize());

   return rhoIn(i);
}

SxDensity SxRhoMixer::getMixedRho ()
{
   SX_CLOCK(Timer::RhoMixing);
   SX_CHECK  (rhoIn.getSize()  >= 1, rhoIn.getSize());
   SX_CHECK  (rhoOut.getSize() >= 1, rhoOut.getSize());
   SX_CHECK (rhoIn.getSize() == rhoOut.getSize(),
             rhoIn.getSize(), rhoOut.getSize());

   int m = (int)rhoIn.getSize() - 1;

   // update preconditioner (some need the density)
   preconditioner.setRho (rhoOut(m));
 
   // --- compute new optimal input density
   SxDensity rhoOpt;
   if (type == Linear)  rhoOpt = getMixedRhoLinear ();
   else                 rhoOpt = getMixedRhoPulay ();

   {
      SxDensity R = rhoOut(m) - rhoIn(m);
      normR = sqrt(R.normSqr ());
      if (R.hasSpin ())
         normRPol = sqrt(R.spin ().normSqr ());
   }

   // --- SxRho stuff
   if (rhoOpt.checkType<SxRho> ())  {
      // --- compute norm of residue
      RhoR &rhoRIn  = rhoIn(m).getRef<SxRho> ().rhoR;
      RhoR &rhoROut = rhoOut(m).getRef<SxRho> ().rhoR;
      int nSpin = (int)rhoRIn.getSize ();

      // --- profile residue for debugging
      if (residueProfile) SX_MPI_MASTER_ONLY
      {
         SxMeshR R = rhoROut(0) - rhoRIn(0), RPol;
         bool spin = nSpin > 1;
         if (spin)  {
            RPol.copy (R);
            R    += rhoROut(1) - rhoRIn(1);
            RPol -= rhoROut(1) - rhoRIn(1);
         }
         const SxRBasis *r 
            = dynamic_cast<const SxRBasis *>(R.getBasisPtr ());
         SX_CHECK(r);
         const SxGBasis &g = r->getGBasis ();
         PsiG RinG = R.to (g), RPolInG;
         if (spin) RPolInG = RPol.to (g);
         FILE *fp = fopen("resprofile.dat","a");
         double norm2RinG = 0.L, norm2R = R.normSqr () * r->dOmega, c2;
         if (fp) {
            PsiG::Iterator itR = RinG.begin ();
            for (int ig = 0; ig < g.ng; ++ig, ++itR)  {
               norm2RinG += ( c2 = (*itR).absSqr () );
               if (spin)  {
                  fprintf(fp,"%f\t%.6e\t%.6e\n", sqrt(g.g2(ig)),
                             c2, RPolInG(ig).absSqr ());
               } else {
                  fprintf(fp,"%f\t%.6e\n", sqrt(g.g2(ig)), c2);
               }
            }
            fprintf(fp, "&\n");
            fclose (fp);
         } else {
            norm2RinG = RinG.normSqr ();
         }
         
         cout << "low  noise: " << (RinG(0).absSqr () / norm2R * 100.) << "%\n";
         cout << "high noise: " << ((1. - norm2RinG / norm2R) * 100.) << "%\n"; 
      }
      
      // --- linear mixing in low-density regions
      if (linearVacMixing && m > 0 && renormModus == RenormOn)  {
         RhoR &rhoROpt = rhoOpt.getRef<SxRho> ().rhoR;
         int ir, nr = (int)rhoROpt(0).getSize ();
         double avgRes = normR ;
         if (avgRes > maxVacDensity) avgRes = maxVacDensity;
         
         double oldMix = 1. - rhoMixing;
         double rho;
         SxMeshR::Iterator rhoIt, mixedIt, oldIt;
         SxMeshR rhoMix, rhoSpin;
         SxMeshR rhoInSpinAvg, rhoOutSpinAvg;

         if (nSpin == 1)  {
            rhoIt = rhoROut(0).begin ();
            oldIt = rhoRIn(0).begin ();
            mixedIt = rhoROpt(0).begin ();
         } else {
            rhoOutSpinAvg = .5 * (rhoROut(UP) + rhoROut(DN));
            rhoInSpinAvg  = .5 * (rhoRIn(UP)  + rhoRIn(DN));
            rhoMix  = .5 * (rhoROpt(UP) + rhoROpt(DN));
            rhoSpin = .5 * (rhoROpt(UP) - rhoROpt(DN));
            rhoIt = rhoOutSpinAvg.begin ();
            oldIt = rhoInSpinAvg.begin ();
            mixedIt = rhoMix.begin ();
         }
         
         const double smoothExp = 1. / log (vacSmoothFactor);
         /*  The following function is a Fermi-like function for the
             ln(minRho). ln(avgRes) plays the role of the Fermi energy, 
             while ln(smoothFactor) plays the role of ekt. The values above
             are not well tested.
         */
         double scaling;
         
         int count = 0;
         for (ir = 0; ir < nr; ++ir, ++rhoIt, ++mixedIt, ++oldIt)  {
            scaling = 1. / (1. + pow(fabs(rho = *rhoIt)/avgRes, smoothExp));
            if (scaling > 0.01)  {
               if (scaling > 0.99) scaling = 1.;
               *mixedIt = (1. - scaling) * *mixedIt
                        + scaling        * (oldMix * *oldIt + rhoMixing * rho);
               ++count;
            }
         }
         if (nSpin > 1)  {
            rhoROpt(UP) = rhoMix + rhoSpin;
            rhoROpt(DN) = rhoMix - rhoSpin;
         }
         cout << "Vacuum: partial linear mixing around " << avgRes << " -> "
              << count << " times (" << (double(count) / double(nr) * 100.) 
              << "%) of all points. " << endl;
      }
   }
   // --- profile residue for debugging
   if (residueProfile && rhoOpt.checkType<SxPAWRho> ()) SX_MPI_MASTER_ONLY
   {
      RhoR &rhoRIn  = rhoIn(m).getRef<SxPAWRho> ().pwRho.rhoR;
      RhoR &rhoROut = rhoOut(m).getRef<SxPAWRho> ().pwRho.rhoR;
      SxMeshR R = rhoROut(0) - rhoRIn(0), RPol;
      int nSpin = (int)rhoRIn.getSize ();
      bool spin = nSpin > 1;
      if (spin)  {
         RPol.copy (R);
         R    += rhoROut(1) - rhoRIn(1);
         RPol -= rhoROut(1) - rhoRIn(1);
      }
      const SxRBasis *r 
         = dynamic_cast<const SxRBasis *>(R.getBasisPtr ());
      SX_CHECK(r);
      const SxGBasis &g = r->getGBasis ();
      PsiG RinG = R.to (g), RPolInG;
      if (spin) RPolInG = RPol.to (g);
      FILE *fp = fopen("resprofile.dat","a");
      double norm2RinG = 0.L, norm2R = R.normSqr () * r->dOmega, c2;
      if (fp) {
         PsiG::Iterator itR = RinG.begin ();
         for (int ig = 0; ig < g.ng; ++ig, ++itR)  {
            norm2RinG += ( c2 = (*itR).absSqr () );
            if (spin)  {
               fprintf(fp,"%f\t%.6e\t%.6e\n", sqrt(g.g2(ig)),
                          c2, RPolInG(ig).absSqr ());
            } else {
               fprintf(fp,"%f\t%.6e\n", sqrt(g.g2(ig)), c2);
            }
         }
         fprintf(fp, "&\n");
         fclose (fp);
      } else {
         norm2RinG = RinG.normSqr ();
      }
      
      cout << "low  noise: " << (RinG(0).absSqr () / norm2R * 100.) << "%\n";
      cout << "high noise: " << ((1. - norm2RinG / norm2R) * 100.) << "%\n"; 
   }
   
   // --- linear mixing in low-density regions
   if (linearVacMixing && m > 0 && rhoOpt.checkType<SxPAWRho> ())  {
      RhoR &rhoRIn  = rhoIn(m).getRef<SxPAWRho> ().pwRho.rhoR;
      RhoR &rhoROut = rhoOut(m).getRef<SxPAWRho> ().pwRho.rhoR;
      RhoR &rhoROpt = rhoOpt.getRef<SxPAWRho> ().pwRho.rhoR;
      int ir, nr = (int)rhoROpt(0).getSize (), nSpin = (int)rhoROpt.getSize ();
      double avgRes = normR;
      if (avgRes > maxVacDensity) avgRes = maxVacDensity;
      
      double oldMix = 1. - rhoMixing;
      double rho;
      SxMeshR::Iterator rhoIt, mixedIt, oldIt;
      SxMeshR rhoMix, rhoSpin;
      SxMeshR rhoInSpinAvg, rhoOutSpinAvg;

      if (nSpin == 1)  {
         rhoIt = rhoROut(0).begin ();
         oldIt = rhoRIn(0).begin ();
         mixedIt = rhoROpt(0).begin ();
      } else {
         rhoOutSpinAvg = .5 * (rhoROut(UP) + rhoROut(DN));
         rhoInSpinAvg  = .5 * (rhoRIn(UP)  + rhoRIn(DN));
         rhoMix  = .5 * (rhoROpt(UP) + rhoROpt(DN));
         rhoSpin = .5 * (rhoROpt(UP) - rhoROpt(DN));
         rhoIt = rhoOutSpinAvg.begin ();
         oldIt = rhoInSpinAvg.begin ();
         mixedIt = rhoMix.begin ();
      }
      
      const double smoothExp = 1. / log (vacSmoothFactor);
      /*  The following function is a Fermi-like function for the
          ln(minRho). ln(avgRes) plays the role of the Fermi energy, 
          while ln(smoothFactor) plays the role of ekt. The values above
          are not well tested.
      */
      double scaling;
      
      int count = 0;
      double nElecOpt = 0., nElecNew = 0.;
      for (ir = 0; ir < nr; ++ir, ++rhoIt, ++mixedIt, ++oldIt)  {
         scaling = 1. / (1. + pow(fabs(rho = *rhoIt)/avgRes, smoothExp));
         nElecOpt += *mixedIt;
         if (scaling > 0.01)  {
            if (scaling > 0.99) scaling = 1.;
            *mixedIt = (1. - scaling) * *mixedIt
                     + scaling        * (oldMix * *oldIt + rhoMixing * rho);
            ++count;
         }
         nElecNew += *mixedIt;
      }
      if (nSpin > 1)  {
         rhoMix *= nElecOpt / nElecNew;
         rhoROpt(UP) = rhoMix + rhoSpin;
         rhoROpt(DN) = rhoMix - rhoSpin;
      } else {
         rhoROpt(0) *= nElecOpt/nElecNew;
      }
      cout << "Vacuum: partial linear mixing around " << avgRes << " -> "
           << count << " times (" << (double(count) / double(nr) * 100.) 
           << "%) of all points. " << endl;
      cout << "rescale: " << (nElecOpt/nElecNew - 1.) << endl;
   }


   // --- renormalize rhoOpt
   if (renormModus == SxRhoMixer::RenormOn)
      rhoOpt.renormalize ();

   // --- clean up
   if ((type == Pulay && m == maxSteps))  {
      if (res2(0) < res2.last ()) {
         cout << "Stagnating convergence detected. Reinitializing Pulay mixer."
              << endl;
         reset ();
      }
   }
   removeFirst ();

#ifdef USE_LOOPMPI
   rhoOpt.syncMPI ();
#endif

   return rhoOpt;
}


SxDensity SxRhoMixer::getMixedRhoLinear ()
{
   int m = (int)rhoIn.getSize()-1;
   SxDensity rhoOpt = rhoIn(m).getCopy (),
             R      = rhoOut(m) - rhoIn(m);

   if (filterPtr)  {
      SxDensity FR = (*filterPtr) | R;
      R += preconditioner * FR - FR;
      rhoOpt.plus_assign_ax (rhoMixing, R);
   } else {
      rhoOpt.plus_assign_ax (rhoMixing, preconditioner * R);
   }
   return rhoOpt;
}

SxDensity SxRhoMixer::getMixedRhoPulay ()
{
   int i, j, m = (int)rhoIn.getSize()-1;
   SxDensity rhoOpt;

   SX_CHECK (rhoIn.getSize() == rhoOut.getSize(),
             rhoIn.getSize(), rhoOut.getSize());

   int bestIdx = m;
   double bestR2 = 0.;
   if (!filterPtr) filterPtr = SxPtr<SxRhoFilter>::create ();
   // --- 1st step: linear mixing
   if (m == 0)  {
      rhoOpt = getMixedRhoLinear ();
   }  else  {
      // --- compute R and dR
      SxArray<SxDensity> R(m+1), dR(m);
      for (i = 0; i < m+1; ++i)  {
         const SxDensity &res = rhoOut(i) - rhoIn(i);
         R(i) = (*filterPtr) | res;
         double R2 = (i < res2.getSize ()) ? res2(i) : res.normSqr ();
         if (i >= res2.getSize ()) res2 << R2;
         if (i < bestIdx || R2 < bestR2) {
            bestIdx = i;
            bestR2 = R2;
         }
         if (i > 0)
            dR(i-1) = R(i) - R(i-1); // ref2 (88d)
      }

      // --- compute Pulay matrix
      SxSymMatrix<TPrecRhoR> A(m);
      SxVector<TPrecRhoR>    B(m), alpha;
      A.set (0.);
      B.set (0.);
      for (i=0; i < m; i++)  {
         for (j=i; j < m; j++)  {
            A(i,j) = (dR(i) | dR(j));   // ref2, (91)
         }
         B(i) = (dR(i) | R(bestIdx)); // ref2, (90) <dR|R_m>
      }

      // TODO: expand is to be removed
      alpha = -(A.inverse().expand() ^ B);     // ref2, (90)
      
      cout << "alpha = " << alpha << endl;
      // cout << "A = " << A.expand() << endl;
      // cout << "B = " << B << endl;
      double r2Now = R(bestIdx) | R(bestIdx),
             r2New = dot(rhoMixing * alpha,
                         (A.expand () ^ (rhoMixing * alpha)) + 2. * B) + r2Now;
      cout << "Pulay R now      : " << sqrt(r2Now) << endl;
      cout << "Pulay R predicted: " << sqrt(r2New) << endl;

      SxDensity mixedR;
      rhoOpt = rhoOut(bestIdx).getCopy ();
      SxDensity dRhoMixed = rhoIn(0) - rhoIn(0); // = zero
      mixedR = R(bestIdx).getCopy ();
      for (i=0; i < m; i++)  {
         // ref2 (92) and (88b)
         dRhoMixed.plus_assign_ax (alpha(i), rhoIn(i+1) - rhoIn(i));
         mixedR.plus_assign_ax (alpha(i), dR(i));
      }
      mixedR = preconditioner * mixedR - R(bestIdx);
      mixedR += (*filterPtr) | dRhoMixed;
      rhoOpt += mixedR;
      /*
      {
         SxDensity change     = rhoOpt - rhoIn(m);
         SxDensity lastChange = rhoIn(m) - rhoIn(m-1);
         double chNorm = sqrt(change.normSqr ());
         double lchNorm = sqrt(lastChange.normSqr ());
         double dotp = (change | lastChange);
         cout << "colinearity = " << fabs(dotp/chNorm / lchNorm) << endl;
         cout << "projection  = " << (dotp/lchNorm / lchNorm) << endl;
      }
      */
   }
   if (fabs(rhoMixing - 1.) > 1e-6 || spinMixing >= 0.)  {
      if (bestIdx != m)  {
         cout << "Mixing from best step (" << (m-bestIdx)
              << " before)" << endl;
      }
      rhoOpt.plus_assign_ax (1. - rhoMixing, rhoIn(bestIdx) - rhoOpt);
      if (spinMixing >= 0. && rhoOpt.hasSpin ())  {
         rhoOpt.plus_assign_aspin (spinMixing - rhoMixing, rhoOpt.spin () - rhoIn(bestIdx).spin ());
      }
   }
   cout << "final Pulay step: " << ((rhoOpt - rhoIn(bestIdx)).normSqr ()) << endl;

   return rhoOpt;
}


void SxRhoMixer::removeFirst ()
{
   SX_CHECK(rhoIn.getSize() == rhoOut.getSize());
   if (rhoIn.getSize() <= maxSteps) return;
   rhoIn.removeFirst ();
   rhoOut.removeFirst ();
   res2.removeFirst ();
}

void SxRhoMixer::reset ()
{
   // --- kick out the oldest density to garantuee some change
   rhoIn.removeFirst ();
   rhoOut.removeFirst ();
   res2.removeFirst ();
   // --- kick out all but the best densities
   while (res2.getSize () > 1)  {
      if (res2(0) < res2.last ())  {
         rhoIn.removeLast ();
         rhoOut.removeLast ();
         res2.removeLast ();
      } else {
         rhoIn.removeFirst ();
         rhoOut.removeFirst ();
         res2.removeFirst ();
      }
   }
   if (res2.getSize ())
      cout << "Keeping density with R=" << sqrt(res2(0)) << endl;;
   /*
   rhoIn.resize  (0);
   rhoOut.resize (0);
   res2.resize   (0);
   */

}

