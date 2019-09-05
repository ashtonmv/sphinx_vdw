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

#include <SxPAWPot.h>
#include <SxDirac.h>
#include <SxProjector.h>
#include <SxBinIO.h>
#include <SxYlm.h>
#include <SxNaturalCubicSpline.h>
#include <SxRadialAtom.h>
#include <SxFileParser.h>
#include <SxRegex.h>
#include <SxFileIO.h>


SxString sciNum = "(-?\\d+.*[e|E][-|+]\\d*)";
SxString intNum = "(-?\\d+)";
SxString realNum = "(-?\\d+\\.\\d*)";

/*
 1. P. E. Bloechl, Phys. Rev. B 50, 17953 (1994).
 2. C. Freysoldt, PAW implementation notes
 3. N. A. W. Holzwarth, "Notes for revised form of atompaw code."
    Wake Forest University, Winston-Salem, NC.  January 28, 2008
    from http://www.wfu.edu/~natalie/papers/pwpaw/notes/atompaw/atompawEqns.pdf
 4. C. Freysoldt, HSE implementation notes
*/


namespace Timer {
   enum PAWPotTimer { ComputeRhoRadial, RadVMatEl, CoreXHSE, CoreXHSEjsb };
}

SX_REGISTER_TIMERS(Timer::PAWPotTimer)
{
   regTimer (Timer::ComputeRhoRadial, "PAW: rad. rho");
   regTimer (Timer::RadVMatEl, "1-center <i|V|j>");
   regTimer (Timer::CoreXHSE, "coreX (HSE)");
   regTimer (Timer::CoreXHSEjsb, "coreX (HSE:jsb)");
}

double toRSpace (double r,  const SxNaturalCubicSpline &fFine, double N,
                 int ngf, double dg, int l = 0)
{
   double sum = 0.;
   for (int ig = 0; ig < ngf; ++ig)  {
      double g = ig * dg;
      sum += g*g * SxYlm::jsb (l, g * r) * fFine.getVal (ig / N)
           * ((ig & 1) ? 4. : 2.);
   }
   sum *= dg / 3.;
   sum /= 2. * sqr(PI); // prefactor
   return sum;
}

SxVector<Double> lowPass (const SxVector<Double> &f, double eps = 1e-4)
{
   if (eps <= 0.) return f;
   int ng = (int)f.getSize ();
   SX_CHECK (eps > 0. && eps < 1e-2, eps);
   double fMax = 0.;
   for (int ig = 0; ig < ng; ++ig)  {
      if (fabs(f(ig)) > fMax) fMax = fabs(f(ig));
   }

   // --- suppress r-space wiggles (low-pass filtering)
   fMax *= eps;
   int igMin = 0;
   for (int ig = 0; ig < ng; ++ig)  {
      if (fabs(f(ig)) > fMax) igMin = ig;  
   }
   if (igMin == ng - 1) igMin = 0;
   double gMin = igMin, gMax = ng;
   //double beta = sqrt(-log(1e-9)) / (gMax - gMin);
   SxVector<Double> fSmooth(ng);
   if (igMin > 0) fSmooth(SxIdx(0,igMin-1)) <<= f(SxIdx(0,igMin-1));
   for (int ig = igMin; ig < ng; ++ig)  {
      //fSmooth(ig) = f(ig) * exp(-sqr(beta * (ig - gMin)));
      double x = (ig - gMin) / (gMax - gMin);
      fSmooth(ig) = f(ig) * (1. - x * x * (3. - 2. * x));
   }
   return fSmooth;
}

SxVector<Double> toRSpace (const SxVector<Double> &f, double gMax,
                           double rMax, double dr, double lowEps = 1e-4,
                           int l = 0)
{
   int ng = (int)f.getSize ();
   double dg = gMax / ng;

   // interpolating function
   SxNaturalCubicSpline fFine(lowPass(f, lowEps), true);

   SxStack<double> res;
   //int N = 40;
   int N = 10;
   for (double r = 0.; r <= rMax; r += dr)  {
      res << toRSpace (r, fFine, N, (ng - 1) * N, dg / N, l);
   }
   return SxVector<Double> (res);
}

SxDiracVec<Double> toRSpace (const SxVector<Double> &f, double gMax,
                           const SxDiracVec<Double> &r, int l = 0)
{
   int ng = (int)f.getSize (), nr = (int)r.getSize ();
   double dg = gMax / ng;

   // interpolating function
   SxNaturalCubicSpline fFine(lowPass(f, 1e-4), true);

   SxDiracVec<Double> res(nr);
   //int N = 40;
   int N = 10;
   for (int ir = 0; ir < nr; ++ir)  {
      res(ir) = toRSpace (r(ir), fFine, N, (ng - 1) * N, dg / N, l);
   }
   return res;
}




SxPAWPot::SxPAWPot ()
   : verbose(false)
{
   // empty
}

SxPAWPot::SxPAWPot (const SxSymbolTable *table)
   //try // catch at end of function
   : SxSpeciesData (table->getGroup("pawPot")),
     verbose(false)
{
   int is=0, nSpecies;
   SxArray<SxString> coreWaveFile;
   SxXC::XCFunctional xc;
   try {
      xc = SxXC::getXCFunctional(table->getGroup("PAWHamiltonian"));
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   bool hybridFunctional = SxXC::isHybrid (xc);
      
   try  {
      SxString potFile, potType;
      SxSymbolTable *pawGroup = table->getGroup("pawPot");
      SxSymbolTable *species;
      nSpecies = pawGroup->getGroup("species")->getNItems("species");

      if (pawGroup->contains ("verbose"))
         verbose = pawGroup->get ("verbose")->toAttribute ();

      // radial mesh
      logDr.resize (nSpecies);
      rc.resize (nSpecies);
      rad.resize (nSpecies);
      rCore.resize(nSpecies);

      // angular mesh for xc
      aGridType.resize (nSpecies);

      // projector's l-values
      lPhi.resize (nSpecies);
      lCore.resize (nSpecies);
      lMax.resize (nSpecies);
      lMaxRho.resize (nSpecies);

      // initialization densities
      rhoInit.resize (nSpecies);
      rhoInitBasis.resize (nSpecies);
      rhoInitQ0.resize (nSpecies);

      // core densities
      rhoCorePS.resize (nSpecies);
      rhoCoreAE.resize (nSpecies);
      // smoothening potential
      vBar.resize (nSpecies);

      // partial waves & projectors
      phiPS.resize (nSpecies);
      phiAE.resize (nSpecies);
      pPS.resize (nSpecies);
      psiCoreAE.resize (nSpecies);
      projBasis.resize (nSpecies);

      // kinetic energy and overlap matrices
      deltaKin.resize (nSpecies);
      deltaS.resize (nSpecies);
      coreX.resize (nSpecies);

      // Core energy
      coreEnergy.resize (nSpecies);
      // Initial guess data
      foccInit.resize (nSpecies);
      DijInit.resize (nSpecies);

      coreWaveFile.resize (nSpecies);

      if (pawGroup->contains ("kjxc")
          && pawGroup->get("kjxc")->toAttribute ())  {
         kjXcShape.resize (nSpecies);
      }

      for (species  = pawGroup->getGroup("species");
           species != NULL;
           species  = species->nextSibling("species"))
      {
         potFile = species->get("potential")->toString();
         potType = species->get("potType")->toString();

         if ( potType == "AbInit" ) {
            readAbInit (potFile, is);
         } else if (potType == "AtomPAW" ) {
            readAtomPAW (potFile, is);
         } else if (potType == "CPPAW" ) {
            readCPPAW (potFile, is);
         } else if (potType == "VASP" ) {
            bool useProjG = species->contains("useProjG")
                          ? species->get("useProjG")->toAttribute ()
                          : false;
            readVasp (potFile, is, useProjG);
         } else {
            cout << "EXITING: potType undefined\n";
            SX_EXIT;
         }
     
         if (species->contains ("lMaxRho"))
            lMaxRho(is) = species->get("lMaxRho")->toInt ();
         else
            lMaxRho(is) = 2 * lMax(is);
         if (verbose)
            cout << "lMaxRho(" << prettyName(is) << ")= " << lMaxRho(is)<< endl;

         if (species->contains ("angularGrid"))
            aGridType(is) = SxSphereGrid::GridType (species->get("angularGrid")
                                                    ->toInt ());
         else
            aGridType(is) = SxSphereGrid::Grid_110; // TODO: better default!!

         if (species->contains ("atomicRhoOcc"))  {
            foccInit(is) = species->get("atomicRhoOcc")->toList ();
            if (foccInit(is).getSize () < lPhi(is).getSize ())  {
               foccInit(is).resize (lPhi(is).getSize (), true);
            }
            // --- delete potential-defined valence density (if any)
            if (rhoInit(is).getSize () > 0)  {
               cout << "WARNING: overriding potential-defined initialisation "
                       "density for " << prettyName(is) << endl;
            }
            rhoInit(is) = SxDiracVec<Double> ();
            rhoInitBasis(is) = SxPtr<SxBasis> ();
         } else if (foccInit(is).getSize () == 0) {
            foccInit(is).resize (lPhi(is).getSize ());
            double n = valenceCharge(is);
            for (int ipt = 0; ipt < lPhi(is).getSize (); ++ipt)  {
               int l = lPhi(is)(ipt);
               double f = min (n, 2.  * (2 * l + 1));
               foccInit(is)(ipt) = f;
               n -= f;
               if (n < 0.) n = 0.;
            }
         }

         if (species->contains ("nRadGrid"))  {
            int nrNew = species->get("nRadGrid")->toInt ();
            double r0 = rad(is)(0), rmax = rad(is)(rad(is).getSize () -1);
            double logDrNew = log(rmax/r0) / (nrNew - 1);
            SxDiracVec<Double> newRad(nrNew);
            for (int i = 0; i < nrNew; ++i)
               newRad(i) = r0 * exp(logDrNew * i);
            if (rhoInit(is).getSize () == vBar(is).getSize ()
                && !dynamic_cast<SxRadGBasis*> (rhoInitBasis(is).getPtr ()))
            {
               SX_CHECK (!rhoInitBasis(is));
               rhoInit(is) = interpolateRad (rhoInit(is),r0,logDr(is),newRad);
            }
            rhoCorePS(is) = interpolateRad (rhoCorePS(is),r0,logDr(is),newRad);
            rhoCoreAE(is) = interpolateRad (rhoCoreAE(is),r0,logDr(is),newRad);
            vBar(is) = interpolateRad (vBar(is),r0,logDr(is),newRad);
            int npt = getNProjType(is);
            SxDiracMat<Double> newPhiAE(nrNew, npt),
                               newPhiPS(nrNew, npt),
                               newProj(nrNew, npt);
            newPhiAE.handle->auxData = phiAE(is).handle->auxData;
            newPhiPS.handle->auxData = phiPS(is).handle->auxData;
            newProj.handle->auxData = pPS(is).handle->auxData;
            for (int ipt = 0; ipt < npt; ++ipt)  {
               newPhiAE.colRef(ipt) <<= interpolateRad (phiAE(is).colRef(ipt),
                                                        r0, logDr(is), newRad);
               newPhiPS.colRef(ipt) <<= interpolateRad (phiPS(is).colRef(ipt),
                                                        r0, logDr(is), newRad);
               if (!dynamic_cast<const SxRadGBasis *>(pPS(is).getBasisPtr ()))
                  newProj.colRef(ipt)  <<= interpolateRad (pPS(is).colRef(ipt),
                                                        r0, logDr(is), newRad);
            }
            if (kjXcShape.getSize () > 0 && kjXcShape(is).getSize () > 0)  {
               kjXcShape(is) = interpolateRad (kjXcShape(is), 
                                               r0, logDr(is), newRad);
            }
            phiAE(is) = newPhiAE;
            phiPS(is) = newPhiPS;
            if (!dynamic_cast<const SxRadGBasis *>(pPS(is).getBasisPtr ()))
               pPS(is)   = newProj;
            rad(is)   = newRad;
            logDr(is) = logDrNew;
         }
         phiAE(is).handle->auxData.is = is;
         phiPS(is).handle->auxData.is = is;
         pPS(is).handle->auxData.is   = is;

         // --- debug output
         if (verbose && SxLoopMPI::me() == 0)  {
            for (int i = 0; i < lPhi(is).getSize (); ++i)  {
               SxString id = chemName(is) + "-" + i;
               SxBinIO io ("paw-ae-"+id+".dat", SxBinIO::ASCII_WRITE_ONLY);
               io.writeXYPlot (rad(is), phiAE(is).colRef(i));
               io.close ();
               io.open ("paw-ps-"+id+".dat", SxBinIO::ASCII_WRITE_ONLY);
               io.writeXYPlot (rad(is), phiPS(is).colRef(i));
               io.close ();
               if (!dynamic_cast<const SxRadGBasis *>(pPS(is).getBasisPtr ())) {
                  io.open ("paw-proj-"+id+".dat", SxBinIO::ASCII_WRITE_ONLY);
                  io.writeXYPlot (rad(is), pPS(is).colRef(i));
                  io.close ();
               }
            }
         }

         if (species->contains ("coreX"))  {
            coreX(is) = SxDiracMat<Double> (species->get("coreX")->toList ());
            int npt = (int)lPhi(is).getSize ();
            if (coreX(is).getSize () != npt * npt)  {
               cout << "Invalid size of core-exchange matrix for "
                    << prettyName(is) << endl;
               cout << "Must be " << npt << "x" << npt << ", but has only "
                    << coreX(is).getSize () << " elements.";
               SX_QUIT;
            }
            coreX(is).reshape (npt, npt);
         }

         if (species->contains ("coreWaves"))  {
            coreWaveFile(is) = species->get ("coreWaves")->toString ();
         }

         bool crashOnError = !species->contains("checkOverlap")
                           || species->get("checkOverlap")->toAttribute ();
         checkOverlap(is, crashOnError);

         ++is;
      }

   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   offset.resize (nSpecies);
   for (is = 0; is < nSpecies; ++is)  {
      int npt = (int)lPhi(is).getSize ();
      offset(is).resize (npt);
      offset(is)(0) = 0;
      for (int ipt = 1; ipt < npt; ++ipt)  {
         offset(is)(ipt) = offset(is)(ipt-1) + 2 * lPhi(is)(ipt-1) + 1;
      }
      SX_CHECK (offset(is)(npt-1)+2*lPhi(is)(npt-1) + 1 == getNProj(is),
                offset(is)(npt-1)+2*lPhi(is)(npt-1) + 1, getNProj(is));
      if (verbose)
         cout << "offset(is=" << is << "): " << offset(is) << endl;
   }

   int maxLRho = lMaxRho.maxval (),
       maxL    = max(lMax.maxval (), maxLRho);
   // clebschGordan in computeU needs highest lmax
   //clebschGordan = getClebschGordan (max(maxL, 2*maxLRho), maxL, maxLRho);
   // --- for forces higher cg coefficients needed
   clebschGordan = SxYlm::getClebschGordan (max(maxL+1,2*maxLRho+1), 
                                            2*maxL, 
                                            max(2*maxL+1,maxLRho));

   // set Radbasis
   setBasis (SxPtr<SxRadBasis>::create(rad, logDr));

   recomputeOverlap ();

   if (verbose)  {
      // --- print kinetic energy matrix elements 
      for (is = 0; is < nSpecies; ++is)  {
         for (int ipt = 0; ipt < lPhi(is).getSize (); ++ipt)  {
            int l = lPhi(is)(ipt);
            for (int jpt = 0; jpt < lPhi(is).getSize (); ++jpt)
               if (lPhi(is)(ipt) == lPhi(is)(jpt))  {
                  cout << "<phi" << ipt << "|T|phi" << jpt << ">="
                       << ((phiPS(is).colRef(ipt) 
                           * SxRadialAtom::laplace (phiPS(is).colRef(jpt), l)
                           -phiAE(is).colRef(ipt) 
                           * SxRadialAtom::laplace (phiAE(is).colRef(jpt), l)
                           )
                           * rad(is).cub ()).integrate (logDr(is)) / 2.
                       << endl;
               }
         }
      }
   }

   // --- compute core energy
   for (is = 0; is < nSpecies; ++is)  {

      coreEnergy(is) = 0.;
      // Hartree part
      SxDiracVec<Double> vH;
      vH = SxRadialAtom::getHartreePotential(rhoCoreAE(is));
      coreEnergy(is) += 0.5 * (rhoCoreAE(is) * vH * rad(is).cub ())
                        .integrate (logDr(is));
      if (verbose) sxprintf ("eCore(Hartree)=%.12f\n", coreEnergy(is));
      // xc part
      double eXc = 0.;
      SxRadialAtom::computeXC(rad(is), logDr(is), 
                                  SQRT_1_4PI * rhoCoreAE(is),
                                  xc, &eXc);
      if (verbose) sxprintf ("eXc(core)=%.12f\n", eXc);
      coreEnergy(is) += eXc;
      // nuclear part
      double eCoreNuc = -sqrt(FOUR_PI) * nuclearCharge(is) 
                        * (rhoCoreAE(is) * rad(is).sqr ())
                          .integrate (logDr(is));
      // --- Estimate integral 0..r0 from rho(r0) and rho'(r0)
      {
         double r0  = rad(is)(0),
                r0cub = r0*r0*r0,
                r01 = r0 - rad(is)(1),
                rho0 = rhoCoreAE(is)(0),
                rho1 = rhoCoreAE(is)(1),
                Z    = nuclearCharge(is);
         eCoreNuc -= sqrt(FOUR_PI)* Z
                   * (rho0 * 0.5 * r0 * r0 - (rho0-rho1) / r01 * r0cub/6.);
      }
      if (verbose) sxprintf ("eCoreNuc=%.12f\n", eCoreNuc);
      coreEnergy(is) += eCoreNuc;

      // --- set up core-valence exchange
      if (coreWaveFile(is).getSize () > 0)  {

         if (SxFileParser (coreWaveFile(is)).reads ("All-electron core wave")) {
            readCoreAbinit (coreWaveFile(is), is);
         } else {
            SxAtomicOrbitals orbs;
            orbs.setBasis (radBasisPtr);
            // --- read orbitals (in atompaw wfn format)
            FILE *fp = fopen (coreWaveFile(is).ascii (), "r");
            if (!fp)  {
               cout << "Cannot open " << coreWaveFile(is) << endl;
               SX_EXIT;
            }
            cout << "Reading core waves from " << coreWaveFile(is) << endl;
            orbs.readOrbitals (fp, is);
            fclose (fp);

            // --- save into phiCoreAE, set up lCore
            int nCore = orbs.getNOrbTypes (is);
            lCore(is).resize (nCore);
            psiCoreAE(is).reformat (rad(is).getSize (), nCore);
            psiCoreAE(is).setBasis (*radBasisPtr);
            psiCoreAE(is).handle->auxData.is = is;
            for (int iCore = 0; iCore < nCore; ++iCore)  {
               psiCoreAE(is).colRef (iCore) << (orbs(is, iCore) / rad(is));
               lCore(is)(iCore) = orbs.getL(is, iCore);
            }
         }

         // --- check core
         SxDiracVec<Double> rhoCore (rad(is).getSize ());
         rhoCore.set (0.);
         for (int ic = 0; ic < lCore(is).getSize (); ++ic)  {
            rhoCore.plus_assign_ax (2. * (2 * lCore(is)(ic) + 1),
                                    psiCoreAE(is).colRef (ic).sqr ());
         }
         cout << "Core from waves: " << (rhoCore * rad(is).cub ())
                                      .integrate (logDr(is)) /* * FOUR_PI */ << endl;
         cout << "Core from file: "  << (rhoCoreAE(is) * rad(is).cub ())
                                      .integrate (logDr(is)) * sqrt(FOUR_PI)
                                   << endl;
         
         writePlot ("rhoCoreWave-" + prettyName(is) + ".dat",
                     rad(is), rhoCoreAE(is)*sqrt(FOUR_PI), rhoCore);

      }

      if (hybridFunctional && psiCoreAE(is).getSize () == 0)  {
         int npt = (int)lPhi(is).getSize ();
         coreX(is) = SxDiracMat<Double> (npt, npt);
         double nCore = sqrt(FOUR_PI) 
                      * (rhoCoreAE(is) * rad(is).cub ()).integrate (logDr(is));
         if (nCore > 1e-6)  {
            cout << SX_SEPARATOR;
            cout << "| WARNING" << endl; 
            cout << "| substituting exact core-valence exchange with PBE for "
                 << prettyName(is) << endl;
            cout << SX_SEPARATOR;

            // --- compute core exchange potential
            SxDiracVec<Double> vXCore 
               = SxRadialAtom::computeXC(rad(is), logDr(is), 
                                         SQRT_1_4PI * rhoCoreAE(is), SxXC::PBE,
                                         NULL, SxXCFunctional::ComputeX);
            // --- compute vXCore matrix elements
            for (int ipt = 0; ipt < npt; ++ipt)  {
               for (int jpt = 0; jpt < npt; ++jpt)  {
                  if (lPhi(is)(ipt) == lPhi(is)(jpt))  {
                     coreX(is)(ipt, jpt) = tr(nij(phiAE(is),ipt,jpt) * vXCore);
                  } else {
                     coreX(is)(ipt, jpt) = 0.;
                  }
               }
            }
         } else {
            coreX(is).set (0.);
         }
      }

      if (verbose) {
         // --- write vBar
         writePlot ("vbar-" + prettyName(is) + ".dat", rad(is), vBar(is));
         writePlot ("rhoCore-" + prettyName(is) + ".dat", rad(is),
                     rhoCoreAE(is), rhoCorePS(is));
      }
   }
   setupQijL ();
   print ();
}
/*
// catch exceptions from initializer
catch (SxException e)  {
   e.print ();
   SX_EXIT;
}
*/


SxPAWPot::~SxPAWPot ()
{
   // empty
}

SxDiracVec<Double> getJsbShape (double rCut, int l,
                                const SxDiracVec<Double> rad, double logdr)
{
   double x0[2], start = 0.;
   // search for zero of spherical Bessel functions
   for (int iq = 0; iq < 2; ++iq)  {
      double step = 0.01;
      x0[iq] = start + step;
      double jLow = SxYlm::jsb(l, x0[iq]), jHigh;
      while ((jHigh = SxYlm::jsb(l, x0[iq] + step)) * jLow > 0.)  {
         jLow = jHigh;
         x0[iq] += step;
      }
      // interval search
      while (step > 1e-12)  {
         step *= 0.5;
         double j = SxYlm::jsb(l, x0[iq] + step);
         if (fabs(j) < 1e-12) {
            x0[iq] += step;
            break;
         }
         if (j * jLow > 0.)  {
            x0[iq] += step;
            jLow = j;
         } else {
            SX_CHECK (j * jHigh >= 0.);
            jHigh = j;
         }
      }
      start = x0[iq];
   }

   SxMatrix<Double> A(2,2);
   int nr = (int)rad.getSize ();
   SxArray<SxDiracVec<Double> > jsbRad(2);
   for (int iq = 0; iq < 2; ++iq)  {
      double q = x0[iq] / rCut;
      jsbRad(iq).resize (nr);
      for (int ir = 0; ir < nr; ++ir)
         jsbRad(iq)(ir) = (rad(ir) < rCut) ?  SxYlm::jsb(l, q * rad(ir)) : 0.;

      A(0, iq) = -q * SxYlm::jsb(l+1, x0[iq]); // jsb'
      A(1, iq) = (jsbRad(iq) * pow(rad, l+3)).integrate (logdr);
   }
   SxMatrix<Double> Ainv = A.inverse ();
   return Ainv(0, 1) * jsbRad(0) + Ainv(1, 1) * jsbRad(1);
}

SxDiracVec<Double> getFunction (const SxString &lines, bool verbose, 
                                int nPts, int nPtsMax,
                                double logDr, bool expCont = false)
{
   SxList<SxString> phiList = lines.tokenize ('\n');

   SxDiracVec<Double> rOut(nPtsMax);
   rOut.set (0.);

   SxList<SxString> rLine;
   SxList<SxString>::Iterator itFu = phiList.begin ();

   int lastFullL=nPts-(nPts%3)-3;
   SxList<SxString> rList;

   {
      double x0,x1,x2;
      while ( sscanf(itFu->ascii (), "%lf%lf%lf", &x0, &x1, &x2) != 3)
         itFu++;
   }

   for (int i=0;i<=lastFullL;i+=3) {
      rList=(*itFu).tokenize (' ');
      rOut(i)=rList(0).toDouble ();
      rOut(i+1)=rList(1).toDouble ();
      rOut(i+2)=rList(2).toDouble ();
      itFu++;
   }

   if ( (nPts%3) != 0 ) {
      rList=(*itFu).tokenize (' ');
      if ( (nPts%3) == 1 ) rOut(nPts-1)=rList(0).toDouble ();
      if ( (nPts%3) == 2) {
         rOut(nPts-2)=rList(0).toDouble ();
         rOut(nPts-1)=rList(1).toDouble ();
      }
   }

   if (nPts < nPtsMax)  {
      if (expCont && fabs(rOut(nPts - 1)) > 1e-6)  {
         if (fabs(rOut(nPts -1)) < fabs(rOut(nPts-2))
             && rOut(nPts -1) * rOut(nPts-2) > 0.)
         {
            // exponential continuation
            // this works also for the shifted logarithmic grid of atompaw
            double bp = -log(rOut(nPts - 1)/rOut(nPts - 2))
                      / (exp((nPts-1)*logDr) - exp((nPts-2)*logDr) );
            double ap = rOut(nPts - 1) / exp (-bp * exp(logDr * (nPts-1)));
            if (verbose)  {
               cout << "ap=" << ap << "; bp = " << bp << endl;
               cout << rOut(nPts-2) << endl;
               cout << ap * exp (-bp * exp(logDr * (nPts-2))) << endl;
               cout << rOut(nPts-1) << endl;
               cout << ap * exp (-bp * exp(logDr * (nPts-1))) << endl;
            }
            for (int ir = nPts; ir < nPtsMax; ++ir)  {
               rOut(ir) = ap * exp (-bp * exp(logDr * ir));
            }
         } else {
            // r(n-2) / r0 
            double rPeak = exp((nPts- 2)*logDr);
            double s = exp(logDr);
            double beta = 0.3 * rPeak;
            double b2 = beta * beta;
            // df/dr / r0
            double fp = (rOut(nPts - 1) - rOut(nPts - 3))
                      / (rPeak * s - rPeak / s);
            if (verbose) cout << "fp=" << fp << endl;
            double a = rOut(nPts - 2);
            double b = fp;
            if (verbose) 
               cout << "a=" << a << "; beta=" << beta << "; b = " << b << endl;
            for (int ir = nPts; ir < nPtsMax; ++ir)  {
               double r = exp(ir*logDr);
               rOut(ir) = (a + b * (r-rPeak)) * exp (-sqr(r-rPeak)/b2);
            }
         }
      }
   }

   // --- interpolation to standard logarithmic grid
   SxNaturalCubicSpline spline(rOut, true);
   // present indices correspond to r0 * (exp(logDr * i) - 1)
   // change to r0 * exp(logDr * i)
   int n = (int)rOut.getSize () - 1;
   for (int i = 0; i < n; ++i)  {
      rOut(i) = spline.getVal(log(1. + exp(logDr * i))/logDr);
   }
   // last point is extrapolated
   rOut(n) = spline.getValYExtra(log(1. + exp(logDr * n))/logDr);

   return rOut;

}

SxDiracMat<Double> readWaves (SxList<SxString>::Iterator &it,
                              const SxList<SxString>::Iterator &end,
                              const SxString &marker,
                              const SxString &stopMarker,
                              bool verbose,
                              int nPtsMax, int nWaves,
                              int meshPts,
                              double logdr,
                              bool expCont = false)
{
   SxDiracMat<Double> res(nPtsMax, nWaves);
   SxList<SxString> valList, cLine;

   while ( (*it).contains(marker) == 0 ) it++;

   int i=0;
   for (; it != end && (*it).contains (stopMarker) == 0; it++) {
      if (i >= nWaves)  {
         cout << "Encountered unexpected part of file" << endl;
         SX_EXIT;
      }
      res.colRef (i) <<= getFunction (*it, verbose, meshPts, nPtsMax, logdr, expCont);
      i++;
   }
   
   VALIDATE_VECTOR(res);
   return res;
}


void SxPAWPot::readAbInit (const SxString &filename, int iSpecies)
//TODO: simplify regexpr by a more pragmatic approach (see readAtomPAW)
{

   SxString potFile;
   SxList<SxString> cLine;
   SxXC::XCFunctional xcFunctional;
   
   if (verbose)
      cout << "reading " << filename << " ...\n";

   try {
      potFile = SxFileIO::readBinary (filename,-1);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   // --- splitting of input + joining HEAD marker with corresponding data
   SxList<SxString> potFileTok =  potFile.substitute("=====","&").tokenize('&');
   for (int i=1;;i++) {

      if (i >= potFileTok.getSize ()) break;

      potFileTok(i) += potFileTok(i+1);
      potFileTok.remove(i+1);

    }

   SxList<SxString> valList;

   // check input format
   int pawFormat = -1;
   if( potFileTok(0).contains ("paw4") == 0 ) pawFormat = 4;
   if( potFileTok(0).contains ("paw5") == 0 ) pawFormat = 5;
   if (pawFormat != 4 && pawFormat != 5)
   {
      cout << "\nEXITING: Cannot read potential " << filename << endl;
      cout << "Format is ";
      int c = 0;
      // skip three lines
      for (int i = 0; i < 3 && c < potFile.getSize (); ++i)
         for (; c < potFile.getSize (); c++) if (potFile(c) == '\n') break;
      // skip white space
      for (c++; c < potFile.getSize (); c++) if (potFile(c) != ' ') break;
      for (; c < potFile.getSize (); c++) {
         if (potFile(c) == ' ') break;
         cout << potFile(c);
      }
      cout << endl << "For this routine, paw format must be 4 or 5\n" << endl;
      SX_EXIT;
   }

   // --- fetch core charge
   cLine = SxRegex("\\s+"+realNum+"\\s+"+realNum+".*: zatom,zion,pspdat")
           .match (potFileTok(0));
   double coreCharge = cLine(1).toDouble () - cLine(2).toDouble ();
   nuclearCharge(iSpecies) = cLine(1).toDouble ();
   valenceCharge(iSpecies) = cLine(2).toDouble ();
   if (verbose)  {
      cout << "Z=" << nuclearCharge(iSpecies) << endl;
      cout << "nv=" << valenceCharge(iSpecies) << endl;
   }

   // --- fetch XC-Functional, lMax
   cLine = SxRegex("\\s+"+intNum+"\\s+" +intNum+"\\s+"+intNum+"\\s+"
                   +intNum+"\\s+"+intNum+"\\s+"+realNum+
                   "\\s+: pspcod,pspxc,lmax,lloc,mmax,r2well")
           .match (potFileTok(0));

   // XC-Functional
   if ( cLine(2) == 7 )  {
      if (verbose) cout << "LDA-PW" << endl;
      xcFunctional = SxXC::LDA_PW;
   } else if ( cLine(2) == 11 || cLine(2) == -101130 /* libXc */)  {
      if (verbose) cout << "GGA-PBE" << endl;
      xcFunctional = SxXC::PBE;
   } else {
      cout << "Unknown xc functional " << cLine(2) 
           << " in " << filename << endl;
      SX_EXIT;
   }

   // lmax
   lMax(iSpecies) = cLine(3).toInt ();
   if (verbose) cout << "lMax: " << lMax(iSpecies) << endl;

   // number of AE/PS/Projectors
   cLine = SxRegex("\\s+"+intNum+"\\s+"+intNum+"\\s+ : basis_size,lmn_size")
           .match (potFileTok(0));
   int nWaves = cLine(1).toInt ();

   // --- setup vector with meshIdx-info
   int nIdx = 3*nWaves+4, nIt=0;
   SxDiracVec<Int> meshIdxs(nIdx);
   bool KresseJoubertXC = false;
   for (int i=0; i<potFileTok.getSize(); i++) {

       if ( potFileTok(i).contains ("VHntZC") == 1) {
          cLine = SxRegex("\\s+"+intNum+"\\s+"+intNum+
                          "\\s+\\: radial mesh index,*").match (potFileTok(i));
          meshIdxs(nIt) = cLine(1).toInt ();
          if ( cLine(2).toInt () == 0 ) {
             cout << "VHntZC format 0 (vBare) is not tested." << endl;
             SX_EXIT;
          }
          if ( cLine(2).toInt () == 1 ) KresseJoubertXC = true;
          nIt++;
      }
      else if ( potFileTok(i).contains("radial mesh index") ) {
         cLine = SxRegex("\\s+"+intNum+"\\s+\\: radial mesh index")
                 .match (potFileTok(i));
         meshIdxs(nIt) = cLine(1).toInt ();
         nIt++;  
      }      

   }

   // --- fetch ordering of wavefunctions with respect to l
   SxList<SxString>::Iterator it;
   SxList<SxString> potFileTokLine = potFileTok(0).tokenize ('\n');

   SxArray<int> lOrdering(nWaves);
   for (it = potFileTokLine.begin();it != potFileTokLine.end();it++) {
      if ( (*it).contains("orbitals") == 1 ) {
         cLine = (((*it).left(':')).tokenize(' '));
         for (int i = 0; i < cLine.getSize(); i++) 
            lOrdering(i) = cLine(i).toInt ();   
         break;   
      }   
   }

   lPhi(iSpecies) = lOrdering;
   if (verbose) cout << "lOrdering: " << lOrdering << endl;

   // --- fetch radial grid parameters
   cLine = SxRegex ("\\s+"+intNum+"\\s+: number_of_meshes")
           .match (potFileTok(0));
   int nMeshes = cLine(1).toInt ();
   SxVector<Int> meshPts(nMeshes);
   
   int i = 0, nPtsMax = -1;
   double r0old = 0., logDrOld = 0.;
 
   for (;it != potFileTokLine.end (); it++) {
      if ( (*it).contains ("rad_step") == 1 ) {
         double r0;
         cLine = SxRegex ("^\\s+"+intNum+"\\s+"+intNum + "\\s+" 
                                +intNum+"\\s+"+sciNum+"\\s+"
                                +sciNum+" : mesh*").match (*it);         
         if ( cLine(2) == 2 ) {
            meshPts(i) = (cLine(3).toInt ());
            r0 = cLine(4).toDouble ();
            logDr(iSpecies) = cLine(5).toDouble ();
            if ( (r0 != r0old) && (r0old != 0.) ) {
               cout << "EXITING: r0 is not the same for all grids!" << endl;
               SX_EXIT;
            }
            if ( (logDr(iSpecies) != logDrOld) && (logDrOld != 0.) ) {
               cout << "EXITING: logDr is not the same for all grids!" << endl;
               SX_EXIT;
            }
            i++;

            // --- compute radial grid for rMax
            if ( i == nMeshes ) {
               nPtsMax = meshPts.maxval ();
               SxDiracVec<Double> rMesh(nPtsMax);
               for (int j = 0; j < nPtsMax; j++) 
                  rMesh(j) = r0*(exp (j*logDr(iSpecies)));
               rad(iSpecies) = rMesh;
               break;
            }

            r0old = r0;
            logDrOld = logDr(iSpecies);

         } else {
            cout << "EXITING: Not a logarithmic grid!" << endl;
            SX_EXIT;
         }
      }
   }

   if (verbose)  {
      cout << "meshPts:" << meshPts << endl;
      cout << "logDr: " << logDr(iSpecies) << endl;
   }
   double logdr = logDr(iSpecies);
   int nr = (int)rad(iSpecies).getSize ();
   
   // hack SxRadBasis for getHartreePotential
   SxRadBasis radB;
   radB.radFunc.resize (iSpecies+1);
   radB.radFunc(iSpecies) = rad(iSpecies);
   radB.logDr.resize (iSpecies+1);
   radB.logDr(iSpecies) = logdr;

   //r_cut(PAW)
   cLine = SxRegex("\\s+"+realNum+"\\s+: r_cut\\(PAW\\)").match (potFileTok(0));
   double rPAW = cLine(1).toDouble ();
   
    // --- fetch shape function parameters
   // TODO: compute rc from rpaw for bessel
   // normalized compensation charge
   SxDiracVec<Double> g0File; 
   for (; it != potFileTokLine.end(); it++) {
      if ( (*it).contains ("shape_type") == 1 ) {
         cLine = (((*it).left (':')).tokenize (' '));
         switch (cLine(0).toInt ()) {
            case 1:
               if (verbose)
                  cout << "Gaussian shape function [exp(-(arg/sigma)^lambda)]" 
                       << endl;
               if ( cLine(2).toInt () != 2 ) {
   cout << "EXITING: only Gaussian shape functions are supported (lambda=2)"
        << endl;
                  SX_EXIT;
               }
               // sigma corresponds to rshape
               rc(iSpecies) = cLine(3).toDouble ();
               g0File = exp(-(rad(iSpecies)/rc(iSpecies)).sqr ());
               break;
            case 2:
               // --- fetch rShape
               double rShape;
               if (verbose)  {
                  cout << "Sinc shape function" <<
                          "[(sin(pi*arg/rc_hat)/(pi*arg/rc_hat))**2]\n" <<
                          "Assuming k(rShape)=1e-4" << endl;
               }

               if ( cLine(1).toDouble () < 1e-5 ) rShape = rPAW;
               else rShape = cLine(1).toDouble ();
               
               rc(iSpecies) = rShape/sqrt (-log(1e-4));
               g0File.resize (nr);
               for (int ir = 0; ir < nr; ++ir)  {
                  double r = rad(iSpecies)(ir);
                  if (r < rShape)  {
                     double x = PI * r / rShape;
                     g0File(ir) = sqr(sin(x)/x);
                  } else {
                     g0File(ir) = 0.;
                  }
               }
               break;
             case 3:            
               cout << "Bessel shape function [al1.jl(ql1.r)+al2.jl(ql2.r)]\n"
                    << "Yet to come..." << endl;              
               SX_EXIT
             default:
                  SX_QUIT;
         }
       }
   }
   // normalize compensation charge
   g0File /= FOUR_PI * (g0File * rad(iSpecies).cub ()).integrate(logdr);

   if (verbose) cout << "rc: " << rc(iSpecies) << endl;
   
   // --- fetch AE partial waves 
   int itMeshIdx = 0;  
   int meshIdx = meshIdxs(itMeshIdx);
   phiAE(iSpecies) = readWaves (it = potFileTok.begin (),
                                potFileTok.end (), "PHI", "TPHI", verbose,
                                nPtsMax, nWaves, meshPts(meshIdx-1), 
                                logdr, true);

   for (int iWave = 0; iWave < nWaves; ++iWave)
      phiAE(iSpecies).colRef (iWave) *= 1./rad(iSpecies);
   
   // --- check normalization 
   if (verbose) {
      cout << "Normalization <phiAE|phiAE>" << endl;

      for (int iWave =0;iWave < nWaves; iWave++) 
         cout << " i=" << iWave << ": " 
              << (phiAE(iSpecies).colRef(iWave).sqr ()*rad(iSpecies).cub ())
                 .integrate(logdr) << endl;
   }

   // --- fetch PS partial waves
   itMeshIdx += nWaves;
   meshIdx = meshIdxs(itMeshIdx);
   phiPS(iSpecies) = readWaves (it, potFileTok.end (), "TPHI",
                                "TPROJECTOR", verbose, nPtsMax, nWaves, 
                                meshPts(meshIdx-1),logdr, true);
   for (int iWave = 0; iWave < nWaves; ++iWave)  {
      phiPS(iSpecies).colRef (iWave) *= 1./rad(iSpecies);
      // --- enforce phiAE=phiPS outside cutoff region
      int ir;
      for (ir = 0; ir < nr; ++ir)  {

         if (rad(iSpecies)(ir) > rPAW)
           phiAE(iSpecies)(ir, iWave) = phiPS(iSpecies)(ir, iWave);
      }
   }

   // --- fetch projector functions
   itMeshIdx += nWaves;
   meshIdx = meshIdxs(itMeshIdx);
   pPS(iSpecies) = readWaves (it, potFileTok.end (), "TPROJECTOR",
                              "CORE_DENSITY", verbose, nPtsMax, nWaves,
                              meshPts(meshIdx-1),logdr);
   for (int iWave = 0; iWave < nWaves; ++iWave)
      pPS(iSpecies).colRef (iWave) *= 1./rad(iSpecies);

   // --- check orthonormality condition 
   if (verbose)  {
      cout << "Orthonormality conditions <phiPS|pPS>" << endl;
      for (int ipt=0;ipt<nWaves;ipt++)  {
         for (int jpt=0; jpt<nWaves; jpt++)  { 
            if (lPhi(iSpecies)(ipt) == lPhi(iSpecies)(jpt))  {
               cout << " i=" << ipt << "/j=" << jpt << ": "
                    << (phiPS(iSpecies).colRef(ipt)*pPS(iSpecies).colRef(jpt)
                        *rad(iSpecies).cub ()).integrate(logDr(iSpecies)) << endl;
            }
         }
      }
   }

   // --- fetch AE core density
   itMeshIdx += nWaves;
   meshIdx = meshIdxs(itMeshIdx);
   rhoCoreAE(iSpecies) = getFunction (*it, verbose, meshPts(meshIdx-1),
                                      nPtsMax, logdr);
   rhoCoreAE(iSpecies).handle->auxData.is = iSpecies;
   rhoCoreAE(iSpecies).handle->auxData.l = 0;
   rhoCoreAE(iSpecies).handle->auxData.m = 0;

   if (verbose)  {
      cout << "core charge: " 
           << "density is defined as sum(|phi|^2/4*PI*r^2): "
           << FOUR_PI*(
                 //rad(iSpecies) * (rad(iSpecies) - rad(iSpecies)(0)).sqr ()
                 rad(iSpecies).cub ()
                 *rhoCoreAE(iSpecies)
                 ).integrate (logDr(iSpecies)) << endl;

      cout << "N (4*Pi*N*rhoCoreAE=corecharge): \n" 
           << coreCharge / 
           (FOUR_PI*(
                     (rad(iSpecies).cub ()*rhoCoreAE(iSpecies)
                      ).integrate (logDr(iSpecies)))) << endl;
   }

   // fetch PS core density
   it++; itMeshIdx++;
   meshIdx = meshIdxs(itMeshIdx);
   rhoCorePS(iSpecies) = getFunction (*it, verbose, meshPts(meshIdx-1),
                                      nPtsMax, logdr);
   rhoCorePS(iSpecies).handle->auxData.is = iSpecies;
   rhoCorePS(iSpecies).handle->auxData.l = 0;
   rhoCorePS(iSpecies).handle->auxData.m = 0;

   // fetch "frozen" values of Dij = Dij0 (?)
   // Part of the Dij term calculated in the psp part
   it++;
   SxList<SxString> matList = (*it).tokenize ('\n');
   matList.removeFirst ();
   int matSize = (int)matList.getSize ();

   SxDiracMat<Double> Dij(matSize,matSize);
   if (verbose) cout << "matSize=" << matSize << endl;

   for (int a = 0; a < matSize; a++) {
      cLine = matList(a).tokenize (' ');
      for (int b = 0; b<cLine.getSize(); b++) {
         Dij(a,b) = cLine(b).toDouble ();
         Dij(b,a) = Dij(a,b);
      }
   }
   if (verbose) cout << "Dij:\n" << Dij << endl; 

   // rhoij0= Atomic initialization of rhoij
   it++;
   matList = (*it).tokenize ('\n');
   matList.removeFirst ();
   matSize = (int)matList.getSize ();
   SxDiracMat<Double> rhoij0(matSize,matSize);
   for (int a = 0;a < matSize; a++) {
      cLine = matList(a).tokenize (' ');
      for (int b = 0; b < cLine.getSize (); b++) {
         rhoij0(a,b) = cLine(b).toDouble ();
         rhoij0(b,a) = rhoij0(a,b);
      }
   }
   if (verbose) cout << "rhoij0:\n" << rhoij0 << endl;

   // --- vbar
   it++; itMeshIdx++;
   meshIdx = meshIdxs(itMeshIdx);
   vBar(iSpecies) = getFunction (*it, verbose, meshPts(meshIdx-1), nPtsMax,
                                 logdr);
   for (ssize_t ir = meshPts(meshIdx-1); ir < nPtsMax; ++ir)
      vBar(iSpecies)(ir) = -valenceCharge(iSpecies) / rad(iSpecies)(ir);
   vBar(iSpecies).handle->auxData.is = iSpecies;
   vBar(iSpecies).handle->auxData.l = 0;
   vBar(iSpecies).handle->auxData.m = 0;

   // --- compute kinetic energy matrix elements
   {
      SxDiracVec<Double> rhoijPS, rhoijAE;
      rhoCoreAE(iSpecies).setBasis (&radB);
      SxDiracVec<Double> vZcAE = (-nuclearCharge(iSpecies))/rad(iSpecies)
             + SxRadialAtom::getHartreePotential (rhoCoreAE(iSpecies)),
                         vZcPS = vBar(iSpecies),
                         r3 = rad(iSpecies).cub ();
      // 0..r0 correction
      // e= 0.5 * Z * r0^2 * rho(r0) -  1/6  Z * r0^3 * rho'(r0)
      // integration weight @ r0=rad(0): logdr r0^3 * w_Simpson(0)
      double rhoPrimePot = nuclearCharge(iSpecies) / 6.
                        / (rad(iSpecies)(0) - rad(iSpecies)(1));
      rhoPrimePot=0.; // strangely, this seems to be better for abInit
      vZcAE(0) -= (0.5 * nuclearCharge(iSpecies)/rad(iSpecies)(0) - rhoPrimePot)
                / (logdr /3.);
      vZcAE(1) -= rhoPrimePot / (logdr * 4./3.);
      rhoCoreAE(iSpecies).handle->auxData.basisPtr = NULL;
      double Xij, Sij, Qij;
      SxVector<Int> map(nWaves);
      for (int ipt = 0, ip = 0; ipt < nWaves; ++ipt)  {
         int l = lPhi(iSpecies)(ipt);
         for (int m = -l; m <= l; ++m, ++ip)
            if (m == 0) map(ipt) = ip;
      }
      deltaKin(iSpecies) = SxDiracMat<Double> (nWaves, nWaves);
      deltaS(iSpecies)   = SxDiracMat<Double> (nWaves, nWaves);

      for (int ipt = 0; ipt < nWaves; ++ipt)  {
         for (int jpt = 0; jpt < nWaves; ++jpt)  {
            if (lPhi(iSpecies)(ipt) == lPhi(iSpecies)(jpt))  { 
               rhoijPS = nij(phiPS(iSpecies),ipt,jpt);
               rhoijAE = nij(phiAE(iSpecies),ipt,jpt);
               Qij = FOUR_PI * ((rhoijAE - rhoijPS) * r3).integrate (logdr);
               SX_CHECK_NUM (Qij);
               // Ref. 3, Eq. 98
               Xij = ((vZcAE * rhoijAE - vZcPS * rhoijPS) * r3)
                          .integrate (logdr);
               // Ref. 3, Eq. 99
               Sij = Qij * (vZcPS * g0File * r3).integrate(logdr);
               // Ref. 3, Eq. 97
               deltaKin(iSpecies)(ipt,jpt) = Dij(map(ipt),map(jpt)) - Xij + Sij;
               deltaS(iSpecies)(ipt,jpt)   = Qij / FOUR_PI;
            } else {
               deltaKin(iSpecies)(ipt,jpt) = 0.;
               deltaS(iSpecies)(ipt,jpt)   = 0.;
            }
         }
      }
      if (verbose)  {
         int oldPrec = (int)cout.precision ();
         cout.precision (12);
         cout << "deltaKin:" << endl;
         cout << deltaKin(iSpecies) << endl;
         cout << "deltaS:" << endl;
         cout << deltaS(iSpecies) << endl;
         cout.precision (oldPrec);
      }
   }

   // --- PS Valence density
   it++; itMeshIdx++;
   meshIdx = meshIdxs(itMeshIdx);
   SxDiracVec<Double> rhoValPS = 
      getFunction (*it, verbose, meshPts(meshIdx-1), nPtsMax, logdr);
   
   rhoInit(iSpecies) = (rhoValPS + rhoCorePS(iSpecies)) * sqrt(FOUR_PI);
   rhoInit(iSpecies).handle->auxData.is = iSpecies;
   rhoInit(iSpecies).handle->auxData.l  = 0;
   rhoInit(iSpecies).handle->auxData.m  = 0;
   SxDiracVec<Double> sumVec(nPtsMax);
   sumVec = 0.;
   
   SxVector<Int> map(getNProj(iSpecies));
   for (int ipt = 0, ip = 0; ipt < nWaves; ++ipt)  {
      int l = lPhi(iSpecies)(ipt);
      for (int m = -l; m <= l; ++m, ++ip)
         map(ip) = ipt;
   }
   foccInit(iSpecies).resize (nWaves);
   foccInit(iSpecies).set (0.);
   // --- Maxmal difference for pseudo valence densities
   for (int ia = 0; ia < rhoij0.nRows (); ia++) {
      for (int ja = 0; ja< rhoij0.nCols (); ja++) {
         //sumVec += rhoij0(ia,ja)*phiPS(iSpecies).colRef(map(ia))
         //          *phiPS(iSpecies).colRef(map(ja));
         sumVec += rhoij0(ia,ja)* nij(phiAE(iSpecies),map(ia),map(ja));
      }
      foccInit(iSpecies)(map(ia)) += rhoij0(ia,ia);
   }
   sumVec /= FOUR_PI;

   if (verbose) {
      int idx;
      double d = (rhoValPS-sumVec).maxval (&idx);
      cout << "Maximal difference for pseudo valence densities\n"
           << "(from File - from PS waves): " << d
           << "@" << idx << ": "
           << rhoValPS(idx) << "/" << sumVec(idx) << endl;
   }

   double SQRT_4PI = sqrt(FOUR_PI);
   {

      // Monopole moments of valence and core
      double Qv = foccInit(iSpecies).sum ()
                - FOUR_PI * (rhoValPS * rad(iSpecies).cub ()).integrate (logdr);
      cout << "Qv=" << Qv << endl;
      rhoInitQ0(iSpecies) = Qv;
      if (verbose) cout << "Qv = " << Qv << endl;

      double QZc = - nuclearCharge(iSpecies)
                 + FOUR_PI * ((rhoCoreAE(iSpecies) - rhoCorePS(iSpecies))
                              * rad(iSpecies).cub ()).integrate(logdr);
      
      SxDiracVec<Double> g0 = getGrl(iSpecies, 0) / FOUR_PI;
      //SxDiracVec<Double> g0 = g0File; // DEBUG only
      g0File.setBasis (&radB);
      g0File.handle->auxData.is = iSpecies;
      g0File.handle->auxData.l = 0;
      g0File.handle->auxData.m = 0;
      g0.setBasis (&radB);
      g0.handle->auxData.is = iSpecies;
      g0.handle->auxData.l = 0;
      g0.handle->auxData.m = 0;
      if (verbose)  {
         cout << "g0 moment: "
              << (FOUR_PI * (g0 * rad(iSpecies).cub ()).integrate (logdr)) << endl;
         cout << "g0File moment: "
              << (FOUR_PI * (g0File * rad(iSpecies).cub ()).integrate (logdr)) << endl;
         writePlot ("rhoVal-"+prettyName(iSpecies), rad(iSpecies),
                    rhoValPS + Qv * g0, rhoValPS, sumVec);
      }


      // pseudo-density
      SxDiracVec<Double> rhoPS = rhoCorePS(iSpecies) + rhoValPS;
      // xc and Hartree potentials
      SxDiracVec<Double> vXcKJ, vXcPS, vHdelta; 
      vXcKJ = SxRadialAtom::computeXC(rad(iSpecies), logdr,
                                          rhoPS + Qv * g0File,
                                          xcFunctional);
      vXcPS = SxRadialAtom::computeXC(rad(iSpecies), logdr,
                                          rhoPS,
                                          xcFunctional);
      // --- cf. Ref 3, Eq. (35) (Bloechl) vs Eq. 89 (Kresse-Joubert)
      vHdelta = SxRadialAtom::getHartreeSmooth (
                (rhoCorePS(iSpecies) - Qv * g0File));
      vHdelta += SxRadialAtom::getHartreePotential ((QZc + Qv) * g0);
      vBar(iSpecies) -= vHdelta;

      if (kjXcShape.getSize () == 0)  {
         if (KresseJoubertXC) vBar(iSpecies) +=  vXcKJ -  vXcPS;
      } else {
         if (!KresseJoubertXC) vBar(iSpecies) -=  vXcKJ -  vXcPS;
         kjXcShape(iSpecies) = FOUR_PI * g0File;
      }
      // writePlot ("vxcKJ-" + prettyName(iSpecies) + ".dat",
      //             rad(iSpecies), vXcKJ);
      // writePlot ("vxcBl-" + prettyName(iSpecies) + ".dat",
      //             rad(iSpecies), vXcPS);
   }
   // writePlot ("vbar-" + prettyName(iSpecies) + ".dat", rad(iSpecies),
   //             vBar(iSpecies));

   // switch to Y_00 expansion coefficient for potentials/densities
   rhoCoreAE(iSpecies) *= SQRT_4PI;
   rhoCorePS(iSpecies) *= SQRT_4PI;
   vBar(iSpecies)      *= SQRT_4PI;

}

SxDiracMat<Double> getMatAtomPAW (const SxList<SxString> &potFileTk,
                                  int nWaves, int matSize,
                                  const SxVector<Int> &lPhi,
                                  bool verbose)
{

   SxDiracVec<Double> matEl(matSize);
   SxList<SxString> cLine;

   int i=0;
   while (!potFileTk(i).contains("MATRIX")) i++;
   i++;

   for (int j=0;;j+=3,i++) {
      cLine = potFileTk(i).tokenize(' ');

      if (cLine.getSize () < 3) { 
         
         if (cLine.getSize () == 1) matEl(j) = cLine(0).toDouble ();
         else if (cLine.getSize () == 2) {

            matEl(j) = cLine(0).toDouble ();
            matEl(j+1) = cLine(1).toDouble ();
         }
         
         break;
      }
      
      matEl(j) = cLine(0).toDouble ();
      matEl(j+1) = cLine(1).toDouble ();
      matEl(j+2) = cLine(2).toDouble ();
      
   }

   if (verbose) cout << matEl << endl;
   
   SxDiracMat<Double> res(nWaves,nWaves);
   res.set(0);
   i = 0;
   for (int ipt = 0; ipt < nWaves; ++ipt)  {
      // TODO: l-check!
      for (int jpt = ipt; jpt < nWaves && lPhi(ipt) == lPhi(jpt); ++jpt, ++i)  {
         res(ipt, jpt) = res(jpt, ipt) = matEl(i);
      }
   }

   return res;
}

void SxPAWPot::readCoreAbinit (const SxString &file, int is)
{
   SxFileParser fp(file);
   fp.verbose = true;

   fp.topic ("header");
   fp.read ("All-electron core wave");
   fp.nextLine (2);
   double Z = fp.getDouble ();
   if (fabs(Z-nuclearCharge(is)) > 1e-2)  {
      fp.where ();
      cout << endl << "Inconsistency: Z for " << prettyName(is)
           << " should be " << nuclearCharge(is) 
           << ", but core waves are for Z=" << Z << endl;
      SX_QUIT;
   }
   double nCore = fp.getDouble ();
   if (fabs(nCore+valenceCharge(is)-nuclearCharge(is)) > 1e-2)  {
      fp.where ();
      cout << endl << "Inconsistency: core for " << prettyName(is)
           << " should have " << (nuclearCharge(is) - valenceCharge(is))
           << " electrons, but core waves have " << nCore << endl;
      SX_QUIT;
   }
   fp.nextLine (3);
   int nCoreWave = fp.getInt ();
   fp.nextLine ();

   // --- read l quantum numbers
   lCore(is).resize (nCoreWave);
   for (int iWave = 0; iWave < nCoreWave; ++iWave)
      fp >> lCore(is)(iWave);
   fp.nextLine ();

   // --- read mesh in file
   int nMeshes = fp.getInt ();
   fp.nextLine ();
   int meshId, type, meshSize;
   fp >> meshId >> type >> meshSize;
   if (type != 2)  {
      cout << "Unknown mesh type " << type << " in " << file << endl;
      cout << "Expected mesh type 2 (abinit logarithmic mesh)" << endl;
      SX_EXIT;
   }
   double rStep, logdr;
   fp >> rStep >> logdr;
   fp.nextLine (nMeshes);

   // --- generate file mesh (abinit log grid)
   double rMax = fp.getDouble ();
   SxVector<Double> radFile (meshSize);
   for (int ir = 0; ir < meshSize; ++ir)  {
      radFile(ir) = rStep * (exp (ir * logdr) - 1.);
   }
   rMax = radFile(meshSize - 1);

   // --- read core waves
   psiCoreAE(is).reformat (rad(is).getSize (), nCoreWave);
   psiCoreAE(is).setBasis (*radBasisPtr);
   psiCoreAE(is).handle->auxData.is = is;

   for (int iWave = 0; iWave < nCoreWave; ++iWave)  {
      fp.topic ("core wave " + SxString (iWave+1));
      fp.nextLine ();
      fp.read ("===== Core wave function");
      fp.nextLine ();
      int iMesh = fp.getInt ();
      if (iMesh != 1)  {
         fp.where ();
         cout << ": unexpected mesh " << iMesh << endl;
         SX_QUIT;
      }
      fp.nextLine ();
      int n, l, iSpin;
      fp >> n >> l >> iSpin;
      if (l != lCore(is)(iWave))  {
         fp.where ();
         cout << endl << "Inconsistency: l should be " << lCore(is)(iWave)
              << " according to header, but here reads " << l << endl;
         SX_QUIT;
      }
      fp.nextLine (2);
      SxVector<Double> psi = fp.getVector (meshSize);

      // --- interpolate to actual grid
      SxNaturalCubicSpline spline (radFile, psi);
      for (int ir = 0; ir < rad(is).getSize (); ++ir)  {
         double r = rad(is)(ir);
         psiCoreAE(is)(ir, iWave) = (r < rMax) ? spline.getVal(r) : 0.;
      }
      psiCoreAE(is).colRef (iWave) /= rad(is);
   }
   fp.close ();
}

void SxPAWPot::readAtomPAW (const SxString &filename, int iSpecies)
{

   SxString potFile;
   SxList<SxString> cLine;
   SxXC::XCFunctional xcFunctional = SxXC::Unknown;
   
   if (verbose) cout << "reading file...\n";

   try {
      potFile = SxFileIO::readBinary (filename,-1);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   SxList<SxString> potFileTok = potFile.tokenize ('\n');
   SxList<SxString> valList;
   if (!potFileTok(0).contains("ATOMTYPE"))  {
      cout << "Expected ATOMTYPE in line 1 of " << filename << endl;
      cout << "but found:" << endl
           << potFileTok(0)
           << endl;
      cout << "Is this really an AtomPAW potential file?" << endl;
      SX_QUIT;
   }

   // XC-Functional
   if ( (potFileTok(1).right("ATOMXCTYPE").stripWhiteSpace ()) == "LDA-PW") {
      cout << "LDA-PW" << endl;
      xcFunctional = SxXC::LDA_PW;
   }
   else if ( 
         (potFileTok(1).right("ATOMXCTYPE").stripWhiteSpace ()) == "GGA-PBE") {
       cout << "GGA-PBE" << endl;
       xcFunctional = SxXC::PBE;
   } else if (kjXcShape.getSize () > 0) {
      // silently ignore unknown generating functional if not needed
      cout << "Unrecognized xc functional in " << filename << endl;
      cout << potFileTok(1).right("ATOMXCTYPE").stripWhiteSpace () << endl;
      SX_EXIT;
   }

   // nuclear/valence charge
   nuclearCharge(iSpecies) = (potFileTok(2).right ("ATOMIC_CHARGE")
                             .stripWhiteSpace ()).toDouble ();   
   // atompaw 4.0 has an extra lines here
   while (!potFileTok(3).contains ("CORE_CHARGE"))  {
      potFileTok.remove (3);
      if (potFileTok.getSize () == 3)  {
         cout << "Illegal atomPAW file " << filename << ": no CORE_CHARGE"
              << endl;
         SX_EXIT;
      }
   }
   double coreCharge = (potFileTok(3).right ("CORE_CHARGE")
                       .stripWhiteSpace ()).toDouble ();
   valenceCharge(iSpecies) = nuclearCharge(iSpecies) - coreCharge;

   // rc
   double rPAW = potFileTok(4).right("RC").stripWhiteSpace ().toDouble ();

   // --- fetch radial grid
   double r0;
   int meshSize   = potFileTok(13).right("MESH_SIZE").toInt ();
   int psCoreSize = potFileTok(16).right("CORETAIL_POINTS").toInt ();
   int lcaoSize   = potFileTok(17).right("LCAO_SIZE").toInt ();
   int nPtsMax, psRhoSize = -1;

   if ( lcaoSize >= meshSize ) {
      nPtsMax = lcaoSize;
      for (SxList<SxString>::ConstIterator it = potFileTok.begin ();
           it != potFileTok.end ();
           ++it)
      {
         if ((*it).contains ("PSEUDO_VALENCE_DENSITY"))  {
            SxList<SxString> splitLine = (*it).tokenize (' ');
            if (splitLine.getSize () == 2)  {
               psRhoSize = splitLine(1).toInt ();
               if (psRhoSize > nPtsMax) {
                  nPtsMax = psRhoSize;
                  cout << "Setting max. points to " << nPtsMax 
                       << " from pseudo valence density." << endl;
               }
            } else {
               psRhoSize = nPtsMax;
            }
            break;
         }
      }
      r0 =  potFileTok (15).right ("LOG_GRID").toDouble ();
      logDr(iSpecies) = potFileTok (18).right ("LCAO_STEP").toDouble ();

      SxDiracVec<Double> rMesh(nPtsMax);
      for (int j=0; j<nPtsMax; j++) rMesh(j) = r0*(exp (j*logDr(iSpecies)));   
      rad(iSpecies) = rMesh;

   } else {
          cout << "\n" <<
          "+------------------------------+\n" <<
          "| EXITING: lcaoSize < meshSize |\n" <<
          "+------------------------------+\n" << endl; 
          SX_QUIT; 
   }

   double logdr = logDr(iSpecies);

   // --- fetch shape function parameters
   SxList <SxString> 
      shapeFuParam = potFileTok(5).right ("SHAPE_TYPE").tokenize (' ');

   SxDiracVec<Double> g0File; 
   if (  shapeFuParam(0) == "sinc2" ) {
      cout << "Sinc shape function[(sin(pi*arg/rc_hat)/(pi*arg/rc_hat))**2]\n"
           << "Assuming k(rShape)=1e-4" << endl;
      double rShape;
      //use optional rShape
      if ( shapeFuParam.getSize () > 1 ) 
         rShape = shapeFuParam(1).toDouble ();
      else
         rShape = rPAW;
      rc(iSpecies) = rShape/sqrt (-log(1e-4));
      int nr = (int)rad(iSpecies).getSize ();
      g0File.resize (nr);
      for (int ir = 0; ir < nr; ++ir)  {
         double r = rad(iSpecies)(ir);
         if (r < rShape)  {
            double x = PI * r / rShape;
            g0File(ir) = sqr(sin(x)/x);
         } else {
            g0File(ir) = 0.;
         }
      }
   } else if ( shapeFuParam(0) == "gaussian" ) {
      cout << "Gaussian shape function [exp(-(arg/sigma)^2)]" << endl;
      rc(iSpecies) = shapeFuParam(1).toDouble ();
      g0File = exp(-(rad(iSpecies)/rc(iSpecies)).sqr ());
   } else {
      cout << "Shape function: " << shapeFuParam(0) <<" not recognized" << endl;
      rc(iSpecies) = rPAW/sqrt (-log(1e-4));
      cout << "Shape will be read from file..." << endl;
   }

   // number of orbitals
   int nWaves = potFileTok(6).right ("BASIS_SIZE").stripWhiteSpace().toInt ();

   // --- l-Ordering & lMax
   cLine = potFileTok(8).tokenize (' ');
   SxArray<int> lOrdering(nWaves);

   for (int i=0; i<nWaves; i++) lOrdering(i) = cLine(i).toInt ();
 
   lPhi(iSpecies) = lOrdering;
   cout << "lOrdering: " << lOrdering << endl;

   for (int ipt = 1; ipt < nWaves; ++ipt)  {
      if (lOrdering(ipt) < lOrdering(ipt-1))  {
         cout << "Unexpected partial wave ordering in '" << filename << "'.\n";
         cout << "PAW basis must be sorted by increasing l" << endl;
         SX_QUIT;
      }
   }
 
   lMax(iSpecies) = lOrdering(nWaves-1);

   // --- initial occupation numbers
   foccInit(iSpecies).resize (nWaves);
   cLine = potFileTok(11).tokenize (' ');
   for (int i=0; i<nWaves; i++)
      foccInit(iSpecies)(i) = cLine(i).toDouble ();
   double nValElec = foccInit(iSpecies).sum ();
   if (fabs(nValElec - valenceCharge(iSpecies)) > 1e-3) {
      cout << "WARNING: occupation numbers sum up to "
           << nValElec << ", but valence charge is " 
           << valenceCharge(iSpecies) << endl;
   }
   

   potFileTok = potFile.substitute("END","&").tokenize('&');

   // --- fetch AE core density
   rhoCoreAE(iSpecies) = 
      getFunction (potFileTok(2),verbose,meshSize,nPtsMax,logdr,true);
   
   rhoCoreAE(iSpecies).handle->auxData.is = iSpecies;
   rhoCoreAE(iSpecies).handle->auxData.l = 0;
   rhoCoreAE(iSpecies).handle->auxData.m = 0;
   
   cout << "core charge (rho = sum(|phi|^2/4*PI*r^2)): " 
        << (rad(iSpecies)*rhoCoreAE(iSpecies)).integrate (logdr)
        << endl;

   // --- fetch PS core density
   rhoCorePS(iSpecies) = 
      getFunction (potFileTok(3),verbose, psCoreSize,nPtsMax,logdr,true);
   for (int ir = meshSize; ir < nPtsMax; ++ir)  {
      rhoCoreAE(iSpecies)(ir) = rhoCorePS(iSpecies)(ir);
   }
   for (int ir = meshSize-1; ir; --ir)  {
      if (rad(iSpecies)(ir) < rPAW) break;
      rhoCoreAE(iSpecies)(ir) = rhoCorePS(iSpecies)(ir);
   }
   
   rhoCorePS(iSpecies).handle->auxData.is = iSpecies;
   rhoCorePS(iSpecies).handle->auxData.l = 0;
   rhoCorePS(iSpecies).handle->auxData.m = 0;

   //redefining to our convention
   rhoCoreAE(iSpecies) /= sqrt (FOUR_PI) * rad(iSpecies).sqr ();
   rhoCorePS(iSpecies) /= sqrt (FOUR_PI) * rad(iSpecies).sqr ();

   // PS valence density
   SxDiracVec<Double> rhoValPS = 
      getFunction (potFileTok(4), verbose, psRhoSize, nPtsMax, logdr);
   rhoValPS /= sqrt(FOUR_PI) * rad(iSpecies).sqr ();
   rhoInit(iSpecies) = rhoValPS + rhoCorePS(iSpecies);
   rhoInit(iSpecies).handle->auxData.is = iSpecies;
   rhoInit(iSpecies).handle->auxData.l  = 0;
   rhoInit(iSpecies).handle->auxData.m  = 0;
 
   if (g0File.getSize () == 0)  {
      g0File = getFunction (potFileTok(5), verbose, meshSize, nPtsMax,logdr);
      cout << "Shape norm from file: "
           << (g0File * rad(iSpecies).cub ()).integrate (logdr)
           << endl;
   }
   // normalize compensation charge
   g0File /= FOUR_PI * (g0File * rad(iSpecies).cub ()).integrate (logdr);
   
   // fetch vloc (in Ry)
   vBar(iSpecies) = getFunction (potFileTok(6),verbose,meshSize,nPtsMax,logdr);
   vBar(iSpecies) *= 0.5;
   vBar(iSpecies).handle->auxData.is = iSpecies;
   vBar(iSpecies).handle->auxData.l = 0;
   vBar(iSpecies).handle->auxData.m = 0;
     
   // --- fetch projector functions
   SxList<SxString>::Iterator it;
   pPS(iSpecies) = readWaves (it=potFileTok.begin (),potFileTok.end (),
                              "TPROJECTOR","PHI", verbose, nPtsMax, nWaves,
                              meshSize, logdr);
   
   for (int iWave = 0; iWave < nWaves; ++iWave)
      pPS(iSpecies).colRef (iWave) /= rad(iSpecies);


   // --- fetch AE partial waves   
   phiAE(iSpecies) = readWaves (it,potFileTok.end (),"PHI","TPHI", verbose,
                                nPtsMax,nWaves,meshSize,logdr,true);
      
   for (int iWave = 0; iWave < nWaves; ++iWave)
      phiAE(iSpecies).colRef (iWave) /= rad(iSpecies);
       
   // --- fetch PS partial waves
   phiPS(iSpecies) = 
      readWaves (it,potFileTok.end (),"TPHI","TPHI_LCAO", verbose,
            nPtsMax,nWaves,meshSize,logdr,true);
      // TODO WRONG OVERLAP ?!?!
      //readWaves (it,potFileTok.end (),"TPHI_LCAO","OVERLAP_SIZE", verbose,
      //           nPtsMax,nWaves,lcaoSize,logdr,true);

   SX_LOOP(iWave) {
      phiPS(iSpecies).colRef (iWave) /= rad(iSpecies);
      // --- enforce phiAE=phiPS outside cutoff region
      SX_LOOP(ir) {
         if (rad(iSpecies)(ir) > rPAW)
           phiAE(iSpecies)(ir, iWave) = phiPS(iSpecies)(ir, iWave);
      }

   }

   {
      writePlot ("rhops.dat", rad(iSpecies), rhoValPS);

      // hack SxRadBasis for getHartreePotential
      SxRadBasis radB;
      radB.radFunc.resize (iSpecies+1);
      radB.radFunc(iSpecies) = rad(iSpecies);
      radB.logDr.resize (iSpecies+1);
      radB.logDr(iSpecies) = logdr;

      int nRad = (int)rad(iSpecies).getSize ();
      if (fabs(rhoValPS(nRad - 1)) > 1e-6)  {
         cout << "Valence charge density for " << prettyName(iSpecies)
              << " in " << filename << " does not drop to 0, but is "
              << rhoValPS(nRad - 1)
              << endl;
         SX_EXIT
      }

      // Monopole moment of pseudo-density
      double Q = - sqrt(FOUR_PI) * ((rhoValPS + rhoCorePS(iSpecies))
                                 * rad(iSpecies).cub ()).integrate (logdr);
      Q += nValElec - valenceCharge(iSpecies);
      sxprintf ("compensation Q= %.12f\n", Q);
      rhoInitQ0(iSpecies) = Q;
      SxDiracVec<Double> g0 = getGrl(iSpecies, 0) / FOUR_PI, vScr;
      //SxDiracVec<Double> g0 = g0File; // DEBUG only
      g0File.setBasis (&radB);
      g0File.handle->auxData.is = iSpecies;
      g0File.handle->auxData.l = 0;
      g0File.handle->auxData.m = 0;
      g0.setBasis (&radB);
      g0.handle->auxData.is = iSpecies;
      g0.handle->auxData.l = 0;
      g0.handle->auxData.m = 0;
      cout << "g0 moment: "
           << (FOUR_PI * (g0 * rad(iSpecies).cub ()).integrate (logdr)) << endl;
      cout << "g0File moment: "
           << (FOUR_PI * (g0File * rad(iSpecies).cub ()).integrate (logdr)) << endl;
      //vScr = SxRadialAtom::getHartreePotential (Q * (g0-g0File));
      //vScr = SxRadialAtom::getHartreeSmooth (Q * (g0-g0File));
      vScr = SxRadialAtom::getHartreePotential (Q * g0);
      vScr -= SxRadialAtom::getHartreeSmooth (Q * g0File);
      vBar(iSpecies) -= vScr; 

      if (kjXcShape.getSize () > 0)  {
         Q = valenceCharge(iSpecies) 
           - sqrt(FOUR_PI)*(rhoValPS * rad(iSpecies).cub ()).integrate (logdr);
         kjXcShape(iSpecies) = g0File * FOUR_PI;
         vBar(iSpecies)
            -= SxRadialAtom::computeXC(rad(iSpecies), logdr,
                                       rhoValPS/sqrt(FOUR_PI)
                                       + rhoCorePS(iSpecies)/sqrt(FOUR_PI)
                                       + Q * g0File, xcFunctional)
             - SxRadialAtom::computeXC(rad(iSpecies), logdr,
                                       rhoValPS/sqrt(FOUR_PI)
                                       + rhoCorePS(iSpecies)/sqrt(FOUR_PI),
                                       xcFunctional);
      }
   }

   writePlot ("vbar.dat", rad(iSpecies), vBar(iSpecies));

   vBar(iSpecies) *= sqrt(FOUR_PI);
   // --- check normalization 
   cout << "Normalization <phiAE|phiAE>" << endl;
   for (int iWave =0;iWave < nWaves; iWave++) 
      cout << " i=" << iWave << ": " 
           << (phiAE(iSpecies).colRef(iWave).sqr ()*rad(iSpecies).cub ())
              .integrate(logdr) << endl;

   // --- check orthonormality condition 
   cout << "Orthonormality conditions <phiPS|pPS>" << endl;
   for (int ipt=0;ipt<nWaves;ipt++)  {
      for (int jpt=0; jpt<nWaves; jpt++)  { 
         if (lPhi(iSpecies)(ipt) == lPhi(iSpecies)(jpt))  {
            double phiP = (  phiPS(iSpecies).colRef(ipt)
                           * pPS(iSpecies).colRef(jpt)
                           * rad(iSpecies).cub ()).integrate(logDr(iSpecies));
            cout << " i=" << ipt << "/j=" << jpt << ": " << phiP << endl;
            /*
            // --- renormalize
            if (ipt == jpt)
               pPS(iSpecies).colRef (ipt) /= phiP;
            else
               pPS(iSpecies).colRef (jpt).plus_assign_ax (-phiP,
                     pPS(iSpecies).colRef (ipt));
            */

         }
      }
   }
   
   // --- overlap-matrix deltaS(iSpecies) from recomputeOverlap?
   SxRegex ovlpMatch("\\s+OVERLAP_SIZE\\s+"+intNum);
   while ( ovlpMatch.match (*it).getSize () == 0 ) it++;
   int matSize = ovlpMatch.match ((*it))(1).toInt ();
   SxList<SxString> potFileTk =(*it).tokenize('\n');
   
   deltaS(iSpecies) = getMatAtomPAW (potFileTk, nWaves, matSize,
                                     lPhi(iSpecies), verbose);

   // --- kinetic energy-matrix
   it++;
   potFileTk =(*it).tokenize('\n');

   deltaKin(iSpecies) = getMatAtomPAW (potFileTk, nWaves, matSize, 
                                       lPhi(iSpecies), verbose);
   deltaKin(iSpecies) *= 0.5; // Ry to H

   int oldPrec = (int)cout.precision ();
   cout.precision (12);
   cout << "deltaKin:" << endl;
   cout << deltaKin(iSpecies) << endl;
   cout << "deltaS:" << endl;
   cout << deltaS(iSpecies) << endl;
   cout.precision (oldPrec);
}

void SxPAWPot::createFineBasis (int is, double r0, double rMax, int nPts)
{
   if (!fineRadBasisPtr)  {
      int nSpecies = SxSpeciesData::getNSpecies ();
      fineRadBasisPtr = SxPtr<SxRadBasis>::create ();
      fineRadBasisPtr->radFunc.resize (nSpecies);
      fineRadBasisPtr->logDr.resize (nSpecies);
      fineRadBasisPtr->logDr.set (-1.);
      pPsFine.resize (nSpecies);
   }
   SX_CHECK (rMax > r0 && r0 > 0., rMax, r0);
   SX_CHECK (nPts > 2, nPts);
   double dex = log(rMax/r0) / (nPts - 1);
   fineRadBasisPtr->logDr(is) = dex;
   SxDiracVec<Double> &r = fineRadBasisPtr->radFunc(is);
   r.resize (nPts);
   for (int ir = 0; ir < nPts; ++ir)  {
      r(ir) = r0 * exp (dex * ir);
   }
}

void SxPAWPot::readVasp (const SxString &filename, int is, bool useProjG)
{
   int oldPrec = (int)cout.precision ();
   if (verbose) cout << "Reading file " << filename << endl;
   SxFileParser fp(filename);
   fp.verbose = verbose;
   SxString line;

   fp.topic ("header");
   // title
   fp.nextLine ();
   // valence charge
   double vOld = valenceCharge(is);
   fp >> valenceCharge(is);
   if (vOld > 0. && fabs(vOld - valenceCharge(is)) > 1e-3)  {
      cout << "Valence charge from input file (" << vOld 
           << ") disagrees with POTCAR value (" << valenceCharge(is) << ").\n";
      SX_QUIT;
   }

   // ignore line
   fp.nextLine (2);

   SxXC::XCFunctional xcf = SxXC::Unknown;
   double rShape = -1.;
   for (int i = 0; i < 1000; ++i)  { // break is inside
      line = fp.getLine ();
      if (line.getSize () < 10) continue;
      if (line (10) != '=') break;
      SxString tag = line.left ("=").trim ();
      if (tag == "LEXCH")  {
         SxString xcType = line.subString (12,13);
         if (xcType == "CA")  {
            xcf = SxXC::LDA;
         } else if (xcType == "WI")  {
            xcf = SxXC::LDA_PW; // our best LDA
         } else if (xcType == "PE")  {
            xcf = SxXC::PBE;
         } else {
            cout << "Cannot interprete LEXCH = " << xcType << endl;
            SX_QUIT;
         }
      } else if (tag == "RDEP")  {
         rShape = line.subString (12,19).toDouble ();
      } else if (tag == "POMASS")  {
         double m = line.right ("=").left (";").toDouble ();
         if (m > 0. && reciprocalMass(is) <= 0.)
            reciprocalMass(is) = m;
         if (m > 0. && ionicMass(is) <= 0.)
            ionicMass(is) = m;
      } else if (tag == "RCORE")  {
         rCore(is) = line.subString (12,19).toDouble ();
      } else {
         if (verbose) cout << "ignoring tag=" << tag << endl;
      }
   }
   // ignore line ("Description", header)
   if (!line.contains ("Description")) fp.find (" Description");
   fp.skipWhite ();
   fp.nextLine ();

   lMax(is) = 0;
   int npt = 0;
   {
      int l = -1, type;
      while (fscanf (fp.fp, " %d %*f %d %*f \n", &l, &type) == 2)  {
         npt++;
      }
   }
   lPhi(is).resize (npt, true);
   lPhi(is).set (-1);

   fp.find ("END of PSCTR-control");
   fp.nextLine (2);

   double gMax = fp.getDouble ();
   gMax /= A2B;

   fp.topic ("local potential");
   int nPts = 1000; // hardcoded in VASP
   SxVector<Double> vLoc = fp.getVector (nPts);
   // \tilde v_{Zc} = vLoc - Z/r
   vLoc *= A2B * A2B * A2B / HA2EV;

   double beta = 1., betaLoc = beta;
   for (int ig = 1; ig < nPts; ++ig)  {
      double g = ig * gMax / nPts;
      /*
      cout << g
           << '\t' << vLoc(ig)
           << '\t' << vLoc(ig) * sqr(g)
           << '\t' << (vLoc(ig) - FOUR_PI * valenceCharge(is) / sqr(g)
                       * (1. - exp(-sqr(0.5 * g) / 4.)))
           << endl;
      */
      vLoc(ig) -= FOUR_PI * valenceCharge(is) / sqr(g) 
                  * (1. - exp(-sqr(betaLoc * g) / 4.));
   }

   /*
   {
      double dr = 0.01;
      SxVector<Double> vLocR = toRSpace (vLoc, gMax, 10., dr);
      FILE *out = fopen (("vLoc"+prettyName(is)+".dat").ascii (), "w");
      for (int i = 1; i < vLocR.getSize (); ++i)
         fprintf (out, "%f\t%.16f\n", i*dr,
                  vLocR(i) - valenceCharge(is) * erf(i * dr / beta) / (i * dr));
      fclose (out);
   }
   */

   // consume gradient correction information (if present)
   if (fp.reads ("gradient corrections used for XC")) fp.nextLine (2);

   fp.topic ("core charge density");
   SxVector<Double> psCoreG;
   if (fp.reads ("core charge-density (partial)")) {
      psCoreG = fp.getVector (nPts);
   } else {
      // --- no core ?!
      if (verbose)
         cout << "No pseudo-core found - initialize to zero." << endl;
      psCoreG.resize (nPts);
      psCoreG.set (0.);
   }
   double normCorePS = psCoreG(0);

   /*
   {
      double dr = 0.01;
      SxVector<Double> psCoreR = toRSpace (psCoreG, gMax, 10., dr);

      beta = 0.5;
      // --- core potential
      SxVector<Double> vCoreG(psCoreG.getSize ());
      vCoreG(0) = -PI * sqr(beta) * normCorePS;
      cout << 0. << '\t' << vCoreG(0) << endl;
      for (int ig = 1; ig < nPts; ++ig)  {
         double g = ig * gMax / nPts;
         vCoreG(ig) = FOUR_PI / sqr(g) 
                    * (psCoreG(ig) - normCorePS * exp(-sqr(beta * g)/4.));
         cout << g << '\t' << vCoreG(ig) << endl;
      }
      SxVector<Double> vCoreR = toRSpace (vCoreG, gMax, 10., dr);

      FILE *out = fopen (("psCore-"+prettyName(is)+".dat").ascii (), "w");
      for (int i = 1; i < psCoreR.getSize (); ++i)
         fprintf (out, "%f\t%.16f\t%.12f\n", i*dr,
                  psCoreR(i), vCoreR(i) + normCorePS * erf(i*dr/beta) / (i * dr));
      fclose (out);
   }
   */
   if (fp.reads ("kinetic energy density (partial)"))
      fp.getVector (nPts);
   if (fp.reads ("Kinetic energy density valence"))
      fp.getVector (nPts);

   fp.topic ("atomic pseudo charge-density");
   fp.read ("atomic pseudo charge-density");
   SxVector<Double> psValG = fp.getVector (nPts);

   double normValPS = psValG(0);
   int nG = (int) psValG.getSize ();
   SxPtr<SxRadGBasis> radGBasis 
      = SxPtr<SxRadGBasis>::create (0., gMax * double(nG-1)/double(nG), nG);
   rhoInitBasis(is) = radGBasis;
   // prefactor: sqrt(4pi) / sqrt(2pi)^3
   rhoInit(is) = toVector(psValG + psCoreG) * sqrt(0.5) / PI;
   rhoInit(is).setBasis (*radGBasis);
   rhoInit(is).handle->auxData.is = is;
   rhoInit(is).handle->auxData.l = 0;
   rhoInit(is).handle->auxData.m = 0;

   //cout << normValPS << endl;
   /*
   {
      double dr = 0.01;
      SxVector<Double> psValR = toRSpace (psValG, gMax, 10., dr);

      beta = 0.5;
      // --- core potential
      SxVector<Double> vValG(psValG.getSize ());
      vValG(0) = -PI * sqr(beta) * normValPS;
      cout << 0. << '\t' << vValG(0) << endl;
      for (int ig = 1; ig < nPts; ++ig)  {
         double g = ig * gMax / nPts;
         vValG(ig) = FOUR_PI / sqr(g) 
                    * (psValG(ig) - normValPS * exp(-sqr(beta * g)/4.));
         cout << g << '\t' << vValG(ig) << endl;
      }
      SxVector<Double> vValR = toRSpace (vValG, gMax, 10., dr);

      FILE *out = fopen (("psVal-"+prettyName(is)+".dat").ascii (), "w");
      for (int i = 1; i < psValR.getSize (); ++i)
         fprintf (out, "%f\t%.16f\t%.12f\n", i*dr,
                  psValR(i), vValR(i) + normValPS * erf(i*dr/beta) / (i * dr));
      fclose (out);
   }
   */

   fp.topic ("Non-local part");
   // get max G for reciprocal space projectors
   double gMaxProj = fp.getDouble ();
   gMaxProj /= A2B;
   fp.nextLine ();

   // --- Read projectors and Dij
   deltaKin(is) = SxDiracMat<Double> (npt, npt);
   deltaKin(is).set (0.);

   int nPtsProjG = 100; // hard coded in VASP
   int nPtsProjR = 100; // hard coded in VASP
   SxList<SxVector<Double> > projG, projR;
   SxList<double> rMaxProj;
   double rMinProj = -1.;
   for (int ipt = 0; ipt < npt; /* empty */)  {
      fp.skipWhite ();
      line = fp.getLine ();
      if (line != "Non local Part\n")  {
         if (line == "PAW radial sets\n")  {
            break;
         }
         cout << "Expected Non local Part, but found:" << endl;
         cout << line << endl;
         SX_EXIT;
      }
      int l, nlpro;
      double rMax;
      fp.topic ("projectors (ipt=" + SxString (ipt+1) + ")");
      fp >> l >> nlpro >> rMax;
      rMax *= A2B;
      if (rMinProj < rMax) rMinProj = rMax;
      if (ipt + nlpro > npt)  {
         fp.where ();
         cout << endl;
         cout << "Inconsistency: found more projectors (" << (ipt + nlpro)
              << ") than expected from header (" << npt << ")." << endl;
         cout << "Unknown variant of VASP potential format "
                 "- please contact developers." << endl;
         SX_QUIT;
      }

      // Read Dij block for current l
      for (int i = 0; i < nlpro; ++i)
         for (int j = 0; j < nlpro; ++j)
            fp >> deltaKin(is)(ipt+i,ipt+j);

      for (int i = 0; i < nlpro; ++i, ++ipt)  {
         lPhi(is)(ipt) = l;
         if (lMax(is) < l) lMax(is) = l;
         rMaxProj << rMax;
         fp.read ("Reciprocal Space Part");
         projG << fp.getVector (nPtsProjG);
         fp.read ("Real Space Part");
         projR << fp.getVector (nPtsProjR);

         projG.last () *= A2B * sqrt(A2B);

         // --- smoothen real-space representation
         {
            projG.last ().resize (2 * nPtsProjG, true);
            double p = projG.last ()(nPtsProjG-1),
                   pp = p - projG.last ()(nPtsProjG-2);
            p += pp;
            for (int ig = 0; ig < nPtsProjG; ++ig)  {
               double x = ig / double(nPtsProjG);
               projG.last ()(nPtsProjG + ig) = (1. - x*x*(3. - 2. * x)) * (p + ig * pp) * exp(-sqr(ig/(0.1 * nPtsProjG)));
            }
         }
         if (verbose) {
            double dr = 0.01;
            SxVector<Double> pR
               = toRSpace (projG.last (), 2. * gMaxProj, 10., dr, 0., l);

            SX_MPI_MASTER_ONLY
            {
               writePlot ("p-rec-"+prettyName(is)+"-" + SxString(ipt) + ".dat",
                           dr, toVector(pR));
               writePlot ("p-G-"+prettyName(is)+"-" + SxString(ipt) + ".dat",
                           gMaxProj / nPtsProjG, toVector(projG.last ()));
            }

         }
      }
   }
   deltaKin(is) /= HA2EV;

   if (projG.getSize () != npt)  {
      npt = (int)projG.getSize ();
      // --- resize deltaKin
      SxDiracMat<Double> deltaKinNew(npt, npt);
      for (int ipt = 0; ipt < npt; ++ipt)
         for (int jpt = 0; jpt < npt; ++jpt)
            deltaKinNew(jpt, ipt) = deltaKin(is)(jpt, ipt);
      deltaKin(is) = deltaKinNew;
      // --- determine lMax
      lMax(is) = 0;
      for (int ipt = 0; ipt < npt; ++ipt)
         if (lMax(is) < lPhi(is)(ipt)) 
            lMax(is) = lPhi(is)(ipt);
      // --- resize lPhi
      lPhi(is).resize (npt, true);
   } else {
      fp.read ("PAW radial sets");
   }
   if (verbose) cout << "lPhi = " << lPhi(is) << endl;

   // --- read grid params
   fp.topic ("grid parameters");
   int nRad = fp.getInt ();
   //double psdmax = fp.getDouble ();

   // --- read occupation numbers
   fp.topic ("occupation numbers"); 
   if (verbose) cout << "Scanning for VASP's 'uccopancies'" << endl;
   fp.find ("uccopancies");
   fp.nextLine ();
   SxMatrix<Double> foccMat (npt, npt);
   for (int i = 0; i < foccMat.getSize (); ++i) fp >> foccMat(i);
   DijInit(is) = foccMat;
   SxVector<Double> focc(npt);
   foccInit(is).resize (npt);
   for (int ipt = 0; ipt < npt; ++ipt)  {
      foccInit(is)(ipt) = focc(ipt) 
                        = foccMat (ipt, ipt) * (2*lPhi(is)(ipt) + 1);
   }
   if (verbose) cout << "Sum of focc:" << focc.sum () << endl;

   // --- try to round foccInit
   if (fabs(round(focc.sum ()) - valenceCharge(is)) < 1e-12
       && fabs(focc.sum () - valenceCharge(is)) > 1e-12)  {
      for (int ipt = 0; ipt < npt; ++ipt)  {
         foccInit(is)(ipt) = round(focc(ipt));
      }
      if (fabs(foccInit(is).sum () - valenceCharge(is))
          > fabs(focc.sum () - valenceCharge(is)))  {
         foccInit(is) <<= focc;
      }
   }


   // --- read grid
   fp.topic ("radial grid");
   fp.read ("grid");
   rad(is)   = fp.getVector (nRad);
   rad(is) *= A2B;
   logDr(is) = log(rad(is)(nRad-1) / rad(is)(0)) / (nRad - 1);

   fp.topic ("ae potential");
   fp.read ("aepotential");
   SxDiracVec<Double> vAE = fp.getVector (nRad);
   vAE /= A2B*A2B*A2B * HA2EV;

   fp.topic ("core charge density");
   fp.read ("core charge-density");
   rhoCoreAE(is) = fp.getVector (nRad);
   rhoCoreAE(is) /= A2B;
   rhoCoreAE(is) /= rad(is).sqr ();
   rhoCoreAE(is).handle->auxData.is = is;
   rhoCoreAE(is).handle->auxData.l = 0;
   rhoCoreAE(is).handle->auxData.m = 0;

   fp.topic ("pspotential");
   fp.find ("pspotential");
   fp.reads (" valence only"); // appears in some recent PAW potentials
   SxDiracVec<Double> vPS = fp.getVector (nRad);
   vPS /= A2B*A2B*A2B * HA2EV;

   fp.topic ("pseudo-core");
   fp.read ("core charge-density (pseudized)");
   rhoCorePS(is) = fp.getVector (nRad);
   rhoCorePS(is) /= A2B;
   rhoCorePS(is) /= rad(is).sqr ();
   rhoCorePS(is).handle->auxData.is = is;
   rhoCorePS(is).handle->auxData.l = 0;
   rhoCorePS(is).handle->auxData.m = 0;

   // writePlot ("rhoCore", rad(is), rhoCoreAE(is), rhoCorePS(is));

   if (verbose) writePlot ("potentials", rad(is), vAE, vPS);

   // --- read wave functions
   phiPS(is) = SxDiracMat<Double> (nRad, npt);
   phiAE(is) = SxDiracMat<Double> (nRad, npt);
   fp.topic ("partial waves");
   for (int ipt = 0; ipt < npt; ++ipt)  {
      fp.read ("pseudo wavefunction");
      phiPS(is).colRef(ipt) <<= fp.getVector (nRad);
      phiPS(is).colRef(ipt) /= rad(is);
      fp.read ("ae wavefunction");
      phiAE(is).colRef(ipt) <<= fp.getVector (nRad);
      phiAE(is).colRef(ipt) /= rad(is);
   }
   phiAE(is) /= sqrt(A2B);
   phiPS(is) /= sqrt(A2B);

   fp.read ("End of Dataset");
   fp.close ();

   // --- make radial grid large enough
   int nRadOld = nRad;
   {
      rMinProj *= 1.3; // needed for smooth cutoff of projectors
      SxNaturalCubicSpline coreG(psCoreG, true);
      int N = 20, ng = (int)psCoreG.getSize ();
      double dg = gMax / (ng * N);
      if (normCorePS > 0.)  {
         double rMax = rad(is)(nRad -1), rho;
         for (int i = 0; i < 20; ++i)  {
            rho = toRSpace (rMax, coreG, N, (ng-1)*N, dg);
            if (rho * sqr(rMax) < 1e-10 * normCorePS)  {
               if (rho < 0.) rMax *= 1.3;
               break;
            } else {
               rMax *= 1.1;
            }
         }
         if (rMinProj < rMax) rMinProj = rMax;
      }

      if (rad(is)(nRad-1) < rMinProj)  {
         int nExtra = int(log (rMinProj / rad(is)(nRad-1)) / logDr(is)) + 1;
         rad(is).resize (nRad + nExtra, true);
         for (int i = nRad; i < nRad+nExtra; ++i)
            rad(is)(i) = rad(is)(0) * exp (logDr(is) * i);
         if (verbose)  {
            cout << "Extending radial grid from " << rad(is)(nRad-1)
                 << " to " << rad(is)(nRad+nExtra-1)
                 << " (+" << nExtra << " points)" << endl;
         }

         // --- extend radial functions
         SxDiracMat<Double> newPS (nRad + nExtra, npt),
                            newAE (nRad + nExtra, npt);
         SxIdx old(0, nRad-1);
         for (int ipt = 0; ipt < npt; ++ipt)  {
            newPS.colRef(ipt)(old) <<= phiPS(is).colRef(ipt);
            newAE.colRef(ipt)(old) <<= phiAE(is).colRef(ipt);
            if (newPS(nRad -1, ipt) * newPS(nRad - 2, ipt) > 1e-16 )  {
               double a, b;
               a = log (   rad(is)(nRad-1) * newPS(nRad - 1, ipt) 
                        / (rad(is)(nRad-2) * newPS(nRad - 2, ipt)))
                   / (rad(is)(nRad-1) - rad(is)(nRad-2));
               if (a < 0.)  {
                  // --- exponential decay
                  b = newPS(nRad - 1, ipt) * rad(is)(nRad - 1) / exp (a * rad(is)(nRad-1));
                  for (int i = nRad; i < nRad + nExtra; ++i)
                     newPS(i, ipt) = newAE(i, ipt) = b * exp(a * rad(is)(i)) / rad(is)(i);
               } else {
                  // --- Gaussian cutoff with matched 1st/2nd derivative
                  const SxDiracVec<Double> &r = rad(is);
                  double rOld = r(nRad-1),
                         bCut = 1., c, d;
                  double r12 = r(nRad-1) - r(nRad-2),
                         r23 = r(nRad-2) - r(nRad-3),
                         r34 = r(nRad-3) - r(nRad-4),
                         r13 = r12 + r23,
                         r24 = r23 + r34;

                  // a = f'(-1.5)
                  a = (newPS(nRad - 1, ipt) - newPS(nRad - 2, ipt)) / r12;
                  // c = f'(-2.5)
                  c = (newPS(nRad - 2, ipt) - newPS(nRad - 3, ipt)) / r23;
                  // d = f'(-3.5)
                  d = (newPS(nRad - 3, ipt) - newPS(nRad - 4, ipt)) / r34;
                  // d = 2f"(-3)
                  d = (c - d) / r24;
                  // c = 2f"(-2)
                  c = (a - c) / r13;
                  // c = 2f"(-1)
                  c += (c - d) * r12 / r23;
                  // a = f'(-1)
                  a += 0.25 * c * r12;

                  c+= 1./sqr(bCut) * newPS(nRad-1, ipt);
                  a -= c * rOld;
                  b = newPS(nRad - 1, ipt) - a * r(nRad-1);
                  for (int i = nRad; i < nRad + nExtra; ++i)
                     newPS(i, ipt) = newAE(i, ipt) 
                                   = (b + (a + c * (r(i) - rOld)) * r(i))
                                   * exp(-sqr((r(i)-rOld)/bCut));
               }
            } else {
               for (int i = nRad; i < nRad + nExtra; ++i)
                  newPS(i, ipt) = newAE(i, ipt) = 0.;
            }
         }
         phiPS(is) = newPS;
         phiAE(is) = newAE;

         // --- extend core 
         rhoCorePS(is).resize (nRad + nExtra, true);
         rhoCoreAE(is).resize (nRad + nExtra, true);
         bool expDecay = (normCorePS < 1e-6);
         double a = 0., b = 0.;
         for (int i = nRad; i < nRad + nExtra; ++i)  {
            if (expDecay)  {
               rhoCorePS(is)(i) = rhoCoreAE(is)(i) = a * exp(b * rad(is)(i));
            } else {
               rhoCorePS(is)(i) = rhoCoreAE(is)(i) 
                                = toRSpace (rad(is)(i), coreG, N, (ng-1)*N, dg)
                                * sqrt(FOUR_PI);
               if (rhoCorePS(is)(i) < 1e-6)  {
                  expDecay = true;
                  if (rhoCorePS(is)(i) > 0. && rhoCorePS(is)(i-1) > 0.)  {
                     b = log(rhoCoreAE(is)(i)/rhoCoreAE(is)(i-1))
                       / (rad(is)(i) - rad(is)(i-1));
                     a = rhoCoreAE(is)(i) / exp (rad(is)(i) * b);
                  }
                  if (b >= 0. || a < 0.)  {
                     cout << "Warning: non-exponential decay in core." << endl;
                     cout << "Core is truncated at r=" << rad(is)(i) << endl;
                     a = b = 0.;
                  }
               }
            }
         }
         nRad += nExtra;
      }
   }

   // --- get projectors on radial grid
   pPS(is)     = SxDiracMat<Double> (nRad, npt);
   for (int ipt = 0; ipt < npt; ++ipt)  {
      SxNaturalCubicSpline spline(projR(ipt), true);
      double rMax = rMaxProj(ipt);
      /*
      // --- interpolate real-space projectors
      double xMax = 1. - 1. / nPtsProjR;
      for (int ir = 0; ir < nRad; ++ir)  {
         double x = rad(is)(ir) / rMax;
         if (x >= xMax)
            pPS(is)(ir, ipt) = 0.;
         else
            pPS(is)(ir, ipt) = spline.getVal (nPtsProjR * x)
                             / (A2B * sqrt(A2B));
      }
      */

      // --- Fourier-transform reciprocal-space projectors
      int N = 10;
      double dg = gMaxProj / (nPtsProjG * N);
      SxNaturalCubicSpline pG(projG(ipt), true);
      for (int ir = 0; ir < nRad; ++ir)  {
         double r = rad(is)(ir);
         if (r > 1.3 * rMax)  {
            pPS(is)(ir, ipt) = 0.;
         } else {
            pPS(is)(ir, ipt) = toRSpace (r, pG , N, (2*nPtsProjG-1)*N,
                                         dg, lPhi(is)(ipt));
            if (r > rMax)  {
               double x = (r - rMax) / (0.3 * rMax);
               pPS(is)(ir, ipt) *= (1. - x*x*(3. - 2. *x));
            }
         }
      }

      // --- output (debug)
      if (verbose)  {
         writePlot ("p-" + prettyName(is) + "-" + ipt, rad(is),
                     pPS(is).colRef (ipt), phiPS(is).colRef (ipt),
                     phiAE(is).colRef (ipt));
      }
   }

   // --- orthonormality check
   if (verbose) {
      cout << lPhi(is) << endl;
      for (int ipt = 0; ipt < npt; ++ipt)
         for (int jpt = 0; jpt < npt; ++jpt)
            if (lPhi(is)(ipt) == lPhi(is)(jpt))  {
               cout << "<phiPS" << ipt << "|pPS" << jpt << ">="
                    << (phiPS(is).colRef(ipt) * pPS(is).colRef(jpt)
                        * rad(is).cub ()).integrate (logDr(is))
                    << endl;
            }
   }

   // --- get projectors on radial G grid
   if (useProjG)  {
      cout << "Using VASP G spline projectors" << endl;
      double gLastProj = gMaxProj * double(2*nPtsProjG-1)/double(nPtsProjG);
      projBasis(is) = SxPtr<SxRadGBasis>::create (0., gLastProj, 2*nPtsProjG);
      pPS(is)     = SxDiracMat<Double> (nPtsProjG * 2, npt);
      pPS(is).setBasis (*projBasis(is));
      SX_LOOP(ipt)
         pPS(is).colRef (ipt) <<= toVector(projG(ipt));
      pPS(is) /= TWO_PI * sqrt(TWO_PI);
   }

   double coreCharge = ((rhoCoreAE(is) - rhoCorePS(is)) * rad(is).cub ())
                       .integrate (logDr(is))* sqrt(FOUR_PI)
                     + normCorePS;
   if (verbose)  {
      cout << "norm core num=" << coreCharge << endl;
      cout << "norm PS core exact = " << normCorePS << endl;
   }
   nuclearCharge(is) = valenceCharge(is) + round(coreCharge);
   if (fabs(round(coreCharge) - coreCharge) > 1e-2)  {
      cout << "Failed to determine consistent charges!" << endl;
      cout << "Core charge numerical: " << coreCharge << " => "
           << round(coreCharge) << " (rounded)" << endl;
      cout << "Valence charge as read from file: " << valenceCharge(is) << endl;
      cout << "Deduced nuclear charge: " << nuclearCharge(is) << endl;
      cout << "Discrepancy: " << (round(coreCharge) - coreCharge) << endl;
      SX_EXIT;
   }

   // --- kinetic energy matrix elements 
   SxDiracMat<Double> kinNum(npt, npt);
   kinNum.set (0.);
   //double alpha = 137.036;
   for (int ipt = 0; ipt < npt; ++ipt)  {
      int l = lPhi(is)(ipt);
      for (int jpt = 0; jpt < npt; ++jpt)
         if (lPhi(is)(ipt) == lPhi(is)(jpt))  {
            kinNum(ipt, jpt) = ((phiPS(is).colRef(ipt) 
                     * SxRadialAtom::laplace (phiPS(is).colRef(jpt), rad(is),
                                              logDr(is), l)
                     -phiAE(is).colRef(ipt) 
                     * SxRadialAtom::laplace (phiAE(is).colRef(jpt), rad(is),
                                              logDr(is), l)
                     )
                     * rad(is).cub ()).integrate (logDr(is)) / 2.;
            /*
            cout << "<phi" << ipt << "|T|phi" << jpt << ">="
                 << kinNum(ipt, jpt)
                 << ':' << ((laplace (phiAE(is).colRef(ipt), rad(is), logDr(is), l) * laplace (phiAE(is).colRef(jpt), rad(is), logDr(is), l) * rad(is).cub ()).integrate (logDr(is)) * 3./8. / sqr(alpha))
                 << ':' << ((phiAE(is).colRef(ipt) * laplace (laplace (phiAE(is).colRef(jpt), rad(is), logDr(is), l), rad(is), logDr(is), l) * rad(is).cub ()).integrate (logDr(is)) * 3./8. / sqr(alpha))
                 << endl;
            */
         }
   }

   double logdr = logDr(is);

   // hack SxRadBasis for getHartreePotential
   SxRadBasis radB;
   radB.radFunc = rad;
   radB.logDr = logDr;
   
   // --- get valence densities
   SxDiracVec<Double> rhoValPsBare(radB(is)), rhoValAE(radB(is));
   rhoValPsBare.set (0.);
   rhoValAE.set (0.);
   for (int ipt = 0; ipt < npt; ++ipt)  {
      for (int jpt = 0; jpt < npt; ++jpt)  {
         if (fabs(foccMat(ipt,jpt)) > 1e-12)  {
            if (lPhi(is)(ipt) != lPhi(is)(jpt))  {
               SX_EXIT;
            }
            rhoValPsBare += foccMat(ipt,jpt) * phiPS(is).colRef(ipt)
                                             * phiPS(is).colRef(jpt)
                                             * double(2*lPhi(is)(ipt) + 1);
            rhoValAE     += foccMat(ipt,jpt) * phiAE(is).colRef(ipt)
                                             * phiAE(is).colRef(jpt)
                                             * double(2*lPhi(is)(ipt) + 1);
         }
      }
   }
   SxDiracVec<Double> rhoValPsAug = FOUR_PI * toRSpace (psValG, gMax, rad(is));
   for (int ir = nRadOld; ir < nRad; ++ir)  {
      rhoValPsBare(ir) = rhoValAE(ir) = rhoValPsAug(ir);
   }
   SxDiracVec<Double> g0File2 = (rhoValPsAug-rhoValPsBare);
   g0File2 /= (g0File2 * rad(is).cub ()).integrate (logdr);

   if (verbose)  {
      cout << "rhoValPS:" << (rhoValPsAug * rad(is).cub ()).integrate (logdr)
           << endl;
      cout << "rhoValAE:" << (rhoValAE * rad(is).cub ()).integrate (logdr)
           << endl;
   }

   /*
   // recompute augmentation charge to exactly match AE part
   double deltaQ = ((rhoValAE-rhoValPsAug) * rad(is).cub ()).integrate(logdr);
   rhoValPsAug.plus_assign_ax (deltaQ, g0File2); 
   */

   SX_MPI_MASTER_ONLY
   {
      double deltaCore = coreCharge - (rhoCoreAE(is) * rad(is).cub ())
                         .integrate (logdr) * sqrt(FOUR_PI);
      if (verbose) cout << "delta rhoCoreAE:" << deltaCore << endl;
      if (fabs(deltaCore) > 0.1)  {
         cout << "Norm of core charge has been lost!" << endl;
         SX_EXIT;
      }
   }
   rhoValPsAug /= FOUR_PI;
   rhoValPsBare /= FOUR_PI;
   rhoValAE /= FOUR_PI;

   rhoValAE.handle->auxData.is = is;
   rhoValAE.handle->auxData.l = 0;
   rhoValAE.handle->auxData.m = 0;
   rhoValPsAug.handle->auxData.is = is;
   rhoValPsAug.handle->auxData.l = 0;
   rhoValPsAug.handle->auxData.m = 0;
   rhoValPsAug.setBasis (&radB);
   rhoValPsBare.handle->auxData.is = is;
   rhoValPsBare.handle->auxData.l = 0;
   rhoValPsBare.handle->auxData.m = 0;

   {
      SxDiracVec<Double> nHat = rhoValPsAug - rhoValPsBare;
      nHat *= sqrt(FOUR_PI);
      nHat.handle->auxData = rhoValPsAug.handle->auxData;
      rhoInit(is) -= (*radGBasis) | nHat;
   }
   rhoInit(is) = radGBasis->toSpline (rhoInit(is));

   if (verbose)
      writePlot ("rhoVal-"+prettyName(is), rad(is), rhoValPsAug,
                  rhoValPsBare, rhoValAE);

   // --- Kresse-Joubert ionic potential in radial space
   SxDiracVec<Double> vZcPS = toRSpace (vLoc, gMax, rad(is));
   // add long-range part to vLoc subtracted above
   vZcPS -= valenceCharge(is) * derf(rad(is) / betaLoc) / rad(is);

   {
      /*
      double r0 = rad(is)(0);
      rShape = r0 * exp(logdr * int(log(rShape/r0)/logdr));
      SxDiracVec<Double> g0File = getJsbShape (rShape, 0, rad(is), logdr);
      g0File.handle->auxData.is = is;
      g0File.handle->auxData.l = 0;
      g0File.handle->auxData.m = 0;
      writePlot ("g00", rad(is), g0File);
      cout << (g0File * r3).integrate (logdr) << endl;
      */

      SxDiracVec<Double> rhoijPS, rhoijAE;
      rhoCoreAE(is).setBasis (&radB);
      rhoCorePS(is).setBasis (&radB);
      SxDiracVec<Double> vZcAE = (-nuclearCharge(is))/rad(is) +
             SxRadialAtom::getHartreePotential (rhoCoreAE(is)/ sqrt(FOUR_PI)),
                         r3 = rad(is).cub ();
      
      // 0..r0 correction
      // e= 0.5 * Z * r0^2 * rho(r0) -  1/6  Z * r0^3 * rho'(r0)
      // integration weight @ r0=rad(0): logdr r0^3 * w_Simpson(0)
      double rhoPrimePot = nuclearCharge(is) / 6.
                        / (rad(is)(0) - rad(is)(1));
      rhoPrimePot=0.; // strangely, this seems to be better for abInit
      vZcAE(0) -= (0.5 * nuclearCharge(is)/rad(is)(0) - rhoPrimePot)
                / (logdr /3.);
      vZcAE(1) -= rhoPrimePot / (logdr * 4./3.);

      if (verbose) writePlot ("vZc", rad(is), vZcAE, vZcPS);

      // --- valence Hartree potential
      SxDiracVec<Double> vPsFile = vPS, vAeFile = vAE;
      vAE = SxRadialAtom::getHartreePotential (rhoValAE) ;
      vPS = SxRadialAtom::getHartreePotential (rhoValPsAug) ;

      // correction for missing charge outside rMax
      double deltaV;
      {
         beta = 0.5;
         // --- core potential
         SxVector<Double> vValG(nPts);
         vValG(0) = -PI * sqr(beta) * normValPS;
         for (int ig = 1; ig < nPts; ++ig)  {
            double g = ig * gMax / nPts;
            vValG(ig) = FOUR_PI / sqr(g) 
                       * (psValG(ig) - normValPS * exp(-sqr(beta * g)/4.));
         }
         double r0 = rad(is)(0);
         deltaV = toRSpace (r0, SxNaturalCubicSpline (vValG, true),
                            10, 10 * (nPts - 1), gMax / (10 * nPts));
         deltaV += normValPS * derf(r0/beta) / r0;
         deltaV -= vPS(0);
      }

      //cout << "deltaV = " << deltaV << endl;
      vAE += deltaV;
      vPS += deltaV;

      // writePlot ("vHVal.dat", rad(is), vAE, vPS);
      // --- add xc potential
      vAE += SxRadialAtom::computeXC (rad(is), logdr, 
                                     rhoCoreAE(is)/sqrt(FOUR_PI)+rhoValAE, xcf);
      vPS += SxRadialAtom::computeXC (rad(is), logdr,
                                      rhoCorePS(is)/sqrt(FOUR_PI)+rhoValPsAug,
                                      xcf);

      // --- add core potential
      vAE += vZcAE;
      vPS += vZcPS;

      // --- make sure that potentials agree beyond PAW radius
      for (int ir = nRadOld; ir < nRad; ++ir)  {
         vAE(ir) = vPS(ir);
      }

      if (verbose) {
         writePlot ("vAtom-"+prettyName(is) + ".dat", rad(is), vAE, vPS);
         writePlot ("rhoAtom-"+prettyName(is) + ".dat",
                     rad(is), rhoValAE, rhoValPsAug);
      }

      // --- compute kinetic energy matrix elements
      double Xij = 0., Sij = 0., Qij = 0.;
      SxVector<Int> map(npt);
      for (int ipt = 0, ip = 0; ipt < npt; ++ipt)  {
         int l = lPhi(is)(ipt);
         for (int m = -l; m <= l; ++m, ++ip)
            if (m == 0) map(ipt) = ip;
      }
      deltaS(is)   = SxDiracMat<Double> (npt, npt);
      double vComp = (vPS * g0File2 * rad(is).cub ()).integrate (logdr);
      for (int ipt = 0; ipt < npt; ++ipt)  {
         for (int jpt = 0; jpt < npt; ++jpt)  {
            if (lPhi(is)(ipt) == lPhi(is)(jpt))  { 
               rhoijPS = nij(phiPS(is), ipt, jpt);
               rhoijAE = nij(phiAE(is), ipt, jpt);
               Qij = ((rhoijAE - rhoijPS) * r3).integrate (logdr);
               Xij = ((vAE * rhoijAE - vPS * rhoijPS) * r3).integrate(logdr);
               Sij = vComp * Qij;
               SX_CHECK_NUM (Qij);
               deltaKin(is)(ipt,jpt) -= Xij - Sij;
               deltaS(is)(ipt,jpt)   = Qij;
            } else {
               deltaKin(is)(ipt,jpt) = 0.;
               deltaS(is)(ipt,jpt)   = 0.;
            }
         }
      }
      cout.precision (12);
      //cout << "deltaKin:" << endl;
      //cout << deltaKin(is) << endl;
      //cout << "deltaS:" << endl;
      //cout << deltaS(is) << endl;
   }

   if (verbose)  {
      cout << "delta T:" << endl;
      for (int ipt = 0; ipt < npt; ++ipt)  {
         for (int jpt = 0; jpt < npt; ++jpt)  {
            if (fabs(deltaS(is)(ipt, jpt)) > 1e-10)  {
               cout << ipt << "," << jpt << ": "
                    << (deltaKin(is)(ipt, jpt) - kinNum(ipt, jpt))
                    << " (" 
                    << ((deltaKin(is)(ipt, jpt) - kinNum(ipt, jpt)) / deltaKin(is)(ipt, jpt) * 100.) << "%)\n";
            }
         }
      }
   }

   bool kjxc = kjXcShape.getSize () > 0;
   if (kjxc)  {
      kjXcShape(is).copy (g0File2);
      kjXcShape(is).handle->auxData.is = is;
      kjXcShape(is).handle->auxData.l = 0;
      kjXcShape(is).handle->auxData.m = 0;
   }

   // --- unscreen the local potential

   // set our hard Gaussian shape parameter
   rc(is) = rShape / sqrt(-log(1e-4));
   if (exp(-sqr(rad(is)(nRad-1)/rc(is))) > 1e-10)
      rc(is) = rad(is)(nRad-1) / sqrt(-log(1e-10));

   // Monopole moments of valence and core
   double Qv = FOUR_PI * ((rhoValAE - rhoValPsBare) * rad(is).cub ())
               .integrate (logdr);
   double QZc = - nuclearCharge(is)
                + sqrt(FOUR_PI) * ((rhoCoreAE(is) - rhoCorePS(is))
                                   * rad(is).cub ()).integrate(logdr);
   rhoInitQ0(is) = Qv;

   SxDiracVec<Double> grl = getGrl (is, 0);
   grl.setBasis (&radB);
   grl /= FOUR_PI;
   g0File2 /= FOUR_PI;

   //cout << FOUR_PI * (grl * rad(is).cub ()).integrate (logdr) << endl;
   //cout << FOUR_PI * (g0File2 * rad(is).cub ()).integrate (logdr) << endl;

   vBar(is).copy(vZcPS);
   vBar(is).handle->auxData.is = is;
   vBar(is).handle->auxData.l = 0;
   vBar(is).handle->auxData.m = 0;

   // xc and Hartree potentials
   SxDiracVec<Double> vXcKJ, vXcPS, vHdelta; 
   vXcKJ = SxRadialAtom::computeXC(rad(is), logdr,
                                   rhoValPsAug + rhoCorePS(is)/sqrt(FOUR_PI),
                                   xcf);
   vXcPS = SxRadialAtom::computeXC(rad(is), logdr,
                                   rhoValPsBare + rhoCorePS(is)/sqrt(FOUR_PI),
                                   xcf);
   vHdelta = SxRadialAtom::getHartreePotential
   //vHdelta = SxRadialAtom::getHartreeSmooth
                                               (rhoCorePS(is)/sqrt(FOUR_PI)
                                                + (Qv + QZc) * grl 
                                                - Qv * g0File2);
   // cf. Ref 3, Eq. (35) (Bloechl) vs Eq. 89 (Kresse-Joubert)
   if (kjxc)
      vBar(is) -=  vHdelta ;
   else
      vBar(is) +=  vXcKJ -  vXcPS - vHdelta ;

   // --- smooth cutoff
   int nCut = int(log (1.3 / logdr));
   if (nCut > nRad - nRadOld) nCut = nRad - nRadOld;
   for (int ir = nRadOld; ir < nRad; ++ir)  {
      double x = (ir - nRadOld) / double(nCut);
      if (x < 1.)
         vBar(is)(ir) *= 1. - x * x * (3. - 2. * x);
      else
         vBar(is)(ir) = 0.;
   }

   //writePlot ("vBar.dat", rad(is), vBar(is));

   vBar(is) *= sqrt(FOUR_PI);

   rhoCoreAE(is).handle->auxData.basisPtr = NULL;
   rhoCorePS(is).handle->auxData.basisPtr = NULL;

   cout.precision (oldPrec);
}

void SxPAWPot::readCPPAW (const SxString &filename, int iSpecies)
{
   sxprintf ("read file\n");
   FILE *fp = fopen (filename.ascii(), "r");
   if (!fp)  {
      sxprintf ("ERROR: Can't open file %s\n", filename.ascii());
      SX_EXIT;
   }
   char buffer[10240];

   const int colWidth = 14;  // CP-PAW format: 14 chars each column
   //double r1, dex, psz, aez, rc;
   double r0, dex, psz, aez;
   int nr, nwave;

   int i, idx, c, nCols;
   SxString line;

   // --- read 1st line: header information
   //fscanf (fp, "%lf %lf %d %d %lf %lf %lf",
   //fscanf (fp, "%15lf%10lf%4d%4d%5lf%5lf%15lf",
   //        &r0, &dex, &nr, &nwave, &psz, &aez, &rc(iSpecies));
   fgets (buffer, 10240, fp); line  = buffer;
   r0 = line.subString(0,14).toDouble ();
   dex = line.subString(15,24).toDouble ();
   nr = line.subString(25,28).toInt ();
   nwave = line.subString(29,32).toInt ();
   psz = line.subString(33,37).toDouble ();
   aez = line.subString(38,42).toDouble ();
   rc(iSpecies) = line.subString(43,57).toDouble ();

   if (fabs(valenceCharge(iSpecies) - psz) > 1e-4)  {
      /*
      cout << "WARNING: valence charge from " << filename << " (" << psz
           << ") overrides value from input file (" 
           << valenceCharge(iSpecies) << ")." << endl;
      valenceCharge(iSpecies) = psz;
      // SX_QUIT??
      */
      cout << "WARNING: valence charge from input file ("
           << valenceCharge(iSpecies) << ") overrides value from "
           << filename << " (" << psz << ")." << endl;
      psz = valenceCharge(iSpecies);
   }
   nuclearCharge(iSpecies) = aez;
   logDr(iSpecies) = dex;
  
   // --- read 2nd line: lPhi array
   //     defines the angular momentum l of the i-th projector
   //     or i-th partial wave
   SxList<int> lPhiList;
   lMax(iSpecies) = 0;
   for (i=0; i < nwave; i++)   {
      fscanf (fp, "%d", &c);
      lPhiList << c;
      if (c > lMax(iSpecies)) lMax(iSpecies) = c;
   }
   fgets (buffer, 10240, fp);  // skip rest of line
   lPhi(iSpecies) = lPhiList;

   // --- read vBar
   SxDiracVec<Double> &vBarRef = vBar(iSpecies);
   vBarRef.resize (nr);
   for (i=0; i < nr; /*empty*/)  {
      fgets (buffer, 10240, fp); line  = buffer;
      nCols = (int)line.getSize() / colWidth;
      for (c=0, idx=0; c < nCols && i < nr; c++, idx+=colWidth)
         vBarRef(i++) = line.subString(idx, idx+colWidth-1).toDouble();
   }
   vBarRef.handle->auxData.is = iSpecies;
   vBarRef.handle->auxData.l = 0;
   vBarRef.handle->auxData.m = 0;

   // --- read aecore
   rhoCoreAE(iSpecies).resize (nr);
   for (i=0; i < nr; /*empty*/)  {
      fgets (buffer, 10240, fp); line  = buffer;
      nCols = (int)line.getSize() / colWidth;
      for (c=0, idx=0; c < nCols && i < nr; c++, idx+=colWidth)
         rhoCoreAE(iSpecies)(i++) = line.subString(idx, idx+colWidth-1)
                                    .toDouble();
   }
   rhoCoreAE(iSpecies).handle->auxData.is = iSpecies;
   rhoCoreAE(iSpecies).handle->auxData.l = 0;
   rhoCoreAE(iSpecies).handle->auxData.m = 0;

   // --- read pscore
   rhoCorePS(iSpecies).resize (nr);
   for (i=0; i < nr; /*empty*/)  {
      fgets (buffer, 10240, fp); line  = buffer;
      nCols = (int)line.getSize() / colWidth;
      for (c=0, idx=0; c < nCols && i < nr; c++, idx+=colWidth)
         rhoCorePS(iSpecies)(i++) = line.subString(idx, idx+colWidth-1)
                                    .toDouble();
   }
   rhoCorePS(iSpecies).handle->auxData.is = iSpecies;
   rhoCorePS(iSpecies).handle->auxData.l = 0;
   rhoCorePS(iSpecies).handle->auxData.m = 0;

   // --- read dtkin = <aephi|-delta/2|aephi> - <psphi|-delta/2|psphi> 
   SxDiracMat<Double> dtKin(nwave,nwave);
   for (i=0; i < dtKin.getSize(); /*empty*/)  {
      fgets (buffer, 10240, fp); line  = buffer;
      nCols = (int)line.getSize() / colWidth;
      for (c=0, idx=0; c < nCols && i < dtKin.getSize(); c++, idx+=colWidth)
         dtKin(i++) = line.subString(idx, idx+colWidth-1).toDouble();
   }
   deltaKin(iSpecies) = dtKin;

   // --- read dover = <aephi|aephi> - <psphi|psphi>
   SxDiracMat<Double> dover(nwave,nwave);
   for (i=0; i < dover.getSize(); /*empty*/)  {
      fgets (buffer, 10240, fp); line  = buffer;
      nCols = (int)line.getSize() / colWidth;
      for (c=0, idx=0; c < nCols && i < dover.getSize(); c++, idx+=colWidth)
         dover(i++) = line.subString(idx, idx+colWidth-1).toDouble();
   }
   deltaS(iSpecies) = dover;
   rhoInitQ0(iSpecies) = 0.;
   SX_LOOP(ipt) rhoInitQ0(iSpecies) += foccInit(iSpecies)(ipt)
                                     * deltaS(iSpecies)(ipt, ipt);

   // --- read in triple blocks (projector,aephi,psphi) with l=lPhi
   SxDiracMat<Double>  proj (nr,nwave);   // :i,:ir
   SxDiracMat<Double>  phiAE_(nr,nwave);   // :i,:ir
   SxDiracMat<Double>  phiPS_(nr,nwave);   // :i,:ir
   int r;
   for (i=0; i < nwave; i++)  {
      // --- projector(lPhi(i))
      for (r=0; r < nr; /* empty */ )  {
         fgets (buffer, 10240, fp); line  = buffer;
         nCols = (int)line.getSize() / colWidth;
         for (c=0, idx=0; c < nCols && r < nr; c++, idx+=colWidth)
            proj(r++,i) = line.subString(idx, idx+colWidth-1).toDouble();
      }
      // --- aephi(lPhi(i))
      for (r=0; r < nr; /* empty */ )  {
         fgets (buffer, 10240, fp); line  = buffer;
         nCols = (int)line.getSize() / colWidth;
         for (c=0, idx=0; c < nCols && r < nr; c++, idx+=colWidth)
            phiAE_(r++,i) = line.subString(idx, idx+colWidth-1).toDouble();
      }
      // --- psphi(lPhi(i))
      for (r=0; r < nr; /* empty */ )  {
         fgets (buffer, 10240, fp); line  = buffer;
         nCols = (int)line.getSize() / colWidth;
         for (c=0, idx=0; c < nCols && r < nr; c++, idx+=colWidth)
            phiPS_(r++,i) = line.subString(idx, idx+colWidth-1).toDouble();
      }
   }
   phiPS(iSpecies) = phiPS_;
   phiAE(iSpecies) = phiAE_;
   pPS(iSpecies) = proj;

   // --- generate radial mesh
   SxDiracVec<Double> rMesh(nr); 
   rMesh(0) = r0;
   for (i=1; i < nr; i++)
      rMesh(i) = r0 * exp(dex*i);
   rad(iSpecies) = rMesh;

}

int SxPAWPot::getNProj (int iSpecies) const
{
   int n = 0;
   const SxArray<int> &l = lPhi(iSpecies);
   for (int ip = 0; ip < l.getSize (); ++ip)
      n += 2 * l(ip) + 1;
   return n;
}


void SxPAWPot::print () const
{
   for (int is = 0; is < getNSpecies (); ++is)  {
      cout << SX_SEPARATOR;
      cout << "| " << chemName(is) << " (" << elementName(is) << ")" << endl;
      cout << "| Z=" << nuclearCharge(is) << "; nv=" << valenceCharge(is)
           << endl;
      cout << "| Number of partial waves: " << lPhi(is).getSize () << endl;
      cout << "| Number of projectors   : " << getNProj(is) << endl;
      cout << "| max. l for partials    : " << lMax(is)<< endl;
      cout << "| max. l for density     : " << lMaxRho(is) << endl;
      cout << "| radial grid size       : " << rad(is).getSize () << endl;
      cout << "| angular grid size      : "
           << SxSphereGrid(aGridType(is)).getSize () << endl;
      cout << "| compensation charge rc : " << rc(is) << " bohr" << endl;

      cout << "| AE core deficit        : ";
      double nc = nuclearCharge(is) - valenceCharge(is);
      cout << (nc - sqrt(FOUR_PI) * (rhoCoreAE(is) * rad(is).cub ())
              .integrate (logDr(is))) << endl;
      cout << "| AE core energy         : ";
      sxprintf ("%.12f Hartree\n", coreEnergy(is));
   }
   cout << SX_SEPARATOR;
}

SxRadialMesh SxPAWPot::computeRhoPS (const SxMatrix<Double> &Dij, 
                                     int is, int nSpin) const
{
   SxRadialMesh res = computeRho(Dij, is, phiPS(is));
   if (nSpin > 0) res(0,0) += rhoCorePS(is) / double(nSpin);
   return res;
}

SxRadialMesh SxPAWPot::computeRhoAE (const SxMatrix<Double> &Dij,
                                     int is, int nSpin) const
{
   SxRadialMesh res = computeRho(Dij, is, phiAE(is));
   if (nSpin > 0) res(0,0) += rhoCoreAE(is) / double(nSpin);
   return res;
}

SxRadialMesh SxPAWPot::computeRho (const SxMatrix<Double> &Dij, int is,
                                   const SxDiracMat<Double> &phi) const
{
   SX_CLOCK (Timer::ComputeRhoRadial);
   SX_CHECK (phi.nCols () == lPhi(is).getSize (),
             phi.nCols (), lPhi(is).getSize ());

   int nr = (int)phi.nRows ();
   SxRadialMesh rho(nr, lMaxRho(is), 0.);

   // --- loop over partial wave types
   for (int ipt = 0; ipt < getNProjType(is); ++ipt)  {
      int l1 = lPhi(is)(ipt), off1 = offset(is)(ipt) + l1;
      for (int jpt = 0; jpt < getNProjType(is); ++jpt)  {
         int l2 = lPhi(is)(jpt), off2 = offset(is)(jpt) + l2;

         // compute phi_i*phi_j radial density
         // Ref. 1, Eq (16) (PS) / (17) (AE)
         SxDiracVec<Double> phiIphiJ = nij(phi, ipt, jpt);

         int maxL = min(l2 + l1, lMaxRho(is));

         // now the angular momentum sums
         for (int m1 = -l1; m1 <= l1; ++m1)  {
            int lm1 = SxYlm::combineLm(l1,m1);
            for (int m2 = -l2; m2 <= l2; ++m2)  {
               int lm2 = SxYlm::combineLm(l2,m2);
               
               // get Dij value for this (n1m1,n2m2) pair
               double dij = Dij(off1 + m1, off2 + m2);

               // loop over radial mesh's lm's
               // possible l from triangular condition & parity
               for (int l = abs(l2-l1); l <= maxL; l+=2)  {
                  for (int m = -l; m <= l; ++m)  {
                     double cg = clebschGordan(lm1, lm2, SxYlm::combineLm(l,m));
                     // no i^l factors for PAW radial ylm
                     // but is l1+l2-l (from cg definition)
                     if ((l1 + l2 - l) & 2) cg = -cg;
                     // note: many cg's are zero
                     if (fabs(cg) > 1e-12)
                        rho(l,m).plus_assign_ax (cg * dij, phiIphiJ);
                  }
               }
            }
         }
      }
   }
   if (radBasisPtr) rho.setBasis (*radBasisPtr);
   rho.setIs (is);
   return rho;
}

SxDiracVec<TPrecCoeffG>
SxPAWPot::getAtomRhoG(const SxGBasis &G, int iSpecies) const
{
   SxDiracVec<TPrecCoeffG> g0 = exp((-0.25 * sqr(rc(iSpecies))) * G.g2)
                              * rhoInitQ0(iSpecies);
   if (rhoInit(iSpecies).getSize () > 0)  {
      return (G | rhoInit(iSpecies)) + g0;
   }
   SxDiracVec<Double> rhoRad (getRadBasis ()((int)iSpecies));
   rhoRad <<= rhoCorePS (iSpecies);
   double Y00 = SQRT_1_4PI;
   SxDiracVec<Double> r3 = rad(iSpecies).cub ();
   for (int ipt = 0; ipt <lPhi(iSpecies).getSize (); ++ipt)  {
      double focc    = foccInit(iSpecies)(ipt);
      SxDiracVec<Double> psRho = phiPS(iSpecies).colRef(ipt).sqr ();

      // --- enforce normalization via pseudo-rho
      double norm = (r3 * psRho).integrate (logDr(iSpecies))
                  + deltaS(iSpecies)(ipt,ipt);
      /*
      if  (fabs(norm - 1.) > 1e-6 && fabs(focc) > 1e-12)  {
         cout << "WARNING: " << elementName(iSpecies)
               << " channel " << (ipt+1) << "(l=" << lPhi(iSpecies)(ipt)
               << ") had to be renormalized (norm=" << norm << ")." << endl;
      }
      */
      // always renormalize
      focc /= norm;

      // radial density
      rhoRad.plus_assign_ax(Y00 * focc, psRho);
   }

   rhoRad.handle->auxData.is = iSpecies;
   rhoRad.handle->auxData.l  = 0;
   rhoRad.handle->auxData.m  = 0;
   return ( G | rhoRad) + g0;
}

SxDiracVec<Double> SxPAWPot::getGrl (int is, int l) const
{
   // normalized generalized Gaussian
   // Ref. 1, Eq. 23, Ref. 2 A1
   double prefac = 0.5 * SQRT_PI * pow (rc(is), double(3 + 2 * l));
   for (double dl = l + 0.5; dl > 0.; dl -= 1.) prefac *= dl;
   // Ref. 1, Eq. 23
   SxDiracVec<Double> res = pow (rad(is), double(l))
                          * exp(-(rad(is)/rc(is)).sqr ()) 
                          / prefac;
#ifndef NDEBUG
   double moment = (pow(rad(is), double(l+3)) * res).integrate (logDr(is));
   SX_CHECK (fabs(moment-1.) < 2e-6, moment - 1.);
#endif
   res.handle->auxData.is = is;
   res.handle->auxData.l = l;
   return res;
}

SxMatrix<Double> SxPAWPot::getVMatEl (const SxRadialMesh &V,
                                      const SxDiracMat<Double> &phi) const
{
   SX_CLOCK (Timer::RadVMatEl);
   SX_CHECK(V.meshData.handle);
   int is = V.meshData.handle->auxData.is;
   /// Only phiPS/phiAE is allowed
   SX_CHECK (&phi == &phiPS(is) || &phi == &phiAE(is));
   int np = getNProj(is);
   SxMatrix<Double> res(np, np);
   res.set (0.);
   SxDiracVec<Double> r3 = rad(is).cub ();

   // --- loop over partial wave types
   for (int ipt = 0; ipt < getNProjType(is); ++ipt)  {
      int l1 = lPhi(is)(ipt), off1 = offset(is)(ipt) + l1;
      for (int jpt = 0; jpt < getNProjType(is); ++jpt)  {
         int l2 = lPhi(is)(jpt), off2 = offset(is)(jpt) + l2;

         // compute phi_i*phi_j *r^3
         SxDiracVec<Double> phiIphiJr3 = phi.colRef(ipt) * phi.colRef(jpt) * r3;

         int maxL = min(l2 + l1, lMaxRho(is));

         // now the angular momentum sums
         for (int m1 = -l1; m1 <= l1; ++m1)  {
            int lm1 = SxYlm::combineLm(l1,m1);
            for (int m2 = -l2; m2 <= l2; ++m2)  {
               int lm2 = SxYlm::combineLm(l2,m2);
               
               // loop over radial mesh's lm's
               // possible l from triangular condition & parity
               for (int l = abs(l2-l1); l <= maxL; l+=2)  {
                  for (int m = -l; m <= l; ++m)  {
                     double cg = clebschGordan(lm1, lm2, SxYlm::combineLm(l,m));
                     // no i^l factors for PAW radial ylm
                     // but is l1+l2-l (from cg definition)
                     if ((l1 + l2 - l) & 2) cg = -cg;
                     // note: many cg's are zero
                     if (fabs(cg) > 1e-12)
                        res(off1 + m1, off2 + m2) += 
                           cg * (phiIphiJr3 * V(l,m)).integrate (logDr(is));
                  }
               }
            }
         }
      }
   }
   return res;
}

void SxPAWPot::setBasis (const SxConstPtr<SxRadBasis> &radPtr)
{
   SX_CHECK(radPtr);
   radBasisPtr = radPtr;
   bool kjxc = kjXcShape.getSize () > 0;
   for (int is = 0; is < getNSpecies (); ++is)  {
      if (rhoInit(is).getSize () > 0 && !rhoInitBasis(is))
         rhoInit(is).setBasis (&*radPtr);
      rhoCorePS(is).setBasis (&*radPtr);
      rhoCoreAE(is).setBasis (&*radPtr);
      vBar(is)     .setBasis (&*radPtr);
      phiPS(is)    .setBasis (&*radPtr);
      phiAE(is)    .setBasis (&*radPtr);
      if (!dynamic_cast<const SxRadGBasis*>(pPS(is).getBasisPtr ()))
         pPS(is).setBasis (&*radPtr);
      if (kjxc && kjXcShape(is).getSize () > 0)
         kjXcShape(is).setBasis (radPtr.getPtr ());
   }
}

void SxPAWPot::recomputeOverlap ()
{
   double ovlp;
   for (int is = 0; is < getNSpecies (); ++is)  {
      deltaS(is) = SxDiracMat<Double> (getNProjType(is), getNProjType(is));
      for (int ipt = 0; ipt < getNProjType(is); ++ipt)  {
         for (int jpt = 0; jpt <= ipt; ++jpt)  {
            if (lPhi(is)(ipt) != lPhi(is)(jpt))
               ovlp = 0.;
            else
               // Ref. 1, Sec. VI E, last paragraph (iii)
               ovlp = tr(nij(phiAE(is), ipt, jpt) - nij(phiPS(is), ipt, jpt));
            deltaS(is)(ipt,jpt) = deltaS(is)(jpt, ipt) = ovlp;
         }
      }
      // Renormalize all-electron charge density
      if (fabs(nuclearCharge(is)-valenceCharge(is)) > 1e-6)  {
         double nCoreNow, nCoreExact;
         nCoreNow = sqrt(FOUR_PI) 
                  * (rhoCoreAE(is) * rad(is).cub ()).integrate (logDr(is));
         nCoreExact = nuclearCharge(is) - valenceCharge(is);
         // ensure we make minor change (otherwise: nCore discrepancy?!)
         if (fabs(nCoreNow - nCoreExact) > 1e-2)  {
            cout << "integrated nCore is " << nCoreNow << endl;
            cout << "nCore should be     " << nCoreExact << endl; 
            SX_EXIT;
         }
         cout << "multiplicative renormalize core: error is " 
              << (nCoreNow - nCoreExact) << endl;
         rhoCoreAE(is) *= nCoreExact / nCoreNow;
         rhoCorePS(is) *= nCoreExact / nCoreNow;
         /*
         cout << "renormalizing nuclear charge to correct core norm: "
              << (nCoreNow - nCoreExact) << endl;
              << endl;
         nuclearCharge(is) += nCoreNow - nCoreExact;
         */
      }
   }
}

void SxPAWPot::checkOverlap (int iSpecies, bool crashOnError)
{
   if (dynamic_cast<const SxRadGBasis *>(pPS(iSpecies).getBasisPtr ()))
      return;
   // ---
   // check that eigenvalues of PAW overlap operator are positive
   // 1 + \sum_ij |p_i>S_ij<p_j|
   // => non-trivial eigenfunctions must be linear combination of |p_k>
   // i.e.
   // \sum_j S_ij <p_j|p_k> c_k = (lambda - 1) c_i
   // => all eigenvalues of S_ij <p_j|p_k> must be > -1
   // Setup SP
   int nProjType = getNProjType(iSpecies);
   SxDiracMat<Double> SP(nProjType,nProjType);
   SxDiracVec<Double> rCube = rad(iSpecies).cub ();
   for (int ipt = 0; ipt < nProjType; ipt++)   {
      for (int jpt = ipt; jpt < nProjType; jpt++)   {
         if (lPhi(iSpecies)(ipt) != lPhi(iSpecies)(jpt))
            SP(ipt,jpt) = SP(jpt,ipt) = 0.;
         else
            SP(ipt,jpt) = SP(jpt,ipt) 
               = (pPS(iSpecies).colRef(ipt) * pPS(iSpecies).colRef(jpt) 
                  * rCube).integrate (logDr(iSpecies));
      }
   }
   // --- eigenvalue problem of PAW correction
   SxDiracMat<Complex16> SSP = deltaS(iSpecies) ^ SP;
   SxDiracMat<Complex16>::Eigensystem eig = SSP.eigensystem ();
   // add the 1 to get the overlap eigenvalues
   eig.vals += 1.;
   cout << SX_SEPARATOR;
   cout << "| PAWPotential Overlap Check. Eigenvalues are:" << endl;
   cout << SX_SEPARATOR;
   eig.vals.real ().print(true);
   bool problem = false;
   for (int iVal = 0; iVal < eig.vals.getSize(); iVal++)   {
      if (eig.vals(iVal).re < 1e-10)   {
         cout << "| WARNING: overlap operator is NOT positive definite.\n"
              << "| At least one eigenvalue is zero or smaller." << endl
              << "| PAW potential for " << prettyName(iSpecies)
              << " might be bad!" << endl;
         cout << "| This is almost always catastrophic:" << endl;
         cout << "| - the code may crash" << endl;
         cout << "| - the results may be wrong" << endl;
         if (crashOnError) {
            cout << "| use 'checkOverlap=false;' to override at own risk!" 
                 << endl;
            SX_QUIT;
         }
         problem = true;
      }
   }
   if (problem) {
      double dg = 0.01;
      int ng = int(10./dg);
      SxDiracMat<Double> SPfilt(ng,nProjType);
      SX_LOOP2(ipt,ig)  {
         double g = double(ig) * dg;
         SxDiracVec<Double> bessel = SxRadBasis::jsb(lPhi(iSpecies)(ipt),
                                                     rad(iSpecies) * g);
         SPfilt(ig, ipt) = (bessel * rCube * pPS(iSpecies).colRef(ipt))
                           .integrate (logDr(iSpecies));
      }
      for (double G2 = 95.; G2 > 9.; G2 -= 5.)  {
         SP.set (0.);
         SX_LOOP2(ipt,jpt)  {
            if (lPhi(iSpecies)(ipt) != lPhi(iSpecies)(jpt)) continue;
            SX_LOOP(ig)  {
               double g = double(ig) * dg, g2 = g*g;
               if (g2 > G2) break;
               SP(ipt,jpt) += g2 * dg * SPfilt(ig,ipt) * SPfilt(ig,jpt);
            }
         }
         SP /= 0.5 * PI; 
         SSP = deltaS(iSpecies) ^ SP;
         eig = SSP.eigensystem ();
         eig.vals += 1.;
         cout << "| Eigenvalues at eCut = " << G2 << " Ry: ";
         SX_LOOP (iEigen)  {
            sxprintf("% .3f (", eig.vals(iEigen).re);
            SX_LOOP (ipt) if (eig.vecs(ipt, iEigen).absSqr () > 0.01)
            {
               int l = lPhi(iSpecies)(ipt);
               if (l < 6)
                  cout << SxString("spdfgh")(l);
               else
                  cout << "l=" << l;
               break;
            }
            cout << ") ";
         }
         cout << endl;
      }
   }
   cout << SX_SEPARATOR;
}


void SxPAWPot::extendRad(double newRMax)
{
   int nSpecies = (int)rad.getSize();

   for(int is = 0; is < nSpecies; is++)
   {
      int newSize = int((log(newRMax) - log(rad(is)(0))) / logDr(is)) + 1;
      int oldSize = (int)rad(is).getSize();
      if (newSize <= oldSize)
      {
         cout << "Error in SxPAWPot::extendRad: new Vector length smaller than old Vector length" << endl;
         SX_QUIT;
      }

      // Build new rad, rhoCorePS, rhoCoreAE, and vBar
      SxDiracVec<Double> oldRad = rad(is).getCopy (),
         oldRhoCorePS = rhoCorePS(is).getCopy (),
         oldRhoCoreAE = rhoCoreAE(is).getCopy (),
         oldVBar = vBar(is).getCopy ();
      rad(is) = SxDiracVec<Double>(newSize);
      rhoCorePS(is).resize(newSize);
      rhoCoreAE(is).resize(newSize);
      vBar(is).resize(newSize);
      for(int i = 0; i < newSize; i++)
      {
         if (i < oldSize) {
            rad(is)(i) = oldRad(i);
            rhoCorePS(is)(i) = oldRhoCorePS(i);
            rhoCoreAE(is)(i) = oldRhoCoreAE(i);
            vBar(is)(i) = oldVBar(i);
         } else {
            rad(is)(i) = rad(is)(0)*exp(logDr(is) * double(i));
            rhoCorePS(is)(i) = 0.0;
            rhoCoreAE(is)(i) = 0.0;
            vBar(is)(i) = 0.0;
         }
      }
      int nOrbTypes = (int)phiPS(is).row(0).getSize();
      // Build new phiPS, phiAE, pPS
      SxDiracMat<Double> oldPhiPS = 1.0 * phiPS(is),
         oldPhiAE = 1.0 * phiAE(is),
         oldpPS = 1.0 * pPS(is);
      pPS(is).reformat(newSize, nOrbTypes);
      phiPS(is).reformat(newSize, nOrbTypes);
      phiAE(is).reformat(newSize, nOrbTypes);
      for (int iot = 0; iot < nOrbTypes; iot++)
      {
         double Phin = oldPhiPS.colRef(iot)(oldSize - 1);
         double Phim = oldPhiPS.colRef(iot)(oldSize - 2);
         double Rn = rad(is)(oldSize - 1);
         double Rm = rad(is)(oldSize - 2);
         double lambdaPhi = 10.0;
         if (Phin > 0 && Phim > 0 )
            lambdaPhi = (log(Phin) - log(Phim)) / (Rn - Rm);
         //Output for control
         cout << "iSpecies: " << is << endl;
         cout << "iOrbitalType: " << iot << endl;
         cout << "(Rn - Rm): " << (Rn - Rm) << endl;
         cout << "lambdaPhi: " << lambdaPhi << endl;
         // if not exponential decreasing set lambda positive to fill vectors with zeros
         if ( (lambdaPhi > 0) || (Phin > Phim) ) {
            cout << "WARNING: no exponential decay for PHI (" << is << "|" << iot << ")" << endl;
            lambdaPhi = 10.0;
         }
         
         for(int i = 0; i < oldSize; i++)  {
            phiPS(is).colRef(iot)(i) = oldPhiPS.colRef(iot)(i);
            phiAE(is).colRef(iot)(i) = oldPhiAE.colRef(iot)(i);
         }
         for(int i = oldSize; i < newSize; i++)  {
            if (lambdaPhi < 0)  {
               phiPS(is).colRef(iot)(i) = phiAE(is).colRef(iot)(i) 
                  = phiPS(is).colRef(iot)(oldSize-1) 
                  * exp(lambdaPhi * (rad(is)(i) - rad(is)(oldSize - 1)));
            } 
            else phiPS(is).colRef(iot)(i) = phiAE(is).colRef(iot)(i) = 0.0;
         }
         for(int i = 0; i < oldSize; i++)  
            pPS(is).colRef(iot)(i) = oldpPS.colRef(iot)(i);
         for(int i = oldSize; i < newSize; i++)  
            pPS(is).colRef(iot)(i) = 0.0;
         SxString file = "RadialExtend";
         file += is;
         file += iot;
         file += ".org";
         SxBinIO out; 
         out.open(file, SxBinIO::ASCII_WRITE_ONLY);
         out.writeXYPlot(toVector(rad(is)),toVector(phiPS(is).colRef(iot)));
         out.close();
      }
      }
   // Set new Basis
   setBasis(SxPtr<SxRadBasis>::create(rad, logDr));
}

void SxPAWPot::setupQijL ()
{
   QijL.resize (getNSpecies ());
   omegaPAW.resize (getNSpecies ());
   for (int is = 0; is < getNSpecies (); ++is)  {
      int npt = getNProjType(is);
      QijL(is).reformat (npt, npt);

      // --- loop over partial wave types
      SX_LOOP2(ipt,jpt) {
         QijL(is)(ipt, jpt).resize (lMaxRho(is)+1);

         // compute phi_i*phi_j radial density
         SxDiracVec<Double> phiIphiJ 
            = nij(phiAE(is), ipt, jpt) - nij(phiPS(is), ipt, jpt);

         for (int l = 0; l <= lMaxRho(is); l++)  {
            SxDiracVec<Double> rl = pow(rad(is),double(l+3));
            QijL(is)(ipt, jpt)(l) = (phiIphiJ * rl).integrate (logDr(is));
         }
      }

      // --- calculate AE norm within PAW radius
      // calculation of radius of the PAW sphere
      double rPAW = rc (is) * sqrt(-log(1e-4));
      omegaPAW(is).reformat (npt, npt);
      SxDiracVec<Double> rW (rad(is).getSize ());
      SX_LOOP(ir)  {
         double r = rad(is)(ir);
         rW(ir) = (r < rPAW) ? (r * r * r) : 0.;
      }
      SX_LOOP2(ipt, jpt)  {
         if (lPhi(is)(ipt) == lPhi(is)(jpt))  {
            omegaPAW(is)(ipt, jpt) = (  phiAE(is).colRef(ipt)
                                      * phiAE(is).colRef(jpt)
                                      * rW).integrate (logDr(is));
         } else {
            omegaPAW(is)(ipt, jpt) = 0.;
         }
      }


   }
}

void SxPAWPot::computeXKernel ()
{
   int maxM = 7; // very precise for Mg, O

   // variables for testing moment expansion 
   double dg = 1e-2;
   int nG = int(13. * SxXCFunctional::omegaHSE / dg) + 10;
   SxArray<SxArray2<SxDiracVec<Double> > > nijL;

   xKernel.resize (getNSpecies ());
   MijL.resize (getNSpecies ());
   if (verbose) nijL.resize (getNSpecies ());

   for (int is = 0; is < getNSpecies (); ++is)  {
      int npt = getNProjType (is);
      int n2 = (npt * (npt + 1))/2;
      int n4 = (n2 * (n2 + 1))/2;
      xKernel(is).reformat (lMaxRho(is) + 1, n4);
      
      if (verbose) nijL(is).reformat (n2, lMaxRho(is) + 1);
      MijL(is).reformat (maxM, lMaxRho(is) + 1, n2);
      MijL(is).set (0.);

      SxArray<SxDiracVec<Double> > gL (lMaxRho(is) + 1);
      for (int l = 0; l <= lMaxRho(is); ++l)
         gL(l) = getGrl (is, l);
      
      for (int ipt = 0; ipt < npt; ++ipt)  {
         int l1 = lPhi(is)(ipt);
         for (int jpt = ipt; jpt < npt; ++jpt)  {
            int l2 = lPhi(is)(jpt);
            int ij = get2Idx (ipt, jpt, npt);
            SxDiracVec<Double> nijAE = nij(phiAE(is), ipt, jpt);
            SxDiracVec<Double> nijPS = nij(phiPS(is), ipt, jpt);
            nijAE.handle->auxData.is = is;
            nijAE.handle->auxData.m = 0;
            nijPS.handle->auxData.is = is;
            nijPS.handle->auxData.m = 0;
            for (int kpt = ipt; kpt < npt; ++kpt)  {
               int l3 = lPhi(is)(kpt);
               for (int lpt = (ipt == kpt) ? jpt : kpt; lpt < npt; ++lpt)  {
                  int l4 = lPhi(is)(lpt);
                  int lmin = max(abs(l1-l2), abs(l3-l4));
                  int lmax = min(l1+l2, l3+l4);
                  int kl = get2Idx (kpt, lpt, npt);
                  SX_CHECK (ij <= kl, ij, kl);
                  int ijkl = get2Idx (ij, kl, n2);
                  SxDiracVec<Double> nklAE = nij(phiAE(is), kpt, lpt);
                  SxDiracVec<Double> nklPS = nij(phiPS(is), kpt, lpt);
                  nklAE.handle->auxData.is = is;
                  nklAE.handle->auxData.m = 0;
                  nklPS.handle->auxData.is = is;
                  nklPS.handle->auxData.m = 0;
                  SxDiracVec<Double> nijSoft, nklSoft;
                  for (int l = 0; l <= lMaxRho(is); ++l)  {
                     if ((l1+l2 + l) & 1 || (l3+l4 + l) & 1
                         || (l < lmin) || (l > lmax))
                     {
                        // (ij) or (kl) do not have a l-component
                        xKernel(is)(l, ijkl) = 0.;
                     } else {
                        nijAE.handle->auxData.l = l;
                        nklAE.handle->auxData.l = l;
                        nijSoft = nijPS + QijL(is)(ipt, jpt)(l) * gL(l);
                        nklSoft = nklPS + QijL(is)(kpt, lpt)(l) * gL(l);
                        nijSoft.handle->auxData.l = l;
                        nklSoft.handle->auxData.l = l;

                        // Ref 4, Eq. (48)
                        xKernel(is)(l, ijkl) = 
        SxRadialAtom::computeScreenedHartree (nijAE, nklAE, 0., dg)
      - SxRadialAtom::computeScreenedHartree (nijSoft, nklSoft, 0., dg);

                        if (ij == kl)  {
                           SxDiracVec<Double> deltaNijR3 = nijAE - nijSoft;
                           deltaNijR3 *= rad(is).cub ();

                           // --- calculate Taylor expansion coefficients of 
                           //     g-shape
                           SxArray2<double> F
                             = SxYlm::getJsbDerivCoeffs (l, l + 2*(maxM+1));
                           for (int m = 0; m < maxM; ++m)  {
                              // Ref. 4, Eq. (66)
                              MijL(is)(m,l,ij) = 
        F(0,l+2*m+2) * (deltaNijR3 * pow(rad(is),l+2*m+2)).integrate(logDr(is));
                           }

                           // --- calculate g-shape directly
                           if (verbose)  {
                              nijL(is)(ij, l).resize (nG);
                              for (int ig = 0; ig < nG; ++ig)  {
                                 SxDiracVec<Double> jsb 
                                    = SxRadBasis::jsb (l, (ig * dg) * rad(is));
                                 nijL(is)(ij, l)(ig) = (deltaNijR3 * jsb)
                                                       .integrate(logDr(is));
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
      if (verbose)  {
         FILE *fp = fopen (("nijL-" + prettyName(is) + ".dat").ascii (), "w");
         if (!fp) { cout << "Cannot open output file!" << endl; SX_EXIT; }
         SX_LOOP2(ij, l)  {
            cout << nijL(is)(ij,l).getSize () << endl;
            if (nijL(is)(ij,l).getSize () == 0) continue;
            fprintf (fp, "# ij = %ld, l = %ld\n", ssize_t(ij), ssize_t(l));
            SX_LOOP(ig)  {
               double g = double(ig) * dg;
               double taylor = 0.;
               for (int m = 0; m < maxM; ++m)  {
                  taylor += MijL(is)(m,l,ij) * pow(g, double(l+2*m+2));
               }
               fprintf (fp, "%.6g %.12f %.12f %.12f\n", g , nijL(is)(ij,l)(ig),
                        taylor,  nijL(is)(ij,l)(ig) - taylor);
            }
            fprintf (fp, "&\n");
         }
         fclose (fp);
      }
   }
}

void SxPAWPot::computeCoreX (int is)
{
   SX_CHECK (is >= 0 && is < getNSpecies (), is, getNSpecies ());
   SX_CHECK (psiCoreAE(is).getSize () > 0, is);
   int npt = getNProjType (is);
   int nCore = (int)lCore(is).getSize ();
   int lCoreMax = 0;
   SX_LOOP(iCore) lCoreMax = max(lCoreMax, lCore(is)(iCore));
   SxDiracVec<Double> r3 = rad(is).cub ();
   coreX(is) = SxDiracMat<Double> (npt, npt);
   coreX(is).set (0.);
   // --- bare exchange only
   for (int jpt = 0; jpt < npt; ++jpt)  {
      int lv = lPhi(is)(jpt);
      for (int iCore = 0; iCore < nCore; ++iCore)  {
         int lc = lCore(is)(iCore);
         SxDiracVec<Double> psiC = psiCoreAE(is).colRef (iCore);
         SxDiracVec<Double> rhoPhiCore = phiAE(is).colRef (jpt) * psiC;
         for (int l = abs(lc - lv); l <= lc + lv; l+=2)  {
            // note: sum_{mc,M} <li mi | lc mc L M><lc mc L M | lj mj>
            //       = delta(li,lj) delta(mi,mj)
            //       * wigner3j (li,lc,L,0,0,0)^2
            //       * (2 lc + 1)(2 L + 1) / (4pi)
            double C = sqr(SxYlm::wigner3j(lv, lc, l, 0, 0, 0))
                     * (2 * lc + 1) * (2 * l + 1) / FOUR_PI;
            
            // --- compute V[l, rhoPhiCore]
            rhoPhiCore.handle->auxData.l = l;
            SxDiracVec<Double> psiVl;
            psiVl = SxRadialAtom::getHartreePotential(rhoPhiCore);

            // multiply with psiCore
            psiVl *= psiC;

            // multiply with r^3
            psiVl *= r3;

            // --- calculate matrix elements
            for (int ipt = 0; ipt < npt; ++ipt)  {
               if (lPhi(is)(ipt) == lv)  {
                  coreX(is)(ipt, jpt) -= C * (phiAE(is).colRef (ipt) * psiVl)
                                             .integrate (logDr(is));
               } 
            }
         }
      }
   }
   if (SxXCFunctional::omegaHSE > 1e-10)  {
      SX_CLOCK (Timer::CoreXHSE);
      // --- screening part: reciprocal space integration
      double gMax = 13. * SxXCFunctional::omegaHSE;
      double dg = 1e-4; // untested
      int ng = int(gMax / dg);
      SX_CHECK (ng > 10 && ng < 1000000, ng);
      SxArray2<SxDiracVec<Double> > rhoPhiCore(npt,nCore);
      for (int iCore = 0; iCore < nCore; ++iCore)  {
         SxDiracVec<Double> psiC = psiCoreAE(is).colRef (iCore);
         SX_LOOP(ipt) rhoPhiCore(ipt,iCore) = phiAE(is).colRef (ipt) * psiC;
      }
      SxArray<double> rhoPhiCoreG(npt);
      for (int l = 0; l <= lCoreMax + lMax(is); ++l)  {
         for (int ig = 0; ig < ng; ++ig)  {
            double g = ig * dg;
            SX_START_TIMER (Timer::CoreXHSEjsb);
            SxDiracVec<Double> jsbR3 = SxRadBasis::jsb(l, g * rad(is));
            jsbR3 *= r3;
            SX_STOP_TIMER (Timer::CoreXHSEjsb);
            for (int iCore = 0; iCore < nCore; ++iCore)  {
               int lc = lCore(is)(iCore);
               SX_LOOP(ipt)
                  rhoPhiCoreG(ipt) = (rhoPhiCore(ipt, iCore) * jsbR3)
                                     .integrate (logDr(is));
               SX_LOOP2(ipt,jpt)  {
                  if (lPhi(is)(ipt) != lPhi(is)(jpt)) continue;
                  int lv = lPhi(is)(ipt);
                  if (l < abs(lv - lc) || l > lv+lc) continue;
                  if (((l + lv + lc) & 1)) continue; // odd sum of l values
                  // for C, see above
                  double C = sqr(SxYlm::wigner3j(lv, lc, l, 0, 0, 0))
                           * (2 * lc + 1) * (2 * l + 1) / FOUR_PI;

                  // HSE implementation notes, Eq. 10
                  double dLR = 8. // combined prefactors Eq. 9+10
                      * rhoPhiCoreG(ipt) * rhoPhiCoreG(jpt)
                      * exp(-g*g/(4. * sqr(SxXCFunctional::omegaHSE)));
                  dLR *= weightSimpson(ig, ng) * dg;
                  coreX(is)(ipt,jpt) += C * dLR;
               }
            }
         }
      }
   }
   cout << "coreX = ["; 
   for (int ipt = 0; ipt < npt; ++ipt)  {
      if (ipt > 0) cout << ",\n         ";
      cout << '[';
      for (int jpt = 0; jpt < npt; ++jpt)  {
         if (jpt > 0) cout << ", ";
         cout << coreX(is)(ipt, jpt);
      }
      cout << "]";
   }
   cout << "];" << endl;
   double eX = 0.;
   for (int ipt = 0; ipt < foccInit(is).getSize (); ++ipt)  {
      for (int jpt = 0; jpt < foccInit(is).getSize (); ++jpt)  {
         eX += foccInit(is)(ipt) * foccInit(is)(jpt)
               * coreX(is)(ipt, jpt);
      }
   }
   cout << "core-valence exact exchange: " << eX << endl;
}

SxArray<SxArray<SxDiracVec<Double> > > SxPAWPot::getPhiPS ()
{
   SX_CHECK (radBasisPtr);
   int nSpecies = getNSpecies();
   SxArray<SxArray<SxDiracVec<Double> > > result(nSpecies);
   for(int is = 0; is < nSpecies; is++)   {
      int nProj = (int)phiPS(is).row(0).getSize();
      result(is).resize(nProj);
      for(int ip = 0; ip < nProj; ip++)   {
         result(is)(ip) = phiPS(is).colRef(ip).getCopy ();
         result(is)(ip).setBasis(&*radBasisPtr);
         result(is)(ip).handle->auxData.is = is;
         result(is)(ip).handle->auxData.ia = -1;
         result(is)(ip).handle->auxData.n = ip;
         result(is)(ip).handle->auxData.l = lPhi(is)(ip);
         result(is)(ip).handle->auxData.m = NONE_M;
      }
   }

   return result;
}

SxDiracVec<Double> SxPAWPot::getPhiPS (int is, int ip)
{

   SxDiracVec<Double> result  = phiPS(is).colRef(ip).getCopy ();
   result.setBasis(&*radBasisPtr);
   result.handle->auxData.is = is;
   result.handle->auxData.ia = -1;
   result.handle->auxData.n = ip;
   result.handle->auxData.l = lPhi(is)(ip);
   result.handle->auxData.m = NONE_M;
   return result;
}

