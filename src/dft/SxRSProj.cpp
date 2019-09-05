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

#include <SxRSProj.h>
#include <SxProjector.h>

enum RsTimer {
   RsInit, RsProject, RsProjectSum,
   RsGrad, RsGradSum, RsYlm, RsMatSumProj, RsCollect,
   RsMatSumGrad, RsDistribute, RsPhase, RsSplines
};

SX_REGISTER_TIMERS(RsTimer)
{
   regTimer (RsInit, "rs init");
   regTimer (RsProject, "rs project");
   regTimer (RsProjectSum, "rs project sum");
   regTimer (RsGrad, "rs grad");
   regTimer (RsGradSum, "rs grad sum");
   regTimer (RsYlm, "rs Ylm");
   regTimer (RsMatSumProj, "rs proj matrix");
   regTimer (RsCollect, "rs collect");
   regTimer (RsMatSumGrad, "rs grad matrix");
   regTimer (RsDistribute, "rs distribute");
   regTimer (RsPhase, "rs phase");
   regTimer (RsSplines, "rs spline collect");
}

SxRSProj::SxRSProj (const SxAtomicStructure &str,
                    const SxArray<SxArray<SxDiracVec<Double> > > pIn,
                    const SxRadGBasis &radG,
                    const SxPtr<SxRBasis> &rBasisIn,
                    double betaIn,
                    double drIn, double rMax, double dPhi,
                    int nrnIn)
  : npg(0), rBasisPtr(rBasisIn), dr(drIn), nrn(nrnIn), beta(betaIn), nPsiMax(8)
{
   setup (str, pIn, radG, rMax, dPhi);
}

SxRSProj::SxRSProj (const SxSymbolTable *table,
                    const SxAtomicStructure &str,
                    const SxArray<SxArray<SxDiracVec<Double> > > pIn)
  : npg(0), dr(1e-3), nrn(20), nPsiMax(8)
{
   double dPhi;
   SxMesh3D mesh;
   double eCut = SxGBasis::getECut(table->topLevel ());
   double rMax, dg;
   // --- read symbol table
   try  {
      // get group
      if (table->getName () != "rsProj")
         table = table->getGroup ("rsProj");

      // --- mesh
      if (table->contains ("mesh"))  {
         mesh = SxVector3<Int>(table->get("mesh")->toIntList ());
         if (!SxGBasis::isCommensurableMesh(mesh, str.cell, true))  {
            cout << "ERROR: The provided nlEES mesh is not commensurable!\n";
            SxGBasis::getCommensurableMesh (str.cell, mesh);
            SX_QUIT;
         }
      } else {
         double meshAccuracy = 1.;
         if (table->contains ("meshAccuracy"))
            meshAccuracy = table->get("meshAccuracy")->toReal ();
         mesh = SxGBasis::getCommensurableMesh (eCut, str.cell, meshAccuracy);
      }

      // Gaussian beta
      beta = table->get ("beta")->toReal ();
      // real-space cutoff
      rMax = table->contains ("rMax")
           ? table->get("rMax")->toReal ()
           : 10.;
      // desired accuracy (for auto-cutoff)
      dPhi = table->get ("dPhi")->toReal ();
      dg = table->get ("dg")->toReal ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   // create real-space basis
   rBasisPtr = SxPtr<SxRBasis>::create (mesh, str.cell);

   // --- create radial G basis
   double cutG = max (sqrt(eCut)*1.0001, sqrt(-log(1e-16)/(0.25*beta*beta)));
   int nG = int(cutG / dg);
   cout << "nG=" << nG << endl;
   SxRadGBasis radG(0., cutG*1.0001, nG, SxRadGBasis::Linear);

   setup (str, pIn, radG, rMax, dPhi, eCut);

}

void SxRSProj::setup (const SxAtomicStructure &str,
                      const SxArray<SxArray<SxDiracVec<Double> > > pIn,
                      const SxRadGBasis &radG,
                      double rMax, double dPhi, double eCut)
{
   SX_CLOCK (RsInit);
   SX_CHECK (rBasisPtr);
   SX_CHECK (dr > 0., dr);
   SX_CHECK (beta > 0., beta);
   SX_CHECK (nrn > 0, nrn);
   int nSpecies = str.getNSpecies ();
   SX_CHECK (nSpecies == pIn.getSize (), nSpecies, pIn.getSize ());
   projId.resize (nSpecies);
   rProj.resize (nSpecies);
   rCut.resize (nSpecies);
   lmax.resize (nSpecies);
   lmax.set (0);

   int nr = int(ceil(rMax / dr));

   SxDiracVec<Double> g = radG.getRadGFunc (),
                      gauss = exp((-beta*beta*0.25) * g.sqr ());
      
   // --- output
   SxMesh3D mesh = rBasisPtr->getMesh ();
   double dOmega = str.cell.volume / double(mesh.getSize ());
   SxCell recCell = str.cell.getReciprocalCell ();
   SxCell fullCell(recCell.basis(0) * mesh(0),
                   recCell.basis(1) * mesh(1),
                   recCell.basis(2) * mesh(2));
   fullCell *= 0.5;
   fullCell.setup ();
   Coord gMax = fullCell.getHeights ();
   cout << SX_SEPARATOR;
   cout << "| Real space projectors" << endl;
   cout << "| beta         :  " << beta << endl;
   cout << "| mesh         :  " << mesh << endl;
   cout << "| dPhi         :  " << dPhi << endl;
   if (eCut > 0.)  {
      cout << "| damp 1(a1)   :  "
           << exp(-sqr(beta) * gMax(0) * ( gMax(0) - sqrt(eCut))) << endl;
      cout << "| damp 1(a2)   :  "
           << exp(-sqr(beta) * gMax(1) * ( gMax(1) - sqrt(eCut))) << endl;
      cout << "| damp 1(a3)   :  "
           << exp(-sqr(beta) * gMax(2) * ( gMax(2) - sqrt(eCut))) << endl;
      cout << "| FFT blowup   :  " 
           << (1e-14 * exp(sqr(0.5 * beta) * eCut)) << endl;
   }
   cout << "| damp 3       :  " << gauss(gauss.getSize () - 1) << endl;

   double gMax3 = min(gMax(0), min(gMax(1), gMax(2)));
   gMax3 = 2 * gMax3 - sqrt(eCut);

   // --- compute soft projectors in spline representation
   npg = 0;
   for (int is = 0; is < nSpecies; ++is)  {
      int npt = (int)pIn(is).getSize ();
      rProj(is).resize (npt);
      rCut(is) = 0.;
      int npl = 0; // number of projectors per atom
      for (int ipt = 0; ipt < npt; ++ipt)  {
         int l = pIn(is)(ipt).handle->auxData.l;
         npl += (2 * l + 1);
         if (lmax(is) < l) lmax(is) = l;
         SxDiracVec<Double> pG = radG | pIn(is)(ipt);

         SxDiracVec<Double> p2Trunc = pG.sqr ();
         for (int ig = 0; ig < pG.getSize (); ++ig)  {
            if (g(ig) < gMax3) p2Trunc(ig) = 0.;
            else               p2Trunc(ig) *= sqr(gauss(ig));
         }
         cout << "| wrap error   <  "
              << sqrt(radG.integrate (p2Trunc)) << endl;

         pG *= gauss;

         // --- set up r-space at r-mesh points
         SxVector<Double> pRad(nr+nrn);
         for (int ir = 0; ir < nr+nrn; ++ir)  {
            double r = (ir - nrn) * dr;
            pRad(ir) = radG.integrate (pG * SxRadBasis::jsb (l, fabs(r) * g));
            if (ir < nrn && (l & 1))
               pRad(ir) = -pRad(ir);
         }
         pRad *= sqrt(2. / PI);

         // setup spline from data points
         rProj(is)(ipt) = SxNaturalCubicSpline (pRad);

         // --- need larger real-space cutoff?
         for (int ir = (int)pRad.getSize () - 1; ir; --ir)  {
            double r = (ir - nrn) * dr;
            if (r < rCut(is)) break;
            if (fabs(pRad(ir) * r * r) > dPhi)  {
               rCut(is) = r + dr;
               break;
            }
         }

      }
      cout << "| rCut  (is=" << (is+1) << ") :  " << rCut(is) << endl;
      cout << "| points(is=" << (is+1) << ") :  " 
           << int(4./3.*PI*rCut(is)*rCut(is)*rCut(is)/dOmega) << endl;

      // --- set up projector map
      projId(is).resize (npl);
      for (int ipt = 0, ipl = 0; ipt < npt; ++ipt)  {
         int l = pIn(is)(ipt).handle->auxData.l;
         for (int m = -l; m <= l; ++m, ++ipl)
            projId(is)(ipl) = ProjQuantumNumbers (ipt, l, m);
      }

      // total number of projectors
      npg += npl * str.getNAtoms (is);
   }
   cout << SX_SEPARATOR;
}

SxDiracVec<Complex16> 
SxRSProj::project (const SxAtomicStructure &str,
                   const SxDiracVec<Complex16> &uk,
                   const Coord &kVec) const
{
   SX_CLOCK (RsProjectSum);
   SX_CHECK (str.getNSpecies () == projId.getSize (),
             str.getNSpecies (), projId.getSize ());

   const SxRBasis &rBasis = uk.getBasis<SxRBasis> ();
   SX_CHECK ((rBasis.cell - str.cell).absSqr ().sum () < 1e-10);
   SxMesh3D mesh = rBasis.getMesh ();
   
   SxCell gridCell (str.cell.basis(0) / mesh(0),
                    str.cell.basis(1) / mesh(1),
                    str.cell.basis(2) / mesh(2));

   SxDiracVec<Complex16> res(npg);

   int offset = 0;
   for (int is = 0; is < projId.getSize (); ++is)  {
      int npl = (int)projId(is).getSize ();
      SxVector<Double> ylm(sqr(lmax(is)+1));
      SxArray<SxComplex16> p2(npl);
      ProjQuantumNumbers *pId = projId(is).elements;

      for (int ia = 0; ia < str.getNAtoms(is); ++ia, offset += npl)  {
         Coord r0Rel = gridCell.carToRel(str.getAtom(is,ia));

         // circumscribing box
         Coord f = gridCell.getBoundingBox (rCut(is));
         f += 0.001; // safety range for numerical noise at grid cell boundaries
         
         // range of grid cells (mesh vector)
         SxVector3<Int> from, to, xyz;
         from = ceil(r0Rel-f);
         to   = floor(r0Rel+f);

         // --- set up 1D phase factors
         SX_START_TIMER (RsPhase);
         SxVector<Complex16> phase1D((to-from).sum () +3);
         int off0 = 0,
             off1 = to(0) - from(0) + 1,
             off2 = off1 + to(1) - from(1) + 1;
         off0 -= from(0);
         off1 -= from(1);
         off2 -= from(2);

         xyz.set (0);
         for (xyz(0) = from(0); xyz(0) <= to(0); ++xyz(0))  {
            Coord rVec = gridCell.relToCar (xyz);
            phase1D(off0+xyz(0)) = exp (I * (rVec ^ kVec));
         }
         xyz(0) = 0;
         for (xyz(1) = from(1); xyz(1) <= to(1); ++xyz(1))  {
            Coord rVec = gridCell.relToCar (xyz);
            phase1D(off1+xyz(1)) = exp (I * (rVec ^ kVec));
         }
         xyz(1) = 0;
         for (xyz(2) = from(2); xyz(2) <= to(2); ++xyz(2))  {
            Coord rVec = gridCell.relToCar (xyz);
            phase1D(off2+xyz(2)) = exp (I * (rVec ^ kVec));
         }
         SX_STOP_TIMER (RsPhase);

         p2.set (0.);
         for (xyz(0) = from(0); xyz(0) <= to(0); ++xyz(0))  {
            for (xyz(1) = from(1); xyz(1) <= to(1); ++xyz(1))  {
               for (xyz(2) = from(2); xyz(2) <= to(2); ++xyz(2))  {
                  Coord rVec = gridCell.relToCar (xyz);
                  Coord rel = rVec - str.constRef (is, ia);
                  double rAbs = rel.norm ();
                  if (rAbs > rCut(is)) continue;

                  if (rAbs > 1e-10)  {
                     //SX_CLOCK (RsYlm);
                     SxYlm::getYlmArray (lmax(is), rel, &ylm);
                  } else  {
                     ylm.set (0.);
                     ylm(0) = 1.;
                  }

                  SxComplex16 psiK = uk(mesh.getMeshIdx(xyz, SxMesh3D::Unknown))
                                   * phase1D(off0 + xyz(0))
                                   * phase1D(off1 + xyz(1))
                                   * phase1D(off2 + xyz(2));

                  double ir = rAbs / dr + nrn;
                  double rpr;
                  /*
                  for (int ipl = 0; ipl < npl; ++ipl)  {
                     const ProjQuantumNumbers &id = projId(is)(ipl);
                     if (id.m == -id.l) rpr = rProj(is)(id.n).getValYExtra (ir);
                     p2(ipl) += ylm(SxYlm::combineLm(id.l,id.m)) * rpr * psiK;
                  }
                  */
                  for (int ipl = 0; ipl < npl; )  {
                     //const ProjQuantumNumbers &id = projId(is)(ipl);
                     const ProjQuantumNumbers &id = pId[ipl];
                     rpr = rProj(is)(id.n).getValYExtra (ir);
                     double *ylmp = ylm.elements + SxYlm::combineLm(id.l,0);
                     for (int m = -id.l; m <= id.l; ++m, ++ipl) {
                        p2(ipl) += ylmp[m] * rpr * psiK;
                     }
                  }
               }
            }
         }
         for (int ipl = 0; ipl < npl; ++ipl)  {
            const ProjQuantumNumbers &id = projId(is)(ipl);
            SxComplex16 p = SxYlm::getYlmNormFactor (id.l,id.m) 
                          * gridCell.volume
                          * p2(ipl);
            // missing (-i)^l factor
            if (id.l & 1) p *= -I;
            if (id.l & 2) p = -p;
            res(offset + ipl) = p;
         }
      }
   }
   SX_CHECK (offset == npg, offset, npg);
   return res;
}

#define SX_RS_MAX_N_R 64;
SxDiracVec<Complex16> 
SxRSProj::project (const SxAtomicStructure &str,
                   const SxArray<SxDiracVec<Complex16> > &uk,
                   const Coord &kVec) const
{
   SX_CLOCK (RsProjectSum);
   SX_CHECK (str.getNSpecies () == projId.getSize (),
             str.getNSpecies (), projId.getSize ());
   int nStates = (int)uk.getSize ();

   const int maxnrs = SX_RS_MAX_N_R;
   const SxRBasis &rBasis = uk(0).getBasis<SxRBasis> ();
   SX_CHECK ((rBasis.cell - str.cell).absSqr ().sum () < 1e-10);
   SxMesh3D mesh = rBasis.getMesh ();
   
   SxCell gridCell (str.cell.basis(0) / mesh(0),
                    str.cell.basis(1) / mesh(1),
                    str.cell.basis(2) / mesh(2));

   SxDiracMat<Complex16> res(npg, nStates);
   SxMatrix<Complex16> psiK(maxnrs, nStates);

   int offset = 0;
   for (int is = 0; is < projId.getSize (); ++is)  {
      int npl = (int)projId(is).getSize ();
      SxVector<Double> ylm(sqr(lmax(is)+1)), work(lmax(is)+1);
      SxMatrix<Complex16> p2(npl, nStates);
      SxMatrix<Complex16> phi(npl, maxnrs);
      ProjQuantumNumbers *pId = projId(is).elements;

      int npt =  (int)rProj(is).getSize ();
      int nRSpline = (int)rProj(is)(0).getCoeff ()(0).getSize ();
      SxVector<Double> sCoeff(npt * 2 * nRSpline);
      for (int ipt = 0; ipt < npt; ++ipt)  {
         const SxArray<SxVector<Double> > &polyCoeff = rProj(is)(ipt).getCoeff ();
         const SxVector<Double> &a0 = polyCoeff(0), 
                                &a2 = polyCoeff(1);
         for (int ir = 0; ir < nRSpline; ++ir)  {
            sCoeff(2 * (ir * npt + ipt)    ) = a0(ir);
            sCoeff(2 * (ir * npt + ipt) + 1) = a2(ir);
         }
      }
      VALIDATE_VECTOR(sCoeff);
      SxArray<double> rprAll(npt);

      for (int ia = 0; ia < str.getNAtoms(is); ++ia, offset += npl)  {
         Coord r0Rel = gridCell.carToRel(str.getAtom(is,ia));

         // circumscribing box
         Coord f = gridCell.getBoundingBox (rCut(is));
         f += 0.001; // safety range for numerical noise at grid cell boundaries
         
         // range of grid cells (mesh vector)
         SxVector3<Int> from, to, xyz;
         from = ceil(r0Rel-f);
         to   = floor(r0Rel+f);

         // --- set up 1D phase factors
         SX_START_TIMER (RsPhase);
         SxVector<Complex16> phase1D((to-from).sum () +3);
         int off0 = 0,
             off1 = to(0) - from(0) + 1,
             off2 = off1 + to(1) - from(1) + 1;
         off0 -= from(0);
         off1 -= from(1);
         off2 -= from(2);

         xyz.set (0);
         for (xyz(0) = from(0); xyz(0) <= to(0); ++xyz(0))  {
            Coord rVec = gridCell.relToCar (xyz);
            phase1D(off0+xyz(0)) = exp (I * (rVec ^ kVec));
         }
         xyz(0) = 0;
         for (xyz(1) = from(1); xyz(1) <= to(1); ++xyz(1))  {
            Coord rVec = gridCell.relToCar (xyz);
            phase1D(off1+xyz(1)) = exp (I * (rVec ^ kVec));
         }
         xyz(1) = 0;
         for (xyz(2) = from(2); xyz(2) <= to(2); ++xyz(2))  {
            Coord rVec = gridCell.relToCar (xyz);
            phase1D(off2+xyz(2)) = exp (I * (rVec ^ kVec));
         }
         SX_STOP_TIMER (RsPhase);

         p2.set (0.);
         int irs = 0, ipls = 0;
         SxArray<int> idxU(maxnrs);
         for (xyz(0) = from(0); xyz(0) <= to(0); ++xyz(0))  {
            for (xyz(1) = from(1); xyz(1) <= to(1); ++xyz(1))  {
               for (xyz(2) = from(2); xyz(2) <= to(2); ++xyz(2))  {
                  Coord rVec = gridCell.relToCar (xyz);
                  Coord rel = rVec - str.constRef (is, ia);
                  double rAbs = rel.norm ();
                  if (rAbs > rCut(is)) continue;

                  if (rAbs > 1e-10)  {
                     //SX_CLOCK (RsYlm);
                     SxYlm::getYlmArray (lmax(is), rel, &ylm);
                  } else  {
                     ylm.set (0.);
                     ylm(0) = 1.;
                  }

                  SxComplex16 phase = phase1D(off0 + xyz(0))
                                    * phase1D(off1 + xyz(1))
                                    * phase1D(off2 + xyz(2));

                  idxU(irs) = (int)mesh.getMeshIdx(xyz, SxMesh3D::Unknown);

                  // --- compute f(r) for all projector types
                  double ir = rAbs / dr + nrn;
                  {
                     //SX_CLOCK (RsSplines);
                     double t = floor(ir);
                     int iir = int(t);
                     t = ir - t;
                     double u = 1. - t, t2 = t*t - 1., u2 = u*u - 1.;
                     int isc = npt * iir * 2;
                     for (int ipt = 0; ipt < npt; ++ipt, isc+=2)
                        rprAll(ipt)  = (sCoeff(isc) + u2 * sCoeff(isc+1)) * u;
                     for (int ipt = 0; ipt < npt; ++ipt, isc+=2)
                        rprAll(ipt) += (sCoeff(isc) + t2 * sCoeff(isc+1)) * t;
                  }
                  // --- compute projector phi_m(r_i) = f(r) * ylm(r)
                  int ipt = 0;
                  for (int ipl = 0; ipl < npl; )  {
                     //const ProjQuantumNumbers &id = projId(is)(ipl);
                     const ProjQuantumNumbers &id = pId[ipl];
                     //rpr = rProj(is)(id.n).getValYExtra (ir);
                     //double rpr = rprAll(id.n);
                     SxComplex16 rpr = rprAll(ipt++) * phase;
                     int lm = SxYlm::combineLm(id.l,-id.l);
                     for (int m = -id.l; m <= id.l; ++m) {
                        //phi(ipls++) = ylm(lm++) * rpr * phase;
                        phi(ipls++) = ylm(lm++) * rpr;
                     }
                     ipl += 2*id.l + 1;
                  }
                  // --- sum psi(n, r) phi(r, m) for several r's
                  if (++irs == maxnrs)  {
                     // --- collect unk
                     SX_START_TIMER (RsCollect);
                     for (int iState = 0, psiKidx = 0; iState < nStates; ++iState)  {
                        const SxDiracVec<Complex16> &uI = uk(iState);
                        for (int jrs = 0; jrs < irs; ++jrs)
                           psiK(psiKidx++) = uI(idxU(jrs));
                     }
                     SX_STOP_TIMER (RsCollect);
                     SX_CLOCK (RsMatSumProj);
                     p2 += phi ^ psiK;
                     irs = 0;
                     ipls = 0;
                  }
               }
            }
         }
         if (irs > 0)  {
            SxMatrix<Complex16> phi2 = phi(SxIdx(0, irs * npl-1));
            phi2.reshape (npl, irs);
            SX_START_TIMER (RsCollect);
            SxMatrix<Complex16> psiK2(irs, nStates);
            for (int iState = 0, idx = 0, idx2 = 0; iState < nStates; ++iState, idx2 += maxnrs-irs)  {
               const SxDiracVec<Complex16> &uI = uk(iState);
               for (int jrs = 0; jrs < irs; ++jrs)
                  //psiK2(idx++) = psiK(idx2++);
                  psiK2(idx++) = uI(idxU(jrs));
            }
            SX_STOP_TIMER (RsCollect);
            SX_CLOCK (RsMatSumProj);
            p2 += phi2 ^ psiK2;
         }
         for (int iState = 0; iState < nStates; ++iState) {
            for (int ipl = 0; ipl < npl; ++ipl)  {
               const ProjQuantumNumbers &id = projId(is)(ipl);
               SxComplex16 p = SxYlm::getYlmNormFactor (id.l,id.m) 
                             * gridCell.volume
                             * p2(ipl, iState);
               // missing (-i)^l factor
               if (id.l & 1) p *= -I;
               if (id.l & 2) p = -p;
               res(offset + ipl, iState) = p;
            }
         }
      }
   }
   SX_CHECK (offset == npg, offset, npg);
   return res;
}


SxDiracVec<Complex16> 
SxRSProj::project (const SxDiracVec<TPrecCoeffG> &psi) const
{
   SX_CLOCK (RsProject);
   const SxGBasis &gk = psi.getBasis<SxGBasis> ();
   SX_CHECK (gk.structPtr);
   SX_CHECK (rBasisPtr);
   if (antiGauss.getSize () == gk.ng && antiGauss.getBasisPtr () == &gk)  {
      // empty
   } else {
      antiGauss = exp(gk.g2 * (beta*beta*0.25));
      antiGauss.setBasis (&gk);
   }
   Coord kVec = gk.getK ();
   /*
   if (psi.nCols () <= 1)
      return project (*gk.structPtr, *rBasisPtr | (antiGauss * psi), kVec);
   */
   int nStates = (int)psi.nCols ();
   if (nStates <= nPsiMax)  {
      SxArray<SxDiracVec<TPrecCoeffG> > psiR(nStates);
      for (int iState = 0; iState < nStates; ++iState)  {
         psiR(iState) = *rBasisPtr | (antiGauss * psi.colRef(iState));
      }
      return project (*gk.structPtr, psiR, kVec);
   } else {
      SxArray<SxDiracVec<TPrecCoeffG> > psiR(nPsiMax);
      SxDiracMat<TPrecCoeffG> res(npg, nStates);
      int iStateBlock = 0, nBlock = nPsiMax;
      for (int iState = 0; iState < nStates; /* done inside */)  {
         psiR(iStateBlock++) = *rBasisPtr | (antiGauss * psi.colRef(iState++));
         if (iStateBlock == nBlock)  {
            SxDiracVec<Complex16> resBlock 
               = res(SxIdx((iState-nBlock) * npg, iState*npg -1));
            resBlock.reshape (npg, nBlock);
            resBlock <<= project(*gk.structPtr, psiR, kVec);
            nBlock = min(nPsiMax, nStates - iState);
            psiR.resize (nBlock);
            iStateBlock = 0;
         }
      }
      return res;
   }
}

SxDiracVec<Complex16> 
SxRSProj::gradient (const SxAtomicStructure &str,
                    const SxDiracVec<Complex16> &pPsi,
                    const Coord &kVec) const
{
   SX_CLOCK (RsGradSum);
   SX_CHECK (str.getNSpecies () == projId.getSize (),
             str.getNSpecies (), projId.getSize ());

   const SxRBasis &rBasis = *rBasisPtr;
   SX_CHECK ((rBasis.cell - str.cell).absSqr ().sum () < 1e-10);
   SxMesh3D mesh = rBasis.getMesh ();
   
   SxCell gridCell (str.cell.basis(0) / mesh(0),
                    str.cell.basis(1) / mesh(1),
                    str.cell.basis(2) / mesh(2));

   SxDiracVec<Complex16> res(*rBasisPtr);
   res.set (0.);

   int offset = 0;
   for (int is = 0; is < projId.getSize (); ++is)  {
      int npl = (int)projId(is).getSize ();
      SxVector<Double> ylm(sqr(lmax(is)+1)), work(lmax(is)+1);
      SxArray<SxComplex16> p2(npl);
      int lMax = lmax(is);
      ProjQuantumNumbers *pId = projId(is).elements;

      for (int ia = 0; ia < str.getNAtoms(is); ++ia, offset += npl)  {
         // --- collect projections with prefactors
         for (int ipl = 0; ipl < npl; ++ipl)  {
            const ProjQuantumNumbers &id = projId(is)(ipl);
            p2(ipl) = pPsi(offset + ipl) * SxYlm::getYlmNormFactor (id.l,id.m);
            // missing (i)^l factor
            if (id.l & 1) p2(ipl) *= I;
            if (id.l & 2) p2(ipl) = -p2(ipl); 
         }

         Coord r0Rel = gridCell.carToRel(str.getAtom(is,ia));

         // circumscribing box
         Coord f = gridCell.getBoundingBox (rCut(is));
         f += 0.001; // safety range for numerical noise at grid cell boundaries
         
         // range of grid cells (mesh vector)
         SxVector3<Int> from, to, xyz;
         from = ceil(r0Rel-f);
         to   = floor(r0Rel+f);

         // --- set up 1D phase factors
         SX_START_TIMER (RsPhase);
         SxVector<Complex16> phase1D((to-from).sum () +3);
         int off0 = 0,
             off1 = to(0) - from(0) + 1,
             off2 = off1 + to(1) - from(1) + 1;
         off0 -= from(0);
         off1 -= from(1);
         off2 -= from(2);

         xyz.set (0);
         for (xyz(0) = from(0); xyz(0) <= to(0); ++xyz(0))  {
            Coord rVec = gridCell.relToCar (xyz);
            phase1D(off0+xyz(0)) = exp (-I * (rVec ^ kVec));
         }
         xyz(0) = 0;
         for (xyz(1) = from(1); xyz(1) <= to(1); ++xyz(1))  {
            Coord rVec = gridCell.relToCar (xyz);
            phase1D(off1+xyz(1)) = exp (-I * (rVec ^ kVec));
         }
         xyz(1) = 0;
         for (xyz(2) = from(2); xyz(2) <= to(2); ++xyz(2))  {
            Coord rVec = gridCell.relToCar (xyz);
            phase1D(off2+xyz(2)) = exp (-I * (rVec ^ kVec));
         }
         SX_STOP_TIMER (RsPhase);

         for (xyz(0) = from(0); xyz(0) <= to(0); ++xyz(0))  {
            for (xyz(1) = from(1); xyz(1) <= to(1); ++xyz(1))  {
               for (xyz(2) = from(2); xyz(2) <= to(2); ++xyz(2))  {
                  Coord rVec = gridCell.relToCar (xyz);
                  Coord rel = rVec - str.constRef (is, ia);
                  double rAbs = rel.norm ();
                  if (rAbs > rCut(is)) continue;

                  if (rAbs > 1e-10)  {
                     //SX_CLOCK (RsYlm);
                     //if (lmax(is) == 1)  {
                     //   double ri = 1. / rAbs;
                     //   ylm(1) = SQRT2 * rel(1) * ri;
                     //   ylm(2) = rel(2) * ri;
                     //   ylm(3) = SQRT2 * rel(0) * ri;
                     //} else {
                        SxYlm::getYlmArray (lMax, rel, &ylm);
                     //}
                  } else  {
                     ylm.set (0.);
                     ylm(0) = 1.;
                  }

                  SxComplex16 psiK(0., 0.), phase;

                  double ir = rAbs / dr + nrn;
                  double rpr = 0.1;
                  /*
                  for (int ipl = 0; ipl < npl; ++ipl)  {
                     //const ProjQuantumNumbers &id = projId(is)(ipl);
                     const ProjQuantumNumbers &id = pId[ipl];
                     if (id.m == -id.l) rpr = rProj(is)(id.n).getVal (ir);
                     psiK += ylm(SxYlm::combineLm(id.l,id.m)) * rpr * p2(ipl);
                  }
                  */
                  for (int ipl = 0; ipl < npl; )  {
                     //const ProjQuantumNumbers &id = projId(is)(ipl);
                     const ProjQuantumNumbers &id = pId[ipl];
                     rpr = rProj(is)(id.n).getValYExtra (ir);
                     double *ylmp = ylm.elements + SxYlm::combineLm(id.l,0);
                     for (int m = -id.l; m <= id.l; ++m, ++ipl) {
                        psiK += ylmp[m] * rpr * p2(ipl);
                     }
                  }
                  phase = phase1D(off0 + xyz(0))
                        * phase1D(off1 + xyz(1))
                        * phase1D(off2 + xyz(2));
                  res(mesh.getMeshIdx(xyz, SxMesh3D::Unknown)) += psiK * phase;
               }
            }
         }
      }
   }
   SX_CHECK (offset == npg, offset, npg);
   return res;
}

SxArray<SxDiracVec<Complex16> >
SxRSProj::gradientN (const SxAtomicStructure &str,
                     const SxDiracVec<Complex16> &pPsi,
                     const Coord &kVec) const
{
   SX_CLOCK (RsGradSum);
   SX_CHECK (str.getNSpecies () == projId.getSize (),
             str.getNSpecies (), projId.getSize ());

   const SxRBasis &rBasis = *rBasisPtr;
   SX_CHECK ((rBasis.cell - str.cell).absSqr ().sum () < 1e-10);
   SxMesh3D mesh = rBasis.getMesh ();
   
   SxCell gridCell (str.cell.basis(0) / mesh(0),
                    str.cell.basis(1) / mesh(1),
                    str.cell.basis(2) / mesh(2));

   int nStates = (int)pPsi.nCols ();
   const int maxnrs = SX_RS_MAX_N_R;

   SxArray<SxDiracVec<Complex16> > res(nStates);
   for (int iState = 0; iState < nStates; ++iState)  {
      res(iState) = SxDiracVec<Complex16> (*rBasisPtr);
      res(iState).set (0.);
   }
   SxMatrix<Complex16> psiK(maxnrs, nStates);

   int offset = 0;
   for (int is = 0; is < projId.getSize (); ++is)  {
      int npl = (int)projId(is).getSize ();
      SxVector<Double> ylm(sqr(lmax(is)+1)), work(lmax(is)+1);
      SxMatrix<Complex16> p2(npl, nStates);
      SxMatrix<Complex16> phi(npl, maxnrs);
      int lMax = lmax(is);
      ProjQuantumNumbers *pId = projId(is).elements;

      for (int ia = 0; ia < str.getNAtoms(is); ++ia, offset += npl)  {
         // --- collect projections with prefactors
         for (int iState = 0; iState < nStates; ++iState)  {
            for (int ipl = 0; ipl < npl; ++ipl)  {
               const ProjQuantumNumbers &id = projId(is)(ipl);
                SxComplex16 p = pPsi(offset + ipl, iState)
                              * SxYlm::getYlmNormFactor (id.l,id.m);
               // missing (i)^l factor
               if (id.l & 1) p *= I;
               if (id.l & 2) p = -p; 
               p2(ipl, iState) = p;
            }
         }

         Coord r0Rel = gridCell.carToRel(str.getAtom(is,ia));

         // circumscribing box
         Coord f = gridCell.getBoundingBox (rCut(is));
         f += 0.001; // safety range for numerical noise at grid cell boundaries
         
         // range of grid cells (mesh vector)
         SxVector3<Int> from, to, xyz;
         from = ceil(r0Rel-f);
         to   = floor(r0Rel+f);

         // --- set up 1D phase factors
         SX_START_TIMER(RsPhase);
         SxVector<Complex16> phase1D((to-from).sum () +3);
         int off0 = 0,
             off1 = to(0) - from(0) + 1,
             off2 = off1 + to(1) - from(1) + 1;
         off0 -= from(0);
         off1 -= from(1);
         off2 -= from(2);

         xyz.set (0);
         for (xyz(0) = from(0); xyz(0) <= to(0); ++xyz(0))  {
            Coord rVec = gridCell.relToCar (xyz);
            phase1D(off0+xyz(0)) = exp (I * (rVec ^ kVec));
         }
         xyz(0) = 0;
         for (xyz(1) = from(1); xyz(1) <= to(1); ++xyz(1))  {
            Coord rVec = gridCell.relToCar (xyz);
            phase1D(off1+xyz(1)) = exp (I * (rVec ^ kVec));
         }
         xyz(1) = 0;
         for (xyz(2) = from(2); xyz(2) <= to(2); ++xyz(2))  {
            Coord rVec = gridCell.relToCar (xyz);
            phase1D(off2+xyz(2)) = exp (I * (rVec ^ kVec));
         }
         SX_STOP_TIMER(RsPhase);

         int irs = 0, ipls = 0;
         SxArray<int> idxU(maxnrs);
         for (xyz(0) = from(0); xyz(0) <= to(0); ++xyz(0))  {
            for (xyz(1) = from(1); xyz(1) <= to(1); ++xyz(1))  {
               for (xyz(2) = from(2); xyz(2) <= to(2); ++xyz(2))  {
                  Coord rVec = gridCell.relToCar (xyz);
                  Coord rel = rVec - str.constRef (is, ia);
                  double rAbs = rel.norm ();
                  if (rAbs > rCut(is)) continue;

                  if (rAbs > 1e-10)  {
                     //SX_CLOCK (RsYlm);
                     //if (lmax(is) == 1)  {
                     //   double ri = 1. / rAbs;
                     //   ylm(1) = SQRT2 * rel(1) * ri;
                     //   ylm(2) = rel(2) * ri;
                     //   ylm(3) = SQRT2 * rel(0) * ri;
                     //} else {
                        SxYlm::getYlmArray (lMax, rel, &ylm);
                     //}
                  } else  {
                     ylm.set (0.);
                     ylm(0) = 1.;
                  }

                  SxComplex16 phase = phase1D(off0 + xyz(0))
                                    * phase1D(off1 + xyz(1))
                                    * phase1D(off2 + xyz(2));

                  int idx = (int)mesh.getMeshIdx(xyz, SxMesh3D::Unknown);
                  idxU(irs) = idx;

                  double ir = rAbs / dr + nrn;
                  double rpr = 0.1;
                  for (int ipl = 0; ipl < npl; )  {
                     //const ProjQuantumNumbers &id = projId(is)(ipl);
                     const ProjQuantumNumbers &id = pId[ipl];
                     rpr = rProj(is)(id.n).getValYExtra (ir);
                     int lm = SxYlm::combineLm(id.l,-id.l);
                     for (int m = -id.l; m <= id.l; ++m, ++ipl) {
                        phi(ipls++) = ylm(lm++) * rpr * phase;
                     }
                  }
                  // --- sum phi(r, m) p(m, n) for several r's
                  if (++irs == maxnrs)  {
                     {
                        SX_CLOCK (RsMatSumGrad);
                        psiK = phi.overlap (p2);
                     }
                     // --- distribute to unk
                     SX_START_TIMER (RsDistribute);
                     for (int iState = 0, psiKidx = 0; iState < nStates; ++iState)  {
                        SxDiracVec<Complex16> &uI = res(iState);
                        for (int jrs = 0; jrs < irs; ++jrs)
                           uI(idxU(jrs)) += psiK(psiKidx++);
                     }
                     SX_STOP_TIMER (RsDistribute);
                     irs = 0;
                     ipls = 0;
                  }

               }
            }
         }
         if (irs > 0)  {
            {
               SxMatrix<Complex16> phi2 = phi(SxIdx(0, irs * npl-1));
               phi2.reshape (npl, irs);
               SX_CLOCK (RsMatSumGrad);
               psiK = phi2.overlap (p2);
            }
            // --- distribute to unk
            SX_START_TIMER (RsDistribute);
            for (int iState = 0, psiKidx = 0; iState < nStates; ++iState)  {
               SxDiracVec<Complex16> &uI = res(iState);
               for (int jrs = 0; jrs < irs; ++jrs)
                  uI(idxU(jrs)) += psiK(psiKidx++);
            }
            SX_STOP_TIMER (RsDistribute);
            irs = 0;
            ipls = 0;
         }

      }
   }
   SX_CHECK (offset == npg, offset, npg);
#ifndef NDEBUG
   for (int iState = 0; iState < nStates; ++iState)  {
      VALIDATE_VECTOR (res(iState));
   }
#endif

   return res;
}

SxDiracVec<TPrecCoeffG> 
SxRSProj::gradient (const SxDiracVec<Complex16> &pPsi,
                    const SxGBasis &gk) const
{
   SX_CLOCK (RsGrad);
   SX_CHECK (gk.structPtr);
   if (antiGauss.getSize () == gk.ng && antiGauss.getBasisPtr () == &gk)  {
      // empty
   } else {
      antiGauss = exp(gk.g2 * (beta*beta*0.25));
      antiGauss.setBasis (&gk);
   }
   Coord kVec = gk.getK ();
   if (pPsi.nCols () <= 1)
      return (gk | gradient(*gk.structPtr, pPsi, kVec)) * antiGauss;
   int nStates = (int)pPsi.nCols ();
   SxDiracMat<TPrecCoeffG> res(gk.ng, nStates);
   res.setBasis (&gk);
   int iStateBlock = 0, nBlock = 0;
   SxArray<SxDiracVec<Complex16> > psiR;
   for (int iState = 0; iState < nStates; ++iState, iStateBlock++)  {
      if (iStateBlock == nBlock)  {
         nBlock = min(nPsiMax, nStates-iState);
         SxDiracVec<Complex16> pPsiBlock 
            = pPsi(SxIdx(iState * npg, (iState+nBlock)*npg -1));
         pPsiBlock.reshape (npg, nBlock);
         psiR = gradientN (*gk.structPtr, pPsiBlock, kVec);
         iStateBlock = 0;
      }
      res.colRef(iState) <<= (gk | psiR(iStateBlock)) * antiGauss;
   }
   return res;
}

