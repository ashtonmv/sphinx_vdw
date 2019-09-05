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

#include <SxHubbardU.h>
#include <SxPAWPot.h>
#include <SxSimpleParser.h>
#include <SxProjector.h>
#include <SxHubbardMO.h>

// Reference 1: SPHInX Hubbard U implementation notes
// Reference 2: M. Cococcioni, S. de Gironcoli, Phys. Rev. B 71, 035105 (2005)

SxHubbardU::SxHubbardU ()
   : energy(0.), eDoubleCounting (0.)
{
   // empty
}

SxVector<Complex16>
SxHubbardU::computeIncr (int iSite, const SxVector<Complex16> &nij)
{
   SX_CHECK (nij.nCols () == nij.nRows (),
             nij.nCols (), nij.nRows ());
   int nOrbital = (int)nij.nRows ();
   SX_CHECK(nOrbital > 0, nOrbital);
   SxVector<Complex16> res;
   res.reformat (nOrbital, nOrbital);

   VALIDATE_VECTOR(nij);
   //SxMatrix<Complex16>::Eigensystem eig
   //   = SxMatrix<Complex16>(nij).eigensystem ();
   //cout << "nij ev=" << eig.vals << endl;

   double alpha = (alphaBias.getSize () > 0) ? alphaBias(iSite) : 0.;

   //cout << "alpha=" << alpha << endl;
   SX_LOOP(i)  {
      // ref 1, Eq. (1), ref. 2, Eq. (9)  (diagonal-only part)
      energy += 0.5 * Ueff(iSite) * nij(i,i).re;
      SX_LOOP(j)  {
         // ref 1, Eq. (1), ref. 2, Eq. (9)  (full matrix part)
         energy -= 0.5 * Ueff(iSite) * (nij(i,j) * nij(j,i)).re;

         // ref 1, Eq. (9)  (full matrix part)
         res(i,j) = -0.5 * Ueff(iSite) * (nij(i,j) + nij(j,i)).re;

         // ref 1, Eq. (13)
         eDoubleCounting += 0.5 * Ueff(iSite) * (nij(i,j) * nij(i,j)).re;
      }
      // ref 1, Eq. (9)  (diagonal-only part)
      res(i,i) += 0.5 * Ueff(iSite);
      // alpha contribution to energy ref.2 , Eq. (17)
      energy += alpha * nij(i,i).re;
      // ref 2, Eq. (18)
      res(i,i) += alpha;
      // ref. 2, Eq.(7)
      siteOccupation(iSite) += nij(i,i).re; 
   }
   //cout << (eig.vecs.adjoint () ^ res ^  eig.vecs) * HA2EV << endl;
   //sxprintf ("Hubbard site occupation(%d)=%.8f\n", iSite+1, siteOccupation);
   return res;
}

void SxHubbardU::resize (ssize_t nSite)
{
   Ueff.resize (nSite, true);
   siteOccupation.resize (nSite, true);
   if (alphaBias.getSize () > 0)
      alphaBias.resize (nSite, true);
   if (nSite == 0)  {
      atomicSiteMap.resize (0);
   }
}

void SxHubbardU::setZero ()
{
   energy = eDoubleCounting = 0.;
   siteOccupation.set (0.);
}

SxVector<Double> SxHubbardU::computeAtom(int iTl, const SxVector<Double> &Dij,
                                           double fFull)
{
   SX_CHECK (Dij.nRows () == Dij.nCols (),
             Dij.nRows (), Dij.nCols ());
   ssize_t npl = Dij.nRows ();
   SxVector<Double> Vij;
   Vij.reformat (npl, npl);
   Vij.set (0.);
   // --- to take into account non-spinpolarized cases, we do the following:
   // - Dij contains the summation over spins -> scale by 1/2 to get nij
   // - Vij is fine (acts on each spin the same way)
   // - energy must include spin summation. This is done by changing for
   //   each atom first to the single-spin (scale energy by 1/2), do the
   //   incremental compute, and then do the spin summation (multiply by 2).
   //   For the spin-polarized case, fFull is 1, and the scaling does not
   //   do anything
   double fInv = 1. / fFull;

   // rescale energy to single spin
   energy *= fInv;
   eDoubleCounting *= fInv;

   // --- run over orbital types with Hubbard U
   for (int iType = 0; iType < atomicSiteMap(iTl).getSize (); ++iType)  {
      const AtomicSite& site = atomicSiteMap(iTl)(iType);
      int nOrb = 2 * site.l + 1;
      SX_CHECK(nOrb > 0, nOrb, site.l);
      SxVector<Complex16> nik, HijSite;
      nik.reformat (nOrb, nOrb);
      // extract the occupation matrix for this orbital
      // Ref. 1, Eq. (6)
      SX_LOOP2(i,j)
         nik(i,j) = fInv * Dij(site.offset + i, site.offset + j) * site.Pnl;
      // do the Hubbard U incremental compute
      HijSite = computeIncr(site.iSite, nik);
      SX_CHECK(HijSite.imag ().normSqr () < 1e-20);
      // Ref. 1, Eq. (12)
      SX_LOOP2(i,j)
         Vij(site.offset + i, site.offset + j) += HijSite(i,j).re * site.Pnl;
   }
   // rescale energy to full occupation
   energy *= fFull;
   eDoubleCounting *= fFull;
   return Vij;
}

/// AtomicSite has empty constructor/destructor
SXSTACK_SIMPLE(SxHubbardU::AtomicSite)

static void inputProblem (const SxSymbolTable *table)
{
   cout << "Hubbard U in " 
        << table->parserFilename
        << " line " << table->parserLineNumber
        << ":";
}

void SxHubbardU::read (const SxSymbolTable *table,
                       const SxPtr<SxPAWPot> &pawPotPtr,
                       const SxAtomicStructure &structure)
{
   SX_CHECK (pawPotPtr);
   const SxPAWPot &pawPot = *pawPotPtr;
   SxArray<SxStack<AtomicSite> > sites(structure.getNAtoms ());
   ssize_t nSites = getSize ();
#ifndef NDEBUG
   if (nSites > 0)  {
      // --- we override the atomic site map
      SX_CHECK(atomicSiteMap.getSize () == 0, atomicSiteMap.getSize ());
      cout << "In " << __FILE__ << ":" << __LINE__ << endl;
      cout << "WARNING: calling SxHubbardU::read with some previous data"
           << endl;
      // this may be ok if you initialize the atomic sites after other
      // sites. If you need to reinitialize, use resize (0) first
   }
#endif
   SxStack<double> Ustack;
   SxStack<double> bias;
   bool useBias = alphaBias.getSize () > 0;
   if (useBias)  {
      SX_LOOP(i) bias << alphaBias(i);
   }
   // for overriding (updating) settings
   SxList<int> updateId;
   SxList<double> updateU, updateBias; 

   SYMBOLPARSE(table)  {
      SYMBOLGROUP("HubbardU")  {
         bool verbose = SYMBOLGET("verbose").toBool ();
         FOREACH_SYMBOLGROUP("site")  {
            double Usite = SYMBOLGET("U");
            Usite /= HA2EV; // units is eV
            int ipt      = SYMBOLGET("projectorType");
            bool update  = SYMBOLGET("update").toBool ();
            ipt--;

            int is = -1;
            SxString label;

            // --- get species
            if (HAVE_SYMBOL("species"))  {
               is = SYMBOLGET("species");
               is--;
            } else if (HAVE_SYMBOL("element"))  {
               is = pawPot.find (SYMBOLGET("element"));
               if (is < 0)  {
                  cout << "Hubbard U: element '" 
                       << SYMBOLGET("element")->toString ()
                       << "' not found." << endl;
               }
            } else if (HAVE_SYMBOL("label"))  {
               label = SYMBOLGET("label");
               if (!structure.hasLabels ())  {
                  inputProblem (SYMBOLGROUP_TABLE);
                  cout << "No labels defined in structure" << endl;
                  SX_QUIT;
               }
               // --- find species id for this label (and check that there's
               //     a unique one)
               const SxArray<SxString> &labels = structure.getLabels ();
               SX_LOOP(iTl)  {
                  if (labels(iTl) == label)  {
                     if (is == -1)
                        is = structure.getISpecies ((int)iTl);
                     else if (structure.getISpecies ((int)iTl) != is)  {
                        inputProblem (SYMBOLGROUP_TABLE);
                        cout << "Multiple species for same label:"  << endl;
                        SX_LOOP(jTl) if (labels(jTl) == label)
                           cout << "Atom " << (jTl + 1) << ": species "
                                << (structure.getISpecies ((int)jTl) + 1) 
                                << endl;
                        SX_QUIT;
                     }
                  }
               }
            }
            // 
            double alpha = SYMBOLGET("shift") || 0.;
            alpha /= HA2EV;
            if (fabs(alpha) > 1e-10)  {
               if (!useBias)
                  for (int i = 0; i < nSites; ++i)
                     bias << 0.;
               useBias = true;
            }
            if (is < 0)  {
               inputProblem (SYMBOLGROUP_TABLE);
               cout << " incomplete site specification" << endl;
               cout << "Use species, element, or label" << endl;
            }

            // --- check ipt 
            if (ipt < 0 || ipt >= pawPot.getNProjType (is))  {
               inputProblem (SYMBOLGROUP_TABLE);
               cout << "Illegal projector type " << (ipt + 1)
                    << " for species " << (is + 1)
                    << ". Valid range is 1.." 
                    << pawPot.getNProjType(is) << endl;
               SX_QUIT;
            }

            AtomicSite site;
            site.l      = pawPot.lPhi(is)(ipt);
            site.offset = pawPot.offset(is)(ipt);
            // Ref. 1, Eq. (7)
            //site.Pnl       = 1./tr(pawPot.pPS(is).colRef(ipt).sqr ());
            site.Pnl       = 1.;
            cout << "Pnl=" << site.Pnl << endl;

            const SxArray<SxString> *labels = NULL;
            if (label.getSize () > 0) labels = &structure.getLabels ();
            SX_LOOP(iAtom)  {
               ssize_t iTl = structure.getIAtom (is, iAtom);
               // in case labels are used: check if label matches 
               if (labels && (*labels)(iTl) != label) continue;
               if (update)  {
                  if (sites(iTl).isEmpty () 
                      || sites(iTl).top ().offset != site.offset)
                  {
                     inputProblem (SYMBOLGROUP_TABLE);
                     cout << endl;
                     cout << "'update' flag is only allowed directly after a "
                             " generic setup for the same projector type."
                          << endl;
                     SX_QUIT
                  }
                  site.iSite = sites(iTl).top ().iSite;
                  updateId << site.iSite;
                  updateU << Usite;
                  updateBias << alpha;
               } else {
                  // add new atomic site 
                  site.iSite = int(nSites++);
                  sites(iTl) << site;
                  Ustack     << Usite;
                  if (useBias) bias << alpha;
               }
               if (verbose)  {
                  cout << "Site " << (site.iSite + 1) << ": projector "
                       << (ipt+1) << " (l=" << site.l << ") at atom "
                       << (iAtom+1) << " of species " << (is+1)
                       << " (total atom id: " << (iTl+1) << ") with U="
                       << (Usite * HA2EV) << " eV" << endl;
               }
            }
         }
         FOREACH_SYMBOLGROUP("MO")  {
            // create new MO site
            ssize_t iHubMO = moSite.getSize ();
            moSite.resize (iHubMO + 1, true);
            moSite(iHubMO) = SxPtr<SxHubbardMO>::create ((int)nSites);
            // read
            moSite(iHubMO)->read (SYMBOLGROUP_TABLE, structure, pawPotPtr);

            // read U and potential bias
            double Usite = SYMBOLGET("U");
            Usite /= HA2EV; // units is eV
            double alpha = SYMBOLGET("shift") || 0.;
            alpha /= HA2EV;
            if (fabs(alpha) > 1e-10)  {
               if (!useBias)
                  for (int i = 0; i < nSites; ++i)
                     bias << 0.;
               useBias = true;
            }

            // add U's to complete list of U's
            int nSitesMO = moSite(iHubMO)->getNSite ();
            for (int iSite = 0; iSite < nSitesMO; ++iSite)  {
               Ustack << Usite;
               if (useBias) bias << alpha;
            }

            nSites += nSitesMO;
         }
         FOREACH_SYMBOLGROUP("AO")  {
            // create new MO site
            ssize_t iHubMO = moSite.getSize ();
            moSite.resize (iHubMO + 1, true);
            moSite(iHubMO) = SxPtr<SxHubbardMO>::create ((int)nSites);
            // read
            moSite(iHubMO)->readAO (SYMBOLGROUP_TABLE, structure, pawPotPtr);

            // read U and potential bias
            double Usite = SYMBOLGET("U");
            Usite /= HA2EV; // units is eV
            double alpha = SYMBOLGET("shift") || 0.;
            alpha /= HA2EV;
            if (fabs(alpha) > 1e-10)  {
               if (!useBias)
                  for (int i = 0; i < nSites; ++i)
                     bias << 0.;
               useBias = true;
            }

            // add U's to complete list of U's
            int nSitesMO = moSite(iHubMO)->getNSite ();
            for (int iSite = 0; iSite < nSitesMO; ++iSite)  {
               Ustack << Usite;
               if (useBias) bias << alpha;
            }

            nSites += nSitesMO;
         }
      }
   }
   resize (nSites);
   atomicSiteMap.resize (structure.getNAtoms ());
   SX_LOOP(iTl) atomicSiteMap(iTl) = sites(iTl);
   // fill in effective U from stack from the end (to keep the
   // previous ones, if any)
   while (!Ustack.isEmpty ()) Ueff(--nSites) = Ustack.pop ();

   // set defined alphas
   if (useBias) alphaBias = bias;

   // implement updates (requires random access) 
   while (updateId.getSize () > 0)  {
      int iSite = updateId(0);
      Ueff(iSite) = updateU(0);
      if (useBias) alphaBias(iSite) = updateBias(0);
      updateId.removeFirst ();
      updateU.removeFirst ();
      updateBias.removeFirst ();
   }

}

void SxHubbardU::computeMO (const SxArray<SxPtr<SxBlockDensityMatrix> > &Pij,
                            const SxAtomicStructure &structure)
{
   for (int i = 0; i <moSite.getSize (); i++)  {
      moSite(i)->compute (this, *Pij(i), structure);
   }
}
