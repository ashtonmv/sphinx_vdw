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

#include <SxFermi.h>
#include <SxConstants.h>
//#include <fhing.h>
#include <stdio.h>
#include <SxLoopMPI.h>

SxFermi::SxFermi ()
   : nElectrons(-1),
     nSpin(-1),
     fFull(0.),
     orderMethfesselPaxton(-1),
     kpPtr (NULL)
{
   init ();
}


SxFermi::SxFermi (int nStates, int nSpin_, int nk)
   : orderMethfesselPaxton(-1)
{
   SX_CHECK (nStates > 0);
   SX_CHECK (nSpin_ > 0);
   nElectrons = -1.;
   nSpin      = nSpin_;
   fFull      = (nSpin == 1) ? 2. : 1.;
   init ();
   kpPtr        = NULL;
   eps          = Eps     (nStates, nSpin, nk);
   focc         = Focc    (nStates, nSpin, nk);
   epsSortList  = EpsSort (nStates, nSpin, nk);
   valBandIdx   = SxBundle3<Int> (nStates, nSpin, nk);
   conBandIdx   = SxBundle3<Int> (nStates, nSpin, nk);

   nValBands.resize (nk); nConBands.resize(nk);
   int i, ik, iSpin;
   for (ik=0; ik < nk; ik++)  {
      nValBands(ik).resize(nSpin);
      nConBands(ik).resize(nSpin);
   }

   for (ik=0; ik < nk; ik++)  {
      for (iSpin=0; iSpin < nSpin; iSpin++)  {
         for (i=0; i < nStates; i++)  {
            epsSortList(i,iSpin,ik) = i;
         }
      }
   }

}

SxFermi::SxFermi (Real8 nElectrons_, int nStates, int nSpin_,
                  const SxKPoints &gk)
   : orderMethfesselPaxton(-1)
{
   SX_CHECK (nStates > 0);
   SX_CHECK (nSpin_ > 0);
   nElectrons = nElectrons_;
   nSpin      = nSpin_;
   fFull      = (nSpin == 1) ? 2. : 1.;
   init ();
   kpPtr        = &gk;
   eps          = Eps     (nStates, nSpin, gk.nk);
   focc         = Focc    (nStates, nSpin, gk.nk);
   epsSortList  = EpsSort (nStates, nSpin, gk.nk);
   valBandIdx   = SxBundle3<Int> (nStates, nSpin, gk.nk);
   conBandIdx   = SxBundle3<Int> (nStates, nSpin, gk.nk);

   nValBands.resize (gk.nk); nConBands.resize(gk.nk);
   int i, ik, iSpin;
   for (ik=0; ik < gk.nk; ik++)  {
      nValBands(ik).resize(nSpin);
      nConBands(ik).resize(nSpin);
   }

   for (ik=0; ik < gk.nk; ik++)  {
      for (iSpin=0; iSpin < nSpin; iSpin++)  {
         for (i=0; i < nStates; i++)  {
            focc(i,iSpin,ik) = nElectrons / nStates / nSpin;
            eps(i,iSpin,ik)  = i;
            epsSortList(i,iSpin,ik) = i;
         }
      }
   }

}


void SxFermi::init ()
{
   printHartree     = false;
   eFermi           = 0.;
   spinMoment       = 0.;
   eqEntropy        = 0.;
   noneqEntropy     = 0.;
   keepSpinMoment   = false;
}


SxFermi::~SxFermi ()
{
   // empty
}

void SxFermi::setOccupencies (const Focc &foccIn)
{
   focc = foccIn;
   int ik, iSpin, i;
   for (ik=0; ik < kpPtr->nk; ik++)
      for (iSpin=0; iSpin < nSpin; iSpin++)
         for (i=0; i < focc.getNStates(ik); i++)
            epsSortList(i,iSpin,ik) = i;
   spinMoment = getSpinMoment ();
}

namespace { // C++ way to hide visibility outside this file
/// utility function for getIdxList
inline int translateIdx(int i, int nMax, const SxString &name)
{
   if (i < 0) i += nMax; // from end-indexing
   else       i--;      // change to C-indexing starting from 0
   if (i < 0 || i >= nMax)  {
      cout << "Illegal " << name << "-index " << i <<  "." << endl;
      SX_QUIT;
   }
   return i;
}
}

SxList<int> SxFermi::getIdxList (const SxSymbolTable *table,
                                 int nMax,
                                 const SxString &name)
{
   SxList<int> res;
   if (table->contains("values", true))  {
      res = table->get("values", true)->toIntList ();
      for (int i = 0; i < res.getSize (); ++i)
         res(i) = translateIdx(res(i), nMax, name);
   }
   if (table->contains("range", true))  {
      SxArray<int> range = table->get("range", true)->toIntList ();
      if (range.getSize () != 2)  {
         cout << name << "-range does not have [from,to] format." << endl;
         SX_QUIT;
      }
      for (int i = 0; i < 2; ++i)
         range(i) = translateIdx(range(i), nMax, name);
      if (range(0) > range(1))  {
         cout << name << "-range=[" << range(0) <<  "," << range(1)
              <<  "] has negative length" << endl;
         SX_QUIT;
      }
      res.resize (range(1)-range(0)+1);
      for (int j = 0, i = range(0); i <= range(1); ++i, ++j) res(j) = i;
   }
   return res;
}


void SxFermi::readOccupations (const SxSymbolTable *table)
{
   SX_CHECK (getNk () > 0);
   SX_CHECK (getNSpin () > 0);
   SX_CHECK (getNStates () > 0);
   SX_CHECK (nElectrons >= 0.);
   SX_CHECK (kpPtr);

   int nk = getNk (), nStates = getNStates ();

   const SxSymbolTable *kp, *spin, *band;
   focc.set (0.);
   try  {
      // get right table
      if (table->getName () != "occupations")
         table = table->getGroup ("occupations");

      // --- parse full-list
      if (table->contains ("values"))  {
         SxArray<int> vals = table->get("values")->toIntList ();
         if (nk * nSpin * nStates != vals.getSize ())  {
            cout << "Wrong number of occupations in input file!" << endl;
            cout << "nk = " << nk << endl;
            cout << "nSpin = " << nSpin << endl;
            cout << "nStates = " << nStates << endl;
            cout << "Expected " << (nk*nSpin*nStates)
                 << "occupation numbers, but found " << vals.getSize ()
                 << "." << endl;
            SX_QUIT;
         }
         int iVal = 0;
         for (int ik = 0; ik < nk; ++ik)
            for (int iSpin = 0; iSpin < nSpin; ++iSpin)
               for (int i = 0; i < nStates; ++i, ++iVal)
                  focc(i,iSpin,ik) = vals(iVal);
      } else {
         // --- parse kPoint-spin-bands groups
         SxList<int> kList, spinList, bandList;
         SxList<int>::Iterator kIt, spinIt, bandIt;

         if (table->containsGroup("kPoints"))  {
            kp = table->getGroup ("kPoints");
         } else {
            kp = table;
            for (int ik = 0; ik < nk; ++ik) kList.append (ik);
         }

         for ( ; kp; kp = (kp == table) ? NULL : kp->nextSibling ("kPoints"))
         {
            // get k indices for kPoints group
            if (kp != table) kList = getIdxList(kp, nk, "k");

            if (kp->containsGroup ("spin"))  {
               spin = kp->getGroup ("spin");
            } else {
               spin = kp;
               // no spin group -> all spins
               spinList.resize (0);
               for (int iSpin = 0; iSpin < nSpin; ++iSpin) spinList << iSpin;
            }

            for (; spin ;
                 spin = (spin == kp) ? NULL : spin->nextSibling ("spin"))
            {

               // get spin indices for spin group
               if (spin != kp) spinList = getIdxList (spin, nSpin, "spin");

               for (band = spin->getGroup ("bands"); band;
                     band = band->nextSibling ("bands"))
               {
                  bandList = getIdxList(band, nStates, "band");

                  // get occupation number
                  double foccVal = band->get ("focc")->toReal ();

                  // --- set specified occupation numbers to foccVal
                  for (kIt = kList.begin (); kIt != kList.end (); ++kIt)
                     for (spinIt = spinList.begin ();
                          spinIt != spinList.end (); ++spinIt)
                        for (bandIt = bandList.begin ();
                             bandIt != bandList.end (); ++bandIt)
                           focc(*bandIt, *spinIt, *kIt) = foccVal;
               }
            }
         }
      }

   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }

   // --- consistency check with nElectrons
   double sum = 0.;
   for (int ik = 0; ik < nk; ++ik)
      for (int iSpin = 0; iSpin < nSpin; ++iSpin)
         for (int i = 0; i < nStates; ++i)
            sum += focc(i,iSpin,ik) * kpPtr->weights(ik);
   if (fabs(sum - nElectrons) > 1e-4)  {
      printOccupation ();
      cout << "Wrong occupation number sum!" << endl
           << "expected " << nElectrons
           << " but found " << sum << ".\n";
      SX_QUIT;
   }

   spinMoment = getSpinMoment ();
}

int SxFermi::getNk () const
{
   SX_CHECK (focc.getNk () == eps.getNk (), focc.getNk (), eps.getNk ());
   return focc.getNk ();
}

int SxFermi::getNStates (int ik) const
{
   SX_CHECK (ik < getNk (), ik, getNk ());
   SX_CHECK (focc.getNStates (ik) == eps.getNStates (ik),
             focc.getNStates (ik), eps.getNStates (ik));
   return focc.getNStates (ik);
}

void SxFermi::setSpinMoment (Real8 spinMoment_, bool keepFixed_)
{
   spinMoment     = spinMoment_;
   keepSpinMoment = keepFixed_;
}


void SxFermi::fermiDistribution (Real8 ekt, Real8 mixFactor)
{
   SX_CHECK (mixFactor >= 0. && mixFactor <= 1.0, mixFactor);
   SX_CHECK (kpPtr);
   SX_CHECK (kpPtr->nk == getNk (), kpPtr->nk, getNk ());
   if (ekt < 1e-10)  ekt = 1e-10;

   Focc foccOld;
   int iSpin, ik;
   Real8 nEl = nElectrons;

   // --- save current occupencies in foccOld
   if (mixFactor < 1.)   foccOld = Focc (focc);

   if (nSpin != 1 && keepSpinMoment)  {
      eSpinFermi[SPIN_UP]   = getFermiEnergy (ekt, nEl/2.+spinMoment/2., SPIN_UP);
      eSpinFermi[SPIN_DOWN] = getFermiEnergy (ekt, nEl/2.-spinMoment/2., SPIN_DOWN);
   }  else  {
      eFermi = getFermiEnergy (ekt, nEl);
   }

   // --- equilibrium entropy
   eqEntropy = getEntropy ();

   // --- mix new and old occupation
   if (mixFactor < 1.)  {
      for (ik=0; ik < kpPtr->nk; ik++)
         for (iSpin=0; iSpin < nSpin; iSpin++)
            focc(iSpin,ik) <<= (1.-mixFactor) * foccOld(iSpin,ik)
               +     mixFactor  * focc(iSpin,ik);
   }

   // --- nonequilibrium entropy
   noneqEntropy = getEntropy ();

   // --- semiconductors: update idx array for valence and conduction bands
   updateValCon ();


   checkOccupation ();

   if (!keepSpinMoment) spinMoment = getSpinMoment ();
}



double SxFermi::getEBand (enum FoccHandling mode) const
{
   SX_CHECK (kpPtr);
   const SxKPoints &gk = *kpPtr;
   SX_CHECK (mode == UseFocc || mode == NoFocc);
   double eBand = 0.;
   int ik, iSpin, i, nStates;
   if (mode == UseFocc)  {
      for (ik=0; ik < gk.nk; ik++)  {
         nStates = focc.getNStates (ik);
         for (iSpin=0; iSpin < nSpin; iSpin++) {
            for (i=0; i < nStates; i++)  {
               eBand += focc(i,iSpin,ik) * gk.weights(ik) * eps(i,iSpin,ik);
            }
         }
      }
   } else { // mode == NoFocc
      for (ik=0; ik < gk.nk; ik++)  {
         nStates = focc.getNStates (ik);
         for (iSpin=0; iSpin < nSpin; iSpin++) {
            for (i=0; i < nStates; i++)  {
               eBand += gk.weights(ik) * eps(i,iSpin,ik);
            }
         }
      }
   }
   return eBand;
}


void SxFermi::checkOccupation ()
{
   SX_CHECK (kpPtr);
   int iSpin, ik;
   Real8 sum = 0.;
   const SxKPoints &gk = *kpPtr;
   for (ik=0; ik < gk.nk; ik++)
      for (iSpin=0; iSpin < nSpin; iSpin++)
         sum += focc(iSpin,ik).sum() * gk.weights(ik);


   if ( fabs (sum - nElectrons) > 0.001 )  {
      printOccupation ();
      cout << "\nThe norm of occupation differs from the number ";
      cout << "of electrons.\n";
      cout << "Norm = " << sum;
      cout << "(should be " << nElectrons << ")\n";
      SX_EXIT;
   }
}


Real8 SxFermi::getSpinMoment ()
{
   SX_CHECK (kpPtr);
   int ik, iSpin=0;
   Real8 sum = 0.;

   const SxKPoints &gk = *kpPtr;
   for (ik=0; ik < gk.nk; ik++)
      sum += focc(iSpin,ik).sum() * gk.weights(ik);  // TODO: weight

   return (2.*sum - nElectrons);
}

Real8 SxFermi::getEntropy ()
{
   SX_CHECK (kpPtr);
   int i, ik, nStates, iSpin;
   Real8 sum = 0., iksum;
   Real8 s = fFull;

   for (ik=0; ik < kpPtr->nk; ik++)  {
      nStates = focc.getNStates (ik);
      iksum = 0.;
      for (iSpin=0; iSpin < nSpin; iSpin++) {
         if (orderMethfesselPaxton < 0)  {
            // --- Fermi-Dirac
            for (i=0; i < nStates; i++)  {
               double f = focc(i,iSpin,ik);
               if ( f > 0. && f < s && fabs(s-f) > 1e-15 )
                  iksum -= (s-f)*log(1.-f/s) + f*log(f/s);
            }
         } else {
            /// --- Methfessel-Paxton
            double eF = keepSpinMoment ? eSpinFermi[iSpin] : eFermi;
            // empirical relation: for this, entropy becomes < 1e-16
            double xLim = 6.5 + 0.5 * orderMethfesselPaxton;
            for (i=0; i < nStates; i++)  {
               double entropy = 0.;
               double x = (eps(i,iSpin,ik) - eF)*beta;
               if ( fabs(x) < xLim )  {
                  dfoccMethfesselPaxton (x, orderMethfesselPaxton, &entropy);
                  iksum += 0.5 * fFull * entropy;
               }
            }
         }
      }
      sum += iksum * kpPtr->weights(ik);
   }

   return sum;
}

//inline
Real8 SxFermi::foccMethfesselPaxton(Real8 x, int order)
{
   // ref 1: Methfessel, Paxton, Phys Rev. B. 40, 3616 (1989)
   double A = 1./SQRT_PI; // A0, ref 1, Eq. for A_n for n=0
   double f = 0.5 * derfc(x); // complementary error function S0 in ref1
   double dn = 1.;
   double H2n = exp(-x*x); // H0(x) * exp(-x^2)
   double H2n1 = 2. * x * H2n;  // H1(x) * exp(-x^2)
   for (int n = 1; n <= order; n++, dn+=1.)  {
      A = -0.25 * A / dn; // => A0 * (-1)^n/ (n! 4^n)
      // H2n1 is H_(2n-1)
      f += A * H2n1;
      // H2n is currently H_(2n-2)
      // now compute H_2n (Hermite polynomial recursion H_n+1
      H2n = 2. * (x * H2n1 - (2. * dn - 1.) * H2n);
      // now compute H_(2n+1)
      H2n1 = 2. * (x * H2n - (2. * dn) * H2n1);
   }
   return f;
}

//inline
Real8 SxFermi::dfoccMethfesselPaxton(Real8 x, int order, double *lastDelta)
{
   // ref 1: Methfessel, Paxton, Phys Rev. B. 40, 3616 (1989)
   double A = 1./SQRT_PI; // A0, ref 1, Eq. for A_n for n=0
   double dn = 1.;
   double H2n = exp(-x*x); // H0(x) * exp(-x^2)
   double H2n1 = 2. * x * H2n;  // H1(x) * exp(-x^2)
   double df = A * H2n;
   for (int n = 1; n <= order; n++, dn+=1.)  {
      A = -0.25 * A / dn; // => A0 * (-1)^n/ (n! 4^n)
      // H2n1 is H_(2n-1)
      // H2n is currently H_(2n-2)
      // now compute H_2n (Hermite polynomial recursion H_n+1
      H2n = 2. * (x * H2n1 - (2. * dn - 1.) * H2n);
      // now compute H_(2n+1)
      H2n1 = 2. * (x * H2n - (2. * dn) * H2n1);
      df += A * H2n;
   }
   if (lastDelta) *lastDelta = A * H2n; // needed for entropy
   return df;
}


Real8 SxFermi::fermiFunction (Real8 energy, Real8 nEl, int spin)
{
   SX_CHECK (kpPtr);
   int i, nStates, ik, iSpin, iSpin0, iSpin1;
   int nk = getNk ();
   Real8 xArg;
   Real8 f;
   Real8 s = fFull;


   if (   nSpin == 1 || !keepSpinMoment )   {
      iSpin0 = 0; iSpin1 = nSpin;
   }  else  {
      iSpin0 = spin; iSpin1 = spin+1; // i.e. no loop over iSpin
   }

   double maxArg = (orderMethfesselPaxton < 0)
                 ? log(1e-16) // Fermi-Dirac
                 : -(sqrt(-log(1e-16))+sqrt(fabs(orderMethfesselPaxton)));
   for (ik=0; ik < nk; ik++)  {
      nStates = focc.getNStates (ik);
      double kSum = 0.;
      for (iSpin=iSpin0; iSpin < iSpin1; iSpin++) {
         for (i=0; i < nStates; i++)  {
            xArg  = (eps(i,iSpin,ik) - energy) * beta;
            if (xArg < maxArg)
               f = s;
            else if (xArg > -maxArg)
               f = 0.;
            else {
               if (orderMethfesselPaxton < 0)  {
                  f = s / (exp(xArg)+1.);
               } else {
                  f = fFull * foccMethfesselPaxton(xArg, orderMethfesselPaxton);
               }
            }
            kSum += f;
            focc(i,iSpin,ik) = f;
         }
      }
      nEl -= kSum * kpPtr->weights(ik);
   }
   return -nEl;
}

Real8 SxFermi::dFoccFermi(int i, int iSpin, int ik) const
{
   if (orderMethfesselPaxton < 0)  {
      double f = focc(i, iSpin, ik);
      return f * (1. - f / fFull) * beta;
   } else {
      double eF = keepSpinMoment ? eSpinFermi[iSpin] : eFermi;
      double x = (eps(i, iSpin, ik) - eF) * beta;
      return fFull * dfoccMethfesselPaxton (x, orderMethfesselPaxton) * beta;
   }
}

Real8 SxFermi::getFermiEnergy (Real8 ekt, Real8 nEl, int spin)
{
   beta = (ekt > 0.) ? 1./ekt : 1e5 /* just a large value */;
   // for large numbers of kPoints numerical noise increases
   double epsF = 1e-12 * getNk ();
   double eMin = -100., eMax = 100.;
   double fMin, fMax;
   // --- find lower bound
   bool fail = false;
   while ((fMin = fermiFunction (eMin, nEl, spin)) > 0.) {
      if (fabs(fMin) < epsF) break;
      eMin *= 2.;
      if (eMin < -1e5) { fail = true; break; }
   }
   // --- find upper bound
   while ((fMax = fermiFunction (eMax, nEl, spin)) < 0.) {
      if (fabs(fMax) < epsF) break;
      eMax *= 2.;
      if (eMax > 1e5) { fail = true; break; }
   }
   // error on failure
   if (fail)  {
      cout << "eMin=" << eMin * HA2EV << " eV, Delta n=" << fMin << endl;
      cout << "eMax=" << eMax * HA2EV << " eV, Delta n=" << fMax << endl;
      cout << "Can't find root of fermi function!" << endl;
      SX_EXIT;
   }

   // --- interval section
   for (int it = 0; it < 200; it++)  {
      double eMiddle = 0.5 * (eMin + eMax);
      double fMiddle = fermiFunction (eMiddle, nEl, spin);
      if (fabs(fMiddle) < epsF) return eMiddle;
      if (fMiddle < 0.)  {
         eMin = eMiddle;
         fMin = fMiddle;
      } else {
         eMax = eMiddle;
         fMax = fMiddle;
      }
      if (fabs(eMin - eMax) < 1e-12 && fabs(fMin - fMax) < epsF)
         return eMiddle;
   }

   // if we get here, the following situation occurred:
   // - no isolator, i.e. not all k-points have same number of electrons
   // - k-points between which electrons are distributed have different weights
   // - xMin and xMax enclose one eigenvalue up to some accuracy delta x, which
   //   is determined by the numerical accuracy
   // - ekt is comparable to or smaller than delta x
   // - for Fermi energy = xMin < eigenvalue the state is (almost completely)
   //   unoccupied
   // - for Fermi energy = xMax > eigenvalue the state is (almost completely)
   //   occupied
   // - the difference in the fermiFunction between the two situations is
   //   weights (ik), so yMin and yMax may be below and above 0, respectively,
   //   because the other electrons are distributed in different portions of
   //   weights (jk)
   // Solution: increase ekt locally

   // --- to be on the safe side: exit if ekt is large
   if (ekt > 1. /* Hartree */)  {
      cout << SX_SEPARATOR;
      cout << "| Can't distribute electrons among k-points correctly" << endl;
      cout << "| even with huge electronic temperature (ekt = " << (ekt * HA2EV);
      cout << " eV)" << endl << SX_SEPARATOR;

      cout << "| remainder: " << fermiFunction (eMin, nEl, spin)
            << "@E=" << (eMin * HA2EV) << " eV" << endl;
      cout << "| remainder: " << fermiFunction (eMax, nEl, spin)
            << "@E=" << (eMax * HA2EV) << " eV" << endl;
      cout << "| Writing spectrum to eps-crash.dat" << endl;
      SX_MPI_MASTER_ONLY writeSpectrum ("eps-crash", "dat");
      SX_EXIT;
   }
   // jump from very low ekt to 1e-5 (just a number)
   // increase ekt by a factor of 10 if its not that tiny
   double newEkt = (ekt < 1e-6 / HA2EV) ? (1e-5 / HA2EV) : (ekt * 10.);
   cout << SX_SEPARATOR;
   cout << "| WARNING: Fermi distribution failed!" << endl;
   cout << "|          Can't distribute electrons among k-points with\n";
   cout << "|          ekt = " << ((ekt > 1e-9) ? (ekt * HA2EV) : 0.);
   cout << " eV. Trying ekt = " << (newEkt * HA2EV) << " eV now.\n";
   cout << SX_SEPARATOR;
   cout.flush ();
   return getFermiEnergy (newEkt, nEl, spin);
}


void SxFermi::resize (int nStates)
{
   SX_CHECK (getNk() > 0, getNk());
   SX_CHECK (nStates > 0, nStates);
   const PrecEps LARGE_EPS = 1e30;  // shift new states energetically away
   int iSpin, ik, nk = getNk ();
   for (ik=0; ik < nk; ik++)  {
      for (iSpin=0; iSpin < nSpin; iSpin++)  {
         eps.bundle(ik)(iSpin).resize         (nStates, true, LARGE_EPS);
         focc.bundle(ik)(iSpin).resize        (nStates, true);
         epsSortList.bundle(ik)(iSpin).resize (nStates, true);
         valBandIdx.bundle(ik)(iSpin).resize  (nStates, true);
         conBandIdx.bundle(ik)(iSpin).resize  (nStates, true);
      }
   }

   // --- sort eigenvalues (for printing)
   for (ik=0; ik < nk; ik++)
      for (iSpin=0; iSpin < nSpin; iSpin++)
         epsSortList(iSpin,ik) = eps(iSpin,ik).getSortIdx();

   updateValCon ();
}


void SxFermi::resize (const SxArray<int> &nPerK)
{
   SX_CHECK (getNk() > 0, getNk());
   int iSpin, ik, n, nk = getNk ();
   const PrecEps LARGE_EPS = 1e30;  // shift new states energetically away
   for (ik=0; ik < nk; ik++)  {
      n = nPerK (ik);
      SX_CHECK (n > 0, n);
      for (iSpin=0; iSpin < nSpin; iSpin++)  {
         eps.bundle(ik)(iSpin).resize         (n, true, LARGE_EPS);
         focc.bundle(ik)(iSpin).resize        (n, true);
         epsSortList.bundle(ik)(iSpin).resize (n, true);
         valBandIdx.bundle(ik)(iSpin).resize  (n, true);
         conBandIdx.bundle(ik)(iSpin).resize  (n, true);
      }
   }

   // --- sort eigenvalues (for printing)
   for (ik=0; ik < nk; ik++)
      for (iSpin=0; iSpin < nSpin; iSpin++)
         epsSortList(iSpin,ik) = eps(iSpin,ik).getSortIdx();

   updateValCon ();
}


void SxFermi::updateValCon (double fOccMin, double fUnoccMax)
{
   int i, nStates, iSpin, ik;
   int nc, nv;
   PrecFocc f;

   fUnoccMax *= 2./fFull;

   for (ik=0; ik < getNk (); ik++)  {
      nStates = focc.getNStates (ik);
      for (iSpin=0; iSpin < nSpin; iSpin++)  {
         valBandIdx(iSpin,ik).set (-1);
         conBandIdx(iSpin,ik).set (-1);
         for (i=0, nc=0, nv=0; i < nStates; i++)  {
            f = focc(i,iSpin,ik);

            if (f > fOccMin)   valBandIdx(nv++,iSpin,ik) = i;
            if (f < fUnoccMax) conBandIdx(nc++,iSpin,ik) = i;
         }
         nValBands(ik)(iSpin) = nv;
         nConBands(ik)(iSpin) = nc;
      }
   }
}


int SxFermi::getNValenceBands (int iSpin, int ik) const
{
   int nv = nValBands(ik)(iSpin);
   SX_CHECK (nv > 0 && nv <= focc.getNStates(ik), nv, focc.getNStates(ik));
   return nv;
}

int SxFermi::getNValenceBands () const
{
   int nv = 0;
   for (int ik = 0; ik < getNk (); ++ik)
      for (int iSpin = 0; iSpin < nSpin; ++iSpin)
         if (nv < nValBands(ik)(iSpin)) nv = nValBands(ik)(iSpin);
   return nv;
}

int SxFermi::getNConductionBands () const
{
   int nc = 0;
   for (int ik = 0; ik < getNk (); ++ik)
      for (int iSpin = 0; iSpin < nSpin; ++iSpin)
         if (nc < nConBands(ik)(iSpin)) nc = nConBands(ik)(iSpin);
   return nc;
}


int SxFermi::getNConductionBands (int iSpin, int ik) const
{
   int nc = nConBands(ik)(iSpin);
   SX_CHECK (nc >= 0 && nc < focc.getNStates(ik), nc, focc.getNStates(ik));
   return nc;
}


int SxFermi::getValenceBandIdx (int i, int iSpin, int ik) const
{
   return valBandIdx (i,iSpin,ik);
}


int SxFermi::getConductionBandIdx (int i, int iSpin, int ik) const
{
   return conBandIdx(i,iSpin,ik);
}


bool SxFermi::isSemiconductor ()
{
   int nVal, ik;
   int nk = getNk ();
   SX_CHECK (nk > 0, nk);
   for (int iSpin = 0; iSpin < getNSpin (); iSpin++)  {
      nVal = nValBands(0)(iSpin);
      for (ik = 1; ik < nk; ik++)  {
         if (nValBands(ik)(iSpin) != nVal) return false;
         if (nValBands(ik)(iSpin) + nConBands(ik)(iSpin) != getNStates(ik))
            return false;
      }
   }
   return true;
}

double SxFermi::getBandGap (int *maxValIk, int *minConIk)
{
   int ik, iSpin, iv, nv, ic, nc;
   double maxValEps, minConEps, maxValEpsTl, minConEpsTl;
   int    maxValEpsIk = 0, minConEpsIk = 0;

   for (ik=0; ik < kpPtr->nk; ik++)  {
      for (iSpin=0; iSpin < nSpin; iSpin++)  {
         nv = getNValenceBands (iSpin, ik);
         nc = getNConductionBands (iSpin, ik);
      }
   }

   maxValEpsTl = eps(valBandIdx(0,0,0),0,0);
   minConEpsTl = eps(conBandIdx(0,0,0),0,0);

   for (ik=0; ik < kpPtr->nk; ik++)  {
      for (iSpin=0; iSpin < nSpin; iSpin++)  {
         nv = getNValenceBands (iSpin, ik);
         nc = getNConductionBands (iSpin, ik);

         // --- is it a semiconductor?
         if (nv + nc != eps.getNStates (ik))
            return -1.;

         maxValEps = eps(valBandIdx(0,iSpin,ik),iSpin,ik);
         minConEps = eps(conBandIdx(0,iSpin,ik),iSpin,ik);

         for (iv = 0; iv < nv; iv++)  {
            maxValEps = maximum (maxValEps, eps(valBandIdx(iv,iSpin,ik),iSpin,ik));
            if (maxValEps > maxValEpsTl)  {
               maxValEpsTl = maxValEps;
               maxValEpsIk = ik;
            }
         }
         for (ic = 0; ic < nc; ic++)  {
            minConEps = minimum (minConEps, eps(conBandIdx(ic,iSpin,ik),iSpin,ik));
            if (minConEps < minConEpsTl)  {
               minConEpsTl = minConEps;
               minConEpsIk = ik;
            }
         }
      }
   }

   double gap = fabs(minConEpsTl - maxValEpsTl);

   // --- return type of gap
   if (maxValIk && minConIk)  {
      if (gap > 1e-5)  {
         *maxValIk = maxValEpsIk;
         *minConIk = minConEpsIk;
      }  else  {
         *maxValIk = -1;
         *minConIk = -1;
      }
   }

   return gap;
}


void SxFermi::printOccupation (bool final) const
{
   SX_CHECK (kpPtr);

   int iSpin, ik;
   // int gapValIk, gapConIk;

   double scale  = (printHartree) ?  1. : HA2EV;
   SxString unit = (printHartree) ? "H" : "eV";
   const char *u = unit.ascii();

   // --- sort eigenvalues (for printing)
   for (ik=0; ik < kpPtr->nk; ik++)
      for (iSpin=0; iSpin < nSpin; iSpin++)
         epsSortList(iSpin,ik) = eps(iSpin,ik).getSortIdx();

   if (nSpin == 1) {
      sxprintf ("| Fermi energy:  %12.6f %s\n", eFermi*scale, u);
   } else  {
      if (!keepSpinMoment) {
         sxprintf ("| Fermi energy:  %12.6f %s\n", eFermi*scale, u);
      } else {
         sxprintf ("| Fermi energy spin Up: %12.6f %s, spin Down: %12.6f %s \n",
               eSpinFermi[SPIN_UP]*scale, u, eSpinFermi[SPIN_DOWN]*scale, u);
      }
      sxprintf ("| Spin moment:  %12.6f B\n", spinMoment);
   }
   sxprintf ("| Entropy equilibrium %12.6f , nonequilibrium  %12.6f \n",
            eqEntropy, noneqEntropy);

   for (ik=0; ik < kpPtr->nk; ik++)  {
      for (iSpin=0; iSpin < nSpin; iSpin++)  {
         printK(iSpin, ik, true, final);
         if (nSpin - iSpin - 1 > 0)  sxprintf ("|\n");
         else                        cout << SX_SEPARATOR;
      }
   }

//   // --- print band gap if this is a bandstructure calculation
//   bool keepRho = elecMinimCtrl->keepRho;
//   if (keepRho && gap > 1e-5)  {
//      if ( gapValIk == gapConIk)
//         sxprintf ("| Direct band gap: %g eV   (see marked k-point (*))\n", gap);
//      else
//         sxprintf ("| Indirect band gap: %g eV   (see marked k-points (*))\n", gap);
//      cout << SX_SEPARATOR;
//   }
}

void SxFermi::printSmearing () const
{
   if (orderMethfesselPaxton < 0)
      cout << "Fermi-Dirac";
   else if (orderMethfesselPaxton == 0)
      cout << "Gaussian";
   else
      cout << "Methfessel-Paxton order " << orderMethfesselPaxton;
   cout << endl;
}

void SxFermi::printK(int iSpin, int ik, bool printOcc, bool final) const
{
   double scale  = (printHartree) ?  1. : HA2EV;
   const char *u = (printHartree) ? "H" : "eV";

   SX_CHECK (kpPtr);
   SxVector3<Double> k = kpPtr->kVec(ik);
   bool useLabels = (kpPtr->kLabels.getSize() == kpPtr->nk);

   sxprintf ("| ik=%d", ik+1);
   if (nSpin != 1)  sxprintf (", iSpin=%d", iSpin);
   if (useLabels && kpPtr->kLabels(ik) != "")
      sxprintf (" *%s* ", kpPtr->kLabels(ik).ascii());

   sxprintf (" @ (%g,%g,%g), ", k(0), k(1), k(2));
// if (ik == gapValIk || ik == gapConIk)  sxprintf (" (*) ");
   sxprintf (" w=%g", kpPtr->weights(ik));
   if (final) sxprintf ("\n| final eig [%s]: ", u); else sxprintf ("\n| eig [%s]: ", u);
   for (int i=0; i < focc.getNStates(ik); i++)
      //sxprintf ("%8.4f ", eps(i,iSpin,ik)*scale);
      sxprintf ("%8.4f ", eps(epsSortList(i,iSpin,ik),iSpin,ik)*scale);
   if (printOcc)  {
      if (final) sxprintf ("\n| final focc:     "); else sxprintf ("\n| focc:     ");
      for (int i=0; i < focc.getNStates(ik); i++)
         //sxprintf ("%8.4f ", focc(i,iSpin,ik));
         sxprintf ("%8.4f ", focc(epsSortList(i,iSpin,ik),iSpin,ik));
   }
   sxprintf ("\n");
}




void SxFermi::writeSpectrum (const SxString &filebase,
                             const SxString &fileext) const
{
   FILE *fp = NULL;
   SxString filename;
   int i, iSpin, ik, nStates;
   for (iSpin=0; iSpin < nSpin; iSpin++)  {

      filename = filebase;
      if (nSpin == 2)  filename += "." + SxString(iSpin);
      if (fileext.getSize () > 0) filename += "." + fileext;

      fp = fopen (filename.ascii(), "w");
      if (!fp)  {
         sxprintf ("Can't open filename %s\n", filename.ascii());
         SX_EXIT;
      }
      // --- find maximum number of states
      int nMaxStates = 0;
      for (ik=0; ik < getNk (); ik++)
         nMaxStates = maximum(nMaxStates, focc.getNStates(ik));
      // --- print header line
      fprintf (fp, "# Eigenspectrum: all energies in eV\n");
      fprintf (fp, "# -ik-  |  i = ");
      for (i=0; i < nMaxStates; i++)  fprintf (fp, "%-4d\t", i + 1);
      fprintf (fp, "\n");
      // --- print data
      double fillVal = 0.;
      for (ik=0; ik < getNk (); ik++)  {
         fprintf (fp, "%d\t", ik + 1);
         nStates = focc.getNStates(ik);
         for (i=0; i < nStates; i++)
            fprintf (fp, "%12.6f\t",
                     eps(epsSortList(i,iSpin,ik),iSpin,ik)*HA2EV);
         for (i=nStates; i < nMaxStates; i++)
            fprintf (fp, "%g\t", fillVal);
         fprintf (fp, "\n");
      }
      fclose (fp);
   }
}


// TODO   Check support for parallel NetCDF4 IO.
// TODO   check if it works w/ SxParallelHierarchy
void SxFermi::write (const SxBinIO &io) const
{
   try  {
      int iSpin, ik, nStates;
      int nElem = 0;
      int nk = getNk ();
      SxVector<Int> nPerK (nk);

      for (ik=0; ik < nk; ik++)  {
         if (SxLoopMPI::myWork(ik))
            nPerK(ik) = focc.getNStates(ik);
         else
            nPerK(ik) = 0;
      }
      SxLoopMPI::sum(nPerK);

      // --- get total number of states
      for (ik=0; ik < nk; ik++)  {
         if (SxLoopMPI::myWork(ik))
            nElem += focc.getNStates(ik) * nSpin;
      }
      nElem = SxLoopMPI::sum(nElem);

      int offset = 0;

      // --- write dimensions
      const SxString nElemName = "nAllStates";
      io.addDimension (nElemName, nElem);
      io.addDimension ("nSpin",   nSpin);
      io.addDimension ("nk",      nk);

      // --- write data

//      if (!io.contains ("nPerK") || (io.ncMode == SxBinIO::WRITE_DATA))
      if ((io.ncMode == SxBinIO::WRITE_HEADER) || (SxLoopMPI::me() == 0))
         io.write ("nPerK", nPerK, "nk");

      for (ik=0; ik < nk; ik++)
      {
         bool flag = false;
         nStates = nPerK(ik);   // focc.getNStates (ik);
         for (iSpin=0; iSpin < nSpin; iSpin++)
         {
            if (SxLoopMPI::myWork(ik)) /* do not care about the spin variable */
            {
               io.write ("epsSortIdx", epsSortList(iSpin,ik), nElemName, offset);
               io.write ("focc", focc(iSpin,ik), nElemName, offset);
               io.write ("eps", eps(iSpin,ik), nElemName, offset);
               if (io.ncMode == SxBinIO::WRITE_HEADER)
               {
                  flag = true;
                  break;
               }
            }
            offset += nStates;
         }
         if (flag) break;
      }

      if ((io.ncMode != SxBinIO::WRITE_HEADER) && (SxLoopMPI::nr() == 1))
      {
         // at the end, all data must be written
         SX_CHECK (offset == nElem, offset, nElem);
      }
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}


#include <fstream>
void SxFermi::readSpectrum (const SxString &file,
                            SxCell *cellPtr,
                            SxKPoints *kp,
                            int iSpin)
{
   SX_CHECK (getNSpin () > 0);
   SX_CHECK (getNStates () > 0);
   SX_CHECK (getNk () > 0);
   SX_CHECK (iSpin >= 0 && iSpin < getNSpin (), iSpin, getNSpin ());

   bool check = cellPtr && kp;
   SxCell bMat;
   if (check) bMat = TWO_PI * cellPtr->inv;

   ifstream epsFile;
   epsFile.open (file.ascii ());
   if (!epsFile)  {
      cout << "Eigenvalue file '" << file << "' can't be opened." << endl;
      SX_QUIT;
   }
   //cout << "Reading eps from '" << file << '\'' << endl;

   // discard first two lines if comment ("#")
   if (epsFile.peek () == '#')
      epsFile.ignore (100000 /* large */,'\n');
   if (epsFile.peek () == '#')
      epsFile.ignore (100000 /* large */,'\n');

   int ik, nk = getNk (), id;
   int i,  n = getNStates ();
   for (ik = 0; ik < nk; ik++)  {
      if (!epsFile)  {
         cout << "Unexpected end of file '" << file << "'" << endl;
         epsFile.close ();
         SX_QUIT;
      }
      if (epsFile.peek () == '#')  {
         Coord k1;
         epsFile.get (); // '#'
         epsFile >> k1(0) >> k1(1) >> k1(2);
         epsFile >> ws; // skip whitespace
         if (check)  {
            Coord k2 = bMat.inv ^ kp->getK(ik);
            if ((k1 - k2).normSqr () > 1e-7)  {
               cout << "k-point mismatch for ik = " << (ik+1);
               cout << ".\nExpected " << k2 << ", but found " << k1 << endl;
               epsFile.close ();
               SX_QUIT;
            }
         }
      }
      epsFile >> id;
      if (id != ik + 1)  {
         cout << "Unexpected item in '" << file << "'. Expected " << (ik + 1);
         cout << ", but found something else." << endl;
         char next[100];
         epsFile.getline (next, 100);
         cout << "Next characters are:" << endl << next << endl;
         epsFile.close ();
         SX_QUIT;
      }
      for (i = 0; i < n; i++)  {
         epsFile >> eps(i,iSpin,ik);
         if (!epsFile.good ())  {
            cout << "Failed reading '" << file << "' at ";
            cout << "ik = " << (ik + 1);
            cout << "; i = " << (i + 1);
            cout << endl;
            epsFile.close ();
            SX_QUIT;
         }
      }
      eps(iSpin, ik) /= HA2EV; // reading in eV, internally Hartree
      epsFile >> ws; // skip whitespace
      epsSortList(iSpin,ik) = eps(iSpin,ik).getSortIdx();
   }
   epsFile.close ();
}

void SxFermi::peekSpectrumFile (const SxString &epsFile,
                                int *nkPtr,
                                int *nStatesPtr)
{
   SX_CHECK (nkPtr);
   SX_CHECK (nStatesPtr);
   int &nk = *nkPtr;
   int &nStates = *nStatesPtr;
   SX_CHECK (epsFile.getSize () > 0);

   // --- open file
   ifstream file (epsFile.ascii ());

   if (!file)  {
      cout << "Can't open file '" << epsFile << "'." << endl;
      SX_QUIT;
   }

   int ik;
   nk = 0;

   // --- get nk
   while (file.good ())  {
      if (file.peek () != '#')  {
         file >> ik;
         if (!file.good ()) break;
         nk++;
         if (ik != nk)  {
            cout << "Invalid file format in '" << epsFile;
            cout << "'.\n Expected " << nk;
            cout << ", but found " << ik << '.' << endl;
            SX_QUIT;
         }
      }
      file.ignore (100000,'\n');
   }
   file.close ();
   SX_CHECK (nk > 0, nk);

   // --- get nStates
   file.clear ();
   file.open (epsFile.ascii ());
   nStates = 0;
   double eps;
   while (file.good ())  {
      if (file.peek () != '#')  {
         file >> ik;
         if (ik == nk)  {
            while (file.good ())  {
               file >> eps;
               if (file.good ()) nStates++;
            }
            break;
         }
      }
      file.ignore (100000,'\n');
   }
   file.close ();
   SX_CHECK (nStates > 0, nStates);

}

void SxFermi::read (const SxBinIO &io, bool keepNStates)
{
   int iSpin, ik, nk, nSpin_, nStates = -1;
   try  {
      nk       = io.getDimension ("nk");
      nSpin_   = io.getDimension ("nSpin");
      if (nSpin_ != nSpin && nSpin > 0)  {
         cout << "SxFermi::readFermi - ";
         cout << "Inconsistency: nSpin = " << nSpin_;
         cout << " (should be = " << nSpin << ") in " << io.filename << endl;
         SX_QUIT;
      } else {
         nSpin = nSpin_;
      }

      if (!keepNStates)  {
         focc.bundle.resize (nk);
         eps.bundle.resize (nk);
      }
      SX_CHECK (nk = getNk (), nk, getNk ());

      epsSortList.bundle.resize (nk);
      valBandIdx.bundle.resize (nk);
      conBandIdx.bundle.resize (nk);
      nValBands.resize (nk);
      nConBands.resize (nk);

      SxDiracVec<Int> nPerK(nk);
      int offset = 0;
      io.read ("nPerK", &nPerK, nk);

      int i;
      PrecEps eps_max;
      for (ik=0; ik < nk; ik++)  {
         if (keepNStates)
            nStates = getNStates(ik);
         else
            nStates = nPerK(ik);
         focc.bundle(ik).resize (nSpin);
         eps.bundle(ik).resize (nSpin);
         epsSortList.bundle(ik).resize (nSpin);
         valBandIdx.bundle(ik).resize(nSpin);
         conBandIdx.bundle(ik).resize(nSpin);
         nValBands(ik).resize (nSpin);
         nConBands(ik).resize (nSpin);
         for (iSpin=0; iSpin < nSpin; iSpin++)  {
            focc.bundle(ik)(iSpin).resize (nStates);
            eps.bundle(ik)(iSpin).resize (nStates);
            epsSortList.bundle(ik)(iSpin).resize (nStates);
            valBandIdx.bundle(ik)(iSpin).resize( nStates);
            conBandIdx.bundle(ik)(iSpin).resize( nStates);


            if (nStates <= nPerK(ik))  {
               io.read ("focc", &focc(iSpin,ik), nStates, offset);
               io.read ("eps", &eps(iSpin,ik),  nStates, offset);
            } else {
               io.read ("focc", &focc(iSpin,ik), nPerK(ik), offset);
               io.read ("eps", &eps(iSpin,ik),  nPerK(ik), offset);

               // --- set unread elements
               eps_max = eps(iSpin,ik).maxval ();
               for (i = nPerK(ik); i < nStates; i++)  {
                  focc(i,iSpin,ik) = 0.;
                  eps(i,iSpin,ik) = eps_max + i;
               }
            }

            offset += nPerK(ik);

            // --- sort eigenvalues for printing
            epsSortList(iSpin,ik) = eps(iSpin,ik).getSortIdx();
         }
      }
   }  catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   updateValCon ();

   // now calculate number of electrons
   double nElectrons_ = 0;
   SxKPoints kPoints;
   SxVector<Double> weights;
   if(kpPtr == NULL) {
      try {
         kPoints.read(io);
      }
      catch (SxException e)   {
         e.print ();
         SX_EXIT;
      }
      weights = kPoints.weights;
   } else {
      weights = kpPtr->weights;
   }

   double fMax = 0.;
   for(ik = 0; ik < nk; ik++)   {
      for(iSpin = 0; iSpin < nSpin; iSpin++)   {
         for(int iband = 0; iband < nStates; iband++)   {
            nElectrons_ += weights(ik) * focc(iband,iSpin,ik);
            if (focc(iband,iSpin,ik) > fMax)
               fMax = focc(iband,iSpin,ik);
         }
      }
   }

   if((nElectrons > 1e-15) && fabs(nElectrons - nElectrons_) > nElectrons*1e-8)   {
      cout << "SxFermi::readFermi - ";
      cout << "Inconsistency: nElectrons differ by "
           << fabs(nElectrons - nElectrons_) << " !\n"
           << "nElectrons = " << nElectrons_
           << " (should be = " << nElectrons << ") in " << io.filename << endl;
      SX_QUIT;
   } else {
      nElectrons = nElectrons_;
   }

   if (fabs(fMax - 1.) < 1e-6) fFull = 1.;
   else if (fabs(fMax - 2.) < 1e-6) fFull = 2.;
   else if (fMax > 1. && fMax < 2.) {
      fFull = 2.;
      cout << "Warning: setting fFull = 2. from max. observed f=" << fMax
           << endl;
   } else {
      fFull = (nSpin == 1) ? 2. : 1.;
      cout << "Warning: setting fFull=" << fFull << " from nSpin=" << nSpin
           << endl;
      cout << "max. observed focc= " << fMax << endl;
   }
}

Focc SxFermi::getFoccByWindow (const double lowEnergy, const double highEnergy, const double ekt)
{

   int nStates = getNStates ();
   int nk      = getNk ();
   Focc result (nStates, nSpin, nk);
   result.set(0.0);

   double s = 3. - double(nSpin);

   for (int iState = 0; iState < nStates; iState++)  {
      for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
         for (int ik = 0; ik < nk; ik++)  {
            double energy = eps(iState,iSpin,ik);
            double low,high;
            if (fabs(ekt) > 1e-10)  {
               low  = 1.0 / (exp((lowEnergy - energy) / ekt) + 1.0);
               high = 1.0 / (exp((energy - highEnergy) / ekt) + 1.0);
            } else {
               low  = energy < lowEnergy ? 0.0 : 1.0;
               high = energy > highEnergy ? 0.0 : 1.0;
            }
            result(iState,iSpin,ik) = s * low * high;
         }
      }
   }

   return result;
}

int SxFermi::getHOMO ()
{
   int result = -1;

   int nStates = getNStates ();
   int nk      = getNk ();

   for (int iState = 0; iState < nStates; iState++)  {
      for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
         for (int ik = 0; ik < nk; ik++)  {
            if ((focc(iState,iSpin,ik) > 1e-12) && iState > result)
               result = iState;
         }
      }
   }

   return result;
}
