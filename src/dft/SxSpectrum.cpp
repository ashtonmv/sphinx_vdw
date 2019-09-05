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

#include <SxSpectrum.h>

SxSpectrum::SxSpectrum ()
{
   nPerPeak = 30;
   computed = false;
   enMin = dE = 0.;
   broad = 0.;
   nPoints = nTime = -1;
   diffusion = 0.1;
}

void SxSpectrum::setDiffusion (double D)  {
   SX_CHECK (D > 0. && D < 1., D);
   diffusion = D;
   resize (0);
}

void SxSpectrum::set(double eMin, double eMax, double broadening, int nSpectra)
{
   broad = broadening;
   SX_CHECK (broad > 0., broad);
   // set up internal constants
   nPoints = int(nPerPeak * (eMax - eMin) / broad);
   if (nPoints < 10) nPoints = 10;
   dE = (eMax-eMin)/(nPoints-1);
   double fac = .25 * (broad*broad)/(dE*dE);
   nTime = int(fac / diffusion);
   diffusion = fac / double(nTime);

   // --- set up energy scale
   enMin = eMin;
   energies.resize (nPoints);
   SxVector<Double>::Iterator eIt = energies.begin ();
   for (int i = 0; i < nPoints; ++i) *eIt++ = eMin + i * dE;

   // resize 
   if (nSpectra >= 0) 
      resize (nSpectra);
   else 
      resize (getNSpectra () == 0 ? 1 : getNSpectra ());
}

void SxSpectrum::resize (int nSpectra)
{
   SX_CHECK (nSpectra >= 0);
   spectra.reformat (nSpectra, nPoints);
   spectra.set (0.);
   computed = false;
}

void SxSpectrum::addPeak (double ePeak, double weight)
{
   SX_CHECK (! computed);
   SX_CHECK (getNSpectra () == 1);
   int pos = int(ceil((ePeak - enMin)/dE));
   if (pos < 1 || pos >= nPoints) return; // outside energy scale
   // distribute over closest energy grid points
   double dW = (energies(pos) - ePeak) / dE;
   spectra(pos)   += weight * (1. - dW) / dE;
   spectra(pos-1) += weight *       dW  / dE;
}

void SxSpectrum::addPeak (double ePeak, const SxVector<Double> &weights)
{
   SX_CHECK (! computed);
   SX_CHECK (getNSpectra () == weights.getSize (),
             getNSpectra () , weights.getSize ());
   int pos = int(ceil((ePeak - enMin)/dE));
   if (pos < 1 || pos >= nPoints) return; // outside energy scale
   
   // distribute over closest energy grid points
   double dW = (energies(pos) - ePeak) / dE;
   spectra.colRef(pos  ).plus_assign_ax ((1. - dW) / dE, weights);
   spectra.colRef(pos-1).plus_assign_ax (      dW  / dE, weights);
}

void SxSpectrum::compute ()
{
   SX_CHECK (! computed);

   int nSpectra = getNSpectra ();
   spectra = spectra.transpose ();

   // Number of points to include outside the visible spectrum
   int nBorder = int (.25 * (broad*broad)/(dE*dE) );

   int iSpec, i, it;
   SxVector<Double>::Iterator dsIt, sIt, specIt;
   int size = nPoints + 2 * nBorder;
   SxVector<Double> s(size), ds(size);
   double last,current;
   double bottom, top;
   double &D = diffusion;
   
   for (iSpec = 0; iSpec < nSpectra; ++iSpec)  {
      // copy spectrum into work array s
      sIt = s.begin ();
      specIt = spectra.colRef (iSpec).begin ();
      for (i = 0; i < nBorder; ++i) *sIt++ = 0.;
      for (i = 0; i < nPoints; ++i) *sIt++ = *specIt++;
      for (i = 0; i < nBorder; ++i) *sIt++ = 0.;

      bottom = top = 0.;
      // --- diffusion process
      //   s(i,it+1) = s(i,it) + D * [s(i+1,it) + s(i-1,it) - 2 * s(i,it)]
      //             =    s    + D * ds                    
      for (it = 0; it < nTime; ++it)  {

         // get ds
         dsIt = ds.begin ();
         sIt = s.begin ();
         last = *sIt++;
         current = *sIt++;
         // ds(0)   = s(1)  -s(0);
         *dsIt++ = last + bottom - 2.*current;
         bottom += D*(current-bottom);
         for (i = 1; i < size-1; ++i, ++dsIt, ++sIt)  {
            // ds(i) = s(i+1) + s(i-1) - 2.*s(i);
            *dsIt = /* next */ *sIt + last - 2. * current;
            last = current;
            current = *sIt;
         }
         // ds(size-1) = s(size-2)-s(size-1);
         *dsIt = last + top - 2. * current;
         top += D*(current-top);

         // s += D * ds;
         s.plus_assign_ax (D,ds);
      }

      // copy spectrum back
      spectra.colRef(iSpec) <<= s(SxIdx(nBorder,nBorder+nPoints-1));
   }
   computed = true;
   
}

void SxSpectrum::fprint (FILE *fp)
{
   int iE, iSpec;
   int nSpectra = getNSpectra ();
   for (iE = 0; iE < nPoints; ++iE)  {
      sxfprintf(fp, "%f", energies(iE));
      for (iSpec = 0; iSpec < nSpectra; ++iSpec)
         sxfprintf(fp, "\t%f", spectra(iE,iSpec));
      sxfprintf(fp, "\n");
   }
}


