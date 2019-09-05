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

#ifndef _SX_SPECTRUM_H_
#define _SX_SPECTRUM_H_

#include <SxVector.h>
#include <SxDFT.h>

/** \brief Gaussian broadened spectra

    \b SxClass = S/PHI/nX ...

    This class provides utilities to create a spectrum from a set
    of peaks.

    \author Christoph Freysoldt */
class SX_EXPORT_DFT SxSpectrum
{
   public:

      /// Number of points per Gauss peak
      int nPerPeak;

      /** The spectra, after computation : (iE, iSpec)
        This contains the spectra.
        */
      SxVector<Double> spectra;
      /// Energy scale
      SxVector<Double> energies;

      /// Empty constructor
      SxSpectrum ();

      /// Destructor
      ~SxSpectrum () { /* empty */} 

      /** \brief Sets the target diffusion constant. 
        \param D target diffusion constant, 0 < D << 1

        \note There is a reasonable default, so this is mainly
              for debugging purposes. It must be called before
              the #set routine
      */ 
      void setDiffusion (double D);

      /** Set up the energy range and create new spectra.

        \note Only this routine considers nPerPeak and diffusion, so any
              changes afterwards have no effect.
        */
      void set(double eMin, double eMax, double broadening, int nSpectra = -1);

      /// Set number of spectra
      void resize (int nSpectra);

      /// Get number of spectra
      int getNSpectra () const { 
         return computed ? (int)spectra.nCols () : (int)spectra.nRows (); 
      }

      /** Add a peak */
      void addPeak (double ePeak, double weight);

      /** Add a peak to one of the spectrums */
      inline void addPeak (double ePeak, double weight, int iSpectrum);

      /** Add a peak to each spectrum (with different weights) */
      void addPeak (double ePeak, const SxVector<Double> &weights);

      /** Compute spectra from peaks by diffusion.
          The addPeak function adds only delta peaks to the spectrum. In
          this function, these are broadened by a simulated diffusion
          process, which gives a Gaussian broadening.
        */
      void compute ();

      /// Print spectra to file
      void fprint (FILE *fp);

      /// Get number of energy points
      int getNPoints () const { return nPoints; }

   private:
      /// Compute must be called only once
      bool computed;
      /// Lower energy boundary
      double enMin;
      /// Energgy spacing
      double dE;
      
      /// Gaussian broadening
      double broad;

      /** \brief Diffusion constant 
          The Gauss peaks are computed by letting the delta peaks diffuse.
          This constant must be lower than 1, it is set to a reasonable 
          default 0.1 .
        */
      double diffusion;

      /// Number of points in spectrum
      int nPoints;

      /// Number of time steps
      int nTime;


};

void SxSpectrum::addPeak (double ePeak, double weight, int iSpectrum)
{
   SX_CHECK (! computed);
   SX_CHECK (iSpectrum >= 0 && iSpectrum < getNSpectra (),
             iSpectrum, getNSpectra ());
   int pos = int(ceil((ePeak - enMin)/dE));
   if (pos < 1 || pos >= nPoints) return; // outside energy scale
   // distribute over closest energy grid points
   double dW = (energies(pos) - ePeak) / dE;
#ifdef USE_OPENMP
#pragma omp atomic update
#endif
   spectra(iSpectrum, pos)   += weight * (1. - dW) / dE;
#ifdef USE_OPENMP
#pragma omp atomic update
#endif
   spectra(iSpectrum, pos-1) += weight *       dW  / dE;
}

#endif /* _SX__H_ */
