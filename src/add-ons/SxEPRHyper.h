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

#ifndef _SX_EPR_HYPER_H_
#define _SX_EPR_HYPER_H_

#include <SxRho.h>
#include <SxMuPW.h>
#include <SxExt.h>
#include <SxPartialWaveBasis.h>

/** \brief EPR hyperfine parameters

    \b SxClass = S/PHI/nX EPR hyperfine parameters

    \author Christoph Freysoldt, freysoldt@mpie.de
    \author Gernot Pfanner, pfanner@mpie.de
 */
class SX_EXPORT_EXT SxEPRHyper
{
   public:
      /// spin density in reciprocal space
      PsiG spinDensityG;
      
      /** \brief Compute \f$\langle r^{-3}\rangle \f$ matrix element
           @param rad   Logarithmic radial grid
           @param up    radial part of p-function
           @param logDr numerical d(ln r)=ln(r_{i+1}/r_i)
           @return the matrix element

           The definition is
           \f[
           \langle \phi_p |r^{-3}|\phi_p\rangle =
           \int_0^{\infty} r^2 dr \phi^*_p(r) \frac{1}{r^3} \phi_p(r)
           \f]
           where \f$\phi_p(r)=u_p(r)/r\f$.
           
        */

      SxConstPtr<SxRadBasis> radBasisPtr;

      SxPseudoPot psPot;

      double rcut;

      /// ae/ps ratio for s-density
      SxArray<double> sEnhancement;

      /// container for partial waves differences
      SxArray<double> deltaRhoPSAE;

      /// <r-3> ae-ps difference
      SxArray<double> rm3diff;

      /// gyromagnetic ratio for different species (in MHz/T)
      SxArray<double> gyromagneticRatio;
      
      /// density at r=0 for high cut-off (improves convergency)
      SxArray<double> rhoNucleus;

      // nuclear charge
      SxArray<int> nucCharge;

      /// Container for rad & psi
      class SxRadPsi  {
         public:
            /// Radial mesh
            SxDiracVec<Double> rad,
            /// u-function (= psi*r)
                               u;
      };

      /// Read a function from the fort.39 file from fhipp
      static SxRadPsi readFort39 (FILE *fp);

      /// Read gyromagnetic ratios from input file
      void readHfData (const SxSymbolTable *table, bool gyrRatioTable);

            /** \brief Read all-electron/pseudo waves from fhi98pp fort.39 files

        This reads the fort.39 files and computes
        1) the ae/ps-ratio of s-waves at the origin
        2) the difference in the <r-3> matrix elements of p-waves
        */
      void readAePsData (const SxArray<SxString> &psiFiles, bool nonRelIso);

      /// Compute hf-parameter (single-projector approach)
      void compute (const SxAtomicStructure &structure,  
                    const SxPtr<SxMuPW> &muPtr = SxPtr<SxMuPW>(), 
                    const SxPW &waves=SxPW(),
                    const SxFermi &fermi=SxFermi());

};

// --- cutoff-function
inline SxDiracVec<Double> fCut(SxDiracVec<Double> rad, double rcut)
{
      int i = 0, nr = int(rad.getSize ());
      SxDiracVec<Double> fCut(nr);
      fCut.set(0);

      while (rad(i) <= rcut)  {
         fCut(i) = 1.;
         // fCut(i) = sqr((sin(PI*rad(i)/1.25)/(PI*rad(i)/1.25)));
         i++;
      }

      // fCut(i) = (rad(i) < rcut) ? (exp(-rad(i)*rad(i)/0.125)) : 0.;
      // fCut(i) = (rad(i) < rcut) ? 1. : (54.5982*exp(-rad(i)/0.25));

      return fCut;
}

#endif /* _SX_EPR_HYPER_H_ */
