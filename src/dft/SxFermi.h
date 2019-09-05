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

#ifndef _SX_FERMI_H_
#define _SX_FERMI_H_

#include <SxTypes.h>
#include <SxKPoints.h>
#include <SxBinIO.h>
#include <SxDFT.h>

/** brief Fermi distribution function

  \b SxFermi = S/PHI/nX Fermi distribution

  In order to treat systems ate finite temperatures a Fermi like 
  function has to be
  
  \ingroup group_dft
  \author  Sixten Boeck
  */
class SX_EXPORT_DFT SxFermi
{
   public:
      
      bool    printHartree;
      Focc    focc;
      Eps     eps;
      mutable EpsSort epsSortList;
      /// Total number of electrons in the distribution
      Real8   nElectrons;
      /// Number of spin channels
      int     nSpin;
      /// Spin prefactor (2 if non-polarized, 1 otherwise)
      Real8   fFull;
      /// Fermi energy
      Real8   eFermi;
      /// Fermi energies for spin up/down channels if spin moment is fixed
      Real8   eSpinFermi[2]; //:iSpin
      /// \f$\beta = 1/kT\f$
      Real8   beta;
      /// Spin moment
      Real8   spinMoment;
      Real8   eqEntropy;
      Real8   noneqEntropy;
      bool    keepSpinMoment;
      /// Order of MethfesselPaxton scheme, -1 for Fermi-Dirac
      int     orderMethfesselPaxton;

      /// Whether to use occupation numbers or not
      enum FoccHandling { NoFocc, UseFocc};

      const SxKPoints *kpPtr;

      
      SxFermi ();
      SxFermi (Real8 nElectrons, int nStates, int nSpin, const SxKPoints &);
      SxFermi (int nStates, int nSpin, int nk);
      SxFermi (const SxBinIO &io, bool keepNStates = false) { 
         kpPtr = NULL; 
         nSpin = -1;
         nElectrons = -1;
         read (io, keepNStates);
         orderMethfesselPaxton = -1;
      }
      ~SxFermi ();

      void setOccupencies (const Focc &initialFocc);

   protected:
      /** \brief Get index-list from symbol table.
          Parses range[from, to] and values=[...]
          Negative indices are interpreted as from-end-counting.

          \note Auxiliary routine for readOccupations
        */
      static SxList<int> getIdxList(const SxSymbolTable *table,
                                     int nMax,
                                     const SxString &name);
   public:
      /** \brief Read occupation numbers from symbol table
          \note: the size must be initialized
        */
      void readOccupations (const SxSymbolTable *);
      void setSpinMoment (Real8 spinMoment, bool keepFixed=false);
      void fermiDistribution (Real8 ekt, Real8 mixFactor=1.0);
      void printOccupation (bool final=false) const;
      /** \brief Print one k-point
           @param iSpin spin index
           @param ik    k-point index
           @param printOcc whether to print the occupation numbers, too
        */
      void printK (int iSpin, int ik, bool printOcc, bool final=false) const;
      /// print the type of smearing, with newline
      void printSmearing () const;

      SxBundle3<Int> valBandIdx;
      SxBundle3<Int> conBandIdx;
      SxList<SxList<int> > nValBands, nConBands;  // :ik,:iSpin

      /// Get number of states 
      int getNStates (int ik = 0) const;
      /// Get number of k-points
      int getNk      () const;
      /// Get number of spin channels
      inline int getNSpin () const { return nSpin; }

      int getNValenceBands     (int iSpin, int ik) const;
      /// Get max. number of valence bands for all k and spin
      int getNValenceBands     () const;
      /// Get max. number of conduction bands for all k and spin
      int getNConductionBands  () const;
      int getNConductionBands  (int iSpin, int ik) const;
      int getValenceBandIdx    (int i, int iSpin, int ik) const;
      int getConductionBandIdx (int i, int iSpin, int ik) const;
      /**
        \brief Is this a semiconductor?

        The criteria for a semiconductor we use here is that the number
        of occupied states is the same for all k-points. It is furthermore
        checked that the number of occupied and unoccupied states is the
        total number of states.
        It is crucial that updateValCon () has been called before this 
        routine.
        */
      bool isSemiconductor ();
      double getBandGap (int *maxValIk=NULL, int *minConIk=NULL);
      double getEBand (enum FoccHandling) const;

      /** \brief Update number of occupied and unoccupied bands
        \param fOccMin   minimum occupation to regard band as occupied
        \param fUnoccMax maximum occupation to regard band as unoccupied
        \note At finite temperature, a band may be 
              considered occupied and unoccupied at the same time.
        \note For nSpin == 2, the fUnoccMax factor is multiplied by 2!
      */
      void updateValCon (double fOccMin   =    1e-7, 
                         double fUnoccMax = 1.-1e-7);
      void checkOccupation ();
      /// Write spectrum to ascii files
      void writeSpectrum (const SxString &filebase, 
                          const SxString &fileext) const;
      /// Read spectrum from ascii file
      void readSpectrum (const SxString &file,
                         SxCell *cellPtr = NULL,
                         SxKPoints *kp = NULL,
                         int iSpin = 0);
      /** \brief Try to extract nk and nStates from a spectrum file (eps.dat)
        \param epsFile    filename
        \param nkPtr      pointer to nk
        \param nStatesPtr pointer to nStates
        */
      static void peekSpectrumFile (const SxString &epsFile,
                                    int *nkPtr, int *nStatesPtr);
                         
      void write (const SxBinIO &) const;
      void read (const SxBinIO &, bool keepNStates = false);
      Real8 getSpinMoment ();
      Real8 getEntropy ();

      /// Methfessel-Paxton occupation number
      static Real8 foccMethfesselPaxton(Real8 x, int order);
      /// Methfessel-Paxton distribution function
      static Real8 dfoccMethfesselPaxton(Real8 x, int order, double *lastDelta=NULL);

      /** \brief Derivative of occupation number with respect to eFermi
          @param i     state index
          @param iSpin spin index
          @param ik    k-point index
          @return \f$\frac{d f}{d E^{\rm Fermi}}\f$ (including spin prefactor)

          @note The spin prefactor is 2 for spin-nonpolarized calculations,
                otherwise 1.
          @note: trivially,
          \f$ \frac{df}{d\epsilon} = -\frac{df}{d E^{\rm Fermi}}\f$

          For the Fermi-Dirac distribution (excluding spin prefactor)
          \f[
             \frac{d f}{d E^{\rm Fermi}} = \frac{1}{k_B T} f (1-f) \;.
          \f]

          For the Methfessel-Paxton distribution (excluding spin prefactor)
          \f[
             \frac{d f}{d E^{\rm Fermi}} =
               \frac{1}{2}\left\[ D_N(x) - D_{N-1}(x)\right\]
          \f]
          with \f$ x = \frac{\epsilon - E^{\rm Fermi}}{kT} \f$.
          \f$D_N\f$ denotes the Methfessel-Paxton approximation to the
          Delta function
          \f[
              D_N(x) = \sum_n A_n H_{2n}(x) e^{-x^2}
          \f]
          with Hermite polynomials H and
          \f[
               A_n = \frac{1}{n! 4^n \sqrt\pi}
          \f]
        */
      Real8 dFoccFermi(int i, int iSpin, int ik) const;


   protected:
      /// Compute the excess number of electrons in a spin channel
      Real8 fermiFunction (Real8 energy, Real8 nEl, int spin);
      /// Compute the Fermi energy (possibly restricted to 1 spin channel)
      Real8 getFermiEnergy (Real8 ekt, Real8 nEl, int spin=0); 
   public:

      /// Resize number of states (keep spin and k-points as is)
      void resize (int nStates);
      /// Resize number of states (variable per k-point)
      void resize (const SxArray<int> &nPerK);

      Focc getFoccByWindow (const double lowEnergy, const double highEnergy, const double ekt);

      int getHOMO ();

   protected:

      void init ();
};


#endif /* _SX_FERMI_H_ */
