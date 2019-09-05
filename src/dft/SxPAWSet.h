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

#ifndef _SX_PAW_SET_H_
#define _SX_PAW_SET_H_

#include <SxDFT.h>
#include <SxPWSet.h>
#include <SxPAWBasis.h>
#include <SxNArray.h>

/** \brief Complete set of PAW wavefunctions


    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DFT SxPAWSet : public SxPWSet
{
   public:
      /** \brief Constructor
        @param gk      |G+k> basis
        @param pBasis  partial wave basis
        @param n       number of states
        @param nSpin   number of spins
        @param tmpDir  if non-zero, keep waves on disk in this directory
        */
      SxPAWSet (SxPtr<SxGkBasis> gkPtrIn,
                const SxPtr<SxPartialWaveBasis> &pBasis,
                int n, int nSpin, const SxString &tmpDir = "");
      /** \brief Constructor
        @param gk          |G+k> basis
        @param pawBasisIn  PAW basis
        @param n           number of states
        @param nSpin       number of spins
        @param tmpDir      if non-zero, keep waves on disk in this directory
        */
      SxPAWSet (SxPtr<SxGkBasis> gkPtrIn,
                const SxArray<SxPtr<SxPAWBasis> > pawBasisIn,
                int n, int nSpin, const SxString &tmpDir = "");
      
      SxPAWSet (SxConstPtr<SxPAWPot> potPtr, 
                const SxAtomicStructure &structure, 
                const SxBinIO &io);

      /// Copy constructor (not implemented)
      SxPAWSet (const SxPAWSet &in) : SxPWSet (in) { SX_EXIT; }

      /// Assignment operator (not implemented)
      void operator= (const SxPAWSet &) { SX_EXIT; }

      /// Destructor
      virtual ~SxPAWSet ();
   protected:
      /// Number of states
      int nStates;

      /// Number of spin channels
      int nSpin;

      /// Wave functions
      mutable SxArray2<PsiGI> waves;

      /// Basis sets (:ik)
      SxArray<SxPtr<SxPAWBasis> > pawBasis;

      //@{ Waves on disk
      /// True if waves are kept on disk
      bool keepWavesOnDisk;

      /// loaded k-point
      mutable int loadedK;

      /// loaded spin channel
      mutable int loadedSpin;

      /// File containing wave function coefficients
      SxBinIO wavesFile;

      /// Create scratch file
      void createScratchFile (const SxString &tmpDir);

   public:
      /// Load waves
      void loadWaves (int iSpin, int ik) const;

      /// Flush waves
      void flushWaves () const;

      /// Whether waves are kept on disk
      virtual bool wavesOnDisk () const { return keepWavesOnDisk; }

      //@}
   public:
      /// Get the basis
      virtual const SxBasis& getBasis (int ik) const
      {
         SX_CHECK (ik >= 0 && ik < getNk (), ik, getNk ());
         SX_CHECK (pawBasis.getSize () == getNk (),
                   pawBasis.getSize (), getNk ());
         return *pawBasis(ik);
      }


      /** \brief Returns the number of available states. */
      virtual int getNStates (int =0) const
      { 
         return nStates; 
      }
      /** \brief Returns the number of available spin channels. */
      virtual int getNSpin () const {
         return nSpin;
      }
      /** \brief Returns the number of \b k points */
      virtual int getNk () const
      { 
         return int(pawBasis.getSize ());
      }
      /** \brief Returns the pBasis */
      SxConstPtr<SxPartialWaveBasis> getPBasis () const
      {
         return pawBasis(0)->pBasis;
      }
      /** \brief Set the GkBasis 
          The implementation must set gkBasisPtr and loop over all its
          array-indices and call the setBasis routine.
       */
      virtual void setGkBasisPtr (SxPtr<SxGkBasis>);

      /** \brief returns a single Bloch-like state
          This function will be used from the SxPWHamiltonian to get
          an Bloch-like state i, with the spin quantum number \f$\sigma\f$, 
          and at the \b k point k

          \param   i      state or band index
          \param   iSpin  spin quantum number \f$\sigma\f$
          \param   k      index of the \b k point

          \return  state \f$\langle\bf{G+k}|\Psi_{i,\sigma,\bf{k}}\rangle\f$
       */
      virtual const SxGBasis::TPsi operator() (int i, int iSpin, int ik) const;

      /** \brief returns a single Bloch-like state
          This function will be used from the SxPWHamiltonian to get
          an Bloch-like state i, with the spin quantum number \f$\sigma\f$, 
          and at the \b k point k

          \param   i      state or band index
          \param   iSpin  spin quantum number \f$\sigma\f$
          \param   k      index of the \b k point

          \return  state \f$\langle\bf{G+k}|\Psi_{i,\sigma,\bf{k}}\rangle\f$
       */
      virtual SxGBasis::TPsi operator() (int i, int iSpin, int ik);

      /** \brief returns Bloch-like states
          This function will be used from the SxPWHamiltonian to get
          the Bloch-like states with the spin quantum number \f$\sigma\f$, 
          and at the \b k point k

          \param   iSpin  spin quantum number \f$\sigma\f$
          \param   k      index of the \b k point

          \return  state \f$\langle\bf{G+k}|\Psi_{i,\sigma,\bf{k}}\rangle\f$
       */
      virtual const SxGBasis::TPsi& operator() (int iSpin, int ik) const;
      /** \brief returns Bloch-like states
          This function will be used from the SxPWHamiltonian to get
          the Bloch-like states with the spin quantum number \f$\sigma\f$, 
          and at the \b k point k

          \param   iSpin  spin quantum number \f$\sigma\f$
          \param   k      index of the \b k point

          \return  state \f$\langle\bf{G+k}|\Psi_{i,\sigma,\bf{k}}\rangle\f$
       */
      virtual SxGBasis::TPsi& operator() (int iSpin, int ik);

      using SxPWSet::operator();

      /** \brief Derefences the wavefunction as a matrix {ng x nBlock}
        \param iStart  starting index for block
        \param nBlock  number of states in the block
        \param iSpin   spin index
        \param ik      k-index
      */
      virtual PsiG getBlock (int iStart, int nBlock, int iSpin, int ik);

      /// Minimize memory usage (if possible)
      virtual void memMinimize ();

      /// Create new SxPAWSet with same size as current one
      virtual SxPtr<SxPWSet> getNew () const;

      virtual void read (const SxBinIO &io, int mode = SxPWSet::KeepGkBasis);
      virtual void write (SxBinIO &io) const;

      void readPAWBasis (const SxBinIO &io, SxConstPtr<SxPAWPot> potPtr, const SxAtomicStructure &structure);
};

#endif /* _SX_PAW_SET_H_ */
