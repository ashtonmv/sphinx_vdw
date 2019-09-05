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

#ifndef _SX_PW_PSI_SET_H_
#define _SX_PW_PSI_SET_H_

#include <SxPsiSet.h>
#include <SxGBasis.h>
#include <SxGkBasis.h>
#include <SxDFT.h>
#include <SxLoopMPI.h>

class SxFermi;

/**
  \brief Abstract class for all kind of derived plane-wave classes.

  This abstract class defines virtual functions called from the corresponding
  Hamiltonian (SxPWHamiltonian). Thus, deriving from SxPWSet allows
  using any kind of basis-set to be used with a plane-wave Hamitlonian.

  \sa      \ref page_dirac
  \sa      SxPWHamiltonian
  \ingroup group_dirac
  \author  Sixten Boeck
  */
class SX_EXPORT_DFT SxPWSet
: public SxPsiSet,
  public SxThis<SxPWSet>
{
   public:
      SxPWSet ();
      SxPWSet (SxPtr<SxGkBasis> inPtr); 

      virtual ~SxPWSet () { }

      /** \brief Returns the number of available states. */
      virtual int getNStates (int =-1) const { return -1; }
      /** \brief Returns the number of available spin channels. */
      virtual int getNSpin () const         { return -1; }
      /** \brief Returns the number of \b k points */
      virtual int getNk () const            { return -1; }
      /** \brief Returns a reference to \f$|\bf{G+k}|\rangle\f$. */
      SxGkBasis &getGkBasis () const;
      /** \brief Returns the SxPtr to \f$|\bf{G+k}|\rangle\f$. */
      SxPtr<SxGkBasis> getGkBasisPtr () const
      {
         return gkBasisPtr;
      };
      /** \brief Set the GkBasis 
          The implementation must set gkBasisPtr and loop over all its
          array-indices and call the setBasis routine.
       */
      virtual void setGkBasisPtr (SxPtr<SxGkBasis> ) = 0;

      /// Get the basis
      virtual const SxBasis& getBasis (int ik) const
      {
         SX_CHECK (gkBasisPtr);
         SX_CHECK (ik >= 0 && ik < gkBasisPtr->getNk (),
                   ik, gkBasisPtr->getNk ());
         return getGkBasis ()(ik);
      }

      /** \brief returns a single Bloch-like state
          This function will be used from the SxPWHamiltonian to get
          an Bloch-like state i, with the spin quantum number \f$\sigma\f$, 
          and at the \b k point k

          \param   i      state or band index
          \param   iSpin  spin quantum number \f$\sigma\f$
          \param   k      index of the \b k point

          \return  state \f$\langle\bf{G+k}|\Psi_{i,\sigma,\bf{k}}\rangle\f$
       */
      virtual const SxGBasis::TPsi operator() (int i, int iSpin, int ik) const
         =0;

      /** \brief returns a single Bloch-like state
          This function will be used from the SxPWHamiltonian to get
          an Bloch-like state i, with the spin quantum number \f$\sigma\f$, 
          and at the \b k point k

          \param   i      state or band index
          \param   iSpin  spin quantum number \f$\sigma\f$
          \param   k      index of the \b k point

          \return  state \f$\langle\bf{G+k}|\Psi_{i,\sigma,\bf{k}}\rangle\f$
       */
      virtual SxGBasis::TPsi operator() (int i, int iSpin, int ik)
         =0;

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

      /** \brief Derefences the wavefunction as a matrix {ng x nBlock}
        \param iStart  starting index for block
        \param nBlock  number of states in the block
        \param iSpin   spin index
        \param ik      k-index
      */
      virtual PsiG getBlock (int iStart, int nBlock, int iSpin, int ik);

      /// Returns a single Bloch state (autoloop)
      template<class T1, class T2, class T3>
      const SxGBasis::TPsi operator() (const T1 &i,
                                       const T2 &iSpin,
                                       const T3 &ik) const
      {
         SxAutoLoop::setLimit (i, getNStates ((int)ik));
         SxAutoLoop::setLimit (iSpin, getNSpin ());
         SxAutoLoop::setLimit (ik, getNk ());
         return operator() ((int)i, (int)iSpin, (int)ik);
      }
      
      /// Returns a single Bloch state (autoloop)
      template<class T1, class T2, class T3>
      SxGBasis::TPsi operator() (const T1 &i,
                                       const T2 &iSpin,
                                       const T3 &ik)
      {
         SxAutoLoop::setLimit (i, getNStates ((int)ik));
         SxAutoLoop::setLimit (iSpin, getNSpin ());
         SxAutoLoop::setLimit (ik, getNk ());
         return operator() ((int)i, (int)iSpin, (int)ik);
      }

      /// returns Bloch-like states (autoloop)
      const SxGBasis::TPsi& operator() (const SxAutoLoop &iSpin, 
                                        const SxAutoLoop &ik) const
      {
         iSpin.setLimit (getNSpin ());
         ik.setLimit (getNk ());
         return operator() ((int)iSpin, (int)ik);
      }
      /// returns Bloch-like states (autoloop)
      SxGBasis::TPsi& operator() (const SxAutoLoop &iSpin, 
                                  const SxAutoLoop &ik)
      {
         iSpin.setLimit (getNSpin ());
         ik.setLimit (getNk ());
         return operator() ((int)iSpin, (int)ik);
      }

      /// Minimize memory usage (if possible)
      virtual void memMinimize ()
      {
         // empty
      }

      /** \brief Whether waves are kept on disk
        @note This information is relevant when creating new wave
        containers (e.g. in SxHamSolver::allStateCG), which should
        behave the same way.
        */
      virtual bool wavesOnDisk () const { return false; }

      /// Create new SxPWSet with same size as current one
      virtual SxPtr<SxPWSet> getNew () const;

      /// What to read from sxb file
      enum ReadFlag {
         /// Default new nStates and read Gk
         Default        = 0x000,
         /// keep nStates
         KeepNStates   =  0x001,
         /// keep GkBasis
         KeepGkBasis   =  0x010,
         /// keep All 
         KeepAll       =  0x011,
         /// save Memory in G+k basis
         SaveMemory    =  0x100
      };


   protected:
      /** \brief Pointer to the \f$ | \bf{G+k} \rangle \f$ basis. */
      SxPtr<SxGkBasis> gkBasisPtr;

   public:
      virtual void read (const SxBinIO &, int = Default) {SX_EXIT;}
      /** \brief Write a complete waves file
          This function writes a complete waves file. It should keep track
          of major format changes (first digit).
          Minor format changes can be implemented in the corresponding
          write functions.
          \param filename          the filename, usually "waves.sxb"
          \param fermi             Fermi object
          \param structure         structure object
          \param chemName          chemical symbols for all species
          \param Gk                G+k basis
          \param writeParserBuffer if true, write SxParser_buffer as "input"
          \param writeFunc         optional callback function for writing
                                   the waves, otherwise SxPW::write is used.
                                   In order to use this, derive a class from
                                   SxPW, write a new write function and give it
                                   here.
        */
      void writeWavesFile (const SxString          &filename,
                           const SxFermi           &fermi,
                           const SxAtomicStructure &structure,
                           bool                    writeParserBuffer) const;

      /** Write wave functions to io file */
      virtual void write (SxBinIO &io) const;
};

#endif /* _SX_PSI_SET_H_ */
