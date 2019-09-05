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

#ifndef _SX_PW_SET_H_
#define _SX_PW_SET_H_


#include <stdio.h>
#include <SxError.h>
#include <SxTypes.h>
#include <SxList.h>
#include <SxDirac.h>
#include <SxPrecision.h>
#include <SxGkBasis.h>
#include <SxString.h>
#include <SxPWSet.h>
#include <SxBinIO.h>
#include <SxTimer.h>
#include <SxOverlap.h>
#include <SxDFT.h>


/** \brief Container of plane-waves \f$ 
              \langle \bf{G+k} | \Psi_{i,\sigma,\bf{k}} \rangle
           \f$

  SxPW = S/PHI/nX Plane Waves
           

  \ingroup DFT
  \ingroup group_dirac
  \ingroup group_dft
  \sa      \ref page_dirac
  \author  Sixten Boeck 
 */
class SxFermi;
class SX_EXPORT_DFT SxPW 
: public SxPWSet
{
   public:

      /** The ways how the wavefunctions are stored.
        */
      enum StorageModel  {
         /// All wavefunctions are stored in memory
         InMemory,
         /// Wavefunctions are stored in a scratchfile
         KeepOnDisk,
         /** \brief Wavefunctions are read from a binary file
             \note The wavefunctions MUST NOT BE modified. All
             modifications are ignored.
           */
         ReadOnDemand,
         /** \brief Wavefunctions are read from a binary file
             \note The wavefunctions MUST NOT BE modified. All
             modifications are ignored.
           */
         ReadOneByOne,
         /// Uninitialized
         Unknown
      };

      mutable SxArray<SxArray<SxGkBasis::TPsiGI> >  waves; // :ik,:iSpin,:{ig x i}

      SxPW ();
      /** \brief Create properly sized SxPW-object
          @param nStates number of states
          @param nSpin   number of spin channels
          @param tmpDir  if non-empty, the waves are stored on disk
                         in directory tmpDir
          \todo  tmpDir  should be automatic, and StorageModel=KeepOnDisk
                         should be input parameter
       */
      SxPW (int nStates, int nSpin, SxPtr<SxGkBasis> ,
                const SxString &tmpDir="");
      SxPW (const SxPW &);
      /** \brief Create properly sized SxPW-object
          @param nStates number of states
          @param nSpin   number of spin channels
          @param tmpDir  if non-empty, the waves are stored on disk
                         in directory tmpDir
          \todo  tmpDir  should be automatic, and StorageModel=KeepOnDisk
                         should be input parameter
       */
      SxPW (int nSpin, const SxList<int> &ngPerK,
                const SxString &tmpDir="");
      /** Initialize waves from binary file.
        \param filename the binaries file name
        \param how   determines how to initialize
                     - InMemory      Just read the waves
                     - KeepOnDisk    not supported
                     - ReadOnDemand  register file, then read it on demand
                     - Unknown       a simple way to crash
        \sa StorageModel
        
        This is the standard constructor for add-ons that only read the
        waves, but never modify them.
        \Example
\code
SxString inFile = cli.option ("-w|--waves","waves file", "the waves to read")
                  .toString ("waves.sxb");
...
cli.finalize ();
...
SxPW waves (inFile, SxPW::ReadOnDemand);

// how to reduce the number of states (fewer than in file, of course)
int differentNStates;
waves.changeNStates (differentNStates);

// if you need a G+k basis, read it from the file!
SxGkBasis gk;
try  {
  SxBinIO io(wavesFile,SxBinIO::BINARY_READ_ONLY);
  gk.read (io);
} catch (SxException e)  {
  e.print ();
  SX_EXIT;
}
waves.setGkBasis (gk);
\endcode
        */
      SxPW (const SxString filename, enum StorageModel how, const SxString &tmpDir="");
      virtual ~SxPW ();

      /** \brief Transfer operator

          \Note This works only as postponed constructor call, i.e.
          \code
waves = SxPW(....);
          \endcode
        */
      void operator= (const SxPW &);

      /** \brief Set the GkBasis */
      virtual void setGkBasisPtr (SxPtr<SxGkBasis> );

      /** \brief Resizes waves with new number of states

          This function reallocates the wavefunction coefficient arrays
          according to the provided number of states specified for each
          k-point.
          \param nPerK  number of states per k-point
          \sa SxHamSolver::setNStates
          \sa SxFermi::setNStates
          \sa SxPWHamiltonian::setNStates */
      void setNStates (const SxArray<int> &nPerK);


      //@{ Derefences the wavefunction as a matrix {ng x i}
      virtual SxGkBasis::TPsi& operator() (int iSpin, int ik);
      virtual const SxGkBasis::TPsi& operator() (int iSpin, int ik) const;
      //@}

      //@{ Derefences the wavefunction as a matrix {ng x (idx.star...idx.end) }
      SxGkBasis::TPsi operator() (const SxIdx &idx, int iSpin, int ik);
      const SxGkBasis::TPsi
      operator() (const SxIdx &idx, int iSpin, int ik) const;
      virtual PsiG getBlock (int iStart, int nBlock, int iSpin, int ik)  {
         return operator() (SxIdx(iStart, iStart+nBlock-1), iSpin, ik);
      }
      //@}

      /// Create new SxPW with same size as current one
      virtual SxPtr<SxPWSet> getNew () const;

      //@{ Extracts one state from the set of wavefunctions
      virtual SxGBasis::TPsi operator() (int i, int iSpin, int ik);
      virtual const SxGBasis::TPsi operator() (int i, int iSpin, int ik) const;
      //@}

      inline int getNk () const
      {
         if (stored == InMemory)  return int(waves.getSize());
         else                     return nkPoints;
      }


      inline int getNSpin () const
      {
         return nSpin; // works for all storage
      }


      virtual int getNStates (int ik=-1) const;

      void changeNStates (int);


      /** \brief Initialize wavefunction coefficients with random numbers

          This function presets the wavefunction coefficients with
          randomized numbers. The randomized wavefunctions are not
          orthonormalized. If required call ::orthonormalize afterwards.

          \par Example:
\code
   SxGkBasis Gk (...);
   SxPW waves (nStates, nSpin, Gk);
   waves.randomize ();
   waves.orthonormalize ();
\endcode
      */
      void randomize ();

      virtual SxPW &normalize ();
      /**
        Performs orthogonalization either according to Gram-Schmidt method
        or L�din scheme. 
        \par (1) Gram-Schmidt method

        According to Gram-Schmidt a orthogonalized set of vectors 
        \f$ \{\tilde{\Psi}\} \f$ can be obtained from a set of 
        non-orthogonal vectors \f$ \{\Psi\} \f$ by
        \f[
           \tilde{\Psi}_i = \Psi_i 
                          - \sum_{j=i}^n \Psi_j^\dagger \Psi_i  \Psi_j.
        \f].

        \par (2) L�din-scheme
        \f[
           \tilde{\Psi}^\dagger = U^\frac{1}{2}  \Psi^\dagger 
        \f]
        where 
        \f[ 
           U^\frac{1}{2} = \Psi^\dagger 1_\varepsilon \Psi, \qquad
           S \Psi = \varepsilon \Psi. 
        \f]
        Here S denotes the overlap matrix. 
        A detailed description you'll find at SxMatrix::getU.
        \brief Orthogonalization of a set of vectors.
        \note In contrast to Gram-Schmidt in L�din method matrices will 
        be used rather than vectors.
        \todo    In L�din part there are still too many copy operations!
        \see     SxPW::orthonormalize.
        \returns Reference to orthogonalized \f$ { \tilde{\Psi} } \f$
        */
      SxPW &orthogonalize  (enum SxOrthoMethod method = GramSchmidt);
      SxPW &orthogonalize  (SxPW &in);
      /**
        This routine orthonormalizes a given set of vectors.
        \brief Orthonormalization of a set of vectors.
        @param uPtr stores the transformation matrices, if required
        \see   SxPW::orthogonalize
        */
      SxPW &orthonormalize (enum SxOrthoMethod method=GramSchmidt,
                       SxArray<SxArray<SxDiracMat<TPrecCoeffG> > > *uPtr=NULL);

      /** \brief Normalization procedures */
      enum NormMethod { DONT_NORMALIZE, NORMALIZE, RENORMALIZE };

      /** \brief Set a vector orthogonal to the first n states

          \param psiPtr     the vector to be orthonormalized
          \param firstN     number of states to be used for orthogonalization
          \param iSpin      spin channel
          \param ik         k-point
          \param norm
                 - DONT_NORMALIZE do not normalize
                 - NORMALIZE      normalize
                 - RENORMALIZE    normalize, assume <psi|psi> = 1.
          \param threshold  orthogonal means  <i|psi> <= threshold

          \author C. Freysoldt, freyso@fhi-berlin.mpg.de
        */
      void setOrthogonal (SxDiracVec<TPrecCoeffG> *psiPtr,
                         int firstN, int iSpin, int ik,
                         enum NormMethod normal,
                         double threshold = 1e-14);

      bool isOrthogonal ();
      bool isOrthonormal ();

      virtual void read (const SxBinIO &io, int mode = SxPWSet::Default);
      void readHDF5 (const SxString &file, int mode = SxPWSet::Default);
      virtual void write (SxBinIO &io) const;

      size_t getNBytes () const;

      void setZero ();


   protected:
      /// Number of states per k-point (:nk)
      SxList<int>  nStatesPerK;   // :ik
      /// Number of coefficients per k-point (:nk)
      SxList<int>  nGIPerK;   // :ik
      int nkPoints;
      int nStates;
      mutable int nSpin;


      /// \name Storing waves on disk
      //{
      /// The currently loaded k-point (-1 if none)
      mutable int  loadedK;
      /// The currently loaded k-point (-1 if none)
      mutable int  loadedSpin;
      /** \brief Internal storage model */
      mutable enum StorageModel stored;
      /// The scratch file directory
      SxString     tmpDir;
      /// The scratch file
      mutable SxBinIO wavesFile;
      //}

      /// \name Storing waves on disk
      //{
      /** \brief Store waves to disk
          \sa StorageModel
        */
      void flushWaves () const;
      /** \brief Read waves from disk
          \sa StorageModel
        */
      void loadWaves (int ik,int iSpin) const;
      /** \brief Read single wave function from disk
          \sa StorageModel
        */
      void loadWaves (int i, int iSpin, int ik) const;
      /** \brief Read single wave function from netcdf file
        */
      static PsiG readPsi (const SxBinIO &io, int i, int iSpin, int ik);

      /** \brief Creates the scratchfile for #StorageModel #KeepOnDisk */
      void createScratchFile ();
      //}

      void registerMemoryObservers ();
   public:

      /// Whether waves are kept on disk
      virtual bool wavesOnDisk () const
      {
         return (stored == KeepOnDisk);
      }
      /** \brief Minimize memory consumption

          If #StorageModel is #KeepOnDisk, unload the current k-point.
       */
      virtual void memMinimize ();
//#define SXPW_FINDLOOP
#ifdef SXPW_FINDLOOP
   protected:
      mutable int findLoopLastIk, findLoopLastIspin;
      void findLoopUpdate (int iSpin, int ik) const;
   public:
      void findLoopReset () const;
#endif
};

#ifdef SXPW_FINDLOOP
inline void SX_NEW_LOOP(const SxPW &pw)
{
   pw.findLoopReset ();
}

inline void SX_NEW_LOOP (const SxPWSet &pw)
{
   const SxPW *ptr = dynamic_cast<const SxPW *>(&pw);
   if (ptr) ptr->findLoopReset ();
}
#else
inline void SX_NEW_LOOP (const SxPWSet &) {}
#endif

//SxBinIO &operator<< (SxBinIO &, const SxPW &);
//SxBinIO &operator>> (SxBinIO &,       SxPW &);


template<>
inline size_t getNBytes<SxPW> (const SxPW &in)
{
   return in.getNBytes ();
}

namespace Timer {
   enum PWTimer {
      ortho,
      loewdinS
   };
}

SX_REGISTER_TIMERS (Timer::PWTimer)
{
   using namespace Timer;
   regTimer (ortho,      "Orthogonalization");
   regTimer (loewdinS,   "Loewdin S");
}




#endif /* _SX_PW_SET_H_ */
