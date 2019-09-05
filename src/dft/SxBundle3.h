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
#ifndef _SX_BUNDLE3_H_
#define _SX_BUNDLE3_H_

#include <SxList.h>
#include <SxDirac.h>
#include <SxTypeMapper.h>
#include <SxError.h>
//#include <SxPrecision.h>
#include <SxBinIO.h>


/** \brief Template container class to store data matching to \f$
              (i,\sigma,\bf{k})
           \f$, e.g. one-particle energies, occupancies, etc.

    \b SxBundle3 = SPHInX Bundle Container for 3 dimensions

    This class actually is a three dimensional SxArray<T>. It will be
    used e.g. for the occupancy numbers or the set wavefunction.

    \author Sixten Boeck
  */
template<class T>
class SxBundle3
{
   public:

      /** The data container */
      SxArray<SxArray<SxDiracVec<T> > >  bundle;   // :ik,:iSpin,:i

      /// \brief Standard constructor
      SxBundle3 () { };
      /** \brief Allocate the bundle array.
       
          The bundle array is allocated with \em nStates states,
          \em nSpin spin channels, and \em ik \b k points. */
      SxBundle3 (int nStates, int nSpin, int nk);

      /// \brief Copy constructor
      SxBundle3 (const SxBundle3<T> &);

      /// \brief Destructor
      ~SxBundle3 ();

      /// \brief Standard assignment operator
      SxBundle3<T> &operator= (const SxBundle3<T> &);
         
      /** \brief Extracts all elements belonging to \f$(\sigma,\bf{k})\f$ */
      inline SxDiracVec<T> &operator() (ssize_t iSpin, ssize_t ik);
      inline SxDiracVec<T> &
      operator() (const SxAutoLoop &iSpin, const SxAutoLoop &ik)
      {
         iSpin.setLimit (getNSpin ());
         ik.setLimit (getNk ());
         return operator()(iSpin.i, ik.i);
      }
      /** \brief Extracts all elements belonging to \f$(\sigma,\bf{k})\f$ */
      inline const SxDiracVec<T> &operator() (ssize_t iSpin, ssize_t ik) const;
      inline const SxDiracVec<T> &
      operator() (const SxAutoLoop &iSpin, const SxAutoLoop &ik) const
      {
         iSpin.setLimit (getNSpin ());
         ik.setLimit (getNk ());
         return operator()(iSpin.i, ik.i);
      }
      /** \brief Extracts the element belonging to \f$(i,\sigma,\bf{k})\f$ */
      inline typename T::Type &
      operator() (ssize_t i, ssize_t iSpin, ssize_t ik);
      /** \brief Extracts the element belonging to \f$(i,\sigma,\bf{k})\f$ */
      inline const typename T::Type &
      operator() (ssize_t i, ssize_t iSpin, ssize_t ik) const;
      /** \brief Extracts the element belonging to \f$(i,\sigma,\bf{k})\f$ */
      inline typename T::Type &
      operator() (SxAutoLoop &i, SxAutoLoop &iSpin, SxAutoLoop &ik)
      {
         i.setLimit (getNStates ((int)ik.i));
         iSpin.setLimit (getNSpin ());
         ik.setLimit (getNk ());
         return operator()(i.i, iSpin.i, ik.i);
      }
      /** \brief Extracts the element belonging to \f$(i,\sigma,\bf{k})\f$ */
      inline const typename T::Type &
      operator() (SxAutoLoop &i, SxAutoLoop &iSpin, SxAutoLoop &ik) const
      {
         i.setLimit (getNStates ((int)ik.i));
         iSpin.setLimit (getNSpin ());
         ik.setLimit (getNk ());
         return operator()(i.i, iSpin.i, ik.i);
      }
      /** \brief Initialize all elements with a specified value */
      void set (const typename T::Type &);

      /** \brief Save element array */
      void write (SxBinIO &) const;
      /** \brief Load element array */
      void read  (SxBinIO &);

      /** \brief Returns the number of states */
      inline int getNStates (int ik=0) const;

      /** \brief Returns the number of spins */
      inline int getNSpin (int ik=0) const;
      
      /** \brief Returns the number of k-points */
      inline int getNk () const;

      inline void synMPI ();  // LoopMPI
};

#include <SxBundle3.hpp>

#endif /* _SX_BUNDLE3_H_ */
