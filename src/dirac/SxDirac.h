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

#ifndef _SX_DIRAC_H_
#define _SX_DIRAC_H_

#include <SxDiracLib.h>
#include <stdio.h>
#include <SxError.h>
#include <atomic>

#define NONE_M -1000

class SxBasis;

/** \brief Auxillary information for S/PHI/nX Dirac Vectors

    \b SxDiracAux = S/PHI/nX Auxiliary Data for Dirac Vectors

    Besides the numerical vector array S/PHI/nX Dirac Vectors need to know
    addition information about the basis-set on which they are defined.
    This class contains those information.
    Beside the basis-set pointer it also contains the quantum numbers to
    specify the Dirac vector. That can be \f$(i,\sigma,\bf{k})\f$ for
    plane-waves or \f$(i_s,i_a,n,l,m)\f$ for atomic orbitals and so on.

    <b>Adding new basis-sets:</b>

    If a new basis-set requires a new type of Dirac vectors changes in this
    class have to be done only if the set of quantum numbers is insufficient
    to describe the new vector. In this case a new (integer) quantum number
    has to be added. In the SxDiracAux::initialize routine the new quantum
    number has to be initialized with a intentionally wrong value (e.g. -1).
    This initialization is used in the DEBUG mode of S/PHI/nX to identify 
    uninitialized S/PHI/nX Dirac vectors.
  
    \ingroup  group_dirac
    \sa       \ref page_dirac
    \author   Sixten Boeck
  */
class SxDiracAux
{
   public:

      /// \brief The band or state index (or -1 if undefined)
      int i;
      /// \brief The index of the spin channel (or -1 if undefined)
      int iSpin;
      /// \brief The index of the \b k point (or -1 if undefined)
      int ik;
      /// \brief The index of the species (or -1 if undefined)
      int is;
      /// \brief The index of the atom of that species (or -1 if undefined)
      int ia;
      /// \brief The main quantum number (or -1 if undefined)
      int n;
      /// \brief The angular quantum number (or -1 if undefined)
      int l;
      /// \brief The magnetic projection number (or -1 if undefined)
      int m;

      /// \brief The pointer to the basis class
      const SxBasis *basisPtr;

      /// \brief Standard constructor, initialize the quantum numbers with -1.
      inline SxDiracAux ()     { init (); }

      /** \brief Used for internal reasons.
        
          This function is necessary for *SxVec::auxData = 0 calls. */
      inline SxDiracAux (int)  { init (); }

      /// \brief Destructor. 
      inline ~SxDiracAux ()    { /* empty */ }

      /** \brief actual constructor

          This function is called from the vector class because "malloc" is
          use instead of "new". Additionally we don't want the assingment 
          operator to be called either. Hence, the standard constructor is 
          not called. */
      inline void init ()  {
         i = iSpin = ik = is = ia = n = l = m = -1; 
         basisPtr = NULL;
      }

      /// Minimal initialization function
      inline void initPtr () { basisPtr = NULL; }
      
      inline SxDiracAux &operator= (const SxDiracAux &in)  {
         //if (!basisPtr)  basisPtr = in.basisPtr;
         if (in.basisPtr)  basisPtr = in.basisPtr;
         i  = in.i; iSpin = in.iSpin; ik = in.ik;
         is = in.is; ia = in.ia; n = in.n; l = in.l; m = in.m;
         return *this;
      }

      /** for SxVec::operator^ */
      inline bool operator== (const SxDiracAux &)  {
         //return basisPtr->getType() == in.basisPtr->getType();
         return true;
      }
};

inline SxDiracAux assign (const SxDiracAux &a, const SxDiracAux &b)
{
//#  ifndef NDEBUG
//      if (a.basisPtr && b.basisPtr)  
//         SX_CHECK (a.basisPtr == b.basisPtr);
//#  endif
   if (a.basisPtr)  return a;
   else             return b;      
}


struct SxDiracHandle 
{
   std::atomic<int> refCounter;
   ssize_t size;
   ssize_t nRows;
   ssize_t nCols;
   char parameters;
   SxDiracAux auxData;
};



#ifdef SXVEC
#   undef SXVEC
#endif
#define SXVEC        SxDiracVec

#ifdef SXMATBASIS
#   undef SXMATBASIS
#endif
#define SXMATBASIS   SxDiracVec

#ifdef SXMAT
#   undef SXMAT
#endif
#define SXMAT        SxDiracMat

#ifdef SXSYMMAT    
#   undef SXSYMMAT   
#endif
#define SXSYMMAT     SxDiracSymMat

#ifdef SXAUX
#   undef SXAUX
#endif
#define SXAUX        SxDiracAux

#define SXAUXHH      <SxDirac.hh>
#undef SXAUXHPP
#define SXAUXHPP     <SxDirac.hpp>

#ifdef MEMHANDLE
#   undef MEMHANDLE
#endif
#define MEMHANDLE    SxDiracHandle

#include <SxVec.h>
#include <SxMat.h>
#include <SxSymMat.h>
#include <SxVector.h>

template<class T>
SxVector<T> toVector (const SxDiracVec<T> &in)
{
   ssize_t n = in.getSize();
   SxVector<T> res (n);

   typename T::Type *srcPtr = in.elements;
   typename T::Type *dstPtr = res.elements;

   for (ssize_t i=0; i < n; i++)  *dstPtr++ = *srcPtr++;

   if (in.nCols () > 0)
      res.reshape (in.nRows (), in.nCols ());
   else
      res.reshape (in.nRows ());

   return res;
}

template<class T>
SxDiracVec<T> toVector (const SxVector<T> &in)
{
   ssize_t n = in.getSize();
   SxDiracVec<T> res (n);

   typename T::Type *srcPtr = in.elements;
   typename T::Type *dstPtr = res.elements;

   for (ssize_t i=0; i < n; i++)  *dstPtr++ = *srcPtr++;

   if (in.nCols () > 0)
      res.reshape (in.nRows (), in.nCols ());
   else
      res.reshape (in.nRows ());

   return res;
}

template <class T>
class SxVecTraits<SxDiracVec<T> >
{
   public:
   typedef T ScalarType;
   typedef SxDiracMat<T> MatType;
};

/// Writes an x-y1-y2... (up to 3 y's) plot to a file
SX_EXPORT_DIRAC void writePlot (const SxString &fileName,
                const SxDiracVec<Double> &x,
                const SxDiracVec<Double> &y1,
                const SxDiracVec<Double> &y2 = SxDiracVec<Double> (),
                const SxDiracVec<Double> &y3 = SxDiracVec<Double> () );

/// Writes an (dx*i)-y1-y2... (up to 2 y's) plot to a file
SX_EXPORT_DIRAC void writePlot (const SxString &fileName,
                double dx,
                const SxDiracVec<Double> &y1,
                const SxDiracVec<Double> &y2 = SxDiracVec<Double> ());

#endif /* _SX_DIRAC_H_ */
