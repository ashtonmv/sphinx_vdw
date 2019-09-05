// ---------------------------------------------------------------------------
//
//           The general purpose cross platform C/C++ framework
//
//                       S x A c c e l e r a t e
//
//           Home:       https://www.sxlib.de
//           License:    Apache 2
//           Authors:    see src/AUTHORS
//
// ---------------------------------------------------------------------------

#ifndef _SX_VECTOR_H_
#define _SX_VECTOR_H_

#include <SxTypeDefs.h> /* ssize_t */
#include <atomic>

/** \brief Empty auxiliary data for purely numerical vectors.

    \b SxAuxVec = SFHIngX Auxiliary Vector Data

    The SFHIngX Algebra Class can host additional information for
    each vector or matrix. For example a SFHIngX Dirac vector needs
    to store its quantum numbers as well as a pointer to its basis-set
    object. This class is an \b empty auxiliary data container used for
    plain vectors and matrices.

  \author  Sixten Boeck
  */
class SxAuxVec
{
   public:

      /// \brief Empty standard constructor.
      inline SxAuxVec ()              { /* empty */ }

      /** \brief for internal reasons only.
        
          This function is needed for *SxVec::auxData = 0 calls */
      inline SxAuxVec (int)           { /* empty */ }

      /// \brief Empty destructor
      inline ~SxAuxVec ()             { /* empty */ }

      /** \brief actual constructor

          This function is called from the vector class because "malloc" is
          use instead of "new". Additionally we don't want the assingment 
          operator to be called either. Hence, the standard constructor is 
          not called. */
      inline void init ()             { /* empty */ }

      /// Minimal initialization function
      inline void initPtr () { /* empty */ }

      /// \brief Empty assignment operator
      inline SxAuxVec &operator= (const SxAuxVec &)  {
         return *this;
      }

      /// \brief Empty equality operator
      inline bool operator== (const SxAuxVec &)  {
         return true;
      }
};

struct SxVectorHandle 
{
   std::atomic<int> refCounter;
   ssize_t size;
   ssize_t nRows;
   ssize_t nCols;
   char parameters;
   SxAuxVec auxData;
};

#ifdef SXVEC
#   undef SXVEC
#endif
#define SXVEC     SxVector


#ifdef SXAUX
#   undef SXAUX
#endif
#define SXAUX SxAuxVec

inline SxAuxVec assign (const SxAuxVec &a, const SxAuxVec &)
{
   return a;
}

#ifdef MEMHANDLE
#   undef MEMHANDLE
#endif
#define MEMHANDLE    SxVectorHandle


#include <SxVec.h>

#endif /* _SX_VECTOR_H_ */
