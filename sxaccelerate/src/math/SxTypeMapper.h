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
#ifndef _SX_TYPE_MAPPER_H_
#define _SX_TYPE_MAPPER_H_

#include <SxComplex.h>

/**
  This class is necessary for implementation of a algebraic template 
  class to serve \URL[BLAS]{http://www.netlib.org/blas} and 
  \URL[LAPACK]{http://www.netlib.org/lapack}. 
  For some driver routines one has to
  know the corresponding real (that is floating point operation type) 
  or complex datatype to the (at compiling time)
  unknown template <class T> type (Eigensolver, absSqr, etc). 
  In order to provide the compiler with knowlegde about "datatype families"
  this pseudoclass has been generated.
  Without such relations between low level datatypes a central precision
  control wouldn't be easily implementable.

  After including SxTypeMapper the following type definitions can be used:
  \begin{verbatim}
  Int          - family of 4-byte integral types (int)
  Float        - family of single precision real types (float)
  Double       - family of double precision real types (double)
  Complex8     - family of single precision complex types
  Complex16    - family of double precision complex types
  \end{verbatim}

  For more information please refer to \Ref{SxArray}, \Ref{SxVector}, and 
  \Ref{SxMatrix}.

  @author  Sixten Boeck
  @ingroup Numerics
  */
template<class T, class R, class C=SxComplex<T> >
struct SxTypeMapper 
{ 
   /** Type refers to the low level datatype, such as 
       float, SxComplex8,... */
   typedef T                   Type; 
   /** The corresponding real value to type <T>. For real
       types like float or double it refers to itself */
   typedef R                   Real; 
   /** The complex type that fits to <T>. For complex
       types it refers to itself */
   typedef C                   Cmplx;
   /** Gives the datatype family for the corresponding
       complex type */
   typedef SxTypeMapper<C,R,C> TCmplx;
   /** Gives the datatype family for the corresponding
       real type */
   typedef SxTypeMapper<R,R,C> TReal;
};


// --- include specializations of SxTypeMapper
#include <SxTypeCaster.h>


#endif /* _SX_TYPE_MAPPER_H_ */
