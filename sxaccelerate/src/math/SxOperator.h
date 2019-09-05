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

#ifndef _SX_OPERATOR_H_
#define _SX_OPERATOR_H_

#include <SxPtr.h>

template<class T> class SxOperator;

/** \brief General operator/transformation or pipeline of 
           operators/transformations

    \b SxOperatorBase = S/PHI/nX general Operator.

    In S/PHI/nX pipelining of operators and transformations is supported.
    It is possible to define a set of operators or transformations and
    concatinate them to form a new operator/transformation.
    A numerical method can now transform its input data or output data by
    calling the newly defined operator. By using this technique the numerical
    algorithm can be defined independently from transformations.

    To define a new operator, you must derive the new operator class from
    SxOperatorBase<T> (or one of the special subclasses mentioned below),
    where T is the type that the operator works on (no 
    type-change operators may be defined). Then, you have to implement
    virtual functions
     - operator*    for out-of-place transformations
     - applyInPlace for in-place transformations data *= op;
     - getCopy      for operator piping

     - getCopy must be implemented if your operator should support
    piping: The resulting composed operator must of course 
    have access to the individual parts. Typically, the single
    operators are copied into SxPtrs via getCopy during
    the concatination (or append).
    If you manage your operator by SxPtr anyway
    you do not need this function.
    
    A standard implementation for getCopy is provided by the
    SXOPERATOR_GETCOPY macro. It assumes your operator class has
    a proper copy constructor.

    \example
    See examples/sxpipes.cpp

    \author Christoph Freysoldt, freysoldt@mpie.de,
            reimplementing the original idea of Sixten Boeck, boeck@mpie.de 
    */
template<class T>
class SxOperatorBase
{
   public:
      /// Apply the operator
      virtual T operator*(const T &in) const = 0;
      /** Apply the operator in place
        @param data the object that the operator is applied to
        */
      virtual void applyInPlace (T &data) const = 0;

      /** \brief Get a autopointer copy of this object
          Implement this function for your operator if you cannot garantuee
          that only SxPtr will manage your operator instances.
          \example
          \code
SxMyOp : public SxOperatorBase<double> {
   public:
      /// Copy constructor
      SxMyOp (const SxMyOp &in);

      // Define standard getCopy
      SXOPERATOR_GETCOPY(SxMyOp, double);

      ...
};
      \endcode
      \note The argument will always be NULL. It is only needed to specify
            the type for multiple-type operators.
      \sa SXOPERATOR_GETCOPY
        */
	  inline virtual SxPtr<SxOperatorBase<T> > getCopy (T *) const;

/** SXOPERATOR_GETCOPY
    \param SxOpClass name of the operator class
    \param T Type on which operator works
  */
#define SXOPERATOR_GETCOPY(SxOpClass, T) \
      inline virtual SxPtr<SxOperatorBase< T > > getCopy (T *) const \
      { return SxPtr< SxOpClass >::create (*this); }

      /// Destructor
      inline virtual ~SxOperatorBase () { /* empty */ }

      /// Apply the operator
      inline T operator|(const T &in) const { return operator*(in); }

      /** \brief Enable type for multi-type operator
      If your operator works for different types, you may
      derive it from SxOperatorBase<?> for each of the types. However,
      you then need to explicitly enable each of the types by this macro.
      \example
      \code
class SxMyOp : public SxMyOperator<int>, public SxMyOperator<double>
{
   public: 
      SXOPERATOR_TYPE (int);
      SXOPERATOR_TYPE (double);
      ...
};
      \endcode
      \note
      This make the operator| of the base class directly visible in the
      derived class. Necessary to resolve ambiguities in multi-type operators

      see B. Stroustrup, The C++ programming language, 
      Addison Wesley, 2000, Chap. 15.2.2, pp393
      */
#define SXOPERATOR_TYPE(T) \
      using SxOperatorBase<T >::operator|

#ifdef MSVC
	  class Indirect : public SxOperatorBase<T> {};
#define SXOP_LINKFIX ::Indirect
#else
#define SXOP_LINKFIX
#endif
   protected:
      /// Constructor
      inline SxOperatorBase () { /* empty */ }

};

template<class T>
T &operator*= (T &data, const SxOperatorBase<T> &op)
{
   op.applyInPlace (data);
   return data;
}


/** \brief General operator/transformation or pipeline of 
           operators/transformations

    \b SxOperator = S/PHI/nX general Operator pipe.

    In S/PHI/nX pipelining of operators and transformations is supported.
    It is possible to define a set of operators or transformations and
    concatinate them to form a new operator/transformation.
    A numerical method can now transform its input data or output data by
    calling the newly defined operator. By using this technique the numerical
    algorithm can be defined independently from transformations.

    \par Introduction: Structure optimization
    Structure optimization is in principle nothing but a multidimensional
    minimization. But, in order to accomplish a good performance or to 
    calculate specific details input forces and/or input geometries are 
    transformed before entering the minimization equations.
    This class allows to deal with abstract operators which are applied
    to e.g. geometries or forces. The actual definition of those operators
    is well seperated from the minimization routines.
    It is even possible to setup the operator pipeline from the input file
    and thus, generate a new operator on-the-fly without recompilation of the
    code.

    To define a new operator, you must derive the new operator class from
    SxOperatorBase<T> (see there for more information).

    \example
    See examples/sxpipes.cpp

    \author Christoph Freysoldt, freysoldt@mpie.de,
            reimplementing the original idea of Sixten Boeck, boeck@mpie.de 
  */
template<class T>
class SxOperator : public SxOperatorBase<T>
{
   private:
      /// The operators (first to be applied)
      SxConstPtr<SxOperatorBase<T> > first;
      /// The rest of the pipe
      SxConstPtr<SxOperatorBase<T> > next;

   public:
      /// Constructor 
      SxOperator() { /* empty */ }
      /// Constructor 
      SxOperator(const SxOperatorBase<T> &in) 
        : first(in.getCopy ((T*)NULL)) { /* empty */ }
      /// Constructor 
      SxOperator(const SxConstPtr<SxOperatorBase<T> > &in) 
        : first(in) { /* empty */ }

      /// Assignment operator
      SxOperator<T> & operator= (const SxOperatorBase<T> &);
      /// Assignment operator (as reference)
      SxOperator<T> & operator= (const SxConstPtr<SxOperatorBase<T> >&);

      /// Get a copy
      SXOPERATOR_GETCOPY (SxOperator<T>, T)

      /** @{
          \name Piping functionality
        */

      /// Append a further operator
      inline void append (const SxOperatorBase<T> &in)
      {
         append (in.getCopy ((T*)NULL));
      }

      /// Append autopointer to end of pipe
      void append (const SxConstPtr<SxOperatorBase<T> >&);
      
      /// Append a further operator
      SxOperator<T> &operator<< (const SxOperatorBase<T> &in) 
      { 
         append(in.getCopy ((T*)NULL)); 
         return *this;
      }
      /// Append a further operator (as reference)
      SxOperator<T> &operator<< (const SxConstPtr<SxOperatorBase<T> >&in)
      { 
         append(in); 
         return *this;
      }

      /// Prepend operator
      void prepend(const SxOperatorBase<T> &pre)
      {
         prepend (pre.getCopy ((T*)NULL));
      }
      /// Prepend operator
      void prepend(const SxConstPtr<SxOperatorBase<T> > &pre);
      
      /// Prepend operator
      SxOperator<T> &operator|=(const SxOperatorBase<T> &pre)
      {
         prepend (pre.getCopy ((T*)NULL));
         return *this;
      }

      /// Prepend operator
      SxOperator<T> &
      operator|= (const SxConstPtr<SxOperatorBase<T> > &pre)
      {
         prepend (pre);
         return *this;
      }
      ///@}

      /// Apply operations out of place
      virtual T operator* (const T &in) const;
      /// Apply operations in place
      virtual void applyInPlace (T &data) const;

      class Identity;

      /// Check whether operator chain is empty
      bool isEmpty () const { return !first; }

};

template<class T>
class SxOperator<T>::Identity : public SxOperatorBase<T>
{
   public:
      /// do nothing
      virtual T operator* (const T& in) const
      {
         return in;
      }

      /// do nothing
      virtual void applyInPlace (T&) const { /* empty */ }

      /// Get another identity
      inline virtual SxPtr<SxOperatorBase<T> > getCopy (T *) const
      {
         return SxPtr<typename SxOperator<T>::Identity>::create ();
      }

};


#include <SxOperator.hpp>

#endif /* _SX_OPERATOR_H_ */
