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

#ifndef _SX_N_ARRAY_H_
#define _SX_N_ARRAY_H_

#include <SxArray.h>

template <unsigned N>
class SxNTuple
{
   public:
      ssize_t val[N];
      inline SxNTuple<N> &operator= (const SxNTuple<N> &in);
      inline SxNTuple<N> &operator= (ssize_t x);
      SxNTuple (const SxNTuple<N> &in);
      inline ssize_t &operator[] (unsigned);
      inline ssize_t operator[] (unsigned) const;

      SxNTuple ();
      SxNTuple (const SxNTuple<N-1> &, ssize_t);
      SxNTuple (ssize_t, const SxNTuple<N-1> &);
      // template <class... Rest>
      // SxNTuple (ssize_t first, Rest... rest)
      // {
      //    val[0] = first;
      //    // use placement new for Rest
      //    new (val+1) SxNTuple<N-1> (rest...);
      // }
      ssize_t product () const;

      bool operator> (ssize_t n)
      {
         for (unsigned int i = 0; i < N; ++i)
            if (val[i] <= n) return false;
         return true;
      }
};

template<unsigned N>
inline ssize_t& SxNTuple<N>::operator[] (unsigned i)
{
   SX_CHECK (i < N, i, N);
   return val[i];
}

template<unsigned N>
inline ssize_t SxNTuple<N>::operator[] (unsigned i) const
{
   SX_CHECK (i < N, i, N);
   return val[i];
}

template<unsigned N>
SxNTuple<N> &SxNTuple<N>::operator= (const SxNTuple<N> &in)
{
   for (unsigned int i = 0; i < N; ++i) val[i] = in.val[i];
   return *this;
}

template<unsigned N>
SxNTuple<N> &SxNTuple<N>::operator= (ssize_t x)
{
   for (unsigned int i = 0; i < N; ++i) val[i] = x;
   return *this;
}

template<unsigned N>
SxNTuple<N>::SxNTuple ()
{
   for (unsigned int i = 0; i < N; ++i) val[i] = -1;
}

template<unsigned N>
SxNTuple<N>::SxNTuple (const SxNTuple<N> &in)
{
   operator= (in);
}

template<unsigned N>
SxNTuple<N>::SxNTuple (const SxNTuple<N-1> &in, ssize_t lastVal)
{
   reinterpret_cast<SxNTuple<N-1>*>(this)->operator= (in);
   val[N-1]=lastVal;
}

template<unsigned N>
SxNTuple<N>::SxNTuple (ssize_t firstVal, const SxNTuple<N-1> &in)
{
   val[0]=firstVal;
   reinterpret_cast<SxNTuple<N-1>&>(val[1])->operator= (in);
}

template<>
class SxNTuple<1>
{
   ssize_t val[1];
   SxNTuple (ssize_t in)
   {
      val[0] = in;
   }

   void operator= (const SxNTuple<1> &in)
   {
      val[0] = in.val[0];
   }

};

template<unsigned N>
ssize_t SxNTuple<N>::product () const
{
   ssize_t res = val[0];
   for (unsigned int i = 1; i < N; ++i) res *= val[i];
   return res;
}

/** \brief Multidimensional array

    \author C. Freysoldt, freysoldt@mpie.de */
template <class T, unsigned N>
class SxNArray : public SxArray<T>
{
   public:
      /// Empty constructor
      SxNArray () { dim = 0; }

      /// Deprecate element access from base class (not implemented)
      inline T& operator()(ssize_t);
      /// Deprecate element access from base class (not implemented)
      inline const T& operator()(ssize_t) const;

      /// reformat (generic)
      void reformat (const SxNTuple<N> &dimIn);
      /// Element access (generic)
      T& operator()(const SxNTuple<N> &idx);
      /// Element access (generic)
      const T& operator()(const SxNTuple<N> &idx) const;

      /// Resize (only 0 allowed)
      inline void resize (ssize_t n)
      {
         SX_CHECK (n==0, n);
         SxArray<T>::resize (0);
         dim = 0;
      }

   protected:
      /// Dimensions
      SxNTuple<N> dim;

      // --- Generic index generation functions (via variadic templates)
      /*
      template <unsigned Dim, class... SizeT>
      ssize_t getIdx (ssize_t i, SizeT... rest) const inline
      {
         return i + dim[N-Dim] * getIdx<N-1>(iDim+1, rest...);
      }

      template <>
      ssize_t getIdx<1> (unsigned iDim, ssize_t i) const inline
      {
         return i;
      }
      */

      // Generic reformat (via variadic templates)
      /*
      template <class... SizeT>
      void reformat (SizeT... dims)
      {
         dim = SxNTuple<N> (dims);
         SxArray<T>::resize (dim.product ());
      }
      */

   public:
      // --- Generic access functions (via variadic templates)
      /*
      template <class... SizeT>
      T& operator () (SizeT... args)
      {
         return SxArray<T>::operator()(getIdx<N> (args...));
      }

      template <class... SizeT>
      const T& operator () (SizeT... args) const
      {
         return SxArray<T>::operator()(getIdx<N> (args...));
      }
      */

      inline ssize_t getDim(unsigned i) const
      {
         SX_CHECK (i < N, i, N);
         return dim[i];
      }

      // --- Generated code (until variadic templates work)
      /*
#!/usr/bin/perl -w

my $n = shift;
exit unless defined $n;

my $parList = "X1";
my ($computeIdx, $closeBracket) = ("i1", "");
my ($dimCode) = "dim[0] = n1;";
my ($prod) = "n1";
my ($check_idx) = "SX_CHECK (i1 >= 0 && i1 < dim[0], i1, dim[0]);";

for my $i (2 .. $n)  {
   $parList .= ", X$i";
   my $j = $i-2;
   $computeIdx = $computeIdx . " + dim[$j] * (i$i";
   $closeBracket .= ")";
   my $k = $i-1;
   $dimCode .= "\n         dim[$k] = n$i;";
   $prod .= " * n$i";
   $check_idx .= "\n         SX_CHECK (i$i >= 0 && i$i < dim[$k], i$i, dim[$k]);"
}
$computeIdx .= $closeBracket;

my ($iParList, $nParList, $iArgList, $nArgList)
    = ($parList, $parList, $parList, $parList);
$iParList =~ s/X/ssize_t i/g;
$nParList =~ s/X/ssize_t n/g;
$iArgList =~ s/X/i/g;
$nArgList =~ s/X/n/g;

print <<END_CODE ;
      // --- generated code for N=$n
      size_t getIdx ($iParList) const
      {
         SX_CHECK (N == $n, N);
         $check_idx
         return $computeIdx;
      }

      T& operator() ($iParList)
      {
         SX_CHECK (N == $n, N);
         return SxArray<T>::operator()(getIdx ($iArgList));
      }

      const T& operator() ($iParList) const
      {
         SX_CHECK (N == $n, N);
         return SxArray<T>::operator()(getIdx ($iArgList));
      }

      void reformat ($nParList)
      {
         SX_CHECK (N == $n, N);
         SxArray<T>::resize ($prod);
         $dimCode
      }

      SxNArray ($nParList)
      {
         SX_CHECK (N == $n, N);
         SxArray<T>::resize ($prod);
         $dimCode
      }

END_CODE

       */
      // --- generated code for N=2
      size_t getIdx (ssize_t i1, ssize_t i2) const
      {
         SX_CHECK (N == 2, N);
         SX_CHECK (i1 >= 0 && i1 < dim[0], i1, dim[0]);
         SX_CHECK (i2 >= 0 && i2 < dim[1], i2, dim[1]);
         return i1 + dim[0] * (i2);
      }

      T& operator() (ssize_t i1, ssize_t i2)
      {
         SX_CHECK (N == 2, N);
         return SxArray<T>::operator()(getIdx (i1, i2));
      }

      const T& operator() (ssize_t i1, ssize_t i2) const
      {
         SX_CHECK (N == 2, N);
         return SxArray<T>::operator()(getIdx (i1, i2));
      }

      void reformat (ssize_t n1, ssize_t n2)
      {
         SX_CHECK (N == 2, N);
         SxArray<T>::resize (n1 * n2);
         dim[0] = n1;
         dim[1] = n2;
      }

      SxNArray (ssize_t n1, ssize_t n2)
      {
         SX_CHECK (N == 2, N);
         SxArray<T>::resize (n1 * n2);
         dim[0] = n1;
         dim[1] = n2;
      }

      // --- generated code for N=3
      size_t getIdx (ssize_t i1, ssize_t i2, ssize_t i3) const
      {
         SX_CHECK (N == 3, N);
         SX_CHECK (i1 >= 0 && i1 < dim[0], i1, dim[0]);
         SX_CHECK (i2 >= 0 && i2 < dim[1], i2, dim[1]);
         SX_CHECK (i3 >= 0 && i3 < dim[2], i3, dim[2]);
         return i1 + dim[0] * (i2 + dim[1] * (i3));
      }

      T& operator() (ssize_t i1, ssize_t i2, ssize_t i3)
      {
         SX_CHECK (N == 3, N);
         return SxArray<T>::operator()(getIdx (i1, i2, i3));
      }

      const T& operator() (ssize_t i1, ssize_t i2, ssize_t i3) const
      {
         SX_CHECK (N == 3, N);
         return SxArray<T>::operator()(getIdx (i1, i2, i3));
      }

      void reformat (ssize_t n1, ssize_t n2, ssize_t n3)
      {
         SX_CHECK (N == 3, N);
         SxArray<T>::resize (n1 * n2 * n3);
         dim[0] = n1;
         dim[1] = n2;
         dim[2] = n3;
      }

      SxNArray (ssize_t n1, ssize_t n2, ssize_t n3)
      {
         SX_CHECK (N == 3, N);
         SxArray<T>::resize (n1 * n2 * n3);
         dim[0] = n1;
         dim[1] = n2;
         dim[2] = n3;
      }

      // --- generated code for N=4
      size_t getIdx (ssize_t i1, ssize_t i2, ssize_t i3, ssize_t i4) const
      {
         SX_CHECK (N == 4, N);
         SX_CHECK (i1 >= 0 && i1 < dim[0], i1, dim[0]);
         SX_CHECK (i2 >= 0 && i2 < dim[1], i2, dim[1]);
         SX_CHECK (i3 >= 0 && i3 < dim[2], i3, dim[2]);
         SX_CHECK (i4 >= 0 && i4 < dim[3], i4, dim[3]);
         return i1 + dim[0] * (i2 + dim[1] * (i3 + dim[2] * (i4)));
      }

      T& operator() (ssize_t i1, ssize_t i2, ssize_t i3, ssize_t i4)
      {
         SX_CHECK (N == 4, N);
         return SxArray<T>::operator()(getIdx (i1, i2, i3, i4));
      }

      const T& operator() (ssize_t i1, ssize_t i2, ssize_t i3, ssize_t i4) const
      {
         SX_CHECK (N == 4, N);
         return SxArray<T>::operator()(getIdx (i1, i2, i3, i4));
      }

      void reformat (ssize_t n1, ssize_t n2, ssize_t n3, ssize_t n4)
      {
         SX_CHECK (N == 4, N);
         SxArray<T>::resize (n1 * n2 * n3 * n4);
         dim[0] = n1;
         dim[1] = n2;
         dim[2] = n3;
         dim[3] = n4;
      }

      SxNArray (ssize_t n1, ssize_t n2, ssize_t n3, ssize_t n4)
      {
         SX_CHECK (N == 4, N);
         SxArray<T>::resize (n1 * n2 * n3 * n4);
         dim[0] = n1;
         dim[1] = n2;
         dim[2] = n3;
         dim[3] = n4;
      }

};

template <class T, unsigned N>
void SxNArray<T,N>::reformat (const SxNTuple<N> &dimIn)
{
   SX_CHECK (dimIn > 0);
   SxArray<T>::resize (dimIn.product ());
   dim = dimIn;
}

template <class T, unsigned N>
T& SxNArray<T,N>::operator()(const SxNTuple<N> &idx)
{
   ssize_t idx1D = 0;
   for (int i = N; i; )  {
      --i;
      idx1D = idx[i] + dim[i] * idx1D;
   }
   return SxArray<T>::operator()(idx1D);
}

template <class T, unsigned N>
const T& SxNArray<T,N>::operator()(const SxNTuple<N> &idx) const
{
   ssize_t idx1D = 0;
   for (int i = N; i; )  {
      --i;
      idx1D = idx[i] + dim[i] * idx1D;
   }
   return SxArray<T>::operator()(idx1D);
}

template <class T>
class SxArray2 : public SxNArray<T, 2> {
   public:
      SxArray2 (int n1, int n2) : SxNArray<T, 2> (n1, n2)
      {
         // empty
      }

      SxArray2 () { /* empty */ }

      T& operator() (const SxAutoLoop &i1, const SxAutoLoop &i2)
      {
         i1.setLimit (this->dim[0]);
         i2.setLimit (this->dim[1]);
         return SxArray<T>::operator()(this->getIdx (i1.i, i2.i));
      }

      const T& operator() (const SxAutoLoop &i1, const SxAutoLoop &i2) const
      {
         i1.setLimit (this->dim[0]);
         i2.setLimit (this->dim[1]);
         return SxArray<T>::operator()(this->getIdx (i1.i, i2.i));
      }

};

template <class T>
class SxArray3 : public SxNArray<T, 3> {
   public:
      SxArray3 (int n1, int n2, int n3)
         : SxNArray<T, 3> (n1, n2, n3)
      {
         // empty
      }
      SxArray3 () { /* empty */ }
};

#endif /* _SX_N_ARRAY_H_ */
