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

#ifndef _SX_DIR_TENSOR3_H_
#define _SX_DIR_TENSOR3_H_

#include <SxMatrix.h>
#include <SxUtil.h>
#include <SxString.h>
#include <SxBinIO.h>

/** \brief This is a container for 3-dimensional matrices.

    \b SxTensor3 = SPHInX 3-dimesional matrix

    Later this class will be modified.

    \author Matthias Wahn, wahn@mpie.de */
template<class T>
class SxDirTensor3
{
   public:
      /** data container */
      SxArray<SxDiracMat<T> > tensor;  // e.g.:  :ik, :n, :m

      /** dimensions */
      int n1, n2, n3;

      /** empty constructor */
      SxDirTensor3 ();
      /** copy constructor */
      SxDirTensor3 (const SxDirTensor3<T> &t);
      /** standard constructor */
      SxDirTensor3 (int n1_, int n2_, int n3_);
      /** destructor */
      ~SxDirTensor3 ()  { /* empty */ };

      /** reformat */
      void reformat (int n1_, int n2_, int n3_);

      /** assign a matrix to a k index */
      void set (const SxDiracMat<T> &mat, int i1);

      /** assign a tensor element to the index triple (i1,i2,i3) */
      void set (const typename T::Type &elem, int i1, int i2, int i3);

      /** set all tensor elements to a certain value */
      void setAllElem (const typename T::Type &in);

      /** set the diagonals of each k-matrix to a certain number, e.g.
          to 1 to obtain a set of unit matrices.
          
          \note This function does not set the tensor's "space diagonal",
          i.e. the elements with k = n = m. */
      void setDiagonals (const typename T::Type &in);

      /** get i1-th matrix */
      SxDiracMat<T> &operator() (int i1);
      const SxDiracMat<T> &operator() (int i1) const;

      /** get (i1,i2,i3)-th element */
      typename T::Type &operator() (int i1, int i2, int i3);
      const typename T::Type &operator() (int i1, int i2, int i3) const;

      /** assign tensor3 to tensor3 */
      inline SxDirTensor3<T> operator= (const SxDirTensor3<T> &in);

      /** print the tensor matrixwise --- for debugging */
      void print ();
      void print (int z);              // prints the first z matrices
      void print (const SxString &s);  // prints the tensor's name
      void print (const SxString &s,   // prints name and dumps
                  int z);              //    first z matrices

      //-----------------------------------------------------------------------
      /**@name I/O */
      //-----------------------------------------------------------------------
      void writeDirTensor3 (const SxString &file) const;
      void writeDirTensor3 (const SxBinIO &io) const;
      void readDirTensor3 (const SxString &file);
      void readDirTensor3 (const SxBinIO &io);
      void readTensor3 (const SxString &file);
      void readTensor3 (const SxBinIO &io);
      /** \brief obsolete -- just to be compatible with the old SPHInX */
      void readUkmn (const SxString &file);
};

/** elementwise sum of two tensors */
template<class A, class B>
inline SxDirTensor3<typename SxTypeCaster<A,B>::TRes>
operator+ (const SxDirTensor3<A> &a, const SxDirTensor3<B> &b)
{
   int n1 = a.n1, n2 = a.n2, n3 = a.n3;

   SX_CHECK (a.n1 == b.n1, a.n1, b.n1);
   SX_CHECK (a.n2 == b.n2, a.n2, b.n2);
   SX_CHECK (a.n3 == b.n3, a.n3, b.n3);

   SxDirTensor3<typename SxTypeCaster<A,B>::TRes> res(n1,n2,n3);

   for (int i1 = 0; i1 < n1; i1++)  {
      res.set ( a(i1) + b(i1), i1 );
   }

   return res;
}

/** &Uuml;berschiebung without sum
    (If anybody knows the English expression, feel free to replace the
    word. ;-)) */
template<class A, class B>
inline SxDirTensor3<typename SxTypeCaster<A,B>::TRes>
operator* (const SxDirTensor3<A> &a, const SxDirTensor3<B> &b)
{
   int n1 = a.n1, n2 = a.n2, n3 = a.n3;

   SX_CHECK (a.n1 == b.n1, a.n1, b.n1);
   SX_CHECK (a.n2 == b.n2, a.n2, b.n2);
   SX_CHECK (a.n3 == b.n3, a.n3, b.n3);

   SxDirTensor3<typename SxTypeCaster<A,B>::TRes> res(n1,n2,n3);
   SxDiracMat<typename SxTypeCaster<A,B>::TRes> mat;

   for (int i1 = 0; i1 < n1; i1++)  {
      mat = (a(i1) ^ b(i1));
      res.set (mat, i1);
   }

   return res;
}

#include <SxDirTensor3.hpp>

#endif /* _SX_DIR_TENSOR3_H_ */
