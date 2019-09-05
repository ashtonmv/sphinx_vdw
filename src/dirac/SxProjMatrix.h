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

#ifndef _SX_PROJ_MATRIX_H_
#define _SX_PROJ_MATRIX_H_

#include <SxDirac.h>

/** \brief Generic projector matrix

    \b SxClass = S/PHI/nX projector matrix base class

    This class implements block algorithms for projector matrices
    defined in the basis |B> by
    \f[
    <B|M|B'> = \sum_i <B|\phi_i> f_i <\phi_i|B'>
    \f]
    where number of projectors \f$\phi_i\f$ (nProj) is considerably
    lower than the number of basis function in |B> (nElements).
    Although it is conceptually very convenient and efficient to
    use these projectors rather than the full matrix, the numerical
    implementation is usually awkward if numerical efficiency plays
    a role, notably when the projectors are not stored in one block.

    This is what this class is for. It collects a number of projectors
    (given by blocksize) into a single block and uses then highly efficient
    matrix-vector or matrix-matrix routines for applying the matrix.

    It does however not know how to produce the required projectors
    \f$\phi_i\f$ and the inner factors \f$ f_i \f$. This is why you
    cannot use this class directly.

    In order to use this class
      - derive a class from it, e.g.
        \code
  SxMyClass : public SxProjMatrix<Complex16>
  {
        \endcode
      - add all the necessary variables to generate the projectors and
        an appropiate constructor
        \code
     public:
        const SGkBasis &g;
        SxArray<PsiGI> &referenceOrbitals;
        SxMyClass (const SxGBasis &gIn, const SxArray<PsiGI> &refOrb)
          : SxProjMatrix<Complex16> (gIn.ng, 64, refOrb.getSize),
            g(gIn), referenceOrbitals(refOrb)
        {
           // empty
        }
        \endcode
      - add a virtual destructor
        \code
        virtual ~SxMyClass () {}
        \endcode
      - overload (at least) getProjector(int) and getFactor(int)
        \code
        virtual SxDiracVec<Complex16> getProjector(int) const;
        virtual SxComplex16 getFactor(int) const;
        \endcode
      - if appropiate, overload getProjector with specified destination
        \code
        virtual void getProjector(int i, SxDiracVec<T> *target) const;
  };
        \endcode

    \author C. Freysoldt, freyso@fhi-berlin.mpg.de */
template <class T>
class SxProjMatrix
{
   protected:
      /// Number of projectors
      int nProj;
      /// Block size for blocked algorithm
      int blocksize;
      /// Number of basis functions
      int nElements;
   public:
      /// \name Dimensions
      //@{
      /// Get number of elements (basis set size)
      inline int getNElements () const { return nElements; }
      /// Get number of projectors
      inline int getNProj () const { return nProj; }

      /** \brief Set the block size.
        \note If the blockSize exceeds the number of projectors, it
              is reduced to the number of projectors. If it is smaller
              than 1, it is set to 1 (unblocked algorithm).
        */
      inline void setBlockSize (int blockSize) {
         if (blockSize < 1)
            blocksize = 1;
         else if (blockSize > nProj)
            blocksize = nProj;
         else
            blocksize = blockSize;
      }
      //@}
   protected:
      /** \brief Constructor
        \param nProjectors    number of projectors
        \param blockSize      maximum blockSize
        \param nBasisElements size of the basis
        */
      SxProjMatrix (int nProjectors, int blockSize, int nBasisElements)
         : nProj(nProjectors),
           nElements(nBasisElements)
      {
         SX_CHECK (nProjectors > 0);
         SX_CHECK (nBasisElements > 0);
         setBlockSize (blockSize);
      }
   public:

      virtual ~SxProjMatrix () { };

      /// \name Using the projector matrix
      //@{
      /** \brief Apply matrix to right hand side vector or matrix
          \f[
          \sum_i |\phi_i> f_i <\phi_i|\psi_n>
          \f]
          */
      inline SxDiracVec<T> apply (const SxDiracVec<T> &right) const
      {
         SX_CHECK (right.nRows () == nElements, right.nRows (), nElements);
         int nCols = right.nCols ();
         if (nCols == 0) nCols = 1;
         SxDiracVec<T> result(right.getSize (), 0.);
         result.reshape (nElements, nCols);
         if (right.nCols () > 1)
            applyAndAddMultiple(&result, right);
         else
            applyAndAddSingle(&result, right);
         return result;
      }
      /** \brief Apply matrix to right hand side vector or matrix and add 
                 to left
           \note This is equivalent to
           \code
*left += projMatrix ^ right;
           \endcode
           but slightly more efficient (avoids one vector addition).
           \note left and right MUST be different.
        */
      inline void applyAndAdd (      SxDiracVec<T> *left,
                               const SxDiracVec<T> &right) const
      {
         SX_CHECK (right.nRows () == nElements, right.nRows (), nElements);
         int nCols = (int)right.nCols ();
         if (nCols == 0) nCols = 1;
         if (right.nCols () > 1)
            applyAndAddMultiple(left, right);
         else
            applyAndAddSingle(left, right);
      }

   protected:
      /** Low-level kernel for doing projections & forming matrix elements
        */
      void projectionKernel (SxDiracVec<T> *res,
                             const SxDiracVec<T> &left,
                             const SxDiracVec<T> &right) const;
   public:

      /** \brief Get matrix elements
          \param left left hand side L, a (nElements x nL) matrix or vector
                      (nL=1)
          \param right left hand side R, a (nElements x nR) matrix or vector
                      (nR=1)
          \return a (nL x nR) matrix
          This calculates the matrix elements of the projector matrix
          \f[
          \sum_i <\psi_n|\phi_i> f_i <\phi_i|\psi_m>
          \f]
          \note left and right may be the same and this is acknowledged
                (projections are calculated only once)
          \note It is by far more efficient to call this routine than
                to work on the right hand side separately.
          \code
matrixElement = left.adjoint () ^ (projMatrix ^ right);    // slow
matrixElement = projMatrix.getMatrixElements(left, right); // fast
          \endcode
        */
      SxDiracVec<T> getMatrixElements (const SxDiracVec<T> &left,
                                       const SxDiracVec<T> &right) const
      {
         SxDiracVec<T> result;
         projectionKernel (&result, left, right);
         return result;
      }


      /** \brief Apply matrix to right hand side vector or matrix
          \f[
          \sum_i |\phi_i> f_i <\phi_i|\psi_n>
          \f]
       */
      inline SxDiracVec<T> operator^ (const SxDiracVec<T> &right) const
      {
         return apply (right);
      }

      /// Apply matrix to right hand side vector and add to left vector
      void applyAndAddSingle   (      SxDiracVec<T> *left,
                                const SxDiracVec<T> &right) const;
      /// Apply matrix to right hand side matrix and add to left matrix
      void applyAndAddMultiple (      SxDiracVec<T> *left,
                                const SxDiracVec<T> &right) const;
      //@}
   protected:
      /// An alias for SxDiracVec<T>
      typedef SxDiracVec<T> VecType;
      /// An alias for T::Type
      typedef typename T::Type CoeffType;

      /// \name Functions to be implemented by derived classes
      //@{
      /** \brief Get projector
          \note
          If the target-version of this function is implemented in
          the derived class, this function is never called. You can
          then implement a pseudo-function, which SX_EXITs if called by mistake,
          with the HAS_TARGET_GETPROJECTOR macro.
        */
      virtual SxDiracVec<T> getProjector (int i) const = 0;

      /// Write projector to target
      virtual void getProjector(int i, SxDiracVec<T> *target) const;
      /// Get inner factor
      virtual typename T::Type getFactor (int i) const = 0;
      /// Whether to use the factors
      virtual bool useFactors () const { return true; }

      /** Implements empty SxDiracVec<T> getProjector(int) const function
        */
#define HAS_TARGET_GETPROJECTOR \
      virtual VecType getProjector(int) const \
      { SX_EXIT; }

      /** Implements getFactor(int) routine that always returns 1
        */
#define NO_GETFACTOR \
      virtual CoeffType getFactor(int) const { return 1.; } \
      virtual bool useFactors () const { return false; }
      //@}

      
   public:
      /** \brief This routine is called before the projections start
          \param right hand side of the coming projection
          \return Number of right hand side vectors
        */
      virtual int prepare (const SxDiracVec<T> &right) const
      {
         int nCols = (int)right.nCols ();
         return (nCols == 0 && right.getSize () == nElements ? 1 : nCols);
      }

      /// \brief This routine is called after the projections are finished
      virtual void finalize () const { /* empty */ }
      

      /** \brief Get projections for a defined block
        @param iStart    first projector index
        @param iEnd      last projector index
        @param projBlock projectors (first...last)
        @param right     this is to be projected
        @return projections
        \note The resulting vector elements will be changed!
        */
      virtual SxDiracVec<T> 
      getProjectedSingle (int iStart, int iEnd,
                          const SxDiracVec<T> &projBlock,
                          const SxDiracVec<T> &right) const;
      /** \brief Get projections for a single projector
        @param iProj     projector index
        @param projector projector
        @param right     this is to be projected
        @return projection
        */
      virtual typename T::Type
      getProjectedSingle (int iProj,
                          const SxDiracVec<T> &projector,
                          const SxDiracVec<T> &right) const;
      /** \brief Get projections for a defined block
        @param iStart    first projector index
        @param iEnd      last projector index
        @param projBlock projectors (first...last)
        @param right     this is to be projected
        @return projections
        \note The resulting vector elements will be changed!
        */
      virtual SxDiracVec<T> 
      getProjectedMultiple (int iStart, int iEnd,
                            const SxDiracVec<T> &projBlock,
                            const SxDiracVec<T> &right) const;
      /** \brief Get projections for a single projector
        @param iProj     projector index
        @param projector projector
        @param right     this is to be projected
        @return projection
        */
      virtual SxDiracVec<T> 
      getProjectedMultiple (int iProj,
                            const SxDiracVec<T> &projector,
                            const SxDiracVec<T> &right) const;
      class SaveProjections;
};

      /** \brief Projector matrix with internal projection caching
      */
template <class T>
class SxProjMatrix<T>::SaveProjections : public SxProjMatrix<T>
{
   public:
      /// Whether the projection should be cached
      bool useCache;
      /// Clear the cache
      void clearCache ();
   protected:
      /// Whether something is cached
      mutable bool cacheUsed;
      /// The cached projections
      mutable SxDiracVec<T> savedProjections;
      /// Constructor
      SaveProjections (int nProjectors,
                       int blockSize,
                       int nBasisElements)
         : SxProjMatrix<T> (nProjectors, blockSize, nBasisElements),
           useCache(true),
           cacheUsed(false)
      {
         // empty
      }
      virtual SxDiracVec<T> getProjector (int i) const = 0;
      virtual typename T::Type getFactor (int i) const = 0;
   public:
      /// Prepare projections
      virtual int prepare (const SxDiracVec<T> &) const;
      /// Finalize projections
      virtual void finalize () const;
      /** \brief Get projections for a defined block
        @param iStart    first projector index
        @param iEnd      last projector index
        @param projBlock projectors (first...last)
        @param right     this is to be projected
        @return projections
        \note The resulting vector elements will be changed!
        */
      virtual SxDiracVec<T>
      getProjectedSingle (int iStart, int iEnd,
                          const SxDiracVec<T> &projBlock,
                          const SxDiracVec<T> &right) const;
      /** \brief Get projections for a single projector
        @param iProj     projector index
        @param projector projector
        @param right     this is to be projected
        @return projection
        */
      virtual typename T::Type
      getProjectedSingle (int iProj,
                          const SxDiracVec<T> &projector,
                          const SxDiracVec<T> &right) const;

      /** \brief Get projections for a defined block
        @param iStart    first projector index
        @param iEnd      last projector index
        @param projBlock projectors (first...last)
        @param right     this is to be projected
        @return projections
        \note The resulting vector elements will be changed!
        */
      virtual SxDiracVec<T>
      getProjectedMultiple (int iStart, int iEnd,
                            const SxDiracVec<T> &projBlock,
                            const SxDiracVec<T> &right) const;
      /** \brief Get projections for a single projector
        @param iProj     projector index
        @param projector projector
        @param right     this is to be projected
        @return projection
        */
      virtual SxDiracVec<T>
      getProjectedMultiple (int iProj,
                            const SxDiracVec<T> &projector,
                            const SxDiracVec<T> &right) const;
      /** \brief Get projections only.
          This calculates the projections of the right hand side \f$\psi_n\f$
          onto the projectors, i.e.
          \f[
          <\phi_i|\psi_n>
          \f]
          \param right a (nElements x nColsRight) matrix or a vector (nCols=1)
          \return a (nProj x nColsRight) matrix
        */
      SxDiracVec<T> getProjection(const SxDiracVec<T> &right);
      /** \brief Get projections only from extended vector.
          This calculates the projections of the right hand side \f$\psi_n\f$
          onto the projectors, i.e.
          \f[
          <\phi_i|\psi_n>
          \f]
          \param right a (nExtended x nColsRight) matrix or a vector (nCols=1)
                       where nExtended >= nElements
          \return a (nProj x nColsRight) matrix
        */
      SxDiracVec<T> getProjectionFromExtended(const SxDiracVec<T> &right);

      /** \brief Apply the matrix with precalculated projections
        \note This routine does multiply with the inner factor!
        */
      SxDiracVec<T> gradient (const SxDiracVec<T> &proj);

      virtual ~SaveProjections () {/* empty */}
};


#include <SxProjMatrix.hpp>

#endif /* _SX_PROJ_MATRIX_H_ */
