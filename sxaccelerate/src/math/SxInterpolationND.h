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


#ifndef _INTERPOLATIONND_H_
#define _INTERPOLATIONND_H_

// Including the header-files
// -- SPHInX-Headers
#include <SxMatrix.h>
#include <SxMath.h>
#include <SxVector.h>
#include <SxArray.h>
#include <SxInterpolation1D.h>

// -- C++-Standardheaders

/** \brief ...

    \b SxClass = SPHInX ...

    This class is used to solve n-dimensional interpolation problems.
    It weights the values of the known grid points that surround the 
    requested point and thus creates a reasonable (not very exact) 
    assumption for the searched value. If you desired a better 
    approximation for unknown points, you should use an interpolation
    that fits with higher-order polynomials (not just linear ones).
    
    \author Thomas Uchdorf, t.uchdorf@mpie.de */
class SX_EXPORT_MATH SxBilinearInterpolationND
{

   public:

      // Constructors
      // -- Constructor
      /** \brief Creates a BilinearInterpolation object */
      SxBilinearInterpolationND ();

      // -- Destructor
      /** \brief Destroys a BilinearInterpolation object */
      ~SxBilinearInterpolationND ();

      // Methods
      /** \brief Returns the assumed value at the requested position x0 */
      double computeBilinearInterpolationND (
         const SxVector<Double> &x0,
         const SxVector<Double> &values,
         const SxVector<Int> &dimension);

      /** \brief Returns recursivly the index of the grid vector
        (for example [0, 0, 0] -> 0, [0, 0, 1] -> 1, ...,
        [0, 0, zDim - 1] -> zDim - 1, [0, 1, 0] -> zDim, ..., 
        [0, yDim - 1, zDim - 1] -> zDim - 1 + (zDim * yDim - 1),
        ...)*/      
      int getIdx (const SxVector<Int> &vec, int i);

      /** \brief Returns the grid vector for the corresponding index
      (for example 0 -> [0, 0, 0] or zDim -> [0, 1, 0]) */
      SxVector<Int> getVector (int idx);

      /** \brief Computes the differences of indices' needed to make
        one step in each coordinate direction */
      void computeStep ();

      /** \brief Computes the weighted addend that helps determining
        the requested value */
      double computeAddend ();

      /** \brief Calls computeAddend () for all surrounding positions */
      void genericForLoop (int iDim);
      
   protected:

      // Methods
      
      // Members
      /** \brief Vector that saves the corresponding values for all grid
        positions */
      SxVector<Double> f;

      /** \brief Contains pieces of information on the dimensions of the
        grid */
      SxVector<Int> dim;

      /** \brief Reflects the step widths needed to move in each direction
        (minimizes the number of flops)  */ 
      SxVector<Int> idxStep;

      /** \brief Specifies the surrounding points (for example 2-dimensional: 
        00 x-----x 11
           |     |
           |     |
           |     |
        01 x-----x 10
           )*/
      SxArray<bool> iCoord;

      /** \brief Describes the sum of the weighted values
        (assumption for the requested value)*/
      double sum;

      /** \brief Index at a certain surrounding point used to simplify 
        movement through the surrounding grid points 
        (for example 2-dimensional: index of the top-left position) */
      int defaultIdx;

      /** \brief Quotient of the differences between the requested position 
        and the surrounding points dimension specific minimum and the total
        length of the dimension ((x - min (dim)) / (max (dim) - min (dim))) */
      SxVector<Double> rel;

};

/** \brief ...

    \b SxClass = SPHInX ...

    This class is used to solve n-dimensional interpolation problems.
    It solves the problem by conducting multiple one-dimensional 
    interpolations with user specified order.
    
    \author Thomas Uchdorf, t.uchdorf@mpie.de */
class SxInterpolationND
{
   public:

      // Constructors
      // -- Constructor
      /** \brief Creates a InterpolationND object */
      SxInterpolationND ();
      

      /** \brief Creates a InterpolationND object */
      SxInterpolationND (const SxVector<Double> &values, const SxVector<Int> &dimension);      

      // -- Destructor
      /** \brief Destroys a InterpolationND object */
      ~SxInterpolationND ();

      // Methods
      /** \brief Returns recursivly the index of the grid vector
        (for example [0, 0, 0] -> 0, [0, 0, 1] -> 1, ...,
        [0, 0, zDim - 1] -> zDim - 1, [0, 1, 0] -> zDim, ..., 
        [0, yDim - 1, zDim - 1] -> zDim - 1 + (zDim * yDim - 1),
        ...)*/      
      int getIdx (const SxVector<Int> &vec, int i);

      /** \brief Returns the grid vector for the corresponding index
      (for example 0 -> [0, 0, 0] or zDim -> [0, 1, 0]) */
      SxVector<Int> getVector (int idx);

      /** \brief Returns the assumed value at the requested position x0.
        Here interpolationDim specifies the order of the interpolation in the
        various coordinate directions (for example [2, 3, 3] means square 
        interpolation in the [1, 0, 0] and cubic interpolation in the [0, 1, 0]
        and [0, 0, 1] direction) */
      double getVal (const SxVector<Double> &x0, const SxVector<Int> &interpolationDim);

   protected:

      // Methods
      /** \brief Function that works recursivly and returns the interpolated 
        value at position x0 */
      double rec (int iDim, const SxVector<Int> &foo,
            const SxVector<Int> &totalMin, const SxVector<Double> &x0, const SxVector<Int> &interpolationDim);
      
      // Members
      /** \brief Vector that saves the corresponding values for all grid
        positions */
      SxVector<Double> f;

      /** \brief Contains pieces of information on the dimensions of the
        grid */
      SxVector<Int> dim;
   
};


#endif /* _INTERPOLATIONND_H_ */
