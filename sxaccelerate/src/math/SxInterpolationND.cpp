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

// Including the header-file
#include <SxInterpolationND.h>

#undef max
#undef min

// Constructors
// -- Constructor
SxBilinearInterpolationND::SxBilinearInterpolationND ()
{
   
   // empty

}

// -- Destructor
SxBilinearInterpolationND::~SxBilinearInterpolationND ()
{

   // empty

}

/*
   3-dimensional index function
int SxBilinearInterpolationND::getIdx (const SxVector3<Int> &vec, const SxVector3<Int> &dim)
{
   return vec(2) + dim(2) * (vec(1) + dim(1) * vec(0));
}
*/

int SxBilinearInterpolationND::getIdx (const SxVector<Int> &vec, int i)
{
   
   SX_CHECK (i >= 0, i);
   
   if (i == 0)  {
      
      return vec(0);
      
   } else  {
      
      return vec(i) + dim(i) * getIdx (vec, i - 1);
      
   }
   
}

SxVector<Int> SxBilinearInterpolationND::getVector (int idx)
{

   SX_CHECK (idx >= 0 && idx < dim.product (), idx, dim.product ());

   SxVector<Int> ret(dim.getSize ());
  
   // Determining the first value for the vector by computing the product
   // from all dimension values without the first one
   int prod (dim.product ());
   prod /= dim(0); 
   ret(0) = idx / prod;

   // Calculating the values for the coordinates that lie in between
   int i (0);
   for (i = 0; i < dim.getSize () - 1; ++i)  {

      prod /= dim(i + 1);
      ret (i+1) = idx / prod - ret(i) * dim(i+1);

   }
   
   // Determining the last value for the vector with modulo
   // mathematical background:
   // idx mod b = vec(i) <=> vec(i) = idx - [idx / b] * b;
   // idx -= vec(i), i = nDim - 1 ,..., 0
   i = static_cast<int>(dim.getSize ()) - 1; 
   ret(i) = idx - (idx / dim(i)) * dim(i);
   
   return ret; 

}

/*
   3-dimensional vector function
SxVector<Int> SxBilinearInterpolationND::getVector (const SxVector<Int> &dim, int idx)
{
   SX_CHECK (idx >= 0 && idx < dim.product (), idx, dim.product ());
   SxVector<Int> ret (3);
   int coeff (idx / dim(2));
   
   ret(0) = coeff / dim(1);
   ret(1) = coeff - ret(0) * dim(1);
   ret(2) = idx - coeff * dim(2);
   
   return ret;

}
*/

void SxBilinearInterpolationND::computeStep ()
{

   idxStep.resize (dim.getSize ());
   idxStep(dim.getSize () - 1) = 1;
   int prod (1);
   for (int i = 1; i < dim.getSize (); ++i)  {

      prod *= dim(dim.getSize () - i);
      idxStep(dim.getSize () - 1 - i) = prod;

   }

}

double SxBilinearInterpolationND::computeAddend ()
{

   SxVector<Double> tmp (dim.getSize ());

   int currentIdx (defaultIdx);

   for (int iDim = 0; iDim < dim.getSize (); ++iDim)  {
      
      // Testing if the current coordinate resembles a minimum or maximum value
      if (iCoord(iDim) == (bool) 0)  {

         tmp(iDim) = 1. - rel(iDim);

      } else  {

         tmp(iDim) = rel(iDim);
         currentIdx += idxStep(iDim);

      }

   }

   return f(currentIdx) * tmp.product ();

}

void SxBilinearInterpolationND::genericForLoop (int iDim) 
{ 

   // loop over minima and maxima (neighbors)
   // 0 respectively false represents a minimum
   // 1 respectively true signs a maximum
   for (int i = 0; i < 2; ++i)  { 

#     ifdef WIN32
         // allow the following bool-cast
#        pragma warning(disable:4800)
#     endif
      iCoord(iDim) = (bool) i; 

      if (iDim < dim.getSize () - 1) 

         genericForLoop (iDim + 1); 

      else {

         sum += computeAddend ();

      }

   } 

} 

double SxBilinearInterpolationND::computeBilinearInterpolationND (
      const SxVector<Double> &x0,
      const SxVector<Double> &values,
      const SxVector<Int> &dimension)
{
   
   SX_CHECK (dimension.getSize () != 0 && x0.getSize () == dimension.getSize (), x0.getSize (), dimension.getSize ());
   
   for (int i = 0; i < dimension.getSize (); ++i) {

      SX_CHECK (x0(i) >= 0 && x0(i) < dimension(i), i, x0(i), dimension(i));

   }
  
   f = values;
   dim = dimension;
   
   rel.resize (dim.getSize ());
   SxVector<Int> min (dim.getSize ()), max (dim.getSize ());

   //	Searching the framing points
   for (int i =0; i < dim.getSize (); ++i)  {

      min(i) = (int) (x0(i));
      max(i) = (int) (x0(i)) + 1;
      rel(i) = 1. * (x0(i) - min(i)) /  (max(i) - min(i));

   }

   // Determining the default index (index of (min[0], ..., min[nDim - 1]))
   defaultIdx = getIdx (min, static_cast<int>(dim.getSize ()) - 1);
     
   // Computing the differences of indices' needed to make one step for each coordinate direction
   computeStep ();

   // Initializing the sum which resembles the return value with 0
   sum = 0.;

   // Initializing the bool array that specifies the surrounding points
   iCoord.resize (dim.getSize ());
   for (int iDim = 0; iDim < dim.getSize (); ++iDim)  {

      iCoord(iDim) = false;

   }

   // Calling a recursive procedure which represents a generic for loop and that
   // determines the assumed value for the requested position
   genericForLoop (0);
  
   return sum;   

}

SxInterpolationND::SxInterpolationND ()
{
   
   // empty

}

SxInterpolationND::SxInterpolationND (const SxVector<Double> &values, const SxVector<Int> &dimension) : f (values), dim (dimension)
{

   // empty
   
}

SxInterpolationND::~SxInterpolationND ()
{

   // empty
 
}

int SxInterpolationND::getIdx (const SxVector<Int> &vec, int i)
{
   
   SX_CHECK (i >= 0, i);
   
   if (i == 0)  {
      int ret (vec(0));
      SX_CHECK (ret >= 0, ret);      
      return ret;
      
   } else  {
      int ret (vec(i) + dim(i) * getIdx (vec, i - 1));
      SX_CHECK (ret >= 0, ret);            
      return ret;
      
   }
   
}

SxVector<Int> SxInterpolationND::getVector (int idx)
{

   SX_CHECK (idx >= 0 && idx < dim.product (), idx, dim.product ());

   SxVector<Int> ret(dim.getSize ());
  
   // Determining the first value for the vector by computing the product
   // from all dimension values without the first one
   int prod (dim.product ());
   prod /= dim(0); 
   ret(0) = idx / prod;

   // Calculating the values for the coordinates that lie in between
   int i (0);
   for (i = 0; i < dim.getSize () - 1; ++i)  {

      prod /= dim(i + 1);
      ret (i+1) = idx / prod - ret(i) * dim(i+1);

   }
   
   // Determining the last value for the vector with modulo
   // mathematical background:
   // idx mod b = vec(i) <=> vec(i) = idx - [idx / b] * b;
   // idx -= vec(i), i = nDim - 1 ,..., 0
   i = static_cast<int>(dim.getSize ()) - 1; 
   ret(i) = idx - (idx / dim(i)) * dim(i);
   
   return ret; 

}


double SxInterpolationND::getVal (const SxVector<Double> &x0, const SxVector<Int> &interpolationDim)
{
   // {(y, f (x0(0), y))} <- interpol1D (x, f(x, y[i]), x0(0)) for all y[i]
   
   // {((y, z), f (x0(0), y, z))} <- interpol1D (x, f(x, y[i], z[j]), x0(0)) for all y[i] and all z[j]
   
   SX_CHECK (dim.getSize () == interpolationDim.getSize (), dim.getSize (), interpolationDim.getSize ());
   for (int iDim = 0; iDim < dim.getSize (); ++iDim)  {

      SX_CHECK (interpolationDim(iDim) + 1 < dim(iDim), interpolationDim(iDim), dim(iDim));

   }

   SxVector<Int> min (dim.getSize ()), totalMin (dim.getSize ());

   //	Searching the framing points
   for (int i = 0; i < dim.getSize (); ++i)  {

      min(i) = (int) (x0(i));
      totalMin(i) = (int) (x0(i));

      // Decrementing the index of the total minimum if there are not enough
      // points available for the desired interpolation order
      while (totalMin(i) + interpolationDim(i) + 1 > dim(i))  {

         --totalMin(i);      

      }
      
   }
  
   SxVector<Int> foo (dim.getSize ());
   foo = 0;
   return rec ((int)dim.getSize () - 1, foo, totalMin, x0, interpolationDim);   
   
}

double SxInterpolationND::rec (int iDim, const SxVector<Int> &foo,
      const SxVector<Int> &totalMin, const SxVector<Double> &x0,
      const SxVector<Int> &interpolationDim)
{

   // Creating a vector that stores the interpolation results of the inner
   // loops
   SxVector<Double> tmp (interpolationDim (iDim) + 1);

   // Looping over all indices of the current dimension
   for (int i = 0; i < interpolationDim (iDim) + 1; ++i)  {

      // Modifying the current indices stored in foo respectively newFoo
      SxVector<Int> newFoo (foo.getSize ());
      newFoo.copy (foo);
      newFoo(iDim) = i;

      // Interpolating (x, f(x, y, z)) at (x0, y, z) if the first plane was reached
      if (iDim == 1) {

         SxVector<Double> x (interpolationDim(0) + 1), y (interpolationDim(0) + 1);
         
         for (int iCol = 0; iCol < interpolationDim(0) + 1; ++iCol)  {

            newFoo(0) = iCol;
            SxVector<Int> vec (dim.getSize ());      
            vec = totalMin + newFoo;
            x(iCol) = getIdx (vec, static_cast<int>(dim.getSize ()) - 1);
            y(iCol) = f((int)x(iCol));
            x(iCol) = totalMin(0) + iCol;

         }
         
         SxPolynomialInterpolation1D interpol (x, y, Lagrange); 
         
         tmp(i) = interpol.getVal (x0(0));

      }
     
      // Stepping one loop deeper and passing the current indices in the vector
      // newFoo
      if (iDim > 1) {

         tmp(i) = rec (iDim - 1, newFoo, totalMin, x0, interpolationDim);

      }

   }

   // Constructing a vector which contains the indices of the current dimension
   SxVector<Int> tmpSrc (interpolationDim(iDim) + 1);
   for (int iLin = 0; iLin < interpolationDim(iDim) + 1; ++iLin)  {

      tmpSrc(iLin) = totalMin(iDim) + iLin;

   }

   // Interpolating (y, interpolation (x, f (x, y, z), x0)) at y0 if the 
   // second plane was reached otherwise interpolating
   // (z, interpolation ([x,y], f(x, y, z), [x0, y0]) at z0
   SxPolynomialInterpolation1D finalInterpol (tmpSrc, tmp, Lagrange);

   // Returning the current interpolation result to the higher loop levels
   return finalInterpol.getVal (x0(iDim));

}
