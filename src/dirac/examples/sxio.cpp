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

#include <SxUtil.h>
#include <SxException.h>
#include <SxBinIO.h>


/** \example sxb.cpp

    This example demonstrates how cross-platform SFHIngX binary files (*.sxb)
    can be created and read in.

    \author Sixten Boeck
*/
int main ()
{
   int i, mode;

   // --- save a simple 1d-vector to disc
   SxDiracVec<Double> vec(10);
   for (i=0; i < vec.getSize(); i++)   vec(i) = sqrt( (double)i );

   // --- save 1d-vector to disc
   try  {
      SxBinIO out ("sxb-vec1d.sxb", SxBinIO::BINARY_WRITE_ONLY);
      for (mode=SxBinIO::WRITE_HEADER; mode != SxBinIO::WRITE_DATA; mode++)  {
         out.setMode (mode);
         out.addDimension ("vecSize", vec.getSize()); // (1) create dimension
         out.write ("myVector", vec, "vecSize");     // (2) write vector
      }
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }


   // --- create some 3d-data
   SxVector3<Int> dim (10, 20, 30);
   SxVector3<Double> t, center;
   SxMatrix3<Double> cell (10,0,0, 0,8,0, 0,0,15);
   SxDiracVec<Double> xyzData (dim.product());
   center = dim / 2;
   i=0;
   for (t(0)=0; t(0) < dim(0); t(0)++)  {
      for (t(1)=0; t(1) < dim(1); t(1)++)  {
         for (t(2) = 0; t(2) < dim(2); t(2)++)  {
            xyzData(i++) = sqrt ( (t - center).absSqr().sum() );
         }
      }
   }

   // --- save 3d-mesh to disc
   try  {
      SxBinIO out ("sxb-mesh3d.sxb", SxBinIO::BINARY_WRITE_ONLY);
      out.writeMesh (xyzData, cell, dim);  // (1) write header information
      out.setMode (SxBinIO::WRITE_DATA);      // (2) switch to data mode
      out.writeMesh (xyzData, cell, dim);  // (3) write data 
      out.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

}
