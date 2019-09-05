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

#include <SxHessian.h>
#include <SxNeighbors.h>

SxHessian::SxHessian (const SxMatrix<Double> &full)
{
   int nAtoms = (int)full.nCols () / 3;
   SX_CHECK (full.nCols () == 3 * nAtoms, full.nCols (), nAtoms);
   SX_CHECK (full.nRows () == 3 * nAtoms, full.nRows (), nAtoms);

   hessian.resize (nAtoms);
   neighbors.resize (nAtoms);
   for (int ia = 0; ia < nAtoms; ++ia)  {
      hessian(ia).resize (nAtoms);
      neighbors.resize (nAtoms);
      for (int ja = 0; ja < nAtoms; ++ja)  {
         neighbors(ia)(ja) = ja;
         SxMatrix3<Double> &D = hessian(ia)(ja);
         for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
               D(i,j) = full(3 * ia + i, 3*ja + j);
      }
   }
}

SxMatrix<Double> SxHessian::getFull () const
{
   int nAtoms = (int)hessian.getSize ();
   SxMatrix<Double> res(3 * nAtoms, 3*nAtoms);
   res.set (0.);

   for (int ia = 0; ia < nAtoms; ++ia)  {
      for (int in = 0; in < nAtoms; ++in)  {
         int ja = neighbors(ia)(in);
         const SxMatrix3<Double> &D = hessian(ia)(ja);
         for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
               res(3 * ia + i, 3*ja + j) = D(i,j);
      }
   }
   return res;
}

void SxHessian::writeFull (SxBinIO &io,
                           const SxMatrix<Double> &full)
{
   SX_CHECK (full.nCols () == full.nRows (), full.nCols (), full.nRows ());
   int nDof = (int)full.nCols ();
   SX_CHECK (nDof > 0, nDof);
   try  {
      io.addDimension ("nDoF", nDof);
      io.write ("hessian", full, "nDoF", "nDoF");
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}

SxMatrix<Double> SxHessian::readFull (const SxBinIO &io)
{
   int nDoF = io.getDimension ("nDoF");
   SxMatrix<Double> res (nDoF, nDoF);
   try {
      io.read ("hessian", &res, nDoF, nDoF);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   return res;
}

void SxHessian::write (const SxString &fileName,
                       const SxMatrix<Double> &hessianFull,
                       const SxAtomicStructure &str)
{

   SX_CHECK (hessianFull.nCols () == hessianFull.nRows (),
             hessianFull.nCols (), hessianFull.nRows ());
   SX_CHECK (3 * str.getNAtoms () == hessianFull.nCols (),
             str.getNAtoms (), hessianFull.nCols ());
   try {
      SxBinIO io(fileName, SxBinIO::BINARY_WRITE_ONLY);
      writeFull (io, hessianFull);
      str.write (io);
      io.setMode (SxBinIO::WRITE_DATA);
      writeFull (io, hessianFull);
      str.write (io);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}

