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

#include <SxBlockDensityMatrix.h>
#include <SxBinIO.h>
#include <SxLoopMPI.h>

/// Assign from a density
void SxBlockDensityMatrix::operator= (const SxDensity &in)
{
   SX_CHECK (in.checkType<SxBlockDensityMatrix> ());
   rho = in.getRef<SxBlockDensityMatrix> ().rho;
}

/// Add a density
void SxBlockDensityMatrix::operator+= (const SxDensity &x)
{
   SX_CHECK (x.checkType<SxBlockDensityMatrix> ());
   const SxBlockDensityMatrix& xRef = x.getRef<SxBlockDensityMatrix> ();
   SX_CHECK (getNSpin () == xRef.getNSpin (),
             getNSpin (), xRef.getNSpin ());
   SX_CHECK (getNSite () == xRef.getNSite (),
             getNSite (), xRef.getNSite ());
   for (int iSpin = 0; iSpin < getNSpin (); ++iSpin)
      for (int iSite = 0; iSite < getNSite (); ++iSite)
         rho(iSpin, iSite) += xRef.rho(iSpin, iSite);
}

/// Subtract a density
void SxBlockDensityMatrix::operator-= (const SxDensity &x)
{
   SX_CHECK (x.checkType<SxBlockDensityMatrix> ());
   const SxBlockDensityMatrix& xRef = x.getRef<SxBlockDensityMatrix> ();
   SX_CHECK (getNSpin () == xRef.getNSpin (),
             getNSpin (), xRef.getNSpin ());
   SX_CHECK (getNSite () == xRef.getNSite (),
             getNSite (), xRef.getNSite ());
   for (int iSpin = 0; iSpin < getNSpin (); ++iSpin)
      for (int iSite = 0; iSite < getNSite (); ++iSite)
         rho(iSpin, iSite) -= xRef.rho(iSpin, iSite);
}

/// axpy-like operation
void SxBlockDensityMatrix::plus_assign_ax (double a, const SxDensity &x)
{
   SX_CHECK (x.checkType<SxBlockDensityMatrix> ());
   const SxBlockDensityMatrix& xRef = x.getRef<SxBlockDensityMatrix> ();
   SX_CHECK (getNSpin () == xRef.getNSpin (),
             getNSpin (), xRef.getNSpin ());
   SX_CHECK (getNSite () == xRef.getNSite (),
             getNSite (), xRef.getNSite ());
   for (int iSpin = 0; iSpin < getNSpin (); ++iSpin)
      for (int iSite = 0; iSite < getNSite (); ++iSite)
         rho(iSpin, iSite).plus_assign_ax(a, xRef.rho(iSpin, iSite));
}

/// axpy-like operation for the spin density
void SxBlockDensityMatrix::plus_assign_aspin (double a, const SxDensity &x)
{
   SX_CHECK (x.checkType<SxBlockDensityMatrix> ());
   const SxBlockDensityMatrix& xRef = x.getRef<SxBlockDensityMatrix> ();
   SX_CHECK (getNSpin () == 2, getNSpin ());
   SX_CHECK (xRef.getNSpin () == 1, xRef.getNSpin ());
   SX_CHECK (getNSite () == xRef.getNSite (),
             getNSite (), xRef.getNSite ());
   for (int iSite = 0; iSite < getNSite (); ++iSite)  {
      rho(0, iSite).plus_assign_ax ( 0.5 * a, xRef.rho(0,iSite));
      rho(1, iSite).plus_assign_ax (-0.5 * a, xRef.rho(0,iSite));
   }
}


/** \brief Difference of two densities
  \note The return value must be an indirect density 
  with a pointer to the specific density type.
  */
SxDensity SxBlockDensityMatrix::operator- (const SxDensity &x) const
{
   SX_CHECK (x.checkType<SxBlockDensityMatrix> ());
   const SxBlockDensityMatrix& xRef = x.getRef<SxBlockDensityMatrix> ();
   SxPtr<SxBlockDensityMatrix> res = res.create ();
   res->computeDiff (*this, xRef);
   return SxDensity(res);
}

/** Get a copy (as a pointer)
  \note The return value must be an indirect density with a pointer to
  the specific density type.
*/
SxDensity SxBlockDensityMatrix::getCopy () const
{
   SxPtr<SxBlockDensityMatrix> res = res.create ();
   res->rho.reformat (getNSpin (), getNSite ());
   for (int iSpin = 0; iSpin < getNSpin (); ++iSpin)
      for (int iSite = 0; iSite < getNSite (); ++iSite)
         res->rho(iSpin, iSite).copy (rho(iSpin, iSite));
   return SxDensity(res);
}

/** Get spin density (as a pointer)
  \note The return value must be an indirect density with a pointer to
  the specific density type.
*/
SxDensity SxBlockDensityMatrix::spin () const
{
   SX_CHECK (getNSpin () == 2, getNSpin ());
   SxPtr<SxBlockDensityMatrix> res = res.create ();
   res->rho.reformat (1, getNSite ());
   for (int iSite = 0; iSite < getNSite (); iSite++)
      res->rho(0,iSite) = rho(0,iSite) - rho(1, iSite);
   return SxDensity(res);
}

/// Check if this is a spin-polarized density
bool SxBlockDensityMatrix::hasSpin () const
{
   return getNSpin () > 1;
}

/// Write density to file
//void SxBlockDensityMatrix::writeRho (const SxBinIO &) const;

/// Read density from file
//void SxBlockDensityMatrix::readRho (const SxBinIO &file);

/// Synchronize across MPI tasks
void SxBlockDensityMatrix::syncMPI ()
{
#ifdef USE_LOOPMPI
   ssize_t nElem = 0;
   // rho as a 1D flattened array
   SxArray<SxMatrix<Double> >& rho1 = rho;
   for (int i = 0; i < rho1.getSize (); ++i)
      nElem += rho1(i).getSize ();
   if (nElem > 1024 * 1024)  {
      cout << "Warning: large MPI sync in " << __FILE__ 
           << ":" << __LINE__ << endl;
      // better algorithm: split data into blocks to overlap
      // in-mem copying with inter-process communication
   }
   SxVector<Double> data(nElem);
   SX_MPI_MASTER_ONLY {
      ssize_t iData = 0;
      for (int i = 0; i < rho1.getSize (); ++i)  {
         size_t n = (size_t)rho1(i).getSize ();
         memcpy (data.elements + iData, rho1(i).elements, n * sizeof(double));
         iData += n;
      }
      SX_CHECK (iData == nElem);
   }
   SxLoopMPI::bcast(data, 0);
   SX_MPI_SLAVE_ONLY {
      size_t iData = 0;
      for (int i = 0; i < rho1.getSize (); ++i)  {
         size_t n = (size_t)rho1(i).getSize ();
         memcpy (rho1(i).elements, data.elements + iData, n * sizeof(double));
         iData += n;
      }
      SX_CHECK (iData == nElem);
   }
#endif
}

void SxBlockDensityMatrix::computeDiff (const SxBlockDensityMatrix &x, 
                                        const SxBlockDensityMatrix &y)
{
   SX_CHECK (x.getNSpin () == y.getNSpin (), x.getNSpin (), y.getNSpin ());
   SX_CHECK (x.getNSite () == y.getNSite (), x.getNSite (), y.getNSite ());
   SX_CHECK (getNSpin () == x.getNSpin (), getNSpin (), x.getNSpin ());
   SX_CHECK (getNSite () == x.getNSite (), getNSite (), x.getNSite ());
   SX_LOOP2(iSpin, iSite)
      rho(iSpin, iSite) = x.rho(iSpin, iSite) - y.rho(iSpin, iSite);
}

void SxBlockDensityMatrix::copyRho (const SxBlockDensityMatrix &x)
{
   SX_CHECK (getNSpin () == 0 || getNSpin () == x.getNSpin (),
             getNSpin (), x.getNSpin ());
   SX_CHECK (getNSite () == 0 || getNSite () == x.getNSite (),
             getNSite (), x.getNSite ());
   if (rho.getSize () == 0)
      rho.reformat (x.getNSpin (), x.getNSite ());
   SX_LOOP2(iSpin, iSite)
      rho(iSpin, iSite).copy (x.rho(iSpin, iSite));
}

void SxBlockDensityMatrix::readRho (const SxBinIO &io, const SxString &varName,
                                    int *offset)
{
   int offsetStart = *offset;
   try {
      for (int iSpin = 0; iSpin < getNSpin (); ++iSpin)  {
         *offset = offsetStart;
         for (int iSite = 0; iSite < getNSite (); ++iSite)  {
            io.readRow (varName, &rho(iSpin, iSite), iSpin, *offset);
            *offset += (int)rho(iSpin,iSite).getSize ();
         }
      }
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
}

void SxBlockDensityMatrix::writeRho (SxBinIO &io, const SxString &varName,
                                     int *offset) const
{
   int offsetStart = *offset;
   try {
      for (int iSpin = 0; iSpin < getNSpin (); ++iSpin)  {
         *offset = offsetStart;
         for (int iSite = 0; iSite < getNSite (); ++iSite)  {
            io.writeRow (varName, rho(iSpin, iSite), iSpin, *offset);
            *offset += (int)rho(iSpin,iSite).getSize ();
         }
      }
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
}

SxDiracVec<Complex16>
SxBlockDensityMatrix::operator^ (const SxDiracVec<Complex16> &p) const
{
   SX_CHECK (p.handle);
   int iSpin = p.handle->auxData.iSpin;
   SxDiracVec<Complex16> res(p.getSize ());
   res.reshape (p.nRows (), p.nCols ());
   res.handle->auxData = p.handle->auxData;
   ssize_t nStates = p.nCols ();
   int offset = 0;
   for (int iSite = 0; iSite < getNSite (); ++iSite)  {
      const SxMatrix<Double> &Aij = rho(iSpin, iSite);
      int npl = (int)Aij.nCols ();
      SX_CHECK (Aij.nCols () == Aij.nRows (),
                Aij.nCols (), Aij.nRows ());
      for (int iState = 0; iState < nStates; ++iState)  {
         // possible improvement: introduce inner loop of few sites here
         // such that larger segments of res/p are processed
         // criteria: all inner loop Aij should fit into L1 cache
         for (int ipl = 0; ipl < npl; ++ipl)  {
            SxComplex16 sum = 0.;
            for (int jpl = 0; jpl < npl; ++jpl)
               sum += Aij(ipl, jpl) * p(offset + jpl, iState);
            res(offset + ipl, iState) = sum;
         }
      }
      offset += npl;
   }
   SX_CHECK (offset == p.nRows (), offset, p.nRows ());
   res.setBasis (p.getBasisPtr ());
   return res;
}

void SxBlockDensityMatrix::sumMPI ()
{
#ifdef USE_LOOPMPI
   for (int iSpin = 0; iSpin < getNSpin (); ++iSpin)
      for (int iSite = 0; iSite < getNSite (); ++iSite)
         SxLoopMPI::sum (rho(iSpin, iSite));
#endif
}
