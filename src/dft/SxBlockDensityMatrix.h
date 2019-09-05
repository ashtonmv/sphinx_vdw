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

#ifndef _SX_BLOCK_DENSITY_MATRIX_H_
#define _SX_BLOCK_DENSITY_MATRIX_H_

#include <SxDFT.h>
#include <SxDensity.h>
#include <SxNArray.h>
#include <SxMatrix.h>
#include <SxDirac.h>

/** \brief Container class for block density matrices

    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DFT SxBlockDensityMatrix : public SxDensity
{
   public:
      /// Density matrices for all sites (iSpin,iSite)
      SxArray2<SxMatrix<Double> > rho;

      /// Assignment operator
      void operator= (const SxBlockDensityMatrix &in)
      {
         if (&in == this) return;
         rho = in.rho;
      }

      /// Get number of spins
      int getNSpin () const { return int(rho.getDim(0)); }
      /// Get number of sites
      int getNSite () const { return int(rho.getDim(1)); }

      /// Resize
      void resize (int nSpin, int nSite)
      {
         rho.reformat (nSpin, nSite);
      }

      ///@{
      /// Get a 1-center matrix
      SxMatrix<Double> &operator() (int iSpin, int iSite)  {
         return rho(iSpin, iSite);
      }

      const SxMatrix<Double> &operator() (int iSpin, int iSite) const  {
         return rho(iSpin, iSite);
      }

      // autoloop wrapper
      template<class T1, class T2>
      SxMatrix<Double> &operator() (const T1 &iSpin, const T2 &iSite)
      {
         SxAutoLoop::setLimit(iSpin, getNSpin ());
         SxAutoLoop::setLimit(iSite, getNSite ());
         return rho((int)iSpin, (int)iSite);
      }

      // autoloop wrapper
      template<class T1, class T2>
      const SxMatrix<Double>& operator()(const T1 &iSpin, const T2 &iSite) const
      {
         SxAutoLoop::setLimit(iSpin, getNSpin ());
         SxAutoLoop::setLimit(iSite, getNSite ());
         return rho((int)iSpin, (int)iSite);
      }
      ///@}

      // --- SxDensity interface ---

      /// Assign from a density
      virtual void operator= (const SxDensity &);

      ///\name Vector-like operations
      ///@{
      /// Add a density
      virtual void operator+= (const SxDensity &x);

      /// Subtract a density
      virtual void operator-= (const SxDensity &x);

      /// axpy-like operation
      virtual void plus_assign_ax (double a, const SxDensity &x);

      /// axpy-like operation for the spin density
      virtual void plus_assign_aspin (double a, const SxDensity &x);

      // not overloaded here:
      //virtual double operator| (const SxDensity &x) const;

      // not overloaded here:
      //virtual double normSqr () const;
      ///@}

      /** \brief Difference of two densities
        \note The return value must be an indirect density 
        with a pointer to the specific density type.
        */
      virtual SxDensity operator- (const SxDensity &x) const;

      /** Get a copy (as a pointer)
        \note The return value must be an indirect density with a pointer to
        the specific density type.
      */
      virtual SxDensity getCopy () const;

      /** Get spin density (as a pointer)
        \note The return value must be an indirect density with a pointer to
        the specific density type.
      */
      virtual SxDensity spin () const;

      /// Check if this is a spin-polarized density
      virtual bool hasSpin () const;

      /// Synchronize across MPI tasks
      virtual void syncMPI ();

      /// Sum across MPI tasks
      void sumMPI ();

      /// Sets rho from x - y
      void computeDiff (const SxBlockDensityMatrix &x, 
                        const SxBlockDensityMatrix &y);
      
      /// Copy all density matrices (real copy)
      void copyRho (const SxBlockDensityMatrix &x);

      /** \brief Read rho
        @param io netcdf file
        @param varName  variable name to write in, must be dimensioned nSpin, nSite
        @param offset   pointer to site offset

        This routine reads the density matrix data. The matrices must have been properly
        dimensioned beforehand.
        */
      void readRho (const SxBinIO &io, const SxString &varName, int *offset);

      /** \brief Write rho
        @param io netcdf file
        @param varName  variable name to write in, must be dimensioned nSpin, nSite
        @param offset   pointer to site offset

        This routine writes the density matrix data into the predefined variable.
        */
      void writeRho (SxBinIO &io, const SxString &varName, int *offset) const;

      /** \brief Multiply matrices with projections
        @param p   \f$\langle p_i|\psi\rangle\f$

        \f[
        \sum_j A_{ij} <p_j | \psi>
        \f]

        Aij are the matrices stored in rho.
        */
      SxDiracVec<Complex16> operator^ (const SxDiracVec<Complex16> &p) const;

};

#endif /* _SX_BLOCK_DENSITY_MATRIX_H_ */
