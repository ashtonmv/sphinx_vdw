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

#ifndef _SX_PAW_RHO_H_
#define _SX_PAW_RHO_H_

#include <SxDFT.h>
#include <SxDensity.h>
#include <SxRho.h>
#include <SxRadMat.h>
#include <SxPAWPot.h>

/** \brief PAW density


    \author C. Freysoldt freysoldt@mpie.de */
class SX_EXPORT_DFT SxPAWRho : public SxDensity
{
   public:
      /// Constructor
      SxPAWRho () { /* empty */ }

      /// Constructor
      SxPAWRho (const SxConstPtr<SxPAWPot> &pawPot);

      /// Destructor
      virtual ~SxPAWRho () {/* empty */}

      /// Assign from a density
      virtual void operator= (const SxDensity &);

      /// Assign from SxPAWRho
      void operator= (const SxPAWRho &);

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

      /// Scalar product
      virtual double operator| (const SxDensity &x) const;

      /** \brief Square norm
          \note The result must be the same as the scalar product with
          itself, i.e.
          \f$
          || \rho ||^2 = (\rho|\rho)
          \f$
        */
      virtual double normSqr () const;
      ///@}

      /** \brief Difference of two densities
        \note The return value must be an indirect density 
        with a pointer to the specific density type.
        */
      virtual SxDensity operator- (const SxDensity &x) const;

      /// Get a copy (as a pointer)
      virtual SxDensity getCopy () const;

      /// Get spin density (as a pointer)
      virtual SxDensity spin () const;

      /// Check if this is a spin-polarized density
      virtual bool hasSpin () const;

      /// Renormalize
      virtual void renormalize ();

      /// Read density from file
      virtual void readRho (const SxString &);

      /// Write density to file
      virtual void writeRho (const SxString &) const;

      /// Synchronize across MPI tasks
      virtual void syncMPI ();

      /// Plane-wave part of density
      SxRho pwRho;

      /// 1-center part of density
      SxRadMat Dij;

      /// Additional AO density matrices
      SxArray<SxPtr<SxBlockDensityMatrix> > blockAO;

      /// Total pseudo-core density (as contained in pwRho)
      //SxMeshR psCore;

      /// Get number of spin channels
      inline int getNSpin () const {
         SX_CHECK (pwRho.getNSpin () == Dij.getNSpin (),
                   pwRho.getNSpin (), Dij.getNSpin ());
         return pwRho.getNSpin ();
      }

      /// Pointer to PAW potential
      SxConstPtr<SxPAWPot> potPtr;

      /// Compute number of electrons
      double getNorm () const;

      /// Compute spin density
      double getSpin () const;
      /// Compute spin moments
      SxArray<double> getSpinMom (const SxAtomicStructure    &str) const;

      /// Set up density from atomic orbitals
      void atomicChargeDensity (const SxAtomicStructure    &str,
                                const SxConstPtr<SxPAWPot> &potPtrIn,
                                const SxRBasis             &R,
                                      int                   nSpin);
      void atomicChargeDensity (const SxAtomicStructure    &str,
                                    const SxConstPtr<SxPAWPot> &potPtrIn,
                                    const SxRBasis             &R,
                                    SxArray<SxArray<double> >  &atomSpin);

      void atomicChargeDensity (const SxAtomicStructure    &str,
                                const SxConstPtr<SxPAWPot> &potPtrIn,
                                const SxRBasis             &R,
                                const SxArray2<SxVector<Double> > &focc);
};

#endif /* _SX_PAW_RHO_H_ */
