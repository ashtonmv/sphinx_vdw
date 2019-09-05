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

#ifndef _SX_DENSITY_H_
#define _SX_DENSITY_H_
#include <SxPtr.h>
#include <SxString.h>
#include <SxDFT.h>
class SxBinIO;

/** \brief Generic density

    This abstract class defines the interface for densities, which can be
    used in the general Hamiltonian solver. The actual implementation is done
    in a derived class, e.g. SxPAWRho.

    This is what you can do with a density:
    -# add, subtract, and scale. 
      You must overload these virtual functions:
      - operator+=, operator-=
      - plus_assign_ax
    -# square norm and scalar product of two densities (for Pulay mixer)
      Virtual functions to overload:
      - operator|
      - optional: normSqr (); (x.normSqr () defaults to (x|x))
    -# copy.
      Virtual functions
      - operator =
      - getCopy
    -# enforce number of electrons.
      Virtual functions:
      - renormalize

    \note Densities are assumed to be real-valued, i.e.,
    no complex scaling or complex scalar products are supported.

    \note
    The computation of rho is not interfaced here, because
    (a) the computation ingredients may vary, and
    (b) updating a density object is usually done INSIDE a Hamiltonian object
        (where you already know the derived type), and
    (c) computing a new density should work without needing an old one
        as auxiliary.
    Densities must be computed by the specific Hamiltonian (which may in turn
    call a computation routine of the specific density).

    \author C. Freysoldt freysoldt@mpie.de */
class SX_EXPORT_DFT SxDensity
{
   private:
      /** \brief Pointer for indirect densities
        Some routines create new densities. These must be of the derived
        type, but this violates virtual function signatures. The solution
        is to create these new densities as a SxPtr. To keep the interface
        intuitive, this indirection is wrapped by the SxDensity class.
        \example
        \code
SxDensity SxMyDensity::operator- (const SxDensity &x)
{
   const SxMyDensity &xRef = x.getRef<SxMyDensity> ();
   SxPtr<SxMyDensity> resPtr = SxPtr<SxMyDensity>::create (...);
   // compute the difference with xRef and resPtr
   ...
   return SxDensity (resPtr);
}
        \endcode
        */
      SxPtr<SxDensity> ptr;

   public:
      /// Constructor
      SxDensity () { /* empty */ }

      /// Constructor for indirect densities
      explicit SxDensity (const SxPtr<SxDensity> &inPtr);

      /// Destructor
      virtual ~SxDensity () { /* empty */ }

      /// Get reference to derived class (if possible)
      template<class T>
      T & getRef ()  {
         T * res = dynamic_cast<T *> (ptr ? ptr.getPtr () : this);
         SX_CHECK (res);
         return *res;
      }

      /// Get reference to derived class (if possible)
      template<class T>
      const T & getRef () const {
         const T * res = dynamic_cast<const T *> (ptr ? ptr.getPtr () : this);
         SX_CHECK (res);
         return *res;
      }

      /// Check if cast to derived class is possible
      template<class T>
      bool checkType () const {
         const T * res = dynamic_cast<const T *> (ptr ? ptr.getPtr () : this);
         return bool(res);
      }

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

      /** \brief Scalar product
        \note The prefactor convention is such that this should
        be the (volume) average, i.e.,
        \f[
        <\rho_1|\rho_2> = \frac{1}{\Omega}
                          \int_{\Omega} d^3 r~ \rho_1(r)\rho_2(r)
        \f]
        This yields a scale-independent, transferable measure.
        */
      virtual double operator| (const SxDensity &x) const;

      /** \brief Square norm
          \note The result must be the same as the scalar product with
          itself, i.e.
          \f$
          || \rho ||^2 = (\rho|\rho)
          \f$
        \note The prefactor convention is such that this should
        be the (volume) average, i.e.,
        \f[
        ||\rho_1|\rho_2||^2 = \frac{1}{\Omega}
                              \int_{\Omega} d^3 r~ |\rho_1(r)|^2
        \f]
        This yields a scale-independent, transferable measure.
        */
      virtual double normSqr () const;
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

      /// Renormalize
      virtual void renormalize ();

      /// Write density to file
      virtual void writeRho (const SxString &) const;

      /// Write density to file
      virtual void writeRho (SxBinIO &) const;

      /// Read density from file
      virtual void readRho (const SxString &);

      /// Read density from file
      virtual void readRho (const SxBinIO &file);

      /// Synchronize across MPI tasks
      virtual void syncMPI ();

};

#endif /* _SX_DENSITY_H_ */
