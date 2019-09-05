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

#ifndef _SX_BASIS_H_
#define _SX_BASIS_H_

#include <SxMemConsumer.h>
#include <SxDiracLib.h>
#include <SxString.h>
#include <SxArray.h>
#include <SxPrecision.h>

class SxGBasis;
class SxRBasis;
class SxRadBasis;
class SxAOBasis;
class SxPartialWaveBasis;
class SxPAWBasis;
class SxRadRBasis;
class SxRadGBasis;
template <class T> class SxDiracVec;

#define VIRTUAL_PROJECT_TO(B)                                                 \
      virtual void projectTo (const B *, const void *, void *,                \
                              const float &, const float &) const             \
      {  projectionFailed ( #B , "(R4->R4)"); }                               \
      virtual void projectTo (const B *, const void *, void *,                \
                              const double &, const float &) const            \
      {  projectionFailed ( #B , "(R8->R4)"); }                               \
      virtual void projectTo (const B *, const void *, void *,                \
                              const SxComplex8 &, const float &) const        \
      {  projectionFailed ( #B , "(C8->R4)"); }                               \
      virtual void projectTo (const B *, const void *, void *,                \
                              const SxComplex16 &, const float &) const       \
      {  projectionFailed ( #B , "(C16->R4)"); }                              \
                                                                              \
      virtual void projectTo (const B *, const void *, void *,                \
                              const float &, const double &) const            \
      {  projectionFailed ( #B , "(R4->R8)"); }                               \
      virtual void projectTo (const B *, const void *, void *,                \
                              const double &, const double &) const           \
      {  projectionFailed ( #B , "(R8->R8)"); }                               \
      virtual void projectTo (const B *, const void *, void *,                \
                              const SxComplex8 &, const double &) const       \
      {  projectionFailed ( #B , "(C8->R8)"); }                               \
      virtual void projectTo (const B *, const void *, void *,                \
                              const SxComplex16 &, const double &) const      \
      {  projectionFailed ( #B , "(C16->R8)"); }                              \
                                                                              \
      virtual void projectTo (const B *, const void *, void *,                \
                              const float &, const SxComplex8 &) const        \
      {  projectionFailed ( #B , "(R4->C8)"); }                               \
      virtual void projectTo (const B *, const void *, void *,                \
                              const double &, const SxComplex8 &) const       \
      {  projectionFailed ( #B , "(R8->C8)"); }                               \
      virtual void projectTo (const B *, const void *, void *,                \
                              const SxComplex8 &, const SxComplex8 &) const   \
      {  projectionFailed ( #B , "(C8->C8)"); }                               \
      virtual void projectTo (const B *, const void *, void *,                \
                              const SxComplex16 &, const SxComplex8 &) const  \
      {  projectionFailed ( #B , "(C16->C8)"); }                              \
                                                                              \
      virtual void projectTo (const B *, const void *, void *,                \
                              const float &, const SxComplex16 &) const       \
      {  projectionFailed ( #B , "(R4->C16)"); }                              \
      virtual void projectTo (const B *, const void *, void *,                \
                              const double &, const SxComplex16 &) const      \
      {  projectionFailed ( #B , "(R8->C16)"); }                              \
      virtual void projectTo (const B *, const void *, void *,                \
                              const SxComplex8 &, const SxComplex16 &) const  \
      {  projectionFailed ( #B , "(C8->C16)"); }                              \
      virtual void projectTo (const B *, const void *, void *,                \
                              const SxComplex16 &, const SxComplex16 &) const \
      {  projectionFailed ( #B , "(C16->C16)"); }                             \


#define REGISTER_PROJECTOR(FROM,TO,PROJECTOR)                                 \
     virtual void projectTo (const TO *b, const void *in, void *o,            \
                             const float &,const float &) const               \
     {                                                                        \
        const SxDiracVec<Float> &vec = *((const SxDiracVec<Float> *)in);      \
        SxDiracVec<Float> &res = *((SxDiracVec<Float> *)o);                   \
        res = PROJECTOR (b, vec);                                             \
     }                                                                        \
     virtual void projectTo (const TO *b, const void *in, void *o,            \
                             const float &, const double &) const             \
     {                                                                        \
        const SxDiracVec<Float>  &vec = *((const SxDiracVec<Float> *)in);     \
        SxDiracVec<Double> &res = *((SxDiracVec<Double> *)o);                 \
        res = PROJECTOR (b, vec);                                             \
     }                                                                        \
     virtual void projectTo (const TO *b, const void *in, void *o,            \
                             const float &, const SxComplex8 &) const         \
     {                                                                        \
        const SxDiracVec<Float>  &vec =*((const SxDiracVec<Float> *)in);      \
        SxDiracVec<Complex8> &res = *((SxDiracVec<Complex8> *)o);             \
        res = PROJECTOR (b, vec);                                             \
     }                                                                        \
     virtual void projectTo (const TO *b, const void *in, void *o,            \
                             const float &, const SxComplex16 &) const        \
     {                                                                        \
        const SxDiracVec<Float>  &vec = *((const SxDiracVec<Float> *)in);     \
        SxDiracVec<Complex16> &res = *((SxDiracVec<Complex16> *)o);           \
        res = PROJECTOR (b, vec);                                             \
     }                                                                        \
                                                                              \
     virtual void projectTo (const TO *b, const void *in, void *o,            \
                             const double &,const float &) const              \
     {                                                                        \
        const SxDiracVec<Double> &vec = *((const SxDiracVec<Double> *)in);    \
        SxDiracVec<Float> &res = *((SxDiracVec<Float> *)o);                   \
        res = PROJECTOR (b, vec);                                             \
     }                                                                        \
     virtual void projectTo (const TO *b, const void *in, void *o,            \
                             const double &, const double &) const            \
     {                                                                        \
        const SxDiracVec<Double>  &vec = *((const SxDiracVec<Double> *)in);   \
        SxDiracVec<Double> &res = *((SxDiracVec<Double> *)o);                 \
        res = PROJECTOR (b, vec);                                             \
     }                                                                        \
     virtual void projectTo (const TO *b, const void *in, void *o,            \
                             const double &, const SxComplex8 &) const        \
     {                                                                        \
        const SxDiracVec<Double>  &vec =*((const SxDiracVec<Double> *)in);    \
        SxDiracVec<Complex8> &res = *((SxDiracVec<Complex8> *)o);             \
        res = PROJECTOR (b, vec);                                             \
     }                                                                        \
     virtual void projectTo (const TO *b, const void *in, void *o,            \
                             const double &, const SxComplex16 &) const       \
     {                                                                        \
        const SxDiracVec<Double>  &vec = *((const SxDiracVec<Double> *)in);   \
        SxDiracVec<Complex16> &res = *((SxDiracVec<Complex16> *)o);           \
        res = PROJECTOR (b, vec);                                             \
     }                                                                        \
                                                                              \
     virtual void projectTo (const TO *b, const void *in, void *o,            \
                             const SxComplex8 &, const float &) const         \
     {                                                                        \
        const SxDiracVec<Complex8> &vec                                       \
                 = *((const SxDiracVec<Complex8> *)in);                       \
        SxDiracVec<Float> &res = *((SxDiracVec<Float> *)o);                   \
        res = PROJECTOR (b, vec);                                             \
     }                                                                        \
     virtual void projectTo (const TO *b, const void *in, void *o,            \
                             const SxComplex8 &, const double &) const        \
     {                                                                        \
        const SxDiracVec<Complex8>  &vec                                      \
                 = *((const SxDiracVec<Complex8> *)in);                       \
        SxDiracVec<Double> &res = *((SxDiracVec<Double> *)o);                 \
        res = PROJECTOR (b, vec);                                             \
     }                                                                        \
     virtual void projectTo (const TO *b, const void *in, void *o,            \
                             const SxComplex8 &, const SxComplex8 &) const    \
     {                                                                        \
        const SxDiracVec<Complex8>  &vec                                      \
                = *((const SxDiracVec<Complex8> *)in);                        \
        SxDiracVec<Complex8> &res = *((SxDiracVec<Complex8> *)o);             \
        res = PROJECTOR (b, vec);                                             \
     }                                                                        \
     virtual void projectTo (const TO *b, const void *in, void *o,            \
                             const SxComplex8 &, const SxComplex16 &) const   \
     {                                                                        \
        const SxDiracVec<Complex8>  &vec                                      \
                = *((const SxDiracVec<Complex8> *)in);                        \
        SxDiracVec<Complex16> &res = *((SxDiracVec<Complex16> *)o);           \
        res = PROJECTOR (b, vec);                                             \
     }                                                                        \
                                                                              \
                                                                              \
                                                                              \
     virtual void projectTo (const TO *b, const void *in, void *o,            \
                             const SxComplex16 &, const float &) const        \
     {                                                                        \
        const SxDiracVec<Complex16> &vec                                      \
                 = *((const SxDiracVec<Complex16> *)in);                      \
        SxDiracVec<Float> &res = *((SxDiracVec<Float> *)o);                   \
        res = PROJECTOR (b, vec);                                             \
     }                                                                        \
     virtual void projectTo (const TO *b, const void *in, void *o,            \
                             const SxComplex16 &, const double &) const       \
     {                                                                        \
        const SxDiracVec<Complex16>  &vec                                     \
                 = *((const SxDiracVec<Complex16> *)in);                      \
        SxDiracVec<Double> &res = *((SxDiracVec<Double> *)o);                 \
        res = PROJECTOR (b, vec);                                             \
     }                                                                        \
     virtual void projectTo (const TO *b, const void *in, void *o,            \
                             const SxComplex16 &, const SxComplex8 &) const   \
     {                                                                        \
        const SxDiracVec<Complex16>  &vec                                     \
                = *((const SxDiracVec<Complex16> *)in);                       \
        SxDiracVec<Complex8> &res = *((SxDiracVec<Complex8> *)o);             \
        res = PROJECTOR (b, vec);                                             \
     }                                                                        \
     virtual void projectTo (const TO *b, const void *in, void *o,            \
                             const SxComplex16 &, const SxComplex16 &) const  \
     {                                                                        \
        const SxDiracVec<Complex16>  &vec                                     \
                = *((const SxDiracVec<Complex16> *)in);                       \
        SxDiracVec<Complex16> &res = *((SxDiracVec<Complex16> *)o);           \
        res = PROJECTOR (b, vec);                                             \
     }                                                                        \
                                                                              \

#define FOR_ALL_VECTYPES(MACRO) \
   MACRO(Float) \
   MACRO(Double) \
   MACRO(Complex8) \
   MACRO(Complex16)


/** \brief Abstract basis-set class used for direct Dirac-like projections.

    \b SxBasis = S/PHI/nX Basis-sets 

    This class is a major class to establish the \ref page_dirac. It is 
    used to perform direct projections onto basis-sets. Direct means that
    the destination basis is explicitly known and Dirac's bra-ket notation
    is used.

    The ::SxBasis class comes with various \b projectTo \em subroutines and 
    \b projectTo \em functions. The \em subroutines are called from 
    ::SxBasis to perform indirect projections (like SxDiracVec<T>::toR, 
    SxDiracVec<T>::toG) whereas the \em functions are called from explicit
    Dirac projectors (see SxProjector.h).

    <b>Adding new basis-sets:</b>

    For each new basis-set the developer has to add:
      - a \b projectTo \em function
      - a \b projectTo \em subroutine 
        (which calls internally its projectTo function)

    \ingroup group_dirac
    \sa      \ref page_dirac
    \sa      SxBasis
    \author  Sixten Boeck
  */
class SX_EXPORT_DIRAC SxBasis
: public SxMemConsumer,
  public SxThis<SxBasis>
{
   private:
      /** \brief Pointers to registered bases

          It is sometimes necessary that basis transformations
          require auxiliary data. These are stored in one of the
          bases (and should be mutable). In order to relate this
          data correctly with the other basis, a pointer to the
          second basis is registered and stored here.

          To add a pointer use the function ::registerBasis.

          If several bases of the same type are registered, the
          id of a specific basis can be obtained via getBasisId.

          \sa \ref page_dirac */
      mutable SxArray<const SxBasis *> registeredBases;

   public:
      /// \brief Number of sampling points of the corresponding Dirac vector.
      virtual ssize_t getNElements () const = 0;

      /// \brief standard destructor
      virtual ~SxBasis ();

      /** \brief Register basis unless already registered
          @param basisPtr basis to be registered
          @return true if basis was new and registered,
                  false if basis was already registered
          @note This function needs rarely to be called outside the basis
                layer, since bases may autoregister in the first projection.
                If newly implemented bases require registration at some point,
                they should autoregister right in place. (That's why we made
                this stuff mutable).
        */
      bool registerBasis (const SxBasis &basis) const;


      /** \brief Register a basis of unknown type

          \note This is a hook-in for derived bases that
                provide special-case registration. The correct type
                must be obtained by dynamic_cast.
          \note The deregistration hook-in is deregister.
        */
      virtual void registerUnknown (const SxBasis &basis) const
      {
         registerBasis (basis);
      }

   protected:
      /** \brief Deregister a basis
        This performs the SxBasis deregistration process.
       */
      void deregisterBasis (const SxBasis *) const;
      
      /** \brief Deregister a basis from derived basis
        This function needs to be overloaded if the deregistration
        should remove some internal data of the derived basis.
        Use dynamic_cast<> to figure out the type of the deregistered 
        basis.

        Remember to call deregisterBasis at the end if you overload
        this.
          
       */
      virtual void deregister (const SxBasis *basis) const
      {
         deregisterBasis (basis);
      }
      
      /** Deregister all bases
          \note This must be called in the derived basis' destructor.
        */
      void deregisterAll ();
      
   public:

      /** \brief Get basis id
          If more than one basis of some type is registered, each
          basis gets an id, starting from 0. In this way, member data
          of the derived basis can be correctly be mapped to the basis
          in question.
          @return the id, or -1 if basis is not registered

          \example
          \code
int iFFT = getBasisId<SxRBasis> (targetBasis);
fft3d(iFFT).fftForward (data, res);
          \endcode
        */
      template <class BasisType>
      int getBasisId (const BasisType* basis) const;

      /// Check if a basis is registered
      bool isRegistered (const SxBasis *basis) const
      {
         for (int i = 0; i < registeredBases.getSize (); ++i)
            if (registeredBases(i) == basis) return true;
         return false;
      }

      /** \brief placeholder for any \f$
                    \mathrm{tr}_\mathbf{B}(\hat{\varrho}\hat{\mathbf{A}})
                 \f$  */
#define VIRTUAL_TR(T)                                                         \
     virtual Real8 tr (const SxDiracVec<T> &) const                           \
     {                                                                        \
        printf ("Trace for SxDiracVec<%s> in %s basis not yet implemented.\n",\
                #T, getType ().ascii ());                                     \
        SX_EXIT; return 0.;                                                   \
     }
      FOR_ALL_VECTYPES (VIRTUAL_TR)
   
   private:
      /// Crash on an unimplemented projection
      void projectionFailed (const char *basis, const char *types) const;
   public:
      /** \brief placeholder for any <??|X> projector */
      VIRTUAL_PROJECT_TO (SxBasis);
      /** \brief placeholder for any <R|X> projector */
      VIRTUAL_PROJECT_TO (SxRBasis);
      /** \brief placeholder for any <G|X> projector */
      VIRTUAL_PROJECT_TO (SxGBasis);
      /** \brief placeholder for any <r|X> projector */
      VIRTUAL_PROJECT_TO (SxRadBasis);
      /** \brief placeholder for any <radR|X> projector */
      VIRTUAL_PROJECT_TO (SxRadRBasis);
      /** \brief placeholder for any <radG|X> projector */
      VIRTUAL_PROJECT_TO (SxRadGBasis);
      /** \brief placeholder for any <mu|X> projector */
      VIRTUAL_PROJECT_TO (SxAOBasis);
      /** \brief placeholder for any <p|X> projector */
      VIRTUAL_PROJECT_TO (SxPartialWaveBasis);
      /** \brief placeholder for any <G,p|X> projector */
      VIRTUAL_PROJECT_TO (SxPAWBasis);

      /** \brief placeholder the Laplacian in any basis-set */
      virtual Real8 laplacian (const void *) const {
         SX_EXIT;  
         return 0.;
      }

      virtual int getNComp () const {
         SX_EXIT;
         return -1;
      }

      /** \brief Print debug information about the basis */
      virtual void print () const;

      /** Very simple description of basis a la "|R>" or "|G>"
        */
      virtual SxString getType () const = 0;

      /// Scalar product
      virtual double scalarProduct (const SxDiracVec<Double> &x,
                                    const SxDiracVec<Double> &y) const;

      /// Scalar product
      virtual SxComplex16 scalarProduct (const SxDiracVec<Complex16> &x,
                                         const SxDiracVec<Complex16> &y) const;

   protected:

      /** \brief Standard constructor */
      SxBasis () {/* empty */}

   public:

      
//    virtual SxDiracVec<T> laplacian   (const SxDiracVec<T> &) const;

};


template <class BasisType>
int SxBasis::getBasisId (const BasisType* basis) const
{
   int i, id = 0;
   // find basis of desired type
   for (i = 0; i < registeredBases.getSize (); ++i)  {
      if (registeredBases(i) == basis) return id;
      // increment id, if other basis is found
      if (dynamic_cast<const BasisType*>(registeredBases(i))) id++;
   }
   // basis is not registered
   SX_EXIT;
}

//template<class T>
//Laplacian (const SxDiracVec<T> &in)
//{
//   return in.laplacian ();
//}

//template<class T>
// SxDiracVec<T> SxBasis::laplacian (const SxDiracVec<T> &) const
//{
//   SX_EXIT;
//}

//template<class T>
// SxDiracVec<T> SxBasis::laplacian (const SxDiracVec<T> &) const
//{
//   SX_EXIT;
//}

#endif /* _SX_BASIS_H_ */
