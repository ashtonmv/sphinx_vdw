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

#ifndef _SX_PROJECTOR_H_
#define _SX_PROJECTOR_H_

#include <SxBasis.h>
#include <SxDirac.h>
#include <SxLaplacian.h>

#include <SxRBasis.h>
#include <SxGBasis.h>
#include <SxRadBasis.h>
#include <SxRadRBasis.h>
#include <SxRadGBasis.h>
#include <SxAOBasis.h>

#define SUM(a,b)  (b)

/** \brief Template to repesent an arbitrary projector \f$
              \langle A | B \rangle
           \f$
  \sa     \ref page_dirac
  \ingroup group_dirac
  \author Sixten Boeck
  */
template<class Bra, class Ket>
class SxProjector
{
   public:
      inline SxProjector (const Bra *a, const Ket *b)  {
         SX_CHECK (a);
         SX_CHECK (b);
         braPtr = a;
         ketPtr = b;
      }

      inline const Bra *bra() const  {
         return braPtr;
      }

      inline const Ket *ket() const  {
         return ketPtr;
      }

   protected:
      const Bra *braPtr;
      const Ket *ketPtr;
};


// --------------------------------------------------------------------------


/** < R | G >
  \sa     \ref page_dirac
  \author Sixten Boeck
  */
inline SxProjector<SxRBasis,SxGBasis> 
operator| (const SxRBasis &a, const SxGBasis &b)
{
   return SxProjector<SxRBasis,SxGBasis> (&a, &b);
}

/** < G | R >
  \sa     \ref page_dirac
  \author Sixten Boeck
  */
inline SxProjector<SxGBasis,SxRBasis> 
operator| (const SxGBasis &a, const SxRBasis &b)
{
   return SxProjector<SxGBasis,SxRBasis> (&a, &b);
}

/** < r | G >
  \sa     \ref page_dirac
  \author Sixten Boeck
  */
inline SxProjector<SxRadBasis,SxGBasis> 
operator| (const SxRadBasis &a, const SxGBasis &b)
{
   return SxProjector<SxRadBasis,SxGBasis> (&a, &b);
}

/** < G | r >
  \sa     \ref page_dirac
  \author Sixten Boeck
  */
inline SxProjector<SxGBasis,SxRadBasis> 
operator| (const SxGBasis &a, const SxRadBasis &b)
{
   return SxProjector<SxGBasis,SxRadBasis> (&a, &b);
}

/** < mu | G >
  \sa     \ref page_dirac
  \author C. Freysoldt (copied from Sixten ;-)
  */
inline SxProjector<SxAOBasis,SxGBasis> 
operator| (const SxAOBasis &a, const SxGBasis &b)
{
   return SxProjector<SxAOBasis,SxGBasis> (&a, &b);
}

/** < G | mu >
  \sa     \ref page_dirac
  \author C. Freysoldt (copied from Sixten ;-)
  */
inline SxProjector<SxGBasis,SxAOBasis> 
operator| (const SxGBasis &a, const SxAOBasis &b)
{
   return SxProjector<SxGBasis,SxAOBasis> (&a, &b);
}

/** < p | G >
  \sa     \ref page_dirac
  \author C. Freysoldt (copied from Sixten ;-)
  */
inline SxProjector<SxPartialWaveBasis,SxGBasis> 
operator| (const SxPartialWaveBasis &a, const SxGBasis &b)
{
   return SxProjector<SxPartialWaveBasis,SxGBasis> (&a, &b);
}

/** < G | p >
  \sa     \ref page_dirac
  \author C. Freysoldt (copied from Sixten ;-)
  */
inline SxProjector<SxGBasis,SxPartialWaveBasis> 
operator| (const SxGBasis &a, const SxPartialWaveBasis &b)
{
   return SxProjector<SxGBasis,SxPartialWaveBasis> (&a, &b);
}

// --------------------------------------------------------------------------

/** < B | vec >
  \sa     \ref page_dirac
  \author Sixten Boeck
  */
template<class T, class B>
SxDiracVec<typename B::TBasisType>
operator| (const B &b, const SxDiracVec<T> &v)
{
   SX_CHECK (static_cast<const SxBasis*> (&b)); // must be a basis
   SX_CHECK (v.getBasisPtr());
   SxDiracVec<typename B::TBasisType> res;
   v.getBasisPtr()->projectTo (&b, (const void *)&v, (void *)&res, 
                               (typename T::Type)0.,
                               (typename B::TBasisType::Type)0.);
   return res;
}

/** < B | vec >
  \sa     \ref page_dirac
  \author Sixten Boeck
  */
template<class T>
SxDiracVec<T>
operator| (const SxBasis &b, const SxDiracVec<T> &v)
{
   SX_CHECK (v.getBasisPtr());
   SxDiracVec<T> res;
   v.getBasisPtr()->projectTo (&b, (const void *)&v, (void *)&res, 
                               (typename T::Type)0.,
                               (typename T::Type)0.);
   return res;
}


// --------------------------------------------------------------------------

/** < vec | R >
  \sa     \ref page_dirac
  \author Sixten Boeck
  */
template<class T>
SxDiracVec<TRBasisType> operator| (const SxDiracVec<T> &v, const SxRBasis &b)
{
   return ( b | v.conj() ).conj();
}

/** < vec | r >
  \sa     \ref page_dirac
  \author Sixten Boeck
  */
template<class T>
SxDiracVec<TRadBasisType> operator| (const SxDiracVec<T> &v, 
                                     const SxRadBasis &b)
{
   return ( b | v ).conj();
}


/** < vec | G >
  \sa     \ref page_dirac
  \author Sixten Boeck
  */
template<class T>
SxDiracVec<TGBasisType> operator| (const SxDiracVec<T> &v, const SxGBasis &b)
{
   return ( b | v.conj() ).conj();
}

/** < vec | mu >
  \sa     \ref page_dirac
  \author C. Freysoldt, copied from Sixten
  */
template<class T>
SxDiracVec<TAOBasisType> operator| (const SxDiracVec<T> &v, const SxAOBasis &b)
{
   return ( b | v ).conj ();
}
/** < vec | g >
  \sa     \ref page_dirac
  \author B. Lange, copied from Sixten
  */
template<class T>
SxDiracVec<TRadGBasisType> operator| (const SxDiracVec<T> &v, 
                                      const SxRadGBasis &b)
{
   return ( b | v ).conj();
}

// --------------------------------------------------------------------------
template<class T>
SxLaplacianPsi<T> operator| (const SxDiracVec<T> &v, const SxLaplacian &)
{
   return SxLaplacianPsi<T> (v);
}

template<class T>
typename T::Real operator| (const SxLaplacianPsi<T> &lv, const SxDiracVec<T> &v)
{
   SX_CHECK (v.getSize() == lv.getVec().getSize(),
             v.getSize(),   lv.getVec().getSize());
   SX_CHECK (v.getBasisPtr() == lv.getVec().getBasisPtr());
   return v.laplacian ();
}

// --------------------------------------------------------------------------


/** Projector * Dirac
  \sa     \ref page_dirac
  \author Sixten Boeck
  */
template<class Bra, class Ket>
SxDiracVec<typename Bra::TBasisType> 
operator* (const SxProjector<Bra,Ket> &proj, 
           const SxDiracVec<typename Ket::TBasisType> &vec)
{
   // --- check that basis sets are matching, e.g.
   //     | A > < B |  is not allowed
   SX_CHECK (proj.ket());
   SX_CHECK (proj.ket() == vec.getBasisPtr());
   SxDiracVec<typename Bra::TBasisType> res;
   proj.ket()->projectTo (proj.bra(), (const void *)&vec, (void *)&res, 
                         (typename Ket::TBasisType::Type)0.,
                         (typename Bra::TBasisType::Type)0.);
   return res;
}

/** Dirac * Projector
  \sa     \ref page_dirac
  \author Sixten Boeck
  */
template<class Bra, class Ket>
SxDiracVec<typename Bra::TBasisType> 
operator* (const SxDiracVec<typename Ket::TBasisType> &vec,
           const SxProjector<Bra,Ket> &proj)
{
   // --- check that basis sets are matching, e.g.
   //     | A > < B |  is not allowed
   SX_CHECK (proj.bra());
   SX_CHECK (proj.bra() == vec.getBasisPtr());
   SxDiracVec<typename Bra::TBasisType> res;
   proj.bra()->projectTo (proj.ket(), (const void *)&vec, (void *)&res, 
                         (typename Ket::TBasisType::Type)0.,
                         (typename Bra::TBasisType::Type)0.);
   return res;
}


// --------------------------------------------------------------------------
/** \brief placeholder for \f$ \mathrm{tr} (\hat{\varrho} \hat{A}) \f$

    Dependent on the basis-set this function returns the expectation
    value of an operator \f$ \langle \hat{A} \rangle \f$.
    Using the notation of a density matrix the expectation value can
    be evaluated according to
    \f[
       \langle \hat{A} \rangle
         =
       \mathrm{tr} (\hat{\varrho} \hat{A})
    \f]
    with the density matrix
    \f[
       \hat{\varrho} = \sum_i | \Psi_i \rangle  \langle \Psi_i |
    \f]

    The actual implementation does not use a density matrix. Instead the
    product of the density with an operator is taken as argument which
    is a Dirac vector again. This Dirac vector is aware of the basis-set
    and calls the corresponding ::tr function of this basis-set. 
    Each basis-set should overload the corresponding ::tr function.
    Since the ::tr function is basis-set dependent we write the
    basis-set as index 
    \f[
       \langle \hat{A_{\mathbf{B}}} \rangle
         =
       \mathrm{tr}_{\mathbf{B}} (\hat{\varrho} \hat{A})
    \f]
    with \b B beeing the basis-set.
    \param  in  the evaluated product of \f$ \hat{\varrho} \hat{A} \f$ 
    \return     \f$ \mathrm{tr} (\hat{\varrho} \hat{A} ) \f$*/
template<class T> 
Real8 tr (const SxDiracVec<T> &in)
{
   SX_CHECK (in.getBasisPtr());
   return in.getBasisPtr()->tr (in);
}

template <class T>
typename T::Type
operator| (const SxDiracVec<T> &x, const SxDiracVec<T> &y)
{
   return x.getBasisPtr ()->scalarProduct (x, y);
}

#endif /* _SX_PROJECTOR_H_ */
