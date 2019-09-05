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

#ifndef _SX_DIRAC_HPP_
#define _SX_DIRAC_HPP_


#include <SxBasis.h>

template<class T>
SxDiracVec<T>::SxDiracVec (const SxBasis &basis)
   : elements(NULL)
{
   init (basis.getNElements ());
   setBasis (&basis);
}

template<class T>
SxDiracVec<T>::SxDiracVec (const SxBasis *basisPtr)
   : elements(NULL)
{
   SX_CHECK (basisPtr);
   init (basisPtr->getNElements ());
   setBasis (basisPtr);
}

template<class T>
SxDiracVec<T>::SxDiracVec (const SxBasis &basis, ssize_t nCol)
   : elements(NULL)
{
   ssize_t nRow = basis.getNElements ();
   init (nRow * nCol);
   reshape (nRow, nCol);
   setBasis (&basis);
}

template<class T>
void SxDiracVec<T>::setBasis (const SxBasis *ptr)
{
   handle->auxData.basisPtr = ptr;
}

template<class T>
void SxDiracVec<T>::setBasis (const SxBasis &basis)
{
   handle->auxData.basisPtr = &basis;
}

template<class T>
const SxBasis *SxDiracVec<T>::getBasisPtr () const
{
   return (const SxBasis *)handle->auxData.basisPtr;
}

template <class T>
template <class B>
const B& SxDiracVec<T>::getBasis () const
{
   const B* ptr = dynamic_cast<const B*>(getBasisPtr ());
   SX_CHECK (ptr);
   return *ptr;
}

//template<class T>
//SxDiracVec<T> SxDiracVec<T>::laplacian () const
//{
//   SX_CHECK (getBasisPtr());
//   return getBasisPtr()->laplacian (*this);
//}

template<class T>
template<class B>
SxDiracVec<typename B::TBasisType> SxDiracVec<T>::to (const B& destBasis) const
{
   // static cast makes sure at compile time that B is derived from SxBasis
   (void)static_cast<const SxBasis*>(&destBasis);
   // --- check whether Dirac vector has been set up
   SX_CHECK (getBasisPtr());
   const SxBasis  *srcBasis = getBasisPtr();

   //return srcBasis->projectTo (destBasis, (const void *)this);
   SxDiracVec<typename B::TBasisType> res;
   srcBasis->projectTo (&destBasis, (const void *)this, (void *)&res,
                        (typename T::Type)0.,
                        (typename B::TBasisType::Type)0.);
   return res;
}

template<class T>
typename T::Real SxDiracVec<T>::laplacian () const
{
   SX_CHECK (getBasisPtr());
   const SxBasis *srcBasis  = getBasisPtr ();

   return (typename T::Real)srcBasis->laplacian ((const void *)this);
}

template<class T>
SxDiracVec<T> SxDiracVec<T>::getComponent(int iComp) const
{
   SX_CHECK (getBasisPtr());
   int nComp = getBasisPtr()->getNComp();
   int size = (int)(*this).getSize() / nComp;
   SxDiracVec out;
   out.resize(size);

   out <<= (*this)(SxIdx(iComp * size, (iComp + 1) * size - 1));

   return out;
}

template<class T>
void SxDiracVec<T>::setComponent(int iComp, const SxDiracVec<T> &component)
{
   SX_CHECK (getBasisPtr());
   int nComp = getBasisPtr()->getNComp();
   int size = (int)(*this).getSize() / nComp;
//   int compSize = component.getSize();

   (*this)(SxIdx(iComp * size, (iComp + 1) * size - 1)) = component;
}

template <class T>
SxDiracVec<T> & SxDiracVec<T>::setAux (const SxDiracAux &aux)
{
    SX_CHECK (handle);
    /// Check that basis size matches vector dimension
    SX_CHECK (!aux.basisPtr ||
              (   getSize () == aux.basisPtr->getNElements ()
              ||   nRows () == aux.basisPtr->getNElements ()),
              aux.basisPtr->getNElements (), getSize (), nRows () );
    handle->auxData = aux;
    return *this;
}

#endif /* _SX_DIRAC_HPP_ */
