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

#ifndef _SX_DIRAC_HH_
#define _SX_DIRAC_HH_


// template<class T>
// class SxDiracVec
// {
//     ...   see (SxVec.h)

       explicit SxDiracVec (const SxBasis&);
       explicit SxDiracVec (const SxBasis &basis, ssize_t nCol);
       explicit SxDiracVec (const SxBasis*);

       void setBasis (const SxBasis*);
       void setBasis (const SxBasis&);
       const SxBasis *getBasisPtr () const;

       template<class B>
       const B& getBasis () const;

       template <class B>
       SxDiracVec<typename B::TBasisType> to (const B&) const;
       /** \brief Get one component of a multicomponent wave (as copy) */
       SxDiracVec<T> getComponent (int) const;
       void setComponent (int, const SxDiracVec<T> &);

       typename T::Real          laplacian () const;

       /// Set aux data
       SxDiracVec<T> &setAux (const SxDiracAux &aux);


//     ...   see (SxVec.h)
// };

#endif /* _SX_DIRAC_HH_ */
