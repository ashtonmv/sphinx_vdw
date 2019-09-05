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

#ifndef _SX_LAPLACIAN_H_
#define _SX_LAPLACIAN_H_

/** \brief ...

    \b SxClass = S/PHI/nX ...

    ....

    \author Sixten Boeck, boeck@mpie.de */
class SxLaplacian
{
   public:
      inline SxLaplacian () { }
      inline ~SxLaplacian () { }

};


/** \brief ...

    \b SxClass = S/PHI/nX ...

    ....

    \author Sixten Boeck, boeck@mpie.de */
template<class T>
class SxLaplacianPsi
{
   public:
      inline SxLaplacianPsi (const SxDiracVec<T> &v_) : v(v_) { }
      inline       SxDiracVec<T> getVec ()       { return v; }
      inline const SxDiracVec<T> getVec () const { return v; }

   protected:
      SxDiracVec<T> v;
};

#endif /* _SX_LAPLACIAN_H_ */
