
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

#ifndef _SX_SYM_OP_H_
#define _SX_SYM_OP_H_
#include <SxMatrix3.h>
#include <SxVector3.h>
#include <SxGeom.h>
#include <SxPrecision.h>

/** \brief Symmetry transformation (affine)

    \b SxSymOp = S/PHI/nX affine symmetry operation

    This class provides an affine transformation, used mainly
    for symmetries:

    \f[
      x' = S ^ x + t
    \f]

    consisting of a rotational part S and a translational part t.
    Note that first the rotation and then the translation is applied.

    \author Christoph Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_GEOM SxSymOp
{
   public:
      /// Rotation
      SxMatrix3<TPrecTauR> rot;
      /// Translation
      SxVector3<TPrecTauR> trans;
      
      /// Constructor
      SxSymOp (const SxMatrix3<TPrecTauR> &rotIn 
               = SxMatrix3<TPrecTauR>(1,0,0,0,1,0,0,0,1),
               const SxVector3<TPrecTauR> transIn
               = SxVector3<TPrecTauR> (0., 0., 0.))
         : rot(rotIn), trans(transIn)
      {
         SX_CHECK(fabs(fabs(rot.determinant ()) - 1.) < 1e-8,
                  rot.determinant ());
      }

      /// Apply to single vector
      inline SxVector3<TPrecTauR> 
      operator^(const SxVector3<TPrecTauR> &x) const
      {
         SxVector3<TPrecTauR> res(rot ^ x);
         res += trans;
         return res;
      }

      /// Self-apply
      inline SxSymOp operator^(const SxSymOp &in) const
      {
         return SxSymOp(rot ^ in.rot, operator^(in.trans)); 
      }

      /// Inverse operator
      inline SxSymOp inverse () const
      {
         SX_CHECK(fabs(fabs(rot.determinant ()) - 1.) < 1e-8,
                  rot.determinant ());
         // inv := (S^-1, -S^-1 t)
         SxSymOp res(rot.transpose ()); 
         res.trans = -(res.rot ^ trans);
         return res;
      }

      bool isSymmorphic (double eps = 1e-8) const
      {
         return    fabs(trans(0)) < eps 
                && fabs(trans(1)) < eps
                && fabs(trans(2)) < eps;
      }

      inline bool operator== (const SxSymOp &op) const {
         SxSymOp res = *this ^ op.inverse ();
         if (res.rot.trace () < 3 || res.trans.norm () > 1e-6)  
         { return false; }   //WARNING trans periodic!
         return true;
      }
};

inline std::ostream &operator<< (std::ostream &out, const SxSymOp &op)
{
   out << "S=" << op.rot;
   if (!op.isSymmorphic ()) out << ",t=" << op.trans;
   return out;
} 

#endif /* _SX_SYM_OP_H_ */
