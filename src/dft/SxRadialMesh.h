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

#ifndef _SX_RADIAL_MESH_H_
#define _SX_RADIAL_MESH_H_
#include <SxDFT.h>
#include <SxRadBasis.h>
#include <SxDirac.h>
#include <SxYlm.h>

/** \brief Radial mesh (multiple l,m) container

    \b SxClass = S/PHI/nX Radial Mesh

    Provides a container for radial meshes. This is needed for
    densities or potentials in radial coordinates.

    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DFT SxRadialMesh
{
   public:
      /// Max l component
      int lmax;

      /// The mesh data
      SxDiracMat<Double> meshData;

      /// Constructor
      SxRadialMesh () : lmax(-1) { /* empty */ }

      /// Sized constructor
      SxRadialMesh (int nr, int lmaxIn)  
      : lmax(lmaxIn), meshData(nr, (lmaxIn+1)*(lmaxIn+1))
      {
         SX_CHECK (nr > 0, nr);
         SX_CHECK(lmax >= 0, lmax);
      }

      /// Sized and initialized constructor
      SxRadialMesh (ssize_t nr, int lmaxIn, double val)  
      : lmax(lmaxIn), meshData(nr, (lmaxIn+1)*(lmaxIn+1))
      {
         SX_CHECK (nr > 0, nr);
         SX_CHECK(lmax >= 0, lmax);
         meshData.set (val);
      }

      /// Return reference to (l,m) channel
      SxDiracVec<Double> operator() (int l, int m)  {
         SX_CHECK (l <= lmax && l >= 0, l, lmax);
         SX_CHECK (abs(m) <= l, l, m);
         SX_CHECK(meshData.handle);
         SxDiracVec<Double> res = meshData.colRef(SxYlm::combineLm(l,m));
         res.handle->auxData.is = meshData.handle->auxData.is;
         res.handle->auxData.ia = meshData.handle->auxData.ia;
         res.handle->auxData.l = l;
         res.handle->auxData.m = m;
         return res;
      }

      /// Return reference to ir-th element of (l,m)-channel
      double &operator() (int ir, int l, int m)  {
         SX_CHECK (l >= 0 && l <= lmax, l, lmax);
         SX_CHECK (abs(m) <= l, l, m);
         SX_CHECK (ir >= 0 && ir < getNRad (), ir, getNRad ());
         return meshData(ir, SxYlm::combineLm (l,m));
      }

      /// Return reference to (l,m) channel
      const SxDiracVec<Double> operator() (int l, int m) const {
         SX_CHECK (l <= lmax && l >= 0, l, lmax);
         SX_CHECK (abs(m) <= l, l, m);
         SX_CHECK(meshData.handle);
         SxDiracVec<Double> res = meshData.colRef(SxYlm::combineLm(l,m));
         res.handle->auxData.is = meshData.handle->auxData.is;
         res.handle->auxData.ia = meshData.handle->auxData.ia;
         res.handle->auxData.l = l;
         res.handle->auxData.m = m;
         return res;
      }

      /// Return ir-th element of (l,m)-channel
      double operator() (int ir, int l, int m) const {
         SX_CHECK (l >= 0 && l <= lmax, l, lmax);
         SX_CHECK (abs(m) <= l, l, m);
         SX_CHECK (ir >= 0 && ir < getNRad (), ir, getNRad ());
         return meshData(ir, SxYlm::combineLm (l,m));
      }

      /// Resize
      void resize (int nr, int lMaxIn)  {
         SX_CHECK (nr > 0 || ((nr == 0) && (lMaxIn==0)), nr);
         SX_CHECK(lMaxIn >= 0, lMaxIn);
         meshData = SxDiracMat<Double> (nr, (lMaxIn+1)*(lMaxIn+1));
         lmax = lMaxIn;
      }

      /// Resize
      void resize (const SxRadBasis &rad, int is, int lMaxIn) {
         lmax = lMaxIn;
         int nr = (int)rad.radFunc(is).getSize ();
         resize (nr, lMaxIn);
         meshData.setBasis (&rad);
         SX_CHECK(meshData.handle);
         meshData.handle->auxData.is = is;
      }

      /// Set with constant value
      void set (double val)  {
         meshData.set (val);
      }

      /// Set atom id (keep species)
      void setAtomId (int ia)  {
         SX_CHECK (meshData.getSize ());
         meshData.handle->auxData.ia = ia;
      }

      /// Set species & atom id 
      void setAtomId (int is, int ia)  {
         SX_CHECK (meshData.getSize ());
         meshData.handle->auxData.is = is;
         meshData.handle->auxData.ia = ia;
      }

      /// Set basis pointer
      void setBasis (const SxRadBasis &rad)  {
         SX_CHECK(meshData.handle);
         meshData.setBasis (&rad);
      }

      const SxRadBasis &getBasis () const {
         const SxRadBasis *rad 
            = dynamic_cast<const SxRadBasis *>(meshData.getBasisPtr ());
         SX_CHECK (rad);
         return *rad;
      }

      /// Resize according to model mesh
      void resize (const SxRadialMesh &in)  {
         if (in.meshData.getSize () == 0)
            resize (0,0);
         else  {
            meshData = SxDiracMat<Double> (in.meshData.nRows (),
                                           in.meshData.nCols ());
            lmax = in.lmax;
            SX_CHECK (meshData.nCols () == sqr(lmax+1),
                      meshData.nCols (), lmax + 1);
            meshData.handle->auxData.is = in.meshData.handle->auxData.is;
            meshData.handle->auxData.ia = in.meshData.handle->auxData.ia;
         }
      }

      /// Get number of radial points
      int getNRad () const {
         SX_CHECK(meshData.getSize () > 0);
         return (int)meshData.nRows ();
      }

      void setIs (int is)
      {
         SX_CHECK(meshData.handle);
         SX_CHECK(!meshData.getBasisPtr () || (is >= 0
                  && is < meshData.getBasis<SxRadBasis> ().getNSpecies ()),
                  is, meshData.getBasis<SxRadBasis> ().getNSpecies ());
         meshData.handle->auxData.is = is;
      }

      int getIs () const
      {
         SX_CHECK (meshData.handle);
         return meshData.handle->auxData.is;
      }

      void operator+= (const SxRadialMesh &in)
      {
         SX_CHECK (lmax == in.lmax, lmax, in.lmax);
         SX_CHECK (meshData.nRows () == in.meshData.nRows (),
                   meshData.nRows (), in.meshData.nRows ());
         meshData += in.meshData;
      }

};

inline SxRadialMesh operator- (const SxRadialMesh &x, const SxRadialMesh &y)
{
   SxRadialMesh res;
   res.meshData = x.meshData - y.meshData;
   res.lmax = x.lmax;
   return res;
}

inline SxRadialMesh operator+ (const SxRadialMesh &x, const SxRadialMesh &y)
{
   SxRadialMesh res;
   res.meshData = x.meshData + y.meshData;
   res.lmax = x.lmax;
   return res;
}

#endif /* _SX_RADIAL_MESH_H_ */
