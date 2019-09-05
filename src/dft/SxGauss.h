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

#ifndef _SX_GAUSS_H_
#define _SX_GAUSS_H_


#include <SxDFT.h>
#include <SxConfig.h>
#include <SxPW.h>
#include <SxString.h>
#include <SxGkBasis.h>
#include <SxRBasis.h>


/** \brief initialises a set of Gaussians from input file and allows
           projections of wavefunctions onto these Gaussians

    \b SxGauss = S/PHI/nX Gaussians

    \author Matthias Wahn, wahn@fhi-berlin.mpg.de */

class SX_EXPORT_DFT SxGauss
{
   public:
      /** Enumeration type speciefies, whether the real space representation
          of the Gaussians shall be given on a radial grid or on the (real
          space) FFT mesh.

          \note default: radial space
       */
      enum RealSpaceMesh { Radial, FFT, None };

      /** constructor */
      SxGauss ();

      /** constructor */
      SxGauss (SxPtr<SxGkBasis> gkPtr_, const SxCell &cell, const SxSymbolTable *top);

      /** destructor */
      ~SxGauss ()  { /* empty */ }

      /** set the |G+0> basis */
      void setGBasis (SxGBasis &gBasis);

      /** set the mesh on which the real space shall be represented */
      void setRealSpaceMesh (RealSpaceMesh realSpaceMesh_);

      /** number of initial guess functions

          \note It must be equal to the number of bands considered.
       */
      int nGauss;

      /** positions of the initial guess functions in Cartesian
          coordinates */
      SxList<Coord> positions;

      /** widths of the Gauss functions */
      SxList<double> widths;

      /** labels, e.g., bond site, off bond site, on Ga atom, ... */
      SxList<SxString> labels;

      /** dump parameters from input file */
      void print () const;

      /** compute Gaussians centered around the origin on the radial mesh */
      void computeGaussiansRadialMesh ();

      /** compute Gaussians on the FFT mesh, considering their periodic
          images within the neighboured unit cells */
      void computeGaussiansFFTMesh ();

      /** transformes the Gaussians created on the radial or the FFT mesh into
          their G space representation */
      void computeGaussiansGSpace ();

      /** transformes the Gaussians created on the radial or the FFT mesh into
          their |G+k> space representation */
      void computeGaussiansGkSpace ();

      /** kept just for comparison with the old code */
      void computeGaussiansOld ();

      /** compute structure factors of the Gaussians in |G> space,
          required to shift them to the specified positions */
      void computeStructureFactorsG ();

      /** compute structure factors of the Gaussians in |G+k> space,
          required to shift them to the specified positions */
      SxDiracVec<TPrecCoeffG> computeStructureFactorsGk (int ik, int iGauss);

      /** write Gaussians to output file */
      void write ();

      /** set of functions being projected onto the Gaussians */
      SxPW                               projections;

      /** projection on a set of wave functions according to Eq. (62) of
          <a href="http://prola.aps.org/abstract/PRB/v56/i20/p12847_1"> Phys.
          Rev. B, 56, 12847 (1997) </a>
       */
      void project (const SxPW &waves, int iBottom, int iTop);
      void projectOld (const SxPW &waves, int iBottom, int iTop);

      /** orthonormalise projections according to Loewdin or Gram-Schmidt
          orthogonalisation scheme */
      void orthonormalizeProjections (SxOrthoMethod method);

      /** Yields the transformation describing the projection from the
          Gaussians onto the Bloch states according to Eq. (62) in
          <a href="http://prola.aps.org/abstract/PRB/v56/i20/p12847_1">
          Phys. Rev. B, 56, p. 12 847 (1997) </a>. */
      SxArray<SxArray<SxDiracMat<TPrecCoeffG> > > getProjectionTrafo () const;

      /** Yields the L&ouml;wdin transformation orthonormalising the results
          of the Gauss-projection. */
      SxArray<SxArray<SxDiracMat<TPrecCoeffG> > > getLoewdinTrafo () const;

      /// TODO: write add-on sxrepeatfunction with template for all gen. types
      //  void repeatFunction (............)

      /// TODO: move this function to a better spot in the code
      static SxArray<Coord>
         computeMeshPoints (const RelVec &dim, const RelVec &repeat,
                            const SxCell &cell_);

   protected:

      /** pointer to the set of |G+k> basises */
      /*const*/ SxPtr<SxGkBasis>                   gkSetPtr;

      /** pointer to the |G+0> basis */
      SxGBasis                                    *gPtr;

      /** the real space basis */
      SxRBasis                                     R;

      /** unit cell */
      SxCell                                       cell;

      /** the FFT mesh */
      SxVector3<Int>                               meshDim;

      // TODO: to be removed -- debugging only
      SxArray<SxDiracVec<TPrecPhi> >               gaussians;

      /** Gaussian "bell curves" centered around the origin on a
          1-dimensional logarithmic mesh (radial mesh)
      
          \note These functions are not created as periodic functions. */
      SxArray<SxDiracVec<TPrecPhi> >               gaussiansRadialMesh;

      /** periodic Gaussian "bell curves" in real space given on the FFT mesh */
      SxArray<SxDiracVec<TPrecPhi> >               gaussiansFFTMesh;

      /** periodic Gaussian "bell curevs" transformed to the G space */
      SxArray<SxDiracVec<TPrecCoeffG> >            gaussiansGSpace;

      /** periodic Gaussians transformed to the set of |G+k> spaces
      
          \note Please note, that the periodicity region does not
                correspond to the unit cell, unless
                \f$ \mathbf{k} = \mathbf{0} \f$.
       */
      SxPW                                         gaussiansGkSpace;

      /** pointer to the radial basis sets for the Gaussians */
      SxPtr<SxRadBasis>                    rPtr;

      RealSpaceMesh                                realSpaceMesh;

      /** structure factors of the Gaussians in |G> space, required to shift
          them to the specified positions */
      SxArray<SxDiracVec<TPrecCoeffG> >            structureFactorsG;

      /** The transformation describing the projection from the Gaussians
          onto the Bloch states according to Eq. (62) in
          <a href="http://prola.aps.org/abstract/PRB/v56/i20/p12847_1">
          Phys. Rev. B, 56, p. 12 847 (1997) </a> */
      SxArray<SxArray<SxDiracMat<TPrecCoeffG> > >  trafoProjection;

      /** The L&ouml;wdin transformation orthonormalising the results of the
          projection */
      SxArray<SxArray<SxDiracMat<TPrecCoeffG> > >  trafoLoewdin;
};

#endif /* _SX_GAUSS_H_ */
