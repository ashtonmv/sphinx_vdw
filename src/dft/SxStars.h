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

#ifndef _SX_STARS_H_
#define _SX_STARS_H_

#include <SxDFT.h>
#include <SxConfig.h>
#include <SxGkBasis.h>


/** \brief Calculates the star(s) of a k-point set.

    \b SxStars = S/PHI/nX Stars of k-points

    \author Matthias Wahnovic, wahn@fhi-berlin.mpg.de */
class SX_EXPORT_DFT SxStars
{
   public:

      SxStars ();
      SxStars (const SxCell                      &cell_,
               const SxArray<SxVector3<TPrecG> > &kVec_,
               const SxVector<TPrecWeights>      &weightsOrig_,
                     bool                         useInvSymmetry_=true);
      ~SxStars ();

      void compute ();

      void computeStar (int ik);

      /// Get a star (in relative coordinates)
      SxArray<SxVector3<TPrecG> > getStar(int ik) const;
      PrecWeights getWeight (int ik) const;
      SxArray<SxArray<int> > getRotations () const;
      int getRotation (int ik, int iRot) const;

      int  getSizeOfMaxStar () const;

   protected:

      SxCell                                  cell;
      SxArray<SxVector3<TPrecG> >             kVec;
      bool                                    useInvSymmetry;

   public:

      int                                     nkOrig;

      /// The stars (in relative coordinates)
      SxArray<SxArray<SxVector3<TPrecG> > >   stars;
      SxArray<PrecWeights>                    weights;
      SxArray<SxArray<int> >                  rotations;

   protected:
      SxVector<TPrecWeights>                weightsOrig;

};

#endif /* _SX_STARS_H_ */
