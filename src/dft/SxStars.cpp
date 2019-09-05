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

#include <SxStars.h>
#include <SxAtomicStructure.h>
#include <SxRotation.h>


SxStars::SxStars ()
   : useInvSymmetry (true)
{
   // empty
}

SxStars::SxStars (const SxCell                      &cell_,
                  const SxArray<SxVector3<TPrecG> > &kVec_,
                  const SxVector<TPrecWeights>      &weightsOrig_,
                        bool                         useInvSymmetry_)
   : cell        (cell_),
     kVec        (kVec_),
     useInvSymmetry (useInvSymmetry_)
{
   nkOrig = int(kVec.getSize());
   stars.resize (nkOrig);
   weights.resize (nkOrig);
   rotations.resize (nkOrig);
   weightsOrig = SxVector<TPrecWeights> (weightsOrig_);

   SX_CHECK (weightsOrig.getSize() == nkOrig,
             weightsOrig.getSize(), nkOrig);
}

SxStars::~SxStars ()
{
   // empty
}

void SxStars::compute ()
{
   SX_CHECK (nkOrig > 0, nkOrig);

   weights.resize (nkOrig);

   for (int ik = 0; ik < nkOrig; ik++)  {
      
      // --- compute stars
      computeStar (ik);

      // --- compute weigths
      weights(ik) = weightsOrig(ik) / double(stars(ik).getSize());
   }

   cout << SX_SEPARATOR;
   cout << "weights: " << weights << endl;

   // --- check sum of k-points
   PrecWeights sum = 0.;
   for (int ik = 0; ik < nkOrig; ik++)
      sum += weights(ik) * double(stars(ik).getSize());
   if ( fabs(sum - 1.) > 1.e-8 )  {
      sxprintf ("Error: Sum of k-point weights is incorrect:\n");
      sxprintf ("       sum = %15.12f (should be 1)\n", sum);
      SX_QUIT;
   }
}

void SxStars::computeStar (int ik)
{
   SX_CHECK (ik >= 0 && ik < nkOrig, ik, nkOrig);

   SxVector3<TPrecG> k       = kVec(ik);  // k in cart. coordinates
   SX_CHECK (cell.symGroupPtr);
   SxArray<SymMat>   symOp   = cell.symGroupPtr->getSymmorphic ();
   SxCell            recCell = cell.getReciprocalCell ();

   cout << SX_SEPARATOR;
   cout << "| star of k = " << recCell.carToRel(k)
        << "   (relative coord.)" << endl;
   cout << SX_SEPARATOR;

   SxMatrix3<TPrecG>            S;
   SxVector3<TPrecG>            kRot;
   int                          iOp, nOp = int(symOp.getSize()), ikRot;
   double                       diff, kNorm, kRotNorm;
   bool                         append;
   SxList<SxVector3<TPrecG > >  star;
   SxList<int>                  rots;
   SxMatrix3<TPrecG>            Id = SxMatrix3<TPrecG> (1., 0., 0.,
                                                        0., 1., 0.,
                                                        0., 0., 1.);
   
   // --- append k itself to the list
   star.resize (0);
   star.append (k);

   rots.resize (0);
   rots.append (-1);  // The -1 is just a placeholder here. It will be replaced
                      // by the index of the identity-symmetry-operator later.

   SxArray<SxSymType> symType (nOp);
   SxArray<int> starFromOp (2 * nOp);

   // --- compute star
   for (iOp = 0; iOp < nOp; iOp++)  {

      S    = symOp(iOp);
      symType (iOp) = SxRotation(S).getType ();
      kRot = S ^ k;                      // k' = S k in cart. coordinates

      append = true;

      for (ikRot = 0; ikRot < star.getSize(); ikRot++)  {

         diff = recCell.getMapped (kRot - star(ikRot), SxCell::Origin)
                .normSqr();

         if (diff < 1.e-6)  {
            append = false;
            break;
         }
      }

      bool selfInverse = symType(iOp).isSelfInverse ();
      if (append)  {
         star.append (kRot);
         rots.append (selfInverse ? iOp : -1);
      } else {
         // self-inverse operation ? => use as generator
         if (selfInverse && rots(ikRot) == -1)
         {
            rots(ikRot) = iOp;
         }
      }

      starFromOp (iOp) = ikRot;

      // --- find out index of identity
      if (S == Id)  rots(0) = iOp;
   }

   // --- points coming out of time reversal symmetry
   if (useInvSymmetry)  {

      for (iOp = 0; iOp < nOp; iOp++)  {

         S    = (-1.) * symOp(iOp);  // -I * S
         kRot = S ^ k;

         append = true;

         for (ikRot = 0; ikRot < star.getSize(); ikRot++)  {

            diff = recCell.getMapped (kRot - star(ikRot), SxCell::Origin)
                   .normSqr();

            if (diff < 1.e-6)  {
               append = false;
               break;
            }
         }

         bool selfInverse = symType(iOp).isSelfInverse ();
         if (append)  {
            star.append (kRot);
            // + nOp marks time reversal symmetry
            rots.append (selfInverse ? (iOp + nOp) : -1);
         } else {
            if (selfInverse && rots(ikRot) == -1)
            {
               // self-inverse operation ? => use as generator
               rots(ikRot) = iOp + nOp;
            }
         }
         starFromOp (iOp+nOp) = ikRot;
      }
   }

   // --- make sure for each star generating operation that the inverse
   //     is included, too (via itself or as a pair of star members)
   // --- self-paired members
   for (ikRot = 0; ikRot < rots.getSize (); ++ikRot)  {
      if (rots(ikRot) >= 0)  {
         cout << "member " << (ikRot+1) << " from self-adjoint symmetry"
              << endl;
      }
   }
   // --- search other pairs
   for (iOp = 0; iOp < 2 * nOp; ++iOp)  {
      ikRot = starFromOp(iOp);
      // --- check whether this star member is symmetry-paired
      if (rots(ikRot) >= 0) continue;

      // --- find the inverse
      S = symOp(iOp % nOp);
      int jOp;
      for (jOp = 0; jOp < nOp; ++jOp) 
         if ((S ^ symOp(jOp)) == Id) break;
      SX_CHECK (jOp < nOp);
      if (iOp >= nOp) jOp += nOp;

      int jkRot = starFromOp(jOp);

      // --- check whether star member from S^{-1} is symmetry-paired
      if (rots(jkRot) >= 0) continue;

      // => this is a new symmetry pair
      rots(ikRot) = iOp;
      rots(jkRot) = jOp;
      cout << "symmetry pair: members " << (ikRot+1) << "<->" << (jkRot+1)
           << " via symmetries " << (iOp+1) << "/" << (jOp+1) << endl;
   }

   // --- print out generating symmetry ops
   bool ok = true;
   for (ikRot = 0; ikRot < rots.getSize (); ++ikRot)  {
      iOp = rots(ikRot);
      if (iOp == -1)  {
         cout << "Star symmetry inconsistency: Missing symmetry pair." << endl;
         ok = false;
      }
      cout << "Star member " << (ikRot + 1) << " is generated from "
           << symType (iOp % nOp).identifier;
      if (iOp >= nOp) cout << " and time reversal";
      cout << endl;
   }
   if (!ok)  { SX_EXIT; }

   int nkStar = int(star.getSize());
   SX_CHECK (rots.getSize() == nkStar, rots.getSize(), nkStar);

   // sxprintf ("nkStar = %d\n", nkStar);

   stars(ik).resize (nkStar);
   rotations(ik).resize (nkStar);

   // --- convert to relative coordinates
   for (ikRot = 0; ikRot < nkStar; ikRot++)  {
      stars(ik)(ikRot) = recCell.carToRel( star(ikRot) );
      rotations(ik)(ikRot) = rots(ikRot);
      sxprintf ("k(%d) = ", ikRot);
      cout << stars(ik)(ikRot) << endl;
   }

   // --- check equal norm of all points
   kNorm = k.norm ();
   for (ikRot = 0; ikRot < nkStar; ikRot++)  {
      kRotNorm = star(ikRot).norm ();

      if ( fabs(kNorm - kRotNorm) > 1.e-6 )  {
         sxprintf ("Error: Each member of the star should have the same "
                   "norm.\n");
         SX_EXIT;
      }
   }
}

SxArray<SxVector3<TPrecG> > SxStars::getStar (int ik) const
{
   SX_CHECK (ik >= 0, ik);
   SX_CHECK (ik < stars.getSize(), ik, stars.getSize());
   SX_CHECK (stars(ik).getSize() > 0, stars(ik).getSize());
   return stars(ik);
}

PrecWeights SxStars::getWeight (int ik) const
{
   SX_CHECK (ik >= 0, ik);
   SX_CHECK (ik < weights.getSize(), ik, weights.getSize());
   return weights(ik);
}

SxArray<SxArray<int> > SxStars::getRotations () const
{
   SX_CHECK (rotations(0).getSize() > 0, rotations(0).getSize());
   return rotations;
}

int SxStars::getRotation (int ik, int iRot) const
{
   SX_CHECK (ik >= 0 && iRot >= 0, ik, iRot);
   SX_CHECK (ik < rotations.getSize(), ik, rotations.getSize());
   SX_CHECK (iRot < rotations(ik).getSize(), iRot, rotations(ik).getSize());
   return rotations(ik)(iRot);
}

int SxStars::getSizeOfMaxStar () const
{
   if (stars(0).getSize() == 0)  return -1;

   int max = -1;
   for (int ik = 0; ik < nkOrig; ik++)
      if (stars(ik).getSize() > max)  max = int(stars(ik).getSize());

   return max;
}
