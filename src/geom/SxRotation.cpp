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

#include <SxRotation.h>
#include <SxMatrix.h>

SxRotation::SxRotation(const SxVector3<Double> &axis, double angle)
{
   SX_CHECK (axis.normSqr () > 1e-10);
   // --- get two orthogonal directions
   SxList<Coord> axes;
   axes << axis << Coord(1.,0.,0.) << Coord(0.,1.,0.) << Coord(0.,0.,1.);
   int i,j;
   double norm;
   for (i = 0; i < 3; /* nothing */)  {
      // orthogonalize
      for (j = 0; j < i; ++j)
         axes(i) -= (axes(i) ^ axes(j)) * axes(j);
      // normalize
      norm = axes(i).norm ();
      if (norm < 1e-6)
         axes.remove(i);
      else
         axes(i++) /= norm;
   }
   // Setup rotation around x-axis
   double cosA = cos(angle), sinA = sin(angle);
   SymMat rot (1., 0.  , 0. ,
               0., cosA,sinA,
               0.,-sinA,cosA);

   // rotate to actual axis
   SymMat U(axes(0), axes(1), axes(2));
   if (U.determinant () < 0.) U = SymMat(axes(0),axes(2),axes(1));
   (*this) = U.transpose () ^ rot ^ U;
}

bool SxRotation::operator== (const SymMat &in) {
   int i,j;
   for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++)
         if (fabs(in.m[i][j] - m[i][j]) >= 1e-5) return false;
   return true;
}


SxString SxRotation::getName (const SymMat& in)
{
   return getType(in).identifier;
}


SxSymType SxRotation::getType (const SymMat& in)
{
   const PrecTauR eps = 1e-3;
   SxSymType type;

   SX_CHECK (fabs(fabs(in.determinant ()) - 1.) < eps, in.determinant ());

   SxMatrix<TPrecTauR>::Eigensystem eig;
   eig = SxMatrix<TPrecTauR>(in).eigensystem ();

   SxList<int> one, minus;
   for (int i = 0; i < 3; i++)  {
      if (fabs(eig.vals(i).re - 1.) < eps) one << i;
      if (fabs(eig.vals(i).re + 1.) < eps) minus << i;
   }
   int nOne = int(one.getSize ());
   int nMinusOne = int(minus.getSize ());
   
   if (nOne == 3)      
      return SxSymType (SxSymType::Identity); // unit transformation 
   if (nMinusOne == 3) 
      return SxSymType(SxSymType::Inversion);

   int id; // nr of eigenvector to be printed

   if (nOne == 2 && nMinusOne == 1)  {
      // --- mirror plane
      type.classification = SxSymType::Mirror;
      id = minus(0);
   } else {
      // --- rotation axis or rotation mirror axis
      SX_CHECK (nOne == 1 || nMinusOne == 1, nOne, nMinusOne);

      if (nOne == 1)  {
         id = one(0); // rotation axis
         type.classification = SxSymType::Rotation;
      } else {
         id = minus(0); // rotation mirror axis
         type.classification = SxSymType::MirrorRotation;
      }

      // --- determine axis count
      SxComplex16 en = eig.vals((id+1)%3);
      SX_CHECK (fabs(en.absSqr () -  1.) < eps, en.re, en.im);
      SX_CHECK (fabs(en.re - 1.) > eps);
      // en = e^{2pi i/n} = (cos 2pi/n, sin 2pi/n) = (re, im)
      // => n = 2pi / atan(im/re)
      double count = fabs(TWO_PI / atan2 (en.im, en.re));
      if (count >= 1000 || fabs(count - round(count)) > eps)  {
         cout << "Count of symmetry axis is " << count << ". Too strange...\n";
         SX_EXIT;
      }
      type.axisCount = int(lround(count));
   }

   type.opCoord = Coord::toVec3Ref(&eig.vecs(3*id));

   return type;
}


size_t SxRotation::hash (const SxRotation &in)
{
   // --- data: double m[3][3];
   return SxHashFunction::hash (const_cast<double*>(in.m[0]), 9);
}
