// ---------------------------------------------------------------------------
//
//           The general purpose cross platform C/C++ framework
//
//                       S x A c c e l e r a t e
//
//           Home:       https://www.sxlib.de
//           License:    Apache 2
//           Authors:    see src/AUTHORS
//
// ---------------------------------------------------------------------------

#ifndef _SX_SPHERE_GRID_H_
#define _SX_SPHERE_GRID_H_

#include <SxMath.h>
#include <SxVector3.h>
#include <SxVector.h>
#include <SxMatrix.h>

/** \brief Spherical grids a la Lebedev-Laikov

  References (according to van Wuellens code on www.ccl.net)

   [1] V.I. Lebedev, and D.N. Laikov
       "A quadrature formula for the sphere of the 131st
        algebraic order of accuracy"
       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.

   [2] V.I. Lebedev
       "A quadrature formula for the sphere of 59th algebraic
        order of accuracy"
       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 

   [3] V.I. Lebedev, and A.L. Skorokhodov
       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 

   [4] V.I. Lebedev
       "Spherical quadrature formulas exact to orders 25-29"
       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 

   [5] V.I. Lebedev
       "Quadratures on a sphere"
       Computational Mathematics and Mathematical Physics, Vol. 16,
       1976, pp. 10-24. 

   [6] V.I. Lebedev
       "Values of the nodes and weights of ninth to seventeenth 
        order Gauss-Markov quadrature formulae invariant under the
        octahedron group with inversion"
       Computational Mathematics and Mathematical Physics, Vol. 15,
       1975, pp. 44-51.

       The grid data is taken from a C code found on the internet, which
       probably is the original code by Dmitrij Laikov, or derived from it. 

    \author Christoph Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_MATH SxSphereGrid
{
   public:
      /// Empty constructor
      SxSphereGrid () { /* empty */ }
      
      /// Available grids
      enum GridType {
         Grid_6, Grid_14, Grid_26, Grid_38, Grid_50, Grid_74, Grid_86,
         Grid_110, Grid_146, Grid_170, Grid_194, Grid_230, Grid_266, Grid_302,
         Grid_350, Grid_434, Grid_590, Grid_770, Grid_974, Grid_1202,
         Grid_1454, Grid_1730, Grid_2030, Grid_2354, Grid_2702, Grid_3074,
         Grid_3470, Grid_3890, Grid_4334, Grid_4802, Grid_5294, Grid_5810
      };

      /// Grid weights
      SxVector<Double> weights;

      /// Grid vectors (:iom, :idir)
      SxMatrix<Double> xyz;

      /// Return grid size
      ssize_t getSize () const
      {
         return weights.getSize ();
      }

      /// Create grid
      void create (GridType type);

      /// Constructor
      SxSphereGrid (GridType type)
      {
         create (type);
      }

      /// Spherical harmonics on angular grid
      SxMatrix<Double> Ylm;
      /// Vector spherical harmonics Psi on angular grid (iOm:idir + 3*ilm)
      SxMatrix<Double> Psilm;

      /// Compute real spherical harmonics
      void setupYlm (int lmax);

      /// Compute real spherical harmonics
      void setupYlmPsilm (int lmax);

      /// Constructor
      SxSphereGrid (GridType type, int lmax)
      {
         create (type);
         setupYlm (lmax);
      }

      /// Get xyz for point i
      SxVector3<Double> getXyz(int i)  {
         SX_CHECK (i >= 0 && i < xyz.nRows (), i, xyz.nRows ());
         return SxVector3<Double> (xyz(i,0), xyz(i,1), xyz(i,2));
      }

   private:
      /// Current index (during setup only)
      int idx;

      /// Resize grid
      void resize (int n)
      {
         weights.resize (n);
         xyz = SxMatrix<Double> (n, 3);
         idx = 0;
      }

      /// Add a new point
      inline void addPoint (const SxVector3<Double> &x)
      {
         SX_CHECK (idx >= 0 && idx < xyz.nRows (), idx, xyz.nRows ());
         xyz(idx, 0) = x(0);
         xyz(idx, 1) = x(1);
         xyz(idx, 2) = x(2);
         ++idx;
      }

      ///@{
      /** The naming is motivated by special points on a cube
          (corner, edge center, diag)
          or of the cartesian coordinate system (axis, plane).
        */
      /// Generate 6 points on the main axes
      void generateAxis (double w);
      /// Generate 12 points on the edge centers
      void generateEdgeCenter (double w);
      /// Generate 8 points on the corners
      void generateCorner (double w);
      /// Generate 24 points on a face diagonal plane
      void generateDiag (double w, double a);
      /// Generate 24 points on axis face
      void generatePlane (double w, double a);
      /// Generate 48 general points
      void generateGeneral (double w, double a, double b);
      //@}

   public:
      /// \name Vector permutation auxiliaries
      /// @{
      /// Change vector component signs according to bitmap "signs" 
      static SxVector3<Double> sign (const SxVector3<Double> &x, int signs);
      /// Generate (1,0,0) and permutations (i < 3)
      static SxVector3<Double> permute1 (int i);
      /// Generate 1/sqrt(2) (1,1,0) and permutations (i < 3)
      static SxVector3<Double> permute2 (int i);
      /// Generate (a,a,b) and permutations (2a^2 + b^2 = 1) (i < 3)
      static SxVector3<Double> permute2 (int i, double a);
      /// Generate (a,b,0) and permutations (a^2 + b^2 = 1) (i < 6)
      static SxVector3<Double> permute3 (int i, double a);
      /// Generate (a,b,c) and permutations (a^2 + b^2 + c^2 = 1) (i < 6)
      static SxVector3<Double> permute3 (int i, double a, double b);
      //@}

};

#endif /* _SX_SPHERE_GRID_H_ */
