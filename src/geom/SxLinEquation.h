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

#ifndef _SX_LIN_EQUATION_H_
#define _SX_LIN_EQUATION_H_

#include <SxAtomicStructure.h>
#include <SxGeom.h>

/**
  @brief Linear equations in 3D space
  */
class SX_EXPORT_GEOM SxLinEquation 
{
   protected:
      /// Whether equation system is soluble
      bool soluble;
      /// Coefficients
      SxList<Coord> leftSide;
      /// Right side vector
      SxList<PrecTauR> rightSide;
   public:
      /// Constructor
      SxLinEquation (const SxMatrix3<TPrecTauR> &mat, const Coord &res);
      /// Constructor
      SxLinEquation () : soluble (true) {}

      /// Add equation
      void addEquation (const Coord &coeff, PrecTauR res);

      /// Get number of equations
      int getSize () const { return int(leftSide.getSize ()); }

      /// Concat two equation systems
      void operator&= (const SxLinEquation& linEq);

      /// Concat two equation systems
      SxLinEquation operator&& (const SxLinEquation &eq2) const;

      /// Assignment
      SxLinEquation& operator= (const SxLinEquation &in);

      /// Comparison
      bool operator== (const SxLinEquation &eq2);

      /// Is system soluble?
      bool isSoluble () const { return soluble; }

      /// A coefficient vector c is considered zero if c.absSqr () < epsZero
      static PrecTauR epsZero;

      /// verbal description of the solution
      SxString getName (const SxString &solName = "point") const;

      /**
        Result depends on dimensionality of the solution space
        - 0 (point) -- get the point
        - 1 (line)  -- get the base point of the origin on that line
        - 2 (plane) -- get the base point of the origin on that plane
        - 3 (space) -- fail (SX_EXIT)
        - insoluble -- fail (SX_EXIT)
        
        @brief Get offset point
      */
      Coord getPoint () const;

      /// Return solution space dimension (-1 for insoluble equations)
      int getDimension () const
      {
         return (soluble ? (3 - getSize ()) : -1);
      }

};
      


#endif /* _SX_LIN_EQUATION_H_ */
