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

#ifndef _SX_LATTICE_SHELLS_H_
#define _SX_LATTICE_SHELLS_H_

#include <SxDFT.h>
#include <SxConfig.h>
#include <SxTypes.h>
#include <SxCell.h>
#include <SxArray.h>

/** \brief Finds and stores shell-wise the lattice vectors belonging to
           neighboured unit cells

    \b SxLatticeShells = S/PHI/nX accessing the different shells of neighboured
                         unit cells

    \author Matthias Wahn, wahn@fhi-berlin.mpg.de */
class SX_EXPORT_DFT SxLatticeShells
{
   public:
      /** lattice vector of the unit cell, relative to which neighbourhood
          is defined here */
      Coord referenceVec;

      /** how many shells (levels of neighbourhood: 1st nearest, 2nd nearest,
          and so on ...) shall be considered */
      int nShells;

      /** empty constructor */
      SxLatticeShells ();

      /** constructor */
      SxLatticeShells (const SxCell &cell_, int nShells_,
                       const Coord &referenceVec_=Coord(0.));
      
      /** destructor */
      ~SxLatticeShells ()  { /* empty */ };

      /** get i-th nearest neighbours */
      SxArray<Coord > getShell (int i) const;

      /** get array of all neighbours */
      SxArray<Coord > getAllNeighbVecs () const;

      /** get a lattice vector from the array of all neighboured lattice
          vectors by index */
      Coord           getNeighbVec (int idx) const;

      /** the i-th element of the return is the shell of the i-th element
          in the return of ::getAllNeighbVecs */
      SxArray<int>    getAllShells () const;

      /** get the shell of the i-th entry in the array of all neighbours */
      int             getShellOfIdx (int idx) const;

      /** get number of all neighboured lattice vectors in total */
      int             getNAllVecs () const;

      /** find all neighboured lattice vectors up to the specified shell */
      void            find ();

   protected:
      /** shell-wise storage of the neighboured lattice vectors */
      SxArray<SxList<Coord > > shells;

      /** all neighboured lattice vectors as a long array */
      SxArray<Coord > allNeighbourVecs;

      /** the shells of all neighboured lattice vectors as a long array */
      SxArray<int>   allShells;

      /** the unit cell */
      SxCell cell;

      /** dump a list of SxVector3<Double>'s */
      void printShell (const SxList<Coord> &list) const;

      /** dump an array of SxVector3<Double>'s */
      void printShell (const SxArray<Coord> &array) const;
};

#endif /* _SX_LATTICE_SHELLS_H_ */
