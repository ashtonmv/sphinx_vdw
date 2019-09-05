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

#ifndef _SX_GRID_H_
#define _SX_GRID_H_

#include <SxCell.h>
#include <SxMesh3D.h>
#include <SxGeom.h>

class SxAtomicStructure;

/** \brief Partitioning into smaller grid cells

    \b SxGrid = S/PHI/nX subcell partitioning for atomic structures

    This class computes a subcell partitioning for periodic
    structures. For this, the (regularized) periodic unit cell
    is partitioned into a regular grid of subcells. For each atom in the
    structure, its grid cell id is computed. For each subcell,
    a list of atoms indices contained in this subcell is set up.
    This index list provides quick access to atoms lying in a
    certain region of space. 
    
    SxGrid is used in the find algorithm of SxAtomicStructure
    as well as for setting up neighbor structures. New applications
    are welcome, but please study the aforementioned algorithms
    carefully to learn how to use SxGrid in practice.

    \note The subcell partitioning is valid only for the
    structure provided at setup time. If you change the structure,
    you must update the grid (by assigning a new grid).

    \author Christoph Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_GEOM SxGrid
{
   public:
      /// The full cell (regularized)
      SxCell regCell;

      /// Number of grid cells 
      SxMesh3D gridMesh;

      /// The grid cells
      SxCell gridCell;

      /// Number of atoms per grid cell
      SxVector<Int> nAtPerGridCell;

      /// A copy of the structure, mapped into the regular cell
      SxArray<Coord> mappedStr;
   protected:
      /** \brief Grid cell Offsets for idx
        \note Grid cell iGrid starts at idx(offset(iGrid))
        */
      SxVector<Int> offset;
      /// Atom indices
      SxVector<Int>  idx;
   public:
      /// Constructor from structure and partitioning mesh
      SxGrid (const SxAtomicStructure &str,
              const SxVector3<Int> &gridMeshIn);

      /** \brief Constructor from structure and grid tuning factor
          @param str              atomic structure
          @param atomsPerGridCell targeted average number of atoms per grid cell

          @note atomsPerGridCell should be roughly the number of atoms
                you will try to retrieve from the grid in your algorithm.
                E.g., for structure matching and finding atoms a
                number close to 1 is appropiate (3 seems a good choice).
                For finding neighbours, 10 or 20 are appropriate.
                These numbers are only approximate and serve as hints
                only. Performance can critically depend on this number.
        */
      SxGrid (const SxAtomicStructure &str, int atomsPerGridCell);

      /** \empty Empty constructor
         \note Later initialization via assignment:
         \code
SxGrid grid;
...
grid = SxGrid(structure, 3);
         \endcode
        */
      SxGrid () : gridMesh(0,0,0) { /* empty */ }
   protected:
      /// Internal setup routine containing the algorithm
      void setup (const SxAtomicStructure &);
   public:

      /** \brief Get iterator for idx for specified grid cell
          \sa computeItEnd
        */
      SxVector<Int>::Iterator getIdxIt(int iGrid) const
      {
         SX_CHECK (idx.getSize () > 0);
         SxVector<Int>::Iterator res = idx.begin ();
         res += offset(iGrid);
         return res;
      }
      /** \brief Get iterator end for specified grid cell
        \example
        \code
SxGrid grid(...);
SxVector<Int>::Iterator idxIt  = grid.getIdxIt(iGrid),
                        idxEnd = grid.computeItEnd(iGrid);
for (; idxIt != idxEnd; idxIt++)  {
   ...
}
        \endcode
      */
      SxVector<Int>::Iterator computeItEnd(int iGrid) const
      {
         SX_CHECK (idx.getSize () > 0);
         SxVector<Int>::Iterator res = idx.begin ();
         res += offset(iGrid) + nAtPerGridCell(iGrid);
         return res;
      }

      /** \brief Suggest a grid mesh with a quite compact gridCell
          @param minSize desired (minimum) size of the mesh
        */
      static SxVector3<Int> suggestMesh(const SxCell &cell, int minSize);
};

#endif /* _SX_GRID_H_ */
