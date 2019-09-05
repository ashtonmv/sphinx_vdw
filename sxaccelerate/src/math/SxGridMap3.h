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

#ifndef _SX_GRID_MAP_3_H_
#define _SX_GRID_MAP_3_H_

#include <SxArray.h>
#include <SxUtil.h>
#include <SxVector3.h>

/** \brief 3D Grid Map

    \b SxGridMap3 = SFHIngX Linear Grid Map

    3D grid map is mapped to the SxArray. It is designed to be
    used: create once - read many times as fast as possible.
    A map needs to be rebuild to add/remove a new element after
    the map is already created.

    The grid structure maps always only points (0D).
    The type of points can vary with parameter T.

    To build a structure from the type Coord (SxVector3<Double>):
\code
   SxGridMap3<Double> atomicStructureGrid;
\endcode

    \par Example:
\code
   // --- the coordinates of the same species
   SxArray<SxVector3<Float> > coords;
   ... fill the array

   // --- create the grid structure with the size of cell 6.8 (1.8A)
   SxGridMap3<Float> grid;
   coords = grid.build (coords, 6.8f);

   // --- 1. primary selection
   // --- get all cells within the distance 6.8 (1.8A)
   // --- from the atom position (0,0,0)
   // --- =1 cell + 26 neighbouring cells
   SxArray<ssize_t> nodes;
   SxVector3<Float> atom (0.0f, 0.0f, 0.0f);
   nodes = grid.intersect (atom);

   // --- 2. secondary selection
   // --- print all bonds within specified distance lower than 6.8 (1.8A)
   // --- from the atom position (0,0,0)
   ssize_t iNode;
   ssize_t nNodes = nodes.getSize ();
   ssize_t iAtom, iAtomFrom, iAtomTo;

   SxVector3<Float> distance;

   // --- covalent radius of the atom
   float colvalentRadius1 = 1.0f;

   // --- covalent radius of all atoms in the grid
   float colvalentRadius2 = 0.5f;

   for (iNode=0; iNode < nNodes; ++iNode)  {
      // --- get the interval from the current node
      iAtomFrom = grid.getFrom (nodes(iNode));
      iAtomTo   = grid.getTo (nodes(iNode));
      for (iAtom=iAtomFrom; iAtom<=iAtomTo; ++iAtom)  {
         // --- check the distance
         // --- specified distance: colvalentRadius1 + colvalentRadius2
         distance = coords(iAtom) - atom;
         if (distance.norm () < colvalentRadius1 + colvalentRadius2)  {
            // --- print the coordinates
            cout << coords(iAtom) << endl;
         }
      }
   }
\endcode

    SxGridMap3 can be seen as the bottom level of SxOctreeMap.
    \sa SxOctreeMap

    \author Vaclav Bubnik, bubnik@mpie.de */
template<class T>
class SX_EXPORT_UTIL SxGridMap3
{
   public:
      /** Constructor. */
      SxGridMap3 ();

      /** Destructor. */
      ~SxGridMap3 ();

      // --- intervals
      ssize_t getFrom (ssize_t iCell_) const;
      ssize_t getTo (ssize_t iCell_) const;
      ssize_t getSize () const;

      // --- cells
      const SxVector3<Int> &getDim () const { return dim; }

      // --- bounding box
      const SxVector3<T> &getMin () const { return min; }
      const SxVector3<T> &getMax () const { return max; }
      
      SxVector3<Int> getCell (const SxVector3<T> &coord_) const;

      // --- creation
      SxArray<ssize_t> getPermutation  (const SxArray<SxVector3<T> > &coords_,
                                        typename T::Type             cell_);
                                        
      SxArray<SxVector3<T> > build (const SxArray<SxVector3<T> > &coords_,
                                    typename T::Type             cell_);

      // --- intersection
      // coordinate
      void intersect (const SxVector3<T> &coord_,
                      SxArray<ssize_t>   *nodes_,
                      ssize_t            *nNodes_) const;
                      
      // min-max box
      void intersect (const SxVector3<T> &min_,
                      const SxVector3<T> &max_,
                      SxArray<ssize_t>   *nodes_,
                      ssize_t            *nNodes_) const;

      // sphere
      void intersect (const SxVector3<T> &coord_,
                      typename T::Type   radius_,
                      SxArray<ssize_t>   *nodes_,
                      ssize_t            *nNodes_) const;


      SxArray<ssize_t> intersect (const SxVector3<T> &coord_) const;
      SxArray<ssize_t> intersect (const SxVector3<T> &min_,
                                  const SxVector3<T> &max_) const;
      SxArray<ssize_t> intersect (const SxVector3<T> &coord_,
                                  typename T::Type   radius_) const;

   protected:
      // --- cells
      typename T::Type cellSize;
      SxVector3<Int> dim;

      // --- intervals;
      SxArray<ssize_t> from;
      SxArray<ssize_t> to;

      // --- bounding box
      SxVector3<T> min;
      SxVector3<T> max;
};

// Default constructor.
template<class T>
SxGridMap3<T>::SxGridMap3 ()
   : cellSize (0)
{
   // empty
}

// Destructor.
template<class T>
SxGridMap3<T>::~SxGridMap3 ()
{
   // empty
}


template<class T>
ssize_t SxGridMap3<T>::getFrom (ssize_t iCell_) const
{
   return from(iCell_);
}


template<class T>
ssize_t SxGridMap3<T>::getTo (ssize_t iCell_) const
{
   return to(iCell_);
}


template<class T>
ssize_t SxGridMap3<T>::getSize () const
{
   return from.getSize ();
}


template<class T>
SxVector3<Int> SxGridMap3<T>::getCell (const SxVector3<T> &coord_) const
{
   SxVector3<Int> result;
   SxVector3<T> coord = coord_ - min;
   
   result(0) = 1 + (ssize_t)(coord(0) / cellSize);
   result(1) = 1 + (ssize_t)(coord(1) / cellSize);
   result(2) = 1 + (ssize_t)(coord(2) / cellSize);
   
   return result;
}


template<class T>
SxArray<ssize_t> SxGridMap3<T>::getPermutation
   (const SxArray<SxVector3<T> >  &coords_,
    typename T::Type              cell_)
{
   ssize_t iCoord;
   const ssize_t nCoords = coords_.getSize ();
   SxArray<ssize_t> lookupTable (nCoords);
   
   if (nCoords < 1)  {
      // --- abortion, an empty grid can not be created
      dim = SxVector3<Int> (0, 0, 0);
      return lookupTable;
   }
   // --- at least one point...

   // --- get bounding box
  // --- the first point is both min and max
   SxVector3<T> pos;
   min = coords_ (0);
   max = coords_ (0);
   for (iCoord=0; iCoord < nCoords; ++iCoord)  {
      pos = coords_ (iCoord);
      if (pos(0) < min(0)) min(0) = pos(0);
      if (pos(1) < min(1)) min(1) = pos(1);
      if (pos(2) < min(2)) min(2) = pos(2);
      if (pos(0) > max(0)) max(0) = pos(0);
      if (pos(1) > max(1)) max(1) = pos(1);
      if (pos(2) > max(2)) max(2) = pos(2);
   }

   cellSize = cell_;
   SX_CHECK (cellSize > 0, cellSize);

   SxVector3<T> d = max - min;
   dim(0) = 3 + (ssize_t)(d(0) / cellSize);
   dim(1) = 3 + (ssize_t)(d(1) / cellSize);
   dim(2) = 3 + (ssize_t)(d(2) / cellSize);
   ssize_t nCell = dim.product ();
   ssize_t dim01 = dim(0)*dim(1);

   SxArray<ssize_t> intervals (nCell);
   intervals.set (0);

   from.resize (nCell);
   to.resize (nCell);

   // --- 1. pass: count occupation
   ssize_t x,y,z;
   SxVector3<T> coord;
   for (iCoord=0; iCoord < nCoords; ++iCoord)  {
      coord = coords_(iCoord) - min;
      x = 1 + (ssize_t)(coord(0) / cellSize);
      y = 1 + (ssize_t)(coord(1) / cellSize);
      z = 1 + (ssize_t)(coord(2) / cellSize);
      intervals(x+y*dim(0)+z*dim01) += 1;
   }

   // --- prepare mapping from grid to linear array
   ssize_t iSeg=0;
   ssize_t iTo=0;
   ssize_t idx=0;
   ssize_t iCell;
   for (iCell=0; iCell < nCell; iCell++)  {
      from(idx)      = iSeg;
      iTo            = intervals(idx);
      to(idx)        = iSeg + iTo-1;
      iSeg          += iTo;
      intervals(idx) = 0;
      idx++;
   }

   // --- 2. pass: build permutation lookup table
   // ---          to sort points by the cells
   for (iCoord=0; iCoord < nCoords; ++iCoord)  {
      // --- find index for the point sorted by the cells
      coord = coords_(iCoord) - min;
      x = 1 + (ssize_t)(coord(0) / cellSize);
      y = 1 + (ssize_t)(coord(1) / cellSize);
      z = 1 + (ssize_t)(coord(2) / cellSize);
      idx = x + y*dim(0) + z*dim01;

      //lookupTable(iCoord) = from(idx) + intervals(idx);
      lookupTable(from(idx) + intervals(idx)) = iCoord;
      intervals(idx) += 1;
   }

   return lookupTable;
}


template<class T>
SxArray<SxVector3<T> > SxGridMap3<T>::build
   (const SxArray<SxVector3<T> > &coords_,
    typename T::Type cell_)
{
   SxArray<ssize_t> lookupTable = getPermutation (coords_, cell_);
   
   // --- permutation   
   SxArray<SxVector3<T> > coordsSorted (coords_);
   if (coordsSorted.getSize () > 1)  {
      coordsSorted.sortByIdx (lookupTable);
   }

   return coordsSorted;
}

// coordinate
template<class T>
void SxGridMap3<T>::intersect (const SxVector3<T> &coord_,
                               SxArray<ssize_t>   *nodes_,
                               ssize_t            *nNodes_) const
{
   SX_CHECK (nodes_);
   SX_CHECK (nNodes_);

   ssize_t nNodes=0;
   ssize_t x,y,z,idx;
   ssize_t dim01;

   SxVector3<T> coord = coord_ - min;
   x = 1 + (ssize_t)(coord(0) / cellSize);
   y = 1 + (ssize_t)(coord(1) / cellSize);
   z = 1 + (ssize_t)(coord(2) / cellSize);

   dim01 = dim(0)*dim(1);
   idx = x + y*dim(0) + z*dim01;

   if (x > 0 && x < dim(0)-1)  {
   if (y > 0 && y < dim(1)-1)  {
   if (z > 0 && z < dim(2)-1)  {
      // --- insert 27 cells
      (*nodes_)(nNodes++) = idx;

      (*nodes_)(nNodes++) = idx + 1;
      (*nodes_)(nNodes++) = idx + 1 - dim01;
      (*nodes_)(nNodes++) = idx + 1 + dim01;
      (*nodes_)(nNodes++) = idx + 1 + dim(0);
      (*nodes_)(nNodes++) = idx + 1 + dim(0) - dim01;
      (*nodes_)(nNodes++) = idx + 1 + dim(0) + dim01;
      (*nodes_)(nNodes++) = idx + 1 - dim(0);
      (*nodes_)(nNodes++) = idx + 1 - dim(0) - dim01;
      (*nodes_)(nNodes++) = idx + 1 - dim(0) + dim01;

      (*nodes_)(nNodes++) = idx - 1;
      (*nodes_)(nNodes++) = idx - 1 - dim01;
      (*nodes_)(nNodes++) = idx - 1 + dim01;
      (*nodes_)(nNodes++) = idx - 1 + dim(0);
      (*nodes_)(nNodes++) = idx - 1 + dim(0) - dim01;
      (*nodes_)(nNodes++) = idx - 1 + dim(0) + dim01;
      (*nodes_)(nNodes++) = idx - 1 - dim(0);
      (*nodes_)(nNodes++) = idx - 1 - dim(0) - dim01;
      (*nodes_)(nNodes++) = idx - 1 - dim(0) + dim01;

      (*nodes_)(nNodes++) = idx - dim01;          // (x) + (y)*dim(0) + (z-1)*dim01;
      (*nodes_)(nNodes++) = idx + dim01;          // (x) + (y)*dim(0) + (z+1)*dim01;
      (*nodes_)(nNodes++) = idx + dim(0);         // (x) + (y+1)*dim(0) + (z)*dim01;
      (*nodes_)(nNodes++) = idx + dim(0) - dim01; // (x) + (y+1)*dim(0) + (z-1)*dim01;
      (*nodes_)(nNodes++) = idx + dim(0) + dim01; // (x) + (y+1)*dim(0) + (z+1)*dim01;
      (*nodes_)(nNodes++) = idx - dim(0);         // (x) + (y-1)*dim(0) + (z)*dim01;
      (*nodes_)(nNodes++) = idx - dim(0) - dim01; // (x) + (y-1)*dim(0) + (z-1)*dim01;
      (*nodes_)(nNodes++) = idx - dim(0) + dim01; // (x) + (y-1)*dim(0) + (z+1)*dim01;
   } // end for z
   } // end for y
   } // end for x

   *nNodes_ = nNodes;
}

// min-max box
template<class T>
void SxGridMap3<T>::intersect (const SxVector3<T> &min_,
                               const SxVector3<T> &max_,
                               SxArray<ssize_t>   *nodes_,
                               ssize_t            *nNodes_) const
{
   SX_CHECK (nodes_);
   SX_CHECK (nNodes_);

   ssize_t nNodes=0;
   ssize_t x, y, z;
   ssize_t idx;
   ssize_t dim01;
   
   // --- get [min, max] interval in cells
   SxVector3<Int> minCell = getCell (min_);
   SxVector3<Int> maxCell = getCell (max_);

   dim01 = dim(0)*dim(1);
      
#if 1
   // --- fast rejection
   (*nNodes_) = 0;
   if (max_(0) < min(0)) return;
   if (max_(1) < min(1)) return;
   if (max_(2) < min(2)) return;
   if (min_(0) > max(0)) return;
   if (min_(1) > max(1)) return;
   if (min_(2) > max(2)) return;

   // --- cut
   if (minCell(0) < 1) minCell(0) = 1;
   if (minCell(1) < 1) minCell(1) = 1;
   if (minCell(2) < 1) minCell(2) = 1;
   if (maxCell(0) > dim(0)-2) maxCell(0) = dim(0)-2;
   if (maxCell(1) > dim(1)-2) maxCell(1) = dim(1)-2;
   if (maxCell(2) > dim(2)-2) maxCell(2) = dim(2)-2;

   // --- append all valid cels (cut) from the interval [minCell, maxCell]
   for (x=minCell(0); x <= maxCell(0); ++x)  {
      for (y=minCell(1); y <= maxCell(1); ++y)  {
         for (z=minCell(2); z <= maxCell(2); ++z)  {
            idx = x + y*dim(0) + z*dim01;
            (*nodes_)(nNodes++) = idx;
         }
      }
   }
#else
   // --- append all valid cels from the interval [minCell, maxCell]
   for (x=minCell(0); x <= maxCell(0); ++x)  {
      if (x > 0 && x < dim(0)-1)  {
         for (y=minCell(1); y <= maxCell(1); ++y)  {
            if (y > 0 && y < dim(1)-1)  {
               for (z=minCell(2); z <= maxCell(2); ++z)  {
                  if (z > 0 && z < dim(2)-1)  {
                     idx = x + y*dim(0) + z*dim01;
                     (*nodes_)(nNodes++) = idx;
                  }
               } // z
            }
         } // y
      }
   } // x
#endif

   *nNodes_ = nNodes;
}

// sphere
template<class T>
void SxGridMap3<T>::intersect (const SxVector3<T> &coord_,
                               typename T::Type   radius_,
                               SxArray<ssize_t>   *nodes_,
                               ssize_t            *nNodes_) const
{
   SxVector3<T> r (radius_, radius_, radius_); // half diagonal
   SxVector3<T> minCoord (coord_ - r);
   SxVector3<T> maxCoord (coord_ + r);
   
   intersect (minCoord, maxCoord, nodes_, nNodes_);
}

// coordinate
template<class T>
SxArray<ssize_t> SxGridMap3<T>::intersect (const SxVector3<T> &coord_) const
{
   SxArray<ssize_t> nodes (from.getSize ());
   ssize_t nNodes=0;

   intersect (coord_, nodes, nNodes);

   nodes.resize (nNodes, true);
   return nodes;
}

// min-max box
template<class T>
SxArray<ssize_t> SxGridMap3<T>::intersect (const SxVector3<T> &min_,
                                           const SxVector3<T> &max_) const
{
   SxArray<ssize_t> nodes (from.getSize ());
   ssize_t nNodes=0;

   intersect (min_, max_, nodes, nNodes);

   nodes.resize (nNodes, true);
   return nodes;
}

// sphere
template<class T>
SxArray<ssize_t> SxGridMap3<T>::intersect (const SxVector3<T> &coord_,
                                           typename T::Type   radius_) const
{
   SxArray<ssize_t> nodes (from.getSize ());
   ssize_t nNodes=0;

   intersect (coord_, radius_, nodes, nNodes);

   nodes.resize (nNodes, true);
   return nodes;
}

#endif /* _SX_GRID_MAP_3_H_ */
