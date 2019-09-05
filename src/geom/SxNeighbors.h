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

#ifndef _SX_NEIGHBORS_H_
#define _SX_NEIGHBORS_H_

#include <SxAtomicStructure.h>
#include <SxGeom.h>
#include <SxGrid.h>

/** \brief Atomic structure next neighbors container

    \b SxNeighbors = S/PHI/nX 

    

    \author Christoph Freysoldt, freysoldt@mpie.de */

class SX_EXPORT_GEOM SxNeighbors
{
   public:
      /// Bitmaps for setup modes
      enum SetupMode {
         /// Store indices
         StoreIdx             = 0x01,
         /// Store absolute positions (no atomInfo)
         StoreAbsNoInfo       = 0x02,
         /// Store absolute positions
         StoreAbs             = 0x03,
         /// Store relative positions (no atomInfo)
         StoreRelNoInfo       = 0x04,
         /// Store relative positions
         StoreRel             = 0x05,
         /// Store absolute and relative positions (no atomInfo)
         StoreAbsAndRelNoInfo = 0x06,
         /// Store absolute and relative positions
         StoreAbsAndRel       = 0x07,
         /// Store squared distances
         StoreDistSqr         = 0x08,
         /// Store all information (positions (abs+rel) + distances)
         StoreAll             = 0x0f,
         /// Include zero distance neighbor
         IncludeZeroDistance  = 0x10
      };
         
      /// Absolute positions of neighbors (if set up)
      SxAtomicStructure absPositions;
      /// Relative positions of neighbors (if set up)
      SxAtomicStructure relPositions;
      /// Atom indices
      SxVector<Int> idx;
      /// Squared distances
      SxVector<TPrecTauR> distSqr;

      /// Empty constructor
      SxNeighbors () {}
      /** \brief Constructor
        @param grid      A subcell partitioning of str
        @param str       the atoms
        @param refPoint  reference point
        @param rcut      cutoff radius
        @param mode      what to do

        */
      SxNeighbors (const SxGrid            &grid,
                   const SxAtomicStructure &str,
                   const Coord             &r0,
                   double                  rcut,
                   int                     mode)
      {
         compute (grid, str, r0, rcut, mode);
      }

      /** \brief Constructor
        @param grid      A subcell partitioning of str
        @param str       the atoms
        @param refPoint  reference point
        @param mode      what to do (see 

        */
      void compute (const SxGrid            &grid,
                    const SxAtomicStructure &str,
                    const Coord             &ref,
                    double                  r0,
                    int                     mode);

      /// Return number of neighbors
      int getSize () const;

};

#endif /* _SX_NEIGHBORS_H_ */
