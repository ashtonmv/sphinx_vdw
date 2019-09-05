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

#ifndef _SX_OCTREE_MAP_H_
#define _SX_OCTREE_MAP_H_

#include <SxArray.h>
#include <SxVector3.h>

/** \brief Linear Octree Map

    \b SxOctreeMap = SFHIngX Linear Octree map

    3D octree map is mapped to the SxArray. It is designed to be
    used: create once - read many times as fast as possible.
    A map needs to be rebuild to add/remove a new element after
    the map is already created.

    The octree structure maps always only points (0D).
    The type of points can vary with parameter T.

    To build a structure from the type Coord (SxVector3<Double>):
\code
   SxOctreeMap<Double> atomicStructureTree;
\endcode

    \par Example:
\code
   // --- the coordinates of the same species
   SxArray<SxVector3<Float> > coords;
   ... fill the array

   // --- create the octree structure with six levels upon the coordinates
   SxOctreeMap<Float> octree;
   coords = octree.build (coords, 6);

   // --- get all nodes within specified distance
   SxArray<ssize_t> nodes;
   SxVector3<Float> atom (0.0f, 0.0f, 0.0f);

   // --- A), B), C) show various selection radii thanks to hierarchical
   // --- structure.

   // --- A)
   // --- get all nodes within the distance 100.0f from the atom position
   // --- if all atoms are inside the selected radius then only a one node will
   // --- be returned with node values: [0, coords.getSize () - 1]
   nodes = octree.intersect (atom, 100.0f);

   // --- B)
   // --- get all nodes within the distance 10.0f from the atom position
   nodes = octree.intersect (atom, 10.0f);

   // --- C)
   // --- get all nodes within the distance 1.0f from the atom position
   nodes = octree.intersect (atom, 1.0f);

\endcode

    SxOctreeMap can be seen as hierarchical SxGridMap3.
    Hence an octree is a good choice only if there are some benefits
    from its hierarchy and to use a grid otherwise.
    \sa SxGridMap3

    \par Intervals:

    The SxOctreeMap structure is not an octree by its definition since
    the values are stored in all nodes and not only in the leaves.
    The values are in the form of intervals to minimize the amount
    of information necessary to describe mapped data. A list of indices
    {1,2,3,4} mapping a part of an array which belongs to a leaf is then
    stored as an interval [1, 4].
    \li {1, 2, 3, ...., n} is stored as [1, n]
    \li {2} is stored as [2, 2]
    \li 0 is stored as [a, b] where a > b

    The intervals in nodes which are not the leaves describe the whole
    interval of their children nodes. If a node A has two subnodes
    B[2,8] and C[9,9] then the node A contain the interval A[2,9].

    \author Vaclav Bubnik, bubnik@mpie.de */
template<class T>
class SxOctreeMap
{
   public:
      /** Constructor. */
      SxOctreeMap ();

      /** Destructor. */
      ~SxOctreeMap ();

      // --- intervals maped to an array
      ssize_t getFrom (ssize_t iNode_);
      ssize_t getTo (ssize_t iNode_);

      // --- boundary
      SxVector3<T> getCenter (ssize_t iNode_);
      typename T::Type getRadius (ssize_t iNode_);
      SxVector3<T> getMin (ssize_t iNode_);
      SxVector3<T> getMax (ssize_t iNode_);

      /** Return the number of nodes.
          \return the size of octree structure. */
      ssize_t getNNodes ();

      SxArray<ssize_t> getPermutation (const SxArray<SxVector3<T> > &coords_,
                                       int nTreeLevels_);

      SxArray<SxVector3<T> > build (const SxArray<SxVector3<T> > &coords_,
                                    int nTreeLevels_);



      // --- to find bonds ...... 31.32s
      SxArray<ssize_t> intersect (const SxVector3<T> &coord_,
                                  typename T::Type radius_);

      // --- to find bonds ...... 27.25s
      void intersect (const SxVector3<T> &coord_,
                      typename T::Type radius_,
                      SxArray<ssize_t> *nodes_,
                      ssize_t *nNodes_);

      void intersect (double frustrumPlanes_[24],
                      SxArray<ssize_t> *nodes_,
                      ssize_t *nNodes_);
                      
      // --- intersection with n planes
      void intersect (const SxArray<SxVector3<Float> > &normals_,
                      const SxArray<float>             &offsets_,
                      SxArray<ssize_t>                 *nodes_,
                      ssize_t                          *nNodes_);

   protected:
      void compileCoord (int iLevel_,
                         ssize_t idx_,
                         const SxVector3<T> &coord_,
                         SxVector3<T> min_,
                         SxVector3<T> max_);

      int nTreeLevels;

      // --- intervals;
      SxArray<ssize_t> from;
      SxArray<ssize_t> to;
      SxArray<ssize_t> length;

      // --- auxiliary stack allows iteration instead of recursion
      SxArray<ssize_t> stack;

      // --- in fact these four arrays can be computed
      // --- any time from just one bounding box
      // --- in the root level
      // --- the idea is then to save CPU time:
      // --- a) min+max/2 to get the center
      // --- b) max-min/2 to get half diagonal
      // --- c) sqrt or ^ to get radius or square radius
      //
      // --- The values are not initialized for empty cells.
      SxArray<SxVector3<T> >    minBox;
      SxArray<SxVector3<T> >    maxBox;
      SxArray<SxVector3<T> >    centers;
      SxArray<typename T::Type> radii;
};

// Default constructor.
template<class T>
SxOctreeMap<T>::SxOctreeMap ()
   : nTreeLevels (0)
{
   // empty
}

// Destructor.
template<class T>
SxOctreeMap<T>::~SxOctreeMap ()
{
   // empty
}


template<class T>
ssize_t SxOctreeMap<T>::getFrom (ssize_t iNode_)
{
   return from(iNode_);
}


template<class T>
ssize_t SxOctreeMap<T>::getTo (ssize_t iNode_)
{
   return to(iNode_);
}


template<class T>
SxVector3<T> SxOctreeMap<T>::getCenter (ssize_t iNode_)
{
   return centers(iNode_);
}


template<class T>
typename T::Type SxOctreeMap<T>::getRadius (ssize_t iNode_)
{
   return radii(iNode_);
}


template<class T>
SxVector3<T> SxOctreeMap<T>::getMin (ssize_t iNode_)
{
   return minBox(iNode_);
}


template<class T>
SxVector3<T> SxOctreeMap<T>::getMax (ssize_t iNode_)
{
   return maxBox(iNode_);
}


template<class T>
ssize_t SxOctreeMap<T>::getNNodes ()
{
   return from.getSize ();
}


template<class T>
SxArray<ssize_t> SxOctreeMap<T>::getPermutation
   (const SxArray<SxVector3<T> > &coords_,
    int nTreeLevels_)
{
   ssize_t iCoord;
   const ssize_t nCoords = coords_.getSize ();

   nTreeLevels = nTreeLevels_;

   if (nTreeLevels < 1)  {
      // --- at least one level == root
      nTreeLevels = 1;
   }

   // --- nTreeLevels:  0  1   2   3    4
   // ---    treeSize:  0  1   9  73   585
   // ---                          v
   // ---     example:            73 = 1 + 8 + 64
   int iTree;
   int treeSize=1;
   int shift=1;
   for (iTree=1; iTree < nTreeLevels; iTree++)  {
      shift = shift<<3;
      treeSize += shift;
   }

   length.resize (treeSize);
   length.set (0);
   from.resize (treeSize);
   to.resize (treeSize);

   centers   = SxArray<SxVector3<T> > (treeSize);
   radii     = SxArray<typename T::Type> (treeSize);
   stack     = SxArray<ssize_t> (treeSize);
   minBox    = SxArray<SxVector3<Float> > (treeSize);
   maxBox    = SxArray<SxVector3<Float> > (treeSize);

   // --- get bounding box
   SxVector3<T> min;
   SxVector3<T> max;
   ssize_t i;

   if (nCoords > 0)  {
      // --- at least one point
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
   }

   // --- 1. pass: count occupation
   for (iCoord=0; iCoord < nCoords; ++iCoord)  {
      compileCoord (0, 0, coords_(iCoord), min, max);
   }

   // --- prepare mapping from octree to linear array
   ssize_t nSeg=1;
   ssize_t iSeg=0;
   ssize_t iTo=0;
   ssize_t idx=0;
   ssize_t offset=0;
   for (iTree=0; iTree < nTreeLevels; iTree++)  {
      iSeg=0;
      for (i=0; i < nSeg; ++i)  {
         idx         = i + offset;
         from(idx)   = iSeg;
         iTo         = length(idx);
         to(idx)     = iSeg + iTo-1;
         iSeg        += iTo;
         length(idx) = 0;
      }
      offset += nSeg;
      nSeg*=8;
   }

   // --- 2. pass: build permutation lookup table
   // ---          to sort points by the octants
   SxArray<ssize_t> lookupTable (coords_.getSize ());
   SxVector3<T> coord;
   SxVector3<T> center;
   ssize_t a, b, c;

   for (iCoord=0; iCoord < nCoords; ++iCoord)  {
      // --- find index sorted by octants for the point
      coord  = coords_(iCoord);

      // --- traverse the octree
      idx = 0;
      for (iTree=0; iTree < nTreeLevels-1; iTree++)  {
         // --- the center of the current octant
         center = centers(idx);

         // --- determine next sub-node octant
         a = (coord(0) > center(0)) ? 4 : 0; // 2^2
         b = (coord(1) > center(1)) ? 2 : 0; // 2^1
         c = (coord(2) > center(2)) ? 1 : 0; // 2^0

         // ---  base  + octant
         idx = 8*idx+1 + a+b+c;
      }

      lookupTable(from(idx) + length(idx)) = iCoord;
      length(idx) += 1;
   }


   return lookupTable;
}


template<class T>
SxArray<SxVector3<T> > SxOctreeMap<T>::build
   (const SxArray<SxVector3<T> > &coords_,
    int nTreeLevels_)
{
   SxArray<ssize_t> lookupTable = getPermutation (coords_, nTreeLevels_);

   // --- permutation
   SxArray<SxVector3<T> > coordsSorted (coords_);
   if (coordsSorted.getSize () > 1)  {
      coordsSorted.sortByIdx (lookupTable);
   }

   return coordsSorted;
}


template<class T>
SxArray<ssize_t> SxOctreeMap<T>::intersect (const SxVector3<T> &coord_,
                                            typename T::Type radius_)
{
   SxArray<ssize_t> nodes (from.getSize ());
   ssize_t nNodes=0;

   intersect (coord_, radius_, &nodes, &nNodes);

   nodes.resize (nNodes, true);
   return nodes;
}


template<class T>
void SxOctreeMap<T>::intersect (const SxVector3<T> &coord_,
                                typename T::Type radius_,
                                SxArray<ssize_t> *nodes_,
                                ssize_t *nNodes_)
{
   SX_CHECK (nodes_);
   SX_CHECK (nNodes_);

   ssize_t nNodes=0;

   int i;
   SxVector3<T> c; // octant center
   SxVector3<T> distance;
   typename T::Type d; // octant radius
   typename T::Type r;

   ssize_t nIndices = from.getSize (); // number of octree nodes
   ssize_t itIdx = 0; // current node, start from root (==0)
   ssize_t nextIdx; // next node
   //SxArray<ssize_t> stack (nIndices); // stack with node indices
   ssize_t *sp = stack.elements; // the top of the stack

   SX_CHECK (sp != NULL);
   if (stack.getSize () < 1)  {
      *nNodes_ = nNodes;
      return;
   }

   // --- push the root node to the stack
   *sp++ = itIdx;

   // --- test for endless loop
   SX_CHECK (sp - 1 == stack.elements);

   while (sp != stack.elements)  {
      // --- recover index from the stack -> recursion
      sp -= 1;
      itIdx = sp[0];

      if (to(itIdx) - from(itIdx) < 0)  {
         // --- skip an empty node
         continue;
      }

      // --- octant center and radius
      c = centers(itIdx);
      r = radii(itIdx);

      distance = centers(itIdx) - coord_;
      d = distance.norm ();

      if (d < radius_ - radii(itIdx))  {
         // --- completely in
         (*nodes_)(nNodes++) = itIdx;
      }  else  {
         if (d < radius_ + radii(itIdx))  {
            // --- the next index = subnode, the first child node
            // --- from the current level, the nodes and leafs are
            // --- mapped to the array so we are using indices
            // --- instead of pointers.
            nextIdx = 8*itIdx + 1;

            // --- simmilar to (level < treeLevels)
            // --- we do not need to save actual levels to the stack
            if (nextIdx < nIndices)  {
               // --- partialy in -> division to suboctants...
               // --- push the suboctants indices to the stack
               for (i=0; i < 8; i++)  {
                  *sp++ = nextIdx+i;
               }
            }  else  {
               // --- partialy inside
               (*nodes_)(nNodes++) = itIdx;
            }
         }
         else  {
            // --- completely out
         }
      }
   }

   *nNodes_ = nNodes;
}


template<class T>
void SxOctreeMap<T>::intersect (double frustrumPlanes_[24],
                                SxArray<ssize_t> *nodes_,
                                ssize_t *nNodes_)
{
   SX_CHECK (nodes_);
   SX_CHECK (nNodes_);

   ssize_t nNodes=0;

   int i;
   SxVector3<T> c; // octant center
   SxVector3<T> distance;
   double d; // the distance between node and frustrum plane
   int positive; // the number of nodes completely inside frustrum
   int partial; // the number of nodes partialy inside the camera frustrum
   typename T::Type r; // octant radius

   ssize_t nIndices = from.getSize (); // number of octree nodes
   ssize_t itIdx = 0; // current node, start from root (==0)
   ssize_t nextIdx; // next node
   //SxArray<ssize_t> stack (nIndices); // stack with node indices
   ssize_t *sp = stack.elements; // the top of the stack

   SX_CHECK (sp != NULL);
   if (stack.getSize () < 1)  {
      *nNodes_ = nNodes;
      return;
   }

   // --- push the root node to the stack
   *sp++ = itIdx;

   // --- test for endless loop
   SX_CHECK (sp - 1 == stack.elements);

   while (sp != stack.elements)  {
      // --- recover index from the stack -> recursion
      sp -= 1;
      itIdx = sp[0];

      if (to(itIdx) - from(itIdx) < 0)  {
         // --- skip an empty node
         continue;
      }

      // --- reset in/out states
      positive = 0;
      partial = 0;

      // --- octant center and radius
      c = centers(itIdx);
      r = radii(itIdx);

      // --- test octant sphere envelop with 6 frustrum planes
      // --- faster then to test all 8 octant points
      for (i=0; i < 6; i++)  {
         // --- distance from the frustrum plane
         d = frustrumPlanes_[i*4]*c(0) +
               frustrumPlanes_[i*4+1]*c(1) +
               frustrumPlanes_[i*4+2]*c(2) +
               frustrumPlanes_[i*4+3];

         if (d > r)  {
            // --- octant is completely in for the current plane
            positive++;
         }

         if (d + r > 0.0)  {
            // --- octant intersects with the current plane
            partial++;
         }
      }

      if (positive == 6)  {
         // --- completely in -> draw
         // --- the vertex culling can be disabled (extension)
         (*nodes_)(nNodes++) = itIdx;
      }  else  {
         if (partial == 6)  {
            // --- the next index = subnode, the first child node
            // --- from the current level, the nodes and leafs are
            // --- mapped to the array so we are using indices
            // --- instead of pointers.
            nextIdx = 8*itIdx + 1;

            // --- simmilar to (level < treeLevels)
            // --- we do not need to save actual levels to the stack
            if (nextIdx < nIndices)  {
               // --- partialy in -> division to suboctants...
               // --- push the suboctants indices to the stack
               for (i=0; i < 8; i++)  {
                  *sp++ = nextIdx+i;
               }
            }  else  {
               // --- partialy inside -> draw
               (*nodes_)(nNodes++) = itIdx;
            }
         }
      }
   }

   *nNodes_ = nNodes;
}


template<class T>
void SxOctreeMap<T>::intersect (const SxArray<SxVector3<Float> > &normals_,
                                const SxArray<float>             &offsets_,
                                SxArray<ssize_t>                 *nodes_,
                                ssize_t                          *nNodes_)
{
   SX_CHECK (nodes_);
   SX_CHECK (nNodes_);
   SX_CHECK (normals_.getSize () == offsets_.getSize (),
             normals_.getSize (), offsets_.getSize ());

   ssize_t nNodes=0;

   int i;
   SxVector3<T> c; // octant center
   SxVector3<T> distance;
   double d; // the distance between node and frustrum plane
   ssize_t positive; // the number of nodes completely inside
   ssize_t partial; // the number of nodes partialy inside
   typename T::Type r; // octant radius
   
   ssize_t nPlanes = normals_.getSize ();

   ssize_t nIndices = from.getSize (); // number of octree nodes
   ssize_t itIdx = 0; // current node, start from root (==0)
   ssize_t nextIdx; // next node
   ssize_t *sp = stack.elements; // the top of the stack

   SX_CHECK (sp != NULL);
   if (stack.getSize () < 1)  {
      *nNodes_ = nNodes;
      return;
   }

   // --- push the root node to the stack
   *sp++ = itIdx;

   // --- test for endless loop
   SX_CHECK (sp - 1 == stack.elements);

   while (sp != stack.elements)  {
      // --- recover index from the stack -> recursion
      sp -= 1;
      itIdx = sp[0];

      if (to(itIdx) - from(itIdx) < 0)  {
         // --- skip an empty node
         continue;
      }

      // --- reset in/out states
      positive = 0;
      partial = 0;

      // --- octant center and radius
      c = centers(itIdx);
      r = radii(itIdx);

      // --- test octant sphere envelop with the planes
      for (i=0; i < nPlanes; i++)  {
         // --- distance from the frustrum plane
         d = normals_(i)(0)*c(0) +
             normals_(i)(1)*c(1) +
             normals_(i)(2)*c(2) +
             offsets_(i);

         if (d > r)  {
            // --- octant is completely in for the current plane
            positive++;
         }

         if (d + r > 0.0)  {
            // --- octant intersects with the current plane
            partial++;
         }
      }

      if (positive == nPlanes)  {
         // --- completely in
         (*nodes_)(nNodes++) = itIdx;
      }  else  {
         if (partial == nPlanes)  {
            nextIdx = 8*itIdx + 1;

            if (nextIdx < nIndices)  {
               // --- partialy in -> division to suboctants...
               for (i=0; i < 8; i++)  {
                  *sp++ = nextIdx+i;
               }
            }  else  {
               // --- partialy inside
               (*nodes_)(nNodes++) = itIdx;
            }
         }
      }
   }

   *nNodes_ = nNodes;
}


template<class T>
void SxOctreeMap<T>::compileCoord (int iLevel_,
                                   ssize_t idx_,
                                   const SxVector3<T> &coord_,
                                   SxVector3<T> min_,
                                   SxVector3<T> max_)
{
   iLevel_++;
   SxVector3<T> r = ((max_ - min_) * 0.5);
   SxVector3<T> center = min_ + r;

   length(idx_) += 1;
   centers(idx_) = center;
   typename T::Type radius = r.norm ();
   radii(idx_) = radius;
   minBox(idx_) = min_;
   maxBox(idx_) = max_;

   if (iLevel_ < nTreeLevels)  {
      ssize_t a=0;
      ssize_t b=0;
      ssize_t c=0;

      if (coord_(0) > center(0))  {
          min_(0) = center(0);
          a = 4;
      }  else  {
          max_(0) = center(0);
      }

      if (coord_(1) > center(1))  {
         min_(1) = center(1);
         b = 2;
      }  else  {
         max_(1) = center(1);
      }

      if (coord_(2) > center(2))  {
         min_(2) = center(2);
         c = 1;
      }  else  {
         max_(2) = center(2);
      }
      ssize_t iSegment = a+b+c;
      idx_ = 8*idx_+1;

      return compileCoord (iLevel_, idx_+iSegment, coord_, min_, max_);
   }
}

#endif /* _SX_OCTREE_MAP_H_ */
