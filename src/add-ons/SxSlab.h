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

#ifndef _SX_SLAB_H_
#define _SX_SLAB_H_

#include <SxMatrix3.h>
#include <SxVector3.h>
#include <SxList.h>
#include <SxTypes.h>
#include <SxSortedList.h>
#include <SxAtomicStructure.h>
#include <SxExt.h>


/** @brief Generate slab geometry from bulk geometry

    \b SxSlab = S/PHI/nX Slab Generator

    ....
    @ingroup  group_addons
    @author   Christoph Freysoldt, freyso@fhi-berlin.mpg.de */
class SX_EXPORT_EXT SxSlab
{
   public:

   /// Surface selection of slab (for saturation)
   enum Direction {None, Upward, Downward, Both};
      
   /// Constructor
   SxSlab ();
   /// Destructor
   ~SxSlab () { /* empty */}

   /// Bulk structure
   SxAtomicStructure bulkStructure;

   /// Compute new slab structure
   SxAtomicStructure compute (Direction dir,
                 int nLayer,
                 int nVacLayer,
                 bool zOrthogonal,
                 double threshold,
                 double bondLength,
                 double minDist,
                 const SxString & chemSymbol);

   /// Miller indices of plane
   RelVec hkl;
   /// Slab cell
   SxCell newCell;
   /// Slab structure
   SxAtomicStructure slabStructure;

   /// Set bulk cell
   void setBulk  (const SxAtomicStructure &);
   /// Set plane (Miller indices)
   void setPlane (const RelVec &);
   /// Set plane (Miller indices)
   void setPlane (int, int, int);

   /// Auxiliary functions
   //@{
   
   /// Divide int vector by greatest common divisor
   static bool shortenVector (SxVector3<Int>*);
   /// Largest common divisor a la Euklid
   static int getGcd (int x, int y);
   /// Integer Gram-Schmidt orthogonalization
   static void integerGramSchmidt(RelVec *v1Ptr, 
                                  RelVec *v2Ptr,
                                  const CellMat &cell);
   /// Set minVec to newVec if newVec is shorter, or minVec has 0. length
   /// @param newVec    vector to be tested
   /// @param minVecPtr pointer to minVec
   /// @brief Store shortest non-zero vector
   static void storeShortest(const SxVector3<Double> &newVec,
                                   SxVector3<Double> *minVecPtr);
   /// Fill cell with lattice points
   static SxList<RelVec> getLatticePointsInCell (const SxMatrix3<Int>&);
   /** 
     The factor is determined such that the resulting matrix is
     integer.
     @param m matrix to be inverted
     @param factor pointer to int, which is set to the multiplication factor
     @brief Get inverse of integer matrix multiplied by factor
     */
   static SxMatrix3<Int> intInverse (const SxMatrix3<Int> m, int *factor);
   //@}
   /** Find a cell where the first two basis vectors are orthogonal to hkl
       and call findElementaryCell afterwards
       @brief Get hkl adapted bulk cell
   */
   void findCell ();
   /** 
     @param relCell a cell, that has the right (b1, b2) plane
     @brief Find elementary cell in cell spanned by relCell, keep (b1,b2) plane
    */
   void findElementaryCell (const SxMatrix3<Int> &relCell);
   /// Rotate newCell such that x || c1 and z _|_ (c1, c2)
   void rotateCell ();
   /**
     @param nLayer number of layers in slab
     @param nVacLayer number of vacuum layers
     @param zOrthogonal if true, set lateral component of c3 axis to 0.
     @brief Set newCell and newTau
     */
   void getSlab (int nLayer, int nVacLayer, bool zOrthogonal = false);
   /**
     @param nLayer number of layers in slab
     @param nVacLayer number of vacuum layers
     @param zOrthogonal if true, set lateral component of c3 axis to 0.
     @param dir which direction to saturate
     @param bondThreshold criterium how to identify broken bonds
     @param newBondLength bond length to saturation atom
     @param identicalThreshold minimum distance between saturation atoms,
            average coordinates of closer atoms
     
     Uses getNextNeighbours to find broken bonds. According to the
     criteria, saturation atoms are introduced to saturate these bonds.
     @brief Set newCell and newTau, saturate
     */
   void getSaturatedSlab (int nLayer, int nVacLayer,
                          bool zOrthogonal,
                          Direction dir,
                          double bondThreshold,
                          double newBondLength,
                          double identicalThreshold);

   /// Rotation matrix that rotates (hkl)-plane into xy-plane
   CellMat axisRotation;

};


#endif /* _SX_SLAB_H_ */
