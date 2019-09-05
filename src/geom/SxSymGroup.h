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

#ifndef _SX_SYM_GROUP_H_
#define _SX_SYM_GROUP_H_
#include <SxSymOp.h>
#include <SxGeom.h>
#include <SxMesh3D.h>
#include <SxNArray.h>


/// A cell basis. The basis vectors are stored columnwise
typedef SxMatrix3<Double>  CellMat;
/// A real space rotation matrix
typedef SxMatrix3<Double>  SymMat;

// --- forward declarations (we still reside quite deep)
class SxCell;
class SxBinIO;
class SxSymbolTable;

/** \brief Symmetry transformation (affine) group

    \b SxSymGroup = S/PHI/nX affine symmetry operation group

    This class provides an symmetry group. Each symmetry is

    \f[
      x' = S ^ x + t
    \f]

    consisting of a rotational part S and a translational part t.
    Note that first the rotation and then the translation is applied.

    When t is zero, this is a symmorphic symmetrie. Symmorphic 
    symmetries are stored as the first nSymmorphic symmetries

    \author Christoph Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_GEOM SxSymGroup
{
   protected:
      /// Generating set of primitive symmetry operations
      SxArray<SxSymOp> primOps;
      /// Generating primitive cell
      CellMat primCell;
      /// Generating mesh
      SxMesh3D mesh;

      /// Number of symmorphic symmetries
      int nSymmorphic;
   public:

      /// Empty constructor
      SxSymGroup () : nSymmorphic(-1) { /* empty */}
      
      /** \brief Symmorphic symmetry only constructor
          @param symmorphic Pure rotational symmetries
          @param cell       primitive cell
          @param meshIn     optional supercell (relative to cell)
        */
      SxSymGroup (const SxArray<SymMat> &symmorphic, 
                  const CellMat &cell = CellMat(0,0,0,0,0,0,0,0,0),
                  const SxMesh3D &meshIn = SxMesh3D (1,1,1));

      /** \brief Constructor from primitive set of operations & supercell
          @param primCellIn primitive cell
          @param primOpsIn  primitive cell operations
          @param cell       optional supercell
      */
      SxSymGroup (const SxCell          &primCellIn,
                  const SxArray<SxSymOp> &primOpsIn, 
                  const SxCell          &cell);//= SxCell ()???
      
      /// leave only primitive symmetries
      void setToPrimitive ()  { mesh = SxMesh3D (1,1,1); }
      
      /// return primCell
      CellMat getPrimCell () const { return primCell; }

      /// Return symmetry operation
      SxSymOp operator()(int iSym) const;

      // get quotient (group); does not include Identity
      SxArray<SxSymOp> operator/ (const SxArray<SxSymOp> &kernel) const;
                       
      /// Return rotational part of symmetry operation:redundant,faster!
      inline const SymMat &getRot(int iSym) const {
         SX_CHECK (iSym >=0 && iSym < getSize ());
         return primOps(iSym % primOps.getSize ()).rot;
      }

      /// Get total number of symmetries
      inline int getSize () const 
      { return int(primOps.getSize () * mesh.getSize ()); } 

      /// Get number of primitive symmetries
      inline int getNPrimitive () const { return int(primOps.getSize ()); } 

      /// Get number of symmorphic symmetries
      inline int getNSymmorphic () const { return nSymmorphic; }
      
      /// Get array of symmorphic rotations
      SxArray<SymMat> getSymmorphic () const;

      /// Get one symmorphic rotation
      const SymMat &getSymmorphic(int iSym) const
      {
         SX_CHECK(iSym >= 0 && iSym < nSymmorphic, iSym, nSymmorphic);
         return primOps(iSym).rot;
      }

      /// \brief Read symmorphic symmetries from file
      void read(const SxBinIO &io, const SxCell &cell);

      /// \brief Write symmorphic symmetries to file
      void write(SxBinIO &io) const;

      /// \brief Printout routine
      void print () const;

      /** \brief Print out in sx input format
          @param output output file stream
          @param nonsymmorphic if false, print only symmorphic symmetries
        */
      void fprintsx (FILE *output, bool nonsymmorphic = false) const;

      /// \brief Check internal consistency
      void check () const;

      /** \brief Set up symmetry multiplication table
        @return table

        @note S(i) * S(j) = S(k) => table(i,j) = k
        */
      static SxArray2<int> getSymMulTable (const SxArray<SymMat> &syms,
                                           double epsSym = 1e-5);

      void read (const SxSymbolTable *table);
};

#endif /* _SX_SYM_GROUP_H_ */
