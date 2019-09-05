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

#ifndef _SX_CELL_H_
#define _SX_CELL_H_

#include <SxMatrix3.h>
#include <SxGeom.h>
#include <SxBinIO.h>
#include <SxSymbolTable.h>
#include <SxPtr.h>
#include <SxSymGroup.h>


/// A real space coordinate (or vector)
typedef SxVector3<Double>  Coord;
/// Lattice vector in relative coordinates
typedef SxVector3<Int>     RelVec;

/** \brief Periodic cell

  This class is derived from a 3x3 matrix, which contains columnwise
  the cell basis. You should not work on this matrix directly, but
  rather use #set to reset the cell.
  Alternatively, whenever you changed the cell, call #setup to update 
  the related variables.

  \ingroup  group_structure
  \author   C. Freysoldt,   freyso@fhi-berlin.mpg.de
  */
class SX_EXPORT_GEOM SxCell : public CellMat
{
   public:
   
      /**
        Note that this is the stored result of Cell::inverse ().
        @brief Inverse of cell
        */
      CellMat inv;
      /// Cell volume
      PrecTauR volume;
      /// Basis vectors
      Coord basis (int i) const { return CellMat::col(i); }
      /** \brief Precision for regularity

        If a basis vector j is as short as possible, the vector
        \f$b_i \pm b_j\f$ must be at least as long:
        \f[ ||b_i \pm b_j||^2 = ||b_i||^2 \pm 2(b_i \cdot b_j) + ||b_j||^2 
            \geq ||b_j||^2 \f]
        \f[ ||b_i||^2 \pm  2(b_i \cdot b_j) \geq 0 \f]

        In order to allow for numerical noise in the left side, we replace
        the 0 on the right side by -#epsRegular.  */
      PrecTauR epsRegular;
#     define EPS_REGULAR_DEFAULT  1e-4
      /** \brief Precision in symmetry searching

        \f[\left|(U^T\cdot U - E)_{ij}\right| < \varepsilon\f]
        */
      PrecTauR epsSym;
#     define EPS_SYM_DEFAULT      1e-5
      
      /** \brief Container for symmetries (NOT set up by SxCell) [cartesian]

        Typically, the atoms in the periodic cell reduce the symmetries
        that the empty cell offers. In order to determine these symmetries,
        the atomic coordinates must be known. Once they have been determined,
        they may be stored in SxCell::symGroup.
        \example
        \code
SxCell cell;
SxAtomicStructure structure;
// --- set up atoms
...
// set periodicity in structure
structure.cell = cell;
// find symmetries compatible with atomic structure and store them in cell
cell.symGroupPtr = SxPtr<SxSymGroup>::create(structure.getSymGroup ());
        \endcode
        
        If you want to store the full symmetry of the cell (ignoring 
        anything inside that may break the symmetry), do it this way
        \code
SxCell cell (...);
cell.symGroupPtr->setup(cell.getLatticeSymmetries (), cell);
        \endcode
        */
      SxPtr<SxSymGroup> symGroupPtr;
      
   public:
      ///@name Constructors and destructor
      //@{
      /// Constructor
      SxCell ();
      
      /// Constructor
      SxCell (const CellMat& inCell);
      /// Constructor
      SxCell (const CellMat& inCell, const SxPtr<SxSymGroup> &symGrp);
      /// Constructor
      SxCell (const SxCell &cell);
       
      /// Constructor from basis vectors
      SxCell (const Coord& a1, const Coord& a2, const Coord &a3);
      /// Constructor from basis vectors
      SxCell (const SxList<Coord> &list);

      /// Constructor from file
      SxCell (const SxBinIO & io);


      /// Constructor from symbol table
      SxCell (const SxSymbolTable *);

      /// Assignment operator
      SxCell& operator= (const SxCell &in);

      /// Destructor
      ~SxCell ();

      /// Sets up internal variables
      void setup ();

      /// Set cell
      void set (const CellMat &cell);
      /// Set cell from basis vectors
      void set (const Coord &a1, const Coord &a2, const Coord &a3);
      /** \brief Read cell from netcdf file
        @param io netcdf file
        */
      void read (const SxBinIO &);
      /** \brief Write to netcdf file
        @param io netcdf file
        */
      void write (SxBinIO &io) const;
      /** \brief Set from set of lattice vectors
        \param points a set of lattice points that should lie around the
               origin. If this is not the case, the routine may fail.
        \return the number of independent basis vectors. If this is less
        than 3, the cell cannot be used as a 3D cell.

        Use this routine if you have obtained a set of lattice points
        in order to derive the regular elementary cell.

        \note experimental feature, may not be reliable
        */
      int setFromGoodVectors(const SxList<Coord> &points);
      /** \brief Set from set of lattice vectors
        \param points a set of lattice points that should lie around the
               origin.
        \return the number of independent basis vectors. If this is less
        than 3, the cell cannot be used as a 3D cell.

        Use this routine if you have obtained a set of lattice points
        in order to derive the regular elementary cell.

        \note This is a wrapper routine to setFromGoodVectors with some
              attempts to improve the success rate.
        */
      int setFromVectors(const SxList<Coord> &points);


      /** \brief Compute generating mesh for supercell
          @param superCell supercell
          @return Mesh dimension

          A generating mesh provides a set of primitive lattice vectors
          that fill a given supercell. While this is simple in case the
          supercell basis vectors are parallel to the primitive cell
          basis vectors, it is non-trivial in the general case.
          
          This routine rotates the cell into one that can be used with
          the returned mesh size \f$m_i\f$ as
          \f[
          \mathbf r = \sum n_i \mathbf a_i ~\textnormal{with}~
          0 \le n_i < m_i
          \f]
          where the resulting vectors are pairwise inequivalent with
          respect to the superlattice.
          \example
          \code
SxCell primCell, supercell, genMesh;
...
genMesh = primCell;
SxVector3<Int> mesh = genMesh.computeGen (supercell);
SxVector3<Int> n;
for (n(0) = 0; n(0) < mesh(0); n(0)++)
   for (n(1) = 0; n(1) < mesh(1); n(1)++)
      for (n(2) = 0; n(2) < mesh(2); n(2)++)
      {
        lVec = genMesh ^ n;
        ...
      }
         \endcode

         \par
         Algorithm:
         The idea is to determine how many times we can go along the different
         primitive cell directions before it becomes equivalent to a supercell
         lattice vector. This is then the mesh size along that 
         direction. For 2nd and 3rd directions, that linear combination with the
         previous one(s) must be chosen that minimizes the multiplicity.
        */
      SxVector3<Int> computeGen(const SxCell &superCell);
      //@}

   private:
      /** \brief Find smallest factor n that makes (n * v) % d == (0,0,0)

          \note This is a simultaneous division of 0 in the d-module
          It is used in computeGen, but will hardly find any other
          reasonable application, which is why I made it private.
        */
      static int moduleDiv3 (const SxVector3<Int> &v, int d);
   public:

      class SX_EXPORT_GEOM CellType  {
         /** This class detects the lattice type and computes the lattice
             constants.

             \author Matthias Wahn, wahn@fhi-berlin.mpg.de
          */
         public:
            /** Enumeration type for the lattice type */
            enum Type  {
               SimpleCubic,
               Tetragonal,
               Orthorhombic,
               FaceCenteredCubic,
               BodyCenteredCubic,
               Hexagonal,
               Unknown
            };

            /**@name contructors and destructor */
            //@{
            /** empty constructor */
            CellType ();

            /** standard constructor */
            CellType (const SxCell &cell);

            /** destructor */
            ~CellType ()  { /* empty */ }
            //@}

            /**@name get functions

               In case of a bcc or an fcc structure, the lattice constant
               \f$ a_{\rm lat} \f$ does not correspond to the length of a
               vector of the unit cell. For this reason, ::getA yields
               the lattice constant \f$ a_{\rm lat} \f$ for fcc or bcc
               structures and not the length of the 1st unit cell vector.

               ::getCosAB, ::getCosBC, and ::getCosAC yield in any case the
               cosines of the angles between the unit cell vectors.
             */
            //@{
            /** yields the type of the lattice structure */
            Type getType () const;

            /** returns the name of the lattice structure */
            SxString getTypeName () const;

            double getA () const;
            double getB () const;
            double getC () const;

            double getCosAB () const;
            double getCosBC () const;
            double getCosAC () const;
            //@}

            /** prints parameters found */
            void print () const;

         protected:
            /**@name the lattice constants and lattice angles */
            //@{
            double a, b, c;
            double bcCos, acCos, abCos;
            //@}

            /** lattice type */
            Type type;

            /** type name, e.g., "face centered cubic" */
            SxString typeName;
      };

      /** \brief Returns the reciprocal cell

          This function computes the reciprocal unit cell vectors
          \f$ \mathbf{b}_1 \f$, \f$ \mathbf{b}_2 \f$, and \f$ \mathbf{b}_3\f$ 
          stored column-wisely in the returned matrix. 
          The reciprocal unit vectors read
          \f[
             \mathbf{b}_1 
                = 2\pi
                   \frac{\mathbf{a}_2 \times \mathbf{a}_3}
                        {\mathbf{a}_1 \cdot (\mathbf{a}_2 \times \mathbf{a}_3)},
          \f]
          \f[
             \mathbf{b}_2 
                = 2\pi
                   \frac{\mathbf{a}_3 \times \mathbf{a}_1}
                        {\mathbf{a}_1 \cdot (\mathbf{a}_2 \times \mathbf{a}_3)},
          \f]
          \f[
             \mathbf{b}_3 
                = 2\pi
                   \frac{\mathbf{a}_1 \times \mathbf{a}_2}
                        {\mathbf{a}_1 \cdot (\mathbf{a}_2 \times \mathbf{a}_3)}.
          \f]
          Here, \f$\mathbf{b}_1\f$, \f$\mathbf{b}_2\f$, and \f$\mathbf{b}_3 \f$ 
          refer to the lattice vectors (columns of this SxCell object).

          Do not mix up the reciprocal cell with the inverse. The
          reciprocal cell \b B is the \e transposed of SxCell::inv
          times \f$ 2\pi \f$! */
      SxCell getReciprocalCell () const;

      /** \brief Get transformation to the conventional cubic cell

        \note: The present algorithm relies on the cubic lattice
        symmetries. Other cells than cubic are currently not covered.
          
        */
      SxMatrix3<Int> getConventional () const;

      /** \brief Get rotation into the standard orientation (if possible)

           @note: this relies on CellType identification, not all
                  variants of a lattice type are correctly identified
        */
      SxMatrix3<Double> getStandardRot () const;

      /// @name Change between relative and absolute
      //@{
      /// Calculate relative coordinate from cartesian
      inline Coord carToRel (const Coord& cartesian) const;
      /// Calculate cartesian coordinate from relative
      inline Coord relToCar (const Coord& relative) const;
      /// Calculate relative transformation from cartesian
      inline SymMat carToRel (const SymMat& cartesian) const;
      /// Calculate cartesian transformation from relative
      inline SymMat relToCar (const SymMat& relative) const;
      
      /// Change from cartesian to relative coordinates
      inline void changeToRel (Coord* cartesian) const;
      /// Change from relative to cartesian coordinates 
      inline void changeToCar (Coord* relative) const;
      /// Change transformation from cartesian to relative coordinates
      inline void changeToRel (SymMat* cartesian) const;
      /// Change transformation from relative to cartesian coordinates 
      inline void changeToCar (SymMat* relative) const;
      //@}

      /// @name Mapping routines
      //@{
      /**
          Mapping means that a lattice vector is added to the coordinate
          such that the results fulfills certain conditions.
          This corresponds to adding integer numbers to the relative
          coordinates (which may be fractional).
          
        
        \brief Possible modes of mapping routines
        */
      enum Mapping { 
         /// make relative coordinates between 0 and 1
         Positive, 
         /// make relative coordinates between -1/2 and 1/2
         Origin, 
         /** \brief map into Wigner-Seitz cell (not yet implemented)

          This yields the shortest vector that can be obtained.
          */
         WignerSeitz };

      /** \brief Return mapped coordinate
          @param in the coordinate to be mapped
          @param mode mapping mode
          @return mapped coordinate

       */
      inline 
      Coord getMapped (const Coord &in, enum Mapping mode = Positive) const;

      /** \brief Return mapped coordinate
          @param in the coordinate to be mapped
          @param mode mapping mode
          @return mapped coordinate

       */
      inline 
      Coord getMappedRel (const Coord &in, enum Mapping mode = Positive) const;
      
      /** \brief Map cartesian coordinate into cell
        @param mode mapping mode, see #Mapping
        */
      inline void map (Coord *, enum Mapping mode = Positive) const;

      /** \brief Map relative coordinate into cell
        @param mode mapping mode, see #Mapping
        */
      inline void mapRel (Coord *, enum Mapping mode = Positive) const;

   protected:
      /** \brief Auxiliary data for Wigner-Seitz mapping
        */
      CellMat reg, regInv;
   public:
      //@}

      

      /// @name Distances and lengthes
      //@{
      /// Real length of a relative vector
      PrecTauR lengthFromRel (const Coord &rel) const;
      /// Real absolute square of a relative vector
      PrecTauR lengthSqFromRel (const Coord &rel) const;
      
      /// Shortest distance between two cartesian coords
      PrecTauR shortestDist (const Coord &cart1, const Coord &cart2) const;

      /** \brief  get heights between parallel planes

          \return The \f$i\f$-th coordinate denotes the height between
                  the plane spanned by the lattice vectors
                  \f$\mathbf{a}_{j}\f$ and \f$\mathbf{a}_{k}\f$ and its
                  parallel plane, with \f$i\f$, \f$j\f$, \f$k\f$ cyclic.

          \author Matthias Wahn
          */
      Coord getHeights () const;

      /** \brief Get size of bounding box in relative coords
        @param rCut radius of sphere inside box
        @return extent of bounding box

        The bounding box has the following property: any vector shorter
        than the cutoff has relative coordinates, whose absolute value
        is smaller than the bounding box.

        This allows to loop over all meshpoints possibly inside the
        cutoff sphere via
        \code
Coord bBox = meshCell.getBoundingBox (rCut);
SxVector3<Int> from = ceil(r0Rel-bBox), to = floor(r0Rel+bBox);
        \endcode
        */
      Coord getBoundingBox (double rCut) const;
      //@}

      /// @name Regular cells
      //@{

      /** \brief Compute a regular (compact) elementary cell

        Some algorithms assume that the chosen basis of the lattice has the
        following properties:
        - the basis vectors are as short as possible
        - the angle between basis vectors is as large as possible
        
        We call this a regular cell. 
       */
      inline SxCell getRegularCell () const
      {
         SxCell res(getRegularLattice(*this, epsRegular));
         res.epsSym = epsSym;
         res.epsRegular = epsRegular;
         return res;
      }
      /// \brief Computational routine for regular cells
      static CellMat getRegularLattice (const CellMat &cell, double eps);
      //@}


      /// @name Symmetry
      //@{
      /**
       \brief Get (rotation and mirror) symmetries in cartesian coordinates

       \note Driver routine

       \note This does NOT set symGroup. This returns the full symmetry of
             the current cell.
       */
      SxList<SymMat> getLatticeSymmetries () const;

      /// Check if sym is a lattice symmetry
      bool isLatticeSymmetry (const SxMatrix3<Double> &sym) const;
   protected:
      /**
       \brief Get (rotation and mirror) symmetries in cartesian coordinates

       \note Computational routine. Cell must be regular.

       \note This does NOT set symGroup. This returns the full symmetry of
             the current cell.
       */
      SxList<SymMat> getLatticeSymmetries (bool regular) const;
      
};

/// Streaming operator
SX_EXPORT_GEOM ostream &operator<< (ostream &, const SxCell &);

Coord SxCell::getMapped (const Coord &in, enum Mapping mode) const
{
   SX_CHECK (mode == Positive || mode == Origin || mode == WignerSeitz);
   Coord res = in;
   map (&res, mode);
   return res;
}

Coord SxCell::getMappedRel (const Coord &in, enum Mapping mode) const
{
   SX_CHECK (mode == Positive || mode == Origin || mode == WignerSeitz);
   Coord res = in;
   mapRel (&res, mode);
   return res;
}

// this is implemented inline so the if (mode..) block can be compiled out
// this is not implemented via mapRel, since the explicit implementation is 
// significantly faster and mapping can be time critical (e.g. neighbor lists)
void SxCell::map (Coord *in, enum Mapping mode) const
{
   SX_CHECK (mode == Positive || mode == Origin || mode == WignerSeitz);
   if (mode == Positive)  {
      Coord rel = carToRel(*in);
      for (int i = 0; i < 3; i++)
        *in -= floor(rel(i)) * basis(i);
      
   } else if (mode == Origin)  {
      Coord rel = carToRel(*in);
      for (int i = 0; i < 3; i++)
        *in -= round(rel(i)) * basis(i);
      
   } else if (mode == WignerSeitz)  {
      // stage 1: map into regular cell
      Coord rel = regInv ^ (*in);
      for (int i = 0; i < 3; ++i) rel(i) = round(rel(i));
      *in -= reg ^ rel;
      RelVec t;
      Coord ins, res = *in;
      double n2, nrm2 = (*in).normSqr ();
      for (t(0) = -1; t(0) <=1; t(0)++)  {
         for (t(1) = -1; t(1) <=1; t(1)++)  {
            for (t(2) = -1; t(2) <=1; t(2)++)  {
               ins = (*in + (reg ^ t));
               n2 = ins.normSqr ();
               if (n2 < nrm2)  {
                  res = ins;
                  nrm2 = n2;
               }
            }
         }
      }
      *in = res;
   }
}

void SxCell::mapRel (Coord *in, enum Mapping mode) const
{
   SX_CHECK (mode == Positive || mode == Origin || mode == WignerSeitz);
   if (mode == Positive)  {
      for (int i = 0; i < 3; i++)
        (*in)(i) -= floor((*in)(i));
      
   } else if (mode == Origin)  {
      for (int i = 0; i < 3; i++)
        (*in)(i) -= round((*in)(i));
      
   } else if (mode == WignerSeitz)  {
      cout << "relative Wigner-Seitz mapping not implemented." << endl;
      // and never will, since one should work with absolute coordinates
      // in general, and in particular for Wigner-Seitz mapping
      SX_EXIT;
   }
}

inline void operator%=(Coord & x, const SxCell &cell)
{
   cell.map (&x);
}

// ----- Change between relative and cartesian ---
Coord SxCell::carToRel (const Coord &cartesian) const
{
   return inv ^ cartesian;
}

Coord SxCell::relToCar (const Coord &relative) const
{
   return (*(const CellMat *)this) ^ relative;
}

void SxCell::changeToRel (Coord* cartesian) const
{
   *cartesian = carToRel (*cartesian);
}

void SxCell::changeToCar (Coord* relative) const
{
   *relative = relToCar (*relative);
}

SymMat SxCell::carToRel (const SymMat& cartesian) const
{
   return (inv ^ cartesian ^ (*(const CellMat*)this));
}

SymMat SxCell::relToCar (const SymMat& relative) const
{
   return ((*(const CellMat*)this) ^ relative ^ inv);
}

void SxCell::changeToRel (SymMat* cartesian) const
{
   *cartesian = carToRel (*cartesian);
}

void SxCell::changeToCar (SymMat* relative) const
{
   *relative = relToCar (*relative);
}

#endif /* _SX_CELL_H_ */
