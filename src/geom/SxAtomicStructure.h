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

#ifndef _SX_ATOMIC_STRUCTURE_H_
#define _SX_ATOMIC_STRUCTURE_H_

#include <SxCell.h>
#include <SxMatrix.h>
#include <SxPrecision.h>
#include <SxSymbolTable.h>
#include <SxSpeciesData.h>
#include <SxAtomInfo.h>
#include <SxSymGroup.h>
#include <SxGeom.h>
#include <SxLoopMPI.h>

class SxGrid;

// say that Coord is a simple type (no destructor, memcpy ok)
SXSTACK_SIMPLE(Coord)

/**
  This class contains the structural information of a periodic system
  - cell (lattice periodicity)
  - atom positions (absolute, cartesian)
  - metainformation about the atoms (species...)

  \author C. Freysoldt, freyso@fhi-berlin.mpg.de
  \brief    Atomic structure class
  \ingroup  group_structure
  */

class SX_EXPORT_GEOM SxAtomicStructure
{
   public:
      /// The coordinates (absolute, cartesian) as 3 x n matrix
      SxMatrix<TPrecTauR> coords;

      /// Number of atoms
      int nTlAtoms;

      /// Number of unique species
      int nSpecies;

      /// The periodic boundary conditions
      SxCell cell;

      /** \brief Pointer to species information.

        This points to the species information, which may be shared
        among different atomic structures.

        If atomInfo is set to NULL (or its SxPtr equivalent), species
        information is not used.

        In most cases (see below for exceptions), binary operations on
        two atomic structures require that they share atomInfo (or have
        no species information at all).

        If two atomic structures
        are generated independently, but refer to the same atoms, the
        atomInfo of the different structures should be set to the same
        value explicitly:
        \code
SxAtomicStructure str1, str2;
// generate str1 and str2 independently
...
// let str2's atomInfo share that of str1
str2.replaceInfo (str1.atomInfo);
       \endcode

       Note that atomInfo may also store relations between different,
       but related structures, e.g. a subset of a structure or a
       different ordering of the atoms.
       Some operators (+=,-=,<<=) exploit this information and allow
       to operate on structures with different atomInfos.

       \sa
       - SxAtomInfo
       - examples/sxatomstructderiv.cpp.
        */
      SxConstPtr<SxAtomInfo> atomInfo;

      /// @name Construction
      //@{
      /// Possible internal states
      enum ConstructionMode {
         /// Uninitialized
         Unknown,
         /// Under construction: coordinates in lists
         ListConstruction,
         /// Ready for usage
         Finished
      };

      enum MetadataTypes {
         Labels,
         StickyFilter,
         Elements
      };

   protected:
      /// Temporary atomic positions in creation mode
      SxList<SxList<Coord > > posList;

      /// Internal status of the construction
      enum ConstructionMode creationStatus;

      //@}
   public:
      /// @name Constructors
      //@{
      /// Possible ways to initialize from existing SxAtomicStructure
      enum DataHandling {
         /// Copy the coordinates
         Copy,
         /// Reference the coordinates (default)
         Reference
      };

      /// Constructor: no setup
      SxAtomicStructure ();
      /// Constructor: set only cell, add atoms later
      explicit SxAtomicStructure (const SxCell &inCell);
      /// Constructor from existing SxAtomicStructure
      SxAtomicStructure (const SxAtomicStructure &in,
                         enum DataHandling = Reference);
      /// Constructor from symbol table
      explicit SxAtomicStructure (const SxSymbolTable* table,
                                  const SxString &coordName = "coords");
      /// Constructor from complete atom list
      SxAtomicStructure (const SxCell &inCell,
                         const SxList<SxList<Coord> > &atoms);
      /// Constructor for pure coordinate lists (no cell, 1 species)
      SxAtomicStructure (const SxArray<Coord> &coord);
      /// Constructor for pure coordinate lists (no cell, 1 species)
      explicit SxAtomicStructure (const SxList<Coord> &coord);
      /// Constructor for pure coordinate lists (no cell, 1 species)
      explicit SxAtomicStructure (const SxVector<TPrecTauR> &coord,
                                  enum DataHandling = Reference);

      /** \brief Constructor for sized, uninitialized coordinates

          The next step must be to initialize the coordinates via
          #coords, #set, #setAtom, or #ref.
          \example
          \code
SxAtomicStructure str(cell, nAtoms);
for (int i = 0; i < nAtoms; ++i) str.ref(i) = Coord(i,i,i);
          \endcode
        */
      SxAtomicStructure (const SxCell &cellIn,
                         int  nAtoms,
                         const SxAtomInfo::ConstPtr &atomInfoIn
                               = SxAtomInfo::ConstPtr ());
      /** \brief Constructor for sized, uninitialized coordinates

          The next step must be to initialize the coordinates via
          #coords, #set, #setAtom, or #ref.
          \example
          \code
SxAtomicStructure str(nAtoms);
for (int i = 0; i < nAtoms; ++i) str.ref(i) = Coord(i,i,i);
          \endcode
        */
      explicit SxAtomicStructure (int  nAtoms,
                         const SxAtomInfo::ConstPtr &atomInfoIn
                               = SxAtomInfo::ConstPtr ());
      /** \brief Coordinate-setting constructor for derived structures
        @param atomInfoIn atom info, must be directly derived from from.atomInfo
        @param from the structure from this is derived

        This constructor collects the coordinates of the new structure from
        the parent according to the mapping defined in atomInfoIn.
        */
      SxAtomicStructure (const SxAtomInfo::ConstPtr &atomInfoIn,
                         const SxAtomicStructure &from);


      /// Constructor from binary file
      explicit SxAtomicStructure (const SxBinIO &);
      /// Destructor
      ~SxAtomicStructure () { /* empty */}
      //@}

      /// @name Setting up the structure
      //@{

      /// Bring to defined size
      void resize (int nAtoms);

      /// \brief Get a new, uninitialized structure with same size,
      /// cell & atomInfo than this
      inline SxAtomicStructure getNewStr () const;

      /// Add new atom (last species)
      void addAtom (const Coord &pos);
      /// Add new atom in given species
      void addAtom (int iSpecies, const Coord &pos);
      /// Remove atom
      void removeAtom (int iSpecies, int iAtom);

      /// Add several atoms
      void append (const SxAtomicStructure &in);

      /// New species
      void newSpecies ();
      /// Set atomic positions from list
      void copy (const SxAtomicStructure &);

      /** \brief Repeat atomic structure along lattice vectors

          The repeat function expects three integers which
          declare the number of repetitions in the direction
          of the three lattice vectors. The number of repetitions
          is 1 smaller than the received parameter, e.g. with
          na1,na2,na3 = 1, the output is the identity!
          \param na1   repetition in the direction of the \f$
                          \mathbf{a}_1
                       \f$
          \param na2   repetition in the direction of the \f$
                          \mathbf{a}_2
                       \f$
          \param na3   repetition in the direction of the \f$
                          \mathbf{a}_3
                       \f$
          \return      the repeated structure
          \author      Blazej Grabowski, blazejgrabowski@gmx.de */
      SxAtomicStructure repeat (int na1, int na2, int na3, bool symmetrize=true) const;
      /** \brief Repeat atomic structure along lattice vectors

          The repeat function expects three integers which
          declare the number of repetitions in the direction
          of the three lattice vectors. The number of repetitions
          is 1 smaller than the received parameter, e.g. with
          \f$ \mathbf{n}_\mathrm{repVec}=(1,1,1) \f$, the output
          is the identity!
          \param repVec   repetition in the direction of the
                          \f$ \mathbf{a}_i \f$ with \f$ i = x,y,z\f$
          \return      the repeated structure
          \author      Blazej Grabowski, blazejgrabowski@gmx.de */
      SxAtomicStructure repeat (const SxVector3<Int> &repVec, bool symmetrize=true ) const;
       /** \brief Repeat atomic structure in a general way

          \param repMatrix  the columns contain the new basis vectors
          \return           the repeated structure

          The repeat function expects a 3x3 matrix, the columns of
          which contain the new basis vectors in relative coordinates.
          This allows to create supercells with rotated axes.
          \example
          \code
SxAtomicStructure c2x2 = structure.repeat (SxMatrix3<Int> ( 1,1,0,
                                                           -1,1,0,
                                                            0,0,1));
          \endcode
          would create a c2x2 supercell if structure is a tetragonal
          slab structure.
          \note The SxVector3 form of the repeat function corresponds
          to a diagonal repMatrix.
          \note This is extremely slow for large repetitions due to the
                naive implementation! The repetition along unit cell
                vectors is faster.
          \todo Find a more efficient way.

         */
      SxAtomicStructure repeat (const SxMatrix3<Int> &repMatrix, bool symmetrize=true) const;

      /// Set structure from complete coordinate list
      void set (const SxList<SxList<Coord> > &atoms);
      /// Set structure from complete coordinate list
      void set (const SxVector<Double> &pos, enum DataHandling = Copy);

      /// Replace atom info
      SxAtomicStructure &replaceInfo (const SxAtomInfo::ConstPtr &newInfo)
      {
         SX_CHECK ((nTlAtoms <= 0) || newInfo);
         SX_CHECK (newInfo->nAtoms.sum () == nTlAtoms,
                   newInfo->nAtoms.sum (), nTlAtoms);
         atomInfo = newInfo;
         return *this;
      }

      /// Check if labels are available
      bool hasLabels () const
      {
         if (!atomInfo) return false;
         return atomInfo->meta.contains (Labels);
      }
      /// Request labels (crash if not available)
      const SxArray<SxString>& getLabels () const;
      /// Read elements from symbol table
      void readElements (const SxSymbolTable *table)
      {
         SX_CHECK(atomInfo);
         SxArray<SxString> chemNames = SxSpeciesData::getElements (table);
         atomInfo->meta.update (Elements, chemNames);
         SX_CHECK (chemNames.getSize () == getNSpecies (),
                   chemNames.getSize (), getNSpecies ());
      }
      /** \brief Return reference to chemical element names

        The function acquires the chemical symbols from the meta data in
        atomInfo, assuming they have been attached to the metadata.
      */
      const SxArray<SxString>& getElements () const
      {
         SX_CHECK(atomInfo);
         SX_CHECK(atomInfo->meta.contains (Elements));
         const SxArray<SxString> &elem = atomInfo->meta.get (Elements);
         SX_CHECK (elem.getSize () == getNSpecies (),
                   elem.getSize (), getNSpecies ());
         return elem;
      }

      /// Start creation mode
      void startCreation ();
      /// End creation mode
      void endCreation ();

      /// Set position of existing atom
      inline void setAtom (int iTlAtom, const Coord &pos)
      {
         ref(iTlAtom) = pos;
      }
      /// Set position of existing atom
      inline void setAtom (int iSpecies, int iAtom, const Coord &pos);
      /// Fill all coordinates with the same value
      SxAtomicStructure& set (const Coord &fillVal);
      /// Fill all coordinates with the same value
      inline void set (double x, double y, double z) { set(Coord(x,y,z)); }
      /** Import atoms from stack
          @param stack to be imported
          @param is species to be imported
        */
      void importStack(const SxStack<Coord> &stack, int is = -1);

      /**  \brief Sort atoms with ascending z, y, x
           @param derive if true, set up parent relation to unsorted structure
           \note The structure must not be in creation mode.
        */
      void sort (bool derive = false);
      //@}

      /// \name Access to coordinates and species
      //@{
      /// Get atom
      inline const Coord operator() (int iTlAtom) const
      {
         return Coord(constRef(iTlAtom));
      }
      /// Get atom
      inline const Coord operator() (int iSpecies, int iAtom) const
      {
         return Coord(constRef(iSpecies, iAtom));
      }
      /// Get atom
      inline const Coord getAtom (int iSpecies, int iAtom) const
      {
         return Coord(constRef(iSpecies, iAtom));
      }
      /// Get atom
      inline const Coord getAtom (int iTlAtom) const
      {
         return Coord(constRef(iTlAtom));
      }
      /// Get atom
      inline const Coord operator() (SxAutoLoop &iTlAtom) const
      {
         iTlAtom.setLimit (nTlAtoms);
         return Coord(constRef((int)iTlAtom.i));
      }
      /// Get atom
      inline const Coord operator() (SxAutoLoop &iSpecies, 
                                     const SxAutoLoop &iAtom) const
      {
         SX_CHECK (atomInfo);
         iSpecies.setLimit (getNSpecies ());
         iAtom.setLimit (getNAtoms ((int)iSpecies.i));
         return Coord(constRef((int)iSpecies.i, (int)iAtom.i));
      }
      /// Get atom
      inline const Coord getAtom (SxAutoLoop &iSpecies,
                                  const SxAutoLoop &iAtom) const
      {
         SX_CHECK (atomInfo);
         iSpecies.setLimit (getNSpecies ());
         iAtom.setLimit (getNAtoms ((int)iSpecies.i));
         return Coord(constRef((int)iSpecies.i, (int)iAtom.i));
      }
      /// Get atom
      inline const Coord getAtom (SxAutoLoop &iTlAtom) const
      {
         iTlAtom.setLimit (nTlAtoms);
         return Coord(constRef((int)iTlAtom.i));
      }
      /// Get total atom number
      inline int getIAtom(int iSpecies, int iAtom) const;

      /// Get total atom number (autoloop)
      template<class T1, class T2>
      inline int getIAtom(const T1 &iSpecies, const T2 &iAtom) const
      {
         SxAutoLoop::setLimit (iSpecies, getNSpecies ());
         SxAutoLoop::setLimit (iAtom, getNAtoms((int)iSpecies));
         return getIAtom((int)iSpecies,(int)iAtom);
      }

      /// Get species number of atom iAtom
      int getISpecies (int iTlAtom, int *iAtomPtr = NULL) const;

      /*
      /// Get reference to coordinates
      SxVector<TPrecTauR> coordRef ()
      {
         int nElem = coords.getSize();
         SxVector<TPrecTauR> res (coords.elements, nElem);
         res.reshape (nElem, 1);
         return res;
      }
      */

      /** \brief Get reference to coordinates from const object
          \note The return vector is still a reference to the
                atomic structure and MUST NOT BE changed. Get a
                true copy first if necessary.
                \code
SxVector<TPrecTauR> onlyRead, res;
onlyRead = structure.coordRef;
res = onlyRead * 3.; // OK, because onlyRead (and hence structure)
                     // is not affected

SxVector<TPrecTauR> degreeOfFreedom;
degreeOfFreedom.copy (structure.coordRef);
degreeOfFreedom *= 3.; // OK, only copy is affected

SxVector<TPrecTauR> buggyExample = structure.coordRef ();

buggyExample *= 3.; // WRONG! structure is changed.
                \endcode
        */
      SxVector<TPrecTauR> coordRef () const
      {
         ssize_t nElem = coords.getSize();
         SxVector<TPrecTauR> res (coords.elements, nElem);
         res.reshape (nElem, 1);
         return res;
      }

      /// Get number of species
      int getNSpecies () const {
        SX_CHECK(atomInfo);
        return (int)atomInfo->nAtoms.getSize ();
      }
      /// Get total number of atoms
      int getSize () const { return nTlAtoms; }
      /// Get total number of atoms
      int getNAtoms () const { return nTlAtoms; }
      /// Get number of atoms for certain species
      int getNAtoms (int iSpecies) const {
         SX_CHECK(atomInfo);
         return atomInfo->nAtoms(iSpecies);
      }

      /// Get atom index range for specified species
      SxIdx getRange (int iSpecies) const
      {
         SX_CHECK (atomInfo);
         return atomInfo->getRange (iSpecies);
      }

      bool isCreationMode () const {
         return (creationStatus == ListConstruction);
      }

      /** @{
          Get reference to atom coordinate
       */
      inline       SxVector3<TPrecTauR>& ref (ssize_t iTlAtom);
      inline       SxVector3<TPrecTauR>& ref (ssize_t iSpecies, ssize_t iAtom);
      inline const SxVector3<TPrecTauR>& constRef (ssize_t iTlAtom) const;
      inline const SxVector3<TPrecTauR>& constRef (ssize_t iSpecies,
                                                   ssize_t iAtom) const;
      inline SxVector3<TPrecTauR>& ref (SxAutoLoop &iTlAtom)
      {
         iTlAtom.setLimit (nTlAtoms);
         return ref(iTlAtom.i);
      }
      inline SxVector3<TPrecTauR>& ref (SxAutoLoop &iSpecies,
                                        SxAutoLoop &iAtom)
      {
         SX_CHECK (atomInfo);
         iSpecies.setLimit (getNSpecies ());
         iAtom.setLimit (getNAtoms ((int)iSpecies.i));
         return ref(iSpecies.i, iAtom.i);
      }
      /** @} */
      /** Get the sum of the coordinate vectors over all atoms

        \f[ \left(\begin{array}{ccc}x&y&z\end{array}\right)
            = \sum_i \left(\begin{array}{ccc}x_i&y_i&z_i\end{array}\right)
        \f]
        \note The center of mass is obtained by
        \code
Coord centerOfMass = structure.sum () / structure.getNAtoms ();
        \endcode
        */
      SxVector3<Double> sum () const;

      /// Get \f$|\mathbf r|^2\f$ of coordinates
      SxVector<TPrecTauR> absSqr () const;


      //@}

      /// @name Map atoms into cell
      //@{
      /** @brief Map atoms back into cell
        \note This modifies the structure. Use getMapped for a mapped copy.

        @sa SxCell::map
        */
      SxAtomicStructure& operator%= (const SxCell &) { return map(cell); }
      SxAtomicStructure& map (const SxCell &,
                               enum SxCell::Mapping mode = SxCell::Positive);
      SxAtomicStructure& map (enum SxCell::Mapping mode = SxCell::Positive)
      {
         return map (cell, mode);
      };

      /// Get mapped version of current structure
      SxAtomicStructure getMapped(enum SxCell::Mapping mode = SxCell::Positive)
      {
         return SxAtomicStructure(*this,Copy).map (mode);
      }

      /// Get mapped version of current structure (mapping mode Positive)
      SxAtomicStructure operator% (const SxCell &) const;

      //@}

      /// @name Rotations
      //@{
      /// Rotate atoms around origin
      inline SxAtomicStructure & operator^= (const SxMatrix<TPrecTauR> &rot);
      SxAtomicStructure & operator^= (const SymMat& rot)
      {
         return operator^= (SxMatrix<TPrecTauR>(rot));
      }
      //@}

      /// @name Vector3 dot products
      //@{
      /// Dot product
      SxVector<TPrecTauR> operator^ (const Coord &x) const;
      /** \brief Dot product of two structures
        @return Matrix of dotproducts
        */
      SxVector<TPrecTauR> operator^(const SxAtomicStructure &x) const;
      //@}


      /** \brief Assign from atomic structure (reference)

        \Note In order to get a copy, use #copy
       */
      SxAtomicStructure &operator= (const SxAtomicStructure &in);

      /** \brief Overwrite atom coordinates in an existing structure; only reorders
       * when one structures is derived from the other
          \note This routines checks for parent/child relations.
        */
      void operator<<= (const SxAtomicStructure &newPos);

      /// @name Translations
      //@{
      SxAtomicStructure operator+ (const SxAtomicStructure &in) const;
      SxAtomicStructure operator- (const SxAtomicStructure &in) const;


      SxAtomicStructure operator+ (const Coord &trans) const;
      SxAtomicStructure operator- (const Coord &trans) const;

      SxAtomicStructure & operator+= (const SxAtomicStructure &in);
      SxAtomicStructure & operator-= (const SxAtomicStructure &in);

      SxAtomicStructure & operator+= (const Coord &trans);
      SxAtomicStructure & operator-= (const Coord &trans);
      //@}

      /// @name Scaling
      //@{
      inline SxAtomicStructure operator- () const;

      /// Scale by constant factor
      template <class T>
      SxAtomicStructure & operator*= (T scale)
      {
         SX_CHECK(!isCreationMode ());
         coords *= PrecTauR(scale);
         return *this;
      }

      /// Divide by constant factor
      template <class T>
      SxAtomicStructure & operator/= (T div)
      {
         SX_CHECK(!isCreationMode ());
         coords /= PrecTauR(div);
         return *this;
      }

      /// Scale by species-dependent factor
      SxAtomicStructure operator* (const SxVector<TPrecTauR> &v) const;

      /// Scale by species-dependent factor
      SxAtomicStructure operator/ (const SxVector<TPrecTauR> &v) const;

      /// Scale by species-dependent factor
      SxAtomicStructure & operator*= (const SxVector<TPrecTauR> &v);

      /// Scale by constant factor
      inline SxAtomicStructure operator*(double s) const;

      /// Divide by constant factor
      inline SxAtomicStructure operator/ (double s) const
      {
         return *this * (1./s);
      }

      /// Divide by constant factor
      inline SxAtomicStructure operator/ (int n) const
      {
         return *this * (1./n);
      }

      /// Scale by species-dependent factor
      SxAtomicStructure & operator/= (const SxVector<TPrecTauR> &v);
      //@}

      /// Comparison operator
      inline bool operator== (const SxAtomicStructure &x) const;
      /**  @brief Comparison; if mappingPtr is given, it contains the atom mapping

        @param x structure to be compared
        @param mappingPtr points to mapping vector, which is set by this
               routine. It may be omitted.
        @return true if mapping exists

        Two structures are equal if the atoms within one species can be
        mapped such that the absolute components of the distance vector
        between the two are less than epsEqual.

        If cell is initialized, the periodicity is taken into account and
        atoms that differ only by a lattice translation are taken as equivalent.

        The mapping vector *mappingPtr will contain the mapping:
        atom i in (*this) corresponds to atom j in x --> mapping(i)=j.
        It is meaningful only if the return value is true. */
      bool isEqual (const SxAtomicStructure &x,
                    SxVector<Int> *mappingPtr = NULL,
                    const bool debug = false) const;

      /** \brief Grid-supported matching
        @param grid        subcell partitioning for this structure
        @param toBeMatched structure to be matched.

        Match a different structure against this, i.e. find out
        which atoms of this are equivalent to those in toBeMatched.

        For large structures, subcell partitioning (see #SxGrid) provides
        a better choice for structure matching.

        The return value is a reordering info that tells how the atoms of
        toBeMatched match against this. It is derived from this atomInfo.

        @return If the structures do not match, a zero info pointer
                is returned.
        */
      SxAtomInfo::ConstPtr match (const SxGrid            &grid,
                                  const SxAtomicStructure &toBeMatched) const;

      /** \brief Grid-supported atom search
        @param x     coordinate to be found
        @param grid  subcell partitioning for this structure
        @return      total atom index of first matching atom, or -1 if no match
        */
      int find (const Coord &x, const SxGrid &grid) const;

      /// Comparison operator
      inline bool operator!= (const SxAtomicStructure &x) const;

      /** \brief Check the periodicity of the structure */
      inline bool isPeriodic () const { return (cell.volume > 0.L); }

      /// Precision of comparisons
      PrecTauR epsEqual;
#     define SX_EPS_STRUCT_DEFAULT 1e-4

      /// \name Symmetries
      //@{

      /** \brief Get symmorphic symmetries (rotations around origin) [cartesian]

        This routine finds the symmetries that are compatible with the
        periodicity and the atomic structure. This does not work if
        the cell is not set.

        \sa updateSymmetries
        */
      SxList<SymMat> getSymmetries () const;

      /// Set up symmorphic symmetries internally
      inline void updateSymmetries (bool symmorphicOnly = true)
      {
         SX_CHECK(symmorphicOnly); // non-symmorphic not implemented yet
         cell.symGroupPtr = SxPtr<SxSymGroup>::create(getSymmetries ());
      }

      /// Possible modes for the symmetry group search
      enum SymSearchMode {
         /// Find also symmetries that are not compatible with the superlattice
         FullSymmetries,
         /// Find only symmetries compatible with the superlattice
         SupercellCompatible
      };
      /** \brief Determine full symmetry group (including non-symmorphic)
          @param mode if set to SupercellCompatible, only symmetries
                 where the rotational part is compatible with the superlattice
                 are determined.
          @param primCell optional primitive cell for speed-up
        */
      SxSymGroup getSymGroup (const SymSearchMode mode = FullSymmetries,
                              const SxCell &primCell = SxCell ()) const;

      /// Set up full symmetry group internally (inCell)
      inline void updateSymGroup (const SymSearchMode mode = FullSymmetries)
      {
         cell.symGroupPtr = SxPtr<SxSymGroup>::create(getSymGroup (mode));
      }

      /** \brief Get primitive cell if atomic structure is a supercell
          This routine returns the primitive cell of a supercell structure.
          The cell may have a strange shape, but the lattice is usually right.
        */
      SxCell getPrimitiveCell () const;

      /** \brief Get primitive structure if atomic structure is a supercell
          This routine returns the primitive structure of a supercell structure.

          \note primCell must be a primitive cell to the structure
        */
      SxAtomicStructure getPrimStr (const SxCell &primCell) const;

      /** \brief Return a new structure where each class of equivalent atoms
       *         forms a separate species.

          Equivalent means equal after mapping into redCell.
          The order of the species is untouched.

          @param redCell reduced cell that determines equivalence of the atoms
          @return structure with species consisting of equivalent atoms

        */
      SxAtomicStructure splitSpecies (const SxCell &redCell) const;

      /** \brief Get inequivalent atoms

        Due to symmetries, different atoms in the supercell can
        be equivalent. This routine filters out one representative
        atom for each equivalence class.

        @param primCell primitive cell with allowed symmetries
        @return a set of inequivalent atoms.
        */
      SxAtomicStructure getInequivalentAtoms (const SxCell &primCell, SxArray<int> *equivalentListPtr = NULL) const;
      SxAtomicStructure getInequivalentAtoms (SxArray<int> *equivalentListPtr = NULL) const;

      //@}

      /// \name Output
      //@{
      /** \brief formatted printout of atomic structure

          This function prints a S/PHI/nX block containing the atomic structure
          and species information. It is useful to call this function right
          after initialization of a structure. This allows the user to verify
          the provided structural data. */
      void print (const SxSpeciesData &) const;

      /// Options for fprint
      enum PrintOptions {
         /// Print no species information (except comments)
         DefaultPrint = 0x0,
         /// Print list of symmetries
         PrintSymmetries = 0x02,
         /// per default, print structure {} as outermost group. CustomGroupName suppresses this.
         CustomGroupName = 0x04
      };
      /** \brief Print to file in S/PHI/nX input format
          \param output            output file handle
          \param printOptions      see PrintOptions
          \param forces            forces

          \note there are driver routines to this function with fewer
                parameters
        */
      void fprint (FILE *output,
                   int printOptions = DefaultPrint,
                   const SxAtomicStructure &forces = SxAtomicStructure ()
                  ) const;
      /** \brief Print to file in S/PHI/nX input format
          \param output            output file handle
          \param forces            forces
       */
      void fprint (FILE *output,
                   const SxAtomicStructure &forces) const
      {
         fprint(output, DefaultPrint, forces);
      }

      /** \brief Print to file in S/PHI/nX input format in relative coordinates
          \param output            output file handle
        */
      void fprintRel (FILE *output) const;

      /** \brief Write to binary file

        \param io The binary file handle
        \param chemName Optional list of chemical symbols for species
        */
      void write (SxBinIO &io) const;
      /** \brief Write a binary file

        \param fileName output file name
        \param chemName Optional list of chemical symbols for species
        */
      void write (const SxString &fileName) const;
      /** \brief Read structure from binary file */
      void read (const SxBinIO &);
      //@}

      SxAtomicStructure cut (const Coord lower, const Coord upper) const;
};

/// Rotate atoms around origin
inline SxAtomicStructure operator^ (const SxMatrix<TPrecTauR> &rot,
                                    const SxAtomicStructure &str);

/// Rotate atoms around origin
inline SxAtomicStructure operator^ (const SymMat &S,
                                    const SxAtomicStructure &str);

/// Symmetry-transform atoms
inline SxAtomicStructure operator^ (const SxSymOp &op,
                                    const SxAtomicStructure &str);

/// Multiply with constant factor
inline SxAtomicStructure operator* (double s,
                                    const SxAtomicStructure &structure)
{
   return structure * s;
}

inline
SxAtomicStructure operator* (const SxVector<TPrecTauR> &v,
                             const SxAtomicStructure &structure)
{
   return structure * v;
}

inline SxAtomicStructure operator+ (const Coord &trans,
                                    const SxAtomicStructure &structure)
{
   return structure + trans;
}

SX_EXPORT_GEOM
SxAtomicStructure operator- (const Coord &trans,
                             const SxAtomicStructure &structure);

inline
SxVector<TPrecTauR> operator^(const Coord &x, const SxAtomicStructure &s)
{
   return s ^ x;
}

SX_EXPORT_GEOM
std::ostream &operator<< (std::ostream &s, const SxAtomicStructure &str);

#include <SxAtomicStructure.hpp>

#ifdef USE_LOOPMPI
template<>
inline void SxLoopMPI::sum (SxAtomicStructure &in)
{
   SxLoopMPI::sum (in.coords.elements, in.coords.elements, in.getNAtoms () * 3);
}
template<>
inline void SxLoopMPI::bcast (SxAtomicStructure &inout, int source)
{
   SxLoopMPI::bcast (inout.coords.elements, inout.getNAtoms () * 3, source);
}
#else
template<> inline void SxLoopMPI::sum (SxAtomicStructure &) {}
template<> inline void SxLoopMPI::bcast (SxAtomicStructure &, int) {}
#endif

#endif /* _SX_ATOMIC_STRUCTURE_H_ */
