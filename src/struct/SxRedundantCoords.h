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

#ifndef _SX_REDUNDANT_COORDS_H_
#define _SX_REDUNDANT_COORDS_H_

#include <SxStruct.h>
#include <SxAtomicStructure.h>
#include <SxSortedList.h>

namespace {
   class SxNearN;
}

/** \brief Redundant internal coordinates

    This class provides redundant internal coordinates for structure
    optimizations. At present, only interatomic distances are implemented.

    The coordinates are automatically set up from an initial atomic structure.
    For this, we first sort interatomic distances into sets of similar cases
    (based on element combination and distance). Then, we select for each
    atom one or more sets such that they form a nearest-neighbor environment.
    The #maxDist, #primarySetLimit, and #rmsThreshold parameters control
    the set generation. The #planeCutLimit control the final selection.

    Each set of interatomic distances is associated with an interatomic
    force constant (saved in #param). This allows to link atomic displacements
    with (very approximate) forces. To get forces from displacements, use
    #applyH. To get the best set of displacements for a given set of forces,
    use #solve. To parameterize, you need the derivative of the forces with
    respect to the parameters for a given displacement from #getParamDeriv.
    Since the interatomic distance vectors vary with the atomic structure,
    these routines always require the current atomic structure.

    \author Christoph Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_STRUCT SxRedundantCoords
{
   public:
      /// Max. distance for initial neighors
      double maxDist;
      /** \brief Parameter to define primary interatomic distance sets
          In a first step, all similar interatomic distances are sorted into
          sets. This parameter defines "similar" as the maximum of
          |log(d1/d2)|. The default is 0.05 (~ 5%)
        */
      double primarySetLimit;
      /** \brief Parameter to control primary set merging
          After the interatomic distance sets have been determined, we
          can determine their root-mean-square variation. If the centers
          of two sets lie close compared to their joint rms, we merge
          them. This parameter determines the limiting center distance
          as a multiple of the joint rms. The default is 3.
        */
      double rmsThreshold;
      /** \brief Parameter to control final selection
          Once the possible candidates for nearest neighbors have been
          determined, we set up the final list "shell by shell". For
          any intermediate list of neighbors, we construct a polyhedron
          by putting planes perpendicular to the interatomic connection
          line, at a particular position relative to the interatomic
          distance.
          New shells are only included if one of their atoms lies within
          this polyhedron. The parameter defines the relative position
          of the polyhedrons boundary planes (1 = through the neighbor,
          0=through the origin). The default is 0.95.
        */
      double planeCutLimit;

      /// Control verbosity
      mutable bool verbose;

      /// Constructor
      SxRedundantCoords ();

      /// Automatic setup
      void setup (const SxAtomicStructure &structure);

      /** \brief Verify that the bond classification is still ok
          @param str current structure
          @return whether or not the classification is OK for str

          This routine reruns the classification algorithm on each
          current type of bond. If any of these types would now be
          classified into two or more independent types, the function
          returns false, suggesting that the setup routine should
          be called again.
        */
      bool verifyClasses (const SxAtomicStructure &str) const;

      /// Classify bond angles
      void getAngles (const SxAtomicStructure &structure);
      
      /// Set Born-von-Karman bond-orthogonal terms for selected atoms
      void getBornVonKarmanAngles (const SxArray<int> &atoms);
   protected:
      /// Auxiliary routine: get neighbor candidates
      SxArray<SxSortedList<SxNearN> >
      getCandidates(const SxAtomicStructure &structure);
      /** \brief Auxiliary routine: classify bonds according to bond length
          @param pairDist sorted(!) list of log (squared distance)
          @param setStartPtr will contain entry points for each class
          @return cycle set: cyclic set of indices of bonds belonging to
                  the same class.
        */
      SxArray<ssize_t> classifyBonds (const SxArray<double> &pairDist,
                                      SxList<ssize_t> *setStartPtr) const;
      /// Total atom id of neighbors (iTl:iNeighbor)
      SxArray<SxArray<ssize_t> > idAtom;
      /// Type of neighbor (iTl:iNeighbor)
      SxArray<SxArray<int> > ricType;
      /** \brief Lattice shift for each neighbor (iTl:iNeighbor)

        In order to determine the relative position of each neighbor,
        we have to add the lattice shift to the reference atom (as specified
        by idAtom), and subtract the absolute position of the current atom.
        This detour is necessary because atomic positions may change, and
        jTl and latShift is sufficient to reconstruct the absolute neighbor
        position.
        */
      SxArray<SxArray<Coord> > latShift;

      /// Number of bond types
      int nBond;
      /// Number of angle types
      int nAngle;

      /// Bond angle type for each pair of neighbors
      SxArray<SxArray<int> > angType;
   public:
      /// Get neighbor ids
      const SxArray<ssize_t>& getNeighbors(int iTl) const {
         return idAtom(iTl);
      }

      /// Get number of parameters for Hessian
      int getNParam () const { return nBond + nAngle; }

      /// Get number of parameters for Hessian
      int getNBond () const { return nBond; }

      /// Hesse Parameters
      SxVector<Double> param;

      /// Classification of bond
      class BondType {
         public:
            int iSpecies1, iSpecies2;
            double dAvg;
      };

      /// Bond type classifications
      SxArray<BondType> bondTypes;

      /// Classification of bond angle
      class AngleType {
         public:
            int iSpeciesCentral, jSpecies1, jSpecies2;
            int iBondType1, iBondType2;
            double cosAvg;
      };

      /// Bond angle type classifications
      SxArray<AngleType> angleTypes;

      /** \brief Get force change
        @param str current atomic structure
        @param x (small) displacement
        @return (approximate) forces
        */
      SxVector<Double> applyH (const SxAtomicStructure &str,
                               const SxVector<Double> &x) const;

      /** \brief Get parameter derivative for bond force constants
          This returns a nDof x nBond matrix
        */
      SxVector<Double> getParamDeriv (const SxAtomicStructure &str,
                                      const SxVector<Double> &x) const;
      /** \brief Get parameter derivative for bond angle constants
        @note This should be used on forces without the contribution of along-bond forces.
          This returns a nDof x nBond matrix
        */
      SxVector<Double> getParamDerivA (const SxAtomicStructure &str,
                                       const SxVector<Double> &x) const;

      /** \brief Apply inverse Hessian
        @param str current atomic structure
        @param dF forces
        @param normWeight regularization parameter for underdetermined cases
        @param accuracy accuracy for iterative solver
        @return displacements

        This routine use a conjugate-gradient iterative solver
        to find the best displacements to reproduce the forces.
        For underdetermined systems, the regularization parameter
        is needed. It adds a normalization term (Tikhonov regularization),
        i.e., the quantity being minimized is
        \f[
        |\Delta \mathbf F - \mathbf H \mathbf x|^2 + \kappa^2 |\mathbf x|^2
        \f]
        For hard modes, the first term dominates. For soft modes, the modes
        are smoothly faded out by the second term.
        This may not work perfectly. It is strongly recommended
        to use an auxiliary (linear) approximation for the
        remainder:
        \code
   double diagInvRest;
   SxRedundantCoords ric;
   SxVector<Double> dF, dxRic, dFrest, dx;
   ...
   dxRic = ric.solve (str, dF, accuracy);
   dFrest = dF - ric.applyH (str, dxRic);
   dx = dxRic + diagInvRest * dFrest;
        \endcode
        */
      SxVector<Double> solve (const SxAtomicStructure &str,
                              const SxVector<Double> dF,
                              double normWeight,
                              double accuracy) const;

      /** \brief Get param values from Schlegel model
        @param elem chemical element symbols

         See B. Schlegel, Theor. Chim. Act. 66, 333 (1984). Table I
         @note: rows > 3 are treated like row=3 :-((
      */
      void paramSchlegel (const SxArray<SxString> &elem);

      /** \brief Get param values from Fischer model
        @param elem chemical element symbols

         See T.H. Fischer, J. Almloef, J. Phys. Chem. 96,9768-9774 (1992),
             Table I
      */
      void paramFischer (const SxArray<SxString> &elem);

};

#endif /* _SX_REDUNDANT_COORDS_H_ */
