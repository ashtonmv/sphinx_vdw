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

#ifndef _SX_HUBBARD_MO_H_
#define _SX_HUBBARD_MO_H_

#include <SxDFT.h>
#include <SxNaturalCubicSpline.h>
#include <SxPAWPot.h>
#include <SxPartialWaveBasis.h>
#include <SxBlockDensityMatrix.h>

class SxHubbardU;

/** \brief MO projectors for DFT + U

    \author Christoph Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DFT SxHubbardMO {
   protected:
      /// dummy atomic structure for setupG (contains setup cell)
      SxAtomicStructure setupStructure;
      /// G-basis for MO setup
      SxGBasis setupG;

      /// Radial G basis for atomic-orbital setup
      SxPtr<SxRadGBasis> radGPtr;

      /// PAW partial wave basis (with projection to setupG)
      SxPtr<SxPartialWaveBasis> pBasis;

      /// Full AO for MO normalization
      SxDiracVec<Double> shapeAO;

      /// Lattice shifts for molecules across cell boundary
      SxArray<Coord> latShift;

   public:
      /// Force contribution from AO->MO transformation
      SxAtomicStructure trafoForce;

   protected:
      /// AO projector for site projection
      SxDiracVec<Double> shapeProj;

      /// Atom indices, molecule by molecule
      SxArray<int> atomIdx;

      /// Number of sites
      int nSite;

   public:
      int getNSite () const { return nSite; }

      /// Offset in the complete list of Hubbard U sites
      int siteOffset;

   public:
      /// Constructor
      SxHubbardMO (int siteOffsetIn = 0);

      /// Real-space truncate a shape provided in radial G space
      SxDiracVec<Double> truncateShapeG (const SxDiracVec<Double> &shape,
                                         double rCut,
                                         double width) const;
   protected:
      /** \brief Setup up atomic orbital and Hubbard projector
          @param ao  Atomic orbital in radial real space
          @param rCut cutoff radius for the truncated projector
          @param truncWidth transition for cutoff
          @param pawPotPtr PAW potential
          @param eCut plane-wave energy cutoff

          The species identity (is) and l quantum number is taken from ao.
          The m value (>=0) of ao determines the rotational quantum number
          of the MO constructed from the ao.
        */
      void setupAO (const SxDiracVec<Double> &ao,
                    double rCut, double truncWidth,
                    const SxPtr<SxPAWPot> &pawPotPtr,
                    double eCut);
      /** \brief Set up the (simple cubic) setup box and its G basis

          @param rSize cubic lattice constant
          @param eCut  plane-wave cutoff
        */
      void setupBox (double rSize, double eCut);

      /** \brief Find the molecules according to MO Group in input file
          @param moGroup HubbardU.MO group
          @param structure structure
          @param speciesInfo needed for chemical symbols, use PAW potential
          @return species id
          */
      int findMolecules (const SxSymbolTable* moGroup,
                         const SxAtomicStructure &structure,
                         const SxSpeciesData &speciesInfo);
      /** \brief Find the atoms according to AO Group in input file
          @param aoGroup HubbardU.AO group
          @param structure structure
          @param speciesInfo needed for chemical symbols, use PAW potential
          @return species id
          */
      int findAtoms(const SxSymbolTable* aoGroup,
                    const SxAtomicStructure &structure,
                    const SxSpeciesData &speciesInfo);
   public:
      /// Read and setup for homonuclear diatomics
      void read (const SxSymbolTable *table,
                 const SxAtomicStructure &structure,
                 const SxPtr<SxPAWPot> &potPtr);
      /// Read and setup for single atoms
      void readAO (const SxSymbolTable *table,
                   const SxAtomicStructure &structure,
                   const SxPtr<SxPAWPot> &potPtr);
      /// Remove intermediate data used for setup
      void finalize ();
   protected:
      /// Minimum distance
      double distFrom;
      /// dr for pNorm interpolation
      double dDist;
      /// Constituent AO's l value
      int l;
      /// MO rotational momentum
      int mMO;
      /// Phase between the atoms
      double sign;

      /// Number of atomic orbitals per site
      int nAoPerSite;

   public:
      /// Return number of atomic orbitals per site (=2 atoms)
      int getNAOperSite () const {
         return nAoPerSite;
      }
      /// Number of interpolation points
      int nInterpolate;
      /// Number of radial points
      int nRad;
      /// Produce debug output
      bool verbose;
      /// Prefix for verbose output files
      SxString prefix;
   protected:
      /// MO projector normalization as function of distance
      SxNaturalCubicSpline pNorm;

      /// Get projector norm
      double getPNorm (double dist)  {
         //return 0.7;
         return pNorm.getValYExtra ((dist - distFrom) / dDist);
      }

      /// Get derivative of projector norm
      double getPNormDeriv (double dist)  {
         //return 0.;
         return pNorm.getDerivYExtra ((dist - distFrom) / dDist) / dDist;
      }

      /** \brief Orbital rotation from the standard z-axis orientation
          to a given molecular axis.

          @param axis the molecular axis (interatomic vector)
          @return a 2x(2L+1) matrix, containing the rotation of the
                  atomic orbital projectors (2L+1 many) to the -M and
                  +M site. The value of M (rotational quantum number)
                  is taken from shapeProj, which gets its value from
                  the #setup routine.
          */
      SxVector<Double> getRot (const Coord &axis) const;

      /// Validate distance
      bool validateDist (double dist) const;

      /** Set up normalization table
          The Hubbard MO projector must be normalized such that
          <MO | P_mm | MO > = 1
          if MO is the normalized orbital for the site index m.


          Since P_mm' = |p_m> N <p_m'|, where |p_m> is a
          modified orbital (shapeProj), the normalization constant
          is
          N = <MO|\hat S|MO>/<MO|p_m> <p_m|MO>

          Note that our shapeProj is set up from two truncated atomic
          orbitals that include the PAW overlap operator for the
          atom they are associated with. In contrast, the MO is set up
          from the full atomic orbital (shapeAO) and uses the
          overlap operator for both atoms.

          For homonuclear diatomic molecules, N depends only on the
          interatomic distance, since the shape perpendicular to the axis
          is dictated by symmetry (the rotational quantum number).

          The number of points is nInterpolate, which defaults to 100.

        */
      void setupNormalization (double distFromIn,
                               double distTo);

   public:
      /// Atomic orbital projectors
      SxPtr<SxAOBasis> aoProj;

      /// Setup AO projectors in G+k space
      void setupProjGk (const SxGkBasis &gk);

      /** \brief Hubbard Hamiltonian in the AO projector space :iSpin, :iSite
        */
      SxBlockDensityMatrix hamProj;

      /** \brief Compute the Hubbard energy and Hamiltonian
        @param hubbardU the Hubbard parent object containing the U parameters,
                        energy, etc.
        @param Pij the AO block density matrix
        @param structure the atomic structure
        */
      void compute (SxHubbardU *hubbardU,
                    const SxBlockDensityMatrix& Pij,
                    const SxAtomicStructure &structure);
   protected:

      /** \brief Compute the Hubbard energy and Hamiltonian for homonuclear
                 diatomics
        @param hubbardU the Hubbard parent object containing the U parameters,
                        energy, etc.
        @param Pij the AO block density matrix
        @param structure the atomic structure
        */
      void computeMO (SxHubbardU *hubbardU,
                      const SxBlockDensityMatrix& Pij,
                      const SxAtomicStructure &structure);
      /** \brief Compute the Hubbard energy and Hamiltonian for single atoms
        @param hubbardU the Hubbard parent object containing the U parameters,
                        energy, etc.
        @param Pij the AO block density matrix
        @param structure the atomic structure
        */
      void computeAO (SxHubbardU *hubbardU,
                      const SxBlockDensityMatrix& Pij,
                      const SxAtomicStructure &structure);
   public:
      /** \brief Add contribution of one k-point to the AO block density matrix
        @param Pij The block density matrix to be computed
        @param waves wave functions for current k-point
        @param weight k-point weight
        @param focc occupation number for each state
        */
      void addToRho (SxBlockDensityMatrix *Pij,
                     const PsiG &waves,
                     double weight,
                     const SxDiracVec<Double> &focc) const;

      /// Symmetrize the AO block density matrix
      void symmetrize (SxBlockDensityMatrix *Pij,
                       const SxAtomicStructure &structure,
                       const SxYlmRotGroup &ylmRot) const;
      /** \brief Compute contribution to gradient of the AO block density matrix

          @param fermi Fermi occupations
          @param P     <AO|Psi> projections
          @param gradP <d/dtau AO | Psi> projection gradients
          @return gradient matrices for all sites. The first ("spin") index is
                  used for the direction.

          @note k-point and spin must be set in the projection's (P) auxData.
        */
      SxBlockDensityMatrix
      computeGradPij (const SxFermi &fermi,
                      const SxDiracVec<Complex16> &P,
                      const SxArray<SxDiracMat<Complex16> > &gradP) const;
      /** \brief symmetrize the gradient of the AO block density matrix

          @param Pij       pointer to the AO block density matrix gradient
          @param structure atomic structure
          @param ylmRot    symmetries in spherical harmonics basis
        */
      void symmetrizeGradPij (SxBlockDensityMatrix *Pij,
                              const SxAtomicStructure &structure,
                              const SxYlmRotGroup &ylmRot) const;

      /** \brief Get forces from the displacement of the AO projectors

          @param gradPij gradient of block density matrix with respect to
                         displacement of the left-hand atom
          @param iSpin   spin-channel this block density matrix refers to
          @return The corresponding forces.
      */
      SxAtomicStructure
      getForce (const SxBlockDensityMatrix &gradPij, ssize_t iSpin) const;
};

#endif /* _SX_HUBBARD_MO_H_ */
