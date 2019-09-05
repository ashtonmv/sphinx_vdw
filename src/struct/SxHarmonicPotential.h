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

#ifndef _SX_HARMONIC_POTENTIAL_H_
#define _SX_HARMONIC_POTENTIAL_H_

#include <SxStruct.h>
#include <SxAtomicStructure.h>
#include <SxPotential.h>
#include <SxSpeciesData.h>
#include <SxGBasis.h>

/** \brief A simple harmonic potential on the displacements


    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_STRUCT SxHarmonicPotential 
   : public SxPotential,
     public SxSpeciesData
{
   public:
      /// List of neighbor ids for each atom (iTlAtom,iNeighbor)
      SxArray<SxArray<int> > neighborIdx;

      /// List of Hesse matrices for each neighbor pair (iTlAtom,iNeighbor)
      SxArray<SxArray<SxMatrix3<Double> > > hesse;

      /** \brief List of equilibrium displacement for each neighbor pair
                 (iTlAtom,iNeighbor)
        */
      SxArray<SxArray<Coord > > dist0;

      /// Dielectric constant (for electrostatic terms)
      double dielecConstant;

      /// @name Anisotropic electrostatics
      ///@{
      /// Born effective charge tensors
      SxArray<SxMatrix3<Double> > Qab;
      /*
      /// Inverse dielectric matrix
      SxMatrix3<Double> epsInv;
      */

      /// Indices of atoms with anisotropic charge tensor
      SxArray<int> idAniso;

      /// Map total atom number to anisotropic id
      SxArray<int> mapAniso;

      /// For each anisotropic atom: reference "frame"
      SxArray<SxArray<int> > frameAtoms;

      /// For each anisotropic atom: reference "frame"
      SxArray<SxArray<Coord> > frameAtomPos;

      /// For each atom: isotropic charge
      SxVector<Double> qIso;

      /// Traceless anisotropy for displacements with respect to reference frame
      SxArray<SxMatrix3<Double> > QFrame;

      /// Position of an anisotropic atom with respect to its reference frame
      SxArray<Coord> framePositions;

      /** \brief Split isotropic and anisotropy of charge tensors
          @param refStruct    reference structure
          @param QabBulk      (for defects) Qab of the bulk system
          @param genId        (for defects) map of defect atoms

          @note
          For defects, corrections are made to the harmonic neighbor
          interactions for those atoms that belong to the bulk region,
          but have a modified anisotropy because they belong to the
          anisotropy reference frame of some defect atom.
        */
      void setupDeltaQ (const SxAtomicStructure &refStruct,
                        const SxArray<SxMatrix3<Double> > &QabBulk,
                        SxArray<int> genId);

      /** \brief Set force ref
        @f forces
        @x if given, the forces refer to these atoms (may be a subset)
        */
      void setExtraForce (const SxAtomicStructure &f,
                          const SxAtomicStructure &x = SxAtomicStructure ());

      /// Whether Coulomb interactions need to be calculated
      bool noCoulomb;

   protected:
      /// Reference structure
      SxAtomicStructure refStr;
      /// External forces on the reference structure
      SxAtomicStructure extraForce;

      /// Drift correction
      SxAtomicStructure driftCorr;

      /// Harmonic potential
      SxMatrix<Double> potEl;
      /// Setup Qab/dielecMatrix from valenceCharge / dielecConstant
      void setupQab ();

      ///@}

      /** \brief Shift interaction to neighbor

        @note This is needed for defects to host all
        defect-bulk interactions at the defect atoms.
        */
      void shiftHesseToNeighbor (int ia, int in);

   public:

      /// Number of generating neighbor sets
      SxVector<Int> nGenSets;

      /** \brief Matrix symmetrizers for generating neighbors
          These matrices map the free parameters of each generating
          Hesse matrix to the 9 entries of the matrix.
        */
      SxArray<SxMatrix<Double> > matSyms;
      /** \brief Matrix symmetrizers for charge tensors
          These matrices map the free parameters of each charge tensor 
          to its 9 entries.
        */
      SxArray<SxMatrix<Double> > matSymQ;


   protected:
      /** \brief Calculate symmetry constraints

        The symmetry constraint is
        \f[
        \mathbf S \mathbf H \mathbf S^{-1} = \mathbf H
        \f]
        This translates to
        \f[
        H(j2,k2) - \sum_{j1,k1} S(j2,j1) H(j1,k1) S(k2, k1) = 0
        \f]
        which can be written as
        \f[
        \sum_{j1,k1} (delta(j1,j2)delta(k1,k2)- S(j2,j1) S(k2, k1)) H(j1,k1)  = 0
        \f]

        The routine automatically includes the inverse symmetries (switch 1 and 2).
        
        The routine determines the free parameter space allowed for the symmetries,
        imposing optionally the constraint that H is a symmetric matrix.
        
        @param matSymId  ids of the symmetries
        @param symmetric  if true, constrain to symmetric matrices
        */
      SxMatrix<Double> getMatrixConstraint (const SxList<int> &matSymId,
                                            bool symmetric = true) const;

   public:
      /// This routines returns the free parameter space for all generating matrices
      SxMatrix<Double> getParamSpace () const;
   protected:
      ///@{
      /// \name Parametrization

      /// Generating neighbor list
      SxArray<Coord>  genNeighbors;
      /// Symmetry list
      SxArray<SymMat> syms;
      /// Map of neighbor origin: which genNeighbor (iTlAtom,iNeighbor)
      SxArray<SxArray<int> >  genNeighborId;
      /// Map of neighbor origin: which sym (iTlAtom,iNeighbor)
      SxArray<SxArray<int> >  genSymId;

   public:
      /// Map of charge tensor origin: which charge tensor (iTlAtom)
      SxArray<int> genQ;
      /// Map of charge tensory origin: which sym (iTlAtom)
      SxArray<int> genQSymId;
      /// Get a symmetry matrix
      const SxMatrix3<Double>& getSym (int iSym) const { return syms(iSym); }


   public:
      Coord getGenPos (int iGen, int iSym = -1) const
      {
         SX_CHECK (iGen >= 0 && iGen < genNeighbors.getSize (),
                   iGen, genNeighbors.getSize ());
         SX_CHECK (iSym < syms.getSize (), iSym, syms.getSize ());
         if (iSym >= 0) 
            return syms(iSym) ^ genNeighbors(iGen);
         else
            return genNeighbors(iGen);
      }
      /// Clear the stored data on neighbor origin
      inline void clearGenInfo ()
      {
         genNeighbors.resize (0);
         syms.resize (0);
         genNeighborId.resize (0);
         genSymId.resize (0);
      }

      /// Number of parameters in the generating Hesse matrices
      inline int getNParam () const
      {
         return int(genNeighbors.getSize ()) * 9;
      }

      /** \brief Gradient of forces with respect to generating
                 Hesse matrices
        return df(a)/d G(b,c)

        */
      SxMatrix<Double> paramGrad (const SxAtomicStructure &tau);

      /** \brief Get parametrization matrix for a charge tensor
        
        The force on an atom due to a local electric field is
        \f[
           f^{(ia)}_\alpha = -Q^{(ia)}_{\beta\alpha} E^{(ia)}_\beta
        \f]

        The charge tensor at atom ia originates from the g'th charge
        tensor in the input file via
        \f[
        \mathbf Q^{(i)} = \mathbf S^{g:ia} \mathbf Q^g 
        \left(\mathbf S^{g:ia} \right)^T
        \f]

        Last, symmetry constraints relate the elements of the charge tensor
        to its degrees of freedom (as determined by #getMatrixConstraint)
        via
        \f[
        Q^g_{\gamma \delta} = \sum_n U^{g}_{\gamma\delta n} p^{g}_n
        \f]

        Minimizing the error between the DFT force and the above
        expressions 
        \f[
        R^2 = \sum_{ia} \left|\mathbf f^{(ia),\textrm{DFT}} + \mathbf (\mathbf E^{(i)})^T Q^{(ia)}\right|^2
        \f]
        with respect to \f$p^g_n\f$ yields
        \f[
        0 = \frac{\partial R^2}{\partial p_n}
        = \sum_{ia,\alpha} \Gamma^{(ia)}_{n\alpha} \sum_{n'}
          \Gamma^{(ia)}_{n'\alpha} 
        + \sum_{ia,\alpha} \Gamma^{(ia)}_{n\alpha} f^{(ia),\textrm{DFT}_\alpha}
        \f]
        with
        \f[
        \Gamma^{(ia)}_{n,\alpha} = \sum_{\beta\gamma\delta}
        U^g_{\gamma\delta n} S^{g:ia}_{\beta \gamma} S^{g:ia}_{\alpha \delta}
        E^{(ia)}_\beta
        \f]
        This routine returns the \f$\Gamma$ matrix for one particular
        atom. The corresponding g is listed in genQ.

        @param iTlAtom  atom id
        @param field  the field at that atom
        @return \f$\Gamma\f$ as a N^g_p x 3 matrix

        */
      SxMatrix<Double> paramQ (int iTlAtom, const Coord &field) const;
      ///@}

      void setHesse (const SxVector<Double> &param);

   public:

      // --- Virtual functions
      virtual ~SxHarmonicPotential ();

      /// SxPotential init function
      virtual void init (const SxAtomicStructure &, const SxSymbolTable *);

      /** \brief Check whether the provided minimization command is known */
      virtual bool isRegistered (const SxSymbolTable *) const
      {
         return true;
      }

      /** \brief Execute the minimization

          \param cmd   minimization command
          \param calc  if \b false the command's input parameter are printed
                       only. */
      virtual void execute (const SxSymbolTable *, bool =true)  {
         return;
      }

      /** \brief compute (unsymmetrized) forces.

          Compute the unsymmetrized forces to each atom according to the
          provided minimization command. Note, usually the forces needs
          to be symmetrized according the the symmetry elements exsiting
          in the current system. In this case SxPotential::getSymForces
          must be used instead.

          \param tau       atomic stucture \f$ \{ \tau_{i_s i_a} \} \f$
          \param cmd       minimization command group */
      virtual SxAtomicStructure getForces (const SxAtomicStructure &tau,
                                           const SxSymbolTable *cmd=NULL);


      /** \brief Returns energy of the current Born-Oppenheimer update

          This function returns the energy computed during the current
          Born-Oppenheimer update, i.e., during the previous call of
          SxPotential::getForces or SxPotential::getSymForces 

          \return energy in Hartree */
      virtual PrecEnergy getEnergy () const;

      /// Get the species data
      virtual SxSpeciesData getSpeciesData () const { return *this; }

      /** \brief Write the current structure
        @param fileName file to write to
        */
      virtual void writeStructure (const SxString &) const {}

   public:
      /// Constructor
      SxHarmonicPotential (const SxSymbolTable *table,
                           const SxAtomicStructure &refStruct);

      /// External electrostatic potential (unscreened)
      PsiG vExt;

      /// Set extra potential
      void setExtraPot (const SxMeshR &vR, bool screened = false);

      /// Get electrostatic forces (and compute energy)
      SxAtomicStructure getElForces (const SxAtomicStructure &tau);

      /// Get 2nd derivative
      SxMatrix<Complex16> getHesse (const Coord &kVec, 
                                    const SxAtomicStructure &tau);
      /// Get 2nd derivative for electrostatic part
      SxMatrix<Complex16> getHesseEl (const Coord &kVec, 
                                    const SxAtomicStructure &tau);
   protected:
      /// Calculate the displacement-induced dipoles
      SxAtomicStructure getMu (const SxAtomicStructure &tau) const;
      /** \brief Auxiliary function for getHesseEl
          Adds 4-th term in Eq. 46 to the Hessian
        */
      void addH4 (int iTl, int jTl,
                  SxMatrix<Complex16> *resPtr,
                  const SxMatrix3<Complex16> &d2EdMu2,
                  const Coord &kVec);
      /** \brief Auxiliary function for getHesseEl
          Adds 2nd and 3rd term in Eq. 46 to the Hessian
        */
      void addH23 (int iTl, int jTl,
                   SxMatrix<Complex16> *resPtr,
                   const PsiG &d2EdRdMu,
                   const Coord &kVec);

      /** \brief Get 2nd derivative for electrostatic part (isotropic charges)
        \note Charge corrections are applied for fields along k-vector

        */
      SxMatrix<Complex16> getHesseElIso (const Coord &kVec, 
                                         const SxAtomicStructure &tau);
   protected:
      /// G-basis for electrostatics
      SxPtr<SxGBasis> gBasis;

      /// Gaussian width
      double beta;

   public:
      /// Gaussian width treated with G-vectors
      static double betaSoft;

   protected:
      /// last energy
      double energy;
};

#endif /* _SX_HARMONIC_POTENTIAL_H_ */
