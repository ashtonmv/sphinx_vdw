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
#ifndef _SX_DIPOLE_CORR_H_
#define _SX_DIPOLE_CORR_H_

#include <SxRho.h>
#include <SxAtomicStructure.h>
#include <SxDFT.h>

/** \brief Dipole correction for slab calculations

     A dipole in a slab calculation gives rise to an artificial
     electric field. This field is corrected for by the dipole
     correction.
    
     The correction field strength E is given by
     
     \f[ E = \frac{4 \pi \mu}{\Omega} \f]
     
     where \f$\mu\f$ is the dipole per unit cell and \f$\Omega\f$ is
     the unit cell volume.

    \author C. Freysoldt, freysoldt@mpie.de */

class SX_EXPORT_DFT SxDipoleCorrZ
{
   public:

      /** \brief Constructor
          For setting up the dipole correction, we need two
          densities: the pure electron density and the
          total charge density from which the Hartree potential
          is computed.

          The minimum of the xy-averaged electron density along z
          defines the boundary between neighbouring slabs. If the 
          density value is too large, the dipole correction is reduced
          to avoid convergence instabilities.

          The Hartree charge density is then used to compute
          the dipole.
          
          The slab must be orthogonal to the z-direction,
          and the cell must be set up correspondingly
          (see #checkCell).
          
        */
      SxDipoleCorrZ (const SxMeshR &rhoRTotal, const SxMeshR &rhoHartreeR,
                     double chargeIn = 0.);

      /** \brief Update the dipole
        @param rhoRTotal total electron density, used to find density minimum
        @param rhoHartreeR total charge density, used to compute dipole
        */
      void update (const SxMeshR &rhoRTotal, const SxMeshR &rhoHartreeR);

      /** \brief Check that the slab unit cell is ok.
          The current implementation of the dipole correction
          requires that the third axis (for which the dipole is corrected)
          is parallel to the z-axis, and the other two are orthogonal
          to it.
          
          @return true if cell is OK.
          
          */
      static bool checkCell (const SxCell &cell);

      /// \brief z-Index of the density minimum
      int rhoMinZ;

      /// Potential difference at z=0 for charged slabs
      /// \brief Minimum density value
      double rhoMin;

      bool fixedDipoleZ;
      /// Position of the dipole layer (between 0 and c)
      double dipolePos;
      /** \brief Absolute position of the dipole layer (close to previous one,
                 can be <0)
      */
      double cutZ;

      int getDipolePos (const SxMeshR &avgRho, double dZ, double accRhoMax);

      /** \brief The dipole moment of the cell

          If \f$z_0\f$ is the boundary between two slabs,
          the dipole moment is given by
          
          \f[
          \mu = \int dx dy \int_{z_0}^{z_0 + a_z} dz z \rho(x,y,z)
          \f]

          where \f$a_z\f$ is the height of the cell in the z-direction
      */
      double dipole;

      /// \brief Additional external field
      double extraField;

      /** \brief The monopole moment of the cell
        */
      double charge;
      
      /** \brief Correct the Hartree potential by a saw-tooth potential
        @param vHartreeRPtr the potential to be corrected

      */
      void correctPotential (SxMeshR *vHartreeRPtr) const;

      /** \brief Correct the forces by the dipole correction field
        @param forcePtr the forces to be corrected
        @param structure original structure (needed for charge corrections)
        @param specData species data containing the valence charges
        */
      void correctForces(SxAtomicStructure *forcePtr,
                         const SxAtomicStructure &str,
                         const SxSpeciesData &specData) const;

      SxDipoleCorrZ ()
         : rhoMinZ(0), rhoMin(0.), fixedDipoleZ (false), cutZ (0.), dipole(0.),
           extraField (0.), charge (0.)
      {
         // empty
      }

      /** \brief Get center of charge for charged slabs */
      inline double getZalign () const
      {
         if (fabs(charge) == 0.) return 0.;
         return dipolePos + dipole/charge;
      }
      
      /** \brief Get right virtual electrode potential at z=0
          @param vCorrected corrected potential in real space 
         */ 
      double getVRZero (const SxMeshR &vCorrected) const;
      /** \brief Get right virtual electrode potential at z=0
          @param vCorrected corrected potential in reciprocal space 
         */ 
      double getVRZero (const PsiG &vCorrected) const;

      /// Get electrode energy for charged slabs
      double getElectrodeEnergy (const SxCell &cell) const;
};

#endif /* _SX_DIPOLE_CORR_H_ */
