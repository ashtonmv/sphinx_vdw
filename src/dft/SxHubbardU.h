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

#ifndef _SX_HUBBARD_U_H_
#define _SX_HUBBARD_U_H_

#include <SxDFT.h>
#include <SxDirac.h>
#include <SxBlockDensityMatrix.h>
//#include <SxHubbardMO.h>
class SxPAWPot;
class SxAtomicStructure;
class SxSymbolTable;
class SxHubbardMO;


/** \brief Hubbard U extension in the "rotationally invariant" formulation

 * Reference 1: SPHInX Hubbard U implementation notes
 * Reference 2: M. Cococcioni, S. de Gironcoli, Phys. Rev. B 71, 035105 (2005)

    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DFT SxHubbardU
{
   public:
      /// Constructor
      SxHubbardU ();

      /// The effective U at each site
      SxArray<double> Ueff;

      /// energy
      double energy;

      /// double counting correction energy
      double eDoubleCounting;

      /// Set energy etc. to zero (before computeIncr)
      void setZero ();
      /** \brief Incremental compute
          @param iSite  which site
          @pram  nij    occupation matrix

          @return \f$\partial E/\partial n_ij
        */
      SxVector<Complex16> computeIncr(int iSite,
                                      const SxVector<Complex16> &nij);

      /** \brief Site bias for fitting/constraints
        */
      SxArray<double> alphaBias;

      /// The total occupation each site
      SxArray<double> siteOccupation;

      /// Get number of sites
      inline ssize_t getSize ()  const
      {
         SX_CHECK (Ueff.getSize () == siteOccupation.getSize (),
                   Ueff.getSize (), siteOccupation.getSize ());
         return Ueff.getSize ();
      }
      /// Resize (keeping current values)
      void resize (ssize_t nSite);

      /** \brief atomic site map

          This container class maps atoms to Hubbard U sites.
          This map must be set up from the outside.

        */
      class AtomicSite {
         public:
            /// Site ID in SxHubbardU
            int iSite;
            /// l quantum number for the orbital type
            int l;
            /** \brief Where m=-l is located within the intra-atomic
                       all-projector index (ipl)
                @note same as SxPAWPot.offset
            */
            int offset;
            /// Inverse Norm of the projector ref 1, Eq. (7)
            double Pnl;
      };

      /// Atomic site map, listing for each atom the Hubbard U sites (:iTlAtom,iSite (at atom)
      SxArray<SxArray<AtomicSite> > atomicSiteMap;

      /** \brief List of MO sites (iTypeMO)

          This list contains the types of MO sites. Each type can consist
          of several sites
        */
      SxArray<SxPtr<SxHubbardMO> > moSite;

      /** \brief Compute contribution of atom iTl
          @param iTl    which atom
          @param Dij    PAW density matrix
          @param fFull  value of full occupation: 1 for spin-polarized
                        calculations, 2 otherwise
          @return contribution to Vij
        */
      SxVector<Double> computeAtom(int iTl, const SxVector<Double> &Dij,
                                     double fFull);

      /// Compute MO sites
      void computeMO (const SxArray<SxPtr<SxBlockDensityMatrix> > &Pij,
                      const SxAtomicStructure &structure);

      /// Read atomic sites from symbol table
      void read (const SxSymbolTable *table,
                 const SxPtr<SxPAWPot> &pawPotPtr,
                 const SxAtomicStructure &structure);
};

#endif /* _SX_HUBBARD_U_H_ */
