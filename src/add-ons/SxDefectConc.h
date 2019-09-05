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

#ifndef _SX_DEFECT_CONC_H_
#define _SX_DEFECT_CONC_H_

#include<SxList.h>
#include<SxString.h>
#include<SxMatrix.h>
#include<SxSymbolTable.h>
#include <SxExt.h>

class SX_EXPORT_EXT SxDefectConc {
   public:
   /// Data describing a single defect
   class Defect {
      public:
         /// Standard formation energy
         SxVector<Double> E0;
         /// Defect charge
         SxVector<Double> charge;
         /// Dependence on chemical potentials
         SxVector<Double> sumFormula;
         /// User description
         SxString name;
         /// \brief Multiplicity (number of configurations)
         int multiplicity;

         /// Diffusion barrier
         double diffusionBarrier;
         /// Diffusion prefactor
         double diffusionPrefactor;
         /// Diffusion charge (for field-driven term)
         double diffusionCharge;

         inline double getDiffusionConstant (double kT) const
         {
            return diffusionPrefactor * exp (-diffusionBarrier / kT);
         }

         /// Empty constructor
         Defect ();
         /** \brief Set defect from input file
           @param table  Table for defect
           @param parent defect container class
           */
         Defect (const SxSymbolTable *table, const SxDefectConc &parent);

         /** \brief Compute defect concentration

           \f[
           c = N e^{-\frac{E}{kT}}
           \f]
           where the formation energy is given by

           \f[
           E = E_0 - \sum_i n_i \mu_i + q E^{\rm Fermi}
           \f]

           \f$n_i\f$ is the composition of the defect (number of extra
           atoms per species), \f$mu_i\f$ the chemical potentials.

           N is the number of configurations.

           \note This concentration is relative to the "site concentration".

           @param kT          \f$kT\f$
           @param mu          the chemical potentials
           @param fermiEnergy the Fermi energy
           @return the concentration
           
           */
         SxVector<Double> getConcentration (double kT,
                                  const SxVector<Double> mu,
                                  double fermiEnergy = 0.) const;
         /** \brief Get charge density for given defect concentration
              @param kT      temperature
              @param eFermi  Fermi energy
              @param cTotal  total defect concentration (all charge states)
          */
         double getChargeDensity (double kT, 
                                  double eFermi, 
                                  double cTotal) const;
         
         /** \brief Compute defect energy

           \f[
           E = E_0 - \sum_i n_i \mu_i + q E^{\rm Fermi}
           \f]
           @param mu          the chemical potentials
           @param fermiEnergy the Fermi energy
           @return the formation energy
           
           */
         inline SxVector<Double> getEnergy (const SxVector<Double> mu,
                                            double fermiEnergy = 0.) const
         {
            return E0 + charge * fermiEnergy - (sumFormula * mu).sum ();
         }

         /// Print
         void print (const SxList<SxString> &elements) const;

   };

   class Band {
      public:
         /// Energy
         double E0;
         /** Effective mass 
             In practice, this is
             \f[ \frac {m_{\rm eff} m_e}{\hbar^2} \f]
             in units of energyUnit^{-1} lengthUnit^-2 
           */
         double mass;
         /// Degeneracy of band
         int degeneracy;
         /// Name
         SxString name;

         /// Constructor
         Band (const SxSymbolTable *table, double eMassByhBarSqr);

         Band () { /* empty */}

         /// Compute the polylogarithm \f$\textrm{Li}_{3/2}(x)\f$
         static double polyLog3_2 (double x)  
         {
            // definition: Li_n(x) = \sum_{i=1}^\infty x^i / i^n
            SX_CHECK (fabs(x) < 1., x);
            if (fabs(x) < 1e-12) return x;
            double res = x, xi = x, d, ii = 2.;
            for (int i = 2; i < 1000; ++i, ii+=1.)  {
               xi *= x; // x^i
               d = xi / (ii * sqrt(ii));
               res += d;
               if (fabs(d) < 1e-12 * fabs(res)) return res;
            }
            SX_EXIT; return 0.;
         }
         /** Get charge density for given Fermi energy
             Integral for electrons (E >> E^{\rm Fermi})
             \f[
             \int_0^{\infty} dE e^{- (E-E^{\rm Fermi}) / kT} 4 \pi 
                                \sqrt{2 m_e m_{\rm eff} \hbar^{-2} E}
             = 
             \f]
           */
         double getChargeDensity (double eFermi, double kT) const;

         void print (double kT) const;
   };

   class Reaction {
      public:
         /// Frequency prefactor
         double prefactor;
         /// Charge of transition state
         double charge;
         /// Energy barrier (absolute scale)
         double barrier;
         /// Reaction coefficients (educts: negative)
         SxVector<Double> coefficients;

         /// Rate constant forward reaction
         double rateConstantForward;
         /// Rate constant backward reaction
         double rateConstantBackward;

         /// Compute rate constants
         void computeConstants (const SxDefectConc &defconc,
                                double eFermi,
                                double ekt);

         /// Compute rate
         SxVector<Double> getRate (const SxVector<Double> &conc);

         Reaction () 
            : prefactor (0.), charge (0.), barrier (0.), 
              rateConstantForward (0), rateConstantBackward(0.)
         { /* empty */ }
                                   
   };

   class ReactionSystem {
      public:
         /// Reactions
         SxArray<Reaction> reactions;

         /// Get size
         ssize_t getSize () const { return reactions.getSize (); }
         /// Concentrations
         SxVector<Double> conc;

         /** \brief Setup conc and rate constants
           @return true, if no concentration is negative and rate
                   constants were computed
         */

         bool setup (const SxVector<Double> &concIn,
                     const SxDefectConc &defconc,
                     double eFermi,
                     double ekt);

      protected:
         /// Time
         double currentTime;
         /// Evolve a few steps
         void evolve (double dt, int nSteps);

      public:
         /// Evolve for some time
         void evolve (double time);

         /// Empty constructor
         ReactionSystem () { /* empty */}
         /// Constructor
         ReactionSystem (const SxSymbolTable *table,
                         const SxDefectConc &defconc)
         {
            read (table, defconc);
         };
         /// Read from symbol table
         void read (const SxSymbolTable *table, const SxDefectConc &defconc);

   };

   /// Possible reactions
   ReactionSystem reactions;

   public:
      /// Energy bands for holes/electrons
      SxList<Band> bands;

      /// Dielectric constant (in units of e^2/[L][E])
      double dielecConstant;

      /// Chemical species
      SxList<SxString> elementNames;
   protected:
      /// List of defects
      SxList<Defect> defects;

      /// Boltzmann constant (in same units as defect energies)
      double kBoltzmann;

   public:

      /// Constructor from symbol table
      SxDefectConc (const SxSymbolTable *);
      
      /// Get number of species
      int getNSpecies () const { return (int)elementNames.getSize (); }
      /// Find a species
      int findSpecies (const SxString &name) const;
      /// Find a defect by name
      int findDefect (const SxString &name) const;
      /// Number of defects
      int getNDefects () const { return (int)defects.getSize (); }

      /// Print
      void print () const;
      /// Print defect concentrations
      void printConcentrations (double eFermi) const;
      
      /// Temperature
      double temperature;
      /// Absolute site concentration
      double siteConc;
      /// Chemical potentials
      SxVector<Double> mu;

      /// List of mobile species
      SxArray<int> mobileSpecies;


      /** \brief Compute thermodynamical equilibrium charge density at 
          given Fermi energy
      */
      double getChargeDensity (double eFermi) const;
      
      /** \brief Compute thermodynamical equilibrium charge density at 
          given Fermi energy
      */
      double getChargeDensity (double eFermi, 
                               const SxVector<Double> &conc) const;
      
      /** \brief Compute thermodynamical equilibrium charge density 
          from free electrons & holes at given Fermi energy
      */
      double getFreeChargeDensity (double eFermi) const;

      /// Class defining a range
      class Gap {
         public: 
            double from;
            double to;
            Gap (double a, double b) : from(a), to(b) {}
      };
      /// Get range of energy gap (for Fermi searches)
      Gap getGapRange () const;
         
      /// Find the thermodynamical equilibrium Fermi energy
      double getEFermi () const;

      /// Get kT
      inline double getkT () const { return kBoltzmann * temperature; }

      /// Get defect iterator begin 
      inline SxList<Defect>::ConstIterator begin () const 
      { 
         return defects.begin (); 
      }
      
      /// Get defect iterator end
      inline SxList<Defect>::ConstIterator end () const 
      { 
         return defects.end (); 
      }

      /** Get concentrations of species 
          @param defectConc concentrations of defects
          @param active  true for chemically active defects, false for
                         inactive ones
          @return total species concentrations from active defects
       */
      SxVector<Double> 
      getConcSpecies (const SxVector<Double> &defectConc,
                      const SxArray<bool> &active) const;


      /** \brief Compute equilibrium concentrations of active defects
          @param eFermi  Fermi energy
          @param chemPot adjusted chemical potentials
          @param active  true for chemically active defects, false for
                         inactive ones
          @param initialConc concentrations for inactive defects
          @return concentrations for all defects

          \note: The concentrations are summed over charge states.
        */
      SxVector<Double> getDefConc(double eFermi,
                                  const SxVector<Double> &chemPot,
                                  const SxArray<bool> &active,
                                  const SxVector<Double> &initialConc) const;

      /** \brief Compute equilibrium concentrations of active defects
          @param eFermi  Fermi energy
          @param active  true for chemically active defects, false for
                         inactive ones
          @param initialConc starting concentrations
          @return concentrations for all defects
          \note: The concentrations are summed over charge states.
        */
      SxVector<Double> 
      equilibrateDefects (double eFermi,
                          const SxArray<bool> &active,
                          const SxVector<Double> &cSIn,
                          bool saveMu = false);

      /** \brief Find chemical potentials for prescribed composition
          @param eFermi  Fermi energy
          @param active  true for chemically active defects, false for
                         inactive ones
          @param cSIn    desired species concentrations from active defects
          @return chemical potentials
        */
      SxVector<Double> 
      chemicalEquilibrium (double eFermi,
                           const SxArray<bool> &active,
                           const SxVector<Double> &cSIn,
                           const SxArray<int> ignoreSpecies = SxArray<int> ())
      const;

      /// \brief Execute equilibrate {} command
      double equilibrate (const SxSymbolTable *cmd);
};

inline double SxDefectConc::Defect::getChargeDensity (double kT, 
                                               double fermiEnergy, 
                                               double cTotal) const
{
   // --- this is a low-speed version (many vector allocations)
   // SxVector<Double> E = (E0 + charge * fermiEnergy);
   // E += -E.minval ();
   // SxVector<Double> conc = exp(-E / kT);
   // conc *= 1. / conc.sum (); // renormalize to 1
   // return (conc * charge).sum () * cTotal;
   
   // --- this is a high-speed version
   double eMin = E0(0) + fermiEnergy * charge(0), E;
   for (int ics = 1; ics < charge.getSize (); ++ics)
      if ((E = E0(ics) + fermiEnergy * charge(ics)) < eMin) eMin = E;

   double cCharge = 0., cSum = 0., conc;
   double beta = 1. / kT;
   for (int ics = 0; ics < charge.getSize (); ++ics)  {
      conc = exp (beta * (eMin - E0(ics) - charge(ics) * fermiEnergy));
      cSum += conc;
      cCharge += charge(ics) * conc;
   }
   return cCharge/cSum * cTotal;
}

#endif /* _SX_DEFECT_CONC_H_ */
