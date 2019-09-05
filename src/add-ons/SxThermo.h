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

#ifndef _SX_THERMO_H_
#define _SX_THERMO_H_

#include <SxString.h>
#include <SxVector.h>
#include <SxMatrix.h>
#include <SxMath.h> 
#include <SxArray.h>
#include <SxAtomicStructure.h>
#include <SxMurn.h>
#include <SxConstants.h>
#include <SxExt.h>


/** \brief Compute thermodynamical properties: (quasi-harmonic) free energies,
 * lattice expansion, bulk modulus & derivative, phononDos, heat capacities
 * (C_V & C_p) from the provided frequencies (in sxb format, multiple volumes). 

 * For the quasi-harmonic approximation murnaghan data must be provided.
 * Phonon extrapolation is supported and two types of quasi-harmonic
 * fitting: numerical and using the murnaghan fits for all T. Extra free energy
 * contributions (electronic / magnetic) can be read in and added to the
 * vibrational free energy.

    \b SxThermo = S/PHI/nX Compute Thermodynamic properties

    \author Matth\'e Uijttewaal, uijttewaal@mpie.de 
 */
      
      //square of number, internal routine: TODO move to SxMath
      inline double sqr (const double &x) {return x*x;}

      //cube of number, internal routine: TODO move to SxMath
      inline double cube (const double &x) {return x*x*x;}

class SX_EXPORT_EXT SxThermo {
   protected:
      //log of vector, also counts and removes negative elements, internal 
      SxVector<Double> logVec (const SxVector<Double> &y) const;

   public:
      // empty constructor
      SxThermo ()  { /*empty*/ }
      
      // destructor
     ~SxThermo ()  { /*empty*/ }

     /**\brief read frequencies from 1 input file

       @param input input file
       @param vol volume per atom to which frequencies correspond (B^3)
       @return vector of the frequencies
       */
      SxMatrix<Double> readFreq (const SxString &input, 
                                       double &vol) const;
      
           /**\brief extrapolate phonons to other volumes (in B^3/atom)

             Either a linear or an exponential extrapolation can be done.
             \note exponential extrapolation is only to positive frequencies and
             identical qPoint grids are assumed!
        @param allFreq frequencies for input volumes (meV)
        @param vol the volumes in B^3/atom (input + to extrapolate)
        @param lin linear interpolation?
        @return frequencies, input + extrapolated (meV)
        */
      SxArray<SxMatrix<Double> > getXtraPhon (
            const SxArray<SxMatrix<Double> > &allFreq,
            const SxArray<double> &vol,
            const bool lin) const;
  
      /**\brief calculate and print the phonon density of states per atom

        \note it rangesfrom 0-40 meV and has on average 10 freq. per interval
       @param vol the volumes (B) to which the frequencies belong
       @param freq all the frequencies
       @param output the output file
       */
      void printDOS (const SxArray<double> &vol, 
                     const SxArray<SxMatrix<Double> > &freq) const;

      /** \brief get the free energy as a quantum statistics sum of
       * the frequencies (meV->H) for 1 volume

        @param freq frequencies (meV) 
        @param maxT maximal temperature (K)
        @param dT temperature step (K)
        @param ignore how to treat imag. freq.
        @return free Energy (H/atom) for T in maxT/dT steps to maxT (K) 
       */
      SxVector<Double> getFreeEn (const SxMatrix<Double> &freq,
                                  const int maxT,
                                  const double dT,
                                  const bool ignore) const;  

      /** \brief read extra free energy contribution (H/atom) when provided

        @param input filename
        @return array of vector with free energies, experimental :(nVol,nT)
        */
      SxArray<SxVector<Double> > readXtraF (const SxString &input) const;

      /**\brief print free energies and C_V for all volumes

        @param allFreeEn free energy (H/atom) for all vol and T
        @param vol the volumes (B^3)
        @param dT temperature step (K)
        @param output the output file
        */
      void printF (const SxArray<SxVector<Double> > &allFreeEn,
                   const SxArray<double> &vol,
                   const double dT,
                         ofstream &output) const;

       /**\brief read murnaghan data input (vol's + Etot's/atom in a.u.)

        @param input the input file
        @param BM bulkModulus to replace
        @param BMD BM derivative to replace
        @return the murnaghan data
        */ 
      SxMurn readMurnDat (const SxString &input,
                          const double &BM,
                          const double &BMD) const;

      /* \brief make an (analytic) 2nd order parametrization of the free energy
       *  at all T

        @param allFreeEn free energies (H/atom, :nVol, nT)
        @param vol volumes/atom (B^3, :nVol)
        @param linF reduce to linear fit
        @return coords with parameters @ all T (H, B)
        \note 2B put in class linear regression
         */
      SxArray<Coord> paramFreeEn (const SxArray<SxVector<Double> > &allFreeEn,
                                  const SxArray<double> &vol,
                                  const bool linF) const; 

      /* \brief take second deriv. of paramF to get param. heat capacity (k_B/B^3)

         @param paramF parametrised free energy for all T (H, B^3)
         @param dT temperature step (K)
         @return parametrised heat capacity (B^3)
       */
      SxArray<Coord> paramHeat (const SxArray<Coord> &paramF,
                                const double dT) const;

      /**\brief make murnaghan fit at all T (in a.u.)

        @param inMurn 0K murnaghan fit 
        @param paramF parametrised vibrational free en. (H/atom, B^3) for all T
        @param vol the volumes (B^3/atom) for which the fit should be made
        @param press the current pressure (H/B^3)
        @return murnaghan fits for all T
        */
      SxArray<SxMurn> fitMurn (SxMurn &inMurn, 
                               const SxArray<Coord> &paramF,
                               const SxArray<double> &vol,
                               const double &press) const;

   protected:
      /* \brief get the equilibrium volume numerically at one temp & pressure,
       * internal routine

         @param inMurn total energy murn fit
         @param paramF parametrised vib. free en. (H, B^3)
         @param press the provided pressure (H/B^3)
         @return equilibrium volume (B^3/atom)
       */
      double getEqVol (SxMurn &inMurn, 
                   const Coord &paramF,
                   const double &press) const;

   public:
        /**\brief make murnaghan parametrisation at all T (in a.u.)

        @param inMurn total E murnaghan fit 
        @param paramF parametrised vib. free en. (H/atom) for all T
        @param press the current pressure (H/B^3)
        @return murnaghan fits for all T
        */
      SxArray<SxMurn> paramMurn (SxMurn &inMurn, 
                   const SxArray<Coord> &paramF,
                   const double &press) const;

      /* \brief get the heat capacity at constant pressure 
       * from that at constant volume and the lattice expansion

         formula: C_p = C_V(p) -T \frac{\partial ^2 F}{\partial T \partial V}
                        \cdot \frac{\partial V}{\partial T}_p
                      = C_V(iT, p) - iT / (K_B dT) \cdot \left(
                      F (iT, V0(iT, p)) + F (iT-1, V0(iT-1, p)) 
                      - F (iT-1, V0(iT, p)) - F (iT, V0(iT-1, p)) \right)
         @param allMurn the total energy murn fit (one pressure)
         @param Cv the parametrised heat capacity (B^3)
         @param dT temperature step (K)
         @return heat capacity at one press (for all T)
       */
      SxArray<double> getCp (SxArray<SxMurn> &allMurn, 
                        const SxArray<Coord> &Cv,
                        const double dT) const;

      /**\brief print heat capacities, internal routine

         general formulae:  C_p=\left(\frac{\partial U + PV}{\partial T}\right)_p
        //==             g(T,p)=f+pv=f-\partial f \partial v_T v.!! 
                    C_V=\left(\frac{\partial U}{\partial T}\right)_V
        
        @param Cv parametrised heat capacity at constant volumes (B^3)
        @param CvXF parametrised extra contribution to the heat capacities
        @param vol the volumes (B^3/atom) for which Cv is printed 
        @param Cp heat capacity at constant pressures (kT/atom)
        @param unitPress the unit of pressure in (H/B^3)
        @param dT the temperature step (K)
        */
      void printHeat (const SxArray<Coord> &Cv,
                      const SxArray<Coord> &CvXF,
                      const SxArray<double> &vol,
                      const SxArray<SxVector<Double> > &Cp,//:nT,nPres
                      const double &unitPress,
                      const double dT) const;
 
      /**\brief print quasi-harm. thermodynamic quantities (free en, vol,
       * linear lattice expansion \alpha, BM, BM')

        @param allMurn murnaghan data for all temperatures
        @param param parametrized free energy :nT
        @param paramXF parametrized extra free energy contribution :nT
        @param dT temperature step (K)
        @param output the output file
        /note the Gr\"uneisen parameter is defined as 
              \gamma = (C_p/C_V-1)/(3\alpha T) 
              and can simply be got from the current output.
        */
      void printThermo (SxArray<SxMurn> &allMurn,
                        SxArray<Coord> &param,
                        SxArray<Coord> &paramXF,
                               const double dT,
                               ofstream &output) const;
};
#endif /* _SX_THERMO_H_ */

