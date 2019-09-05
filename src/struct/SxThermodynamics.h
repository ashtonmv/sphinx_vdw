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

#ifndef _SX_THERMODYNAMICS_H_
#define _SX_THERMODYNAMICS_H_

#include <SxArray.h>
#include <SxConstants.h>
#include <SxPrecision.h>
#include <SxVector.h>
#include <SxSecondaryStructure.h>
#include <SxHessianOps.h>
#include <SxStruct.h>


/**
  @brief A Class to deal with thermodynamic properties as derived from 
         partition functions for ionic movement
         Up to now: Only constant Volume, zero pressure harmonic internal 
                    vibrational partition function considered
         Units (up to now): input frequencies in 1/cm, Kelvin
                            output quantities in kcal/mol, Kelvin 
         under development: quasi-harmonic approximation 
         References:  L. Ismer, Diploma Thesis, TUB (2002)     
                      J. Xie et al. PRB 59, 965 (1999)         
  */
class SX_EXPORT_STRUCT SxThermodynamics 
{
      
   public:
      SxThermodynamics ();
      SxThermodynamics (const SxVector<Complex16> &, int, int, double);
      SxThermodynamics (SxSecondaryStructure &, const SxString &, double);
      SxThermodynamics (SxHessianOps &, int, double);
      ~SxThermodynamics ();
      /**\brief sets harmonic frequency spectrum  
                 (to get the partition function )*/
      void setFrequencies (const SxVector<Complex16> &);
      /**\brief sets the number of external degrees of freedom*/
      void setNExternalDoF (int);
      /**\brief sets the threshold 
               (frequencies below this value are excluded)*/
      void setThreshold (double);
      /**\brief gets the quantum harmonic vibrational entropy */
      double getSVib (double, double);
      /**\brief gets the classical harmonic vibrational entropy */
      double getSVibCl (double, double);
      /**\brief gets the quantum harmonic vibrational heat capacity*/
      double getCVVib (double, double);
      /**\brief gets the classical harmonic vibrational heat capacity*/
      double getCVVibCl (double, double);
      /**\brief gets the quantum harmonic vibrational internal energy*/
      double getUVib (double, double);
      /**\brief gets the classical harmonic vibrational internal energy*/
      double getUVibCl (double, double);
      /**\brief gets the quantum harmonic vibrational free energy*/
      double getFVib (double, double);
      /**\brief gets the classical harmonic vibrational free energy*/
      double getFVibCl (double, double);
      /**\brief gets zero point vibrational energy*/
      double getZPV (double);
      /**\brief (approximately) transforms classical temperature into 
                quantum meachinacal Temperature*/
      double getRescaledTemp (double, double);
      /**\brief gets total energy for a given volume from Murnaghan equ. */ 
      double getEMurn (double);
      /**\brief gets first derivative of total energy 
                for a given volume from Murnaghan equ. */ 
      double getdEdVMurn (double);
      /**\brief gets second derivative of total energy 
                for a given volume from Murnaghan equ. */ 
      double getd2EdV2Murn (double);
      /**\brief gets first derivative of vibrational free energy 
                for a given volume in  grueneisen appr. */ 
      double getdFvibdV (double, double);
      /**\brief gets second derivative of vibrational free energy 
                for a given volume in grueneisen appr. */ 
      double getd2FvibdV2 (double, double);
      /**\brief gets derivative of free energy 
                for a given volume in grueneisen appr. */ 
      double getdFdV (double, double);
      /**\brief solves the equation of state, i.e. returns the
                relaxed volume for a given temperature and pressure*/
      double solveEquOfState (double, double);
     /**\brief sets up parameters for the grueneisen formalism*/
      void setGrueneisen (const SxSymbolTable *);
      /**\brief prints the thermodynamic data into a file*/
      void print (const SxString &, double, double, int);
      

   protected:
      /**\brief contains frequencies in 1\cm*/
      SxVector<Complex16> freqs;
      /**\brief contains permutation of freq-vector */
      SxVector<Int> perm;
      /**\brief ...*/
      int nDoF;
      /**\brief frequencies below threshold are not counted in the
                partition function */
      double threshold;
      /**\brief number of external degrees of freedom (are zero, not counted)*/
      int externalDoF;
      /**\brief norm (the thermodynamic data is always normed to to the 
                supercell*/
      double norm;
      /**\brief sets the norm */
      void setNorm (double);
      /**\brief defines whether the grueneisen-theory is applied or not*/
      bool grueneisen;
      /**\brief static equilibrium volume */
      double V0;
      /**\brief Bulk modulus at T=0 K */
      double B0;
      /**\brief first derivative of Bulk modulus at T=0 K */
      double B0Prime;
      /**\brief Grueneisen coefficients (1/w*dw/dV) */
      SxVector<Double> GECoeff;
      
        
};
      
#endif /* _SX_THERMODYNAMICS_H_ */
