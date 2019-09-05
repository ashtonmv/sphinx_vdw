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

#ifndef _SX_QUAMOL_H_
#define _SX_QUAMOL_H_

#include <SxUtil.h>
#include <SxDirac.h>
#include <SxVector.h>
#include <SxMatrix.h>
#include <SxParser.h>
#include <SxAtomicStructure.h>
#include <SxPseudoPot.h>
#include <SxGkBasis.h>
#include <SxRadBasis.h>
#include <SxRadRBasis.h>
#include <SxRadGBasis.h>
#include <SxAtomicOrbitals.h>
#include <SxAtomicOrbitalsR.h>
#include <SxAtomicOrbitalsG.h>
#include <SxQuantumNumbers.h>
#include <SxProjector.h>
#include <SxConstants.h>
#include <SxTimer.h>
#include <SxError.h>
#include <SxCLI.h>
#include <SxFermi.h>
#include <SxPW.h>
#include <SxYlm.h>
#include <SxPWOverlap.h>
#include <SxPAWOverlap.h>
#include <SxExt.h>
#include <SxCubicSpline.h>

typedef SxDiracMat<Complex16> SxOrbitals; //(iOrbital:ig)
typedef SxDiracMat<Complex16> SxDMatC16; //(iOrbital:ig)
typedef SxDiracVec<Complex16> SxWaveState; //(ig)

class SX_EXPORT_EXT SxQuamol   {
   public:
      /// \brief Standart Constructor
      SxQuamol ();
      /// \brief Standart Destructor
      ~SxQuamol ();
      /// \brief set Function
      void set ( SxAtomicOrbitals &radialsIn,
                 SxConstPtr<SxRadGBasis> radGBasisPtr,
                 SxConstPtr<SxPW> wavePtrIn,
                 SxConstPtr<SxFermi> fermiPtrIn,
                 SxConstPtr<SxOverlapBase> SPtrIn);

      // Variables
      /// \brief QUasi AtoMic OrbitaLS (QUAMOLS) in radGBasis representation
      SxAtomicOrbitalsG radials;

      double functionalValue;
      SxAtomicOrbitalsG grad;

      double spillage;
      SxAtomicOrbitalsG spillageGrad;
      double sigma;

      double eKin;
      SxAtomicOrbitalsG eKinGrad;
      double zeta;
      bool adaptiveZeta;

      double locVal;
      SxAtomicOrbitalsG locGrad;
      double kappa;
      double rStart;
      bool adaptiveKappa;

      /// \brief maximal radial extension from Radialbasis of initial guess
      double rMax;
      void computeLocRGGprime (SxConstPtr<SxRadRBasis> radRBasisPtr);
      SxList<SxDiracMat<Double> > locRGGPrime;

      /// \brief Waves
      SxConstPtr<SxPW> wavesPtr;
      /// \brief Fermi
      SxConstPtr<SxFermi> fermiPtr;
      /// \brief Overlap Operator
      SxConstPtr<SxOverlapBase> SPtr;
      /// \brief Convergence criterium for functional change
      double dF;
      /// \brief Convergence criterium for gradient length
      double dRes;
      /// \brief Maximum number of steps for spillage optimization
      int maxSteps;
      double relLineMin;
      /** \brief Flag to print radials gradients and directions 
           in every optimization step
      **/
      bool print;

      void computeFunctional (const SxAtomicOrbitalsG &functions, 
                              bool calcGrad = true, bool calcKinLoc = true);
      void computeSpillage (const SxAtomicOrbitalsG &functions, 
                            bool calcGrad = true);
      void computeKineticEnergy (const SxAtomicOrbitalsG &functions, 
                            bool calcGrad = true);
      void computeLocalization (const SxAtomicOrbitalsG &functions, 
                            bool calcGrad = true);
      SxAtomicOrbitalsG calcNumericGrad (const SxAtomicOrbitalsG &functions, 
            double h = 1e-6);
      double lineMin (const SxAtomicOrbitalsG &dir,
                      const SxVector3<Double> &x,
                      SxVector3<Double> &y);
      SxComplex<double> parabelFit (double x0, double x1, double x2,
                                        double N0, double N1, double N2);
      /// \brief Expand radialsG into Gk Basis
      SxOrbitals expandRadialsG (const SxAtomicOrbitalsG &functions, int ik);
      SxAtomicOrbitals getOrbitals (const SxAtomicOrbitalsG &functions,
            SxConstPtr<SxRadBasis> radBasisPtr);
      SxAtomicOrbitalsR getOrbitals (const SxAtomicOrbitalsG &functions,
            SxConstPtr<SxRadRBasis> radRBasisPtr);

      void checkTrafo (const SxAtomicOrbitalsG &functionsIn);

      void refine (SxArray<int> &factor);

      SxAtomicOrbitalsG setOrthogonal (
            const SxAtomicOrbitalsG &functions, 
            const SxAtomicOrbitalsG &referenz);

      SxAtomicOrbitalsG setOrthogonal (
            const SxAtomicOrbitalsG &functions);

      void completenessProfile(const SxAtomicOrbitalsG &functions);
      void getNormContribution(const SxAtomicOrbitalsG &functions,
            SxConstPtr<SxRadBasis> radBasisPtr);

      /// \brief optimization routine
      void compute ();

      int printStep;

      bool printLine;

      bool checkGrad;

      SxArray<SxArray<SxDiracMat<Double> > > dRgdRGk;

      void setDRDR (int ik);

      SxArray<SxArray<bool> > fixedOrbitals;

      void setFixedList(const SxSymbolTable *table);

      void printSpillagePerState (const SxAtomicOrbitalsG &functions);

      void updateGuess (SxVector3<Double> &x, SxVector3<Double> &y, double val);
      
      void calcResidues (const SxAtomicOrbitalsG &functions);

};

namespace Timer {
   enum QuamolTimer {
      setup,
      computeFunctional,
      computeSpillage,
      computeKineticEnergy,
      computeLocalization,
      calcNumericGrad,
      stepClock,
      lineMin,
      radG2mu,
      checkTrafo,
      getROrbitals,
      DRDataDRGrid
   };
}

SX_REGISTER_TIMERS(Timer::QuamolTimer)
{
   using namespace Timer;
   regTimer (setup,                "Setup");
   regTimer (computeFunctional,    "compute F");
   regTimer (computeSpillage,      "compute Spillage");
   regTimer (computeKineticEnergy, "compute Ekin");
   regTimer (computeLocalization,  "compute Loc");
   regTimer (calcNumericGrad,      "numeric gradient");
   regTimer (stepClock,            "Step Clock");
   regTimer (lineMin,              "Line Minimization");
   regTimer (radG2mu,              "RadG -> orbitalsG");
   regTimer (checkTrafo,           "trafo check");
   regTimer (getROrbitals,         "get orbitals in R");
   regTimer (DRDataDRGrid,         "DRData-Grid");
}

#endif /* _SX_QUAMOL_H_ */
