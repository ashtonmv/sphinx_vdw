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
#include <SxRadGBasis.h>
#include <SxAtomicOrbitals.h>
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

typedef SxArray<SxDiracMat<Complex16> > SxOrbitals; //(ik,iOrbital:ig)

class SX_EXPORT_EXT SxQuamolRedSubSpace   {
   public:
      /// \brief Standart Constructor
      SxQuamolRedSubSpace ();
      /// \brief Standart Destructor
      ~SxQuamolRedSubSpace ();
      /// \brief set function for jsbBasis
      void set ( SxString fileName,
                 SxConstPtr<SxRadBasis> radBasisPtrIn,
                 SxPW &waveIn,
                 SxGkBasis *GkPtrIn,
                 SxConstPtr<SxRadGBasis> radGBasisPtrIn,
                 SxPtr<SxOverlapBase> SPtrIn,
                 SxAtomicStructure structureIn,
                 double rCut);
      /// \brief set Function for Basis reduction
      void set ( SxString fileName,
                 SxAtomicOrbitals initGuess,
                 SxConstPtr<SxRadBasis> radBasisPtrIn,
                 SxPW &waveIn,
                 SxGkBasis *GkPtrIn,
                 SxConstPtr<SxRadGBasis> radGBasisPtrIn,
                 SxPtr<SxOverlapBase> SPtrIn,
                 SxAtomicStructure structureIn);
      /// \brief QUasi AtoMic OrbitaLS (QUAMOLS) in radBasis
      SxAtomicOrbitals basis,improvedBasis;

      SxOrbitals basisGk;

      /// \brief Waves
      SxPW waves;
      /// \brief Plane-Wave Quamol Projection \f$\langle\Psi|b\rangle \f$
      SxArray<SxDiracMat<Complex16> > beta,
      /** \brief Expansion coefficients matrix 
         \f$ \beta_{m,\mu} = \sum\limits_n{[S^b_{m,n}]^{-1}
                       \langle b_n|c_\mu \rangle} \f$ 
      */             
                          Sb;
      SxArray<SxArray<SxDiracMat<Complex16> > > P;

      bool restartCG;

                          

      SxAtomicStructure structure;

      /// \brief RadBasisPtr
      SxConstPtr<SxRadBasis> radBasisPtr;
      /// \brief RadGBasisPtr;
      SxConstPtr<SxRadGBasis> radGBasisPtr;
      /// \brief GkBasisPtr
      SxGkBasis *GkBasisPtr;

      /// \brief Overlap Operator
      SxPtr<SxOverlapBase> SPtr;

      void info (SxString fileName, SxConstPtr<SxRadBasis> radBasisPtrIn);

      void checkInitialGuess ();

      SxArray<SxArray<SxDiracMat<Complex16> > > getProjections ();

      SxArray<SxDiracMat<Complex16> > getCoefficients ();

      SxArray<SxDiracMat<Complex16> > getOverlapC ();

      SxArray<SxDiracMat<Complex16> > getOverlapB (); 
      
      /** \brief Calculate the squarenorm: 
          \f$\sum\limits_{states,k}{w_k (P\beta[S^c]^{-1}\beta^{\dagger}P^{\dagger})_{states,k;states,k} 
      */   
      
      double calcNormB ();
      double calcNormC ();


      /** \brief calculate Gradient for SxOrbitals orbitals and transform back into SxRadials:
        \f$\sum\limits_{states,k}   {
         (P\beta[S^c]^{-1})_{states,k;\pi})P^{\dagger}_{p;staes,k}
         -(P\beta[S^c]^{-1})_{staes,k;\pi}
         (S^b\beta[S^c]^{-1}\beta P^{\dagger})_{p;states,k}  
         }\f$
      */   
      SxArray<SxDiracMat<Complex16> > calcGradient ();

      SxDiracMat<Complex16> keepSymmetry (const SxDiracMat<Complex16> &MatIN);

      /// \brief Update orbitals with (Gk|radials)
      SxOrbitals radialsToOrbitals (const SxAtomicOrbitals &radialsIN);

      //void betaToRadial ();

      int getOrbital(SxArray<SxQuantumNumbers> &map, int is, int n, int l) const;

      SxDiracMat<Complex16> reduce (SxArray<SxDiracMat<Complex16> > &matIn) const;

      SxArray<SxDiracMat<Complex16> > expand (SxDiracMat<Complex16> &matIn) const;

      void redBetaToRadial ();

      /*
      /// \brief Update radials with (rad|orbitals)
      SxAtomicOrbitals orbitalsToRadials (
            const SxOrbitals &orbitalsIN,
            const SxAtomicOrbitals &radialsIN);
            */

      /// \brief \f$\langle\mu|\nu\rangle\f$
      double dotproductRad (
            const SxDiracVec<Double> &vec1, 
            const SxDiracVec<Double> &vec2);

      double error;

      int maxSteps;

      void compute ();

      SxArray<SxDiracMat<Complex16> > getCopy (
            const SxArray<SxDiracMat<Complex16> > &dataIN);
      
      SxArray<SxDiracMat<Complex16> > getBeta ();

      void setBeta (const SxArray<SxDiracMat<Complex16> > &betaIN);

      SxArray<SxDiracMat<Complex16> > addBeta (
            const SxArray<SxDiracMat<Complex16> > &beta1,
            const SxArray<SxDiracMat<Complex16> > &beta2);

      SxArray<SxDiracMat<Complex16> > skalarMultBeta (
            const double skalar,
            const SxArray<SxDiracMat<Complex16> > &betaIN);
      
      SxComplex<double> lineMin (const SxArray<SxDiracMat<Complex16> > &dir, double sw);

      SxComplex<double> parabelFit (double,double,double,double,double,double);

      SxAtomicOrbitals getJSB (double rCut, int lMax, int nZeros);

      SxArray<SxVector<Double> > getJSBZeros (int lMax, int nZeros);

      double findRootJSB (int l, double guess);

      SxDiracVec<Double> jsb (int l, const SxDiracVec<Double> &z) const;

};

namespace Timer {
   enum QuamolTimer {
      normCalc,
      gradCalc,
      rad2Orb,
      setup,
      overlapB,
      projections
   };
}

SX_REGISTER_TIMERS (Timer::QuamolTimer)
{
   using namespace Timer;
   regTimer (setup,      "Basis setup");
   regTimer (rad2Orb,    "Radial to Orbitals");
   regTimer (overlapB,   "<b_i|b_j>");
   regTimer (projections,"<Psi|b_i>");
   regTimer (normCalc,   "Norm Calculation");
   regTimer (gradCalc,   "Gradient Calculation");
}

#endif /* _SX_QUAMOL_H_ */
