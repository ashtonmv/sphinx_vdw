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

#ifndef _SX_TIME_CORR_H_
#define _SX_TIME_CORR_H_

#include <SxAtomicStructure.h>
#include <SxSpeciesData.h>
#include <SxElemDB.h>
#include <fstream>
#include <SxFFT1d.h>
#include <SxStruct.h>

/** \brief SxTimeCorr 
    a tool to analyses time correlation functions and covariance 
    to extract thermodynamic averages etc.
    \ingroup   group_structure
    \author    Lars Ismer, ismer@fhi-berlin.mpg.de */
class SX_EXPORT_STRUCT SxTimeCorr
{
   public:
      SxList<SxAtomicStructure> trajectoryL;
      SxArray<SxAtomicStructure> trajectory;
      SxArray<SxAtomicStructure> trajectoryN;
      SxArray<SxVector<Complex16> > trajectoryNF;
      SxMatrix<Double> dynamicalMatrix, covarianceMatrix;
      SxMatrix<Double> normalModes;
      SxVector<Complex16> normalFrequencies;
      SxList<double> tL;
      SxArray<double> time;
      double startTime, endTime, dt;
      int nSteps, startIndex, endIndex, length, nDoF;
      SxString scheme, type;
      bool needsUpdateRC, needsUpdateRN, needsUpdateNF,
           needsUpdateDM, needsUpdateCM, needsUpdateModes, 
           needsUpdateFreqs;
      SxFFT1d fft;
      SxArray<double> masses;

      //---------------------------------------------------------------------
      /**@name Constructors and Destructors */
      //---------------------------------------------------------------------
      //@{
      SxTimeCorr ();
      ~SxTimeCorr ();
      /**\brief gets pair correlation function for a given pair of DoF's */
      SxArray<double> getPairCorrelationFunctionDoF (int, int);
      /**\brief gets pair covariance for a given pair of DoF's*/
      double getPairCovarianceDoF (int, int);
      /**\brief gets pair correlation function for a given pair of Atoms */
      SxArray<double> getPairCorrelationFunctionAtom (int, int);
      /**\brief gets auto-correlation function for a given Atom */
      SxArray<double> getAutoCorrelationFunction (int);
      /**\brief gets averaged (over all degrees of freedom) auto-correlation*/
      SxArray<double> getAveragedAutoCorrelationFunction ();
      /**\brief gets generalized frequency spectrum (quantity 
        should be a velocity)*/
      SxArray<double> getGeneralizedFrequencySpectrum ();
      /**\brief gets cross correlation matrix 
        (dynamical matrix for velocities)*/
      SxMatrix<Double> getCrossCorrelationMatrix ();
      /**\brief gets frequency resolution in 1/cm 
        (determined by observation time)*/
      double getDFreq (); 
      /**\brief pushes a new step to the trajectory*/
      void push (const SxAtomicStructure &, double);
      /**\brief loads trajectory from a moldynHist.dat file*/
      void loadTrajectory (const SxString &, const SxString &);
      /**\brief gets index in trajectory for a given real time (in au's)*/
      int getIndex (double);
      /**\brief sets up the type of quantity (structure, velocity, etc.)*/
      void setType (SxString);
      /**\brief sets up the observation range (in time units)*/
      void setObservationRange (int, int); 
      /**\brief gets Covariance Matrix*/
      SxMatrix<Double> getCovarianceMatrix ();
      /**\brief gets normal modes (quasi harmonic approximation) */
      SxMatrix<Double> getNormalModes ();

      /**\brief checks hermiticity (if the input hermitian)
          TODO: this function should be transfered to class SxMatrix*/
      double checkHermiticity (const SxMatrix<Complex16> &);
      /**\brief sorts the eigensystem2 using the similarity to eigensystem 2
          TODO: should be transferred to another class (which one ??)*/
      void sortByEigenvectors 
       (const SxMatrix<Double>::Eigensystem &, SxMatrix<Double>::Eigensystem *);
      
      /**\brief returns effective Hessian, requires as input a 
           velocity trajectory  
          (first argument) and a force trajectory (second argument)
          Ref: To be placed (URGENT, L.Ismer 03/17/05)*/
      SxMatrix<Double> getHessian (SxTimeCorr &, SxTimeCorr &);
      /**\brief sets up the masses*/
      void setMasses (const SxVector<Double> &); 
  
   protected:
      //--- checks current fft mesh and renews it if necessary
      void renewFFT (int, double);
      //--- updates array storage of data
      void updateTrajectory ();
      
      
      //--- NON-TESTED ROUTINES (PROJECT UNDER DEVELOPMENT)
      //   L.Ismer 03/17/05
      //--- update trajectory in normal mode projection
      ///void updateTrajectoryN ();
      //--- update fourier normal mode projection
      //void updateTrajectoryNF ();
      //--- gets fourier transformed normal mode projection of trajectory
      //SxArray<double> getNFProjection (int);
      //--- gets effective Dynamical Matrix (quasi harmonic approximation)
      //SxMatrix<Double> getDynamicalMatrix ();
      //--- update normal Frequencies (quasi harmonic approximation)
      //SxVector<Complex16> getNormalFrequencies ();

      //@}
};


#endif /* _SX_TIME_CORR_H_ */
