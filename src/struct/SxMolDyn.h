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

#ifndef _SX_MOL_DYN_H_
#define _SX_MOL_DYN_H_

#include <SxAtomicStructure.h>
#include <SxSymbolTable.h>
#include <SxPotential.h>
#include <SxSpeciesData.h>
#include <SxTimeCorr.h>
#include <SxArtifactFilter.h>
#include <SxTaylorExpPotential.h>
#include <SxHamSolver.h>
#include <SxStruct.h>

/** \brief Molecular dynamics

    \b SxMolDyn = S/PHI/nX Molecular Dynamics

    ...

    \ingroup   group_structure
    \author    Lars Ismer, ismer@fhi-berlin.mpg.de */
class SX_EXPORT_STRUCT SxMolDyn
{
   public:

      //---------------------------------------------------------------------
      /**@name Constructors and Destructors */
      //---------------------------------------------------------------------
      //@{
      SxMolDyn (const SxAtomicStructure &, SxPotential *);
      SxMolDyn ();
      ~SxMolDyn ();
      //@}
		class SxIntegrator;
		friend class SxMolDyn::SxIntegrator;
		class SxInit;
		friend class SxMolDyn::SxInit;
		class SxOutput;
		friend class SxMolDyn::SxOutput;
		class SxPerformance;
		friend class SxMolDyn::SxPerformance;
		class SxPerturbation;
		friend class SxMolDyn::SxPerturbation;
		class SxState;
		friend class SxMolDyn::SxState;
      class SxGuess;
      friend class SxMolDyn::SxGuess;
      class SxImpSampler;
      friend class SxMolDyn::SxImpSampler;
		
		class SxIntegrator {
			public:
				SxIntegrator ();
				SxIntegrator (const SxSymbolTable *, SxMolDyn *);
				~SxIntegrator ();
				void printParameters ();
				void validateInput ();
				void resetGearPredictorMatrix (double);
				void initialize (SxState *, SxState *);
				void initPCABM ();
				void initVerlet ();
				void predictNewState(SxState *, SxState *);
				void correctNewState(SxState *, SxState *);
				
				SxMolDyn* parent;
				// --- control parameters as readed in from input.sx
				int order, timeSteps;
				double dt;
				SxString scheme;
				// --- variables needed for Gear Predictor-Corrector scheme
				SxMatrix<Double> gearPredictorMatrix;
				SxVector<Double> gearCoeff;
				SxVector<Double> correctionGear;
				SxMatrix<Double> gearDerivatives;
				SxVector<Double> thermoDerivatives;
            // --- variable needed for the BBK scheme
		};


		class SxTDIntegration {
			public:
				SxTDIntegration ();
				SxTDIntegration (const SxSymbolTable *, SxMolDyn *);
				~SxTDIntegration ();
				void printParameters ();
				void validateInput ();
				void initialize ();
				void apply (SxGuess *);
				void takeMeasurement (const SxState &);
				void goOn ();
				void switchOn ();
				SxMolDyn* parent;
				bool switchedOff;
				double sumSteps, dUdLAvg;
		};
      
		class SxThermostat {
			public:
				SxThermostat ();
				SxThermostat (const SxSymbolTable *, SxMolDyn *);
				~SxThermostat ();
				void printParameters ();
				void validateInput ();
				void initialize ();
				void apply (SxState *);
				void goOn ();
				void switchOn ();
            double getInitialEkin ();
				SxMolDyn* parent;
				double T, startTemp, endTemp, mass, storeLambda;
				int tempSteps, L, iSteps, nHarmEqu, iHarmEqu;
				SxString type, increment, scheme;
            double alpha, dT;
				bool switchedOff;
            //--- parameters for Langevin dynamics;
            double gammaLD, gammaLDStore, CLD, ELD, GLD, HLD;
            SxVector<Double> sigmaSqrPlu, sigmaSqrMin;
		};

		class SxPerformance {
			public:
				SxPerformance ();
				SxPerformance (const SxSymbolTable *, SxMolDyn *);
				~SxPerformance ();
				void printParameters ();
				void validateInput ();
            void initialize ();
            void update (const SxState &);
            void extrapolateWavefunction(const SxAtomicStructure &);
            int meditateRestart(int, SxState *, SxState *);
            double getAlpha (const SxAtomicStructure &);

            SxPW actWaves, oldWaves;
            SxMolDyn* parent;
				SxString freezeMode;
				SxMatrix3<Double> constraints;
				bool applyConstraints, extrapolateWaves, restart;
            SxAtomicStructure xOld, xAct, xNew;
		};
		
		class SxPerturbation {
			public:
				SxPerturbation ();
				SxPerturbation (const SxSymbolTable *, SxMolDyn *);
				~SxPerturbation ();
				void printParameters ();
				void validateInput ();

				SxMolDyn* parent;
				bool applyPerturbation;
				int dof, interval;
				SxString gas;
		};
		
		class SxOutput {
			public:
				SxOutput ();
				SxOutput (const SxSymbolTable *, SxMolDyn *);
				~SxOutput ();
				void printParameters ();
				void validateInput ();
				void writeEnergiesVsTime (SxState*);
				void writeState (SxState*);
				void writeRestartInfo (SxState *, SxState*, int);
				void printGFS ();
            void printHessian ();

				SxMolDyn* parent;
				FILE *fpMoldyn;
				FILE *fpMoldynHist;
            FILE *fpFreq;
            FILE *fpHessian;
				int saveWaves;
				bool saveHistory, saveRestartInfo;
            bool inSectionsFreq;
            int printSteps;
            int printStepsFreq;
            bool inSectionsHessian;
            bool plotDih;
            int printStepsHessian;
            int C_one, C_two, N_one, N_two, C_alpha;
		};

		class SxInit {
			public:
				SxInit ();
				SxInit (const SxSymbolTable *, SxMolDyn *);
				~SxInit ();
				void printParameters ();
				void validateInput ();
				SxMolDyn::SxState getInitialState ();
				SxMolDyn* parent;
				SxString scheme, initHistoryFN, initStructureFN, gas, initHessianFN;
				int dof, seed;
            bool fromTemperature, noRandom;
				double initEkin, deltaE;
		};
		
		class SxState {
			public:
				SxState ();
				SxState (SxMolDyn *);
				~SxState ();
				void updateVariables ();
				void plotDihedrals ();
            void set (const SxState &);

				SxMolDyn* parent;
				SxAtomicStructure X, V, F, A, XRand, XRandInv;
				double eKin, ePot, eTot, t, T, hamNH, dUdL, dUdLRef, E0Ref;
				double thermoX, thermoV, thermoA, thermoF;
            int iSteps;
		};
		
      class SxImpSampler {
			public:
				SxImpSampler ();
				SxImpSampler (const SxSymbolTable *, SxMolDyn *);
				~SxImpSampler ();
            
            void inc ();
            void apply ();
            void initialize ();
            
            SxMolDyn* parent;
            SxMolDyn::SxState startState, startOldState;
            SxMolDyn::SxIntegrator intOld;
            int iSteps;
            int nTrust;
            bool isOn;
            double ePotRS, ePotTS, ePotRF, ePotTF, lambda;
      };

		class SxGuess {
			public:
				SxGuess ();
				SxGuess (const SxSymbolTable *, SxMolDyn *);
				~SxGuess ();

				SxMolDyn *parent;
            SxPotential *potential;
            bool isActive;
            double lambda;
		};
      //---------------------------------------------------------------------
      /**@name Interface to input file
         Controlling the class by the \ref tutor_parser. */
      //---------------------------------------------------------------------
      //@{
      void print (const SxSymbolTable *);
		void execute (const SxSymbolTable *, bool calc=true);
		//@}
      static bool isRegistered (const SxSymbolTable *);
   
	protected:
	
		SxIntegrator integrator;
		SxThermostat thermostat;
		SxPerturbation perturbation;
		SxInit init;
		SxOutput output;
		SxPerformance performance;
		SxState state, oldState;
      SxGuess guess;
      SxTDIntegration TDInt;
      SxImpSampler impSampler;
			

      // --- variables needed by the MD
      SxTimeCorr timeCorrVel, timeCorrForce;
		SxAtomicStructure initialTau;
		SxSpeciesData speciesData;
		SxPotential *potential;  // TODO: should be const
      SxHamSolver *hamSolver;
		SxArray<const SxSymbolTable *> elMinimCmds;
		int nDoF;
		SxAtomicStructure mass; // SIXTEN: masses need to implemented
	   SxVector<TPrecTauR> massVec;
      bool useTimeCorr;
      SxArtifactFilter forceFilter, velocityFilter, randomFilter;
      SxPtr<SxTaylorExpPotential> tPotential, tPotential2;
		
      //--- routines needed by the MD
		double getEKin (SxAtomicStructure);
		double getEPot ();
      SxAtomicStructure getForces (const SxAtomicStructure &);
		
      int mdStep;
		
		//-- utility routines
		int getFaculty(int);
		double gaussDev   (double);
};


#endif /* _SX_MOL_DYN_H_ */
