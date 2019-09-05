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

#ifndef _SX_VIBRATIONS_H_
#define _SX_VIBRATIONS_H_

#include <SxString.h>
#include <SxSpeciesData.h>
#include <SxExt.h>

/** \brief sxvibrations
    This "add-on" is a tool to study harmonic vibrations and more
    - is in a pre-alpha-version state yet
    \ingroup  group_addons
    \author   Lars Ismer, ismer@fhi-berlin.mpg.de */
class SX_EXPORT_EXT SxVibrations
{
   public:
      SxVibrations ();
      ~SxVibrations ();
      SxSpeciesData speciesData;
      SxParser::Table table;

      /** \brief  user-defined filenames */
      SxString inputFile, inputHessian,
               inputStructure, inputStructure2, inputPeriodicity,
               inputPeriodicity2, inputEigenvalues, inputHistory,
               moldenFilename, anharmonicFilename, thermodynamicsFilename,
               eigenvaluesFilename, hessianFilename,
               dispersionFilename, couplingsFilename,
               inputBasisHessian, forceConstantsFilename,
               nuMuFilename;

      /** \brief user-defined class of the system (i.e. helix, molecule etc.)*/
      SxString systemClass;

      /** \brief user-defined: which degrees of freedom shall be considered
                 as free translations/rotations*/
      SxString freezeMode;
      /** \brief interpolation type for the phonon-dispersion relation*/
      SxString interpolation;
      /** \brief interpolation type for the pd relation tp calc td -props*/
      SxString tdInterpolation;
      /** \brief rotor-type (for anharmonic treatment of the methyl side-chain
           in the alanine peptide-unit*/
      SxString rotorType;
      /** \brief type of empirical van der Waals correction (ElstnerI,II or
           WTYang I,II, s.a. SxVDW*/
      SxString vdwCorrection;

      /**\brief hessian matrix: is loaded in and then manipulated according
                to user input*/
      SxMatrix<Double> hessian;
      /**\brief basis hessian matrix: in case that the input hessian matrix
        has been calculated in a refinement calculation, this matrix contains
        the according phonon basis */
      SxMatrix<Double> basisHessian;
      /**\brief contains atomiv species in a list*/
      SxList<SxString> speciesNameList;
      /**\brief contains atomic labels in crystal-order, is needed for helices*/
      SxList<SxVector<Double> > periodicity;
      /**\brief contains id's of eigenfrequencies to be replaced
                and replacements*/
      SxList<SxList<double> > replaceList;
      /**\brief atomic structure*/
      SxAtomicStructure tau;
      /**\brief ionic masses*/
      SxAtomicStructure masses;
      /**\brief control flags*/
      bool printAnharmonic, printThermodynamics, printEigenvalues, printMolden,
      overwriteFreezeMode, overwriteVibRots, applyPeriodicity, applyAveraging,
      printHessian, printDispersion, fromDispersion, printCouplings, fileBasis,
      printForceConstants, printNuMu, incorporateRotor, oldCalc,
      lowCoop, isDynamical, replaceCurvature, grueneisen;
      /**\brief all hessian-entries with abosulte values below this
                treshhold are set to zero*/
      double setZero;
      /**\brief pressure (for thermodynamic analysis, feature to be
                transfered from rel-10)*/
      double   p;
      /**\brief initial ans final temperature for thermodyn. analysis*/
      double startT, endT;
      /**\brief z-lattice constant*/
      double zExt;
      /**\brief rotational barrier of methyl-group in case of poly-alanine*/
      double rotorBarrier;
      /**\brief control values*/
      int sigma, replicZ, ahId, nSteps, nDoF, vibRots,
      res, supercellCutoff, tempSamplePoints, branchToReplace, nTurns;

      /**\brief object for treating secondary structure*/
      SxSecondaryStructure peptideChain;
      /**\brief object for treating hessians, dynamical matrices*/
      SxHessianOps hessianOps;
      /**\brief s.a. SxElemDB*/
      SxElemDB elemDB;

      /**\brief Symbol Table for force constant output*/
      SxSymbolTable *FCGroup;
      /**\brief Symbol Table for grueneisen parameters*/
      SxSymbolTable *grueneisenGroup;

      void parseInputFile();
      void printParameters ();
      static void flushTemplate (const SxString &filename);
      void validateInput ();
      void initSystem ();
      void manipulateInput ();
      void printOutput ();
      void printAnharmonicCorrections ();
      void printPhononDispersion ();
      void printPhononCouplings ();
      void printFreqVsMass ();

      SxAtomicStructure loadStructureFHI98 ();
      void getMasses ();
      void setElementsToZero (SxMatrix<Double> *matrix, double);
      void symmetrizeHessian (SxMatrix<Double> *matrix);
   protected:
};


#endif /* _SXVIBRATIONS_H_ */
