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

#ifndef _SX_EXX_SOLVER_H_
#define _SX_EXX_SOLVER_H_
#include <SxHamSolver.h>
#include <SxExx.h>

/** \brief EXX solver

    \b SxClass = S/PHI/nX EXX solver class


    \author M. Wahn
    Reintegrated by C. Freysoldt, freysoldt@mpie.de
*/
class SX_EXPORT_EXX SxExxSolver : public SxHamSolver
{
   public:
      /** \brief Constructor */
      SxExxSolver (const SxAtomicStructure &structureIn,
                   const SxSymbolTable *table);
      /** \brief Check if a command is registered */
      virtual bool isRegistered (const SxSymbolTable *cmd) const;
      /** \brief Check if command is an EXX solver command */
      static bool isExx (const SxSymbolTable *cmd);
      /** \brief Execute specified command */
      virtual void execute (const SxSymbolTable *cmd, bool calc=true);
   public:
      void writeR111Plot (const SxString &filename,
                          const SxDiracVec<Double> &vecR);
      void writeGPlot (const SxString &filename,
                       const SxDiracVec<Double> &vecR);
      void writeGPlot (const SxString &filename,
                       const SxDiracVec<Complex16> &vecG);

   public:

      /** \brief Electronic minimisation using orbital dependend
                 functionals, such as KLI, Slater potential. */
      void odp (const SxSymbolTable *, bool calc=true);  

      /** \brief standard exact exchange formalism */
      void stEXX (const SxSymbolTable *, bool calc=true);

      //---------------------------------------------------------------------
      /**@name Perturbation methods */
      //---------------------------------------------------------------------
      //@{
      /** \brief ground state exact exchange formalism */
      void gsEXX (const SxSymbolTable *, bool calc=true);
      //@}

      /** Solves the Sternheimer equation
          \f[   (H^{0} - \varepsilon_{i}^{0}) | \psi_{i}^{1} \rangle
              + (H^{1} - \varepsilon_{i}^{1}) | \psi_{i}^{0} \rangle
              = 0\,
          \f] using the steepest descent scheme.
     
          \author Matthias Wahn, wahn@fhi-berlin.mpg.de
       */
      PsiG sternheimer (const PsiG &h1psi0, const PsiG &psi1Start);

      PsiG sternheimerStDesc (const PsiG &h1psi0, const PsiG &psi1Start);

      PsiG sternheimerStDesc_write (const PsiG &h1psi0, const PsiG &psi1Start);

      PsiG sternheimerConjGrad (const PsiG &h1psi0, const PsiG &psi1Start);

      PsiG sternheimerConjGrad_write (const PsiG &h1psi0,
                                      const PsiG &psi1Start);

      PsiG sternheimerPrecCG (const PsiG &h1psi0, const PsiG &psi1Start);

      PsiG sternheimerPrecCG_write (const PsiG &h1psi0, const PsiG &psi1Start);

      int minSternSteps, maxSternSteps;
      PrecEnergy epsResidue;

      enum SternMethod  {
         NoMethod = 0,
         StDesc   = 1,
         ConjGrad = 2,
         PrecCG   = 3
      };

      SternMethod  sternMethod;
      int  precondType;
      SxPW waves1, waves1amp;

      double checkSternheimer (const PsiG &psi1,
                               const PsiG &h1psi0,
                               bool        dump=true);

      double sternFunctional (const PsiG &psi1, const PsiG &h1psi0);

      /// Apply conduction band Green function to h1psi0 (which can be anything)
      PsiG greensfunction (const PsiG &h1psi0);

      void updateEXXCPotential (const SxSymbolTable *, int step);

      /** \brief Compute EXX potential for fixed exchange operator (linear mix)

        This computes the induced density and updates the EXX potential
        to make that zero. That's similar to a SCF density mixing.

        This routine uses only linear mixing.
        */
      void relaxEXXPotentialLin (const SxSymbolTable *, int step);

      /** \brief Compute EXX potential for fixed exchange operator (linear mix)

        This computes the induced density and updates the EXX potential
        to make that zero. That's similar to a SCF density mixing.

        This routine uses the mixer from the input file.
        */
      void relaxEXXPotential (const SxSymbolTable *, int step);

      SxMeshG computeRhoInducedVogl (const SxMeshG &vExG,
                                     bool writeConvergence=false);

      SxMeshG computeRhoInduced (const SxMeshG &vExG,
                                 bool writeConvergence=false);

      bool            useGreensfunction;
      // true: do not symmetrize chi
      bool            mimicVoglChi;
      // true: do not symmetrize 
      // (delta E_x / delta psi * delta psi/delta v_x +cc)
      // = Green(i) | V_fock |psi_i>
      bool            useAbdullahSymmetry; 
      // true: do not symmetrize V_Fock
      bool            useNonsymFock;
      // perfect: mimicVoglChi=false, useAbdullahSymmetry=false, useNonSymFock=true
      // acceptable: everything true
      bool            fockStorage;
      bool            checkSternheimerEq;

      SxMeshR rhoOld, rhoDiff;

};

namespace Timer {
   enum ExxSolverTimer {
      ExxSolver,
      ExxStartup,
      RelaxEXX,
      OrbitalDependentPotential,
      Sternheimer,
      GreensFunction
   };
}

SX_REGISTER_TIMERS (Timer::ExxSolverTimer)
{
   using namespace Timer;
   regTimer (ExxSolver,                 "EXX solver");
   regTimer (ExxStartup,                "EXX startup");
   regTimer (RelaxEXX,                  "EXX pot. relaxation");
   regTimer (OrbitalDependentPotential, "odp solver");
   regTimer (Sternheimer,               "Sternheimer");
   regTimer (GreensFunction,            "Green's function");
}

#endif /* _SX_EXX_SOLVER_H_ */
