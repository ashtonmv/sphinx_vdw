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

#include <SxStructOpt.h>
#include <SxProjPairMatrix.h>
#include <SxPDBFast.h>
#include <SxStickyFilter.h>
#include <SxDriftFilter.h>
#include <SxLineFilter.h>
#include <SxHessian.h>
#include <SxTimer.h>
#include <SxSimpleParser.h>
#include <SxProcess.h>
#include <SxHamSolver.h>
#include <SxRedundantCoords.h>
#include <SxFileIO.h>

//---------------------------------------------------------------------------
// Bibliography:
//   ref1 - P. v. Rague Schleyer, 
//          Enc. of Comp. Chemistry, 1, pp.1136-1156
//   ref2 - fh96md paper


SxStructOpt::SxStructOpt () : potentialPtr(NULL)
{
   // empty
}


SxStructOpt::SxStructOpt (SxPotential *P, SxAtomicStructure &str)
   : potentialPtr(NULL)
{
   set (P, str);
}

SxStructOpt::~SxStructOpt ()
{
   // empty
}


//---------------------------------------------------------------------------
// Controlling class by the Input file parser
//---------------------------------------------------------------------------
bool SxStructOpt::isRegistered (const SxSymbolTable *cmd)
{
   SX_CHECK (cmd);
   SxString str = cmd->getName();
   SX_MPI_MASTER_ONLY if (str == "extControl") {
      SxString ctlFile = SxProcess::getenv ("SX_EXT_CTRL");
      SxString resFile = SxProcess::getenv ("SX_EXT_RES");
      bool stop = false;
      if (ctlFile.getSize () == 0) {
         cout << "extControl needs control file name via environment SX_EXT_CTRL" << endl;
         stop = true;
      }
      if (resFile.getSize () == 0) {
         cout << "extControl needs results file name via environment SX_EXT_RES" << endl;
         stop = true;
      }
      if (stop) { SX_QUIT; }
   }
   return (str == "DN" || str == "QN" || str == "linQN" || str == "ricQN"
           || str == "extControl");
}

void SxStructOpt::print (const SxSymbolTable *cmd)
{
   execute (cmd, false);
}

void SxStructOpt::execute (const SxSymbolTable *cmd, bool calc)
{
   SX_CHECK (cmd);
   SX_CHECK (potentialPtr);

   SxString str = cmd->getName();

   // --- clean up files produced by SxStructOpt
   SX_OUTPUT_FILE ("energy-structOpt.dat");
   SX_OUTPUT_FILE ("relaxHist.sx");

   // --- initialize electronic minimizer commands
   elMinimCmds = potentialPtr->getMinimCmds (cmd);
   if (elMinimCmds.getSize() == 0 && !potentialPtr->isRegistered (NULL))  {
      sxprintf ("No valid command found in group 'bornOppenheimer'\n");
      SX_QUIT;
   }

   if      (str == "DN")     dampedNewton       (cmd, calc);
   else if (str == "QN")     quasiNewton        (cmd, calc);
   else if (str == "linQN")  linearQuasiNewton  (cmd, calc);
   else if (str == "ricQN")  ricQuasiNewton     (cmd, calc);
   else if (str == "extControl")  extControl  (cmd, calc);
   else  {
      cout << "ERROR: Unknown command " << str << endl;
      SX_QUIT;
   }
}




//---------------------------------------------------------------------------
// Structure optimization routines
//---------------------------------------------------------------------------
void SxStructOpt::dampedNewton (const SxSymbolTable *cmd, bool calc)
{  
   //TODO revised
   SX_EXIT;
   
   SX_CHECK (potentialPtr);

   maxSteps = cmd->get("maxSteps")->toInt();

   cout << SX_SEPARATOR;
   cout << "| Damped Newton\n";
   cout << SX_SEPARATOR;
   cout << "|  max. steps:    " << maxSteps << endl;
   cout << SX_SEPARATOR;
   if (!calc)  return;


   // --- entities depending on :iSpecies and :iAtom
   SxAtomicStructure tau    = structure;
   SxAtomicStructure f      = getForces(tau);
   SxAtomicStructure tauOld;

   // --- entities depending on :iSpecies only
   SxVector<Double> lambda = species.dampingMass;
   SxVector<Double> mu     = species.reciprocalMass;

   // --- optimization loop
   for (int it=0; it < maxSteps; it++)  {
      f = /*F **/ getForces(tau);
      //tau = (1. + lambda)*tau - lambda*tauOld + mu*f;     // ref2 (??)
      tau = lambda*tau;
      tauOld.copy (tau);

      cout << SX_SEPARATOR;
      cout << "it=" << it+1 << endl;
   }

   // --- synchronize structure and tau;
   //structure.set (tau); TODO
}

/** \brief Filter for atomic structure optimization steps

    This class filters cartesian steps (in degree-of-freedom space) by
    wrapping SxAtomicStructure-based filter-chain.

    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_STRUCT SxStepFilter : public SxOperatorBase<SxVector<Double> >
{
   protected:
      SxConstPtr<SxAtomInfo> atomInfo;
      SxCell cell;
   public:
      /// Constructor
      SxStepFilter (const SxAtomicStructure &str)
         : atomInfo(str.atomInfo),
           cell(str.cell)
      { /* empty */ }

      /// Empty destructor
      virtual ~SxStepFilter () {/* empty */}

      /// The actual filter chain (based on SxAtomicStructure)
      SxOperator<SxAtomicStructure> F;

      /// Apply the center-of-mass-filter out-of-place
      virtual SxVector<Double> operator*(const SxVector<Double> &step) const
      {
         SxAtomicStructure stepAsStr(step);
         stepAsStr.cell = cell;
         stepAsStr.replaceInfo (atomInfo);
         SxVector<Double> res = (F | stepAsStr).coordRef ().getCopy ();
         return res;
      }

      /// Apply the center-of-mass-filter in-place
      virtual void applyInPlace (SxVector<Double> &step) const
      {
         step = operator* (step);
      }

      /** \brief Get SxPtr copy
          \note Standard copy constructor is OK, no dangerous members
      */
      SXOPERATOR_GETCOPY(SxStepFilter,SxAtomicStructure);

};



SxOperator<SxAtomicStructure>
SxStructOpt::getFilter (const SxSymbolTable *table)
{
   SxOperator<SxAtomicStructure> F;
   SxPtr<SxStepFilter> stepF = SxPtr<SxStepFilter>::create (structure);
   double stepMayDrift = false;
   // --- use drift filter unless disabled
   bool useDriftFilter = table->contains ("driftFilter")
                       ? table->get("driftFilter")->toAttribute ()
                       : potentialPtr->noNetForce;
   if (useDriftFilter)  {
      F << SxPtr<SxDriftFilter>::create ();
   } else {
      cout << "Disabling drift filter" << endl;
      stepMayDrift = true;
   }
   table = table->topLevel();

   if (table->containsGroup ("structure"))  {
      table = table->getGroup ("structure");
      SxPtr<SxStickyFilter> stickyFilter = SxPtr<SxStickyFilter>::create(table);

      // --- validate the sticky filter
      SX_CHECK (structure.cell.symGroupPtr);
      const SxSymGroup &S = *structure.cell.symGroupPtr;
      SxVector<Int> equivalentIdx (structure.getNAtoms(), 0);
      int           nSym = S.getSize ();
      cout << "| Validating sticky filter ....\n\n";
      for (int iSym = 0; iSym < nSym; iSym++) {
         structure.isEqual (S(iSym) ^ structure, &equivalentIdx);
         stickyFilter->validate (S(iSym).rot, equivalentIdx);
      }
      structure.atomInfo->meta.attach (SxAtomicStructure::StickyFilter,
                                       stickyFilter->getStickyArray ());

      if (!stickyFilter->isUnity ())  {
         F << stickyFilter;
         stepF->F << stickyFilter;
         stepMayDrift = true;
      }

      SxPtr<SxLineFilter> lineFilter
         = SxPtr<SxLineFilter>::create(table, structure);

      if (!lineFilter->isEmpty ())  {
         F << lineFilter;
         stepF->F << lineFilter;
         stepMayDrift = true;
      }
   }

   if (!stepMayDrift)
      stepF->F << SxDriftFilter ();
   // set up the step filter
   stepFilter = stepF;

   if (F.isEmpty ()) F = SxOperator<SxAtomicStructure>::Identity ();
   return F;
}

namespace {
   inline bool cmdMatch (const char *cmd, const char *pattern)  {
      // run until end of either string, or character mismatch
      while (*cmd == *pattern && (*pattern) && (*cmd)) { cmd++; pattern++ ;}
      return (*pattern == 0); // pattern matches if end of pattern is reached
   }
   template<class T>
   inline const char* name(const SxPtr<T> &)  { return typeid(T).name (); }
   template<class T>
   void requires(const SxPtr<T> &ptr, FILE *extRes)
   {
      if (!ptr)  {
         cout << "No " << name(ptr) << " available!" << endl;
         fprintf (extRes, "ERROR");
         SX_QUIT;
      }
   }
// define template specialization for name(...)
#define NAME_PTR(Type,Name) \
   template <> inline const char* name(const SxPtr<Type> &) { return Name ; }\
   class AbsorbingTheTrailingSemicolon

   NAME_PTR(SxPAWHamiltonian, "PAW Hamiltonian");
   NAME_PTR(SxHubbardU, "Hubbard U");
   NAME_PTR(SxSpinConstraint, "spin constraint");

   // disguised void to silence warning about unused result from fscanf etc.
   class Void { public: template <class T> Void(const T&) {/* empty */} };
}


void SxStructOpt::extControl (const SxSymbolTable *cmd, bool calc)
{
   enum {
      EndControl            = -1,
      GetForces             = -2,
      EnableSpinConstraint  = -3,
      DisableSpinConstraint = -4,
      RecomputeHamiltonian  = -5,
      ReadWaves             = -6,
      ReadRho               = -7,
      WriteWaves            = -8, // not yet implemented
      WriteRho              = -9  // not yet implemented
   };
   SX_CHECK (potentialPtr);
   SxList<SxArray<const SxSymbolTable *> > minimizers;
   SxList<SxString> minName;

   cout << SX_SEPARATOR;
   cout << "| External control\n";
   cout << SX_SEPARATOR;
   bool noForces = true;
   SYMBOLPARSE(cmd)  {
      ssize_t minimId = 1;
      FOREACH_SYMBOLGROUP("bornOppenheimer")  {
         SxString id = SYMBOLGET("id") || (SxString() + minimId);
         minName << id;
         minimizers << potentialPtr->getMinimCmds (SYMBOLGROUP_TABLE);
         cout << "| offering minimizer id = " << id
               << " (" << minimizers.last ().getSize () << "-step)" << endl;
         minimId++;
      }
      noForces = SYMBOLGET("noForces").toBool ();
   }

   cout << SX_SEPARATOR;
   if (!calc) return;

   // get rid of abstraction: need hamSolver and (maybe) PAW Hamiltonian
   SxHamSolver* hamSolver = dynamic_cast<SxHamSolver*>(potentialPtr);
   SX_CHECK (hamSolver);
   SxPtr<SxPAWHamiltonian> pawHam;
   if (dynamic_cast<SxPAWHamiltonian*>(hamSolver->hamPtr.getPtr ()))  {
      pawHam = hamSolver->hamPtr;
   }

   // --- handle slave MPI tasks (the master counterpiece is done in "run")
   SX_MPI_SLAVE_ONLY {
      int minimId = EndControl;
      while(true /* exit is inside */)
      {
         // --- get the minimizer id from master
         minimId = SxLoopMPI::bcast (minimId, 0 /* master */);
         if (minimId == EndControl) break;

         if (minimId == EnableSpinConstraint)  {
            hamSolver->spinConstraint = SxPtr<SxSpinConstraint>::create (hamSolver->potPtr);
            hamSolver->spinConstraint->targetSpin.resize (structure.getNAtoms ());
            hamSolver->spinConstraint->constrainedAtom.resize (structure.getNAtoms ());
            pawHam->nuA.resize (structure.getNAtoms ());
            pawHam->nuA.set (0.);
            continue;
         }
         if (minimId == DisableSpinConstraint)  {
            hamSolver->spinConstraint = SxPtr<SxSpinConstraint> ();
            pawHam->nuA.resize (0);
            continue;
         }
         if (minimId == ReadWaves)  {
            continue;
            int size = SxLoopMPI::bcast(-1, 0);
            SX_CHECK (size > 0, size);
            SxString fileName;
            fileName.resize (size);
            SxLoopMPI::bcast(const_cast<char*>(fileName.ascii ()), size, 0);
            try {
               SxBinIO io(fileName, SxBinIO::BINARY_READ_ONLY);
               hamSolver->wavesPtr->read (io, SxPWSet::KeepAll);
               hamSolver->fermi.read (io);
               io.close ();
            } catch (SxException e)  {
               e.print ();
               SX_EXIT;
            }
         }
         if (minimId == ReadRho)  {
            // reading is done by MPI master
            hamSolver->hamPtr->getRho ().syncMPI ();
            continue;
         }


         // --- broadcast everything that can be set in extControl
         // structure
         SxLoopMPI::bcast (structure, 0);
         // spin constraints
         if (hamSolver->spinConstraint)  {
            SxLoopMPI::bcast (hamSolver->spinConstraint->targetSpin, 0);
            SxLoopMPI::bcast (hamSolver->spinConstraint->constrainedAtom, 0);
         }
         if (pawHam && pawHam->hubbardU)  {
            SxArray<double> &alpha = pawHam->hubbardU->alphaBias;
            alpha.resize (pawHam->hubbardU->getSize ());
            SxLoopMPI::bcast (alpha.elements, alpha.getSize (), 0);
         }

         // --- run the electronic loop
         if (minimId == GetForces)  {
            // special case: "get forces" before "run"
            getForces ();
         } else if (minimId == RecomputeHamiltonian)  {
            hamSolver->fermi.fermiDistribution (hamSolver->ekt);
            pawHam->compute(*hamSolver->wavesPtr, hamSolver->fermi);
         } else if (noForces)  {
            SX_LOOP(iMinim) hamSolver->execute (minimizers(minimId)(iMinim));
         } else {
            potentialPtr->getSymForces (structure, minimizers(minimId));
         }
      }
      return;
   }

   // --- open control and results file
   SxString ctlFile = SxProcess::getenv ("SX_EXT_CTRL");
   SxString resFile = SxProcess::getenv ("SX_EXT_RES");
   (cout << "Opening control file '" << ctlFile << "'...").flush ();
   FILE *extCtl = fopen (ctlFile.ascii (), "r");
   cout << endl;
   if (!extCtl) {
      cout << "Failed to open control file '" << ctlFile
           << "' for reading.";
   }
   (cout << "Opening result file '" << resFile << "'...").flush ();
   FILE *extRes = fopen (resFile.ascii (), "w");
   cout << endl;
   if (!extRes) {
      cout << "Failed to open result file '" << resFile
           << "' for writing.";
   }
   if (! (extRes && extCtl))  {
      if (extCtl) fclose (extCtl);
      if (extRes) fclose (extRes);
      SX_QUIT;
   }

   // --- main control loop
   enum { Crash, Stop, Ignore } onProblem = Stop;
   bool confirm = false;
   SxAtomicStructure forces;
   cout << SX_SEPARATOR;
   while (!feof(extCtl))  {

      (cout << "| extControl: external command=").flush ();
      // ignore white space
      (Void)fscanf (extCtl, " ");

      // read next command
      char extCmd[1024];
      bool okCtl = fgets(extCmd, 1024, extCtl);
      if (!okCtl) {
         cout << endl;
         cout << "ERROR: end of file in '" << ctlFile << "'." << endl;
         fprintf(extRes, "ERROR: reached end of file in control file '%s'.\n"
                 "Send 'end' for a proper shutdown of extControl.\n",
                 ctlFile.ascii ());
         fflush (extRes);
         if (onProblem == Crash) { SX_QUIT; }
         return;
      }
      for (char *p = extCmd; *p != 0; p++)  {
         // chomp the newline
         if (*p == '\n') { *p = 0; break; }
         // lowercase
         if (*p >= 'A' && *p <= 'Z') *p = char(*p  - 'A' + 'a');
      }

      cout << extCmd << endl;
      cout.flush ();
      // --- end of commands
      if (cmdMatch(extCmd, "end")) break;
      // --- confirmation
      else if (cmdMatch(extCmd, "confirm"))  {
         for (int ic = 8; ic < 1024; ic++)  {
            if (extCmd[ic] != ' ')  {
               confirm = (extCmd[ic] != 'n' && extCmd[ic] != 'N');
               break;
            }
         }
      }
      // --- total energy
      else if (cmdMatch(extCmd, "get energy"))  {
         double E = potentialPtr->getPotentialEnergy ();
         fprintf (extRes, "%.16f\n", E);
      }
      // --- atomic forces
      else if (cmdMatch(extCmd, "get forces"))  {
         if (forces.getNAtoms () == 0)  {
            if (noForces) {
               cout << "ERROR: forces have been disabled in input file." << endl;
               fprintf (extRes, "ERROR\n");
               fflush (extRes);
               if (onProblem == Crash) { SX_EXIT; }
               if (onProblem == Stop) { SX_QUIT; }
               break;
            }

            // --- send the minimizer id to slaves
            SxLoopMPI::bcast (GetForces, 0 /* master */);
            // --- broadcast everything that can be set in extControl
            // structure
            SxLoopMPI::bcast (structure, 0);
            // spin constraints
            if (hamSolver->spinConstraint)  {
               SxLoopMPI::bcast (hamSolver->spinConstraint->targetSpin, 0);
            }
            // alpha for Hubbard U
            if (pawHam && pawHam->hubbardU)  {
               SxArray<double> &alpha = pawHam->hubbardU->alphaBias;
               if (alpha.getSize () == 0)  {
                  alpha.resize (pawHam->hubbardU->getSize ());
                  alpha.set (0.);
               }
               SxLoopMPI::bcast (alpha.elements, alpha.getSize (), 0);
            }
            forces = getForces ();
         }
         SX_LOOP(ia)  {
            fprintf (extRes, "%.16f %.16f %.16f\n",
                     forces(ia)(0), forces(ia)(1), forces(ia)(2));
         }
      }
      // --- get number of atoms
      else if (cmdMatch(extCmd, "get natoms"))  {
         fprintf(extRes, "%d\n", structure.getNAtoms ());
      }
      // --- get structure
      else if (cmdMatch(extCmd, "get structure"))  {
         SX_LOOP(ia)  {
            fprintf(extRes, "%.8f %.8f %.8f\n", structure.ref(ia)(0),
                                                structure.ref(ia)(1),
                                                structure.ref(ia)(2));
         }
      }
      // --- get structure in Angstrom (atomic structure algorithm protocol)
      else if (cmdMatch(extCmd, "get positions"))  {
         const SxArray<SxString> &elem = structure.getElements ();
         SX_LOOP2(is, ia)  {
            fprintf(extRes, "%.8f %.8f %.8f %2s\n",
                    structure.ref(is,ia)(0),
                    structure.ref(is,ia)(1),
                    structure.ref(is,ia)(2),
                    elem(is).ascii ());
         }
      }
      // --- get cell in Angstrom (atomic structure algorithm protocol)
      else if (cmdMatch(extCmd, "get cell"))  {
         for (int a = 0; a < 3; a++)  {
            fprintf(extRes, "%.8f %.8f %.8f\n",
                        structure.cell(0,a),
                        structure.cell(1,a),
                        structure.cell(2,a));
         }
      }
      // --- set structure
      else if (   cmdMatch(extCmd, "set structure")
               || cmdMatch(extCmd, "set positions"))
      {
         if (noForces) {
            cout << "ERROR: forces have been disabled in input file." << endl;
            cout << "ERROR: structure setting is not possible without forces"
                 << endl;
            fprintf (extRes, "ERROR\n");
            fflush (extRes);
            if (onProblem == Crash) { SX_EXIT; }
            if (onProblem == Stop) { SX_QUIT; }
            break;
         }
         SX_LOOP(ia)  {
            // read new structure
            int ok = fscanf(extCtl, "%lf %lf %lf",
                            &structure.ref(ia)(0),
                            &structure.ref(ia)(1),
                            &structure.ref(ia)(2));
            if (ok != 3)  {
               cout << "Failed to coordinates for atom no. "
                    << (ia + 1) << endl;
               fprintf (extRes, "ERROR\n");
               SX_QUIT;
            }
         }
      }
      else if (cmdMatch(extCmd, "shift atom"))  {
         int ia = -1;
         if (sscanf(extCmd + 10, "%d", &ia) != 1
             || (ia < 1) || (ia > structure.getNAtoms ()))  {
            cout << "ERROR: atom id is invalid" << endl;
            fprintf (extRes, "ERROR\n");
            fflush (extRes);
            if (onProblem == Crash) { SX_EXIT; }
            if (onProblem == Stop) { SX_QUIT; }
         } else {
            ia--;
            Coord xyz;
            if (fscanf(extCtl, "%lf %lf %lf", &xyz(0), &xyz(1), &xyz(2)) != 3) {
               cout << "ERROR: Failed to read shift" << endl;
               fprintf (extRes, "ERROR\n");
               fflush (extRes);
               if (onProblem == Crash) { SX_EXIT; }
               if (onProblem == Stop) { SX_QUIT; }
            }
            structure.ref(ia) += xyz;
         }
      }
      // --- unsupported commands in atomic structure algorithm protocol
      else if (   cmdMatch(extCmd, "set cell")
               || cmdMatch(extCmd, "get stress"))  {
         if (confirm) fprintf(extRes, "no\n");
         else {
            cout << "Unsupported command " << extCmd << endl;
            fprintf (extRes, "ERROR\n");
            fflush (extRes);
            if (onProblem == Crash) { SX_EXIT; }
            if (onProblem == Stop)  { SX_QUIT; }
         }
      }
      // --- get number of species
      else if (cmdMatch(extCmd, "get nspecies"))  {
         if (confirm) fprintf(extRes, "ok\n");
         fprintf (extRes, "%u\n", structure.getNSpecies ());
         for (int is = 0; is < structure.getNSpecies (); ++is)
            fprintf (extRes, "%u\n", structure.getNAtoms (is));
      }
      // --- get elements
      else if (cmdMatch(extCmd, "get elements"))  {
         if (confirm) fprintf(extRes, "ok\n");
         const SxArray<SxString> &elem = structure.getElements ();
         SX_LOOP(is) fprintf (extRes, "%s\n", elem(is).ascii ());
      }
      // --- get number of spin constraint atoms
      else if (cmdMatch(extCmd, "get nspinconstraints"))  {
         requires (pawHam, extRes);
         fprintf(extRes, "%d\n", (int)pawHam->nuA.getSize ());
      }
      // --- spin constraint forces
      else if (cmdMatch(extCmd, "get nu"))  {
         requires (pawHam, extRes);
         for(ssize_t iNu = 0; iNu < pawHam->nuA.getSize (); iNu++)
            fprintf(extRes, "%.16f\n", pawHam->nuA(iNu));
      }
      // --- set spin constraint
      else if (cmdMatch(extCmd, "set spinconstraint"))  {
         requires(hamSolver->spinConstraint, extRes);
         SxVector<Double> &spins =  hamSolver->spinConstraint->targetSpin;
         SX_CHECK (spins.getSize () > 0);
         SX_LOOP(iSpin)  {
            int ok = fscanf(extCtl, " ");
            int c = fgetc (extCtl);
            if (c == 'x' || c == 'X')  {
               hamSolver->spinConstraint->constrainedAtom(iSpin) = false;
               spins(iSpin) = 0;
               continue;
            } else {
               hamSolver->spinConstraint->constrainedAtom(iSpin) = true;
               ungetc(c, extCtl);
            }
            // read spin constraint
            ok = fscanf(extCtl, "%lf", & spins(iSpin));
            if (ok != 1)  {
               cout << "Failed to read spin for constrained atom no. "
                    << (iSpin+1) << endl;
               fprintf (extRes, "ERROR\n");
               SX_QUIT;
            }
         }
         hamSolver->spinConstraint->checkSym (structure);
      }
      // --- get spin constraint
      else if (cmdMatch(extCmd, "get spinconstraint"))  {
         requires(hamSolver->spinConstraint, extRes);
         SX_CHECK (hamSolver->spinConstraint->targetSpin.getSize () > 0);
         SX_LOOP (ia) {
            if (hamSolver->spinConstraint->constrainedAtom(ia))
               fprintf(extRes, "%.16f\n",
                       hamSolver->spinConstraint->targetSpin(ia));
            else
               fprintf(extRes, "X\n");
         }
      }
      // --- get atomic spins
      else if (cmdMatch(extCmd, "get atomspin"))  {
         if (hamSolver->wavesPtr->getNSpin () != 2)  {
            cout << "No spin in calculation!" << endl;
            fprintf (extRes, "ERROR");
            SX_QUIT;
         }
         SxArray<double> spins = pawHam->pawRho.getSpinMom (structure);
         SX_LOOP (ia) fprintf(extRes, "%.16f\n", spins(ia));
      }
      // --- enable spin constraints
      else if (cmdMatch(extCmd, "enable spinconstraint"))  {
         if (hamSolver->wavesPtr->getNSpin () != 2)  {
            cout << "No spin in calculation!" << endl;
            fprintf (extRes, "ERROR");
            SX_QUIT;
         }
         hamSolver->spinConstraint = SxPtr<SxSpinConstraint>::create (hamSolver->potPtr);
         hamSolver->spinConstraint->targetSpin = pawHam->pawRho.getSpinMom (structure);
         hamSolver->spinConstraint->constrainedAtom.resize (structure.getNAtoms ());
         hamSolver->spinConstraint->constrainedAtom.set (true);
         pawHam->nuA.resize (structure.getNAtoms ());
         pawHam->nuA.set (0.);
         SxLoopMPI::bcast (EnableSpinConstraint, 0);
      }
      // --- disable spin constraints
      else if (cmdMatch(extCmd, "disable spinconstraint"))  {
         hamSolver->spinConstraint = SxPtr<SxSpinConstraint> ();
         pawHam->nuA.resize (0);
         SxLoopMPI::bcast (DisableSpinConstraint, 0);
      }
      // --- get number of hubbard sites
      else if (cmdMatch(extCmd, "get hubbard nsite"))  {
         requires(pawHam, extRes);
         requires(pawHam->hubbardU, extRes);
         fprintf (extRes, "%ld\n", pawHam->hubbardU->getSize ());
      }
      // --- get non-Hubbard-U on-site shifts 
      else if (cmdMatch(extCmd, "get hubbard alpha"))  {
         requires(pawHam, extRes);
         requires(pawHam->hubbardU, extRes);
         if (pawHam->hubbardU->alphaBias.getSize () == 0)  {
            pawHam->hubbardU->alphaBias.resize (pawHam->hubbardU->getSize ());
            pawHam->hubbardU->alphaBias.set (0.);
         }
         SX_LOOP(iSite)
            fprintf(extRes, "%.16f\n", pawHam->hubbardU->alphaBias(iSite));
      }
      // --- get hubbard site occupations
      else if (cmdMatch(extCmd, "get hubbard occ"))  {
         requires(pawHam, extRes);
         requires(pawHam->hubbardU, extRes);
         SX_LOOP(iSite)
            fprintf(extRes, "%.16f\n", pawHam->hubbardU->siteOccupation(iSite));
      }
      // --- recompute Hamiltonian
      else if (cmdMatch(extCmd, "compute hamiltonian"))  {
         requires(pawHam, extRes);
         // --- send the minimizer id to slaves
         SxLoopMPI::bcast (RecomputeHamiltonian, 0 /* master */);
         // --- broadcast everything that can be set in extControl
         // structure
         SxLoopMPI::bcast (structure, 0);
         // spin constraints
         if (hamSolver->spinConstraint)  {
            SxLoopMPI::bcast (hamSolver->spinConstraint->targetSpin, 0);
            SxLoopMPI::bcast (hamSolver->spinConstraint->constrainedAtom, 0);
         }
         // alpha for Hubbard U
         if (pawHam->hubbardU)  {
            SxArray<double> &alpha = pawHam->hubbardU->alphaBias;
            if (alpha.getSize () == 0)  {
               alpha.resize (pawHam->hubbardU->getSize ());
               alpha.set (0.);
            }
            SxLoopMPI::bcast (alpha.elements, alpha.getSize (), 0);
         }
         hamSolver->fermi.fermiDistribution (hamSolver->ekt);
         pawHam->compute(*hamSolver->wavesPtr, hamSolver->fermi);
      }
      // --- set non-Hubbard-U on-site shifts 
      else if (cmdMatch(extCmd, "set hubbard alpha"))  {
         requires(pawHam, extRes);
         requires(pawHam->hubbardU, extRes);
         if (pawHam->hubbardU->alphaBias.getSize () == 0)  {
            pawHam->hubbardU->alphaBias.resize (pawHam->hubbardU->getSize ());
            pawHam->hubbardU->alphaBias.set (0.);
         }
         SX_LOOP(iSite)  {
            int ok = fscanf(extCtl, "%lf", &pawHam->hubbardU->alphaBias(iSite));
            if (!ok)  {
               cout << "Failed to read alpha bias for site " << (iSite+1)
                    << endl;
               fprintf (extRes, "ERROR");
               SX_QUIT;
            }
         }
      }
      // --- read waves
      else if (cmdMatch(extCmd, "read waves"))  {
         SxString fileName(extCmd + 10);
         fileName = fileName.trim ();
         if (fileName == "-")  {
            char buffer[4096];
            buffer[0] = 0;
            fgets (buffer, 4096, extCtl);
            fileName = SxString(buffer);
            if (fileName(fileName.getSize () - 1) == '\n')
               fileName(fileName.getSize () -1 ) = ' ';
            fileName = fileName.trim ();
         }
         if (fileName.getSize () == 0) fileName="waves.sxb";
         cout << "Reading waves and fermi from '" << fileName << "'." << endl;
         SxLoopMPI::bcast (ReadWaves, 0);
         int size = (int)fileName.getSize ();
         SxLoopMPI::bcast (size, 0);
         SxLoopMPI::bcast (const_cast<char*>(fileName.ascii ()),
                           fileName.getSize (), 0);
         try {
            SxBinIO io(fileName, SxBinIO::BINARY_READ_ONLY);
            hamSolver->wavesPtr->read (io, SxPWSet::KeepAll);
            hamSolver->fermi.read (io);
            io.close ();
         } catch (SxException e)  {
            e.print ();
            SX_EXIT;
         }
      }
      // --- read density
      else if (cmdMatch(extCmd, "read rho"))  {
         SxString fileName(extCmd + 8);
         fileName = fileName.trim ();
         if (fileName.getSize () == 0) fileName="rho.sxb";
         cout << "Reading rho from '" << fileName << "'." << endl;
         SxLoopMPI::bcast (ReadRho, 0);
         hamSolver->hamPtr->getRho ().readRho (fileName);
         hamSolver->hamPtr->getRho ().syncMPI ();
      }
      // --- run electronic loop
      else if (cmdMatch(extCmd, "run"))  {
         const char* id = extCmd + 3;
         while (*id == ' ' || *id == '\t') id++;
         int minimId = 0;
         if (*id == 0)  {
            cout << "Standard minimizer" << endl;
         } else {
            minimId = (int)minName.findPos (id);
            if (minimId < 0)  {
               cout << "ERROR: Unknown bornOppenheimer '" << id << "'." << endl;
               fprintf (extRes, "ERROR\n");
               fflush (extRes);
               if (onProblem == Crash) { SX_EXIT; }
               if (onProblem == Stop) { SX_QUIT; }
               continue;
            }
            cout << "Minimizer " << id <<  endl;
         }

         // --- send the minimizer id to slaves
         minimId = SxLoopMPI::bcast (minimId, 0 /* master */);
         // --- broadcast everything that can be set in extControl
         // structure
         SxLoopMPI::bcast (structure, 0);
         // spin constraints
         if (hamSolver->spinConstraint)  {
            SxLoopMPI::bcast (hamSolver->spinConstraint->targetSpin, 0);
            SxLoopMPI::bcast (hamSolver->spinConstraint->constrainedAtom, 0);
         }
         // alpha for Hubbard U
         if (pawHam && pawHam->hubbardU)  {
            SxArray<double> &alpha = pawHam->hubbardU->alphaBias;
            if (alpha.getSize () == 0)  {
               alpha.resize (pawHam->hubbardU->getSize ());
               alpha.set (0.);
            }
            SxLoopMPI::bcast (alpha.elements, alpha.getSize (), 0);
         }

         // --- run the electronic loop
         cout << SX_SEPARATOR;
         if (noForces)  {
            SX_LOOP(iMinim) hamSolver->execute (minimizers(minimId)(iMinim));
         } else {
            forces = potentialPtr->getSymForces (structure, minimizers(minimId));
         }
         cout << SX_SEPARATOR;
      }
      else if (cmdMatch(extCmd, "onproblem=crash"))  { onProblem = Crash; }
      else if (cmdMatch(extCmd, "onproblem=stop"))  { onProblem = Stop; }
      else if (cmdMatch(extCmd, "onproblem=ignore")) { onProblem = Ignore; }
      else {
         cout << "Unknown command" << endl;
         fprintf (extRes, "ERROR\n");
         fflush (extRes);
         if (onProblem == Crash)  {
            SX_EXIT;
         }
      }
      fflush (extRes);
   }
   SX_MPI_MASTER_ONLY {
      // send "end" to slave tasks
      SxLoopMPI::bcast (EndControl, 0 /* master */);
   }
   fclose (extRes);
   fclose (extCtl);

}

// khr: set up a timer to measure the total time spent in the QN loop
namespace Timer
{
   enum qnTimer { qnLoop };
}
SX_REGISTER_TIMERS (Timer::qnTimer)
{
   using namespace Timer;
   regTimer (qnLoop, "QN optimization");
}


void SxStructOpt::quasiNewton (const SxSymbolTable *cmd, bool calc)
{
   maxSteps = cmd->contains ("maxSteps") ? cmd->get("maxSteps")->toInt() : 50;
   double dX    = cmd->contains("dX")    ? cmd->get("dX")->toReal()
                                         : 1e-2;
   double dF    = cmd->contains("dF")    ? cmd->get("dF")->toReal()
                                         : 1e-3;
   double dXavg = cmd->contains("dXavg") ? cmd->get("dXavg")->toReal ()
                                         : dX;
   double dFavg = cmd->contains("dFavg") ? cmd->get("dFavg")->toReal ()
                                         : dF;
   bool   saveRelaxHist = cmd->contains("saveRelaxHist") ? 
                          cmd->get("saveRelaxHist")->toAttribute () : false;

   double dEnergy = cmd->contains("dEnergy")
                  ? cmd->get("dEnergy")->toReal () : 1e-4;
   double maxStepLength = cmd->contains("maxStepLength")
                        ? cmd->get("maxStepLength")->toReal ()
                        : 0.3;
   
   double E  = 0.;
   bool thirdOrderCor = false;
   SxSecondaryStructure peptideChain;
   bool secondaryStructure = false;
   
   if (cmd -> contains ("thirdOrderCor"))
      thirdOrderCor = cmd->get ("thirdOrderCor")->toAttribute ();

   // --- user-defined: if the system is a helix, helical symmetry
   //     is enforced during symmetrisation
   
   if (cmd -> containsGroup("secondaryStructure")) {
      SxSymbolTable *helixGroup = cmd -> getGroup("secondaryStructure");
      SxString motif
         (helixGroup -> get("motif") -> toString ());
      peptideChain = SxSecondaryStructure (structure, species,  
            helixGroup -> get("periodicity") -> toString (),
            motif,
            helixGroup -> get("nTurns") -> toInt (), true);
      secondaryStructure = true;
      structure = peptideChain.getCoords();
   }
   
   SxVector3<Double> comForce;
   
   cout << SX_SEPARATOR;
   cout << "| Quasi Newton (BFGS scheme)\n";
   cout << SX_SEPARATOR;
   cout << "|  max. steps:    " << maxSteps << endl;
   cout << "|  dX (max)  :    " << dX << endl;
   cout << "|  dX (avg)  :    " << dXavg << endl;
   cout << "|  dF (max)  :    " << dF << endl;
   cout << "|  dF (avg)  :    " << dFavg << endl;
   cout << "|  dEnergy   :    " << dEnergy << endl;
   cout << SX_SEPARATOR;

   // --- setup force filter
   SxOperator<SxAtomicStructure> F = getFilter (cmd);

   if (!calc)  return;
   FILE *relaxHistFile = NULL;
   SX_MPI_MASTER_ONLY
   {
      relaxHistFile = sxfopen("relaxHist.sx", "a");
   }

   // --- initialize forces and Hessian
   SxAtomicStructure g, gOld, trialG;
   SxVector<TPrecTauR> s, y, xOld, l,  x = structure.coordRef();
   SxVector<TPrecTauR> trialX, s2, g2;
   SxMatrix<TPrecTauR> B;
   if (cmd->containsGroup ("hessian"))  {
      SxSymbolTable *hessianGroup = cmd->getGroup ("hessian");
      try {
         SxString fileName = hessianGroup->get ("file")->toString ();
         cout << "Reading initial Hessian from " << fileName << endl;
         SxBinIO io(fileName, SxBinIO::BINARY_READ_ONLY);
         B = SxHessian::readFull (io);
         cout << B.eigenvalues () << endl;
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
   } else if (cmd->containsGroup ("ric"))  {
      SxRedundantCoords ric;
      bool withAngle = false;
      bool SchlegelParam = false, FischerParam = false;
      double bondForceDefault = 0., angleForceDefault = 0.;
      SYMBOLPARSE(cmd)  {
         SYMBOLGROUP("ric")  {
            ric.maxDist         = SYMBOLGET("maxDist")           || 10.;
            ric.primarySetLimit = SYMBOLGET("typifyThreshold")   || 0.05;
            ric.rmsThreshold    = SYMBOLGET("rmsThreshold") || 3.;
            ric.planeCutLimit   = SYMBOLGET("planeCutLimit")     || 0.95;
            ric.verbose = true;
            withAngle = SYMBOLGET("withAngles").toBool ();
            SchlegelParam = SYMBOLGET("Schlegel").toBool ();
            FischerParam  = SYMBOLGET("Fischer").toBool ();
            bondForceDefault  = SYMBOLGET("bondConstant") || 0.02;
            angleForceDefault = SYMBOLGET("angleConstant") || 0.1;
         }
      }
      ric.setup (structure);
      if (withAngle) ric.getAngles (structure);
      ric.param.resize (ric.getNParam ());
      if (SchlegelParam)  {
         ric.paramSchlegel (structure.getElements ());
         cout << "param(Schlegel)=" << ric.param << endl;
      } else if (FischerParam) {
         ric.paramFischer (structure.getElements ());
         cout << "param(Fischer)=" << ric.param << endl;
      } else {
         ric.param.set (bondForceDefault);
         if (withAngle)  {
            ssize_t nBond = ric.getNBond ();
            for (ssize_t i = nBond; i < ric.getNParam (); i++)  {
               ric.param(i) = angleForceDefault;
            }
         }
      }
      ssize_t nDof = 3 * structure.getNAtoms ();
      B.reformat (nDof, nDof);
      SxVector<Double> one(nDof);
      for (int i = 0; i < nDof; ++i)  {
         one.set (0.);
         one(i) = -1.;
         B.colRef(i) <<= (ric.applyH (structure, one));
      }
      B *= 0.5;
      cout << (B - B.transpose ()).normSqr () << endl;
      //SX_LOOP2(i,j) { if (i != j) B(i,j) = 0.; }
      B += B.transpose ();
      SxMatrix<Complex16>::Eigensystem eig
         = SxMatrix<Complex16>(B).eigensystem ();
      cout << eig.vals.real () << endl;
      SX_LOOP(i)  {
         if (eig.vals(i).re < 1e-6)  {
            SX_LOOP2(j,k) B(j,k) += (1. - eig.vals(i).re)
                                  * (eig.vecs(j,i) * eig.vecs(k,i).conj ()).re;
         }
      }
      eig = SxMatrix<Complex16>(B).eigensystem ();
      cout << eig.vals.real () << endl;
   } else {
      cout << "Initialising Hessian with 1" << endl;
      B = SxMatrix<TPrecTauR> (nDoF,nDoF).identity ();
   }
   // --- variables for third order correction
   double /*toCorrect, */deltaE, a, b, c, sol1, sol2, curv1, curv2;
   double sMax, gMax, sAvg, gAvg;
   int is, ia;
   FILE *fp;

   {
      SxAtomicStructure force = getForces (x);
      g = F | -force;
      if (relaxHistFile)  {
         fprintf (relaxHistFile, "// --- step 0\n");
         structure.fprint (relaxHistFile, force);
         fflush (relaxHistFile);
      }
   }

   //cout << "forces :\n" << g << endl;
   // calculate convergence parameters
   g2 = g.absSqr ();
   gMax = sqrt(g2.maxval ());
   gAvg = sqrt(g2.sum () / (double)g2.getSize ());


   E = potentialPtr->getPotentialEnergy ();

   sxprintf ("QN[0]: max(|f|)=%12.8f, max(|dx|)=%12.8f, E=%12.8f H\n", 
                gMax, 0., E);
   sxprintf ("rms(dx)=%12.8f, rms(f)=%12.8f\n", 0., gAvg);
   

   SX_MPI_MASTER_ONLY
   {
      SxPDBFast ("tau_ini.pdb").write (structure, species.chemName);
      if (saveRelaxHist)  {
         fp = fopen("relaxHist.dat", "w");
         if (fp)  {
            fprintf (fp, "# Position x y z    Force x y z\n");
            fprintf (fp, "# Iteration 0\n");
            for (is = 0; is < structure.getNSpecies (); is++)  {
               for (ia = 0; ia < structure.getNAtoms (is); ia++)  {
                  fprintf (fp, "% f  % f  % f\t\t % f  % f  % f\n",
                        structure(is,ia)(0), structure(is,ia)(1),
                        structure(is,ia)(2), -g(is,ia)(0),
                        -g(is,ia)(1), -g(is,ia)(2));
               }
            }
            fclose(fp);
         }
      }
   }


   SX_START_TIMER(Timer::qnLoop);
   // --- optimization loop  (khr: implicitly parallel by means of LoopMPI)
   for (int it=0; it < maxSteps; it++)  {

      cout << SX_SEPARATOR;
      //cout << "it=" << it+1 << endl;

      xOld.copy (x);
      gOld = g;
      {
         SxVector<Double> step;
         step = stepFilter | (B.inverse() ^ g.coordRef ()); // ref1 (28)
         double step2Max = 0.;
         for (int ja = 0; ja < structure.getNAtoms (); ++ja)  {
            double step2Atom = sqr(step(3*ja))
                             + sqr(step(3*ja + 1))
                             + sqr(step(3*ja + 2));
            if (step2Atom > step2Max) step2Max = step2Atom;
         }
         if (step2Max > sqr(maxStepLength))  {
            double scale = maxStepLength / sqrt(step2Max);
            cout << "limit step length: scaling step by " << scale << endl;
            step *= scale;
         }
         x -= step;
      }
      //cout << "X(" << it << ")\n" << x << endl;

      // --- print structure after every structOpt step
      SX_MPI_MASTER_ONLY
      {
         fp = fopen("relaxedStr.sx", "w");
         if (fp)  { structure.fprint(fp); fclose(fp); }
      }
      
      if (secondaryStructure) {
         SxSymbolTable *helixGroup = cmd -> getGroup("secondaryStructure");
         SxString motif
             (helixGroup -> get("motif") -> toString ());
         peptideChain = SxSecondaryStructure (structure, species,  
               helixGroup -> get("periodicity") -> toString (),
               motif,
               helixGroup -> get("nTurns") -> toInt (), true);
         structure = peptideChain.getCoords();
         x = SxVector3<Double> (); // forget old reference
         x = structure.coordRef ();
      }

      //cout << "STRUCTURE\n" << structure << endl;

      {
         SxAtomicStructure force = getForces (x);
         g  = F | -force;
         if (relaxHistFile)  {
            fprintf (relaxHistFile, "// --- step %d\n", it + 1);
            structure.fprint (relaxHistFile, force);
            fflush (relaxHistFile);
         }
      }

      deltaE = potentialPtr->getPotentialEnergy () - E;

      E  = potentialPtr->getPotentialEnergy ();


      //--- Third Order Correction:
      //    Third Order Polynmom f(x) = a*x^3 + b*x^2 + c*x 
      //    is fitted along the actual search direction
      if (thirdOrderCor) {
         //--- Coefficients a,b and c are evaluated from projected gradients and energies
         //    along the actual search direction
         l = B.inverse () ^ gOld.coordRef ();
         // f(x) = a*x^4 + b*x^3 + c*x^2 + d
         //toCorrect = sqrt(l.absSqr ().sum ());
         c = -(gOld.coordRef ()^l).chop ();
         b = 3.*deltaE - (2.*c - (g.coordRef ()^l).chop ());
         a = -2.*deltaE + c - (g.coordRef ()^l).chop ();

         if (b*b > 3.*a*c) {
         
            //--- minimum of the polynomial fit is evaluated
            sol1 = (-b - sqrt(b*b - 3.*a*c))/3./a;
            sol2 = (-b + sqrt(b*b - 3.*a*c))/3./a;
            curv1 = 6.*a*sol1 + 2.*b;
            curv2 = 6.*a*sol2 + 2.*b;
           

            if (curv1 > 0) trialX = xOld - sol1*l;
            else {
               if (curv2 > 0) trialX = xOld - sol2*l;
            }
            
            //--- on the evaluated minum forces and energies are measured 
            //   (additional call of the electronic loop necessary)
            SxVector<Double> currentX;
            currentX.copy (x); // note: x is reference to structure and
            trialG = F | -getForces (trialX);

            //--- if the energy on the evaluated is lower then the energy from
            //    the harmonic guess, forces, structure and energy are updated
            //    accordingly
            if (potentialPtr->getPotentialEnergy () < E) {
               g = trialG;
               E = potentialPtr->getPotentialEnergy ();
            } else {
               x <<= currentX;

               // --- print structure after every structOpt step
               SX_MPI_MASTER_ONLY
               {
                  fp = fopen("relaxedStr.sx", "w");
                  if (fp)  { structure.fprint(fp); fclose(fp); }
               }
            }
         } 
      } 

      s  = x - xOld;                             // ref1 (30), 
      y  = g.coordRef () - gOld.coordRef ();     // ref1 (31)

      // calculate convergence parameters
      s2 = SxAtomicStructure(s,SxAtomicStructure::Reference).absSqr ();
      g2 = g.absSqr ();
      sMax = sqrt(s2.maxval ());
      gMax = sqrt(g2.maxval ());
      sAvg = sqrt(s2.sum () / (double)s2.getSize ());
      gAvg = sqrt(g2.sum () / (double)g2.getSize ());

      {
         if (s.normSqr () > 1e-10)
            cout << "| Residual force along optimization line: "
                 << (dot(s,g.coordRef ()) / s.norm ())
                 << " from " << (dot(s,gOld.coordRef ()) / s.norm ()) << " ("
                 << (100. * dot(s,g.coordRef ()) / dot(s, gOld.coordRef ()))
                 << "%)" << endl;
         double dEEstimate = -0.5 * dot(s,B ^ s);
         double dEHarmonic = 0.5 * dot(s,(g + gOld).coordRef ());
         cout << "| estimated energy gain: " << dEEstimate << endl;
         cout << "| obtained energy gain:  " << deltaE << endl;
         cout << "| difference:            " << (deltaE - dEEstimate) << endl;
         cout << "| harmonic energy gain:  " << dEHarmonic << endl;
         cout << "| difference:            " << (deltaE - dEHarmonic) << endl;
      }
      
      // --- print out, and to file 
      sxprintf ("QN[%d]: max(|f|)=%12.8f, max(|dx|)=%12.8f, E=%12.8f H\n", 
               it+1, gMax, sMax, E);
      sxprintf ("rms(dx)=%12.8f, rms(f)=%12.8f\n", sAvg, gAvg);

      SX_MPI_MASTER_ONLY
      {

         SxFileIO::appendToFile
         (  SxString(it+1)                  + "\t"  // iteration
               + SxString(E,"%15.12f")           + "\n"  // E_tot  [H]
         , "energy-structOpt.dat");

         // --- print structure to pdb file
         SxPDBFast ("tau_tmp.pdb").write (structure, species.chemName);

         if (saveRelaxHist)  {
            fp = fopen("relaxHist.dat", "a");
            if (fp)  {
               fprintf (fp, "# Iteration %d\n", it+1);
               for (is = 0; is < structure.getNSpecies (); is++)  {
                  for (ia = 0; ia < structure.getNAtoms (is); ia++)  {
                     fprintf (fp, "% f  % f  % f\t\t % f  % f  % f\n",
                           structure(is,ia)(0), structure(is,ia)(1),
                           structure(is,ia)(2), -g(is,ia)(0),
                           -g(is,ia)(1), -g(is,ia)(2));
                  }
               }
               fclose(fp);
            }
         }
      } // LoopMPI

      // --- ref1 (34)
      double ys = dot(y,s);
      cout << "y ^ s = " << ys << endl;
      //if (ys > 0.) {
         B -= (B^s^s.transpose()^B.transpose()) / (s.transpose()^B^s).chop() 
            - (y^y.transpose())                 / (y.transpose()^s).chop();
      //}


      SX_MPI_MASTER_ONLY
      {
         SxHessian::write ("hessian.sxb", B, structure);
      }

      // --- convergence check
      if (    sMax         < dX      // maximum structural step
           && gMax         < dF      // maximum force
           && sAvg         < dXavg   // square-average structural step
           && gAvg         < dFavg   // square-average force
           && fabs(deltaE) < dEnergy )  {
         sxprintf ("Convergence reached\n"); break;
      }
      
   }
   SX_STOP_TIMER(Timer::qnLoop);
   
   // --- print relaxed structure into a file
   SX_MPI_MASTER_ONLY
   {
      fp = fopen("relaxedStr.sx", "w");
      if (fp)  { structure.fprint(fp); fclose(fp); }
      // --- print structure to pdb file
      SxPDBFast ("tau_end.pdb").write (structure, species.chemName);
   }

   if (relaxHistFile) fclose (relaxHistFile);
}



void SxStructOpt::linearQuasiNewton (const SxSymbolTable *cmd, bool calc)
{
   double dX = 0., dF = 0., dEnergy = 0., diag = 0., maxStepLength = 0.;
   int nProj = 0;
   try {
      maxSteps = cmd->contains ("maxSteps") ? cmd->get("maxSteps")->toInt(): 50;
      dX    = cmd->contains("dX") ? cmd->get("dX")->toReal() : 1e-2;
      dF    = cmd->contains("dF") ? cmd->get("dF")->toReal() : 1e-3;
      nProj = cmd->contains("nProjectors")
            ? cmd->get("nProjectors")->toInt ()
            : 10;
      dEnergy = cmd->contains("dEnergy")
              ? cmd->get("dEnergy")->toReal ()
              : 1e-4;
      diag = cmd->contains("diag") ? cmd->get("diag")->toReal () : 1.;
      maxStepLength = cmd->contains("maxStepLength")
                    ? cmd->get("maxStepLength")->toReal ()
                    : 0.3;
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }

   cout << SX_SEPARATOR;
   cout << "| Quasi Newton (linear scaling BFGS scheme)\n";
   cout << SX_SEPARATOR;
   cout << "|  max. steps                 :    " << maxSteps << endl;
   cout << "|  dX (max)                   :    " << dX << endl;
   cout << "|  dF (max)                   :    " << dF << endl;
   cout << "|  dEnergy                    :    " << dEnergy << endl;
   cout << "|  projector pairs            :    " << nProj << endl;
   cout << "|  diagonal of inverse Hessian:    " << diag << endl;
   cout << "|  max. step length           :    " << maxStepLength << endl;
   cout << SX_SEPARATOR;

   if (!calc)  return;
   FILE *relaxHistFile = NULL;
   SX_MPI_MASTER_ONLY
   {
      relaxHistFile = fopen("relaxHist.sx", "a");
      if (!relaxHistFile)  {
         cout << "WARNING: Failed to open relaxHist.sx for writing" << endl;
      }
   }

   // --- setup force filter
   SxOperator<SxAtomicStructure> F = getFilter (cmd);

   // --- variables for convergence
   SxVector<TPrecTauR> s2, g2;
   double sMax, gMax, sAvg, gAvg, E, deltaE;

   // --- initialize forces and inverse Hessian
   SxVector<TPrecTauR> x = structure.coordRef(), s, y;
   SxAtomicStructure g, gOld;

   {
      SxAtomicStructure force = getForces ();
      g  = F | -force;
      if (relaxHistFile)  {
         fprintf (relaxHistFile, "// --- step 0\n");
         structure.fprint (relaxHistFile, force);
         fflush (relaxHistFile);
      }
   }

   g2 = g.absSqr ();
   gMax = sqrt(g2.maxval ());
   gAvg = sqrt(g2.sum () / (double)g2.getSize ());

   E = potentialPtr->getPotentialEnergy ();
   sxprintf ("linQN[0]: max(|f|)=%12.8f, max(|dx|)=%12.8f, E=%12.8f H\n",
           gMax, 0., E);
   sxprintf ("rms(f)=%12.8f, rms(dx)=%12.8f\n", gAvg, 0.);

   SX_MPI_MASTER_ONLY
   {
      SxFileIO::appendToFile
      (  SxString(0)                     + "\t"  // iteration
       + SxString(E,"%15.12f")           + "\n"  // E_tot  [H]
      , "energy-structOpt.dat");
   }

   SxList<SxVector<TPrecTauR> > sk,yk;
   SxProjPairMatrix invB(diag, nProj); // inverse Hessian

   // --- optimization loop
   for (int it=0; it < maxSteps; it++)  {
      cout << SX_SEPARATOR;

      s = - (invB ^ g.coordRef ()); // ref1 (30) and (28)

      // calculate convergence parameters
      s2 = SxAtomicStructure(s,SxAtomicStructure::Reference).absSqr ();
      sMax = sqrt(s2.maxval ());
      sAvg = sqrt(s2.sum () / (double)s2.getSize ());

      if (sMax > maxStepLength)  {
         // rescale step
         double scale =  maxStepLength / sMax;
         cout << "limit step length: scaling step by " << scale << endl;
         s    *= scale;
         sMax *= scale;
         sAvg *= scale;
      }

      // --- get new structure and forces
      x += s;                      // ref1 (28) and (30)
      //cout << "X(" << it << ")\n" << x << endl;

      // --- print structure after every structOpt step
      SX_MPI_MASTER_ONLY {
         FILE *fp = fopen("relaxedStr.sx", "w");
         if (fp)  { structure.fprint(fp); fclose(fp); }
      }

      gOld = g;
      {
         SxAtomicStructure force = getForces ();
         g  = F | -force;
         if (relaxHistFile)  {
            fprintf (relaxHistFile, "// --- step %d\n", it + 1);
            structure.fprint (relaxHistFile, force);
            fflush (relaxHistFile);
         }
      }

      //cout << "Gradient:" << g << endl;

      cout << SX_SEPARATOR;
      if (s.normSqr () > 1e-10)
         cout << "| Residual force along optimization line: "
              << (dot(s,g.coordRef ()) / s.norm ())
              << " from " << (dot(s,gOld.coordRef ()) / s.norm ()) << " ("
              << (100. * dot(s,g.coordRef ()) / dot(s, gOld.coordRef ()))
              << "%)" << endl;

      y  = g.coordRef () - gOld.coordRef ();  // ref1 (31)

      deltaE = potentialPtr->getPotentialEnergy () - E;
      E  = potentialPtr->getPotentialEnergy ();
      {
         double dEEstimate = 0.5 * dot(s,gOld.coordRef ());
         double dEHarmonic = 0.5 * dot(s,(g + gOld).coordRef ());
         cout << "| estimated energy gain: " << dEEstimate << endl;
         cout << "| obtained energy gain:  " << deltaE << endl;
         cout << "| difference:            " << (deltaE - dEEstimate) << endl;
         cout << "| harmonic energy gain:  " << dEHarmonic << endl;
         cout << "| difference:            " << (deltaE - dEHarmonic) << endl;
      }
      cout << SX_SEPARATOR;

      // calculate convergence parameters
      g2 = g.absSqr ();
      gMax = sqrt(g2.maxval ());
      gAvg = sqrt(g2.sum () / (double)g2.getSize ());


      // --- print out, and to file
      sxprintf ("linQN[%d]: max(|f|)=%12.8f, max(|dx|)=%12.8f, E=%12.8f H\n",
               it+1, gMax, sMax, E);
      sxprintf ("rms(f)=%12.8f, rms(dx)=%12.8f\n", gAvg, sAvg);

      SX_MPI_MASTER_ONLY {
         SxFileIO::appendToFile
         (  SxString(it+1)                  + "\t"  // iteration
          + SxString(E,"%15.12f")           + "\n"  // E_tot  [H]
         , "energy-structOpt.dat");
      }


      // --- convergence check
      if (    sMax         < dX      // maximum structural step
           && gMax         < dF      // maximum force
           && fabs(deltaE) < dEnergy )  {
         sxprintf ("Convergence reached\n"); break;
      }

      // --- ref1 (34)
      {
         double ys = dot(y,s);
         //cout << "y ^ s = " << ys << endl;
         if (ys <= 0.) {
            // update doesn't preserve positive definiteness
            cout << "Suppressing update of approximate inverse Hessian ..."
                 << endl;
         } else {
            SxVector<TPrecTauR> t = invB ^ y;
            double yt = dot(y,t);
            yk.append (y);
            sk.append (s);

            invB.append (s, t, 1./ys,
                         yt/ys+1., -1.,
                         -1.,      0.);

            // --- check for and remove out-of-date data
            bool reinit = false;
            // --- does current B still predict s within reasonable accuracy?
            for (int k = 0; k < yk.getSize (); ++k)  {
               double relErr = ((invB ^ yk(k))-sk(k)).norm () / sk(k).norm ();
               //cout << "k=" << k << ": " << relErr << endl;
               if (relErr > 1.)  {
                  sk.remove(k);
                  yk.remove(k);
                  --k;
                  reinit=true;
               }
            }
            // --- number of projectors exceeded?
            if (sk.getSize () > nProj)  {
               sk.removeFirst ();
               yk.removeFirst ();
               reinit = true;
            }

            if (reinit)  {
               // --- reinitialize B
               invB = SxProjPairMatrix(diag, nProj);
               SxList<SxVector<TPrecTauR> >::Iterator sIt, yIt;
               sIt = sk.begin ();
               yIt = yk.begin ();
               for (; yIt != yk.end (); ++yIt, ++sIt)  {
                  t = invB ^ (*yIt);
                  yt = dot(*yIt,t);
                  ys = dot (*yIt, *sIt);
                  invB.append (*sIt, t, 1./ys,
                               yt/ys+1., -1.,
                               -1.,      0.);
               }
            }

            /*
            {
               int nDof = s.getSize ();
               SxMatrix<Double> B(nDof, nDof);
               SxVector<Double> one(nDof);
               for (int i = 0; i < nDof; ++i)  {
                  one.set (0.);
                  one(i) = 1.;
                  B.colRef(i) <<= (invB ^ one);
               }
               cout << '*' << B << endl;
               cout << B.eigensystem ().vals.real () << endl;
            }
            */
         }
      }
   }

   // --- print relaxed structure into a file
   SX_MPI_MASTER_ONLY {
      FILE *fp = fopen("relaxedStr.sx", "w");
      if (fp)  { structure.fprint (fp); fclose(fp); }
   }
   if (relaxHistFile) fclose (relaxHistFile);
}

enum RicTimer { RicUpdate, RicParam, RicApply, RicBtest };

SX_REGISTER_TIMERS(RicTimer)
{
   regTimer (RicUpdate, "ric update");
   regTimer (RicParam, "ric param");
   regTimer (RicApply, "ric apply");
   regTimer (RicBtest, "ric B test");
}

class SxRicPrecond
{
   public:
      SxVector<Double> operator^ (const SxVector<Double> &f) const;

      SxRedundantCoords ric;

      SxAtomicStructure &structure;

      double accuracy;
      double diagRic;
      double softModeDamping;

      SxRicPrecond (SxAtomicStructure &str)
         : structure(str),
           diagRic(1.)
      {
         accuracy = 1e-8;
         softModeDamping=1e-2;
      }

      /// Setup ric
      void setup (double bondForceDefault, bool withAngle,
                  double angleForceDefault,
                  const SxArray<int> &bvkAtoms);
      void reparamRic (const SxList<SxAtomicStructure> &configs,
                       const SxList<SxAtomicStructure> &forces);
      SxMatrix<Double> getHcart () const;

      double maxDistFromSetup (const SxAtomicStructure &str) const;
   protected:
      SxAtomicStructure setupStr;
};


void SxRicPrecond::setup (double bondForceDefault, bool withAngle,
                       double angleForceDefault,
                       const SxArray<int> &bvkAtoms)
{
   ric.setup (structure);
   if (withAngle) ric.getAngles (structure);
   if (bvkAtoms.getSize () > 0)
      ric.getBornVonKarmanAngles (bvkAtoms);
   ric.param.resize (ric.getNParam ());
   ric.param.set (bondForceDefault);
   if (withAngle)  {
      ssize_t nBond = ric.getNBond ();
      for (ssize_t i = nBond; i < ric.getNParam (); i++)  {
         ric.param(i) = angleForceDefault;
      }
   }
   setupStr.copy (structure);
}

double SxRicPrecond::maxDistFromSetup (const SxAtomicStructure &str) const
{
   double m2 = 0.;
   SX_CHECK (str.getNAtoms () == setupStr.getNAtoms (),
             str.getNAtoms (), setupStr.getNAtoms ());
   for (int ia = 0; ia < str.getNAtoms (); ++ia)  {
      double d2 = (str.constRef(ia) - setupStr.constRef(ia)).normSqr ();
      if (d2 > m2) m2 = d2;
   }
   return sqrt(m2);
}

SxVector<Double> SxRicPrecond::operator^ (const SxVector<Double> &f) const
{
   SX_CLOCK(RicApply);
   SxVector<Double> xRic = ric.solve (structure, -f, softModeDamping, accuracy);
   //cout << "|xRic|=" << xRic.norm () << endl;
   SxVector<Double> fRest = f + ric.applyH(structure, xRic);
   return xRic + diagRic * fRest;
}

void SxRicPrecond::reparamRic (const SxList<SxAtomicStructure> &configs,
                 const SxList<SxAtomicStructure> &forces)
{
   SX_CLOCK(RicParam);
   if (configs.getSize () < 2) return;
   int nBond = ric.getNBond ();
   SxMatrix<Double> A(nBond, nBond);
   SxVector<Double> rhs(nBond);
   A.set (0.);
   rhs.set (0.);
   for (int iConfig = 1; iConfig < configs.getSize (); iConfig++)  {
      SxVector<Double> dx = configs(iConfig  ).coordRef ()
                            - configs(iConfig-1).coordRef ();
      SxAtomicStructure avgStruct = 0.5*(configs(iConfig) + configs(iConfig-1));
      SxVector<Double> B = ric.getParamDeriv (avgStruct, dx);

      SxVector<Double> Bt = B.transpose ();
      A += (1./dx.normSqr ()) * (Bt ^ B);
      rhs += Bt ^ (  forces(iConfig  ).coordRef ()
                   - forces(iConfig-1).coordRef ()) / dx.normSqr ();
   }
   SxVector<Double> paramNew(nBond);
   int nNeg = 0, dofRicParam;
   do {
      SxMatrix<Double>::Eigensystem eig = A.eigensystem ();
      paramNew.set (0.);
      dofRicParam = 0;
      SX_LOOP (i)  {
         SxVector<Double> vec = eig.vecs.colRef(i);
         double val = eig.vals(i).re;
         if (fabs (val) > 1e-10)  {
            // inverse of A
            paramNew.plus_assign_ax (dot(vec, rhs)/val, vec);
            dofRicParam++;
         } else {
            // projected old parameter
            paramNew.plus_assign_ax (dot(vec, ric.param(SxIdx(0,nBond-1))), vec);
         }
      }
      bool anyNeg = false;
      SX_LOOP(ip)  {
         if (paramNew(ip) < -1e-12) {
            //double paramChangeTo = 0.2 * ric.param(ip);
            double paramChangeTo = 0.;
            if (paramChangeTo < 0.) paramChangeTo = 0;
            cout << "Changing param " << (ip+1) << " from "
                 << paramNew(ip) << " to " << paramChangeTo << endl;
            rhs(ip) = paramChangeTo;
            SX_LOOP(jp) A(ip,jp) = A(jp,ip) = 0.;
            A(ip,ip) = 1.;
            paramNew(ip) = paramChangeTo;
            anyNeg = true;
            nNeg++;
            break;
         }
      }
      if (anyNeg) {
         cout << paramNew << endl;
         cout << "Refitting..." << endl;
      } else {
         break;
      }
   } while (nNeg < nBond);
   ric.param(SxIdx(0, nBond-1)) <<= paramNew;
   if (ric.getNParam () > nBond)  {
      SxIdx idxAng(nBond, (int)ric.param.getSize () - 1);
      SxVector<Double> paramOld = ric.param(idxAng).getCopy ();
      ric.param(idxAng).set (0.);
      int nAng = ric.getNParam () - nBond;
      A.reformat (nAng, nAng);
      A.set (0.);
      rhs.resize (nAng);
      rhs.set (0.);
      for (int iConfig = 1; iConfig < configs.getSize (); iConfig++)  {
         SxVector<Double> dx = configs(iConfig  ).coordRef ()
                               - configs(iConfig-1).coordRef ();
         SxAtomicStructure avgStruct = 0.5*(configs(iConfig) + configs(iConfig-1));
         SxVector<Double> B = ric.getParamDerivA (avgStruct, dx);

         SxVector<Double> Bt = B.transpose ();
         A += (1./dx.normSqr ()) * (Bt ^ B);
         SxVector<Double> df = forces(iConfig  ).coordRef ()
                             - forces(iConfig-1).coordRef ();
         //cout << "df=" << df << endl;
         cout << "|df|=" << df.norm ();
         df -= ric.applyH(avgStruct, dx); // subtract bond forces
         //cout << "df2=" << df << endl;
         cout << " |df2|=" << df.norm () << endl;
         rhs += (1./dx.normSqr ()) * (Bt ^ df);
      }
      paramNew.resize (nAng);
      int dofRicParamA;
      int nNegA = 0;
      do {
         SxMatrix<Double>::Eigensystem eig = A.eigensystem ();
         paramNew.set (0.);
         dofRicParamA = 0;
         SX_LOOP (i)  {
            SxVector<Double> vec = eig.vecs.colRef(i);
            double val = eig.vals(i).re;
            if (fabs (val) > 1e-10)  {
               // inverse of A
               paramNew.plus_assign_ax (dot(vec, rhs)/val, vec);
               dofRicParamA++;
            } else {
               // projected old parameter
               paramNew.plus_assign_ax (dot(vec, paramOld), vec);
               cout << eig.vals(i) << endl;
               cout << eig.vecs.colRef (i) << endl;
            }
         }
         bool anyNeg = false;
         SX_LOOP(ip)  {
            if (paramNew(ip) < -1e-12) {
               cout << "Changing param " << (nBond + ip+1) << " from "
                    << paramNew(ip) << " to " << 0 << endl;
               rhs(ip) = 0.;
               SX_LOOP(jp) A(ip,jp) = A(jp,ip) = 0.;
               A(ip,ip) = 1.;
               paramNew(ip) = 0.;
               anyNeg = true;
               nNegA++;
               break;
            }
         }
         if (anyNeg) {
            cout << paramNew << endl;
            cout << "Refitting..." << endl;
         } else {
            break;
         }
      } while (nNegA < nAng);
      dofRicParam += dofRicParamA;
      ric.param(idxAng) <<= paramNew;
      nNeg += nNegA;
   }
   cout << "ric parameters optimized: " << (dofRicParam - nNeg);
   if (nNeg > 0) cout << '+' << nNeg << "x(=0)";
   cout << "/" << ric.getNParam () << endl;
   cout << "ric param: " << ric.param << endl;
   double optShift = 0.;
   for (int iConfig = 1; iConfig < configs.getSize (); iConfig++)  {
      SxVector<Double> dx = configs(iConfig  ).coordRef ()
                            - configs(iConfig-1).coordRef ();
      SxAtomicStructure avgStruct = 0.5*(configs(iConfig) + configs(iConfig-1));
      SxVector<Double> df = forces(iConfig  ).coordRef ()
                          - forces(iConfig-1).coordRef ();
      //cout << "df=" << df << endl;
      cout << "|df|=" << df.norm ();
      df -= ric.applyH(avgStruct, dx); // subtract bond forces
      //cout << "df2=" << df << endl;
      cout << " |df2|=" << df.norm () << endl;
      optShift += dot (dx, df) / dx.normSqr ();
   }
   optShift /= -double (configs.getSize () - 1);
   cout << "optimum shift=" << optShift << endl;
   if (optShift > 0.) softModeDamping = optShift;
}

SxMatrix<Double> SxRicPrecond::getHcart () const
{
   int nDof = structure.getNAtoms () * 3;
   SxMatrix<Double> H(nDof, nDof);
   SxVector<Double> one(nDof);
   for (int i = 0; i < nDof; ++i)  {
      one.set (0.);
      one(i) = -1.;
      H.colRef(i) <<= (ric.applyH (structure, one));
   }
   H *= 0.5;
   cout << (H - H.transpose ()).normSqr () << endl;
   H += H.transpose ();
   SxMatrix<Complex16>::Eigensystem eig
      = SxMatrix<Complex16>(H).eigensystem ();
   cout << eig.vals.real () << endl;
   SX_LOOP(i)  {
      if (eig.vals(i).re < 1e-6)  {
         SX_LOOP2(j,k) H(j,k) += (1. - eig.vals(i).re)
                               * (eig.vecs(j,i) * eig.vecs(k,i).conj ()).re;
      }
   }
   return H;
}

void SxStructOpt::ricQuasiNewton (const SxSymbolTable *cmd, bool calc)
{
   double dX = 0., dF = 0., dEnergy = 0., maxStepLength = 0.;
   int nProj = 0;
   double bondForceDefault  = 0.05;
   double angleForceDefault = 0.1;
   SxRicPrecond invBric(structure); // inverse Hessian
   invBric.ric.verbose = false;
   bool withAngle = false;
   SxArray<int> bvkAtoms;
   SYMBOLPARSE(cmd)  {
      maxSteps      = SYMBOLGET("maxSteps") || 50;
      dX            = SYMBOLGET("dX") || 1e-2;
      dF            = SYMBOLGET("dF") || 1e-3;
      nProj         = SYMBOLGET("nProjectors") || 10;
      dEnergy       = SYMBOLGET("dEnergy") || 1e-4;
      maxStepLength = SYMBOLGET("maxStepLength") || 0.3;
      SYMBOLGROUP("ric")  {
         invBric.ric.maxDist         = SYMBOLGET("maxDist")           || 10.;
         invBric.ric.primarySetLimit = SYMBOLGET("typifyThreshold")   || 0.05;
         invBric.ric.rmsThreshold    = SYMBOLGET("rmsThreshold") || 3.;
         invBric.ric.planeCutLimit   = SYMBOLGET("planeCutLimit")     || 0.95;
         invBric.ric.verbose = true;
         withAngle = SYMBOLGET("withAngles").toBool ();
         bvkAtoms << SYMBOLGET("bvkAtoms");
      }
      invBric.softModeDamping = SYMBOLGET("softModeDamping") || 1e-2;
   }
   for (int i = 0; i < bvkAtoms.getSize (); i++) bvkAtoms(i)--;

   cout << SX_SEPARATOR;
   cout << "| Quasi Newton (RIC + BFGS scheme)\n";
   cout << SX_SEPARATOR;
   cout << "|  max. steps                 :    " << maxSteps << endl;
   cout << "|  dX (max)                   :    " << dX << endl;
   cout << "|  dF (max)                   :    " << dF << endl;
   cout << "|  dEnergy                    :    " << dEnergy << endl;
   cout << "|  projector pairs            :    " << nProj << endl;
   cout << "|  max. step length           :    " << maxStepLength << endl;
   cout << SX_SEPARATOR;
   cout << "|  soft mode damping          :    " << invBric.softModeDamping << endl;
   cout << SX_SEPARATOR;
   cout << "|  RIC setup parameters" << endl;
   cout << "|  max. distance        :    " << invBric.ric.maxDist << endl;
   cout << "|  typify threshold     :    " << invBric.ric.primarySetLimit<<endl;
   cout << "|  type join rms factor :    " << invBric.ric.rmsThreshold << endl;
   cout << "|  polyhedron scaling   :    " << invBric.ric.planeCutLimit << endl;
   cout << SX_SEPARATOR;

   if (!calc)  return;

   invBric.setup (bondForceDefault, withAngle,angleForceDefault, bvkAtoms);

   // --- setup force filter
   SxOperator<SxAtomicStructure> F = getFilter (cmd);

   // --- variables for convergence
   SxVector<TPrecTauR> s2, g2;
   double sMax, gMax, sAvg, gAvg, E, deltaE;

   // --- initialize forces and inverse Hessian
   SxVector<TPrecTauR> x = structure.coordRef(), s, y;
   SxAtomicStructure g, gOld, force;
   SxMatrix<TPrecTauR> H = invBric.getHcart ();

   FILE *relaxHistFile = NULL;
   
   force = getForces ();
   g  = F | -force;

   g2 = g.absSqr ();
   gMax = sqrt(g2.maxval ());
   gAvg = sqrt(g2.sum () / (double)g2.getSize ());

   E = potentialPtr->getPotentialEnergy ();
   sxprintf ("ricQN[0]: max(|f|)=%12.8f, max(|dx|)=%12.8f, E=%12.8f H\n",
           gMax, 0., E);
   sxprintf ("rms(f)=%12.8f, rms(dx)=%12.8f\n", gAvg, 0.);

   SX_MPI_MASTER_ONLY
   {
      SxFileIO::appendToFile
      (  SxString(0)                     + "\t"  // iteration
       + SxString(E,"%15.12f")           + "\n"  // E_tot  [H]
      , "energy-structOpt.dat");

      relaxHistFile = fopen("relaxHist.sx", "a");
      if (!relaxHistFile)  {
         cout << "WARNING: Failed to open relaxHist.sx for writing" << endl;
      } else {
         fprintf (relaxHistFile, "// --- step 0\n");
         structure.fprint (relaxHistFile, force);
         fflush (relaxHistFile);
      }
   }

   SxList<SxVector<TPrecTauR> > sk,yk;

   SxList<SxAtomicStructure> configs, forces;

   configs << SxAtomicStructure(structure, SxAtomicStructure::Copy);
   forces << SxAtomicStructure(force, SxAtomicStructure::Copy);

   // --- optimization loop
   for (int it=0; it < maxSteps; it++)  {

      cout << SX_SEPARATOR;

      // ref1 (30) and (28)
      s = stepFilter | -(H.inverse () ^ g.coordRef ());

      // calculate convergence parameters
      s2 = SxAtomicStructure(s,SxAtomicStructure::Reference).absSqr ();
      sMax = sqrt(s2.maxval ());
      sAvg = sqrt(s2.sum () / (double)s2.getSize ());

      if (sMax > maxStepLength)  {
         // solve s = (H + lambda E)^{-1} g
         // chose lambda such that |s| = sMax * (99%..100%)
         double lambdaMin = 0., lambdaMax = 1.;
         // find upper bound of lambda
         do {
            SxMatrix<Double> Hl = H.getCopy ();
            SX_LOOP(i) Hl(i,i) += lambdaMax;
            s = stepFilter | -(Hl.inverse () ^ g.coordRef ());
            s2 = SxAtomicStructure(s,SxAtomicStructure::Reference).absSqr ();
            sMax = sqrt(s2.maxval ());
            cout << lambdaMax << " " << sMax << endl;
            if (sMax <= maxStepLength) break;
            lambdaMin = lambdaMax;
            lambdaMax *= 2.;
         } while (sMax > maxStepLength);
         // interval section
         double lambda;
         while (sMax > maxStepLength || sMax < 0.99 * maxStepLength)  {
            lambda = 0.5 * (lambdaMin + lambdaMax);
            SxMatrix<Double> Hl = H.getCopy ();
            SX_LOOP(i) Hl(i,i) += lambda;
            s = stepFilter | -(Hl.inverse () ^ g.coordRef ());
            s2 = SxAtomicStructure(s,SxAtomicStructure::Reference).absSqr ();
            sMax = sqrt(s2.maxval ());
            cout << lambda << " " << sMax << endl;
            if (sMax > maxStepLength)  {
               lambdaMin = lambda;
            } else {
               lambdaMax = lambda;
            }
         }
         sAvg = sqrt(s2.sum () / (double)s2.getSize ());
         cout << "limit step length: damping parameter is " << lambda << endl;
      }

      // --- get new structure and forces
      x += s;                      // ref1 (28) and (30)
      //cout << "X(" << it << ")\n" << x << endl;

      // --- print structure after every structOpt step
      SX_MPI_MASTER_ONLY {
         FILE *fp = fopen("relaxedStr.sx", "w");
         if (fp)  { structure.fprint(fp); fclose(fp); }
      }

      gOld = g;
      force = getForces ();
      g  = F | -force;

      //cout << "Gradient:" << g << endl;

      if (relaxHistFile)  {
         fprintf (relaxHistFile, "// --- step %d\n", it + 1);
         structure.fprint (relaxHistFile, force);
         fflush (relaxHistFile);
      }

      cout << SX_SEPARATOR;
      if (s.normSqr () > 1e-10)
         cout << "| Residual force along optimization line: "
              << (dot(s,g.coordRef ()) / s.norm ())
              << " from " << (dot(s,gOld.coordRef ()) / s.norm ()) << " ("
              << (100. * dot(s,g.coordRef ()) / dot(s, gOld.coordRef ()))
              << "%)" << endl;

      configs << SxAtomicStructure(structure, SxAtomicStructure::Copy);
      forces << SxAtomicStructure(force, SxAtomicStructure::Copy);
      if (configs.getSize () > nProj)  {
         configs.removeFirst ();
         forces.removeFirst ();
      }

      y  = g.coordRef () - gOld.coordRef ();  // ref1 (31)

      deltaE = potentialPtr->getPotentialEnergy () - E;
      E  = potentialPtr->getPotentialEnergy ();
      {
         double dEEstimate = 0.5 * dot(s,gOld.coordRef ());
         double dEHarmonic = 0.5 * dot(s,(g + gOld).coordRef ());
         cout << "| estimated energy gain: " << dEEstimate << endl;
         cout << "| obtained energy gain:  " << deltaE << endl;
         cout << "| difference:            " << (deltaE - dEEstimate) << endl;
         cout << "| harmonic energy gain:  " << dEHarmonic << endl;
         cout << "| difference:            " << (deltaE - dEHarmonic) << endl;
      }
      cout << SX_SEPARATOR;

      // calculate convergence parameters
      g2 = g.absSqr ();
      gMax = sqrt(g2.maxval ());
      gAvg = sqrt(g2.sum () / (double)g2.getSize ());


      // --- print out, and to file
      sxprintf ("ricQN[%d]: max(|f|)=%12.8f, max(|dx|)=%12.8f, E=%12.8f H\n",
               it+1, gMax, sMax, E);
      sxprintf ("rms(f)=%12.8f, rms(dx)=%12.8f\n", gAvg, sAvg);

      SX_MPI_MASTER_ONLY {
         SxFileIO::appendToFile
         (  SxString(it+1)                  + "\t"  // iteration
          + SxString(E,"%15.12f")           + "\n"  // E_tot  [H]
         , "energy-structOpt.dat");
      }


      // --- convergence check
      if (    sMax         < dX      // maximum structural step
           && gMax         < dF      // maximum force
           && fabs(deltaE) < dEnergy )  {
         sxprintf ("Convergence reached\n"); break;
      }

      // --- ref1 (34)
      {
         SX_CLOCK (RicUpdate);
         if (dot(y,s) <= 0.) {
            // update doesn't preserve positive definiteness
            cout << "Suppressing update of approximate inverse Hessian ..."
                 << endl;
         } else {
            yk.append (y);
            sk.append (s);

            // --- check for and remove out-of-date data
            /*
            // --- does current B still predict s within reasonable accuracy?
            for (int k = 0; k < yk.getSize (); ++k)  {
               double relErr = ((invB ^ yk(k))-sk(k)).norm () / sk(k).norm ();
               //cout << "k=" << k << ": " << relErr << endl;
               if (relErr > 1.)  {
                  sk.remove(k);
                  yk.remove(k);
                  --k;
                  reinit=true;
               }
            }
            */
            // --- number of projectors exceeded?
            if (sk.getSize () > nProj)  {
               sk.removeFirst ();
               yk.removeFirst ();
            }

            {
               if (invBric.maxDistFromSetup (structure) > 1
                   ||  !invBric.ric.verifyClasses (structure))
               {
                  // redefine internal coordinates
                  invBric.setup(bondForceDefault, withAngle, angleForceDefault,
                                bvkAtoms);
                  /*
                  while (configs.getSize () > 5)  {
                     configs.removeFirst ();
                     forces.removeFirst ();
                  }
                  */
               }
               invBric.reparamRic (configs, forces);
               // --- reinitialize H
               H = invBric.getHcart ();
               SX_LOOP(i) H(i,i) += invBric.softModeDamping;
               SxList<SxVector<TPrecTauR> >::Iterator sIt, yIt;
               sIt = sk.begin ();
               yIt = yk.begin ();
               for (; yIt != yk.end (); ++yIt, ++sIt)  {
                  double ysInv = 1./dot (*yIt, *sIt);
                  SxVector<Double> t = H ^ (*sIt);
                  double sHsInv = 1./dot(*sIt,t);
                  SX_LOOP2(j,i)  {
                     H(i,j) += (*yIt)(i) * (*yIt)(j) * ysInv 
                             -     t(i)  *    t(j)   * sHsInv;
                  }
               }
            }
         }
      }
   }

   // --- print relaxed structure into a file
   SX_MPI_MASTER_ONLY {
      FILE *fp = fopen("relaxedStr.sx", "w");
      if (fp)  { structure.fprint (fp); fclose(fp); }
   }
   if (relaxHistFile) fclose (relaxHistFile);
}

void SxStructOpt::set (SxPotential *P, SxAtomicStructure &str)
{
   potentialPtr = P;
   structure    = str;
   species      = P->getSpeciesData ();

   nDoF = 3 * structure.getNAtoms();
}


//---------------------------------------------------------------------------
// Service routines
//---------------------------------------------------------------------------
SxAtomicStructure SxStructOpt::getForces ()
{
   SX_CHECK (potentialPtr);
   return potentialPtr->getSymForces (structure, elMinimCmds);
}

SxAtomicStructure SxStructOpt::getForces (const SxMatrix<TPrecTauR> &tau)
{
   SX_CHECK (tau.getSize() == structure.coords.getSize(),
             tau.getSize(),   structure.coords.getSize());
   // TODO: ugly!!!
   structure.coords <<= tau;
   return getForces ();
}

SxAtomicStructure SxStructOpt::getForces (const SxAtomicStructure &tau)
{
   SX_CHECK (tau.getNAtoms() == structure.getNAtoms(),
             tau.getNAtoms(),   structure.getNAtoms());

   structure = SxAtomicStructure (tau, SxAtomicStructure::Copy);

   return getForces ();
}
