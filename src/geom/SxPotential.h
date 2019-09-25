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

#ifndef _SX_POTENTIAL_H_
#define _SX_POTENTIAL_H_

#include <SxAtomicStructure.h>
#include <SxVDW.h>
#include <SxSpeciesData.h>
#include <SxSymbolTable.h>
#include <SxForceSym.h>
#include <SxGeom.h>
#include <SxHamiltonian.h>

/** \brief Interface between energy/forces and structure-related algrotihms

    \b SxPotential = SPHInX Potential Interface

    This interface is required in order to use structure related algorithms,
    such as structure optimization or molecular dynamics, together with
    different potentials (e.g. ab-initio potential = SxHamSolver or
    empirical potentials).
    In general a potential returns the energy and forces at a given atomic
    structure.

    The main routines are SxPotential::getForces and SxPotential::getSymForces.

    \par Introducing of a new potential

    A new potential can easily be introduced by deriving from SxPotential.
    The purely abstract functions SxPotential::isRegistered and
    SxPotential::execute have to be overloaded in the derived class.
    Input parameters can be arranged in one or many subgroup(s) of
    main.bornOppenheimer.
    Note also, that the destructor of the derived class must be virtual due
    to the virtually declared functions in SxPotential.

    \par Example:

\code
SxMyPotential : public SxPotential
{
   public:
      SxMyPotential () : SxPotential ()  {  }
      virtual ~SxMyPotential () { }

      virtual isRegistered (const SxSymbolTable *cmd) const
      {
         SxString str = cmd->getName ();
         return (  str == "SteepestDescent" || str == "ConjugateGradient" );
      }


      virtual execute (const SxSymbolTable *cmd, bool calc)
      {
         SxString str = cmd->getName ();
         if      (str == "SteepestDescent")    steepestDescent (cmd, calc);
         else if (str == "ConjugateGradient")  conjugateGradient (cmd, calc);
         else  {
            SX_EXIT;  // unknown command
         }
      }


      void steepestDescent (const SxSymbolTable *cmd, bool calc)
      {
          try  {
             // read in all input parameters from cmd
          } catch (SxException e)  {
             e.print ();
             return;
          }

          // --- print input parameters
          cout << SX_SEPARATOR;
          cout << "| Steepest Descent\n";
          cout << SX_SEPARATOR;
          ...
          cout << SX_SEPARATOR;

          if (!calc)  return

          // --- perform actual calculation
      }

      ...


};
\endcode

    \ingroup group_structure
    \author  Sixten Boeck, boeck@mpie.de */
class SX_EXPORT_GEOM SxPotential
{
   public:
      SxPotential ();
      virtual ~SxPotential ();

      /** \brief Check whether the provided minimization command is known

          Potentials typically will get the content of the BornOppenheimer group
          for the force calculation. For example, electronic minimizers for DFT
          will be in that group. Other potentials such as empirical potentials
          do not need that.

          To produce better error messages, potentials that do not need any
          symbol table for force calculation must return "true" even for NULL.
          Potentials requiring a symbol table must return "false" for NULL.
          Potentials that optionally use a symbol table must return "true"
          for NULL.
     */
      virtual bool isRegistered (const SxSymbolTable *) const=0;

      /** \brief Execute the minimization

          \param cmd   minimization command
          \param calc  if \b false the command's input parameter are printed
                       only. */
      virtual void execute (const SxSymbolTable *cmd, bool calc=true)=0;

      /** \brief compute (unsymmetrized) forces.

          Compute the unsymmetrized forces to each atom according to the
          provided minimization command. Note, usually the forces needs
          to be symmetrized according the the symmetry elements exsiting
          in the current system. In this case SxPotential::getSymForces
          must be used instead.

          \param tau       atomic stucture \f$ \{ \tau_{i_s i_a} \} \f$
          \param cmd       minimization command group */
      virtual SxAtomicStructure getForces (const SxAtomicStructure &tau,
                                           const SxSymbolTable *cmd=NULL)=0;


      /** \brief Returns energy of the current Born-Oppenheimer update

          This function returns the energy computed during the current
          Born-Oppenheimer update, i.e., during the previous call of
          SxPotential::getForces or SxPotential::getSymForces

          \return energy in Hartree */
      virtual PrecEnergy getEnergy () const=0;


      /** \brief compute (unsymmetrized) forces.

          Compute the unsymmetrized forces to each atom according to the
          provided minimization command list and optionaly apply
          van-der-Waals correction on top. Note, usually the forces needs
          to be symmetrized according the the symmetry elements exsiting
          in the current system. In this case SxPotential::getSymForces
          must be used instead.

          \param tau       atomic stucture \f$ \{ \tau_{i_s i_a} \} \f$
          \param cmds      list of minimization command group */
      SxAtomicStructure getForces (const SxAtomicStructure &tau,
                                   const SxArray<const SxSymbolTable *> &cmds);

      /** \brief getting potential energy, including optional van-der-waals
                 correction*/

      double getPotentialEnergy ();


      /** \brief Force symmetrize filter
        */
      SxForceSym forceSymmetrizer;
      /** \brief compute (symmetrized) forces.

          Compute the unsymmetrized forces to each atom according to the
          provided minimization command using SxPotential::getForces.
          The unsymmetrized forces are symmetrized according to the
          symmetry elements existing in the system.

          \param tau       atomic stucture \f$ \{ \tau_{i_s i_a} \} \f$
          \param cmd       minimization command group */
      SxAtomicStructure getSymForces (const SxAtomicStructure &tau,
                                      const SxSymbolTable *cmd);

      /** \brief compute (symmetrized) forces.

          Compute the unsymmetrized forces to each atom according to the
          provided minimization command list using SxPotential::getForces.
          The unsymmetrized forces are symmetrized according to the
          symmetry elements existing in the system.

          \param tau       atomic stucture \f$ \{ \tau_{i_s i_a} \} \f$
          \param cmds      list of minimization command group */
      SxAtomicStructure
      getSymForces (const SxAtomicStructure &tau,
                    const SxArray<const SxSymbolTable *> &cmds);


      /** \brief Return species information

          Every potential type needs specific data about the contributing
          species. Hence, the specific potential (derived class from
          SxPotential) creates its own species data (a derived class of
          SxSpeciesData) on start-up, e.g., by reading these data from
          the input file.
          This function gives access to the calling routine (usually a
          member class related to structure optimization or molecular
          dynamics) */
      virtual SxSpeciesData getSpeciesData () const=0;


      /** \brief Extract all registered commands of a minimizer

          This function extracts all registered command groups of the
          derived class. It allows to perform a series of commands in order
          to reach the Born-Oppenheimer surface. */
      SxArray<const SxSymbolTable *> getMinimCmds (const SxSymbolTable *) const;

      /** \brief Print algorithm data only.

          This function calls the corresponding minimizer without performing
          the actual computation. It only reads the provided input parameters
          and prints them.

          \param cmd  minimization command
          \sa SxPotential::execute */
      void print (const SxSymbolTable *cmd) { execute (cmd, false); }

//       /**\brief apply van-der-Waals correcion */
     bool applyVDWCorrection;
//      /** \brief additive van-der-Waals correction*/
     SxVDW vdwCorrection;
//      /** \brief sets the van-der-Waals correction*/
//     void setVDWCorrection (const SxVDW &);

      bool dEnergyLow;

      /// Set to true if net force should always be zero
      bool noNetForce;

   protected:
      /** \brief checks if a structure is consistent with its symmetry elements

          This function is to check sure if the provided atomic structure
          is consistent with the existing symmetries in the crystal.
          It is called in the DEBUG mode in all SxPotential::getForces and
          SxPotential::getSymForces \b before the actual minimization is performed.

          If it fails it indicates that there is a problem in the calling
          structure optimization or molecular dynamics routine.  */
      bool isSymmetrizedStructure (const SxAtomicStructure &tau) const;



};

#endif /* _SX_POTENTIAL_H_ */
