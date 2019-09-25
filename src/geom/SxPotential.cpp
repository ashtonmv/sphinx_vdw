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

#include <SxPotential.h>
#include <SxGrid.h>

SxPotential::SxPotential ()
   : dEnergyLow(false), noNetForce(true)
{
   // empty
}

SxPotential::~SxPotential ()
{
   // empty
}

SxAtomicStructure
SxPotential::getForces (const SxAtomicStructure &tau,
                        const SxArray<const SxSymbolTable *> &cmds)
{
   SxAtomicStructure f;
   for (int i=0; i < cmds.getSize(); i++)
         f = getForces (tau, cmds(i));

//   // --- extremely ugly, to be cleaned up

   const SxSymbolTable *hamGroup = SxHamiltonian::getHamiltonianGroup (cmds(0) -> topLevel());

/*
    if (SxHamiltonian::applyVDWCorrection) {
      SxArray<SxVector3<Double> > tauArray (tau.nTlAtoms);
      SxAtomicStructure help;
      SxAtomicStructure fVDW;

      vdwCorrection = SxVDW (tau, hamGroup);
      //setVDWCorrection(VDW);

      help.copy(tau);
      //TODO (urgent): switch to SxAtomicStructure
      for (int i = 0; i < tau.nTlAtoms; i++)
         tauArray(i) = help.ref(i);
      vdwCorrection.update (tauArray);
      vdwCorrection.compute ();
      tauArray = vdwCorrection.Forces;
      fVDW.copy(tau);
      for (int i = 0; i < tau.nTlAtoms; i++)
         fVDW.ref(i) = (tauArray(i));
      fVDW.nSpecies = f.nSpecies = tau.nSpecies;
      fVDW.nTlAtoms = f.nTlAtoms = tau.nTlAtoms;
     f = f + fVDW;
   }
*/
   return f;
}

double SxPotential::getPotentialEnergy ()
{
   double ePot = getEnergy ();
   if (applyVDWCorrection) {
      ePot = ePot + vdwCorrection.totalEnergy;
      /*
      cout << VDWCorrection.potentialType << endl;
      cout << VDWCorrection.getTotalEnergy () << endl;
      SX_EXIT;
      */
   }
   return ePot;
}


SxAtomicStructure SxPotential::getSymForces (const SxAtomicStructure  &tau,
                                             const SxSymbolTable      *table)
{
   SX_CHECK (isSymmetrizedStructure (tau));

   // --- compute unsymmetrized forces
   SxAtomicStructure f = getForces (tau, table);
   const SxSymbolTable *hamGroup = SxHamiltonian::getHamiltonianGroup( table -> topLevel() );

//   // --- extremely ugly, to be cleaned up
//   if (applyVDWCorrection) {
   if (hamGroup -> containsGroup("vdwCorrection")) {
      applyVDWCorrection = true;

      SxArray<SxVector3<Double> > tauArray (tau.nTlAtoms);
      SxAtomicStructure help;
      SxAtomicStructure fVDW;
      vdwCorrection = SxVDW (tau, hamGroup);
      //setVDWCorrection(VDW);
      help.copy(tau);

      //TODO (urgent): switch to SxAtomicStructure
      for (int i = 0; i < tau.nTlAtoms; i++)
         tauArray(i) = help.ref(i);

      vdwCorrection.update (tauArray);
      vdwCorrection.compute ();
      tauArray = vdwCorrection.Forces;
      fVDW.copy(tau);

      for (int i = 0; i < tau.nTlAtoms; i++)
         fVDW.ref(i) = (tauArray(i));

      fVDW.nSpecies = f.nSpecies = tau.nSpecies;
      fVDW.nTlAtoms = f.nTlAtoms = tau.nTlAtoms;
      f = f + fVDW;
   }

   // symmetrize forces according to existing symmetry elements

   if (forceSymmetrizer.getNSymmetries () == 0) {
      forceSymmetrizer.setup (tau);
   }

   return forceSymmetrizer | f;
}

SxAtomicStructure
SxPotential::getSymForces (const SxAtomicStructure &tau,
                           const SxArray<const SxSymbolTable *> &cmds)
{
   if (cmds.getSize () == 0) {
      SX_CHECK (isRegistered (NULL));
      return getSymForces (tau, NULL);
   }
   SxAtomicStructure f;
   for (int i=0; i < cmds.getSize(); i++) {
      f = getSymForces (tau, cmds(i));
   }
   return f;
}


SxArray<const SxSymbolTable *>
SxPotential::getMinimCmds (const SxSymbolTable *table) const
{
   SX_CHECK (table);
   SxList<const SxSymbolTable *> cmdList;
   try  {
      if (table->containsGroup ("bornOppenheimer"))
         table = table->getGroup("bornOppenheimer");

      SxSymbolTable *cmd   = NULL;
      for (cmd  = table->begin();
           cmd != NULL;
           cmd  = cmd->nextSibling())
      {
         if (isRegistered (cmd))  cmdList << cmd;
         else  {
            if (!table->useDeprecateFormat())  {
               sxprintf ("Unknown command.\n");
               cmd->print ();
               SX_EXIT;
            }
         }
      }

   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   return SxArray<const SxSymbolTable *> (cmdList);
}


bool SxPotential::isSymmetrizedStructure (const SxAtomicStructure &tau) const
{
   SX_CHECK (tau.cell.symGroupPtr);
   const SxSymGroup &S = *tau.cell.symGroupPtr;

   int nSym = S.getSize();

   // --- verify existing symmetry elements
   if (tau.getNAtoms () < 100)  {
      for (int iSym=0; iSym < nSym; iSym++)
         if (tau != (S(iSym) ^ tau) )
            return false;
   } else {
      SxGrid grid(tau, 3);
      for (int iSym=0; iSym < nSym; iSym++)
         if (! tau.match (grid, S(iSym) ^ tau) )
            return false;
   }

   if (forceSymmetrizer.getNSymmetries () != 0)
      if (!forceSymmetrizer.checkStr (tau)) return false;
   return true;

}

//void SxPotential::setVDWCorrection (const SxVDW &vdw)
//{
//  VDWCorrection = vdw;
//  applyVDWCorrection = true;
//}
