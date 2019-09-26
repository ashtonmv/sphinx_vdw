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

#include <SxHamiltonian.h>
#include <SxPWHamiltonian.h>
#include <SxPAWHamiltonian.h>
#include <SxKdotP8x8.h>
#include <SxKdotP.h>
#include <SxStrain.h>
#include <SxConstants.h>

const SxSymbolTable *
SxHamiltonian::getHamiltonianGroup (const SxSymbolTable *table)
{
   const SxSymbolTable *hGrp = NULL;
   if ( (hGrp = SxPWHamiltonian::getHamiltonianGroup(table)) ) return hGrp;
   else if ( (hGrp = SxPAWHamiltonian::getHamiltonianGroup(table)) ) 
      return hGrp;
   else if ( (hGrp = SxKdotP8x8::getHamiltonianGroup(table)) ) return hGrp;
   else if ( (hGrp = SxKdotP::getHamiltonianGroup(table)) ) return hGrp;
   else if ( (hGrp = SxStrain::getHamiltonianGroup(table)) )   return hGrp;
   cout << "Hamiltonian group not found." << endl;
   cout << "Perhaps new and deprecate input tags are used." << endl;
   SX_QUIT;
   return NULL;
}

int SxHamiltonian::getNSpin (const SxSymbolTable *table)
{
   int nSpin = 1;
   try  {
      bool spinPol = false;
      const SxSymbolTable *hamiltonian = getHamiltonianGroup (table);
      if (hamiltonian->contains("spinPolarized"))
         spinPol = hamiltonian->get("spinPolarized")->toAttribute();
      nSpin = (spinPol) ? 2 : 1;
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   return nSpin;
}

int SxHamiltonian::getNEmptyStates (const SxSymbolTable *table)
{
   int nEmptyStates = 0;
   try {
      const SxSymbolTable *hamiltonian = getHamiltonianGroup (table);
      if (hamiltonian->contains("nEmptyStates"))
         nEmptyStates = hamiltonian->get("nEmptyStates")->toInt();
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
   return nEmptyStates;
}

int SxHamiltonian::getNEmptyStatesChi (const SxSymbolTable *table)
{
   int nEmptyStatesChi = -1;
   try  {
      const SxSymbolTable *hamiltonian = getHamiltonianGroup (table);
      if (hamiltonian->contains("nEmptyStatesChi"))
         nEmptyStatesChi = hamiltonian->get("nEmptyStatesChi")->toInt();
   }  catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   return nEmptyStatesChi;
}

double SxHamiltonian::getNExcessElectrons (const SxSymbolTable *table)
{
   double nExcessElectrons = 0.;
   try {
      const SxSymbolTable *hamiltonian = getHamiltonianGroup (table);
      if (hamiltonian->contains("nExcessElectrons"))
         nExcessElectrons = hamiltonian->get("nExcessElectrons")->toReal();
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
   return nExcessElectrons;
}

double SxHamiltonian::getNElectrons (const SxSymbolTable *table,
                                            const SxArray<double> &valCharge,
                                            const SxAtomicStructure &str)
{
   if (valCharge.getSize () != str.getNSpecies ())  {
      cout << "Inconsistency:" << endl;
      cout << "Number of species mismatch  between structure and potential:"
           << endl;
      cout << "in structure:" << str.getNSpecies () << endl;
      cout << "in potential:" << valCharge.getSize () << endl;
      SX_QUIT;
   }
   // get extra electrons
   double nElectrons = SxHamiltonian::getNExcessElectrons (table);
   // add valence electrons
   for (int iSpecies=0; iSpecies < str.getNSpecies (); iSpecies++)
      nElectrons += valCharge(iSpecies) * str.getNAtoms(iSpecies);
   return nElectrons;
}

int SxHamiltonian::getNStates (const SxSymbolTable *table, double nElectrons,
                               double spinMoment)
{
   int nSpin            = SxHamiltonian::getNSpin (table);
   int nEmptyStates     = SxHamiltonian::getNEmptyStates (table);
   if (nSpin == 2 && fabs(spinMoment) < 1e-8)
      spinMoment=min(10., nElectrons); // 10 = just a large number
   int nStates = nEmptyStates 
               + int (ceil (nElectrons / 2. + fabs(spinMoment) / 2.));
   return nStates;
}

double SxHamiltonian::getEkt (const SxSymbolTable *table)
{
   double ekt = 0.;
   try {
      const SxSymbolTable *hamiltonian = getHamiltonianGroup (table);
      if (hamiltonian->contains("ekt"))
         ekt = hamiltonian->get("ekt")->toReal() / HA2EV;
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
   return ekt;
}

bool SxHamiltonian::withDipoleCorrection (const SxSymbolTable *table)
{
   bool res = false;
   try {
      const SxSymbolTable *hamiltonian = getHamiltonianGroup (table);
      if (hamiltonian->contains("dipoleCorrection"))
         res = hamiltonian->get("dipoleCorrection")->toAttribute();
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
   return res;
}

SxConstPtr<SxOverlapBase> SxHamiltonian::getS () const
{
   return SxPtr<SxPWOverlap>::create ();
}

bool SxHamiltonian::rereadTable (const SxSymbolTable *)
{
   return false;
}

void SxHamiltonian::backToDefault (const SxSymbolTable *)
{
   SX_EXIT;
}

