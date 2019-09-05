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


#include <SxForst.h>

SxForst::SxForst () 
   : SxPotential (),
     speciesData ()
{
   printf ("SxForst::SxForst\n");
}


SxForst::SxForst (const SxAtomicStructure &/*tau*/, const SxSymbolTable *table) 
   : SxPotential (),
     speciesData (table)
{
   printf ("SxForst::SxForst\n");
}


SxForst::~SxForst ()
{
   // empty
}


bool SxForst::isRegistered (const SxSymbolTable *cmd) const
{
   SxString str = cmd->getName ();
   printf (">>>>>%s<<<<<\n", str.ascii());
   return ( str == "FORST" );
}


void SxForst::execute (const SxSymbolTable *cmd, bool calc)
{
   SX_CHECK (cmd);

   // --- parsing input data
   int    someInt    = cmd->get ("someInt")->toInt ();
   double someDouble = cmd->get ("someDouble")->toReal ();
   SxString someStr  = cmd->get ("someString")->toString ();

   // --- print input data
   cout << SX_SEPARATOR;
   cout << "| Forst\n";
   cout << SX_SEPARATOR;
   cout << "| Parameters:\n";
   cout << "|    int:         " << someInt    << endl;
   cout << "|    double:      " << someDouble << endl;
   cout << "|    string:      " << someStr    << endl;
   cout << SX_SEPARATOR;

   if (!calc)  return;

   // --- calculation starts here

   // ...
}


SxAtomicStructure SxForst::getForces (const SxAtomicStructure &tau,
                                      const SxSymbolTable *cmd)
{
   // atomic coordinates -> tau
   execute (cmd);

   // update tau

   return tau;
}


PrecEnergy SxForst::getEnergy () const
{
   return 0.;
}


SxSpeciesData SxForst::getSpeciesData () const
{
   return speciesData;
}
