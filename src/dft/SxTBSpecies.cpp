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


#include <SxTBSpecies.h>

SxTBSpecies::SxTBSpecies ()
   : SxSpeciesData ()
{
   // empty
}

SxTBSpecies::~SxTBSpecies ()
{
   // empty
}

SxTBSpecies::SxTBSpecies (const SxSymbolTable *tableIn)
   : SxSpeciesData ()
{
   int iSpecies, nSpecies;
   SxSymbolTable *species;
   try {
      const SxSymbolTable *table = tableIn->useDeprecateFormat() 
                                 ? tableIn->getGroup("structure")
                                 : tableIn->getGroup("skData");

      SxSpeciesData::readSpecies (table);

      nSpecies = table->getGroup("species")->getNItems ("species");
      skElementName.resize (nSpecies);
      lMax.resize (nSpecies);

      if (!table->useDeprecateFormat())  {
         for (species  = table->getGroup("species"), iSpecies=0;
               species != NULL;
               species  = species->nextSibling("species"), iSpecies++)
         {
            skElementName(iSpecies) = species->get("skElementName")->toString();
            lMax(iSpecies)          = species->get("lMax")->toInt();
            //valence charge is already in the constructor of SxSpeciesData 
            //valenceCharge(iSpecies)  = species->get("valenceCharge")->toInt();
            //get occupation numbers?
         }
         skFilesPath   = table->get("skFilesPath")->toString();
         // Ask for fileformat
         if (table->contains("skFilesFormat"))
            skFilesFormat = table->get("skFilesFormat")->toInt();
         // Standard is S/PHI/nX file format 
         else skFilesFormat = 2;
      }  else {
         sxprintf(" deprecate format is not supported in tight-binding\n");
         sxprintf(" please change your input file to the new format\n");
         SX_QUIT;
      }
      
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
}


