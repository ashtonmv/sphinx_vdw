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

#include <SxForstData.h>

SxForstData::SxForstData () : SxSpeciesData ()
{
   // empty
}


SxForstData::SxForstData (const SxSymbolTable *table) : SxSpeciesData (table)
{
   try  {
      // someValue = table->get("someIdentifier")->toDouble ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}


SxForstData::~SxForstData ()
{
   // empty
}



