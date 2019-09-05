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

#include <SxSpeciesData.h>
#include <SxElemDB.h>

SxSpeciesData::SxSpeciesData ()
{
   // empty
}

SxSpeciesData::SxSpeciesData (const SxSymbolTable *table)
{
   readSpecies (table);
}

SxSpeciesData::SxSpeciesData (const SxList<SxString> &chemNames_)
{
   int i;
   ssize_t nSpecies = chemNames_.getSize();
   chemName.resize (nSpecies);
   elementName.resize (nSpecies);
   nuclearCharge.resize (nSpecies);
   valenceCharge.resize (nSpecies);
   dampingMass.resize (nSpecies);
   reciprocalMass.resize (nSpecies);
   ionicMass.resize (nSpecies);

   SxList<SxString>::ConstIterator it = chemNames_.begin();
   for (i=0; i < nSpecies; i++, it++)  {
      chemName(i)       = *it;
      elementName(i)    = *it;
      nuclearCharge(i)  = -1.;
      valenceCharge(i)  = -1.;
      dampingMass(i)    = -1.;
      reciprocalMass(i) = -1.;
      ionicMass(i)      = -1.;
   }
}

SxSpeciesData::~SxSpeciesData ()
{
   // empty
}


const SxSymbolTable* SxSpeciesData::findTable (const SxSymbolTable *tableIn)
{
   SxString tableName = tableIn->getName ();
   SxList<SxString> possibleName;
   possibleName << "pseudoPot" // pseudopotential
                << "pawPot"    // PAW potential
                << "eamPot"    // EAM potential
                << "skData"    // DFTB potential (slater-koster files)
                ;

   // --- check if tableIn is or contains a potential group
   for (int i = 0; i < possibleName.getSize (); ++i)  {
      if (possibleName(i) == tableName) return tableIn;
      if (tableIn->containsGroup(possibleName(i)))
         return tableIn->getGroup(possibleName(i));
   }
   if (tableIn == tableIn->topLevel ())  {
      SxString validNames = possibleName(0);
      for (int i = 1; i < possibleName.getSize (); ++i)
         validNames += ", " + possibleName(i);

      throw SxException (("No potential group (" + validNames +") found."
                         ).ascii(),
                         __FILE__, __LINE__);
   }
   // try top level
   return findTable(tableIn->topLevel ());
}


void SxSpeciesData::readSpecies (const SxSymbolTable *tableIn)
{
   SX_CHECK (tableIn);
   try {
      const SxSymbolTable *table = tableIn->containsGroup("species")
                                 ? tableIn
                                 : findTable(tableIn);
      const SxSymbolTable *strSpecies = NULL;
      if (tableIn->topLevel ()->containsGroup ("structure"))  {
         strSpecies = tableIn->topLevel ()->getGroup ("structure")
                      ->getGroup ("species");
      }

      SxSymbolTable *species;
      int i, nSpecies = table->getGroup("species")->getNItems ("species");
      if (strSpecies && strSpecies->getNItems ("species") != nSpecies)  {
         cout << "Species mismatch between " << table->getName ()
              << " (" << nSpecies << " species) and structure group ("
              << strSpecies->getNItems ("species") << ")." << endl;
         SX_QUIT;
      }
      // --- allocate memory
      chemName.resize (nSpecies);
      elementName.resize (nSpecies);
      nuclearCharge.resize (nSpecies);
      valenceCharge.resize (nSpecies);
      dampingMass.resize (nSpecies);
      reciprocalMass.resize (nSpecies);
      ionicMass.resize (nSpecies);

      // --- read species data
      for (species  = table->getGroup("species"), i=0;
           species != NULL;
           species  = species->nextSibling("species"), i++)
      {
         if (species->contains("element"))  {
            chemName(i) = species->get("element")->toString().trim();
            // --- check that structure has same species ordering
            if (strSpecies && strSpecies->contains("element"))  {
               SxString strElement = strSpecies->get("element")->toString ();
               if (strElement != chemName(i))  {
                  cout << "Species mismatch between " << table->getName ()
                       << " and structure group for species " << (i+1) 
                       << "!" << endl;
                  cout << table->getName () << ": " << chemName(i) << endl;
                  cout << "structure: " << strElement << endl;
                  SX_QUIT;
               }
            }
         } else {
            chemName(i) = "unknown-" +  SxString(i);
         }
         if (species->contains("name"))
            elementName(i) = species->get("name")->toString().trim();
         else
            elementName(i) = chemName(i);
         nuclearCharge(i)  = -1; // TODO: to be implemented
         valenceCharge(i)  = species->contains ("valenceCharge")
                           ? species->get("valenceCharge")->toReal()
                           : 0.;
         dampingMass(i)    = species->contains("dampingMass")
                           ? species->get("dampingMass")->toReal()
                           : 0.7;
         reciprocalMass(i) = species->contains("reciprocalMass")
                           ? species->get("reciprocalMass")->toReal()
                           : -1.;
         ionicMass(i)      = species->contains ("ionicMass")
                           ?  species->get("ionicMass")->toReal()
                           : reciprocalMass(i);
         if (    species->contains("ionicMass")
             && !species->contains ("reciprocalMass"))
            reciprocalMass(i) = ionicMass(i);
         if (!species->contains("ionicMass")
             && !species->contains("reciprocalMass"))
         {
            // inspect the database
            SxString symbol = chemName(i);
            if (symbol.getSize () > 2) symbol.resize (2, true);
            if (symbol.getSize () == 2
                && (symbol(1) < 'a' || symbol(1) > 'z'))
                symbol.resize (1, true);
            const SxElemDB database = SxElemDB::getElemDB ();
            double mass = database.getAtomicWeight (symbol);
            if (fabs(mass) < 1e-10)  {
               cout << "WARNING: failed to look up element '"
                    << chemName(i) << "'";
               if (chemName(i) != symbol)
                  cout << " as chemical symbol '" << symbol << "'";
               cout << " in element database" << endl;
               cout << "Ionic mass is not known" << endl;
            } else {
               cout << "Lookup ionic mass of " << chemName(i);
               if (chemName(i) != symbol)
                  cout << " (=" << symbol << ")";
               cout << ": " << mass << endl;
               ionicMass(i) = reciprocalMass(i) = mass;
               // also set chemical name, if we are here anyway
               if (elementName(i) == chemName(i))
                  elementName(i) = database.getName (database.getIdx(symbol));
            }
         }
         if (strSpecies) strSpecies = strSpecies->nextSibling ("species");
      }
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
}


int SxSpeciesData::getNSpecies () const
{
   return int(elementName.getSize());
}

int SxSpeciesData::find (const SxString &name) const
{
   if (name.getSize () == 0) return -1;
   int is;
   // try chemName
   for (is = 0; is < chemName.getSize (); ++is)
      if (chemName(is) == name) return is;
   // try elementName
   for (is = 0; is < elementName.getSize (); ++is)
      if (elementName(is) == name) return is;
   // try number
   try {
      is = name.toInt () - 1;
      if (is < 0 || is >= getNSpecies ()) is=-1;
   } catch (SxException e)  {
      // name is not a number
      is = -1;
   }
   return is;
}

SxString SxSpeciesData::prettyName (int is) const
{
   SX_CHECK (is >= 0 && is < getNSpecies ());
   if (elementName(is).getSize () > 0) return elementName(is);
   if (chemName(is).getSize () > 0)    return chemName(is);
   return SxString ("species ") + (is+1);
}

void SxSpeciesData::PrintInfo::setup (const SxParser &parser,
                                      const SxSymbolTable *table)
{
   const SxSymbolTable *species;
   try  {
      if (table->getName () == "species") {
         species = table;
         groupName = table->parent->getName ();
      } else if (table->containsGroup("species"))  {
         species = table->getGroup("species");
         groupName = table->getName ();
      } else {
         const SxSymbolTable *potGroup = findTable(table);
         species = potGroup->getGroup("species");
         groupName = potGroup->getName ();
      }

      int nSpecies = species->getNItems ("species");
      
      includeFiles.resize (nSpecies);
      
      for (int is = 0; is < nSpecies; 
           ++is, species = species->nextSibling ("species"))
      {
         includeFiles(is) = parser.locateFile (species->get("valenceCharge"))
                           .getIncludeName (parser.getSearchPath ());
      }
   } catch (SxException e)  {
      e.print ();
      cout << "Cannot guess species include files." << endl;
   }
}


SxArray<SxString>
SxSpeciesData::guessIncludeFiles(const SxParser &parser, SxParser::Table table)
{
   cout << "WARNING: use of deprecate SxSpeciesData::guessIncludeFiles" << endl;
   return PrintInfo(parser, &*table).includeFiles;
}


void SxSpeciesData::PrintInfo::fprint (FILE *output) const 
{
   SX_CHECK(output);
   if (!isValid ()) return;
   if (groupName == "structure")  {
      cout << "WARNING: suppressing deprecate-style potential output "
              "nested into structure {}. Use pseudoPot {} instead!" 
           << endl;
      return;
   }
   sxfprintf (output, "%s {\n", groupName.ascii ());
   int nSpecies = int(includeFiles.getSize ());
   for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies)  {
      sxfprintf (output, "   species { include ");
      if (includeFiles(iSpecies).getSize () > 0)
         sxfprintf(output, "<%s>; }\n", includeFiles(iSpecies).ascii ());
      else
         sxfprintf(output, "   species_%d; }\n", iSpecies+1);
   }
   sxfprintf(output, "}\n\n");
   
}

SxArray<SxString> SxSpeciesData::getElements (const SxSymbolTable *table)
{
   const SxSymbolTable *group = NULL;
   try {
      group = findTable (table);
   } catch (SxException e) {
      const SxSymbolTable *top = table->topLevel ();
      if (top->containsGroup ("structure"))
         group = top->getGroup ("structure");
      else  {
         e.print ();
         cout << "Couldn't find group to read element names from ..." << endl;
         SX_QUIT;
      }
   }
   SxArray<SxString> elements;
   try {
      const SxSymbolTable *species = group->getGroup("species");
      elements.resize (species->getNItems ("species"));
      // --- read species data
      char dummy ='a';
      for (int i=0;
           species != NULL;
           species  = species->nextSibling("species"), i++)
      {
         if (species->contains("element"))
            elements(i) = species->get("element")->toString().trim();
         else  {
            elements(i) = "Xa"; // dummy name
            elements(i)(1) = dummy; // Xa .. Xz (and then the rest of ASCII)
            dummy++;
         }
      }
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
   return elements;
}

void SxSpeciesData::read (const SxBinIO &io)
{
   SX_CHECK (io.mode == SxBinIO::BINARY_READ_ONLY);
   try {
      SxString chemNameList;
      io.read("chemNames", &chemNameList);
      chemName = chemNameList.tokenize (',');
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   elementName = chemName;
   ssize_t nSpecies = chemName.getSize ();
   nuclearCharge.resize (nSpecies); nuclearCharge.set (0.);
   valenceCharge.resize (nSpecies); valenceCharge.set (0.);
   dampingMass.resize (nSpecies);   dampingMass.set (0.);
   ionicMass.resize (nSpecies);     ionicMass.set (0.);
}
