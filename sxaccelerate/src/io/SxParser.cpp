// ---------------------------------------------------------------------------
//
//           The general purpose cross platform C/C++ framework
//
//                       S x A c c e l e r a t e
//
//           Home:       https://www.sxlib.de
//           License:    Apache 2
//           Authors:    see src/AUTHORS
//
// ---------------------------------------------------------------------------

#include <SxParser.h>  // should become <SxFileParser.h>
#include <SxUtil.h>
#include <SxException.h>
#include <SxDir.h>
#ifdef WIN32
#  include <windows.h>
#  ifndef _MAX_PATH
#    define _MAX_PATH 10240
#  endif /* _MAX_PATH */
#endif /* WIN32 */

void initParser ()
{
   SxParser_init ();
}




SxParser::SxParser ()
   : freeFormat(false)
{
   // empty
}

SxParser::~SxParser ()
{
   // empty
}


void SxParser::setSearchPath (const SxString &path_)
{
   searchPath = path_;
}


SxConstPtr<SxSymbolTable> SxParser::read (const SxString &file, 
                                          const SxString &validator_)
{
   SX_TRACE ();
   SxParser_nestedFiles.removeAll ();
   
   SxString path = searchPath;
   if (searchPath == "")  path = getDefaultSearchPath();

   SxSymbolTable *&symTab = SxSymbolTable::getGlobalPtr ();
   if (symTab)  delete symTab; 

   SxPtr<SxSymbolTable> table;
   table = SxPtr<SxSymbolTable>::create ();
   symTab = &(*table);

   SxParser_filename = file;
   SxParser_init ();
   SxParser_parseFile (SxParser_findInPath(file, path)(0), path, excludeGroups);

   try {
      if (table->contains ("handleExceptions"))  {
         if (table->get("handleExceptions")->toAttribute() == false)  {
            SxException::causeSegFault ();
            sxprintf ("Exceptions will cause segmentation faults!\n");
         }
      }
   } catch (SxException)  {
      // empty
   }

   if (!freeFormat)  {
      try  {
         SxSymbolTable  valTable;
         symTab = &valTable;
         SxString validator;

         // --- user provided validator has higher priority
         if (table->contains("validator"))  {
            validator = table->get("validator")->toString();
         }  else  {
            if (validator_ != "")  
               validator = validator_;
            else  {
               sxprintf ("\n");
               sxprintf ("Error: No format statement found in file %s.\n",
                        file.ascii());
               sxprintf ("       Please provide the format validator in the\n");
               sxprintf ("       input file, e.g.\n");
               sxprintf ("          format myGrammar;\n");
               sxprintf ("       This statement would refer to the validator\n");
               sxprintf ("          $SPHINX/std/myGrammar.std\n");
               SX_QUIT;
            }
         }

         SxParser_filename  = SxParser_findInPath(validator, path)(0);
         SxParser_init ();
         SxString origBuffer = SxParser_buffer; // save original buffer
         SxParser_parseFile (SxParser_filename, path, excludeGroups);
         if (!valTable.contains("ignoreValidation"))  {
            table->validateTable (valTable);
         }
         SxParser_buffer = origBuffer;  // restore original buffer
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      } 
   }

   symTab  = NULL;
   // checkParams (getParam(), __lics);

   return table;
}

SxList<SxFile> SxParser::getNestedFiles () const
{
   return SxParser_nestedFiles;
}

SxFile SxParser::locateFile (const SxSymbol *symbol) const
{
   return (symbol) ? symbol->parserFilename : SxFile();
}

SxFile SxParser::locateFile (const SxSymbolTable *symbol) const
{
   return (symbol) ? symbol->parserFilename : SxFile();
}

int SxParser::locateLine (const SxSymbol *symbol) const
{
   return (symbol) ? symbol->parserLineNumber : 0;
}

int SxParser::locateLine (const SxSymbolTable *symbol) const
{
   return (symbol) ? symbol->parserLineNumber : 0;
}


void SxParser::setVerbose (bool verbose)
{
   SxParser_verbose = verbose;
}


SxString SxParser::getDefaultSearchPath ()
{
   SxString userPath = getenv ("SX_INCLUDE_PATH");
#  ifdef WIN32
      SxString sep = ";";
#  else
      SxString sep = ":";
#  endif /* WIN32 */
   if (!userPath)  {
      // --- try to figure out top directory of package
      try {
         SxString installPath = SxDir::getInstallPath();
         userPath  = installPath + "/share/sphinx";
         userPath += sep + installPath + "/share" + sep;
      } catch (SxException e) {
#ifndef SX_NOINSTALLPATH
         e.print ();
         SX_EXIT;
#endif
      }
      userPath += SxString(getenv ("HOME")) + "/SPHInX";
   }
   return userPath;
}


SxString SxParser::getSearchPath () const
{
   if (searchPath != "")  return searchPath;
   else                   return getDefaultSearchPath ();
}


void SxParser::setValidation (bool state)
{
   freeFormat = !state;
}


void SxParser::excludeGroup (const SxString &name)
{
   excludeGroups << name.tokenize (',');
}

