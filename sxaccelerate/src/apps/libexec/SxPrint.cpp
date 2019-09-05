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

#include <SxCLI.h>
#include <SxParser.h>
#include <SxUtil.h>

int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   SxString input = cli.option ("-i|--input", "file",
                                "SFHIngX input file").toString("input.sx");
   bool printPath   = cli.option ("-p|--path", "Print search path").toBool ();
   bool printTable  = cli.option ("-t|--table", 
                                  "Print input file in table form").toBool ();
// bool printBuffer = cli.option ("-a|--ASCII", 
//                                "Print file in ASCII form with all includes "
//                                "being unrolled."
//                               ).toBool ();
   bool printHash   = cli.option ("--perlhash",
                                  "Print a hash, which can be included from "
                                  "perl using the eval() function.").toBool ();
   bool printNested = cli.option ("-n|--nestedFiles",
                                  "Print a list of nested files").toBool ();
   bool skipCheck   = cli.option ("-s|--skipValidation", 
                                  "Skip grammar validation of input file"
                                 ).toBool ();
   bool printValidator 
      = cli.option ("--print-validator|-V",
                    "print validator file").toBool ();
   cli.finalize ();

   SxParser parser;
   if (printPath)  {
      cout << parser.getSearchPath () << endl;
      return 0;
   }

   parser.setValidation (!skipCheck && !printValidator);
   SxParser::Table table = parser.read (input);
   if (printValidator)  {
      SxString fileName = table->get("validator")->toString ();
      cout << "Validator file: " << SxParser_findInPath(fileName) << endl;
      return 0;
   }

   if (printNested)  {
      SxList<SxFile> nestedFiles = parser.getNestedFiles ();
      SxList<SxFile>::Iterator it;
      for (it = nestedFiles.begin(); it != nestedFiles.end(); ++it)  {
         cout << (*it).getName() << endl;
      }
   }

   if (printTable)  table->print ();
   if (printHash)   table->printHash ();
// if (printBuffer) cout << SxParser_buffer << endl;

   return 0;
}

