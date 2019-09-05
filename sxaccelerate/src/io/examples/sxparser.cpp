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

#include <SxUtil.h>
#include <SxSymbolTable.h>
#include <SxMatrix3.h>
#include <SxParser.h>
#include <SxCLI.h>
#include <stdio.h>
#include <iostream>
#include <cstdio>
#include <ctime>
/** \example parser.cpp

    This example demonstrates how to use the input file parser. Please 
    inspect to file parser.sx that will be read now.
    In particulary you learn how to combine the SFHIngX utility with
    the io class.
    \brief  Demonstration of the input file parser
    \sa     SxList
    \author Sixten Boeck
 */
int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   cli.finalize ();


   SxParser parser;
   SxParser::Table table = parser.read ("parser.sx");

   try {
      // --- read structure {...} group
      SxSymbolTable *structure = table->getGroup ("structure");
//ss
      // --- read scalar values
      cout << "i = " << structure->get("i")->toInt() << endl;
      cout << "d = " << structure->get("d")->toReal() << endl;
      cout << "s = " << structure->get("s")->toString() << endl;

      // --- read vectors/matrices 
      //cout << "v = " << SxMatrix3<Double> (structure->get("cell")->toList())
     //      << endl;
     // cout << "m = " << SxMatrix3<Double> (structure->get("cell")->toList())
     //      << endl;

      // --- dealing with optional values (here an attribute)
      bool doSomething = false;
      if (structure->contains("doSomething"))
         doSomething = structure->get("doSomething")->toAttribute();

      size_t nAtoms = 0;

      SxSymbolTable *species = structure->getGroup ("species");
   std::clock_t start;
   double duration;

   start = std::clock ();


      // --- iterate over all structure.atom {...} groups
      SxSymbolTable *atom = NULL;
      for (atom  = species->getGroup("atom");
           atom != NULL;
           atom  = atom->nextSibling ("atom"))
      {
         //cout << "ATOM: " << SxVector3<Double>(atom->get("coords")->toList())
           //   << "\t" << atom->get("label") << endl;
         nAtoms++;
      }
      std::cout << "Total atoms: " << nAtoms << std::endl;

      duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;

      std::cout << "time duration: " << duration << std::endl;


   } catch (SxException e)  {
      e.print();
   }
}

