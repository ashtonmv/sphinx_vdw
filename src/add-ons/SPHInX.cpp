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




#include <SPHInX.h>
#include <SxMemMonitor.h>
#include <SxAtomicStructure.h>
#include <SxSpeciesData.h>
#include <SxSymFinder.h>
#include <SxStructOpt.h>
#include <SxHamSolver.h>
#include <SxExxSolver.h>
#include <SxTBHamSolver.h>
#include <SxMolDyn.h>
#include <SxFrozenPhonon.h>
#include <SxUtil.h>
#include <SxEAM.h>
#include <SxMathLib.h>
#include <SxParser.h>
#include <SxConfig.h>
#include <SxLoopMPI.h>
#include <SxParallelHierarchy.h>


#ifndef SX_STANDALONE

SPHInX::SPHInX ()
{
   // empty
}

SPHInX::SPHInX (const SxString &inFile)
{
   file = inFile;
}


SPHInX::~SPHInX ()
{
   // empty
}


int SPHInX::compute ()
{
   // --- read from input file
   SxParser parser;
   SxParser::Table table = parser.read (file);

   // --- initialize atomic struture from input file
   SxAtomicStructure str (&*table);

   // --- find symmetry elements in current atomic structure
   //SxSymFinder symFinder;
   //symFinder.compute (str);
   //cout << "SYMMETRIES\n" << symFinder.rot << endl;

   // --- choose the Born-Oppenheimer solver
   SxPotential *potential = NULL;
   try  {
      if (table->containsGroup ("basis"))  {
         SxSymbolTable *cmd = table->getGroup ("main")->begin ();
         if (cmd == NULL) 
            potential = new SxHamSolver (str, &*table);
         else if (SxExxSolver::isExx(cmd))
            potential = new SxExxSolver (str, &*table);
         else
            potential = new SxHamSolver (str, &*table);
      } else if (table->containsGroup ("eamPot")) 
      {
         potential = new SxEAM (str, &*table);
      } else if (table->containsGroup ("tbBasis")) 
         potential = new SxTBHamSolver (str, &*table);
      else
         SX_EXIT;
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   try {
      // --- initial guess for ham solver
      if (table->containsGroup ("basis"))  {  
         if (table->containsGroup ("initialGuess"))  {
            SxSymbolTable *initialGuess = table->getGroup("initialGuess");
            potential->execute (initialGuess);
         } else {
            SX_CHECK (dynamic_cast<SxHamSolver*>(potential));
            dynamic_cast<SxHamSolver*>(potential)->initialGuess (&*table, true);
         }
      } else if (table->containsGroup("eamPot")) {

      } else if (table->containsGroup ("tbBasis")) {
         SxSymbolTable *tbInitialGuess = table->getGroup("tbInitialGuess");
         potential->execute (tbInitialGuess);
      } else
         SX_EXIT;
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   
   // --- enter main loop
   cout << endl;
   cout << SX_SEPARATOR;
   cout << "| Enter Main Loop\n";
   cout << SX_SEPARATOR;
   cout << endl;
   //cout << "TODO: SHOULD BE IN TRY BLOCK\n";
   SxSymbolTable *main = table->getGroup("main");
   SxSymbolTable *cmd  = NULL;
   for (cmd  = main->begin();
        cmd != NULL;
        cmd  = cmd->nextSibling())
   {
      if       (SxStructOpt::isRegistered (cmd))  {  
         SxStructOpt structOpt (potential, str);
         structOpt.execute (cmd);
      } else if (SxMolDyn::isRegistered (cmd))  {
         SxMolDyn    molDyn    (str, potential);
         molDyn.execute (cmd);
      } else if (SxFrozenPhonon::isRegistered (cmd))  {
         SxFrozenPhonon   frozenPhonon (str, potential);
         frozenPhonon.execute (cmd);
      } else if (potential->isRegistered (cmd))  {
         potential->execute (cmd);
      } else  {
         sxprintf ("Unknown command.\n");
         cmd->print ();
         SX_EXIT;
      }
   }

   SxMemMonitor::print ();
   
   // --- clean up
   if (potential)  delete potential;
   
   return 0;
}


#else /* SX_STANDALONE */
#include <SxCLI.h>
SxString getLogo ()
{
   SxString res =
   // --- created with JavE, http://www.jave.de/
   "\n\n\n\n\n"
   "                               ....................,.\n"
   "                             .JGp`````````````` `++dMMa,\n"
   "                          .JYjdMN.``````````````` i+dpdNm.\n"
   "                        .d9.JMVdNBY?`` ```````` `!!++?0dJW,\n"
   "                       .B'.MEOdMM` ```````` `` ```?!?;?+Nc?h,\n"
   "                     .J^ .MEuMMY5` ` ..uuY= ``.TMMMNSxc++z++?o\n"
   "                    .D  .MBdMBYYh,` ````````.qHMH9MMIO3++jg+&Jx\n"
   "                   d%` .MEzM8O.+? `` ..++!`?M8ttttOZZzI`jMSuVHJ;\n"
   "                  J^ `.M8OdNOz .M`````````   `..,?TMNOv.MMKOHNv`\n"
   "                 J! `.M8ltwNOdMH%``````` ``?!`...`?77? .MNMNmg2\n"
   "                J!``.JEOttlOZMMO.  ` ```..J&ggsO?++I.   `MMHY98\n"
   "               J!`` JNOOtttOUm:..... +vTMHY993.++1zOzdxMNmV;` .T,\n"
   "              `` ` .M6OtlltOvMNMNMNv` ``` `....?+OOtzdldU9zO.+,tX.`\n"
   "             ,``  .MBttltttttOdMM0O:` ``` ``````++tOv^ ````?1JAsd,\n"
   "            J``` .MEtlltltlltZwMMZtO.```` +??;;++1?!  ` ...,;zTMMd,\n"
   "           .^ ` .dMOtttttltlttwMMZtO:`....&+HHHBBB5.`TYY=`??lqgMNI4.\n"
   "          .%...?+vHSxOOtltttltwMMROdTYUY=!``.```` ..````7WQmgVHNS&.b\n"
   "         .KJ+?+;tMMadTNsltttltdNMEti. .Iz&++;;?;??;+..``..1VMMMNaJJJr\n"
   "         d6gMMD?zdNWMMSxHmylttMM9``!!``.Otltttlzzzzzzzzzz&&dgMNNJ.J+r\n"
   "        JNMM! .jOtMNOZTMNNWNmdMM!``` ...tlttltttltlltttOUTMMMN8+..,7r\n"
   "       .MMMS+.ttttzMNOOI?TMNTMM5 ````1OOZOtlttttttlttttOttttHMgH?7T3r\n"
   "       .WMMktltltltdMIv+;;+?T1?++: ``` tOtZUUUU9U9TMMMSxgAQgMNAx;++.u\n"
   "         WMNtOgmttttdI?+??` .JMNe+,````?tOtttltltttttZwNMNMMMItO7=TQg,\n"
   "         JMMMMMNt;?tOlvz+;.```MNM!```````?1lttttllttlOMNMM8Oltttzz+..r\n"
   "         .MB1MMEv! .!jzzzz! ` ,BY````````` ?+1tlttOOGdMMMmAAggeAyO+.,`\n"
   "          ! JM8!```` mztOt;````` ``````````.?;ztlttdMMMMNmgAOttwONm+J\n"
   "           .M3``````,M6tzz``.`` Jg,.``````` i+1ttttdMMMMSxggAsyttzqF\n"
   "          .M! ````` JBT=^````` ?N+J.,```````.+?zttzdMMMMNNHmgmdW9TH;\n"
   // in next line: splitting trigraph ??! -> |
   "         .M3?"/**/"?!``  .N2``````` .,HMNMF ````````?+ztttdNMMMNzlttltVVM;"
   "\n"
   "        .H6++!!``` ^` ``````` ...+J+JJ+..``````?;1OtOVHMMMMNxlttzzdF\n"
   "       ,=  ``````` ..........   ?17HMMMNm+ `````???zttlttVHMMNmOttdF\n"
   "    dY```...+ggmmWMWYYY7=````?7777TYYT=7zC``` ??;;?;1tttltttZTMMNwdF\n"
   "    j+J+QgHYY=!`..`` ````` ?7MNmo.. ..J++... ```?+;??+zttltttOzOZWMH,.\n"
   "    MMB=` ....?++! ```` .uXUHMHMMMMMMMMMMMMSxo. `` `!+;1zOtlltttttOOVWH..\n"
   "\n"
   "+--------------------------------------------------------------------------"
   "---\n"
   "|                            S / P H I / n X\n"
   "+--------------------------------------------------------------------------"
   "---\n"
   "|\n"
   "|                 The ab-initio based multiscale library\n"
   "|\n"
   "|       www:      https://sxrepo.mpie.de\n"
   "|       contact:  Christoph Freysoldt, freysoldt@mpie.de\n"
   "|\n"
   "+--------------------------------------------------------------------------"
   "---\n\n";
   return res;
}


int main (int argc, char **argv)
{
   // --- initialize the SPHInX libraries and print logo
   initSPHInXMath ();
   SX_ALLOC_CACHE;

   SxLoopMPI::init (argc, argv);

   SxParallelHierarchy ();

   bool error = false;
   {
      // --- parse command line
      SxCLI cli (argc, argv);

      SxString input = cli.option ("-i|--input", "file", 
                                   "SPHInX input file").toString("input.sx");
#  ifndef NDEBUG                                
      bool segFaults = cli.option ("--noexceptions", 
                                   "Handle exceptions").toBool ();
      if (segFaults)  SxException::causeSegFault ();
#  endif                                
      SxFFT::plannerCLI (cli);
      cli.finalize ();

      cout << getLogo ();

      if (SxUtil::getGlobalObj ().nProcs != 1)  {
         cout << "Number of processors used for openMP: "
              << SxUtil::getGlobalObj ().nProcs
              << endl;
      }
      
      SPHInX application (input);
      error = application.compute ();
      if (!error)  {
         cout << SX_SEPARATOR;
         cout << "| Program exited normally.\n";
         cout << SX_SEPARATOR;
      }
   }


#ifdef USE_LOOPMPI
   SxParallelHierarchy::write ("parallelHierarchy.sx.actual");
   SxLoopMPI::finalize ();  // LoopMPI
#endif
   return error;
}

#endif /* SX_STANDALONE */
