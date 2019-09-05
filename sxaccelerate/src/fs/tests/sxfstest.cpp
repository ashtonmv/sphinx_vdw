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

// --- Including header-files
#include <SxUtil.h>
#include <SxString.h>
#include <SxError.h>
#include <SxCLI.h>

#include <SxFS.h>
#include <SxFSTest.h>
#include <SxFSAction.h>
#include <SxFile.h>
#include <SxDir.h>

//------------------------------------------------------------------------------
// mv and cp_a
//------------------------------------------------------------------------------
   
enum ElemTypes {File = 0, Dir = 1, Nothing = 2, FileSymLink = 3,  
                DirSymLink = 4, EmptySymLink = 5, Self = 6};
enum Operation {Copy = 0, Move = 1};

void createElemOfType  (ElemTypes type, SxString name)
{
   if (type == File)  {
      TRY_N_PRINT_EX (SxFSAction::touch (name));
   } else if (type == Dir)  {
      TRY_N_PRINT_EX (SxFSAction::mkdir (name));
   } else if (type == Nothing)  {
      // empty
   } else if (type == FileSymLink)  {
      SxString target = name + "Target";
      TRY_N_PRINT_EX (SxFSAction::touch (target));
      TRY_N_PRINT_EX (SxFSAction::ln_sf (target, name));
   } else if (type == DirSymLink)  {
      SxString target = name + "Target";
      TRY_N_PRINT_EX (SxFSAction::mkdir (target));
      TRY_N_PRINT_EX (SxFSAction::ln_sf (target, name));
   } else if (type == EmptySymLink)  {
      SxString target = name + "Target";
      // Willingly omitting target creation
      TRY_N_PRINT_EX (SxFSAction::ln_sf (target, name));
   }      
}

bool checkPostcondition (ElemTypes type, SxString name)
{
   if (type == File)  {
      return !SxSymLink (name).exists () && SxFile (name).exists ();
   } else if (type == Dir)  {
      return !SxSymLink (name).exists () && SxDir (name).exists ();
   } else if (type == Nothing)  {
      return !SxSymLink (name).exists () 
                       && !SxDir (name).exists () 
                       && !SxFile (name).exists ();
   } else if (type == FileSymLink)  {
      return SxSymLink (name).exists () && SxSymLink (name).isFile ();
   } else if (type == DirSymLink)  {
      return SxSymLink (name).exists () && SxSymLink (name).isDir ();
   } else if (type == EmptySymLink)  {
      return SxSymLink (name).exists () 
                       && !SxSymLink (name).isFile () 
                       && !SxSymLink (name).isDir ();
   } else  {
      SX_EXIT;
   }
}

void printTestResults ()
{
   globTestSuite->printResults ();
}

void newTest (Operation op)
{
   // Checking whether a file system element with the name of the test 
   // directory already exists
   SxString testDir;
   if (op == Copy)  {
        testDir = "cpDir";
   } else  {
        testDir = "mvDir";
   }

   SX_CHECK (!SxDir (testDir).exists ());
   // Creating the test directory
   TRY_N_PRINT_EX (SxFSAction::mkdir (testDir));
   SxString srcDir ((SxFileInfo (testDir)/"srcDir").getAbsPath ());
   SxString destDir ((SxFileInfo (testDir)/"srcDir").getAbsPath ());
   TRY_N_PRINT_EX (SxFSAction::mkdir (srcDir));
   if (srcDir != destDir)  {
      TRY_N_PRINT_EX (SxFSAction::mkdir (destDir));
   }
   SxList<SxString> fsElem;
   fsElem << "File";
   fsElem << "Dir";
   fsElem << "Nothing";
   fsElem << "FileSymLink";
   fsElem << "DirSymLink";
   fsElem << "EmptySymLink";
   fsElem << "Self";

   // --- Creating some lists
   SxList<bool> correctResult;
   SxList<int> postSrcTypes, postDestTypes;

   // --- Reading in the correct results and postconditions from a file
   SxString referenceFile;
   if (op == Copy)  {
      referenceFile = "cp_aReferenceFile";
   } else  {
      referenceFile = "mvReferenceFile";
   }
   SxString fileContent;
   fileContent = SxString::read (referenceFile);
   // Fetching the lines
   SxList<SxString> lines = fileContent.tokenize ('\n');
   // Traversing all lines apart from the first one
   SxList<SxString>::Iterator itLines, itTokens;
   size_t iLine;
   SX_CHECK (lines.getSize () > 1);
   for (itLines = ++lines.begin (), iLine = 0; 
        itLines != lines.end (); 
        ++itLines, ++iLine)
   {
      cout << iLine << ": " << (*itLines) << endl;
      // Getting the blank-separated tokens
      SxList<SxString> tokens = (*itLines).tokenize (' ');
      // --- Running through the tokens
      size_t iToken;
      for (itTokens = tokens.begin (), iToken = 0; 
           itTokens != tokens.end (); 
           ++itTokens, ++iToken)
      {
         // Willingly ignoring the first two tokens
         // Appending the tokens to the list of correct results and 
         // postconditions
         if (iToken == 2)  {
            correctResult << (!(*itTokens).toInt ());
         } else if (iToken == 3) {
            postSrcTypes << ((ElemTypes) (*itTokens).toInt ());         
         } else if (iToken == 4) {
            postDestTypes << ((ElemTypes) (*itTokens).toInt ());         
         }
      }
   }

   // --- Performing all tests
   SX_INIT_TESTS ();
   SxList<SxString>::Iterator itSrc, itDest;
   SxList<bool>::Iterator itCorrectResult;
   size_t id, srcId (0), destId (6*7);
   int iSrc, iDest;
   SxList<int>::Iterator itPostSrc, itPostDest;
   for (itSrc = fsElem.begin (), itCorrectResult = correctResult.begin (), 
        itPostSrc = postSrcTypes.begin (), 
        itPostDest = postDestTypes.begin (), 
        id = 0, iSrc = 0; 
        itSrc != fsElem.end (); 
        ++itSrc, ++iSrc)
   {
      //if (id <2)  {//TODO Continue testing iteratively 
      if (iSrc != Self)  {
         for (itDest = fsElem.begin (), iDest = 0; 
              itDest != fsElem.end (); 
              ++itDest, ++iDest)
         {
            SxString const &curSrc = (*itSrc);
            SxString const &curDest = (*itDest);
            SxString srcName;
            srcName = (SxFileInfo(srcDir) / 
                       (curSrc + SxString((int)(id + srcId)) + 
                        ((iSrc >= 3 && iSrc <= 5)? ".lnk" : "") )
                       ).getAbsPath ();
            createElemOfType ((ElemTypes)iSrc, srcName);

            SxString destName;
            if (iDest != Self)  {
               destName = (SxFileInfo(destDir) / 
                        (curDest + SxString((int)(id + destId)) +
                         ((iDest >= 3 && iDest <= 5)? ".lnk" : "") 
                        )).getAbsPath ();
               createElemOfType ((ElemTypes)iDest, destName);
            } else  {
               destName = srcName;
            }
            
            // --- Performing the current test
            SX_CREATE_TEST ("'" + curSrc + "' to '" + curDest + "'");
            //SX_CHECK_PRE ();//TODO Implement
            //SX_CHECK_PRE ();//TODO Implement
            if (op == Copy)  {
               SX_MAKE_TEST (SxFSAction::cp (srcName, destName), 
                             (*itCorrectResult));
            } else  {
               SX_MAKE_TEST (SxFSAction::mv (srcName, destName), 
                             (*itCorrectResult));            
            }
            SX_CHECK_POST (checkPostcondition ((ElemTypes)(*itPostSrc), 
                                               srcName));
            SX_CHECK_POST (checkPostcondition ((ElemTypes)(*itPostDest), 
                                               destName));
            SX_QUIT_TEST ();
            
            ++itCorrectResult;
            ++itPostSrc;
            ++itPostDest;
            ++id;
         }
      }
      //}
   }
   // Printing the evaluation of the results
   printTestResults ();
}

//------------------------------------------------------------------------------
// main
//------------------------------------------------------------------------------

/** \example fstest.cpp

  In this example the usage of the SPHInX file info classes is demonstrated.

  \brief  Basic usage of SxFileInfo, SxFile, SxDir and SxSymLink
  \author Thomas Uchdorf
  */
int main (int argc, char **argv)
{
   // --- Using the command line interface to retrieve the program arguments
   SxCLI cli (argc, argv);
   bool enableCp_a = cli.option ("--cp_a",
                               "Enables \"cp_a\"-tests"
                              ).toBool ();
   bool enableMv = cli.option ("--mv",
                               "Enables \"mv\"-tests"
                              ).toBool ();
   cli.finalize (); 

   if (enableCp_a) 
        newTest (Copy);
   if (enableMv)
        newTest (Move);
   return 0;
}
