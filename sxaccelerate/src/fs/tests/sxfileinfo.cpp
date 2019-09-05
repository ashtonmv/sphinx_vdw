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
// mv
//------------------------------------------------------------------------------
void mvDirToExisting (const SxString &srcRootDir,
                      const SxString &destRootDir)
{

   SX_INIT_TESTS ();

   // --- Moving a directory to another directory
   SX_CREATE_TEST ("");
   // Preconditions:
   // srcRootDir + "/dir1" exists and is a directory
   // destRootDir + "/dir2" exists and is a directory
   SX_CHECK_PRE (SxDir (srcRootDir + "/dir1").exists ());
   SX_CHECK_PRE (SxDir (destRootDir + "/dir2").exists ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/dir1",
                                 destRootDir + "/dir2"), true);
   // Postconditions:
   // srcRootDir + "/dir1" does not exist
   // destRootDir + "/dir2/dir1" exists and is a directory that corresponds
   // with the by now deleted directory srcRootDir + "/dir1"
   SX_CHECK_POST (!SxDir (srcRootDir + "/dir1").exists ());
   SX_CHECK_POST (SxDir (destRootDir + "/dir2/dir1").exists ());
   SX_QUIT_TEST ();

   // --- Moving a directory to a symbolic link which points to another
   //     directory
   SX_CREATE_TEST ("");
   // Preconditions:
   // srcRootDir + "/dir3" exists and is a directory
   // destRootDir + "/symLinkToDir4" exists and is a symbolic link to a
   // directory
   SX_CHECK_PRE (SxDir (srcRootDir + "/dir3").exists ());
   SX_CHECK_PRE (SxSymLink (destRootDir + "/symLinkToDir4").exists () &&
                 SxSymLink (destRootDir + "/symLinkToDir4").isValid () &&
                 SxSymLink (destRootDir + "/symLinkToDir4").isDir ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/dir3",
                                 destRootDir + "/symLinkToDir4"),
                 true);
   // Postconditions:
   // srcRootDir + "/dir3" does not exist
   // destRootDir + "/symLinkToDir4/dir3" exists and is a directory that
   // corresponds with the by now deleted directory srcRootDir + "/dir3"
   SX_CHECK_POST (!SxDir (srcRootDir + "/dir3").exists ());
   SX_CHECK_POST (SxDir (destRootDir + "/symLinkToDir4/dir3").exists ());
   SX_QUIT_TEST ();

   // --- Moving a directory to a file
   SX_CREATE_TEST ("");
   // Preconditions:
   // srcRootDir + "/dir5" exists and is a directory
   // destRootDir + "/file1" exists and is a file
   SX_CHECK_PRE (SxDir (srcRootDir + "/dir5").exists ());
   SX_CHECK_PRE (SxFile (destRootDir + "/file1").exists ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/dir5",
                                 destRootDir + "/file1"),
                 false);
   // Postconditions:
   // operation is not permitted
   SX_CHECK_POST (SxDir (srcRootDir + "/dir5").exists ());
   SX_CHECK_POST (SxFile (destRootDir + "/file1").exists ());
   SX_QUIT_TEST ();

   // --- Moving a directory to a symbolic link which points to a file
   SX_CREATE_TEST ("");
   // Preconditions:
   // srcRootDir + "/dir6" exists and is a directory
   // destRootDir + "/symLinkToFile2" exists and is a symbolic link to a file
   SX_CHECK_PRE (SxDir (srcRootDir + "/dir6").exists ());
   SX_CHECK_PRE (SxSymLink (destRootDir + "/symLinkToFile2").exists () &&
                 SxSymLink (destRootDir + "/symLinkToFile2").isValid () &&
                 SxSymLink (destRootDir + "/symLinkToFile2").isFile ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/dir6",
                                 destRootDir + "/symLinkToFile2"),
                 false);
   // Postconditions:
   // operation is not permitted
   SX_CHECK_POST (SxDir (srcRootDir + "/dir6").exists ());
   SX_CHECK_POST (SxSymLink (destRootDir + "/symLinkToFile2").exists () &&
                  SxSymLink (destRootDir + "/symLinkToFile2").isValid () &&
                  SxSymLink (destRootDir + "/symLinkToFile2").isFile ());
   SX_QUIT_TEST ();

   // --- Moving a directory to itself
   SX_CREATE_TEST ("");
   // Preconditions:
   // srcRootDir + "/dir7" exists and is a directory
   SX_CHECK_PRE (SxDir (srcRootDir + "/dir7").exists ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/dir7",
                                 destRootDir + "/dir7"), false);
   // Postconditions:
   // operation is not permitted
   SX_CHECK_POST (SxDir (srcRootDir + "/dir7").exists ());
   SX_QUIT_TEST ();

   // --- Moving a directory to an empty symbolic link
   SX_CREATE_TEST ("");
   // Preconditions:
   // srcRootDir + "/dir19" exists and is a directory
   // destRootDir + "/symLinkToNothing1" is a symbolic link that points to
   // nowhere
   SX_CHECK_PRE (SxDir (srcRootDir + "/dir19").exists ());
   SX_CHECK_PRE (SxSymLink (destRootDir + "/symLinkToNothing1").exists () &&
                 !SxSymLink (destRootDir + "/symLinkToNothing1").isFile () &&
                 !SxSymLink (destRootDir + "/symLinkToNothing1").isDir ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/dir19",
                                 destRootDir + "/symLinkToNothing1"),
                 false);
   // Postconditions:
   // operation is not permitted
   SX_CHECK_POST (SxDir (srcRootDir + "/dir19").exists ());
   SX_CHECK_POST (SxSymLink (destRootDir + "/symLinkToNothing1").exists () &&
                  !SxSymLink (destRootDir + "/symLinkToNothing1").isFile () &&
                  !SxSymLink (destRootDir + "/symLinkToNothing1").isDir ());
   SX_QUIT_TEST ();
}

void mvFileToExisting (const SxString &srcRootDir,
                       const SxString &destRootDir)
{

   SX_INIT_TESTS ();

   // --- Moving a file to a directory
   SX_CREATE_TEST ("");
   // Preconditions:
   // srcRootDir + "/file3" exists and is a file
   // destRootDir + "/dir8" exists and is directory
   SX_CHECK_PRE (SxFile (srcRootDir + "/file3").exists ());
   SX_CHECK_PRE (SxDir (destRootDir + "/dir8").exists ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/file3",
                                 destRootDir + "/dir8"), true);
   // Postconditions:
   // srcRootDir + "/file3" does not exist
   // destRootDir + "/dir8/file3" exists and is a file that corresponds with
   // the by now deleted file srcRootDir + "/file3"
   SX_CHECK_POST (!SxFile (srcRootDir + "/file3").exists ());
   SX_CHECK_POST (SxFile (destRootDir + "/dir8/file3").exists ());
   SX_QUIT_TEST ();

   // --- Moving a file to a symbolic link that points to a directory
   SX_CREATE_TEST ("");
   // Preconditions:
   // srcRootDir + "/file4" exists and is a file
   // destRootDir + "/symLinkToDir9" is a valid symbolic link to a directory
   SX_CHECK_PRE (SxFile (srcRootDir + "/file4").exists ());
   SX_CHECK_PRE (SxSymLink (destRootDir + "/symLinkToDir9").exists () &&
                  SxSymLink (destRootDir + "/symLinkToDir9").isValid () &&
                  SxSymLink (destRootDir + "/symLinkToDir9").isDir ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/file4",
                                 destRootDir + "/symLinkToDir9"),
                 true);
   // Postconditions:
   // srcRootDir + "/file4" does not exist
   // destRootDir + "/symLinkToDir9/file4" exists and is a file that
   // corresponds with the by now deleted file srcRootDir + "/file4"
   SX_CHECK_POST (!SxFile (srcRootDir + "/file4").exists ());
   SX_CHECK_POST (SxFile (destRootDir + "/symLinkToDir9/file4").exists ());
   SX_QUIT_TEST ();

   // --- Moving a file to another file
   SX_CREATE_TEST ("");
   // Preconditions:
   // srcRootDir + "/file5" exists and is a file
   // destRootDir + "/file6" exists and is a file
   SX_CHECK_PRE (SxFile (srcRootDir + "/file5").exists ());
   SX_CHECK_PRE (SxFile (destRootDir + "/file6").exists ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/file5",
                                 destRootDir + "/file6"),
                 true);
   // Postconditions:
   // srcRootDir + "/file5" does not exist
   // destRootDir + "/file6" exists and is a file that corresponds with the
   // by now deleted file srcRootDir + "/file5"
   SX_CHECK_POST (!SxFile (srcRootDir + "/file5").exists ());
   SX_CHECK_POST (SxFile (destRootDir + "/file6").exists ());
   SX_QUIT_TEST ();

   // --- Moving a file to a symbolic link that points to a file
   SX_CREATE_TEST ("");
   // Preconditions:
   // srcRootDir + "/file7" exists and is a file
   // destRootDir + "/symLinkToFile8" is a valid symbolic link to a file
   SX_CHECK_PRE (SxFile (srcRootDir + "/file7").exists ());
   SX_CHECK_PRE (SxSymLink (destRootDir + "/symLinkToFile8").exists () &&
                  SxSymLink (destRootDir + "/symLinkToFile8").isValid () &&
                  SxSymLink (destRootDir + "/symLinkToFile8").isFile ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/file7",
                                 destRootDir + "/symLinkToFile8"),
                 true);
   // Postconditions:
   // srcRootDir + "/file7" does not exist
   // destRootDir + "/symLinkToFile8" exists and is a file that corresponds
   // with the by now deleted file srcRootDir + "/file7"
   SX_CHECK_POST (!SxFile (srcRootDir + "/file7").exists ());
   SX_CHECK_POST (!SxFileInfo (destRootDir + "/symLinkToFile8").isSymLink () &&
                  SxFile (destRootDir + "/symLinkToFile8").exists ());
   SX_QUIT_TEST ();

   // --- Moving a file to itself
   SX_CREATE_TEST ("");
   // Preconditions:
   // srcRootDir + "/file9" exists and is a file
   SX_CHECK_PRE (SxFile (srcRootDir + "/file9").exists ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/file9",
                                 destRootDir + "/file9"),
                 false);
   // Postconditions:
   // operation is not permitted
   SX_CHECK_POST (SxFile (srcRootDir + "/file9").exists ());
   SX_QUIT_TEST ();

   // --- Moving a file to an empty symbolic link
   //TODO Investigate why this does not work when moving from a directory on
   //     scratch to raid!
   SX_CREATE_TEST ("");
   SX_CHECK (SxFile (srcRootDir + "/file19").exists ());
   SX_CHECK (SxSymLink (destRootDir + "/symLinkToNothing2").exists () &&
                 !SxSymLink (destRootDir + "/symLinkToNothing2").isFile () &&
                 !SxSymLink (destRootDir + "/symLinkToNothing2").isDir ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/file19",
                                 destRootDir + "/symLinkToNothing2"), true);
   // Postconditions:
   // srcRootDir + "/file19" does not exist
   // destRootDir + "/symLinkToNothing2" exists and is a file that corresponds
   // with the by now removed file srcRootDir + "/file19"
   SX_CHECK (!SxFile (srcRootDir + "/file19").exists ());
   SX_CHECK (!SxSymLink (destRootDir + "/symLinkToNothing2").exists () &&
                  SxFile (destRootDir + "/symLinkToNothing2").exists ());
   SX_QUIT_TEST ();
}

void mvDirSymLinkToExisting (const SxString &srcRootDir,
                             const SxString &destRootDir)
{

   SX_INIT_TESTS ();

   // --- Moving a symbolic link that addresses a directory to another one
   SX_CREATE_TEST ("");
   // Preconditions:
   // srcRootDir + "/symLinkToDir10" is a valid symbolic link to a directory
   // destRootDir + "/symLinkToDir11" is a valid symbolic link to a directory
   SX_CHECK_PRE (SxSymLink (srcRootDir + "/symLinkToDir10").exists () &&
                 SxSymLink (srcRootDir + "/symLinkToDir10").isValid () &&
                 SxSymLink (srcRootDir + "/symLinkToDir10").isDir ());
   SX_CHECK_PRE (SxSymLink (destRootDir + "/symLinkToDir11").exists () &&
                 SxSymLink (destRootDir + "/symLinkToDir11").isValid () &&
                 SxSymLink (destRootDir + "/symLinkToDir11").isDir ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/symLinkToDir10",
                                 destRootDir + "/symLinkToDir11"),
                 true);
   // Postconditions:
   // srcRootDir + "/symLinkToDir10" does not exist
   // destRootDir + "/symLinkToDir11" is a valid symbolic link to a directory
   // and corresponds with the by now eliminated srcRootDir + "/symLinkToDir10"
   SX_CHECK_POST (!SxSymLink (srcRootDir + "/symLinkToDir10").exists ());
   SX_CHECK_POST (SxSymLink (destRootDir + "/symLinkToDir11").exists () &&
                  SxSymLink (destRootDir + "/symLinkToDir11").isValid () &&
                  SxSymLink (destRootDir + "/symLinkToDir11").isDir ());
   SX_QUIT_TEST ();

   // --- Moving a symbolic link that addresses a directory to a directory
   SX_CREATE_TEST ("");
   // Preconditions:
   // srcRootDir + "/symLinkToDir12" is a valid symbolic link to a directory
   // destRootDir + "/dir13" exists and is a directory
   SxString oldTarget = SxSymLink (srcRootDir + "/symLinkToDir12").getTarget ();
   SX_CHECK_PRE (SxSymLink (srcRootDir + "/symLinkToDir12").exists () &&
                 SxSymLink (srcRootDir + "/symLinkToDir12").isValid () &&
                 SxSymLink (srcRootDir + "/symLinkToDir12").isDir ());
   SX_CHECK_PRE (SxDir (destRootDir + "/dir13").exists ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/symLinkToDir12",
                                 destRootDir + "/dir13"),
                 true);
   // Postconditions:
   // srcRootDir + "/symLinkToDir12" does not exist
   // destRootDir + "/dir13/symLinkToDir12" is a symbolic to a directory that
   // corresponds with the by now removed symbolic link
   // srcRootDir + "/symLinkToDir12"
   SX_CHECK_POST (!SxSymLink (srcRootDir + "/symLinkToDir12").exists ());
   SX_CHECK_POST (SxSymLink (destRootDir + "/dir13/symLinkToDir12").exists () &&
                  SxSymLink (destRootDir + "/dir13/symLinkToDir12").getTarget ()
                  == oldTarget);
   SX_QUIT_TEST ();

   // --- Moving a symbolic link that addresses a directory to a file
   SX_CREATE_TEST ("");
   // Preconditions:
   // srcRootDir + "/symLinkToDir14" is a valid symbolic link to a directory
   // destRootDir + "/file10" exists and is a file
   SX_CHECK_PRE (SxSymLink (srcRootDir + "/symLinkToDir14").exists () &&
                 SxSymLink (srcRootDir + "/symLinkToDir14").isValid () &&
                 SxSymLink (srcRootDir + "/symLinkToDir14").isDir ());
   SX_CHECK_PRE (SxFile (destRootDir + "/file10").exists ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/symLinkToDir14",
                                 destRootDir + "/file10"),
                 true);
   // Postconditions:
   // srcRootDir + "/symLinkToDir14" does not exist
   // destRootDir + "/file10" is a valid symbolic link to a directory that
   // corresponds with the by now removed symbolic link
   // srcRootDir + "/symLinkToDir14"
   SX_CHECK_POST (!SxSymLink (srcRootDir + "/symLinkToDir14").exists ());
   SX_CHECK_POST (SxSymLink (destRootDir + "/file10").exists () &&
                  SxSymLink (destRootDir + "/file10").isValid () &&
                  SxSymLink (destRootDir + "/file10").isDir ());
   SX_QUIT_TEST ();

   // --- Moving a symbolic link that addresses a directory to a symbolic link
   //     that points to a file
   SX_CREATE_TEST ("");
   // Preconditions:
   // srcRootDir + "/symLinkToDir15" is a valid symbolic link to a directory
   // destRootDir + "/symLinkToFile11" is a valid symbolic link to a file
   SX_CHECK_PRE (SxSymLink (srcRootDir + "/symLinkToDir15").exists () &&
                 SxSymLink (srcRootDir + "/symLinkToDir15").isValid () &&
                 SxSymLink (srcRootDir + "/symLinkToDir15").isDir ());
   SX_CHECK_PRE (SxSymLink (destRootDir + "/symLinkToFile11").exists () &&
                 SxSymLink (destRootDir + "/symLinkToFile11").isValid () &&
                 SxSymLink (destRootDir + "/symLinkToFile11").isFile ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/symLinkToDir15",
                                 destRootDir + "/symLinkToFile11"),
                 true);
   // Postconditions:
   // srcRootDir + "/symLinkToDir15" does not exist
   // destRootDir + "/symLinkToFile11" is a valid symbolic link to a directory
   // that corresponds with the by now removed symbolic link
   // srcRootDir + "/symLinkToDir15"
   SX_CHECK_POST (!SxSymLink (srcRootDir + "/symLinkToDir15").exists ());
   SX_CHECK_POST (SxSymLink (destRootDir + "/symLinkToFile11").exists () &&
                  SxSymLink (destRootDir + "/symLinkToFile11").isValid () &&
                  SxSymLink (destRootDir + "/symLinkToFile11").isDir ());
   SX_QUIT_TEST ();

   // --- Moving a symbolic link that addresses a directory to itself
   SX_CREATE_TEST ("");
   // Preconditions:
   // srcRootDir + "/symLinkToDir16" is a valid symbolic link to a directory
   SX_CHECK_PRE (SxSymLink (srcRootDir + "/symLinkToDir16").exists () &&
                 SxSymLink (srcRootDir + "/symLinkToDir16").isValid () &&
                 SxSymLink (srcRootDir + "/symLinkToDir16").isDir ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/symLinkToDir16",
                                 destRootDir + "/symLinkToDir16"),
                 false);
   // Postconditions:
   // operation is not permitted
   SX_CHECK_POST (SxSymLink (srcRootDir + "/symLinkToDir16").exists () &&
                  SxSymLink (srcRootDir + "/symLinkToDir16").isValid () &&
                  SxSymLink (srcRootDir + "/symLinkToDir16").isDir ());
   SX_QUIT_TEST ();

   // --- Moving a symbolic link that addresses a directory to an empty
   //     symbolic link
   SX_CREATE_TEST ("");
   // Preconditions:
   // srcRootDir + "/symLinkToDir20" is a valid symbolic link that points to a
   // directory
   // destRootDir + "/symLinkToNothing3" is a symbolic link that points to
   // nowhere
   oldTarget = SxSymLink (srcRootDir + "/symLinkToDir20").getTarget ();
   SX_CHECK_PRE (SxSymLink (srcRootDir + "/symLinkToDir20").exists () &&
                 SxSymLink (srcRootDir + "/symLinkToDir20").isValid () &&
                 SxSymLink (srcRootDir + "/symLinkToDir20").isDir ());
   SX_CHECK_PRE (SxSymLink (destRootDir + "/symLinkToNothing3").exists () &&
                 !SxSymLink (destRootDir + "/symLinkToNothing3").isFile () &&
                 !SxSymLink (destRootDir + "/symLinkToNothing3").isDir ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/symLinkToDir20",
                                 destRootDir + "/symLinkToNothing3"), true);
   // Postconditions:
   // srcRootDir + "/symLinkToDir20" does not exist
   // destRootDir + "/symLinkToNothing" is a valid symbolic link that points
   // to a directory that corresponds with the target of the by now deleted
   // symbolic link srcRootDir + "/symLinkToDir20"
   SX_CHECK_POST (!SxSymLink (srcRootDir + "/symLinkToDir20").exists ());
   SX_CHECK_POST (SxSymLink (destRootDir + "/symLinkToNothing3").exists () &&
                  SxSymLink (destRootDir + "/symLinkToNothing3").isValid () &&
                  SxSymLink (destRootDir + "/symLinkToNothing3").isDir () &&
                  SxSymLink (destRootDir + "/symLinkToNothing3").getTarget ()
                  == oldTarget);
   SX_QUIT_TEST ();
}

void mvFileSymLinkToExisting (const SxString &srcRootDir,
                              const SxString &destRootDir)
{

   SX_INIT_TESTS ();

   // --- Moving a symbolic link that addresses a file to itself
   SX_CREATE_TEST ("");
   // Preconditions:
   // srcRootDir + "/symLinkToFile12" is a valid symbolic link that points to a
   // file
   SX_CHECK_PRE (SxSymLink (srcRootDir + "/symLinkToFile12").exists () &&
                 SxSymLink (srcRootDir + "/symLinkToFile12").isValid () &&
                 SxSymLink (srcRootDir + "/symLinkToFile12").isFile ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/symLinkToFile12",
                                 destRootDir + "/symLinkToFile12"),
                 false);
   // Postconditions:
   // operation is not permitted
   SX_CHECK_POST (SxSymLink (srcRootDir + "/symLinkToFile12").exists () &&
                  SxSymLink (srcRootDir + "/symLinkToFile12").isValid () &&
                  SxSymLink (srcRootDir + "/symLinkToFile12").isFile ());
   SX_QUIT_TEST ();

   // --- Moving a symbolic link that addresses a file to a symbolic link that
   //     addresses a directory
   SX_CREATE_TEST ("");
   // Preconditions:
   // srcRootDir + "/symLinkToFile13" is a valid symbolic link that points to a
   // file
   // destRootDir + "/symLinkToDir17" is a valid symbolic link that points to a
   // directory
   SxString oldTarget = SxSymLink (srcRootDir + "/symLinkToFile13").getTarget ();
   SX_CHECK_PRE (SxSymLink (srcRootDir + "/symLinkToFile13").exists () &&
                 SxSymLink (srcRootDir + "/symLinkToFile13").isValid () &&
                 SxSymLink (srcRootDir + "/symLinkToFile13").isFile ());
   SX_CHECK_PRE (SxSymLink (destRootDir + "/symLinkToDir17").exists () &&
                 SxSymLink (destRootDir + "/symLinkToDir17").isValid () &&
                 SxSymLink (destRootDir + "/symLinkToDir17").isDir ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/symLinkToFile13",
                                 destRootDir + "/symLinkToDir17"),
                 true);
   // Postconditions:
   // srcRootDir + "/symLinkToFile13" does not exist
   // destRootDir + "/symLinkToDir17/symLinkToFile13" is a symbolic link that
   // corresponds with the by now deleted symbolic link
   // srcRootDir + "/symLinkToFile13"
   SX_CHECK_POST (!SxSymLink (srcRootDir + "/symLinkToFile13").exists ());
   SX_CHECK_POST (SxSymLink (destRootDir + "/symLinkToDir17/symLinkToFile13").exists () &&
                 SxSymLink (destRootDir + "/symLinkToDir17/symLinkToFile13").getTarget ()
                 == oldTarget);
   SX_QUIT_TEST ();

   // --- Moving a symbolic link that addresses a file to a directory
   SX_CREATE_TEST ("");
   // Preconditions:
   // srcRootDir + "/symLinkToFile14" is a valid symbolic link that points to a
   // file
   // destRootDir + "/dir18" exists and is a directory
   oldTarget = SxSymLink (srcRootDir + "/symLinkToFile14").getTarget ();
   SX_CHECK_PRE (SxSymLink (srcRootDir + "/symLinkToFile14").exists () &&
                 SxSymLink (srcRootDir + "/symLinkToFile14").isValid () &&
                 SxSymLink (srcRootDir + "/symLinkToFile14").isFile ());
   SX_CHECK_PRE (SxDir (destRootDir + "/dir18").exists ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/symLinkToFile14",
                                 destRootDir + "/dir18"),
                 true);
   // Postconditions:
   // srcRootDir + "/symLinkToFile14" does not exist
   // destRootDir + "/dir18/symLinkToFile14" is a symbolic link that
   // corresponds with the by now deleted symbolic link srcRootDir + "/symLinkToFile14"
   SX_CHECK_POST (!SxSymLink (srcRootDir + "/symLinkToFile14").exists ());
   SX_CHECK_POST (SxSymLink (destRootDir + "/dir18/symLinkToFile14").exists () &&
                  SxSymLink (destRootDir + "/dir18/symLinkToFile14").getTarget ()
                  == oldTarget);
   SX_QUIT_TEST ();

   // --- Moving a symbolic link that addresses a file to a file
   SX_CREATE_TEST ("");
   // Preconditions:
   // srcRootDir + "/symLinkToFile15" is a valid symbolic link that points to a
   // file
   // destRootDir + "/file16" exists and is a file
   SX_CHECK_PRE (SxSymLink (srcRootDir + "/symLinkToFile15").exists () &&
                 SxSymLink (srcRootDir + "/symLinkToFile15").isValid () &&
                 SxSymLink (srcRootDir + "/symLinkToFile15").isFile ());
   SX_CHECK_PRE (SxFile (destRootDir + "/file16").exists ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/symLinkToFile15",
                                 destRootDir + "/file16"),
                 true);
   // Postconditions:
   // srcRootDir + "/symLinkToFile15" does not exist
   // destRootDir + "/file16" is a valid symbolic link that points to a
   // file and corresponds with the by now deleted symbolic link
   // srcRootDir + "/symLinkToFile15"
   SX_CHECK_POST (!SxSymLink (srcRootDir + "/symLinkToFile15").exists ());
   SX_CHECK_POST (SxSymLink (destRootDir + "/file16").exists () &&
                  SxSymLink (destRootDir + "/file16").isValid () &&
                  SxSymLink (destRootDir + "/file16").isFile ());
   SX_QUIT_TEST ();

   // --- Moving a symbolic link that addresses a file to another one
   SX_CREATE_TEST ("");
   // Preconditions:
   // srcRootDir + "/symLinkToFile17" is a valid symbolic link that points to a
   // file
   // destRootDir + "/symLinkToFile18" is a valid symbolic link that points
   // to a file
   SX_CHECK_PRE (SxSymLink (srcRootDir + "/symLinkToFile17").exists () &&
                 SxSymLink (srcRootDir + "/symLinkToFile17").isValid () &&
                 SxSymLink (srcRootDir + "/symLinkToFile17").isFile ());
   SX_CHECK_PRE (SxSymLink (destRootDir + "/symLinkToFile18").exists () &&
                 SxSymLink (destRootDir + "/symLinkToFile18").isValid () &&
                 SxSymLink (destRootDir + "/symLinkToFile18").isFile ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/symLinkToFile17",
                                 destRootDir + "/symLinkToFile18"),
                 true);
   // Postconditions:
   // srcRootDir + "/symLinkToFile17" does not exist
   // destRootDir + "/symLinkToFile18" is a valid symbolic link that points to
   // the same file srcRootDir + "/symLinkToFile17" addressed before it was
   // deleted
   SX_CHECK_POST (!SxSymLink (srcRootDir + "/symLinkToFile17").exists ());
   SX_CHECK_POST (SxSymLink (destRootDir + "/symLinkToFile18").exists () &&
                  SxSymLink (destRootDir + "/symLinkToFile18").isValid () &&
                  SxSymLink (destRootDir + "/symLinkToFile18").isFile ());
   SX_QUIT_TEST ();

   // --- Moving a symbolic link that addresses a file to an empty symbolic
   //     link
   SX_CREATE_TEST ("");
   oldTarget = SxSymLink (srcRootDir + "/symLinkToFile20").getTarget ();
   // Preconditions:
   // srcRootDir + "/symLinkToFile20" is a valid symbolic link that addresses a
   // file
   // destRootDir + "/symLinkToNothing4" is a symbolic link that points to
   // nowhere
   SX_CHECK_PRE (SxSymLink (srcRootDir + "/symLinkToFile20").exists () &&
                 SxSymLink (srcRootDir + "/symLinkToFile20").isValid () &&
                 SxSymLink (srcRootDir + "/symLinkToFile20").isFile ());
   SX_CHECK_PRE (SxSymLink (destRootDir + "/symLinkToNothing4").exists () &&
                 !SxSymLink (destRootDir + "/symLinkToNothing4").isFile () &&
                 !SxSymLink (destRootDir + "/symLinkToNothing4").isDir ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/symLinkToFile20",
                                 destRootDir + "/symLinkToNothing4"),
                 true);
   // Postconditions:
   // srcRootDir + "/symLinkToFile20" does not exist
   // destRootDir + "/symLinkToNothing4" is a valid symbolic link that
   // addresses a file which represents the target of the by now deleted
   // symbolic link srcRootDir + "/symLinkToFile20"
   SX_CHECK_POST (!SxSymLink (srcRootDir + "/symLinkToFile20").exists ());
   SX_CHECK_POST (SxSymLink (destRootDir + "/symLinkToNothing4").exists () &&
                  SxSymLink (destRootDir + "/symLinkToNothing4").isValid () &&
                  SxSymLink (destRootDir + "/symLinkToNothing4").isFile () &&
                  SxSymLink (destRootDir + "/symLinkToNothing4").getTarget () ==
                  oldTarget);
   SX_QUIT_TEST ();
}

void mvEmptySymLinkToExisting (const SxString &srcRootDir,
                               const SxString &destRootDir)
{

   SX_INIT_TESTS ();

   // --- Moving a symbolic link that addresses nothing to itself
   SX_CREATE_TEST ("");
   // Preconditions:
   // srcRootDir + "/symLinkToNothing5" is a symbolic link that points to nowhere
   SX_CHECK_PRE (SxSymLink (srcRootDir + "/symLinkToNothing5").exists () &&
                 !SxSymLink (srcRootDir + "/symLinkToNothing5").isFile () &&
                 !SxSymLink (srcRootDir + "/symLinkToNothing5").isDir ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/symLinkToNothing5",
                                 destRootDir + "/symLinkToNothing5"),
                 false);
   // Postconditions:
   // operation is not permitted
   SX_CHECK_POST (SxSymLink (srcRootDir + "/symLinkToNothing5").exists () &&
                  !SxSymLink (srcRootDir + "/symLinkToNothing5").isFile () &&
                  !SxSymLink (srcRootDir + "/symLinkToNothing5").isDir ());
   SX_QUIT_TEST ();

   // --- Moving a symbolic link that addresses nothing to a file
   SX_CREATE_TEST ("");
   SxString oldTarget = SxSymLink (srcRootDir + "/symLinkToNothing6").getTarget ();
   // Preconditions:
   // srcRootDir + "/symLinkToNothing6" is a symbolic link that points to
   // nowhere
   // destRootDir + "/file21" exists and is a file
   SX_CHECK_PRE (SxSymLink (srcRootDir + "/symLinkToNothing6").exists () &&
                 !SxSymLink (srcRootDir + "/symLinkToNothing6").isFile () &&
                 !SxSymLink (srcRootDir + "/symLinkToNothing6").isDir ());
   SX_CHECK_PRE (SxFile (destRootDir + "/file21").exists ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/symLinkToNothing6",
                                 destRootDir + "/file21"),
                 true);
   // Postconditions:
   // srcRootDir + "/symLinkToNothing6" does not exist
   // destRootDir + "/file21" is a valid symbolic link that points to the same
   // target as the by now deleted symbolic link srcRootDir + "/symLinkToNothing6"
   SX_CHECK_POST (!SxSymLink (srcRootDir + "/symLinkToNothing6").exists ());
   SX_CHECK_POST (SxSymLink (destRootDir + "/file21").exists () &&
                  !SxSymLink (destRootDir + "/file21").isFile () &&
                  !SxSymLink (destRootDir + "/file21").isDir () &&
                  SxSymLink (destRootDir + "/file21").getTarget () == oldTarget);
   SX_QUIT_TEST ();

   // --- Moving a symbolic link that addresses nothing to a directory
   SX_CREATE_TEST ("");
   oldTarget = SxSymLink (srcRootDir + "/symLinkToNothing7").getTarget ();
   // Preconditions:
   // srcRootDir + "/symLinkToNothing7" is a symbolic link that points to
   // nowhere
   // destRootDir + "/dir21" exists and is a file
   SX_CHECK_PRE (SxSymLink (srcRootDir + "/symLinkToNothing7").exists () &&
                 !SxSymLink (srcRootDir + "/symLinkToNothing7").isFile () &&
                 !SxSymLink (srcRootDir + "/symLinkToNothing7").isDir ());
   SX_CHECK_PRE (SxDir (destRootDir + "/dir21").exists ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/symLinkToNothing7",
                                 destRootDir + "/dir21"),
                 true);
   // Postconditions:
   // srcRootDir + "/symLinkToNothing7" does not exist
   // destRootDir + "/dir21/symLinkToNothing7" is a symbolic link that points
   // to the same target as the by now deleted symbolic link
   // srcRootDir + "/symLinkToNothing7"
   SX_CHECK_POST (!SxSymLink (srcRootDir + "/symLinkToNothing7").exists ());
   SX_CHECK_POST (SxSymLink (destRootDir + "/dir21/symLinkToNothing7").exists () &&
                  !SxSymLink (destRootDir + "/dir21/symLinkToNothing7").isFile () &&
                  !SxSymLink (destRootDir + "/dir21/symLinkToNothing7").isDir () &&
                  SxSymLink (destRootDir + "/dir21/symLinkToNothing7").getTarget ()
                  == oldTarget);
   SX_QUIT_TEST ();

   // --- Moving a symbolic link that addresses nothing to a symbolic link
   //     which points to a file
   SX_CREATE_TEST ("");
   oldTarget = SxSymLink (srcRootDir + "/symLinkToNothing8").getTarget ();
   // Preconditions:
   // srcRootDir + "/symLinkToNothing8" is a symbolic link that points to
   // nowhere
   // destRootDir + "/symLinkToFile22" is a valid symbolic link that points
   // to a file
   SX_CHECK_PRE (SxSymLink (srcRootDir + "/symLinkToNothing8").exists () &&
                 !SxSymLink (srcRootDir + "/symLinkToNothing8").isFile () &&
                 !SxSymLink (srcRootDir + "/symLinkToNothing8").isDir ());
   SX_CHECK_PRE (SxSymLink (destRootDir + "/symLinkToFile22").exists () &&
                 SxSymLink (destRootDir + "/symLinkToFile22").isValid () &&
                 SxSymLink (destRootDir + "/symLinkToFile22").isFile ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/symLinkToNothing8",
                                 destRootDir + "/symLinkToFile22"),
                 true);
   // Postconditions:
   // srcRootDir + "/symLinkToNothing8" does not exist
   // destRootDir + "/symLinkToFile22" is a symbolic link that points to
   // the same target as the by now deleted symbolic link
   // srcRootDir + "/symLinkToNothing8"
   SX_CHECK_POST (!SxSymLink (srcRootDir + "/symLinkToNothing8").exists ());
   SX_CHECK_POST (SxSymLink (destRootDir + "/symLinkToFile22").exists () &&
                  !SxSymLink (destRootDir + "/symLinkToFile22").isFile () &&
                  !SxSymLink (destRootDir + "/symLinkToFile22").isDir () &&
                  SxSymLink (destRootDir + "/symLinkToFile22").getTarget ()
                  == oldTarget);
   SX_QUIT_TEST ();

   // --- Moving a symbolic link that addresses nothing to a symbolic link
   //     which points to a directory
   SX_CREATE_TEST ("");
   oldTarget = SxSymLink (srcRootDir + "/symLinkToNothing9").getTarget ();
   // Preconditions:
   // srcRootDir + "/symLinkToNothing9" is a symbolic link that points to nowhere
   // destRootDir + "/symLinkToDir22" is a valid symbolic link that points to a
   // directory
   SX_CHECK_PRE (SxSymLink (srcRootDir + "/symLinkToNothing9").exists () &&
                 !SxSymLink (srcRootDir + "/symLinkToNothing9").isFile () &&
                 !SxSymLink (srcRootDir + "/symLinkToNothing9").isDir ());
   SX_CHECK_PRE (SxSymLink (destRootDir + "/symLinkToDir22").exists () &&
                 SxSymLink (destRootDir + "/symLinkToDir22").isValid () &&
                 SxSymLink (destRootDir + "/symLinkToDir22").isDir ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/symLinkToNothing9",
                                 destRootDir + "/symLinkToDir22"),
                 true);
   // Postconditions:
   // srcRootDir + "/symLinkToNothing9" does not exist
   // destRootDir + "/symLinkToDir22/symLinkToNothing9" is a symbolic link that
   // points to the same target as the by now deleted symbolic link
   // srcRootDir + "/symLinkToNothing9"
   SX_CHECK_POST (!SxSymLink (srcRootDir + "/symLinkToNothing9").exists ());
   SX_CHECK_POST (SxSymLink (destRootDir + "/symLinkToDir22/symLinkToNothing9").exists () &&
                  !SxSymLink (destRootDir + "/symLinkToDir22/symLinkToNothing9").isFile () &&
                  !SxSymLink (destRootDir + "/symLinkToDir22/symLinkToNothing9").isDir () &&
                  SxSymLink (destRootDir + "/symLinkToDir22/symLinkToNothing9").getTarget ()
                  == oldTarget);
   SX_QUIT_TEST ();

   // --- Moving a symbolic link that addresses nothing to another empty
   //     symbolic link
   SX_CREATE_TEST ("");
   oldTarget = SxSymLink (srcRootDir + "/symLinkToNothing10").getTarget ();
   // Preconditions:
   // srcRootDir + "/symLinkToNothing10" is a symbolic link that points to
   // nowhere
   // destRootDir + "/symLinkToNothing11" is a symbolic link that points to
   // nowhere
   SX_CHECK_PRE (SxSymLink (srcRootDir + "/symLinkToNothing10").exists () &&
                 !SxSymLink (srcRootDir + "/symLinkToNothing10").isFile () &&
                 !SxSymLink (srcRootDir + "/symLinkToNothing10").isDir ());
   SX_CHECK_PRE (SxSymLink (destRootDir + "/symLinkToNothing11").exists () &&
                 !SxSymLink (destRootDir + "/symLinkToNothing11").isFile () &&
                 !SxSymLink (destRootDir + "/symLinkToNothing11").isDir ());
   SX_MAKE_TEST (SxFSAction::mv (srcRootDir + "/symLinkToNothing10",
                                 destRootDir + "/symLinkToNothing11"),
                 true);
   // Postconditions:
   // srcRootDir + "/symLinkToNothing10" does not exist
   // destRootDir + "/symLinkToNothing11" is a symbolic link that points to
   // the same target as the by now deleted symbolic link
   // srcRootDir + "/symLinkToNothing10"
   SX_CHECK_POST (!SxSymLink (srcRootDir + "/symLinkToNothing10").exists ());
   SX_CHECK_POST (SxSymLink (destRootDir + "/symLinkToNothing11").exists () &&
                  !SxSymLink (destRootDir + "/symLinkToNothing11").isFile () &&
                  !SxSymLink (destRootDir + "/symLinkToNothing11").isDir () &&
                  SxSymLink (destRootDir + "/symLinkToNothing11").getTarget ()
                  == oldTarget);
   SX_QUIT_TEST ();

}

void mvTests (const SxString &srcRootDir, const SxString &destRootDir)
{
   for (int iRoot = 0; iRoot < 2; ++iRoot)  {
      SxString rootDir;
      if (iRoot == 0)  {
         rootDir = srcRootDir;
      } else  {
         if (srcRootDir == destRootDir)  {
            break;
         } else  {
            rootDir = destRootDir;
         }
      }
      std::cout << "Just checking \"" << rootDir << "\"" << std::endl;//TEST
      SX_CHECK (!SxDir (rootDir).exists ());

      // --- Setting up the directory structure
      try  {
         SxFSAction::mkdir (rootDir + "");
      } catch (SxException ex)  {
         ex.print ();
      }
      const int nDirs = 22;
      const int nFiles = 22;
      for (int iDir = 0; iDir < nDirs; ++iDir)  {
         SxString dirName ("dir");
         dirName += SxString (iDir + 1);
         const SxString symLinkName = "symLinkToDir" + SxString (iDir + 1);
         try  {
            SxFSAction::mkdir (SxString (rootDir + "/") + dirName);
         } catch (SxException ex)  {
            ex.print ();
         }
         try  {
            SxFSAction::ln_sf (dirName, SxString (rootDir + "/") + symLinkName);
         } catch (SxException ex)  {
            ex.print ();
         }
      }
      for (int iFile = 0; iFile < nFiles; ++iFile)  {
         SxString fileName ("file");
         fileName += SxString (iFile + 1);
         SxString imaginaryFileName ("nowhere");
         imaginaryFileName += SxString (iFile + 1);
         const SxString symLinkName = "symLinkToFile" + SxString (iFile + 1);
         const SxString imaginarySymLinkName = "symLinkToNothing"
                          + SxString (iFile + 1);
         try  {
            SxFSAction::touch (SxString (rootDir + "/") + fileName);
         } catch (SxException ex)  {
            ex.print ();
         }
         SxString (SxString ("I was born in '/mvTestDir/") + fileName + "'"
                  ).appendToFile (SxString (rootDir + "/") + fileName);
         try  {
            SxFSAction::ln_sf (fileName, SxString (rootDir + "/") + symLinkName);
         } catch (SxException ex)  {
            ex.print ();
         }
         try  {
            SxFSAction::ln_sf (imaginaryFileName, SxString (rootDir + "/") +
                               imaginarySymLinkName);
         } catch (SxException ex)  {
            ex.print ();
         }

      }
   }
   // --- Moving directories
   mvDirToExisting (srcRootDir, destRootDir);

   // --- Moving files
   mvFileToExisting (srcRootDir, destRootDir);

   // --- Moving symbolic links to directories
   mvDirSymLinkToExisting (srcRootDir,
                           destRootDir);

   // --- Moving symbolic links to files
   mvFileSymLinkToExisting (srcRootDir,
                            destRootDir);

   mvEmptySymLinkToExisting (srcRootDir,
                             destRootDir);

   // Printing the evaluation of the results
   globTestSuite->printResults ();
}

void createTmpFileTest ()
{
   SxString hello ("Hello");
   SxArray<unsigned char> arr (hello.getSize () + 1);
   memcpy (arr.elements, hello.getElems (), (size_t)(hello.getSize () + 1));
   SxFile file;
   try  {
      file = SxFSAction::createTmpFile ("/tmp", arr);
   } catch (SxException ex)  {
      ex.print ();
   }
   cout << "File name of the generated temporary file: "<< file << endl;
}

//------------------------------------------------------------------------------
// main
//------------------------------------------------------------------------------

/** \example fileinfo.cpp

  In this example the usage of the SPHInX file info classes is demonstrated.

  \brief  Basic usage of SxFileInfo, SxFile, SxDir and SxSymLink
  \author Thomas Uchdorf
  */
int main (int argc, char **argv)
{
   // --- Using the command line interface to retrieve the program arguments
   SxCLI cli (argc, argv);
   bool enableMvRaid = cli.option ("--mvRaid",
                                   "Enables \"mv\"-tests for /raid"
                                  ).toBool ();
   bool enableTmp = cli.option ("--tmp",
                                "Enables \"tmp\""
                               ).toBool ();
   cli.finalize ();

   if (enableMvRaid) mvTests ((SxFSAction::pwd ()/"mvTestDir").getAbsPath (),
                              "/raid/uchdorf/mvTestDir");
   if (enableTmp) createTmpFileTest ();

   return 0;
}
