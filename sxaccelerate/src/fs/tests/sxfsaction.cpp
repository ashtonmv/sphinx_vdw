// --- Including header-files
#include <iostream>
#include <SxCLI.h>
#include <SxString.h>
#include <SxFile.h>
#include <SxFSAction.h>
#include <SxFileIO.h>
#include <SxTime.h>
#include <SxSymLink.h>

void SxFS_createRootFolder (const SxString &rootPath);

void SxFS_touchFile (const SxString &rootPath,
                     const SxString &fileName);

void SxFS_tmpFile (const SxString &rootPath);

void SxFS_mkrmdir (const SxString &rootPath,
                   const SxString &dirName);

void SxFS_cpmvFile (const SxString &rootPath);

void SxFS_symlink (const SxString &rootPath,
                   const SxString &fileName);

void SxFS_cdDir (const SxString &rootPath);

void SxFS_lsDir (const SxString &rootPath);

int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   cli.finalize ();

   try  {

      // --- Langauges that will be a good test scenario.

      // English  - 1 byte  - [example]     [word]
      // Greek    - 2 bytes - [παράδειγμα]  [λέξη]
      // Russian  - 2 bytes - [пример]      [слово]
      // Chinease - 3 bytes - [例]          [字]
      // Thai     - 3 bytes - [ตัวอย่าง]      [คำ]

      // --------- Basic tests

      // GET TMP
      if (SxFileInfo::getTmpStr () == "")
         SX_THROW ("[TMPDIR] or [TMP] environment variable is was not found."
                   "Please verify that the environment variable is set.");
      // GET TMP
      SxDir tmpDir = SxDir(SxFileInfo::getTmpStr ());

      if (!tmpDir.exists ())
         SX_THROW ("[TMPDIR] or [TMP] directory not found. Please verify "
                   "that the directory exists");

      // GET HOME
      if (SxFileInfo::getHomeStr () == "")
         SX_THROW ("'HOME' environment variable is was not found."
                   "Please verify that the environment variable is set.");
      SxDir homeDir = SxDir (SxFileInfo::getHomeStr ());
      if (!homeDir.exists ())
         SX_THROW ("'HOME' directory not found. Please verify that the"
                   "directory exists");

      // Unicode folder paths
      SxList<SxString> unicodePaths;
      unicodePaths << "example"
                   << SxString::unicodeFromUtf8 (u8"παράδειγμα")
                   << SxString::unicodeFromUtf8 (u8"пример")
                   << SxString::unicodeFromUtf8 (u8"例")
                   << SxString::unicodeFromUtf8 (u8"ตัวอย่าง");

      if (unicodePaths.getSize () == 0)  {
        SX_THROW ("No language test cases specified.");
     }

     // Folder in which all tests will be performed
     tmpDir = tmpDir.getAbsPath () + SxString ("/SxFSActionTest/");

      // Clean old tests if they exist
      if (SxFSAction::test_d (tmpDir.getAbsPath ())) {
         SxFSAction::rm_r (tmpDir);
      }

      // Check combinations of Unicode, with different encoding size
      for (const SxString &path1 : unicodePaths)  {
         for (const SxString &path2 : unicodePaths)  {

            // Skip repeating cases, except full ASCII test
            if (path1 == path2 && path1 != unicodePaths(0)) continue;

            // All the tests go here.
            SxString rootPath = tmpDir.getAbsPath () + "/" + path1 + "/"
                              + path2 + "/";
            // Just in case the rootPath is already unicode
            rootPath = SxString::unicodeFromUtf8 (rootPath.getElems ());

            SxFS_createRootFolder (rootPath);

            SxFS_touchFile (rootPath, path1);

            SxFS_tmpFile (rootPath);

            SxFS_mkrmdir (rootPath, path1);

            SxFS_cpmvFile (rootPath);

            SxFS_symlink (rootPath, path1);

            SxFS_cdDir (rootPath);

            // This operation may take some time because of FS flush
            SxFS_lsDir (rootPath);

            // --- Remove the root test folder
            SxString cleanPath = tmpDir.getAbsPath () + "/" + path1 + "/";
            SxFSAction::rm_r (cleanPath);

            if (SxFSAction::test_d (rootPath))
               SX_THROW ("Failed to remove test folder '" + rootPath + "'");

            if (SxFSAction::test_d (cleanPath))
               SX_THROW ("Failed to remove root folder '" + cleanPath + "'");
         }
      }

   } catch (SxException e)  {
      cout << e.toString ().wrap ("ERROR: ") << endl;
      cerr << "Exception stack:" << endl;
      cerr << e.toString (SxException::DebugStack) << endl;
      SX_EXIT;
   }

   return 0;
}


void SxFS_createRootFolder (const SxString &rootPath)
{
   SX_TRACE ();
   try {
      // --- Create root test folder
      //  test_d, mkdir_p

      // Folder should not be present
      if (SxFSAction::test_d (rootPath))
         SX_THROW ("Testing folder '" + rootPath + "' already exists");

      // Create test folder
      SxFSAction::mkdir_p (rootPath);

      // Folder should be represent
      if (!SxFSAction::test_d (rootPath))
         SX_THROW ("SxFSAction::mkdir_p failed to create '"
                  + rootPath + "' folder");

   } catch (SxException e) {
      SX_THROW (e, "SxFS_createRootFolder Failed.");
   }
}

void SxFS_touchFile (const SxString &rootPath, const SxString &fileName)
{
   SX_TRACE ();
   try {
      // --- Create and remove file
      //  test_f, test_f, touch, rm

      SxString touchFile = rootPath + fileName + ".txt";

      // File should not be present
      if (SxFSAction::test_f (touchFile))
         SX_THROW ("File '" + touchFile + "' already exists");

      if (SxFSAction::exists (touchFile))
         SX_THROW ("Object '" + touchFile + "' already exits");

      // Create empty file
      SxFSAction::touch (touchFile);

      // File should be represent
      if (!SxFSAction::test_f (touchFile))
         SX_THROW ("File '" + touchFile + "' not found");

      if (!SxFSAction::exists (touchFile))
         SX_THROW ("Object '" + touchFile + "' not found");

#     ifndef WIN32

         SxFileInfo touchInfo(touchFile);
         int origMode = touchInfo.getMode ();

         // Everyone should be able to read it.
         if (origMode != 0644)
            SX_THROW ("Default access rights of touch file are incorrect");

         // Make it readyonly, and visible only for owner
         SxFSAction::chmod (0400, touchInfo);

         if (touchInfo.getMode () != 0400)
            SX_THROW ("SxFSAction::chmod failed to change access rights.");
#     endif
      SxFSAction::rm (touchFile);

      // Check that file was deleted
      if (SxFSAction::test_f (touchFile))
         SX_THROW ("SxFSAction::rm failed to remove  file'" + touchFile + "'");

   } catch (SxException e) {
      SX_THROW (e, "SxFS_touchFile Failed");
   }
}

void SxFS_tmpFile (const SxString &rootPath)
{
   SX_TRACE ();
   try {

      // --- Create temp file
      // createTmpFile, rm
      SxString tempFile = SxFSAction::createTmpFile ().getAbsPath ();

      if (!SxFSAction::test_f (tempFile))
         SX_THROW ("SxFSAction::createTmpFile failed to create a temporary"
                   " file '" + tempFile + "'");

      // Remove temp file
      SxFSAction::rm (tempFile);

      // Create new temp file in unicode path
      tempFile = SxFSAction::createTmpFile (rootPath).getAbsPath ();

      if (!SxFSAction::test_f (tempFile))
         SX_THROW ("SxFSAction::rm failed to remove file '" + tempFile + "'");

      // Remove temp file
      SxFSAction::rm (tempFile);

   } catch (SxException e) {
      SX_THROW (e, "SxFS_tmpFile failed");
   }
}

void SxFS_mkrmdir (const SxString &rootPath, const SxString &dirName)
{
   SX_TRACE ();
   try {

      // --- Create and remove folder
      //  test_d, mkdir, rmdir
      SxString tempDir = rootPath + dirName;

      // File should not be present
      if (SxFSAction::test_d (tempDir))
         SX_THROW ("Folder '" + tempDir + "' already exists");

      if (SxFSAction::exists (tempDir))
         SX_THROW ("Object '" + tempDir + "' already exits");

      // Create empty file
      SxFSAction::mkdir (tempDir);

      // File should be represent
      if (!SxFSAction::test_d (tempDir))
         SX_THROW ("Folder '" + tempDir + "' not found");

      if (!SxFSAction::exists (tempDir))
         SX_THROW ("Object '" + tempDir + "' already exits");

      SxFSAction::rmdir (tempDir);

      // Check that file was deleted
      if (SxFSAction::test_d (tempDir))
         SX_THROW ("SxFSAction::rmdir failed to remove folder '"
                  + tempDir + "'");

   } catch (SxException e) {
      SX_THROW (e, "SxFS_mkrmdir failed");
   }
}

void SxFS_cpmvFile (const SxString &rootPath)
{
   SX_TRACE ();
   try {

      // ---  Copy and move file
      SxString origFile = SxFSAction::createTmpFile (rootPath).getAbsPath ();
      SxString origFolder = rootPath + "folder";
      SxFSAction::mkdir (origFolder);

      SxList<SxString> objects;
      objects << origFile << origFolder;

#     ifndef WIN32
         SxString origFileLnk = origFile + "_symlnk";
         SxFSAction::ln_sf (origFile, origFileLnk);

         SxString origFolderLnk = origFolder + "_symlnk";
         SxFSAction::ln_sf (origFolder, origFolderLnk);

         objects << origFileLnk << origFolderLnk;
#     endif


      if (!SxFSAction::test_f (origFile)) {
         SX_THROW ("SxFSAction::createTmpFile failed to create test file ''"
                  + origFile + "'");
      }

      for (const SxString &obj : objects) {

         SxFileInfo objFile = SxFileInfo (obj);
         SxString newObj = obj + "_mv";

         if (objFile.isFile ())  {

            if (SxFSAction::test_f (newObj)) {
               SX_THROW ("Test folder not is not empty");
            }

         } else if (objFile.isDir ())  {

            if (SxFSAction::test_d (newObj)) {
               SX_THROW ("Test folder not is not empty");
            }
         }

         int origMode = objFile.getMode ();

         // Move file
         SxFSAction::mv (obj, newObj);

         if (objFile.isFile ())  {

            if (   SxFSAction::test_f (obj)
               || !SxFSAction::test_f (newObj))
            {
               SX_THROW ("SxFSAction::mv failed to move file from '" + obj
                        + "' to '" + newObj + "'");
            }

         } else if (objFile.isDir ())  {

            if (   SxFSAction::test_d (obj)
               || !SxFSAction::test_d (newObj))
            {
               SX_THROW ("SxFSAction::mv failed to move file from '" + obj
                        + "' to '" + newObj + "'");
            }

         }

         if (origMode != SxFileInfo (newObj).getMode ())
            SX_THROW ("SxFSAction::mv file access rights are not the same.");

         // Copy file
         SxFSAction::cp (newObj, obj);

         if (origMode != SxFileInfo (obj).getMode ())
            SX_THROW ("SxFSAction::cp file access rights are not the same.");

         if (objFile.isFile ())  {
            if (!SxFSAction::test_f (obj))  {
               SX_THROW ("SxFSAction::cp failed to copy file from '" + newObj
                        + "' to '" + obj + "'");
            }

         } else if (objFile.isDir ())  {
            if (!SxFSAction::test_d (obj))  {
               SX_THROW ("SxFSAction::cp failed to copy file from '" + newObj
                        + "' to '" + obj + "'");
            }
         }

         SxFSAction::rm_r (obj);
         SxFSAction::rm_r (newObj);
         if (SxFileInfo(obj).exists ())
            SX_THROW (obj + " was not deleted");
         if (SxFileInfo(newObj).exists ())
            SX_THROW (newObj + " was not deleted");
      }

   } catch (SxException e) {
      SX_THROW (e, "SxFS_cpmvFile failed");
   }
}


void SxFS_symlink (const SxString &rootPath,
                   const SxString &fileName)
{
   SX_TRACE ();
   try {
#     ifndef WIN32
         // ---  Copy and move file
         SxString origFile = rootPath + fileName;
         SxString lnkFile = origFile + ".lnk";

         if (  SxFSAction::test_f (origFile)
            || SxFSAction::test_f (lnkFile))
         {
            SX_THROW ("Testing folder not empty.");
         }

         SxFSAction::touch (origFile);

         if (!SxFSAction::test_f (origFile)) {
            SX_THROW ("SxFSAction::touch failed to create file '" + origFile + "'");
         }

         SxFSAction::ln_sf (origFile, lnkFile);

         if (!SxFSAction::test_f (lnkFile))
            SX_THROW ("SxFSAction::ln_sf failed to create sym link file '"
                     + lnkFile + "' to '" + origFile + "'");

         if (  !SxFSAction::test_L (lnkFile)
            || SxFSAction::test_L (origFile))
         {
            SX_THROW ("SxFSAction::test_L works incorrectly.");
         }

         SxList<SxSymLink> symlinks = SxFSAction::getSymLinks (rootPath);
         if (symlinks.getSize () == 0)
            SX_THROW ("SxFSAction::getSymLinks returned empty array.");

         SxString resOrig = symlinks.first ().getTarget ();
         resOrig = SxString::unicodeFromUtf8 (resOrig.getElems ());
         SxString resLnk = symlinks.first ().getAbsPath ();
         resLnk = SxString::unicodeFromUtf8 (resLnk.getElems ());

         if (  resLnk  != lnkFile
            || resOrig != origFile)
         {
            SX_THROW ("SxFSAction::getSymLinks returned unexpected result");
         }

         SxString cpLnkStr = lnkFile + "_copy";
         SxFSAction::cp (lnkFile, cpLnkStr);
         SxSymLink cpLnk (cpLnkStr);
         if (!cpLnk.isValid ())
            SX_THROW ("SxFSAction::cp did not create a valid lnk file");

         SxString cpLnkTarget = SxString::unicodeFromUtf8 (cpLnk.getTarget ().getElems ());
         if (cpLnkTarget != origFile)
            SX_THROW ("SxFSAction::cp of symlink file has invalid target file");

         SxString mvLnkStr = lnkFile + "_move";
         SxFSAction::mv (cpLnkStr, mvLnkStr);
         SxSymLink mvLnk (mvLnkStr);
         if (!mvLnk.isValid ())
            SX_THROW ("SxFSAction::mv did not create a valid lnk file");

         SxString mvLnkTarget = SxString::unicodeFromUtf8 (mvLnk.getTarget ().getElems ());
         if (mvLnkTarget != origFile)
            SX_THROW ("SxFSAction::mv of symlink file has invalid target file");

         SxFSAction::rm (origFile);
         SxFSAction::rm (lnkFile);
         SxFSAction::rm (mvLnkStr);
#     endif

   } catch (SxException e) {
      SX_THROW (e, "SxFS_symlink failed");
   }
}
void SxFS_cdDir (const SxString &rootPath)
{
   SX_TRACE ();

   try {

      // Original directory
      SxString origPath = SxFSAction::pwd ().getAbsPath ();
      SxFileInfo origPathFile (origPath);
      // Check that the current directory is not the same as
      // the testing one
      if (  !origPathFile.exists ()
         || !origPathFile.isDir ())
         SX_THROW ("SxFSAction::pwd returned invalid path '" + origPath + "'");

      // Change to testing directory
      SxFSAction::cd (rootPath);

      // Verify that the current directory has been changed
      SxString currentPath = SxFSAction::pwd ().getAbsPath () + "/";
      if (currentPath != rootPath)
         SX_THROW ("SxFSAction::cd failed to change current directory.");

      // Push current directory to stack and change directory
      SxFSAction::pushd (origPath);

      // Check that current directory has been changed
      if (SxFSAction::pwd ().getAbsPath () != origPath)
         SX_THROW ("SxFSAction::pushd did not change the current directory");

      // Return the directory from the stack
      SxFSAction::popd ();

      // Verify that the returned directory is the same
      currentPath = SxFSAction::pwd ().getAbsPath () + "/";
      currentPath = SxString::unicodeFromUtf8 (currentPath.getElems ());
      if (currentPath != rootPath)
         SX_THROW ("SxFSAction::popd did not returned '" + currentPath
                  + "', while the expected result is '" + rootPath + "'");

      // Return to original path
      SxFSAction::cd (origPath);

   } catch (SxException e) {
      SX_THROW (e, "SxFS_cdDir failed");
   }
}

void SxFS_lsDir (const SxString &rootPath)
{
   SX_TRACE ();
   try {

      // --- Create a file structure which will be tested.
      SxString dataSmall  = "123";
      SxString dataMedium = "12345";
      SxString dataBig    = "1234567890";
      SxFileInfo fileSmall  = rootPath + "small.txt";
      SxFileInfo fileMedium = rootPath + "medium.txt";
      SxFileInfo fileBig    = rootPath + "big.txt";

      // --- Timestamp precision on different file systems may be different
      // ext3  - 1 second
      // ext4  - 1 nanosecond
      // ntfs  - 100 nanoseconds
      // fat   - 2 seconds

      // --- The major problem is that if several writing request are called
      //     close to each other, they can be packaged as single write operation
      //     which causes filesystem to assign same timestamps to them.
      // Milisecond precision is not accessible by default on Windows
#     if defined(LINUX)
         int fsFlushDelayMS = 15;
#     else
         int fsFlushDelayMS = 1000;
#     endif

      //--- There are 7 sleep statments.

      dataSmall.write (fileSmall.getAbsPath ());
      SxTime::msleep (fsFlushDelayMS);

      dataMedium.write (fileMedium.getAbsPath ());
      SxTime::msleep (fsFlushDelayMS);

      dataBig.write (fileBig.getAbsPath ());
      SxTime::msleep (fsFlushDelayMS);

      // An empty folder
      SxFileInfo emptyFolder = rootPath + "emptyFolder";
      SxFSAction::mkdir (emptyFolder);
      SxTime::msleep (fsFlushDelayMS);

      // Folder with another folder inside
      SxFileInfo folder1 = rootPath + "folder1";
      SxFSAction::mkdir (folder1);
      SxTime::msleep (fsFlushDelayMS);

      SxFileInfo recursiveFolder = folder1.getAbsPath () + "/folder";
      SxFSAction::mkdir_p (recursiveFolder);
      SxTime::msleep (fsFlushDelayMS);

      // Folder with a file inside
      SxFileInfo folder2 = rootPath + "folder2";
      SxFSAction::mkdir (folder2);
      SxTime::msleep (fsFlushDelayMS);

      SxFileInfo recursiveFile = folder2.getAbsPath () + "/file.txt";
      dataMedium.write (recursiveFile.getAbsPath ());

      // LS
      SxList<SxFileInfo> resultLs = SxFSAction::ls (rootPath + "/*");
      // 3 files + 2 folders
      if (resultLs.getSize () != 6)
         SX_THROW ("SxFSAction::ls returned incorrect number of elements.");
      if (  !resultLs.contains (fileSmall)
         || !resultLs.contains (fileMedium)
         || !resultLs.contains (fileBig)
         || !resultLs.contains (emptyFolder)
         || !resultLs.contains (folder1)
         || !resultLs.contains (folder2))
      {
         SX_THROW ("SxFSAction::ls returned unexpected content.");
      }

      // Get files
      SxList<SxFile> resultGF = SxFSAction::getFiles (SxFileInfo(rootPath));
      // 3 files
      if (resultGF.getSize () != 3)
         SX_THROW ("SxFSAction::getFiles returned "
                   "incorrect number of elements.");

      if (  !resultGF.contains (fileSmall)
         || !resultGF.contains (fileMedium)
         || !resultGF.contains (fileBig))
      {
         SX_THROW ("SxFSAction::getFiles returned unexpected content.");
      }

      // Get dirs
      SxList<SxDir> resultGD = SxFSAction::getDirs (SxFileInfo(rootPath));
      // 2 folders
      if (resultGD.getSize () != 3)
         SX_THROW ("SxFSAction::getDirs returned incorrect number of elements");

      if (  !resultLs.contains (emptyFolder)
         || !resultGD.contains (folder1)
         || !resultGD.contains (folder2))
      {
         SX_THROW ("SxFSAction::getDirs returned unexpected content.");
      }

      // LS recursively
      SxList<SxString> resultListDir = SxFSAction::listDir (rootPath);
      if (resultListDir.getSize () != 8)
         SX_THROW ("SxFSAction::listDir returned incorrect number of elements");

      // Should contain the LS content + few additional element
      if (  !resultListDir.contains (fileSmall.getAbsPath ())
         || !resultListDir.contains (fileMedium.getAbsPath ())
         || !resultListDir.contains (fileBig.getAbsPath ())
         || !resultListDir.contains (emptyFolder.getAbsPath ())
         || !resultListDir.contains (folder1.getAbsPath ())
         || !resultListDir.contains (folder2.getAbsPath ())
         || !resultListDir.contains (recursiveFolder.getAbsPath ())  // Recursive folder
         || !resultListDir.contains (recursiveFile.getAbsPath ())    // Recursive file
         )
      {
         SX_THROW ("SxFSAction::listDir returned unexpected content");
      }

      // LS_S
      SxList<SxFISortedBySize> resultlss = SxFSAction::ls_S (rootPath + "/*");
      if (resultlss.getSize () != 6)
         SX_THROW ("SxFSAction::ls_S returned incorrect number of elements");

      // On Windows folder size is always 0
      // On Linux ext4 it is 4096
      // On MacOS it is 64 for empty and 96 for non empty folder
      ssize_t idxEmpty = resultlss.findPos (emptyFolder);
      ssize_t idxdir1  = resultlss.findPos (folder1);
      ssize_t idxdir2  = resultlss.findPos (folder2);
      ssize_t idxSmall = resultlss.findPos (fileSmall);
      ssize_t idxMed   = resultlss.findPos (fileMedium);
      ssize_t idxBig   = resultlss.findPos (fileBig);

      if (
#     if defined(WIN32)
            (idxEmpty != 0 && idxEmpty != 1 && idxEmpty != 2)
         || (idxdir1  != 0 && idxdir1  != 1 && idxdir1  != 2)
         || (idxdir2  != 0 && idxdir2  != 1 && idxdir2  != 2)
         || idxSmall  != 3
         || idxMed    != 4
         || idxBig    != 5
#     elif defined(LINUX)
            idxSmall  != 0
         || idxMed    != 1
         || idxBig    != 2
         || (idxEmpty != 3 && idxEmpty != 4 && idxEmpty != 5)
         || (idxdir1  != 3 && idxdir1  != 4 && idxdir1  != 5)
         || (idxdir2  != 3 && idxdir2  != 4 && idxdir2  != 5)
#     elif defined(MACOSX)
            idxSmall != 0
         || idxMed   != 1
         || idxBig   != 2
         || idxEmpty != 3
         || (idxdir1 != 4 && idxdir1 != 5)
         || (idxdir2 != 4 && idxdir2 != 5)
#     endif

         )
      {
         SX_THROW ("SxFSAction::ls_S returned unexpected content.");
      }

      // LS_T
      SxList<SxFISortedByTime> resultlst = SxFSAction::ls_t (rootPath + "/*");
      if (resultlst.getSize () != 6)
         SX_THROW ("SxFSAction::ls_t returned incorrect number of elements.");
      if (  (SxFileInfo)resultlst(0) != fileSmall
         || (SxFileInfo)resultlst(1) != fileMedium
         || (SxFileInfo)resultlst(2) != fileBig
         || (SxFileInfo)resultlst(3) != emptyFolder // Contains nothing
         || (SxFileInfo)resultlst(4) != folder1     // Contains an empty folder
         || (SxFileInfo)resultlst(5) != folder2     // Contains a small file
         )
      {
         SX_THROW ("SxFSAction::ls_t returned unexpected content.");
      }

      SxList<SxFileInfo> resultFind = SxFSAction::find (rootPath + "*.txt");
      if (resultFind.getSize () != 4)
         SX_THROW ("SxFSAction::find returned incorrect number of elements");
      if (  !resultFind.contains (recursiveFile)
         || !resultFind.contains (fileSmall)
         || !resultFind.contains (fileMedium)
         || !resultFind.contains (fileBig)
         )
      {
         SX_THROW ("SxFSAction::find returned unexpected content.");
      }

      // SxFSAction::which returns the first item of SxFSAction::where
      SxString sxunicodeFile;
#     ifdef WIN32
         sxunicodeFile = SxFSAction::which ("dxtest.exe");
#     else
         sxunicodeFile = SxFSAction::which ("sxfsaction");
#     endif
      if (sxunicodeFile == "")
         SX_THROW ("SxFSAction::which returned empty string");

   } catch (SxException e) {
      SX_THROW (e, "SxFS_lsDir failed");
   }
}
