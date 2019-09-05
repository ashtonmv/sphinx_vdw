// --- Including header-files
#include <iostream>
#include <SxCLI.h>
#include <SxString.h>
#include <SxFile.h>
#include <SxFSAction.h>
#include <SxFileIO.h>
#include <SxTime.h>

#ifdef WIN32
#  include <io.h>
#  include <fcntl.h>
#  include <iosfwd>
#  include <string>
#  include <sstream>
#  include <iostream>
#endif

void SxFS_FileInfo (const SxString &rootPath);
void SxFS_File     (const SxString &rootPath);
void SxFS_Dir      (const SxString &rootPath);
void SxFS_FileInfoUnit (const SxFileInfo &obj1,
                        const SxFileInfo &obj2);

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

      if (SxFileInfo::getTmpStr () == "")
         SX_THROW ("[TMPDIR] or [TMP] environment variable is was not found. "
                   "Please verify that the environment variable is set.");
      // GET TMP
      SxDir tmpDir = SxDir (SxFileInfo::getTmpStr ());
      if (!tmpDir.exists ())
         SX_THROW ("[TMPDIR] or [TMP] directory not found. Please verify "
                   "that the directory exists");

      // GET HOME
      if (SxFileInfo::getHomeStr () == "")
         SX_THROW ("'HOME' environment variable is was not found. "
                   "Please verify that the environment variable is set.");
      SxDir homeDir = SxDir (SxFileInfo::getHomeStr ());
      if (!homeDir.exists ())
         SX_THROW ("'HOME' directory not found. Please verify that the"
                   "directory exists");

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
     tmpDir = tmpDir.getAbsPath () + SxString("/SxFileTest/");

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
            SxString rootPath = tmpDir.getAbsPath ()
                              + "/" + path1 + "/" + path2 + "/";
            // Just in case the rootPath is already unicode
            rootPath = SxString::unicodeFromUtf8 (rootPath.getElems ());

            // --- Assume that this function call is working.
            //     Chicken and egg problem. SxFSAction should be tested
            //     only after SxFile, but we need folders....
            SxFSAction::mkdir_p (rootPath);

            SxFS_FileInfo (rootPath);

            SxFS_File (rootPath);

            SxFS_Dir (rootPath);

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
      cout << "Exception stack:" << endl;
      cout << e.toString (SxException::DebugStack) << endl;
      SX_EXIT;
   }

   return 0;
}

void SxFS_File (const SxString &rootPath)
{
   SX_TRACE ();
   try {

      SxString filePath = rootPath + "SxFile.txt";
      SxFile file(filePath);
      SxFile pathFile(rootPath);

      if (file.exists ())
         SX_THROW ("SxFile::exists unexpected result");

      if (pathFile.exists ())
         SX_THROW ("SxFile::exists unexpected result");

      SxString ("12345").write (filePath);
      file.setDirty ();

      if (!file.exists ())
         SX_THROW ("SxFile::exists unexpected result");

      if (file.getIncludeName (rootPath) != file.getName ())
         SX_THROW ("SxFile::getIncludeName unexpected result");

   } catch (SxException e) {
      SX_THROW (e, "SxFS_File Failed.");
   }
}

void SxFS_Dir (const SxString &rootPath)
{
   SX_TRACE ();
   try {

      SxString dirPath = rootPath + "newfolder";
      SxDir dir (dirPath);
      SxDir rootDir (rootPath);

      if (dir.exists ())
         SX_THROW ("SxDir::exists unexpected result");

      if (!rootDir.exists ())
         SX_THROW ("SxDir::exists unexpected result");

      SxFSAction::mkdir (dirPath);
      dir.setDirty ();

      if (!dir.exists ())
         SX_THROW ("SxDir::exists unexpected result");

#     ifdef WIN32
         SxFileInfo exeFile (SxDir::getExecPath () + "/dxtest.exe");
#     else
         SxFileInfo exeFile (SxDir::getExecPath () + "/sxfile");
#     endif

      if (  !exeFile.exists ()
         || !exeFile.isExecutable ())
         SX_THROW ("SxDir::getExecPath unexpected result");

      // Searches for [share] folder in exePath.
      // But in 'examples' and higher, such folder does not exist.
      //SxDir::getInstallPath () << endl;

   } catch (SxException e) {
      SX_THROW (e, "SxFS_Dir Failed.");
   }
}

void SxFS_FileInfoUnit (const SxFileInfo &obj1,
                        const SxFileInfo &obj2)
{
   SX_TRACE ();
   try {

      // --- Even though file exists, SxFileInfo object was created before
      //     object existed. We need to force 'update' to change values.
      if (obj1.exists ())
         SX_THROW ("SxFileInfo::exists unexpected result.");

      // Set the object to 'Dirty' and force the update.
      obj1.setDirty ();

      // Object should already exist exists yet
      if (!obj1.exists ())
         SX_THROW ("SxFileInfo::exists unexpected result.");

      // Different files
      if (SxFileInfo::equals (obj1, obj2))
         SX_THROW ("SxFileInfo::equals unexpected result");

      // Same file, different objects
      if (!SxFileInfo::equals (obj1, SxFileInfo (obj1.getAbsPath ())))
         SX_THROW ("SxFileInfo::equals unexpected result");

      if (!obj1.lastAccessed ())
         SX_THROW ("SxFileInfo::lastAccessed unexpected return value");
      if (!obj1.lastModified ())
         SX_THROW ("SxFileInfo::lastModified unexpected return value");


      if (obj1.lastAccessedMS () / 1000 != obj1.lastAccessed ()) {
         SX_THROW ("SxFileInfo::lastAccessed/lastAccessedMS have "
                     "inconsistent return values");
      }

      if (obj1.lastModifiedMS () / 1000 != obj1.lastModified ())
         SX_THROW ("SxFileInfo::lastModified/lastModifiedMS have "
                  "inconsistent return values");

      if (obj1.lastModifiedMS () == obj2.lastModifiedMS ())
         SX_THROW ("SxFileInfo::lastModifiedMS unexpected result.");

      if (obj1.lastAccessedMS () == obj2.lastAccessedMS ())
         SX_THROW ("SxFileInfo::lastAccessedMS unexpected result.");

      // If symlink, check that the timing is same as the original object
      if (obj1.isSymLink ()) {

         SxFileInfo target(SxSymLink (obj1.getAbsPath ()).getTarget ());
         if (target.lastModifiedMS () != obj1.lastModifiedMS ())
            SX_THROW ("SxFileInfo::lastModifiedMS unexpected result.");
         if (target.lastAccessedMS () != obj1.lastAccessedMS ())
            SX_THROW ("SxFileInfo::lastAccessedMS unexpected result.");
      }

      if (SxFileInfo::compareModificationTime (obj1, obj2) != -1) {
         SX_THROW ("SxFileInfo::compareModificationTime unexpected result");
      }
      if (SxFileInfo::compareModificationTime (obj2, obj1) != 1) {
         SX_THROW ("SxFileInfo::compareModificationTime unexpected result");
      }
      if (SxFileInfo::compareModificationTime (obj1,
                                               SxFileInfo (obj1.getAbsPath ())) != 0) {
         SX_THROW ("SxFileInfo::compareModificationTime unexpected result");
      }

      if (obj1.isFile ()) {

         if (SxFileInfo::compareSize (obj1,
                                      obj2) != -1) {
            SX_THROW ("SxFileInfo::compareSize unexpected result");
         }
         if (SxFileInfo::compareSize (obj2,
                                      obj1) != 1) {
            SX_THROW ("SxFileInfo::compareSize unexpected result");
         }

         SxFileInfo cpObj = obj1.getAbsPath () + "_copy";
         SxFSAction::cp (obj1.getAbsPath (), cpObj.getAbsPath ());
         cpObj.setDirty ();
         // Same content
         if (SxFileInfo::compareSize (obj1,
                                      cpObj) != 0) {
            SX_THROW ("SxFileInfo::compareSize unexpected result");
         }

      }

#     ifndef WIN32

         int origMod = obj1.getMode ();

         // Can not set executable to folder
         if (obj1.isFile ())  {

            // Exe and read only
            SxFSAction::chmod (0555, obj1);
            if (  !obj1.isExecutable (SxFileInfo::AccessGroup::Owner)
               || !obj1.isExecutable (SxFileInfo::AccessGroup::Group)
               || !obj1.isExecutable (SxFileInfo::AccessGroup::Other)
               || !obj1.isExecutable (SxFileInfo::AccessGroup::All))
            {
               SX_THROW ("SxFileInfo::isExecutable unexpected result");
            }
         }

         // Make writeonly
         SxFSAction::chmod (0222, obj1);
         if (  obj1.isExecutable (SxFileInfo::AccessGroup::Owner)
            || obj1.isExecutable (SxFileInfo::AccessGroup::Group)
            || obj1.isExecutable (SxFileInfo::AccessGroup::Other)
            || obj1.isExecutable (SxFileInfo::AccessGroup::All))
         {
            SX_THROW ("SxFileInfo::isExecutable unexpected result");
         }

         // --- Check writable
         if (  !obj1.isWritable (SxFileInfo::AccessGroup::Owner)
            || !obj1.isWritable (SxFileInfo::AccessGroup::Group)
            || !obj1.isWritable (SxFileInfo::AccessGroup::Other)
            || !obj1.isWritable (SxFileInfo::AccessGroup::All))
         {
            SX_THROW ("SxFileInfo::isWritable unexpected result");
         }

         // Make readonly
         SxFSAction::chmod (0444, obj1);

         if (  obj1.isWritable (SxFileInfo::AccessGroup::Owner)
            || obj1.isWritable (SxFileInfo::AccessGroup::Group)
            || obj1.isWritable (SxFileInfo::AccessGroup::Other)
            || obj1.isWritable (SxFileInfo::AccessGroup::All))
         {
            SX_THROW ("SxFileInfo::isWritable unexpected result");
         }

         // --- Check readable
         if (  !obj1.isReadable (SxFileInfo::AccessGroup::Owner)
            || !obj1.isReadable (SxFileInfo::AccessGroup::Group)
            || !obj1.isReadable (SxFileInfo::AccessGroup::Other)
            || !obj1.isReadable (SxFileInfo::AccessGroup::All))
         {
            SX_THROW ("SxFileInfo::isReadable unexpected result");
         }

          SxFSAction::chmod (0000, obj1);
          if (  obj1.isReadable (SxFileInfo::AccessGroup::Owner)
             || obj1.isReadable (SxFileInfo::AccessGroup::Group)
             || obj1.isReadable (SxFileInfo::AccessGroup::Other)
             || obj1.isReadable (SxFileInfo::AccessGroup::All))
          {
             SX_THROW ("SxFileInfo::isReadable unexpected result");
          }

          SxFSAction::chmod (0644, obj1);

         // Default mode
         if (obj1.getMode () != 0644)
            SX_THROW ("SxFileInfo::getMode unexpected result");
         if (obj1.getPerms () != "rw-r--r--")
               SX_THROW ("SxFileInfo::getPerms unexpected result");

         if (obj1.getMode () != SxFileInfo::getMode (obj1.getPerms ()))
            SX_THROW ("static SxFileInfo::getPerms/getMode unexpected result");

         if (obj1.getPerms () != SxFileInfo::getPerms (obj1.getMode ()))
            SX_THROW ("static SxFileInfo::getPerms/getMode unexpected result");

         SxFSAction::chmod (origMod, obj1.getAbsPath ());
         SxFSAction::chmod (origMod, obj2.getAbsPath ());
#     endif

   } catch (SxException e) {
      SX_THROW (e, "SxFS_FileInfoUnit Failed.");
   }
}
void SxFS_FileInfo (const SxString &rootPath)
{
   SX_TRACE ();
   try {

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
         unsigned long fsFlushDelayMS = 15;
#     else
         unsigned long fsFlushDelayMS = 1000;
#     endif

      // --- File

      SxString file1PathName = "file1.txt";
      SxString file1PathStr = rootPath + file1PathName;
      SxFileInfo file1Path = SxFileInfo (file1PathStr);

      // File should not exists yet
      if (file1Path.exists ())
         SX_THROW ("SxFileInfo::exists unexpected result.");

      // Write file
      SxString file1Data ("12345");
      file1Data.write (file1PathStr);
      SxTime::msleep (fsFlushDelayMS);
      SxString file2PathStr = rootPath + "file2.txt";

      SxString file2Data ("123456");
      file2Data.write (file2PathStr);
      SxFileInfo file2Path = SxFileInfo (file2PathStr);

      // Some general path tests
      if (file2Path.getSize () != file2Data.getSize ())
         SX_THROW ("SxFileInfo::getSize unexpected result");

      if (!SxFileInfo::equals (file1Path.getAbsPath (), file1PathStr))
         SX_THROW ("SxFileInfo::getAbsPath unexpected result");

      if (!SxFileInfo::equals (rootPath, SxFileInfo::getPath (file1PathStr) + "/"))
         SX_THROW ("SxFileInfo::getPath unexpected result");

      if (file1PathName != SxFileInfo::getName (file1PathStr))
         SX_THROW ("SxFileInfo::getName unexpected result");

      if (file1Path.getRelPath (rootPath) != SxFileInfo::getName (file1PathStr))
         SX_THROW ("SxFileInfo::getRelPath unexpected result");

      if (!SxFileInfo::isAbs (file1PathStr))
         SX_THROW ("SxFileInfo::isAbs unexpected result");
      if (SxFileInfo::isAbs (file1Path.getName ()))
         SX_THROW ("SxFileInfo::isAbs unexpected result");

      if (file1Path.isSymLink ())
         SX_THROW ("SxFileInfo::isSymLink unexpected result");

      // Unified test
      SxFS_FileInfoUnit (file1Path, file2Path);

      // --- Directory

      SxString dir1PathName = "folder1";
      SxString dir1PathStr = rootPath + dir1PathName;
      SxFileInfo dir1Path = SxFileInfo (dir1PathStr);
      // File should not exists yet
      if (dir1Path.exists ())
         SX_THROW ("SxFileInfo::exists unexpected result.");

      // Create folder1
      SxFSAction::mkdir (dir1PathStr);
      SxFSAction::cp (file1PathStr, dir1PathStr + "/file1");
      SxTime::msleep (fsFlushDelayMS);

      SxString dir2PathStr = rootPath + "folder2";
      SxFSAction::mkdir (dir2PathStr);
      SxFSAction::cp (file2PathStr, dir2PathStr + "/file2");
      SxFileInfo dir2Path = SxFileInfo (dir2PathStr);

      // Some general path tests
      if (!SxFileInfo::equals (dir1Path.getAbsPath (), dir1PathStr))
         SX_THROW ("SxFileInfo::getAbsPath unexpected result");

      if (rootPath != SxFileInfo::getPath (dir1PathStr) + "/")
         SX_THROW ("SxFileInfo::getPath unexpected result");

      if (dir1PathName != SxFileInfo::getName (dir1PathStr))
         SX_THROW ("SxFileInfo::getName unexpected result");

      if (dir1Path.getRelPath (rootPath) != SxFileInfo::getName (dir1PathStr))
         SX_THROW ("SxFileInfo::getRelPath unexpected result");

      if (!SxFileInfo::isAbs (dir1PathStr))
         SX_THROW ("SxFileInfo::isAbs unexpected result");
      if (SxFileInfo::isAbs (dir1Path.getName ()))
         SX_THROW ("SxFileInfo::isAbs unexpected result");

      if (dir1Path.isSymLink ())
         SX_THROW ("SxFileInfo::isSymLink unexpected result");

      // Unified test
      SxFS_FileInfoUnit (dir1Path, dir2Path);

#     ifndef WIN32

         // --- SymLink to file
         SxString lnk1Str = rootPath + "lnk1";
         SxFileInfo lnk1File (lnk1Str);
         if (lnk1File.exists ())
            SX_THROW ("SxFileInfo::exists unexpected result.");
         SxFSAction::ln_sf (file1Path, lnk1Str);
         SxTime::msleep (fsFlushDelayMS);

         SxString lnk2Str = rootPath + "lnk2";
         SxFSAction::ln_sf (file2Path, lnk2Str);
         SxFileInfo lnk2File (lnk2Str);
         if (!lnk2File.exists ())
            SX_THROW ("SxFileInfo::exists unexpected result.");
         if (!lnk2File.isSymLink ())
            SX_THROW ("SxFileInfo::isSymLink unexpected result");

         // Reset access timestamp on the files before the test
         SxString::read (file1PathStr);
         SxTime::msleep (fsFlushDelayMS);
         SxString::read (file2PathStr);

         SxFS_FileInfoUnit (lnk1File, lnk2File);


         // --- SymLink to folder
         SxString lnkDir1Str = rootPath + "lnkdir1";
         SxFileInfo lnkDir1File (lnkDir1Str);
         if (lnkDir1File.exists ())
            SX_THROW ("SxFileInfo::exists unexpected result.");
         SxFSAction::ln_sf (dir1PathStr, lnkDir1Str);
         SxFileInfo (dir1PathStr).getSize ();
         SxTime::msleep (fsFlushDelayMS);

         SxString lnkDir2Str = rootPath + "lnkdir2";
         SxFSAction::ln_sf (dir2PathStr, lnkDir2Str);
         SxFileInfo lnkDir2File (lnkDir2Str);
         if (!lnkDir2File.isSymLink ())
            SX_THROW ("SxFileInfo::isSymLink unexpected result");

         SxFS_FileInfoUnit (lnkDir1File, lnkDir2File);
#     endif

      if (!SxFileInfo::equals (SxFileInfo::resolvePath (file1PathStr), file1PathStr))
         SX_THROW ("SxFileInfo::resolvePath unexpected result");

      // Get name of root folder.
      SxString currentFolder = SxFileInfo::getName (rootPath);
      SxString resolvePathStr = SxFileInfo::resolvePath (rootPath + "../" + currentFolder) + "/";

      if (!SxFileInfo::equals (resolvePathStr, rootPath))
         SX_THROW ("SxFileInfo::resolvePath unexpected result");

      // Replace / with //, and \\ with '\\\\'
      SxString sepStr = SxFileInfo::getSeparator ();
      SxString dirtyPath = file1PathStr.substitute (sepStr, sepStr + sepStr);
      if (!SxFileInfo::equals (SxFileInfo::simplifyPath (dirtyPath), file1PathStr))
         SX_THROW ("SxFileInfo::simplifyPath unexpected result");

      SxString origPath = SxFSAction::pwd ().getAbsPath ();
      SxFSAction::cd (rootPath);
      SxFileInfo localFile (file1PathName);

      if (!SxFileInfo::equals (file1Path, localFile))
         SX_THROW ("SxFileInfo local file is not recognized correctly");

      if (  !SxFileInfo::equals (localFile.getOrig (), file1PathName)
         || !SxFileInfo::equals (file1Path.getOrig (), file1PathStr))
         SX_THROW ("SxFileInfo::getOrig unexpected result");

      SxFSAction::cd (origPath);

      // Not implemented. Returns ["file://" + getAbsPath ();]
      //SxString getAbsURI () const;

      // Only one user and group in docker
      //sxuid_t getUID () const;
      //sxgid_t getGID () const

   } catch (SxException e) {
      SX_THROW (e, "SxFS_FileInfo Failed.");
   }
}
