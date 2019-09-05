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

// Including header files
#include <SxFSDeleteAction.h>
#include <SxFSExploreAction.h>
#include <SxFSError.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#ifdef WIN32
#  include <shortcut.h>
#  include <windows.h>
#  include <direct.h>
#  include <unknwn.h>
#  include <shlobj.h>
#  include <objbase.h>
#  include <olectl.h>
#  include <initguid.h>
#endif /* WIN32 */

// --- Constructors
// Constructor
SxFSDeleteAction::SxFSDeleteAction ()
{
   SX_DBG_TRACE("");
   // empty
}
// Destructor
SxFSDeleteAction::~SxFSDeleteAction ()
{
   SX_DBG_TRACE("");
   // empty
}

// --- Class-functions
void SxFSDeleteAction::rm (const SxFile &target)
{
   SX_DBG_TRACE("");
   SxString const & targetStr = target.getAbsPath ();
   if (!target.exists ())  {
      SX_THROW ("Can't remove file '" + targetStr + "'. "
                "It does not exist.");
   }
#ifdef WIN32
      int error = 1;
      if (targetStr.isUnicode ())
            error = _wunlink ((LPCWSTR)targetStr.utf16 ().elements);
      else  error = _unlink  (targetStr.getElems ());
#else
      int error = unlink (targetStr.getElems ());
#endif

   if (error)  {
      SX_THROW ("Can't remove file '" + targetStr
                + "'. " + SxFSError::getUnlinkErrMsg ());
   }
   target.setDirty ();
}

void SxFSDeleteAction::rmdir (const SxDir &path)
{
   SX_DBG_TRACE("");
   SxString const &pathStr (path.getAbsPath ());
   // Testing if the directory that is to be deleted exists
   if (path.exists ())  {
      // Trying to remove the directory
#     ifdef WIN32
         int ret = 0;
         if (pathStr.isUnicode ())
            ret = _wrmdir ((LPCWSTR)pathStr.utf16 ().elements);
         else
            ret = _rmdir (pathStr.ascii ());
#     else
         int ret = ::rmdir (pathStr.getElems ());
#     endif
      if (ret)  {
         SX_THROW ("Can't remove the directory '"
                   + pathStr + "'. " + SxFSError::getRmdirErrMsg ());
      }
   }  else  {
      SX_THROW ("Can't remove the directory '"
                + pathStr + "'. It does not exist.");
   }
   path.setDirty ();
}

void SxFSDeleteAction::rmSymLink (const SxSymLink &target)
{
   SX_DBG_TRACE("");
   SxString const & targetStr = target.getAbsPath ();

   // Checking whether the target that is to be removed exists and is a
   // symbolic link
   if (!target.exists ())  {
      SX_THROW ("Can't remove '" + targetStr + "'. "
                "There is no symbolic link with such a name.");
   }
   // Initializing the error state so that an exception is thrown if the
   // platform is windows
   int error = -1;
#ifndef WIN32
   // Unlinking the symbolic link
   error = ::unlink (targetStr.getElems ());
#endif /* not WIN32 */
   // Testing if either unlink failed or the platform is windows
   if (error)  {
      SX_THROW ("Can't remove symbolic link '"
                + targetStr + "'. " + SxFSError::getUnlinkErrMsg ());
   }
   target.setDirty ();
}


// --- rm_r ("gcc-481") folder with 84391 files in 4 seconds (CentOS 6.3 vm)
void SxFSDeleteAction::rm_r (const SxFileInfo &target)
{
   SX_DBG_TRACE("");
   if (target.isSymLink ())  {
      SxFSDeleteAction::removeSymLink (target);
   } else if (target.isFile ())  {
      SxFSDeleteAction::removeFile (target);
   } else if (target.isDir ())  {
#ifdef WIN32
      SxList<SxFileInfo> fileList;
      // --- find() 9 seconds for 84391 files (gcc-481 folder)
      fileList << SxFSExploreAction::find ("*", SxDir (target.getAbsPath ()));
      SxList<SxFileInfo>::ConstIterator it;
      for (it = fileList.fromLast (); it != fileList.end (); --it)  {
         if (it->isSymLink ())  {
            SxFSDeleteAction::rmSymLink (*it);
         } else if (it->isFile ())  {
            SxFSDeleteAction::rm (*it);
         } else if (it->isDir ())  {
            SxFSDeleteAction::rmdir (*it);
         }
      }
      SxFSDeleteAction::rmdir (target);
#else

      SxList<SxString> fileList;
      // --- listDir() less than 1 second for 84391 files (gcc-481 folder)
      SxFSExploreAction::listDir (target.getAbsPath (), &fileList);
      SxList<SxString>::ConstIterator it;
      for (it = fileList.fromLast (); it != fileList.end (); --it)  {
         SxFileInfo file(*it);
         if (file.isSymLink ())  {
            SxFSDeleteAction::rmSymLink (file);
         } else if (file.isFile ())  {
            SxFSDeleteAction::rm (file);
         } else if (file.isDir ())  {
            SxFSDeleteAction::rmdir (file);
         }
      }
      // Remove the root folder
      SxFSDeleteAction::rmdir (target);
#endif /* WIN32 */
   } else  {
      SX_THROW ("Can't remove '" + target.getAbsPath () +
                "'. A file, symbolic link or directory with "
                "this name does not exist.");
   }
}

#if 0
// --- rm_r "gcc-481" folder with 84391 files in 36 seconds
//     many calls for regexp(pattern) in SxFSExploreAction
void SxFSDeleteAction::rm_r (const SxFileInfo &target)
{
   SX_DBG_TRACE("");
   SxString const &targetStr = target.getAbsPath ();
   if (target.isSymLink ())  {
      SxFSDeleteAction::removeSymLink (target);
   } else  {
      if (target.isFile ())  {
         SxFSDeleteAction::removeFile (target);
      } else if (target.isDir ())  {
         SxFSDeleteAction::removeDir (target);
      } else  {
         SX_THROW ("Can't remove '" + targetStr +
                   "'. A file, symbolic link or directory with "
                   "this name does not exist.");
      }
   }
}

void SxFSDeleteAction::removeDir (const SxFileInfo &target)
{
   SX_DBG_TRACE("");
   SxList<SxFileInfo> fileInfos = SxFSExploreAction::getFileInfos (target);
   SxList<SxFileInfo>::Iterator fileInfoIt;
   SxSymLink targetLink;
   SxFile targetFile;
   SxDir targetDir;
   for (fileInfoIt = fileInfos.begin ();
        fileInfoIt != fileInfos.end ();
        ++fileInfoIt)
   {
      SxFileInfo &curFileInfo = *fileInfoIt;
      SX_CHECK (curFileInfo.getName () != "." &&
                curFileInfo.getName () != "..");
      if (curFileInfo.isSymLink ())  {
         // the found SxFileInfo-object resembles a symbolic link

         // --- Removing symbolic links that belong to this directory
         targetLink = (target / curFileInfo.getName());
         SxFSDeleteAction::removeSymLink (targetLink);

      } else if (curFileInfo.isFile ())  {
         // the found SxFileInfo-object represents a file

         // --- Removing files which are contained in this directory
         targetFile = (target / curFileInfo.getName());
         SxFSDeleteAction::removeFile (targetFile);

      } else if (curFileInfo.isDir ())  {
         // the found SxFileInfo-object is a directory

         // --- Removing directories that are part of this directory
         targetDir = (target / curFileInfo.getName());
         SxFSDeleteAction::removeDir (targetDir);

      } else  {
         SX_EXIT;
      }
   }

   // --- Deleting the by now empty directory
   SxFSDeleteAction::rmdir (target);
}
#endif /* 0 */

void SxFSDeleteAction::removeFile (const SxFileInfo &target)
{
   SX_DBG_TRACE("");
   SxFSDeleteAction::rm (target);
}

void SxFSDeleteAction::removeSymLink (const SxFileInfo &target)
{
   SX_DBG_TRACE("");
   SxFSDeleteAction::rmSymLink (target);
}
