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
#include <SxFSCreateAction.h>
#include <SxFSExploreAction.h>
#include <SxFSError.h>
#include <SxTime.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#ifdef WIN32
#  include <shortcut.h>
#  include <windows.h>
#  include <unknwn.h>
#  include <shlobj.h>
#  include <objbase.h>
#  include <olectl.h>
#  include <initguid.h>
#else
#  include <dirent.h>
#endif /* WIN32 */


// --- Constructors
// Constructor
SxFSCreateAction::SxFSCreateAction ()
{
   SX_DBG_TRACE("");
   // empty
}
// Destructor
SxFSCreateAction::~SxFSCreateAction ()
{
   SX_DBG_TRACE("");
   // empty
}

void SxFSCreateAction::touch (const SxFileInfo &target, int mode)
{
   SX_DBG_TRACE("");
   SxString const &targetStr (target.getAbsPath ());

   // Checking if the target element is a directory or a file
   if (target.isDir ())  {
      // Updating the timestamp of the directory
      SxFSCreateAction::updateDir (SxFileInfo (targetStr), mode);
   } else  {
      // Updating the timestamp of the file or even creating it if it does not
      // exist
      SxFSCreateAction::updateFile (SxFileInfo (targetStr), mode);
   }
   target.setDirty ();
}

void SxFSCreateAction::mkdir (const SxDir &dir, int umask)
{
   SX_DBG_TRACE("");
   SxString dirStr = dir.getAbsPath ();
   if (dir.exists ())  {
      SX_THROW ("Can't create directory '" + dirStr
                 + "'. File already exists.");
   }

#  ifdef WIN32
   int error = true;
   if (dirStr.isUnicode ())
        error = !CreateDirectoryW ((LPCWSTR)dirStr.utf16 ().elements, NULL);
   else error = !CreateDirectoryA (dirStr.getElems (), NULL);
#  else
   int error = ::mkdir (dirStr.getElems (), (mode_t)umask);
#  endif /* WIN32 */

   if (error)  {
      SX_THROW ("Can't create directory '"+ dirStr +"': "+ sxstrerror());
   }
   dir.setDirty ();
}

void SxFSCreateAction::mkdir_p (const SxDir &destDir, int umask)
{
   SX_DBG_TRACE("");
   SxString destDirStr = destDir.getAbsPath ();
   SxArray<SxString> subfolders = destDirStr.tokenize ('/');

   if (subfolders.getSize () > 0)  {
      ssize_t i=0;
#      ifdef WIN32
          SxDir path(subfolders(i++)); // 'C:'
#      else
          SxDir path("/");
#      endif /* WIN32 */

      for (; i < subfolders.getSize (); i++)  {
         path = path / subfolders(i);
         if (!path.exists ())
            SxFSCreateAction::mkdir (path, umask);
      }
   }
}

SxFileInfo SxFSCreateAction::ln_sf (const SxFileInfo &path,
                                    const SxFileInfo &newLink)
{
   SX_DBG_TRACE("");
   // Willingly not testing whether the path exists or not to enable relative
   // and virtual paths
   SxString const & pathStr = path.getOrig ();
   SxString const & newLinkStr = newLink.getAbsPath ();
#  ifndef WIN32
      int error = -1;
      error = symlink (pathStr.getElems (), newLinkStr.getElems ());
      if (error)  {
         SX_THROW ("Can't create symlink '"
                   + newLinkStr + "' -> '"
                   + pathStr + "'. " +
                   SxFSError::getSymlinkErrMsg ());
      }
#  else
      //createShortcut (newLinkStr, pathStr);
      //DWORD flag = 0;
      //if (CreateSymbolicLink (newLinkStr.getElems (), pathStr.getElems (), flag) == 0)  {
      //   SX_THROW ("Can't create link '" + newLinkStr + "' -> '"
      //             + pathStr + "': " + sxstrerror ());
      //}
      // error: requires user CreateSymbolicLinkPrivilege
      // error: The file or directory is not a reparse point.
      SX_EXIT; // not yet implemented
#  endif /* not WIN32 */

   newLink.setDirty ();
   return SxFileInfo (newLink);
}

SxFile SxFSCreateAction::createTmpFile (const SxString &tmpDir,
                                        const SxArray<unsigned char> &buffer)
{
   SX_DBG_TRACE("");
   // --- Trying to find the directory devoted to temporary files
   SxDir folder;
#  ifdef WIN32
      if (tmpDir == SxString ())  folder = SxString(getenv ("TMP"));
#  else
      if (tmpDir == SxString ())  folder = SxString(getenv ("TMPDIR"));
#  endif
   else                        folder = tmpDir;
   if (!folder.exists ())      folder = SxString("/tmp");  // try fallback
   if (!folder.exists ())  {
      SX_THROW ("No valid temporary folder found. Please "
                "provide the location in the TMPDIR variable.");
   }

   // Defining a template for creation of a template file name (the six
   // trailing X's are mandatory!)
   SxFile res = folder.getAbsPath () + "/sxtmpXXXXXX";
   SxString const & resInitStr = res.getAbsPath ();

   // --- Creating a C-like string to enable the mkstemp ()-function to modify
   //     the content of the proposed template
#  ifdef WIN32
      SxArray<uint16_t> str = resInitStr.utf16 ();
#  else
      SxArray<char> str(resInitStr.getNBytes () + 1);
      memcpy (str.elements,
              resInitStr.getElems (),
              static_cast<size_t>(resInitStr.getNBytes () + 1));
#  endif

   // Initializing the file descriptor with an invalid value so that an
   // exception is thrown in case that no mkstemp ()-function exists
   int fd = -1;
#  ifdef HAVE_MKSTEMP
      fd = ::mkstemp (str.elements);
#  elif defined(WIN32)
      if (_wmktemp_s ((LPWSTR)str.elements, str.getSize ()) != 0)  {
         SX_THROW ("Can't create a temporary file '"
                   + res.getAbsPath () + "': _mktemp_s() failed.");
      }
      int mode = 0600;
      fd = ::_wopen ((LPCWSTR)str.elements, O_WRONLY | O_CREAT, mode);

#  endif /* HAVE_MKSTEMP */

   // --- Checking whether a valid temporary file could be generated
   if (fd < 0)  {
      SX_THROW("Can't create a temporary file according "
               "to the following template '"+ res.getAbsPath ()
               + "'. "+ SxFSError::getOpenErrMsg ());
   } else  {

      // Storing the actually used path
      // But maintaing unicode flag if one was set
#     ifdef WIN32
         SxString tempStr = SxString::fromUtf16 (str.elements);
#     else
         SxString tempStr (str.elements);
#     endif
      if (folder.getAbsPath ().isUnicode ())
         tempStr.setUnicode ();
      res = SxFile (tempStr);

      SxString const &resStr = res.getAbsPath ();

      // --- Writing the passed content to the file in one chunk
      ssize_t error = write (fd, buffer.elements, static_cast<size_t>(buffer.getSize ()));
      if (error == -1)  {
         SX_THROW ("Can't write to temporary file '"
                   + resStr + "'. " + SxFSError::getWriteErrMsg ());
      }
      // Closing the file descriptor that was opened by mkstemp ()
      close (fd);
   }

   res.setDirty ();

   // Returning the temporary file that was just created
   return res;
}

void SxFSCreateAction::updateFile (const SxFileInfo &target, int mode)
{
   SX_DBG_TRACE("");
   SxString const &targetStr = target.getAbsPath ();

   // Renewing the time stamp of the file or even creating it if it does not
   // exist
   int fd = -1;
#  ifdef WIN32
      if (targetStr.isUnicode ()) {
         fd = _wopen ((LPCWSTR)targetStr.utf16 ().elements,
                      O_WRONLY | O_CREAT, mode);
      } else {
         fd = _open (targetStr.ascii (), O_WRONLY | O_CREAT, mode);
      }
#  else
      fd = open (targetStr.getElems (), O_WRONLY | O_CREAT, mode);
#  endif

   if (fd == -1)  {
      SX_THROW ("Can't write to file '" +
                targetStr + "'. " + SxFSError::getOpenErrMsg ());
   }
   close (fd);

#  if defined(LINUX)
      // --- update both atime and mtime
      long t = static_cast<long>(SxTime::getRealTime ());
      struct timeval times[2];
      times[0].tv_sec = t;
      times[0].tv_usec = 0;
      times[1].tv_sec = t;
      times[1].tv_usec = 0;
      if (utimes (targetStr.getElems (), times) != 0)  {
         SX_THROW ("utimes('"+targetStr+"', "+SxString(t)+") error: "
                   + sxstrerror());
      }
#  endif /* LINUX */
}

void SxFSCreateAction::updateDir (const SxFileInfo &/*target*/, int)
{
   SX_DBG_TRACE("");
   return;
//   //TODO Check why the timestamp is not updated correctly or implement a
//   //     working solution using utime
//   SxString const &targetStr = target.getAbsPath ();
//   DIR *dir = opendir (targetStr.getElems ());
//   if (!dir)  {
//      SX_THROW (("Can't open directory '" + targetStr
//                     + "'.").getElems (),
//                    __FILE__, __LINE__);
//   }
//   // Closing the opened directory
//   closedir (dir);
}

//#ifdef WIN32
//void SxFSCreateAction::createShortcut (const SxString &lnkSrc,
//                                       const SxString &lnkDest)
//{
//   SX_DBG_TRACE("");
//   int err = ::sxnumlibs_createShortcut (lnkSrc.getElems (), lnkDest.getElems ());
//   SxFSError::throwShortcutException (err, lnkSrc, lnkDest);
//}
//#endif /* WIN32 */
