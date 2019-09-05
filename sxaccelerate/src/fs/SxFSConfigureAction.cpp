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
#include <SxFSConfigureAction.h>

#include <SxFSAuthorAction.h>
#include <SxFSError.h>

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
#  include <utime.h>
#endif /* WIN32 */

// --- Constructors
// Constructor
SxFSConfigureAction::SxFSConfigureAction ()
{
   // empty
}
// Destructor
SxFSConfigureAction::~SxFSConfigureAction ()
{
   // empty
}

// --- Class-functions
void SxFSConfigureAction::syncUIDAndGIDAndTimes (const SxFileInfo &src,
                                                 const SxFileInfo &dest)
{
#ifndef WIN32
   // Assuring that the user and group IDs match
   unsigned int srcUID;
   if ((srcUID = src.getUID ()) != dest.getUID ())  {
      try {
         SxFSAuthorAction::chown (dest, srcUID);
      } catch (SxException)  {
         // Willingly ignoring thrown exceptions to circumvent unimportant,
         // irritating and disturbing error messages from the user's point of
         // view
      }
   }
   unsigned int srcGID;
   if ((srcGID = src.getGID ()) != dest.getGID ())  {
      try {
         SxFSAuthorAction::chgrp (dest, srcGID);
      } catch (SxException)  {
         // Willingly ignoring thrown exceptions to circumvent unimportant,
         // irritating and disturbing error messages from the user's point of
         // view
      }
   }
#endif /* not WIN32 */
   // Synchronizing the access and modification times
   SxFSConfigureAction::syncTimes (src, dest);
}

void SxFSConfigureAction::syncTimes (const SxFileInfo &src,
                                     const SxFileInfo &dest)
{
   SxString const &srcStr = src.getAbsPath ();
   SxString const &destStr = dest.getAbsPath ();
#ifdef WIN32
   // --- Fetching the times of the source
   FILETIME srcTime;
   // Retrieving the file handle for source
   HANDLE srcFileHandle = INVALID_HANDLE_VALUE;

   if (srcStr.isUnicode ())
        srcFileHandle = CreateFileW ((LPCWSTR)srcStr.utf16 ().elements,
                                     0, 0, NULL,
                                     OPEN_EXISTING,
                                     ((src.isDir ())?
                                     FILE_FLAG_BACKUP_SEMANTICS :
                                     FILE_ATTRIBUTE_NORMAL),
                                     NULL);
   else srcFileHandle = CreateFileA (srcStr.getElems (),
                                     0, 0, NULL,
                                     OPEN_EXISTING,
                                     ((src.isDir ())?
                                     FILE_FLAG_BACKUP_SEMANTICS :
                                     FILE_ATTRIBUTE_NORMAL),
                                     NULL);
   // Checking whether the file handle could be opened
   if (srcFileHandle != INVALID_HANDLE_VALUE)  {

      if (GetFileTime (srcFileHandle, &srcTime, &srcTime, &srcTime))  {
         // Closing the file handle
         CloseHandle (srcFileHandle);

         // Retrieving the file handle for destination
         HANDLE destFileHandle = INVALID_HANDLE_VALUE;
         if (destStr.isUnicode ())
              destFileHandle = CreateFileW ((LPCWSTR)destStr.utf16 ().elements,
                                            0, 0, NULL,
                                            OPEN_EXISTING,
                                            ((dest.isDir ())?
                                            FILE_FLAG_BACKUP_SEMANTICS :
                                            FILE_ATTRIBUTE_NORMAL),
                                            NULL);
         else destFileHandle = CreateFileA (destStr.getElems (),
                                            0, 0, NULL,
                                            OPEN_EXISTING,
                                            ((dest.isDir ())?
                                            FILE_FLAG_BACKUP_SEMANTICS :
                                            FILE_ATTRIBUTE_NORMAL),
                                            NULL);

         // Checking whether the file handle could be opened
         if (destFileHandle != INVALID_HANDLE_VALUE)  {
            // --- Writing the times to the destination
            SetFileTime (destFileHandle, &srcTime, &srcTime, &srcTime);
            // Closing the file handle
            CloseHandle (destFileHandle);
         } else  {
            // Closing the file handle
            CloseHandle (srcFileHandle);
            SX_THROW ("Can't synchronize times of '"
                      + srcStr + "' and '" + destStr + "'. "
                        "Can't get the file handle of '" +destStr+ "'.");
         }
      } else  {
         // Closing the file handle
         CloseHandle (srcFileHandle);
         SX_THROW ("Can't synchronize times of '"
                   + srcStr + "' and '" + destStr + "'. "
                     "Can't retrieve the times of '" + srcStr + "'.");
      }
   } else  {
      SX_THROW ("Can't synchronize times of '"
                + srcStr + "' and '" + destStr + "'. "
                  "Can't get the file handle of '" + srcStr + "'.");
   }
#else
   // --- Fetching the times of the source
/*   struct stat statbuf;
   if (stat (srcStr.getElems (), &statbuf) != 0)  {
      SX_THROW (("Error: Can't synchronize times of '"
                     + srcStr + "' and '" + destStr + "'."
                     + SxFSError::getUtimeErrMsg ()
                    ).getElems (),
                    __FILE__, __LINE__);
   }*/
   // --- Writing the times to the destination
   struct utimbuf srcTime;
   srcTime.actime = (time_t)src.lastAccessed ();//statbuf.st_atime;
   srcTime.modtime = (time_t)src.lastModified ();//statbuf.st_mtime;
   int error = utime (destStr.getElems (), &srcTime);
   if (error)  {
      SX_THROW ("Can't synchronize times of '"
                + srcStr + "' and '" + destStr + "'. "
                + SxFSError::getUtimeErrMsg ());
   }
#endif /* WIN32 */
   dest.setDirty ();
}
