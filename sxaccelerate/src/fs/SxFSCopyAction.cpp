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
#include <SxFS.h>
#include <SxFSCopyAction.h>
#include <SxFSCreateAction.h>
#include <SxFSDeleteAction.h>
#include <SxFSConfigureAction.h>
#include <SxFSExploreAction.h>
#include <SxFSAuthorAction.h>

#include <SxFile.h>
#include <SxDir.h>

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
#endif /* WIN32 */

// --- Constructors
// Constructor
SxFSCopyAction::SxFSCopyAction ()
{
   SX_DBG_TRACE("");
   // empty
}
// Destructor
SxFSCopyAction::~SxFSCopyAction ()
{
   SX_DBG_TRACE("");
   // empty
}

// --- Methods
void SxFSCopyAction::convFileToFile (const SxFileInfo &src,
                                     const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxFile (src).exists ());
   SX_CHECK (SxFile (dest).exists ());
   SxString const &srcStr = src.getAbsPath ();
   SxString const &destStr = dest.getAbsPath ();

   // Performing a self-test
   if (SxFileInfo::equals (src, dest))  {
      SX_THROW ("Can't copy file '" + srcStr + "' to '" + destStr
                + "'. It is not possible to copy a file to itself.");
   } else  {
      // Performing the actual copy-operation
      SxFSCopyAction::convFileToNothing (src, dest, true);

   }
}

void SxFSCopyAction::convFileToDir (const SxFileInfo &src,
                                    const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxFile (src).exists ());
   SX_CHECK (SxDir (dest).exists ());
   SxString const &srcStr = src.getAbsPath ();
   SxString const &destStr = dest.getAbsPath ();

   // Determining the new destination
   SxFileInfo newDest (dest / src.getName ());

   // --- Checking whether the new destination is
   if (newDest.isSymLink ())  {
      // new destination is a symbolic link
      if (newDest.isFile ())  {
         // that addresses a file
         SxFSCopyAction::convFileToFile (src, newDest);
      } else if (newDest.isDir ())  {
         // that addresses a directory
         SX_THROW ("Can't copy file '" + srcStr + "' to '" + destStr
                   + "'. It is not possible to copy a file to this "
                     "directory as it already contains a directory "
                     "with that name.");
      } else  {
         // that has no valid target
         SX_THROW ("Can't copy file '" + srcStr + "' to '" + destStr
                   + "'. It is not possible to copy a file to this "
                     "directory as it already contains an empty "
                     "symbolic link with that name.");
      }
   } else if (newDest.isFile ())  {
      // new destination is a file
      SxFSCopyAction::convFileToFile (src, newDest);
   } else if (newDest.isDir ())  {
      // new destination is a directory
      SX_THROW ("Can't copy file '" + srcStr + "' to '" + destStr
                + "'. It is not possible to copy a file to this "
                  "directory as it already contains a directory "
                  "with that name.");
   } else  {
      // new destination does not exist
      SxFSCopyAction::convFileToNothing (src, newDest);
   }
}

void SxFSCopyAction::convFileToNothing (const SxFileInfo &src,
                                        const SxFileInfo &dest,
                                        bool force)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxFile (src).exists ());
   if (!force)  {
      SX_CHECK (!SxFile (dest).exists () && !SxDir (dest).exists () &&
                !SxSymLink (dest).exists ());
   }
   SxString const &srcStr = src.getAbsPath ();
   SxString const &destStr = dest.getAbsPath ();
   bool destAlreadyExisted (dest.isFile ());
   // Temporarily setting the umask to 0 so that the mask is copied
   // correctly
   mode_t oldUmask = ::umask (0);
   int const srcMode = src.getMode ();
   // --- Opening the source file for reading
   int fdIn, fdOut;
#  ifdef WIN32
      if (srcStr.isUnicode ())
           fdIn = _wopen ((const wchar_t *)srcStr.utf16 ().elements,
                          _O_RDONLY | _O_BINARY);
      else fdIn =  _open (srcStr.ascii (), _O_RDONLY | _O_BINARY);
#  else
      fdIn = open (srcStr.getElems (), O_RDONLY);
#  endif /* WIN32 */
   if (fdIn == -1)  {
      SX_THROW ("Can't read from file '" + srcStr + "'.");
   }

   // --- Opening the destination file for writing
#  ifdef WIN32
      if (destStr.isUnicode ())
           fdOut = _wopen ((const wchar_t *)destStr.utf16 ().elements,
                           _O_WRONLY | _O_CREAT | _O_TRUNC | _O_BINARY,
                           srcMode);
      else fdOut =  _open (destStr.getElems (),
                           _O_WRONLY | _O_CREAT | _O_TRUNC | _O_BINARY,
                           srcMode);
#  else
      fdOut = open (destStr.getElems (), O_WRONLY | O_CREAT | O_TRUNC,
                    srcMode);
#  endif /* WIN32 */
   if (fdOut == -1)  {
      SX_THROW ("Can't write to file '" + destStr + "'.");
   }

   // --- Creating a buffer for the copy-operation

   SxArray<char> buffer (SxFileInfo::BUFFER_SIZE);

   // --- Copying data
#  ifdef WIN32
      unsigned int bufLen = static_cast<unsigned int>(buffer.getSize ());
      int nr = 0;
      int nw = 0;
      while ( (nr = _read (fdIn, buffer.elements, bufLen)) > 0)  {
         nw = _write (fdOut, buffer.elements,
                      static_cast<unsigned int>(nr));
         if (nw != nr)  {
            SX_THROW ("Can't write to file '" + destStr + "'.");
         }
      }
#  else
      ssize_t n;
      while ( (n = read (fdIn, buffer.elements, SxFileInfo::BUFFER_SIZE)) > 0)  {
         write (fdOut, buffer.elements, static_cast<size_t>(n));
      }
#  endif /* WIN32 */

   // --- Closing the file descriptors again
   close (fdOut);
   close (fdIn);

   // Explicitly setting the mask provided the file was not created
   if (destAlreadyExisted)  {
      SxFSAuthorAction::chmod (srcMode, dest);
   }

   // Resetting the umask
   ::umask (oldUmask);

   // Synchronizing the user ID, the group ID as well as the access
   // and modification times
   SxFSConfigureAction::syncUIDAndGIDAndTimes (src, dest);
}

void SxFSCopyAction::convFileToFileSymLink (const SxFileInfo &src,
                                            const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxFile (src).exists ());
   SX_CHECK (SxSymLink (dest).exists () && SxSymLink (dest).isFile ());
   // Treating the destination as if it was a regular file
   SxFSCopyAction::convFileToFile (src, dest);
}

void SxFSCopyAction::convFileToDirSymLink (const SxFileInfo &src,
                                           const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxFile (src).exists ());
   SX_CHECK (SxSymLink (dest).exists () && SxSymLink (dest).isDir ());
   // Treating the destination as if it was a directory
   SxFSCopyAction::convFileToDir (src, dest);
}

void SxFSCopyAction::convFileToEmptySymLink (const SxFileInfo &src,
                                             const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxFile (src).exists ());
   SX_CHECK (SxSymLink (dest).exists () && !SxSymLink (dest).isValid ());
   SxString const &srcStr = src.getAbsPath ();
   SxString const & destStr = dest.getAbsPath ();
   SX_THROW ("Can't copy file '" + srcStr + "' to '" + destStr
             + "'. It is not possible to copy a file to an "
               "empty symbolic link.");
}


void SxFSCopyAction::convDirToFile (const SxFileInfo &src,
                                    const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxDir (src).exists ());
   SX_CHECK (SxFile (dest).exists ());
   SxString const &srcStr = src.getAbsPath ();
   SxString const & destStr = dest.getAbsPath ();
   SX_THROW ("Can't copy directory '" + srcStr + "' to '" + destStr
             + "'. It is not possible to copy a directory to a file.");
}

void SxFSCopyAction::convDirToDir (const SxFileInfo &src,
                                   const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxDir (src).exists ());
   SX_CHECK (SxDir (dest).exists ());
   SxString const & srcDirStr = src.getAbsPath ();
   SxString const & destDirStr = dest.getAbsPath ();
   //std::cout << "void SxFSCopyAction::convDirToDir ('" << srcDirStr;//TEST
   //std::cout << "', '" <<  destDirStr << "')";//TEST
   //std::cout << std::endl;//TEST
   if (SxFileInfo::equals (src, dest))  {
      SX_THROW ("Can't copy directory '"
                + srcDirStr + "' to '" + destDirStr
                + "'. It is not possible to copy a directory to itself.");
   } else  {
      SxFSCopyAction::convDirToNothing (src, (dest / src.getName()), true);
   }
}

void SxFSCopyAction::convDirToNothing (const SxFileInfo &src,
                                       const SxFileInfo &dest,
                                       bool force)
{
   SX_DBG_TRACE("");

   //std::cout << "void SxFSCopyAction::convDirToNothing ('" << src.getAbsPath ();//TEST
   //std::cout << "', '" <<  dest.getAbsPath () << "')";//TEST
   //std::cout << std::endl;//TEST

   SX_CHECK (SxDir (src).exists ());
   if (!force)  {
      SX_CHECK (!SxFile (dest).exists () && !SxDir (dest).exists () &&
                !SxSymLink (dest).exists ());
   }
   SxDir srcDir (src);
   SxDir destDir (dest);

   // Temporarily setting the umask to 0 so that the mask is copied
   // correctly
   mode_t oldUmask = ::umask (0);

   if (!destDir.exists ())  {
      SxFSCreateAction::mkdir (destDir, srcDir.getMode ());
   }

   SxFileInfo curDestLink;

   // --- Copying invalid symbolic links contained in this directory
   SxList<SxSymLink> symLinks = SxFSExploreAction::getSymLinks (src);
   SxList<SxSymLink>::Iterator srcSymLinkIt;
   for (srcSymLinkIt = symLinks.begin ();
        srcSymLinkIt != symLinks.end ();
        ++srcSymLinkIt)
   {
      SxSymLink &curSrcSymLink = *srcSymLinkIt;
      if (!curSrcSymLink.isFile () && !curSrcSymLink.isDir ())  {
         //std::cout << "OH NO 1" << std::endl;//TEST
         curDestLink = (destDir / curSrcSymLink.getName());
         SxFSCopyAction::convEmptySymLinkToNothing (curSrcSymLink,
                                                    curDestLink,
                                                    true);
      }
   }

   // --- Copying files contained in this directory
   SxFile curDestFile;
   SxList<SxFile> files = SxFSExploreAction::getFiles (src);
   SxList<SxFile>::Iterator srcFileIt;
   for (srcFileIt = files.begin ();
        srcFileIt != files.end ();
        ++srcFileIt)
   {
      SxFile &curSrcFile = *srcFileIt;
      // Testing if the current file is really a symbolic link
      if (curSrcFile.isSymLink ())   {
         curDestLink = (destDir / curSrcFile.getName());
         //std::cout << "OH NO 2" << std::endl;//TEST
         SxFSCopyAction::convFileSymLinkToNothing (curSrcFile, curDestLink,
                                                   true);
      } else  {
         curDestFile = SxFile (destDir / curSrcFile.getName());
         //std::cout << "OH NO 3" << std::endl;//TEST
         SxFSCopyAction::convFileToNothing (curSrcFile, curDestFile, true);
      }
   }

   // --- Copying directories that are part of this directory
   SxDir curDestDir;
   //std::cout << "WHERE AM I LOOKING '" << src.getAbsPath ();//TEST
   //std::cout << "'" << std::endl;//TEST
   SxList<SxDir> dirs = SxFSExploreAction::getDirs (src);
   SxList<SxDir>::Iterator srcDirIt;
   for (srcDirIt = dirs.begin (); srcDirIt != dirs.end (); ++srcDirIt)
   {
      SxDir &curSrcDir = *srcDirIt;
         //std::cout << "VERY BAD '" << curSrcDir.getAbsPath ();//TEST
         //std::cout << "'" << std::endl;//TEST
      // Testing if the current directory is really a symbolic link
      if (curSrcDir.isSymLink ())  {
         //std::cout << "OH NO 4" << std::endl;//TEST
         curDestDir = destDir;
         SxFSCopyAction::convDirSymLinkToDir (curSrcDir, curDestDir);
      } else  {
         //std::cout << "OH NO 5" << std::endl;//TEST
         curDestDir = SxDir (destDir / curSrcDir.getName());
         SxFSCopyAction::convDirToNothing (curSrcDir, curDestDir, true);
      }
   }

   // Resetting the umask
   ::umask (oldUmask);

   // Synchronizing the user ID, the group ID as well as the access
   // and modification times
   SxFSConfigureAction::syncUIDAndGIDAndTimes (src, dest);
}

void SxFSCopyAction::convDirToFileSymLink (const SxFileInfo &src,
                                           const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxDir (src).exists ());
   SX_CHECK (SxSymLink (dest).exists () && SxSymLink (dest).isFile ());
   SxString const &srcStr = src.getAbsPath ();
   SxString const & destStr = dest.getAbsPath ();
   SX_THROW ("Can't copy directory '"
             + srcStr + "' to '" + destStr
             + "'. It is not possible to copy a directory to a "
               "symbolic link that addresses a file.");
}

void SxFSCopyAction::convDirToDirSymLink (const SxFileInfo &src,
                                          const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxDir (src).exists ());
   SX_CHECK (SxSymLink (dest).exists () && SxSymLink (dest).isDir ());
   // Determining the new destination
   SxFileInfo newDest (dest / src.getName ());
   // Performing the actual copy-operation
   SxFSCopyAction::convDirToNothing (src, newDest);
}

void SxFSCopyAction::convDirToEmptySymLink (const SxFileInfo &src,
                                            const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxDir (src).exists ());
   SX_CHECK (SxSymLink (dest).exists () && !SxSymLink (dest).isValid ());
   SxString const &srcStr = src.getAbsPath ();
   SxString const & destStr = dest.getAbsPath ();
   SX_THROW ("Can't copy directory '"
             + srcStr + "' to '" + destStr
             + "'. It is not possible to copy a directory to a "
               "symbolic link that addresses nothing.");
}


void SxFSCopyAction::convFileSymLinkToFile (const SxFileInfo &src,
                                            const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && SxSymLink (src).isFile ());
   SX_CHECK (SxFile (dest).exists ());
   // Removing the existing destination file
   SxFSDeleteAction::rm (SxFile (dest));
   // Performing the actual copy-operation on the link
   SxFSCopyAction::convFileSymLinkToNothing (src, dest);
}

void SxFSCopyAction::convFileSymLinkToDir (const SxFileInfo &src,
                                           const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && SxSymLink (src).isFile ());
   SX_CHECK (SxDir (dest).exists ());

   // Determining the new destination
   SxFileInfo newDest (dest/src.getName ());
   SxString const &srcStr = src.getAbsPath ();
   SxString const & newDestStr = newDest.getAbsPath ();

   // Testing how and whether to create the new empty symbolic link
   if (newDest.isSymLink ())  {
      // the new destination is a symbolic link

      if (newDest.isFile ())  {
         // that addresses a file

         SxFSCopyAction::convFileSymLinkToFileSymLink (src, newDest);
      } else if (newDest.isDir ())  {
         // that addresses a directory

         // --- Removing the symbolic link to a directory that exists at the
         //     new destination and inserting a link with the same target as
         //     the source symbolic link
         SxFSDeleteAction::removeSymLink (newDest);
         convFileSymLinkToNothing (src, newDest);
      } else {
         // that addresses nothing

         SxFSCopyAction::convFileSymLinkToEmptySymLink (src, newDest);
      }
   } else if (newDest.isFile ())  {
      // the new destination is a file

      // Replacing the existing file at the new destination by a symbolic link
      // that shares its target with src
      SxFSCopyAction::convFileSymLinkToFile (src, newDest);
   } else if (newDest.isDir ())  {
      // the new destination is a directory
      SX_THROW ("Can't copy '"
                + srcStr + "' to '" + newDestStr
                + "'. It is not possible to copy a symbolic "
                  "link that addresses a file to a directory that "
                  "contains a directory with the same name.");
   } else {
      // the new destination does not exist

      // Performing the actual copy-operation on the link
      SxFSCopyAction::convFileSymLinkToNothing (src, newDest);
   }
}

void SxFSCopyAction::convFileSymLinkToNothing (const SxFileInfo &src,
                                               const SxFileInfo &dest,
                                               bool force)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && SxSymLink (src).isFile ());
   if (!force)  {
      SX_CHECK (!SxFile (dest).exists () && !SxDir (dest).exists () &&
                !SxSymLink (dest).exists ());
   }
   SxSymLink srcLink (src.getOrig ());
   // Getting the link target
   SxString target;
   target = srcLink.getTarget ();
   // Creating a link that has the same target as the source link
   SxFSCreateAction::ln_sf (SxFileInfo (target), dest);
}

void SxFSCopyAction::convFileSymLinkToFileSymLink (const SxFileInfo &src,
                                                   const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && SxSymLink (src).isFile ());
   SX_CHECK (SxSymLink (dest).exists () && SxSymLink (dest).isFile ());
   SxString const &srcStr = src.getAbsPath ();
   SxString const &destStr = dest.getAbsPath ();
   // Performing a self-test
   if (!SxFileInfo::equals (src, dest))  {
      // Removing the existing destination link
      SxFSDeleteAction::removeSymLink (dest);
      SxFSCopyAction::convFileSymLinkToNothing (src, dest);
   } else  {
      SX_THROW ("Can't copy the symbolic link '"
                + srcStr + "' to '" + destStr
                + "'. It is not possible to copy a symbolic link "
                  "that points to a file to itself.");
   }
}

void SxFSCopyAction::convFileSymLinkToDirSymLink (const SxFileInfo &src,
                                                  const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && SxSymLink (src).isFile ());
   SX_CHECK (SxSymLink (dest).exists () && SxSymLink (dest).isDir ());

   // Treating the destination as if it was a directory
   SxFSCopyAction::convFileSymLinkToDir (src, dest);
}

void SxFSCopyAction::convFileSymLinkToEmptySymLink (const SxFileInfo &src,
                                                    const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && SxSymLink (src).isFile ());
   SX_CHECK (SxSymLink (dest).exists () && !SxSymLink (dest).isValid ());
   // Removing the existing destination link
   SxFSDeleteAction::removeSymLink (dest);
   SxFSCopyAction::convFileSymLinkToNothing (src, dest);
}

void SxFSCopyAction::convDirSymLinkToFile (const SxFileInfo &src,
                                           const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && SxSymLink (src).isDir ());
   SX_CHECK (SxFile (dest).exists ());
   // Removing the existing destination file
   SxFSDeleteAction::removeFile (dest);
   // Performing the actual copy-operation on the link
   SxFSCopyAction::convDirSymLinkToNothing (src, dest);
}

void SxFSCopyAction::convDirSymLinkToDir (const SxFileInfo &src,
                                          const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && SxSymLink (src).isDir ());
   SX_CHECK (SxDir (dest).exists ());
   // Determining the new destination
   SxFileInfo newDest (dest/src.getName ());
   // --- Removing the new destination if there is one
   // Performing the actual copy-operation on the link
   SxFSCopyAction::convDirSymLinkToNothing (src, newDest);
}

void SxFSCopyAction::convDirSymLinkToNothing (const SxFileInfo &src,
                                              const SxFileInfo &dest,
                                              bool force)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && SxSymLink (src).isDir ());
   if (!force)  {
      SX_CHECK (!SxFile (dest).exists () && !SxDir (dest).exists () &&
                !SxSymLink (dest).exists ());
   }
   SxSymLink srcLink (src.getOrig ());
   // Getting the link target
   SxString target;
   target = srcLink.getTarget ();
   // Creating a link that has the same target as the source link
   SxFSCreateAction::ln_sf (SxFileInfo (target), dest);
}

void SxFSCopyAction::convDirSymLinkToFileSymLink (const SxFileInfo &src,
                                                  const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && SxSymLink (src).isDir ());
   SX_CHECK (SxSymLink (dest).exists () && SxSymLink (dest).isFile ());
   // Removing the destination link
   SxFSDeleteAction::removeSymLink (dest);
   // Performing the actual copy-operation on the link
   SxFSCopyAction::convDirSymLinkToNothing (src, dest);
}

void SxFSCopyAction::convDirSymLinkToDirSymLink (const SxFileInfo &src,
                                                 const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && SxSymLink (src).isDir ());
   SX_CHECK (SxSymLink (dest).exists () && SxSymLink (dest).isDir ());
   SxFSCopyAction::convDirSymLinkToDir (src, dest);
}

void SxFSCopyAction::convDirSymLinkToEmptySymLink (const SxFileInfo &src,
                                                   const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && SxSymLink (src).isDir ());
   SX_CHECK (SxSymLink (dest).exists () && !SxSymLink (dest).isValid ());
   // Removing the destination link
   SxFSDeleteAction::removeSymLink (dest);
   SxFSCopyAction::convDirSymLinkToNothing (src, dest);
}

void SxFSCopyAction::convEmptySymLinkToFile (const SxFileInfo &src,
                                             const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   // Checking whether the preconditions are fullfilled
   SX_CHECK (SxSymLink (src).exists () && !SxSymLink (src).isValid ());
   SX_CHECK (SxFile (dest).exists ());
   // Removing the existing destination file
   SxFSDeleteAction::removeFile (dest);
   // Performing the actual copy-operation on the link
   SxFSCopyAction::convEmptySymLinkToNothing (src, dest);
}

void SxFSCopyAction::convEmptySymLinkToDir (const SxFileInfo &src,
                                            const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   // Checking whether the preconditions are fullfilled
   SX_CHECK (SxSymLink (src).exists () && !SxSymLink (src).isValid ());
   SX_CHECK (SxDir (dest).exists ());
   // Determining the new destination
   SxFileInfo newDest (dest/src.getName ());
   SxString const &srcStr = src.getAbsPath ();
   SxString const & newDestStr = newDest.getAbsPath ();

   // Testing how and whether to create the new empty symbolic link
   if (newDest.isSymLink ())  {
      // the new destination is a symbolic link

      if (newDest.isFile ())  {
         // that addresses a file

         // Performing the actual copy-operation on the link
         SxFSCopyAction::convEmptySymLinkToFileSymLink (src, newDest);
      } else if (newDest.isDir ())  {
         // that addresses a directory

         // Removing the destination link
         SxFSDeleteAction::removeSymLink (newDest);
         // Performing the actual copy-operation on the link
         SxFSCopyAction::convEmptySymLinkToNothing (src, newDest);

      } else {
         // that addresses nothing

         // Performing the actual copy-operation on the link
         SxFSCopyAction::convEmptySymLinkToEmptySymLink (src, newDest);
      }
   } else if (newDest.isFile ())  {
      // the new destination is a file

      // Performing the actual copy-operation on the link
      SxFSCopyAction::convEmptySymLinkToFile (src, newDest);
   } else if (newDest.isDir ())  {
      // the new destination is a directory

      SX_THROW ("Can't copy '" + srcStr + "' to '" + newDestStr
                + "'. It is not possible to copy an empty "
                  "symbolic link to a directory that contains "
                  "a directory with the same name.");
   } else {
      // the new destination does not exist

      // Performing the actual copy-operation on the link
      SxFSCopyAction::convEmptySymLinkToNothing (src, newDest);
   }
}

void SxFSCopyAction::convEmptySymLinkToNothing (const SxFileInfo &src,
                                                const SxFileInfo &dest,
                                                bool force)
{
   SX_DBG_TRACE("");
   // Checking whether the preconditions are fullfilled
   SX_CHECK (SxSymLink (src).exists () && !SxSymLink (src).isValid ());
   if (!force)  {
      SX_CHECK (!SxFile (dest).exists () && !SxDir (dest).exists () &&
                !SxSymLink (dest).exists ());
   }
   SxSymLink srcLink (src.getOrig ());
   // Getting the link target
   SxString target;
   target = srcLink.getTarget ();

   // Creating a link that has the same target as the source link
   SxFSCreateAction::ln_sf (SxFileInfo (target), dest);

}

void SxFSCopyAction::convEmptySymLinkToFileSymLink (const SxFileInfo &src,
                                                    const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   // Checking whether the preconditions are fullfilled
   SX_CHECK (SxSymLink (src).exists () && !SxSymLink (src).isValid ());
   SX_CHECK (SxSymLink (dest).exists () && SxSymLink (dest).isFile ());
   rmDestSymLinkAndCpEmptySymLink (src, dest);
}

void SxFSCopyAction::convEmptySymLinkToDirSymLink (const SxFileInfo &src,
                                                   const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   // Checking whether the preconditions are fullfilled
   SX_CHECK (SxSymLink (src).exists () && !SxSymLink (src).isValid ());
   SX_CHECK (SxSymLink (dest).exists () && SxSymLink (dest).isDir ());
   // Treating the destination as if it was a directory
   SxFSCopyAction::convEmptySymLinkToDir (src, dest);
}

void SxFSCopyAction::convEmptySymLinkToEmptySymLink (const SxFileInfo &src,
                                                     const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   // Checking whether the preconditions are fullfilled
   SX_CHECK (SxSymLink (src).exists () && !SxSymLink (src).isValid ());
   SX_CHECK (SxSymLink (dest).exists () && !SxSymLink (dest).isValid ());
   SxString const & srcStr = src.getAbsPath ();
   SxString const & destStr = dest.getAbsPath ();
   // Performing a self-test
   if (!SxFileInfo::equals (src, dest))  {
      rmDestSymLinkAndCpEmptySymLink (src, dest);
   } else  {
      SX_THROW ("Can't copy the symbolic link '"
                + srcStr + "' to '" + destStr
                + "'. It is not possible to copy a symbolic link "
                  "that points to nowhere to itself.");
   }
}

void SxFSCopyAction::rmDestSymLinkAndCpEmptySymLink (const SxFileInfo &src,
                                                     const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   // Removing the destination link
   SxFSDeleteAction::removeSymLink (dest);
   // Performing the actual copy-operation on the link
   convEmptySymLinkToNothing (src, dest);
}

SxString SxFSCopyAction::getName () const
{
   SX_DBG_TRACE("");
   return "copy";
}
