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
#include <SxFSMoveAction.h>
#include <SxFSDeleteAction.h>
#include <SxFSCopyAction.h>
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
SxFSMoveAction::SxFSMoveAction ()
{
   SX_DBG_TRACE("");
   // empty
}
// Destructor
SxFSMoveAction::~SxFSMoveAction ()
{
   SX_DBG_TRACE("");
   // empty
}

// --- Methods

void SxFSMoveAction::convFileToFile (const SxFileInfo &src,
                                     const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxFile(src).exists ());
   SX_CHECK (SxFile(dest).exists ());
   SxString const &srcStr = src.getAbsPath ();
   SxString const &destStr = dest.getAbsPath ();

   // Performing a self-test
   if (SxFileInfo::equals (src, dest))  {
      SX_THROW ("Can't move file '" + srcStr + "' to '" + destStr +
                "'. It is not possible to copy a file to itself.");
   } else  {
      SxFSMoveAction::rename (src, dest);
   }
}

void SxFSMoveAction::convFileToDir (const SxFileInfo &src,
                                    const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxFile(src).exists ());
   SX_CHECK (SxDir(dest).exists ());
   SxFSMoveAction::convFileToNothing (src, (dest / src.getName ()), true);
}

void SxFSMoveAction::convFileToNothing (const SxFileInfo &src,
                                        const SxFileInfo &dest,
                                        bool force)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxFile (src).exists ());
   // Ignoring the destination if force is set, i.e., removing destination if
   // necessary
   if (!force)  {
      SX_CHECK (!SxFile (dest).exists () && !SxDir(dest).exists () &&
                !SxSymLink (dest).exists ());
   }
   SxFSMoveAction::rename (src, dest);
}

void SxFSMoveAction::convFileToFileSymLink (const SxFileInfo &src,
                                            const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxFile(src).exists ());
   SX_CHECK (SxSymLink(dest).exists () && SxSymLink(dest).isFile ());

   // Calling a function that performs a special treatment
   SxFSMoveAction::specialMvFileToFileSymLink (src, dest);
}

void SxFSMoveAction::specialMvFileToFileSymLink (const SxFileInfo &src,
                                                 const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   //TODO Substitute the following part by a more sophisticated solution
   //     (implement this special rename function)
   // --- Using a special treatment since "cp -a" behaves different from "conv"
   //     concerning files as source and symbolic links that address files as
   //     destination
   //     Be aware that the wrong behavior only occurs on fall backs to copy
   //     operations instead of rename that means on renaming between different
   //     devices!
   SxFSDeleteAction::rmSymLink (dest);
   SxFSMoveAction::convFileToNothing (src, dest);
}

void SxFSMoveAction::convFileToDirSymLink (const SxFileInfo &src,
                                           const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxFile(src).exists ());
   SX_CHECK (SxSymLink(dest).exists () && SxSymLink(dest).isDir ());
   SxFSMoveAction::convFileToNothing (src, (dest / src.getName ()), true);
}

void SxFSMoveAction::convFileToEmptySymLink (const SxFileInfo &src,
                                             const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxFile(src).exists ());
   SX_CHECK (SxSymLink(dest).exists () && !SxSymLink(dest).isValid ());
   // Calling a function that performs a special treatment
   SxFSMoveAction::specialMvFileToEmptySymLink (src, dest);
}

void SxFSMoveAction::specialMvFileToEmptySymLink (const SxFileInfo &src,
                                                  const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   //TODO Substitute the following part by a more sophisticated solution
   //     (implement this special rename function)
   // --- Using a special treatment since "cp -a" behaves different from "conv"
   //     concerning files as source and symbolic links that address nothing as
   //     destination
   //     Be aware that the wrong behavior only occurs on fall backs to copy
   //     operations instead of rename that means on renaming between different
   //     devices!
   SxFSDeleteAction::rmSymLink (dest);
   SxFSMoveAction::convFileToNothing (src, dest);
}

void SxFSMoveAction::convDirToFile (const SxFileInfo &src,
                                    const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxDir(src).exists ());
   SX_CHECK (SxFile(dest).exists ());
   SxFSMoveAction::rename (src, dest);
}

void SxFSMoveAction::convDirToDir (const SxFileInfo &src,
                                   const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxDir(src).exists ());
   SX_CHECK (SxDir(dest).exists ());
   SxString const &srcStr = src.getAbsPath ();
   SxString const &destStr = dest.getAbsPath ();
   // Performing a self-test
   if (SxFileInfo::equals (src, dest))  {
      SX_THROW ("Can't move directory '" + srcStr + "' to '" + destStr +
                "'. It is not possible to move a directory to itself.");
   } else  {
      SxFSMoveAction::convDirToNothing (src, (dest / src.getName ()), true);
   }
}

void SxFSMoveAction::convDirToNothing (const SxFileInfo &src,
                                       const SxFileInfo &dest,
                                       bool force)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxDir(src).exists ());
   // Ignoring the destination if force is set, i.e., removing destination if
   // necessary
   if (!force)  {
      SX_CHECK (!SxFile(dest).exists () && !SxDir(dest).exists () &&
                !SxSymLink(dest).exists ());
   }
   SxFSMoveAction::rename (src, dest);
}

void SxFSMoveAction::convDirToFileSymLink (const SxFileInfo &src,
                                           const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxDir(src).exists ());
   SX_CHECK (SxSymLink(dest).exists () && SxSymLink(dest).isFile ());
   SxFSMoveAction::rename (src, dest);
}

void SxFSMoveAction::convDirToDirSymLink (const SxFileInfo &src,
                                          const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxDir(src).exists ());
   SX_CHECK (SxSymLink(dest).exists () && SxSymLink(dest).isDir ());
   SxFSMoveAction::convDirToNothing (src, (dest / src.getName ()), true);
}

void SxFSMoveAction::convDirToEmptySymLink (const SxFileInfo &src,
                                            const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxDir(src).exists ());
   SX_CHECK (SxSymLink(dest).exists () && !SxSymLink(dest).isValid ());
   // Deleting the destination
   SxFSDeleteAction::rmSymLink (dest);
   SxFSMoveAction::rename (src, dest);
}

void SxFSMoveAction::convFileSymLinkToFile (const SxFileInfo &src,
                                            const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && SxSymLink (src).isFile ());
   SX_CHECK (SxFile (dest).exists ());
   SxFSMoveAction::rename (src, dest);
}

void SxFSMoveAction::convFileSymLinkToDir (const SxFileInfo &src,
                                           const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && SxSymLink (src).isFile ());
   SX_CHECK (SxDir (dest).exists ());
   SxFSMoveAction::convFileSymLinkToNothing (src, (dest/src.getName ()), true);
}

void SxFSMoveAction::convFileSymLinkToNothing (const SxFileInfo &src,
                                               const SxFileInfo &dest,
                                               bool force)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && SxSymLink (src).isFile ());
   // Ignoring the destination if force is set, i.e., removing destination if
   // necessary
   if (!force)  {
      SX_CHECK (!SxFile (dest).exists () && !SxDir (dest).exists () &&
                !SxSymLink (dest).exists ());
   }
   SxFSMoveAction::rename (src, dest);
}

void SxFSMoveAction::convFileSymLinkToFileSymLink (const SxFileInfo &src,
                                                   const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && SxSymLink (src).isFile ());
   SX_CHECK (SxSymLink (dest).exists () && SxSymLink (dest).isFile ());
   if (SxFileInfo::equals (src, dest)) {
      SX_THROW ("Can't move '" + src.getAbsPath ()
               + "' to '" + dest.getAbsPath ()
               + "'. The symbolic links correspond.");
   } else  {
      SxFSMoveAction::rename (src, dest);
   }
}

void SxFSMoveAction::convFileSymLinkToDirSymLink (const SxFileInfo &src,
                                                  const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && SxSymLink (src).isFile ());
   SX_CHECK (SxSymLink (dest).exists () && SxSymLink (dest).isDir ());
   SxFSMoveAction::convFileSymLinkToNothing (src, (dest/src.getName ()), true);
}

void SxFSMoveAction::convFileSymLinkToEmptySymLink (const SxFileInfo &src,
                                                    const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && SxSymLink (src).isFile ());
   SX_CHECK (SxSymLink (dest).exists () && !SxSymLink (dest).isValid ());
   // Deleting the destination
   SxFSDeleteAction::rmSymLink (dest);
   SxFSMoveAction::rename (src, dest);
}

void SxFSMoveAction::convDirSymLinkToFile (const SxFileInfo &src,
                                           const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && SxSymLink (src).isDir ());
   SX_CHECK (SxFile (dest).exists ());
   SxFSMoveAction::rename (src, dest);
}

void SxFSMoveAction::convDirSymLinkToDir (const SxFileInfo &src,
                                          const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && SxSymLink (src).isDir ());
   SX_CHECK (SxDir (dest).exists ());
   SxFSMoveAction::convDirSymLinkToNothing (src, (dest/src.getName ()), true);
}

void SxFSMoveAction::convDirSymLinkToNothing (const SxFileInfo &src,
                                              const SxFileInfo &dest,
                                              bool force)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && SxSymLink (src).isDir ());
   // Ignoring the destination if force is set, i.e., removing destination if
   // necessary
   if (!force)  {
      SX_CHECK (!SxFile (dest).exists () && !SxDir (dest).exists () &&
                !SxSymLink (dest).exists ());
   }
   SxFSMoveAction::rename (src, dest);
}

void SxFSMoveAction::convDirSymLinkToFileSymLink (const SxFileInfo &src,
                                                  const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && SxSymLink (src).isDir ());
   SX_CHECK (SxSymLink (dest).exists () && SxSymLink (dest).isFile ());
   SxFSMoveAction::rename (src, dest);
}

void SxFSMoveAction::convDirSymLinkToDirSymLink (const SxFileInfo &src,
                                                 const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && SxSymLink (src).isDir ());
   SX_CHECK (SxSymLink (dest).exists () && SxSymLink (dest).isDir ());
   SxFSMoveAction::convDirSymLinkToNothing (src, (dest/src.getName ()), true);
}

void SxFSMoveAction::convDirSymLinkToEmptySymLink (const SxFileInfo &src,
                                                   const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && SxSymLink (src).isDir ());
   SX_CHECK (SxSymLink (dest).exists () && !SxSymLink (dest).isValid ());
   // Deleting the destination
   SxFSDeleteAction::rmSymLink (dest);
   SxFSMoveAction::rename (src, dest);
}

void SxFSMoveAction::convEmptySymLinkToFile (const SxFileInfo &src,
                                             const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && !SxSymLink (src).isValid ());
   SX_CHECK (SxFile (dest).exists ());
   SxFSMoveAction::rename (src, dest);
}

void SxFSMoveAction::convEmptySymLinkToDir (const SxFileInfo &src,
                                            const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && !SxSymLink (src).isValid ());
   SX_CHECK (SxDir (dest).exists ());
   SxFSMoveAction::convEmptySymLinkToNothing (src, (dest/src.getName ()), true);
}

void SxFSMoveAction::convEmptySymLinkToNothing (const SxFileInfo &src,
                                                const SxFileInfo &dest,
                                                bool force)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && !SxSymLink (src).isValid ());
   // Ignoring the destination if force is set, i.e., removing destination if
   // necessary
   if (!force)  {
      SX_CHECK (!SxFile (dest).exists () && !SxDir (dest).exists () &&
                !SxSymLink (dest).exists ());
   }
   SxFSMoveAction::rename (src, dest);
}

void SxFSMoveAction::convEmptySymLinkToFileSymLink (const SxFileInfo &src,
                                                    const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && !SxSymLink (src).isValid ());
   SX_CHECK (SxSymLink (dest).exists () && SxSymLink (dest).isFile ());
   SxFSMoveAction::rename (src, dest);
}

void SxFSMoveAction::convEmptySymLinkToDirSymLink (const SxFileInfo &src,
                                                   const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && !SxSymLink (src).isValid ());
   SX_CHECK (SxSymLink (dest).exists () && SxSymLink (dest).isDir ());
   SxFSMoveAction::convEmptySymLinkToNothing (src, (dest/src.getName ()), true);
}

void SxFSMoveAction::convEmptySymLinkToEmptySymLink (const SxFileInfo &src,
                                                     const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SX_CHECK (SxSymLink (src).exists () && !SxSymLink (src).isValid ());
   SX_CHECK (SxSymLink (dest).exists () && !SxSymLink (dest).isValid ());
   if (SxFileInfo::equals (src, dest)) {
      SX_THROW ("Can't move '" + src.getAbsPath ()
                + "' to '" + dest.getAbsPath ()
                + "'. The symbolic links correspond.");
   } else  {
      // Deleting the destination
      SxFSDeleteAction::rmSymLink (dest);
      SxFSMoveAction::rename (src, dest);
   }
}

void SxFSMoveAction::rename (const SxFileInfo &src,
                             const SxFileInfo &dest)
{
   SX_DBG_TRACE("");
   SxString const &srcStr = src.getAbsPath ();
   SxString const & destStr = dest.getAbsPath ();
   int error = ::rename (srcStr.getElems (),
                         destStr.getElems ());
   if (error)  {
      try {
         SxFSMoveAction::cpAndRm (src, dest);
      } catch (SxException e)  {
         SX_THROW ("Can't rename '" + srcStr + "' to '" + destStr
                   + "': " + SxString(e.getMessage ()));
      }
   }
}

void SxFSMoveAction::cpAndRm (const SxFileInfo &src,
                              const SxFileInfo &dest)
{
   SX_DBG_TRACE("");

   bool existedBefore (false);
   // Checking whether there is already a valid link, file or directory at
   // the destination
   if (dest.exists ())  {
      existedBefore = true;
   }
   // Trying to copy source to destination
   try  {
      SxFSCopyAction cp;
      cp.execute (src, dest);
   } catch (SxException ex)  {
      // Cleaning up rudimentarily if possible
      if (!existedBefore)  {
         SxFSDeleteAction::rm_r (dest);
      }
      throw ex;
   }
   // Removing the source if the copying succeeded
   SxFSDeleteAction::rm_r (src);
}

SxString SxFSMoveAction::getName () const
{
   SX_DBG_TRACE("");
   return "move";
}
