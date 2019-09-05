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

#ifndef _SX_FS_DELETE_ACTION_H_
#define _SX_FS_DELETE_ACTION_H_

// Including header-files
#include <SxFS.h>
#include <SxFileInfo.h>
#include <SxFile.h>
#include <SxDir.h>
#include <SxSymLink.h>

/**
  \brief Enables deletion of files, directories and symbolic links.

  \author Thomas Uchdorf, t.uchdorf@mpie.de
 */
class SX_EXPORT_FS SxFSDeleteAction
{
   public:
      // --- Constructors
      // Constructor
      SxFSDeleteAction ();
      // Destructor
      ~SxFSDeleteAction ();

      // --- Methods
      static void rm (const SxFile &target);

      /**
        \brief Removes an empty directory.

        This function corresponds with the UNIX command "rmdir".
       */
      static void rmdir (const SxDir &);

      /**
        \brief Removes the passed symbolic link
       */
      static void rmSymLink (const SxSymLink &);

      static void rm_r (const SxFileInfo &);

      //static void removeDir (const SxFileInfo &);

      static void removeFile (const SxFileInfo &);

      static void removeSymLink (const SxFileInfo &);
};

#endif /* _SX_FS_DELETE_ACTION_H_ */
