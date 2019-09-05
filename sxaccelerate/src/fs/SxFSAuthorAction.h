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

#ifndef _SX_FS_AUTHOR_ACTION_H_
#define _SX_FS_AUTHOR_ACTION_H_

#include <SxFS.h>
#include <SxFileInfo.h>

/**
  \brief This class contains class-functions related to authorization.

  \author Thomas Uchdorf, t.uchdorf@mpie.de
 */
class SX_EXPORT_FS SxFSAuthorAction
{
   public:
      SxFSAuthorAction ();
      ~SxFSAuthorAction ();

      static void chmod (int perm, const SxFileInfo& target);

      static void chmod (const SxString &modStr, const SxFileInfo& target);

      static void chmod (const SxFileInfo &,
                         SxFileInfo::AccessGroup,
                         SxFileInfo::AccessType);

      // This method is intended for internal usage of the FSAction-classes
      // only
      static void chmod (const SxFileInfo &, int);
      static void chgrp (const SxFileInfo &, sxgid_t);
      static void chown (const SxFileInfo &, sxuid_t);
};

#endif /* _SX_FS_AUTHOR_ACTION_H_ */
