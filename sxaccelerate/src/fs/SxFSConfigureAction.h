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

#ifndef _SX_FS_CONFIGURE_ACTION_H_
#define _SX_FS_CONFIGURE_ACTION_H_

// Including header-files
#include <SxFS.h>
#include <SxFileInfo.h>
#include <SxException.h>

/**
  \brief Contains class-functions for synchronization.

  \author Thomas Uchdorf, t.uchdorf@mpie.de
 */
class SX_EXPORT_FS SxFSConfigureAction
{
   public:
      SxFSConfigureAction ();
      ~SxFSConfigureAction ();

      static void syncUIDAndGIDAndTimes (const SxFileInfo &src,
                                         const SxFileInfo &dest);

   protected:

      static void syncTimes (const SxFileInfo &src,
                             const SxFileInfo &dest);
};

#endif /* _SX_FS_CONFIGURE_ACTION_H_ */
