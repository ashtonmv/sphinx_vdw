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

#ifndef _SX_FS_ERROR_H_
#define _SX_FS_ERROR_H_

// Including header-files
#include <SxString.h>
#include <SxFS.h>
#include <SxException.h>

/**
  \brief Helper-class to create meaningful error messages.

  \author Thomas Uchdorf, t.uchdorf@mpie.de
 */
class SX_EXPORT_FS SxFSError
{
   public:
      // --- Constructors
      // Constructor
      SxFSError ();
      // Destructor
      ~SxFSError ();

      // --- Class-functions
      /** \brief Returns a detailed description of the error that occurred on
        a previous utime ()-call */
      static SxString getUtimeErrMsg ();

      /** \brief Returns a detailed description of the error that occurred on
        a previous chmod ()-call */
      static SxString getChmodErrMsg ();

      /** \brief Returns a detailed description of the error that occurred on
        a previous chown () or chgrp ()-call */
      static SxString getChownErrMsg ();

      /** \brief Returns a detailed description of the error that occurred on
        a previous symlink ()-call */
      static SxString getSymlinkErrMsg ();

      /** \brief Returns a detailed description of the error that occurred on
        a previous rename ()-call */
      static SxString getRenameErrMsg ();

      static SxString getOpenErrMsg ();

      static SxString getWriteErrMsg ();

      static SxString getChdirErrMsg ();

      static SxString getGetcwdErrMsg ();

      static SxString getOpendirErrMsg ();

      static SxString getRmdirErrMsg ();

      /**
        \brief Returns a descriptive error message when errors occur on
        calling readlink
       */
      static SxString getReadlinkErrMsg ();

      /**
        \brief Returns a descriptive error message when errors occur on
        calling unlink
       */
      static SxString getUnlinkErrMsg ();

#ifdef WIN32
      static void throwShortcutException (int,
                                          const SxString &,
                                          const SxString &lnkDest
                                             = SxString ());
#endif /* WIN32 */
};

#endif /* _SX_FS_ERROR_H_ */
