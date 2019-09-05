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

// Defining the header-guard
#ifndef _SX_SYM_LINK_H_
#define _SX_SYM_LINK_H_

// --- Including header-files
#include <SxString.h>
#include <SxFileInfo.h>
#include <SxException.h>
#include <SxFS.h>

/**
  \brief Symbolic link representation in SPHInX

  \b SxSymLink = SPHInX Symbolic link handler

  Use this class for performing elementary functions related to symbolic links,
  such as checking for existence, validity and so on. Take into account that
  files and directories can also be symbolic links and that only a call to
  SxSymLink::exists () determines if the element is a symbolic link.
  SxDir::getFiles () for example also returns valid symbolic links to files
  instead of just pure files.

  \sa       SxDir
  \sa       SxFile
  \ingroup  group_os
  \author   Thomas Uchdorf, t.uchdorf@mpie.de
 */
class SX_EXPORT_FS SxSymLink : public SxFileInfo
{
   public:

      // --- Constructors
      // Constructor
      /**
        \brief Creates a SxSymLink-object
       */
      SxSymLink (const SxString &path_ = SxString ());

      /**
        \brief Turns a SxFileInfo-object into a SxSymLink-object
       */
      SxSymLink (const SxFileInfo &);

      // Destructor
      /**
        \brief Destroys a SxSymLink-object
       */
      virtual ~SxSymLink ();

      // --- Methods
      /**
        \briefs Returns the target, i.e., the file system element the
        symbolic link references, or an empty string if the symbolic link is a
        dead link
       */
      SxString const getTarget () const;

      /**
        \brief Tests whether the file system element exists as well as it is
        a symbolic link. Note, that this function does not check the validity
        of the symbolic link. Thus, dead links also exist. To check whether the
        link is a dead link use SxSymLink::isValid ().
       */
      virtual bool exists () const;

      /**
        \brief Checks if the symbolic link points to a target that is an
        existing file system element.
       */
      bool isValid () const;

#ifdef WIN32
      static SxString resolveShortcutsInPath (const SxString &lnkSrc);

      static SxString resolveShortcutsInAbsPath (const SxString &lnkSrc);
#endif /* WIN32 */

   protected:
#ifdef WIN32
      static SxString resolveShortcuts (const SxString &lnkSrc,
                                        bool keepName = false);

      static SxString getShortcutTarget (const SxString &lnkSrc);
#endif /* WIN32 */
};

#endif /* _SX_SYM_LINK_H_ */

