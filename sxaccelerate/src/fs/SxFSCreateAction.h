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

#ifndef _SX_FS_CREATE_ACTION_H_
#define _SX_FS_CREATE_ACTION_H_

// Including header-files
#include <SxFS.h>
#include <SxFileInfo.h>
#include <SxFile.h>
#include <SxDir.h>
#include <SxSymLink.h>

/**
  \brief Enables timestamp updates as well as creation of files, directories and
  symbolic links.

  \author Thomas Uchdorf, t.uchdorf@mpie.de
 */
class SX_EXPORT_FS SxFSCreateAction
{
   public:
      // --- Constructors
      // Constructor
      SxFSCreateAction ();
      // Destructor
      virtual ~SxFSCreateAction ();

      // --- Methods
      //virtual void create (const SxFileInfo &, int = 0777);
      virtual void update (const SxFileInfo &, int = 0777) = 0;

      static void touch (const SxFileInfo &, int mode = 0644);

      /**
        \brief Creates a new directory.

        This command is equivalent to the UNIX command "mkdir". It creates
        a new empty subfolder inside the specified folder. All directories
        specified in the passed path have to exist. If they are to be created
        as well, SxDir::mkdir_p (const SxDir &, int) should be used instead.
        \par Example:
        Create a subdirectory abc/ under $HOME/def.
        \code
      // Fails if $HOME/def does not exists
      try  {
      SxDir::mkdir ((getHome ()/"def")/"abc");
      } catch (SxException ex)  {
      ex.print ();
      }
      \endcode
       */
      static void mkdir (const SxDir &, int umask = 0777);

      /**
        \brief Creates a new directory and also  non-existing parent
        directories.

        This command is equivalent to the UNIX command "mkdir -p". It creates
        all directories that do not exist and contribute to the path of the
        target directory.
        \par Example:
        Create a subdirectory abc/ under $HOME/def.
        \code
      // Works provided proper permission no matter if $HOME/def exists or
      // not
      try  {
      SxDir::mkdir_p ((getHome ()/"def")/"abc");
      } catch (SxException ex)  {
      ex.print ();
      }
      \endcode
       */
      static void mkdir_p (const SxDir &, int umask = 0777);

      static SxFileInfo ln_sf (const SxFileInfo &,
                               const SxFileInfo &);

      /**
        \brief Create a unique file on the file system.

        For creation of temporary files you have to ensure to use a
        filename, that is unique. This function returns such a
        filename.
        \param tmpDir  location of temporary directory, e.g. "/temp"
        if not provided the content of the
        environment variable TMPDIR is taken.
        \setenv TMPDIR
       */
      static SxFile createTmpFile (const SxString &tmpDir="",
                                   const SxArray<unsigned char> &buffer =
                                   SxArray<unsigned char> ());

   protected:
      static void updateFile (const SxFileInfo &, int = 0777);

      static void updateDir (const SxFileInfo &, int = 0777);

#ifdef WIN32
      static void createShortcut (const SxString &lnkSrc,
                                  const SxString &lnkDest);
#endif /* WIN32 */

};

#endif /* _SX_FS_CREATE_ACTION_H_ */
