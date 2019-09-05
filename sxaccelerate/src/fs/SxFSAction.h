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

#ifndef _SX_FS_ACTION_H_
#define _SX_FS_ACTION_H_

// Including header-files
#include <SxString.h>
#include <SxArray.h>
#include <SxList.h>
#include <SxFileInfo.h>
#include <SxFile.h>
#include <SxDir.h>
#include <SxSymLink.h>
#include <SxFS.h>

/**
  \brief Provides static functions to operate on files, directories and
  symbolic links.

  SxFSAction #exception

  \ingroup  group_os
  \author Thomas Uchdorf, t.uchdorf@mpie.de
 */
class SX_EXPORT_FS SxFSAction
{
   public:
      // --- Constructors
      // Constructor
      /**
        \brief Creates an SxFSAction-object.
       */
      SxFSAction ();
      // Destructor
      /**
        \brief Destroys an SxFSAction-object.
       */
      ~SxFSAction ();

      // --- Class-functions
      // --- Authorize-functions
      static void chmod (const SxString &path, int mode);
      static void chmod (int, const SxString &);
      static void chmod (int, const SxFileInfo &);

      static void chmod (const SxString &, const SxString &);
      static void chmod (const SxString &, const SxFileInfo &);

      /**
        \brief Convenience function
        \sa SxFSAction::chmod (const SxFileInfo &,
        SxFileInfo::AccessGroup,
        SxFileInfo::AccessType)
       */
      static void chmod (const SxString &,
                         SxFileInfo::AccessGroup,
                         SxFileInfo::AccessType);
      /**
        \brief Changes the permissions of a file system element.

        This function can be used to change the file permissions.
        The values of the last two arguments are bitwisely OR'ed.

        \par Example:
        Make a file readable, writable, and executable for the user
        and the group.
        \code
        SxFile file ("myFile.dat");
        SxFSAction::chmod (file, SxFile::Owner | SxFile::Group,
        SxFile::Readable | SxFile::Writable | SxFile::Executable);
        \endcode
       */
      static void chmod (const SxFileInfo &,
                         SxFileInfo::AccessGroup,
                         SxFileInfo::AccessType);

      /**
        \brief Convenience function
        \sa SxFSAction::chgrp (const SxFileInfo &, unsigned int)
       */
      static void chgrp (unsigned int, const SxString &);

      /**
        \brief Changes the group of a file system element.
       */
      static void chgrp (unsigned int, const SxFileInfo &);

      /**
        \brief Convenience function
        \sa SxFSAction::chown (const SxFileInfo &, unsigned int)
       */
      static void chown (unsigned int, const SxString &);

      /**
        \brief Changes the owner of a file system element.
       */
      static void chown (unsigned int, const SxFileInfo &);


      // --- Convert-functions
      /**
        \brief Convenience function
        \sa SxFSAction::mv (const SxFileInfo &, const SxFileInfo &)
       */
      static void mv (const SxString &, const SxString &);
      /**
        \brief Moves and renames file system elements within the file system.
       */
      static void mv (const SxFileInfo &, const SxFileInfo &);

      /**
        \brief Convenience function
        \sa SxFSAction::cp (const SxFileInfo &, const SxFileInfo &)
       */
      static void cp (const SxString &, const SxString &);
      /** \brief Copies files, directories or symbolic links.

        This function is similar to the "cp -a"-command known from UNIX-based systems.
       */
      static void cp (const SxFileInfo &, const SxFileInfo &);

      // --- Create-functions
      /**
        \brief Convenience function
        \sa SxFSAction::touch (const SxFileInfo &, int mode)
       */
      static void touch (const SxString &, int mode = 0644);
      /**
        \brief Creates a file or updates its time stamps.

        This function corresponds with the UNIX-command "touch".
        \param file that is to be created or updated concerning its time stamp
        \param permissions that should be applied on file creation
       */
      static void touch (const SxFileInfo &, int mode = 0644);

      /**
        \brief Convenience function
        UNIX: ln -sf /etc/profile .myprofile
        SX:   ln_sf ("/etc/profile", ".myprofile");
        \sa SxFSAction::ln_sf (const SxFileInfo &, const SxFileInfo &)
       */
      static SxFileInfo ln_sf (const SxString &, const SxString &);

      /** \brief Creates a symbolic link to a file system element.

        This function creates a symbolic link and returns it.
        It is similar to the UNIX-command "ln -sf <original> <link>" */
      static SxFileInfo ln_sf (const SxFileInfo &, const SxFileInfo &);

      /**
        \brief Creates a unique file on the file system.

        For creation of temporary files you have to ensure to use a
        filename, that is unique. This function returns such a
        filename.
        \param tmpDir location of temporary directory, e.g. "/temp"
        if not provided the content of the
        environment variable TMPDIR is taken.
        \param buffer content in bytes that should be used to initialize
        the new temporary file
        \setenv TMPDIR
       */
      static SxFile createTmpFile (const SxString &tmpDir="",
                                   const SxArray<unsigned char> &buffer =
                                   SxArray<unsigned char> ());

      /**
        \brief Convenience function
        \sa SxFSAction::mkdir (const SxDir &, int)
       */
      static void mkdir (const SxString &, int umask = 0777);
      /**
        \brief Creates a new directory.

        This command is equivalent to the UNIX-command "mkdir". It creates
        a new empty subfolder inside the specified folder. All directories
        specified in the passed path have to exist. If they are to be created
        as well, SxFSAction::mkdir_p (const SxDir &, int) should be used instead.
        \par Example:
        Create a subdirectory abc/ under $HOME/def.
        \code
      // Fails if $HOME/def does not exists
      try  {
      SxFSAction::mkdir ((SxFSAction::getHome ()/"def")/"abc");
      } catch (SxException ex)  {
      ex.print ();
      }
      \endcode
       */
      static void mkdir (const SxDir &, int umask = 0777);


      /**
        \brief Convenience function
        \sa SxFSAction::mkdir_p (const SxDir &, int)
       */
      static void mkdir_p (const SxString &, int umask = 0777);

      /**
        \brief Creates a new directory and also  non-existing parent
        directories.

        This command is equivalent to the UNIX-command "mkdir -p". It creates
        all directories that do not exist and contribute to the path of the
        target directory.
        \par Example:
        Create a subdirectory abc/ under $HOME/def.
        \code
      // Works provided proper permission no matter if $HOME/def exists or
      // not
      try  {
      SxFSAction::mkdir_p ((SxFSAction::getHome ()/"def")/"abc");
      } catch (SxException ex)  {
      ex.print ();
      }
      \endcode
       */
      static void mkdir_p (const SxDir &, int umask = 0777);

      // --- Delete-functions
      /**
        \brief Convenience function
        \sa SxFSAction::rm (const SxFile &)
       */
      static void rm (const SxString &);
      /**
        \brief Removes a file.

        This function corresponds with the UNIX command "rm".
        \param file that is to be removed
       */
      static void rm (const SxFile &);

      /**
        \brief Convenience function
        \sa SxFSAction::rmdir (const SxDir &)
       */
      static void rmdir (const SxString &);
      /**
        \brief Removes an empty directory.

        This function corresponds with the UNIX command "rmdir".
       */
      static void rmdir (const SxDir &);

      /**
        \brief Convenience function
        \sa SxFSAction::rm_r (const SxFileInfo &)
       */
      static void rm_r (const SxString &);
      /**
        \brief Removes directories, files or symbolic links (recursively in
        case of directories)
       */
      static void rm_r (const SxFileInfo &);

      // --- Navigate-functions
      /**
        \brief Convenience function
        \sa SxFSAction::cd (const SxDir &)
       */
      static void cd (const SxString &);
      /**
        \brief Makes the directory specified via the passed parameter the
        current working directory.

        \par Example
        Set the current working directory to the user's home directory
        \code
        try  {
        SxFSAction::cd (SxFSAction::getHome ());
        } catch (SxException ex)  {
        ex.print ();
        }
        \endcode
       */
      static void cd (const SxDir &);

      /**
        \brief Convenience function
        \sa SxFSAction::pushd (const SxDir &)
       */
      static void pushd (const SxString &);
      /**
        \brief Changes to the specified new directory and adds the current
        working directory to the stack.

        This function corresponds with the UNIX-command "pushd".
       */
      static void pushd (const SxDir &);

      /**
        \brief Changes to the previous directory and removes the entry from
        the directory stack

        This function corresponds with the UNIX-command "popd".
       */
      static void popd ();

      // --- Explore-functions
      /**
        \brief Returns the user's home directory
       */
      static SxDir getHome ();

      /**
        \brief Returns the directory for temporary files
       */
      static SxDir getTmp ();

      /**
        \brief Returns the current working directory.
       */
      static SxDir pwd ();

      /**
        \brief Recursively traverses all directories starting at the current
        one and searches for a file name that corresponds with the passed
        string. Here, wildcards like '*' are also supported.
       */
      static SxList<SxFileInfo> find (const SxString &, bool chdir=false);

      /**
        \brief Convenience function
        \sa SxFSAction::ls_t (const SxFileInfo &)
       */
      static SxList<SxFISortedByTime> ls_t (const SxString &path);

      /**
        \brief Returns all file system elements corresponding with the
        specified path sorted by time.
        \par Example
        The usage is analogous with SxFSAction::ls ()
       */
      static SxList<SxFISortedByTime> ls_t (const SxFileInfo &path);

      /**
        \brief Convenience function
        \sa SxFSAction::ls_S (const SxFileInfo &)
       */
      static SxList<SxFISortedBySize> ls_S (const SxString &path);

      /**
        \brief Returns all file system elements corresponding with the
        specified path sorted by size.
        \par Example
        The usage is analogous with SxFSAction::ls ()
       */
      static SxList<SxFISortedBySize> ls_S (const SxFileInfo &path);

      /**
        \brief Convenience function
        \sa SxFSAction::ls (const SxFileInfo &)
       */
      static SxList<SxFileInfo> ls (const SxString &path);
      /**
        \brief Returns a list of file system elements like files, directories
        or symbolic links that is contained in the determined directory and
        maybe matches a defined pattern.
        \par Example
        Get all *.dat or *.txt files, directories or symbolic links in the
        current working directory
        \code
        SxFileInfo dir (SxDir (".") / "*.{dat,txt}");
        SxList<SxFileInfo> elements;
        try  {
        elements = SxFSAction::ls (dir);
        } catch (SxException ex)  {
        ex.print ();
        }
        cout << "elements:" << endl;
        SxList<SxFileInfo>::Iterator itElem;
        for (itElem = elements.begin (); itElem != elements.end (); ++itElem)  {
        if (itElem == elements.begin ())  {
        cout << itElem->getAbsPath ();
        } else  {
        cout << endl << itElem->getAbsPath ();
        }
        }
        \endcode
       */
      static SxList<SxFileInfo> ls (const SxFileInfo &path);

      /** Search for executable in binary search path */
      static SxString which (const SxString &file, const SxString &path="");

      /** Search for executable in binary search path */
      static SxList<SxString> where (const SxString &file,
                                     const SxString &path="");

      /** Collect all filenames in given directory. Recursive. */
      static SxList<SxString> listDir (const SxString &path);

      /**
        \brief Convenience function
        \sa SxFSAction::getFiles (const SxFileInfo &)
       */
      static SxList<SxFile> getFiles (const SxString &path);
      /**
        \brief Returns all files in a directory matching a name pattern.

        This function returns all files in a directory matching a
        given pattern, such as "*abc" or "*def".
        \par Example
        Get all files in subdirectories of the current working directory that
        match "*abc" or "*def".
        \code
        SxFileInfo dir (SxDir (".") / "*{abc,def}");
        SxList<SxFile> files;
        try  {
        files = SxFSAction::getFiles (dir);
        } catch (SxException ex)  {
        ex.print ();
        }
        cout << "files:" << endl;
        SxList<SxFile>::Iterator itFiles;
        for (itFiles = files.begin (); itFiles != files.end (); ++itFiles)  {
        if (itFiles == files.begin ())  {
        cout << itFiles->getAbsPath ();
        } else  {
        cout << endl << itFiles->getAbsPath ();
        }
        }
        \endcode
       */
      static SxList<SxFile> getFiles (const SxFileInfo &path);
      /**
        \brief Convenience function
        \sa SxFSAction::getDirs (const SxFileInfo &)
       */
      static SxList<SxDir>  getDirs (const SxString &path);
      /**
        \brief Returns all directories which are listed in the directory
        and, if defined, correspond with a specified pattern.

        Every directory may contain files and/or directories and/or symbolic
        links as child nodes. This function returns only the subdirectories
        which are listed in this directory. The files can be retrieved with
        SxFSAction::getFiles () and the symbolic links with
        SxFSAction::getSymLinks ().
        \par Example
        The usage is analogous with SxFSAction::getFiles ()
       */
      static SxList<SxDir> getDirs (const SxFileInfo &path);

      /**
        \brief Convenience function
        \sa SxFSAction::getSymLinks (const SxFileInfo &)
       */
      static SxList<SxSymLink>  getSymLinks (const SxString &path);
      /**
        \brief Returns all symbolic links which are listed in the directory
        and, if defined, correspond with a specified pattern.
        \par Example
        The usage is analogous with SxFSAction::getFiles ()
       */
      static SxList<SxSymLink>  getSymLinks (const SxFileInfo &path);

      // --- Meta-functions (for convenience)
      /**
        \brief Checks whether the file system element exists.
       */
      static bool exists (const SxString &str);

      /**
        \brief Tests if the file system element exists and is a file.
       */
      static bool test_f (const SxString &str);

      /**
        \brief Checks whether the file system element exists and is a directory.
       */
      static bool test_d (const SxString &str);

      /**
        \brief Tests if the file system element exists and is a symbolic link.
       */
      static bool test_L (const SxString &str);
};

#endif /* _SX_FS_ACTION_H_ */
