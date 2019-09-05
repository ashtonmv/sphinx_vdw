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

#ifndef _SX_NAVIGATE_ACTION_H_
#define _SX_NAVIGATE_ACTION_H_

// Including header-files
#include <SxDir.h>
#include <SxFS.h>
#include <SxException.h>

/**
  \brief Provides functions for the navigation through the directory tree.

  \author Thomas Uchdorf, t.uchdorf@mpie.de
 */
class SX_EXPORT_FS SxFSNavigateAction
{
   public:
      // --- Constructors
      // Constructor
      SxFSNavigateAction ();
      // Destructor
      ~SxFSNavigateAction ();

      // --- Class-functions
      /**
        \brief Convenience function
        \sa SxDir::cd (const SxDir &)
       */
      static void cd (const SxString &);
      /**
        \brief Makes the directory specified via the passed parameter the
        current working directory.

        \par Example
        Set the current working directory to the user's home directory
        \code
        try  {
        SxDir::cd (SxDir::getHome ());
        } catch (SxException ex)  {
        ex.print ();
        }
        \endcode
       */
      static void cd (const SxDir &);

      /**
        \brief Convenience function
        \sa SxDir::pushd (const SxDir &)
       */
      static void pushd (const SxString &);
      /**
        \brief Changes to the specified new directory and adds the current
        working directory to the stack.

        This function corresponds with the UNIX command "pushd".
       */
      static void pushd (const SxDir &);

      /**
        \brief Changes to the previous directory and removes the entry from
        the directory stack

        This function corresponds with the UNIX command "popd".
       */
      static void popd ();


      // --- Methods

   protected:
      // --- Class-members
      /**
        \brief Used as directory-stack by popd and pushd
       */
      static SxList<SxDir> dirStack;

      //TODO Introduce the following member
      //     SxMap<SxFileInfoMgr, SxList<SxDir> > homeAndCwdAndDirStack;
      //     It should form a relation between file info managers e.g. file
      //     managers like konqueror or apple's finder and currently fix static
      //     members like the directory stack, the current working directory or
      //     the home directory

};

#endif /* _SX_NAVIGATE_ACTION_H_ */
