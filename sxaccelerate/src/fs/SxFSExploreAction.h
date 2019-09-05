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

#ifndef _SX_FS_EXPLORE_ACTION_H_
#define _SX_FS_EXPLORE_ACTION_H_

// Including header-files
#include <SxFileInfo.h>
#include <SxFile.h>
#include <SxDir.h>
#include <SxSymLink.h>
#include <SxFS.h>

/**
  \brief Provides all sorts of operations related to exploring/searching files,
  directories and symbolic links.

  \author Thomas Uchdorf, t.uchdorf@mpie.de
 */
class SX_EXPORT_FS SxFSExploreAction
{
   public:
      // --- Constructors
      // Constructor
      SxFSExploreAction ();
      // Destructor
      ~SxFSExploreAction ();

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

      static SxList<SxFISortedByTime> ls_t (const SxFileInfo &);

      static SxList<SxFISortedBySize> ls_S (const SxFileInfo &);

      static SxList<SxFileInfo> find (const SxString &target,
                                      const SxDir &curDir = SxDir ("."));

      static SxList<SxFileInfo> ls (const SxFileInfo &);

      static void listDir (const SxString path, SxList<SxString> *fileList);

      static SxList<SxFileInfo> getFileInfos (const SxFileInfo &);

      static SxList<SxFile> getFiles (const SxFileInfo &);

      static SxList<SxDir>  getDirs (const SxFileInfo &);

      static SxList<SxSymLink>  getSymLinks (const SxFileInfo &);

   protected:

      // --- Class-functions
      static void prepareRegExpPattern (SxString *);

      static void finishRegExpPattern (SxString *);

      static SxList<SxFileInfo> searchNode (const SxFileInfo &curDir,
                                            const SxString   &pattern,
                                            SxList<SxDir>    *foundDirs,
                                            SxList<SxDir>    *subDirs = NULL);

      static SxList<SxFileInfo> processNode (const SxString   &name,
                                             const SxFileInfo &curDir,
                                             const SxString   &pattern,
                                             SxList<SxDir>    *foundDirs,
                                             SxList<SxDir>    *subDirs);


      static bool contains (const SxList<SxFileInfo> &visited,
                            const SxFileInfo &current);

      static void findName (const SxDir &, const SxString &,
                            SxList<SxFileInfo> *, SxList<SxFileInfo> *,
                            bool = false);
};

#endif /* _SX_FS_EXPLORE_ACTION_H_ */
