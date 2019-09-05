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
#ifndef _SX_FILE_H_
#define _SX_FILE_H_

// --- Including header-files
#include <SxFileInfo.h>
#include <SxString.h>
#include <SxList.h>
#include <SxFS.h>

/**
  \brief File representation in SPHInX

  \b SxFile = SPHInX File handler

  Use this class for performing elementary functions related to files,
  such as checking for existence, creating temporary files, copying
  files, creating symbolic links and so on.

  \sa       SxDir
  \sa       SxSymLink
  \ingroup  group_os
  \author   Sixten Boeck, boeck@mpie.de
  \author   Thomas Uchdorf, t.uchdorf@mpie.de
 */
class SX_EXPORT_FS SxFile : public SxFileInfo
{
   public:

      // --- Constructors
      // Constructor
      /**
        \brief Initializes a file with an absolute or relative path.

        Note, that this function will not open or create a file. It
        just initializes the local member variables.
       */
      SxFile (const SxString & path_ = SxString ());

      /**
        \brief Initialize file from its base class object.
       */
      SxFile (const SxFileInfo &);

      // Destructor
      /**
        \brief Destroys the SxFile-object.
       */
      virtual ~SxFile ();

      // --- Methods
      /**
        \brief Get that part of filename that would be sufficient to find
        the file in the given paths. Deprecated!

        \param path ':'-separated list of paths
       */
      SxString getIncludeName (const SxString &path_);

      /**
        \brief Get that part of filename that would be sufficient to find
        the file in the given paths. Deprecated!

        \param list of paths
       */
      SxString getIncludeName (const SxList<SxString> &path_);

      /**
        \brief Tests if a file exists.
       */
      virtual bool exists () const;
};

#endif /* _SX_FILE_H_ */
