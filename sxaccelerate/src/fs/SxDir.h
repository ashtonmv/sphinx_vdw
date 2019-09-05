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
#ifndef _SX_DIR_H_
#define _SX_DIR_H_

// --- Including header-files
#include <SxError.h>
#include <SxFileInfo.h>
#include <SxString.h>
#include <SxFS.h>

/**
  \brief Directory representation in SPHInX

  \b SxDir = SPHInX Directory handler

  Use this class for performing elementary functions related to directory
  handling.

  \sa       SxFile
  \sa       SxSymLink
  \ingroup  group_os
  \author   Sixten Boeck, boeck@mpie.de
  \author   Thomas Uchdorf, t.uchdorf@mpie.de
 */
class SX_EXPORT_FS SxDir : public SxFileInfo
{
   public:

      // --- Constructors
      // Constructor
      /**
        \brief Creates a SxDir-object using the passed relative or absolute
        path. If no path is specified the current directory is used.
        Note, that this constructor will not open or create the folder.
        It just initializes the local member variables.
       */
      SxDir (const SxString &path_ = ".");

      /** Copies a SxFileInfo-object to a SxDir-object.
       */
      SxDir (const SxFileInfo &);

      // Destructor
      /**
        \brief Destroys a SxDir-object.
       */
      virtual ~SxDir ();

      // --- Methods
      /**
        \brief Checks whether a certain directory exists.

        This function evaluates SxFileInfo::exists () and checks additionally
        that the given directory is of the type 'directory'.
       */
      virtual bool exists () const;

      /// return path in which the executable is located
      static SxString getExecPath (); // #exception
      static SxString getInstallPath (); // #exception

};

#endif /* _SX_DIR_H_ */
