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

#ifndef _SX_FS_CONVERT_ACTION_H_
#define _SX_FS_CONVERT_ACTION_H_

// Including header-files
#include <SxException.h>
#include <SxFileInfo.h>
#include <SxFS.h>

/**
  \brief Enables conversion-operations like move or copy.

  \author Thomas Uchdorf, t.uchdorf@mpie.de
 */
class SX_EXPORT_FS SxFSConvertAction
{
   public:
      SxFSConvertAction ();
      virtual ~SxFSConvertAction ();

      void execute (const SxFileInfo &src, const SxFileInfo &dest);

   protected:

      virtual void convFileToFile (const SxFileInfo &,
                                   const SxFileInfo &) = 0;
      virtual void convFileToDir (const SxFileInfo &,
                                  const SxFileInfo &) = 0;
      virtual void convFileToNothing (const SxFileInfo &,
                                      const SxFileInfo &,
                                      bool force = false) = 0;
      virtual void convFileToFileSymLink (const SxFileInfo &,
                                          const SxFileInfo &) = 0;
      virtual void convFileToDirSymLink (const SxFileInfo &,
                                         const SxFileInfo &) = 0;
      virtual void convFileToEmptySymLink (const SxFileInfo &,
                                           const SxFileInfo &) = 0;

      virtual void convDirToFile (const SxFileInfo &,
                                  const SxFileInfo &) = 0;
      virtual void convDirToDir (const SxFileInfo &,
                                 const SxFileInfo &) = 0;
      virtual void convDirToNothing (const SxFileInfo &,
                                     const SxFileInfo &,
                                     bool force = false) = 0;
      virtual void convDirToFileSymLink (const SxFileInfo &,
                                         const SxFileInfo &) = 0;
      virtual void convDirToDirSymLink (const SxFileInfo &,
                                        const SxFileInfo &) = 0;
      virtual void convDirToEmptySymLink (const SxFileInfo &,
                                          const SxFileInfo &) = 0;

      virtual void convFileSymLinkToFile (const SxFileInfo &,
                                          const SxFileInfo &) = 0;
      virtual void convFileSymLinkToDir (const SxFileInfo &,
                                         const SxFileInfo &) = 0;
      virtual void convFileSymLinkToNothing (const SxFileInfo &,
                                             const SxFileInfo &,
                                             bool force = false) = 0;
      virtual void convFileSymLinkToFileSymLink (const SxFileInfo &,
                                                 const SxFileInfo &) = 0;
      virtual void convFileSymLinkToDirSymLink (const SxFileInfo &,
                                                const SxFileInfo &) = 0;
      virtual void convFileSymLinkToEmptySymLink (const SxFileInfo &,
                                                  const SxFileInfo &) = 0;

      virtual void convDirSymLinkToFile (const SxFileInfo &,
                                         const SxFileInfo &) = 0;
      virtual void convDirSymLinkToDir (const SxFileInfo &,
                                        const SxFileInfo &) = 0;
      virtual void convDirSymLinkToNothing (const SxFileInfo &,
                                            const SxFileInfo &,
                                            bool force = false) = 0;
      virtual void convDirSymLinkToFileSymLink (const SxFileInfo &,
                                                const SxFileInfo &) = 0;
      virtual void convDirSymLinkToDirSymLink (const SxFileInfo &,
                                               const SxFileInfo &) = 0;
      virtual void convDirSymLinkToEmptySymLink (const SxFileInfo &src,
                                                 const SxFileInfo &dest) = 0;

      virtual void convEmptySymLinkToFile (const SxFileInfo &src,
                                           const SxFileInfo &dest) = 0;
      virtual void convEmptySymLinkToDir (const SxFileInfo &src,
                                          const SxFileInfo &dest) = 0;
      virtual void convEmptySymLinkToNothing (const SxFileInfo &src,
                                              const SxFileInfo &dest,
                                              bool force = false) = 0;
      virtual void convEmptySymLinkToFileSymLink (const SxFileInfo &src,
                                                  const SxFileInfo &dest) = 0;
      virtual void convEmptySymLinkToDirSymLink (const SxFileInfo &src,
                                                 const SxFileInfo &dest) = 0;
      virtual void convEmptySymLinkToEmptySymLink (const SxFileInfo &src,
                                                   const SxFileInfo &dest) = 0;

      virtual SxString getName () const = 0;
};

#endif /* _SX_FS_CONVERT_ACTION_H_ */
