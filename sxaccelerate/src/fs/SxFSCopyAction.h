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

#ifndef _SX_FS_COPY_ACTION_H_
#define _SX_FS_COPY_ACTION_H_

// Including header-files
#include <SxFSConvertAction.h>
#include <SxFileInfo.h>
#include <SxFS.h>

/**
  \brief Can be used to perform copy-operations

  \author Thomas Uchdorf, t.uchdorf@mpie.de
 */
class SX_EXPORT_FS SxFSCopyAction : public SxFSConvertAction
{
   public:
      SxFSCopyAction ();
      virtual ~SxFSCopyAction ();

   protected:

      virtual void convFileToFile (const SxFileInfo &,
                                   const SxFileInfo &);
      virtual void convFileToDir (const SxFileInfo &,
                                  const SxFileInfo &);
      virtual void convFileToNothing (const SxFileInfo &,
                                      const SxFileInfo &,
                                      bool force = false);
      virtual void convFileToFileSymLink (const SxFileInfo &,
                                          const SxFileInfo &);
      virtual void convFileToDirSymLink (const SxFileInfo &,
                                         const SxFileInfo &);
      virtual void convFileToEmptySymLink (const SxFileInfo &,
                                           const SxFileInfo &);

      virtual void convDirToFile (const SxFileInfo &,
                                  const SxFileInfo &);
      virtual void convDirToDir (const SxFileInfo &,
                                 const SxFileInfo &);
      virtual void convDirToNothing (const SxFileInfo &,
                                     const SxFileInfo &,
                                     bool force = false);
      virtual void convDirToFileSymLink (const SxFileInfo &,
                                         const SxFileInfo &);
      virtual void convDirToDirSymLink (const SxFileInfo &,
                                        const SxFileInfo &);
      virtual void convDirToEmptySymLink (const SxFileInfo &,
                                          const SxFileInfo &);

      virtual void convFileSymLinkToFile (const SxFileInfo &,
                                          const SxFileInfo &);
      virtual void convFileSymLinkToDir (const SxFileInfo &,
                                         const SxFileInfo &);
      virtual void convFileSymLinkToNothing (const SxFileInfo &,
                                             const SxFileInfo &,
                                             bool force = false);
      virtual void convFileSymLinkToFileSymLink (const SxFileInfo &,
                                                 const SxFileInfo &);
      virtual void convFileSymLinkToDirSymLink (const SxFileInfo &,
                                                const SxFileInfo &);
      virtual void convFileSymLinkToEmptySymLink (const SxFileInfo &,
                                                  const SxFileInfo &);

      virtual void convDirSymLinkToFile (const SxFileInfo &,
                                         const SxFileInfo &);
      virtual void convDirSymLinkToDir (const SxFileInfo &,
                                        const SxFileInfo &);
      virtual void convDirSymLinkToNothing (const SxFileInfo &,
                                            const SxFileInfo &,
                                            bool force = false);
      virtual void convDirSymLinkToFileSymLink (const SxFileInfo &,
                                                const SxFileInfo &);
      virtual void convDirSymLinkToDirSymLink (const SxFileInfo &,
                                               const SxFileInfo &);
      virtual void convDirSymLinkToEmptySymLink (const SxFileInfo &src,
                                                 const SxFileInfo &dest);

      virtual void convEmptySymLinkToFile (const SxFileInfo &src,
                                           const SxFileInfo &dest);
      virtual void convEmptySymLinkToDir (const SxFileInfo &src,
                                          const SxFileInfo &dest);
      virtual void convEmptySymLinkToNothing (const SxFileInfo &src,
                                              const SxFileInfo &dest,
                                              bool force = false);
      virtual void convEmptySymLinkToFileSymLink (const SxFileInfo &src,
                                                  const SxFileInfo &dest);
      virtual void convEmptySymLinkToDirSymLink (const SxFileInfo &src,
                                                 const SxFileInfo &dest);
      virtual void convEmptySymLinkToEmptySymLink (const SxFileInfo &src,
                                                   const SxFileInfo &dest);

      void rmDestSymLinkAndCpEmptySymLink (const SxFileInfo &src,
                                           const SxFileInfo &dest);

      virtual SxString getName () const;
};

#endif /* _SX_FS_COPY_ACTION_H_ */
