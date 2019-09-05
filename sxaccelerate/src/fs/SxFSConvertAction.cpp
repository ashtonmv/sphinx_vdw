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

// Including header files
#include <SxFSConvertAction.h>

// --- Constructors
// Constructor
SxFSConvertAction::SxFSConvertAction ()
{
   // empty
}
// Destructor
SxFSConvertAction::~SxFSConvertAction ()
{
   // empty
}

// --- Methods
void SxFSConvertAction::execute (const SxFileInfo &src,
                                 const SxFileInfo &dest)
{
   // Distributing the task of moving source to destination according to the
   // type of the passed objects
   if (src.isSymLink ())  {
      // source is a symbolic link

      if (src.isFile ())  {
         // source is a symbolic link that addresses a file

         if (dest.isSymLink ())  {
            // destination is a symbolic link


            if (dest.isFile ())  {
               // destination is a symbolic link that addresses a file

               convFileSymLinkToFileSymLink (src, dest);
            } else if (dest.isDir ())  {
               // destination is a symbolic link that addresses a directory

               convFileSymLinkToDirSymLink (src, dest);
            } else  {
               // destination is a symbolic link that addresses nothing

               convFileSymLinkToEmptySymLink (src, dest);
            }
         } else if (dest.isFile ())  {
            // destination is a file

            convFileSymLinkToFile (src, dest);
         } else if (dest.isDir ())  {
            // destination is a directory

            convFileSymLinkToDir (src, dest);
         } else  {
            // destination does not exist

            convFileSymLinkToNothing (src, dest);
         }

      } else if (src.isDir ())  {
         // source is a symbolic link that addresses a directory
         if (dest.isSymLink ())  {
            // destination is a symbolic link


            if (dest.isFile ())  {
               // destination is a symbolic link that addresses a file

               convDirSymLinkToFileSymLink (src, dest);
            } else if (dest.isDir ())  {
               // destination is a symbolic link that addresses a directory

               convDirSymLinkToDirSymLink (src, dest);
            } else  {
               // destination is a symbolic link that addresses nothing

               convDirSymLinkToEmptySymLink (src, dest);
            }
         } else if (dest.isFile ())  {
            // destination is a file

            convDirSymLinkToFile (src, dest);
         } else if (dest.isDir ())  {
            // destination is a directory

            convDirSymLinkToDir (src, dest);
         } else  {
            // destination does not exist

            convDirSymLinkToNothing (src, dest);
         }
      } else  {
         // source is a symbolic link that addresses nothing
         if (dest.isSymLink ())  {
            // destination is a symbolic link

            if (dest.isFile ())  {
               // destination is a symbolic link that addresses a file

               convEmptySymLinkToFileSymLink (src, dest);
           } else if (dest.isDir ())  {
               // destination is a symbolic link that addresses a directory

              convEmptySymLinkToDirSymLink (src, dest);
            } else  {
               // destination is a symbolic link that addresses nothing

               convEmptySymLinkToEmptySymLink (src, dest);
            }
         } else if (dest.isFile ())  {
            // destination is a file

            convEmptySymLinkToFile (src, dest);
         } else if (dest.isDir ())  {
            // destination is a directory

            convEmptySymLinkToDir (src, dest);
         } else  {
            // destination does not exist

            convEmptySymLinkToNothing (src, dest);
         }
      }
   } else if (src.isFile ())  {
      // source is a file
      if (dest.isSymLink ())  {
         // destination is a symbolic link


         if (dest.isFile ())  {
            // destination is a symbolic link that addresses a file

            convFileToFileSymLink (src, dest);

         } else if (dest.isDir ())  {
            // destination is a symbolic link that addresses a directory

            convFileToDirSymLink (src, dest);

         } else  {
            // destination is a symbolic link that addresses nothing

            convFileToEmptySymLink (src, dest);
         }
      } else if (dest.isFile ())  {
         // destination is a file

         convFileToFile (src, dest);
      } else if (dest.isDir ())  {
         // destination is a directory

         convFileToDir (src, dest);
      } else  {
         // destination does not exist

         convFileToNothing (src, dest);
      }
   } else if (src.isDir ())  {
      // source is a directory
      if (dest.isSymLink ())  {
         // destination is a symbolic link

         if (dest.isFile ())  {
            // destination is a symbolic link that addresses a file

            convDirToFileSymLink (src, dest);
         } else if (dest.isDir ())  {
            // destination is a symbolic link that addresses a directory

            convDirToDirSymLink (src, dest);
         } else  {
            // destination is a symbolic link that addresses nothing

            convDirToEmptySymLink (src, dest);
         }
      } else if (dest.isFile ())  {
         // destination is a file

         convDirToFile (src, dest);
      } else if (dest.isDir ())  {
         // destination is a directory

         convDirToDir (src, dest);
      } else  {
         // destination does not exist

         convDirToNothing (src, dest);
      }
   } else  {
      // source does not exist, i.e., there is nothing to be moved
      SX_THROW ("Can't " + getName ()
                + " '" + src.getAbsPath () + "' to '"
                + dest.getAbsPath ()
                + "'. There is no file, directory or "
                  "symbolic link with that name.");
   }
}
