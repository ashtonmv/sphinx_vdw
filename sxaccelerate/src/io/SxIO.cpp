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

#include <SxUniqueList.h>
#include <SxIO.h>
#include <unistd.h>
#include <SxString.h>
#include <SxMutex.h>

void deleteFileOnce (const SxString &filename)
{
   static SxMutex mutex;
   static SxUniqueList<SxString> allNames;
   mutex.lock ();
   if (! allNames.contains (filename))  {
      allNames.append (filename);
      int err = unlink (filename.ascii ());
      if (err == EISDIR)  {
         cout << "WARNING: output file " << filename << " exists as a directory"
              << endl;
      } else if (err == EPERM)  {
         // strange. According to man page, this may happen if the OS
         // does not allow to unlink files. Let's try to truncate file...
         FILE *fp = fopen(filename.ascii (), "w");
         if (fp)  {
            fclose (fp);
         } else {
            std::cout << "WARNING: Failed to clean " << filename << endl;
         }
      }
      // ignore all other errors.
   }
   mutex.unlock ();
}

void sxfopenError(const char *name, const char *mode)
{
   cout << "Failed to open '" << name << "' for ";
   switch (mode[0])  {
      case 'r': cout << "reading"; break;
      case 'w': cout << "writing"; break;
      case 'a': cout << "appending"; break;
      default: SX_EXIT; // illegal file mode
   }
   cout << endl;
   SX_QUIT;
}
