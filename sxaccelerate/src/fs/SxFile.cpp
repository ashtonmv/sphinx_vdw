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

// --- Including header-files
#include <SxFile.h>

// --- Constructors
// Constructor
SxFile::SxFile (const SxString &path_)
   : SxFileInfo (path_)
{
   // empty
}

SxFile::SxFile (const SxFileInfo &in)
   : SxFileInfo (in)
{
   // empty
}

// Destructor
SxFile::~SxFile ()
{
   // empty
}

// --- Methods
SxString SxFile::getIncludeName (const SxString &path_)
{
   if (path_.getSize () == 0) return getName ();
#  ifdef WIN32
      return getIncludeName (path_.tokenize (';'));
#  else
      return getIncludeName (path_.tokenize (':'));
#  endif
}

SxString SxFile::getIncludeName (const SxList<SxString> &path_)
{
   const SxString &tmpAbsPath = getAbsPath ().substitute ('\\', '/');
   SxList<SxString> partialPath = tmpAbsPath.tokenize ('/');
   SxList<SxString>::ConstIterator pathIt, pPathIt = partialPath.fromLast ();
   SxString fileName;
   bool found = false;

   for (pPathIt = partialPath.fromLast (); pPathIt.isValid (); --pPathIt)  {
      if (fileName.getSize () > 0)
         fileName = *pPathIt + "/" + fileName;
      else
         fileName = *pPathIt;

      for (pathIt = path_.begin (); pathIt != path_.end (); ++pathIt)  {
         if (SxFile (*pathIt + "/" + fileName).exists ())  {
            found = true;
            break;
         }
      }
      if (found) break;
   }

   if (!found) return SxString ();
   return fileName;
}

bool SxFile::exists () const
{
   // Checking whether the absolute path leads to an existing file
   return isFile () || isSymLink ();
}
