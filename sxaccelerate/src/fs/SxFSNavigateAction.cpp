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
#include <SxFSNavigateAction.h>
#include <SxFSError.h>
#include <SxFSAction.h>

#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#ifdef WIN32
#  include <windows.h>
#  include <direct.h>
#else
#  include <dirent.h>
#endif

// Class-members
SxList<SxDir> SxFSNavigateAction::dirStack;  // used for pushDir and popDir

// --- Constructors
// Constructor
SxFSNavigateAction::SxFSNavigateAction ()
{
   SX_DBG_TRACE("");
   // empty
}
// Destructor
SxFSNavigateAction::~SxFSNavigateAction ()
{
   SX_DBG_TRACE("");
   // empty
}

// Class-functions
void SxFSNavigateAction::cd (const SxString &newCWD)
{
   SX_DBG_TRACE("");
   SxFSNavigateAction::cd (SxDir (newCWD));
}

void SxFSNavigateAction::cd (const SxDir &newCWD)
{
   SX_DBG_TRACE("");
   SxString const &newCWDStr = newCWD.getAbsPath ();
   // Performing the directory change
   int error = -1;
#  ifdef WIN32
      if (newCWDStr.isUnicode ())
         error = ::_wchdir ((LPCWSTR)newCWDStr.utf16 ().elements);
      else
         error = ::_chdir  (newCWDStr.ascii ());
#  else
      error = ::chdir (newCWDStr.getElems ());
#  endif

   if (error)  {
      SX_THROW ("Can't change the current directory to '"
                + newCWDStr + "'. " + SxFSError::getChdirErrMsg ());
   }
}

void SxFSNavigateAction::pushd (const SxString &dir)
{
   SX_DBG_TRACE("");
   SxFSNavigateAction::pushd (SxDir (dir));
}

void SxFSNavigateAction::pushd (const SxDir &dir)
{
   SX_DBG_TRACE("");
   SxFSNavigateAction::dirStack << SxFSAction::pwd ();
   SxFSNavigateAction::cd (dir);
}

void SxFSNavigateAction::popd ()
{
   SX_DBG_TRACE("");
   if (SxFSNavigateAction::dirStack.getSize () > 0)  {
      SxFSNavigateAction::cd (SxFSNavigateAction::dirStack.last());
      SxFSNavigateAction::dirStack.removeLast ();
   }
}
