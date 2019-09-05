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

// Including header-files
#include <SxDir.h>
#include <SxFileIO.h>
#include <SxFSAction.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <SxConfig.h>
#include <SxUtil.h>

#ifdef MACOSX
#  include <mach-o/dyld.h>
#endif

#ifdef WIN32
#  include <windows.h>
#else
#  include <dirent.h>
#endif

// --- Constructors
// Constructor
SxDir::SxDir (const SxString &name_)
   : SxFileInfo (name_)
{
   // empty
}

SxDir::SxDir (const SxFileInfo &rhs)
   : SxFileInfo (rhs)
{
   // empty
}

// Destructor
SxDir::~SxDir ()
{
   // empty
}

// --- Methods
bool SxDir::exists () const
{
   // Checking whether the absolute path leads to an existing directory
   return isDir ();
}


SxString SxDir::getExecPath ()
{
   const size_t BUFLEN = 10240;
   char cPath[BUFLEN];
   sxGetExecPath (cPath, BUFLEN);
   SxString res  (cPath);
   res = res.substitute ('\\', '/');

   // remove name of executable from path
   SxList<SxString> tokens = res.tokenize ('/');
   tokens.removeLast ();

#  ifdef WIN32
      res = "";
#  else
      res = "/";
#  endif /* WIN32 */
   res += SxString::join (tokens, '/');

   SX_CHECK (res != "");

   return res;
}

SxString SxDir::getInstallPath ()
{
   SX_TRACE ();
   // --- try to figure out top directory of package
   SxString installPath, parent;
   SxString execPath = getExecPath ();
   SX_DBG_MSG ("execPath=" << execPath);
#  ifdef WIN32
      // --- cmd.exe environment
      if (SxFSAction::test_d (execPath + "share"))
         installPath = execPath + "/..";
	   else
		   installPath = SxString(SX_SOLUTION_PATH) + "/src";
      return installPath;

#  endif /* WIN32 */
   for (int lvl=0; lvl < 6; ++lvl)  {
      parent += "/..";
      SX_DBG_MSG ("scanning " << execPath << parent << "/share");
      if (SxFSAction::test_d (execPath + parent + "/share"))  {
         installPath = execPath + parent;
         break;
      }
   }
   if (installPath.getSize () == 0 && execPath.contains ("sxaccelerate"))  {
      SxString sxaccelPath = execPath.left ("sxaccelerate");
      if (SxFSAction::test_d (sxaccelPath + "src/share"))  {
         installPath = sxaccelPath + "src";
      }
   }

   if (installPath == "")  {
      SX_THROW("Can't determine program installPath from execPath '"
               + execPath + "'.");
   }

   SxString srcdir = installPath + "/share/.srcdir";
   if (SxFSAction::test_f (srcdir))  {
      installPath = SxFileIO::readLines (srcdir);
      if (installPath.getSize () == 0)  {
         SX_THROW("Can't read program installPath from file '"
                  + srcdir + "': File is empty.");
      }
      if (installPath.contains ("\n"))
         installPath = installPath.left ("\n");
   }
   return installPath;
}
