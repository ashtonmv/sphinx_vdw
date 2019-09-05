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

// --- Including header files
#include <SxFSAction.h>

#include <SxFSMoveAction.h>
#include <SxFSCopyAction.h>
#include <SxFSNavigateAction.h>
#include <SxFSExploreAction.h>
#include <SxFSAuthorAction.h>
#include <SxFSCreateAction.h>
#include <SxFSDeleteAction.h>

#include <stdlib.h>

// --- Constructors
// Constructor
SxFSAction::SxFSAction ()
{
   // empty
}
// Destructor
SxFSAction::~SxFSAction ()
{
   // empty
}

// --- Class-functions
// --- Authorize-functions
void SxFSAction::chmod (const SxString& target,
                        SxFileInfo::AccessGroup group,
                        SxFileInfo::AccessType type)
{
   SxFSAction::chmod (SxFileInfo (target), group, type);
}

void SxFSAction::chmod (const SxFileInfo& target,
                        SxFileInfo::AccessGroup group,
                        SxFileInfo::AccessType type)
{
   SxFSAuthorAction::chmod (target, group, type);
}

void SxFSAction::chmod (const SxString &path, int mode)
{
   SxFSAuthorAction::chmod (SxFileInfo (path), mode);
}

void SxFSAction::chmod (int perm,
                        const SxString& targetStr)
{
   SxFSAction::chmod (perm, SxFileInfo (targetStr));
}

void SxFSAction::chmod (int perm,
                        const SxFileInfo &target)
{
   SxFSAuthorAction::chmod (perm, target);
}

void SxFSAction::chmod (const SxString &modStr,
                        const SxString &targetStr)
{
   SxFSAction::chmod (modStr, SxFileInfo (targetStr));
}

void SxFSAction::chmod (const SxString   &modStr,
                        const SxFileInfo &target)
{
   SxFSAuthorAction::chmod (modStr, target);
}

void SxFSAction::chgrp (unsigned int gid,
                        const SxString &target)
{
   SxFSAction::chgrp (gid, SxFileInfo (target));
}

void SxFSAction::chgrp (unsigned int gid,
                        const SxFileInfo& target)
{
   SxFSAuthorAction::chgrp (target, gid);
}

void SxFSAction::chown (unsigned int uid,
                        const SxString& target)
{
   SxFSAction::chown (uid, SxFileInfo (target));
}

void SxFSAction::chown (unsigned int uid,
                        const SxFileInfo& target)
{
   SxFSAuthorAction::chown (target, uid);
}

// --- Convert-functions
void SxFSAction::mv (const SxString &src,
                     const SxString &dest)
{
   SxFSAction::mv (SxFileInfo (src), SxFileInfo (dest));
}

void SxFSAction::mv (const SxFileInfo &src,
                     const SxFileInfo &dest)
{
   SxPtr<SxFSMoveAction> mvAct = SxPtr<SxFSMoveAction>::create ();
   mvAct->execute (src, dest);
}

void SxFSAction::cp (const SxString &src,
                     const SxString &dest)
{
   SxFSAction::cp (SxFileInfo (src), SxFileInfo (dest));
}

void SxFSAction::cp (const SxFileInfo &src,
                     const SxFileInfo &dest)
{
   SxPtr<SxFSCopyAction> cpAct = SxPtr<SxFSCopyAction>::create ();
   cpAct->execute (src, dest);
}

// --- Create-functions
void SxFSAction::touch (const SxString &target, int mode)
{
   touch (SxFileInfo (target), mode);
}

void SxFSAction::touch (const SxFileInfo &target, int mode)
{
  SxFSCreateAction::touch (target, mode);
}

SxFileInfo SxFSAction::ln_sf (const SxString &path,
                              const SxString &newLink)
{
   return SxFSAction::ln_sf (SxFileInfo (path), SxSymLink (newLink));
}

SxFileInfo SxFSAction::ln_sf (const SxFileInfo &path,
                              const SxFileInfo &newLink)
{
   if (SxFSAction::test_L (newLink.getAbsPath ()))  {
      SxFSAction::rm_r (newLink.getAbsPath ());
   }
   return SxFSCreateAction::ln_sf (path, newLink);
}

SxFile SxFSAction::createTmpFile (const SxString &tmpDir,
                                  const SxArray<unsigned char> &
                                  buffer)
{
   return SxFSCreateAction::createTmpFile (tmpDir, buffer);
}

void SxFSAction::mkdir (const SxString &subFolder,
                        int umask)
{
   SxFSAction::mkdir (SxDir (subFolder), umask);
}

void SxFSAction::mkdir (const SxDir &subFolder,
                        int umask)
{
   SxFSCreateAction::mkdir (subFolder, umask);
}

void SxFSAction::mkdir_p (const SxString &destDir,
                          int umask)
{
   SX_CHECK (destDir.removeWhiteSpace() != "");
   SX_CHECK (!destDir.contains ("***"), destDir); // avoid markers
   SxFSAction::mkdir_p (SxDir (destDir), umask);
}

void SxFSAction::mkdir_p (const SxDir &destDir,
                          int umask)
{
   SX_CHECK (destDir.getPath ().removeWhiteSpace () != "");
   SX_CHECK (!destDir.getPath ().contains ("***"), destDir); // avoid markers
   SxFSCreateAction::mkdir_p (destDir, umask);
}

// --- Delete-functions
void SxFSAction::rm (const SxString &target)
{
   if (target == "")  {
      SX_THROW ("rm(pathname): Pathname can't be empty string.");
   }
   SxFSAction::rm (SxFile (target));
}

void SxFSAction::rm (const SxFile &target)
{
   SxFSDeleteAction::rm (target);
}

void SxFSAction::rmdir (const SxString &path)
{
   if (path == "")  {
      SX_THROW ("rmdir(pathname): Pathname can't be empty string.");
   }
   SxFSAction::rmdir (SxDir (path));
}

void SxFSAction::rmdir (const SxDir &path)
{
   SxFSDeleteAction::rmdir (path);
}

void SxFSAction::rm_r (const SxString &target)
{
   if (target == "")  {
      SX_THROW("rm_r(pathname): Pathname can't be empty string.");
   }
   SxFSAction::rm_r (SxFileInfo (target));
}

void SxFSAction::rm_r (const SxFileInfo &target)
{
   SxFSDeleteAction::rm_r (target);
}

// --- Navigate-functions
void SxFSAction::cd (const SxString &newCWD)
{
   SxFSAction::cd (SxDir (newCWD));
}

void SxFSAction::cd (const SxDir &newCWD)
{
   SxFSNavigateAction::cd (newCWD);
}

void SxFSAction::pushd (const SxString &dir)
{
   SxFSAction::pushd (SxDir (dir));
}

void SxFSAction::pushd (const SxDir &dir)
{
   SxFSNavigateAction::pushd (SxDir (dir));
}

void SxFSAction::popd ()
{
   SxFSNavigateAction::popd ();
}

// --- Explore-functions
SxDir SxFSAction::getHome ()
{
   return SxFSExploreAction::getHome ();
}

SxDir SxFSAction::getTmp ()
{
   return SxFSExploreAction::getTmp ();
}

SxDir SxFSAction::pwd ()
{
   return SxFSExploreAction::pwd ();
}

SxList<SxFileInfo> SxFSAction::find (const SxString &name, bool chdir_)
{
   SX_DBG_TRACE("");

   // --- support of "dirA/dirB/*.ext"
   SxList<SxString> patternList = name.tokenize ('/');
   SxString filePattern = patternList.last (); patternList.removeLast ();
   SxString dir = SxString::join (patternList, '/');
   SX_CHECK (!dir.contains ("*"));  // patterns in folders not yet supported
   if (name(0) == '/')  dir = "/" + dir;

   SxList<SxFileInfo> res;

   if (chdir_)  {
      pushd (dir);
      res = SxFSExploreAction::find (filePattern);
      popd ();
   }  else  {
      res = SxFSExploreAction::find (filePattern, SxDir(dir));
   }

   return res;
}

SxList<SxFISortedByTime> SxFSAction::ls_t (const SxString &path)
{
   SX_DBG_TRACE("");
   return SxFSAction::ls_t (SxFileInfo (path));
}

SxList<SxFISortedByTime> SxFSAction::ls_t (const SxFileInfo &path)
{
   SX_DBG_TRACE("");
   return SxFSExploreAction::ls_t (path);
}

SxList<SxFISortedBySize> SxFSAction::ls_S (const SxString &path)
{
   SX_DBG_TRACE("");
   return SxFSAction::ls_S (SxFileInfo (path));
}

SxList<SxFISortedBySize> SxFSAction::ls_S (const SxFileInfo &path)
{
   SX_DBG_TRACE("");
   return SxFSExploreAction::ls_S (path);
}

SxList<SxFileInfo> SxFSAction::ls (const SxString &path)
{
   return SxFSAction::ls (SxFileInfo (path));
}

SxList<SxFileInfo> SxFSAction::ls (const SxFileInfo &src)
{
   return SxFSExploreAction::ls (src);
}

SxString SxFSAction::which (const SxString &file, const SxString &path_)
{
   SxList<SxString> files = where (file, path_);
   if (files.getSize() == 0)  return "";
   else                       return files.first ();
}

SxList<SxString> SxFSAction::where (const SxString &file,
                                    const SxString &path_)
{
   SxList<SxString> res;
   SxFile f(file);
   if (f.exists ())  res << f.getAbsPath();

   SxString path = path_;
   if (path == "")  path = ::getenv ("PATH");
#  ifdef WIN32
      // unquote sub paths in windows
      path = path.substitute ("\"", "");
#  endif

#  ifdef WIN32
      char sep = ';';
#  else
      char sep = ':';
#  endif
   path.tokenize (sep).foreach ([&](auto it) {
      f = SxFile (*it + '/' + file);
      if (f.exists ())  {
         res << f.getAbsPath ();
      } else {
#        ifdef WIN32
            f = SxFile (*it + '/' + file + ".exe");
            if (f.exists ())   res << f.getAbsPath ();
#        endif
      }
   });
   return res;
}


// Collect all filenames in given directory. Recursive.
SxList<SxString> SxFSAction::listDir (const SxString &path_)
{
   SxList<SxString> result;

   SxFSExploreAction::listDir (path_, &result);

   return result;
}

SxList<SxFile> SxFSAction::getFiles (const SxString &path)
{
   return SxFSAction::getFiles (SxFileInfo (path));
}

SxList<SxFile> SxFSAction::getFiles (const SxFileInfo &path)
{
   return SxFSExploreAction::getFiles (path);
}

SxList<SxDir> SxFSAction::getDirs (const SxString &path)
{
   return SxFSAction::getDirs (SxFileInfo (path));
}

SxList<SxDir> SxFSAction::getDirs (const SxFileInfo &path)
{
   return SxFSExploreAction::getDirs (path);
}

SxList<SxSymLink> SxFSAction::getSymLinks (const SxString &path)
{
   return SxFSAction::getSymLinks (SxFileInfo (path));
}

SxList<SxSymLink> SxFSAction::getSymLinks (const SxFileInfo &path)
{
   return SxFSExploreAction::getSymLinks (path);
}

// --- Meta-functions (for convenience)
bool SxFSAction::exists (const SxString &str)
{
   return SxFileInfo (str).exists ();
}

bool SxFSAction::test_f (const SxString &str)
{
   return SxFile (str).exists ();
}

bool SxFSAction::test_d (const SxString &str)
{
   return SxDir (str).exists ();
}

bool SxFSAction::test_L (const SxString &str)
{
   return SxSymLink (str).exists ();
}
