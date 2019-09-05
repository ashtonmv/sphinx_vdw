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

// --- Including header-filesf
#include <SxFileInfo.h>
#include <SxSymLink.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#ifdef WIN32
#  include <shortcut.h>
#  include <windows.h>
#  include <direct.h>
#endif /* WIN32 */
#include <SxDir.h>

// --- Constructors
// Constructor
SxFileInfo::SxFileInfo (const SxString &path_)
{
   SX_DBG_TRACE("");
   cache = SxPtr<Cache>::create ();

   if (path_ == "")  {
      // --- Fetching the current working directory and storing it as absolute
      //     path
      char buffer[MAX_PATH_LEN];
      if(!getcwd (buffer, MAX_PATH_LEN))  {
         SX_THROW ("Cannot get current working directory. "
                   "getcwd() failed: " + sxstrerror ());
      }
      SX_CHECK (cache);
      cache->setAbsPath (buffer);
      orig = getAbsPath ();
   } else  {
      SxString specificPath;

#ifdef WIN32
      specificPath = path_.substitute ('\\', '/');
#else
      specificPath = path_;
#endif /* WIN32 */

      if (path_.isUnicode ()) specificPath.setUnicode ();

      // The path parameter was defined
      orig = specificPath;
	   //std::cout << "Checking if \"" << specificPath << "\" is absolute..." << std::endl;//TEST
      // Testing if the path is relative or absolute
      if (SxFileInfo::isAbs (specificPath))  {
         // An absolute path has been specified

		   //std::cout << "TRUE" << std::endl;//TEST
         SX_CHECK (cache);
         // Storing the path
         cache->setAbsPath (specificPath);
      } else  {
         // The passed path is relative

		   //std::cout << "FALSE" << std::endl;//TEST
         // --- Getting the current working directory and appending the
         //     transmitted path

#        ifdef WIN32
            uint16_t buffer[MAX_PATH_LEN];
            if (!_wgetcwd ((wchar_t *)buffer, MAX_PATH_LEN))  {
               SX_THROW ("Cannot get current working directory. "
                         "_wgetcwd() failed: " + sxstrerror ());
            }

            SX_CHECK (cache);
            if (specificPath == ".")  {
               cache->setAbsPath (SxString::fromUtf16 (buffer));
            } else  {
               cache->setAbsPath (SxString::fromUtf16 (buffer) +
                                  SxString ('/')
                                  + specificPath);
           
            }
#        else
            char buffer[MAX_PATH_LEN];
            if (!getcwd (buffer, MAX_PATH_LEN))  {
               SX_THROW ("Cannot get current working directory. "
                         "getcwd() failed: " + sxstrerror ());
            }

            SX_CHECK (cache);
            if (specificPath == ".")  {
               cache->setAbsPath (SxString (buffer));
            } else  {
               cache->setAbsPath (SxString (buffer) +
                                  SxString ('/')
                                  + specificPath);

            }
#        endif

      }
   }
}

SxFileInfo::SxFileInfo (const SxFileInfo &rhs)
  : orig (rhs.orig), cache (rhs.cache)
{
   SX_DBG_TRACE("");
   // empty
}

// Destructor
SxFileInfo::~SxFileInfo ()
{
   SX_DBG_TRACE("");
   // empty
}

// --- Methods
void SxFileInfo::setDirty () const
{
   SX_DBG_TRACE("");
   cache->setDirty ();
}

char SxFileInfo::getSeparator ()
{
   SX_DBG_TRACE("");
#ifdef WIN32
   return '\\';
#else
   return '/';
#endif /* WIN32 */
}

bool SxFileInfo::exists () const
{
   SX_DBG_TRACE("");
   // --- Testing if the file or directory specified by the stored absolute
   //     path exists
   SX_CHECK (cache);
   return cache->isExisting ();
}

bool SxFileInfo::isFile () const
{
   SX_DBG_TRACE("");
   // Checking whether the absolute path leads to a file or directory at all
   // and testing if the element defined by the absolute path is a regular file
   SX_CHECK (cache);
   return (cache->isExisting () && (cache->getType () & Cache::File));
}

bool SxFileInfo::isSymLink () const
{
   SX_DBG_TRACE("");
   // Willingly not checking whether the absolute path leads to a file or
   // directory at all, as we would like to find all symbolic links (even
   // those that actually point to nowhere)!
   SX_CHECK (cache);
   return (cache->getType () & Cache::SymLink) != 0;
}

bool SxFileInfo::isDir () const
{
   SX_DBG_TRACE("");
   if (!SxFileInfo::exists ())  return false;

   // Testing if the element defined by the absolute path is a directory
   SX_CHECK (cache);
   return (cache->isExisting () && (cache->getType () & Cache::Dir));
}

bool SxFileInfo::isReadable (AccessGroup group) const
{
   SX_DBG_TRACE("");
   SX_CHECK (cache);
   int mode = cache->getMode ();
   if (cache->getStatus ())  return false;

   int mask = (group & Owner) ? S_IRUSR : 0;
#ifndef WIN32
   mask |= (group & Group) ? S_IRGRP : 0 | (group & Other) ? S_IROTH : 0;
#endif /* not WIN32 */

   // Testing if the element defined by the absolute path is readable
   return (mode & mask) != 0;
}

bool SxFileInfo::isWritable (AccessGroup group) const
{
   SX_DBG_TRACE("");
   SX_CHECK (cache);
   int mode = cache->getMode ();
   if (cache->getStatus ())  return false;

   int mask = (group & Owner) ? S_IWUSR : 0;
#ifndef WIN32
   mask |= (group & Group) ? S_IWGRP : 0 | (group & Other) ? S_IWOTH : 0;
#endif /* not WIN32 */

   // Testing if the element defined by the absolute path is writable
   return (mode & mask) != 0;
}

bool SxFileInfo::isExecutable (AccessGroup group) const
{
   SX_DBG_TRACE("");
   SX_CHECK (cache);
#  ifdef WIN32
      return true;
#  else
      int mode = cache->getMode ();
      if (cache->getStatus ())  return false;
      int mask = (group & Owner) ? S_IXUSR : 0;
      mask |= (group & Group) ? S_IXGRP : 0 | (group & Other) ? S_IXOTH : 0;

      // Testing if the element defined by the absolute path is executable
      return (mode & mask);
#  endif /* WIN32 */
}


SxString SxFileInfo::getPerms () const
{
   SxString res;
   res += isReadable   (Owner) ? "r" : "-";
   res += isWritable   (Owner) ? "w" : "-";
   res += isExecutable (Owner) ? "x" : "-";
   res += isReadable   (Group) ? "r" : "-";
   res += isWritable   (Group) ? "w" : "-";
   res += isExecutable (Group) ? "x" : "-";
   res += isReadable   (Other) ? "r" : "-";
   res += isWritable   (Other) ? "w" : "-";
   res += isExecutable (Other) ? "x" : "-";
   return res;
}

int SxFileInfo::getMode () const
{
   int res = 0;
   if (isReadable   (Owner))  res += 0400;
   if (isWritable   (Owner))  res += 0200;
   if (isExecutable (Owner))  res += 0100;
   if (isReadable   (Group))  res += 0040;
   if (isWritable   (Group))  res += 0020;
   if (isExecutable (Group))  res += 0010;
   if (isReadable   (Other))  res += 0004;
   if (isWritable   (Other))  res += 0002;
   if (isExecutable (Other))  res += 0001;
   return res;
}

int SxFileInfo::getMode (const SxString &perms)
{
   SX_CHECK (perms.getSize () == 9, perms);
   int res = 0;
   if (perms(0) == 'r')  res += 0400;
   if (perms(1) == 'w')  res += 0200;
   if (perms(2) == 'x')  res += 0100;
   if (perms(3) == 'r')  res += 0040;
   if (perms(4) == 'w')  res += 0020;
   if (perms(5) == 'x')  res += 0010;
   if (perms(6) == 'r')  res += 0004;
   if (perms(7) == 'w')  res += 0002;
   if (perms(8) == 'x')  res += 0001;
   return res;
}

SxString SxFileInfo::getPerms (int m)
{
   SxString res = "---------";
   if (m & 0400)  res(0) = 'r';
   if (m & 0200)  res(1) = 'w';
   if (m & 0100)  res(2) = 'x';
   if (m & 0040)  res(3) = 'r';
   if (m & 0020)  res(4) = 'w';
   if (m & 0010)  res(5) = 'x';
   if (m & 0004)  res(6) = 'r';
   if (m & 0002)  res(7) = 'w';
   if (m & 0001)  res(8) = 'x';
   return res;
}


unsigned long SxFileInfo::lastAccessed () const
{
   SX_DBG_TRACE("");
   SX_CHECK (cache);
   return cache->getAccessTime ();
}

unsigned long long SxFileInfo::lastAccessedMS () const
{
   SX_DBG_TRACE("");
   SX_CHECK (cache);
   return cache->getAccessTimeMS ();
}

unsigned long SxFileInfo::lastModified () const
{
   SX_DBG_TRACE("");
   SX_CHECK (cache);
   return cache->getModificationTime ();
}

unsigned long long SxFileInfo::lastModifiedMS () const
{
   SX_DBG_TRACE("");
   SX_CHECK (cache);
   return cache->getModificationTimeMS ();
}


SxString SxFileInfo::getName (const SxString &path_)
{
   SxString path = path_.substitute ('\\', '/');
   if (path.isUnicode ()) {

      if (path.getSize () == 0) return SxString ();

      ssize_t searchStart = path.getSize () - 1;
      if (path.charAt (searchStart) == (uint32_t)'/') {
         searchStart--;
      }

      ssize_t lastSeparator = path.findLast ('/',
                                             searchStart);

      if (lastSeparator > 1)
         return path.subString (lastSeparator + 1, searchStart);

   } else {
      ssize_t i = path.getSize () - 1;
      const char *str = path.getElems ();
      // --- skip slashes at the end '...b///'
      while (i >= 0 && (str[i] == '\\' || str[i] == '/')) i--;

      // --- windows drive name is not file name
      if (i == 1 && str[0] != '/' && str[1] == ':') return SxString();

      // --- find index of the last slash '.../b'
      ssize_t j = i;
      while (i >= 0 && str[i] != '\\' && str[i] != '/') i--;

      if (i >= 0) return path.subStringByBytes (i + 1, j);
      if (j >= 0) return path.subStringByBytes (0, j);
   }

   return SxString ();
}

SxString SxFileInfo::getPath (const SxString &path_)
{
   SxString path = path_.substitute ('\\', '/');
   if (path.isUnicode ()) {

      if (path.getSize () == 0) return SxString ();

      ssize_t searchStart = path.getSize () - 1;
      if (path.charAt (searchStart) == (uint32_t) '/') {
         searchStart--;
      }

      ssize_t lastSeparator = path.findLast ('/',
                                             searchStart);
      // Linux '/' root path
      if (lastSeparator == 0) return SxString ('/');

      // Windows 'C:\\' root path
      if (  lastSeparator == 2
         && path(1) == ':') return path.subString (0, lastSeparator);

      if (lastSeparator > 1)
         return path.subString (0, lastSeparator - 1);

   } else {

      ssize_t i = path.getSize () - 1;
      const char *str = path.getElems ();
      // --- skip slashes at the end '...b///'
      while (i >= 0 && (str[i] == '\\' || str[i] == '/')) i--;

      // --- '/' or 'C:\' keep the last slash in path
      if (i < 0 && path.getSize () > 0) return path(0);
      if (i == 1 && str[0] != '/' && str[1] == ':')  {
         if (str[2] == '\\' || str[2] == '/') return path.subStringByBytes (0, 2);
         else                                 return SxString(str[0]) + ":\\";
      }

      // --- find index of the last slash '.../b'
      while (i >= 0 && str[i] != '\\' && str[i] != '/') i--;

      // --- '/' or 'C:\' keep the last slash in path
      if (i == 2 && str[0] != '/' && str[1] == ':')
         return path.subStringByBytes (0, i);
      if (i == 0) return path(0);

      // --- return path without the last slash
      if (i > 0) return path.subStringByBytes (0, i - 1);
   }

   return SxString ();
}

SxString SxFileInfo::getName () const
{
   return SxFileInfo::getName (getAbsPath ());
}

SxString SxFileInfo::getPath () const
{
   return SxFileInfo::getPath (getAbsPath ());
}

SxString SxFileInfo::getAbsURI () const
{
   SX_DBG_TRACE("");
   //TODO Implement correctly
   return "file://" + getAbsPath ();
}

SxString const & SxFileInfo::getAbsPath () const
{
   SX_DBG_TRACE("");
   SX_CHECK (cache);
   // Returning the absolute path
   return cache->getAbsPath ();
}


SxString SxFileInfo::getRelPath (const SxString &path) const
{
   SX_DBG_TRACE("");
   SxString ret;
   if (path != "")  {
      ret = getRelPath (SxDir (path));
   }
   return ret;
}

SxString SxFileInfo::getRelPath (const SxDir &root_) const
{
   SX_DBG_TRACE("");
   const SxString root = root_.getAbsPath ();
   const SxString path = getAbsPath ();
   SxString res;

   if (path.find (root) == 0)  {
      res = path.right (root);
      if (res.head(1) == '/')
         res = res.right ('/');
   }
   return res;
}

int64_t SxFileInfo::getSize () const
{
   SX_DBG_TRACE("");
   SX_CHECK (cache);
   return cache->getSize ();
}

sxuid_t SxFileInfo::getUID () const
{
   SX_DBG_TRACE("");
   SX_CHECK (cache);
   return cache->getUID ();
}

sxgid_t SxFileInfo::getGID () const
{
   SX_DBG_TRACE("");
   SX_CHECK (cache);
   return cache->getGID ();
}

//SxString SxFileInfo::getModeStr () const
//{
//   SX_DBG_TRACE("");
//   SxString ret = "---------";
//   const ssize_t & s = ret.getSize ();
//   int mode = getMode ();
//   char c;
//   int mod;
//   for (int i = 0; i < s; ++i)  {
//      if ((mod = (i%3)) == 0)  {
//         c = 'r';
//      } else if (mod == 1)  {
//         c = 'w';
//      } else  {
//         c = 'x';
//      }
//      if (mode&(1<<(s-1-i))) ret(i) = c;
//   }
//   return ret;
//}

// --- Protected class-functions
bool SxFileInfo::equals (const SxString &lhs, const SxString &rhs)
{
   SX_DBG_TRACE("");
   return SxFileInfo::equals (SxFileInfo (lhs), SxFileInfo (rhs));
}

// --- Protected class-functions
bool SxFileInfo::equals (const SxFileInfo &lhs, const SxFileInfo &rhs)
{
   SX_DBG_TRACE("");
   SxString const& lhsStr = lhs.getAbsPath ();
   SxString const& rhsStr = rhs.getAbsPath ();

#ifdef WIN32
   // Opening the file passed as left hand parameter to retrieve pieces of
   // information on it
   bool result = false;
   SxString errorMsg;
   HANDLE rhsFileHandle = INVALID_HANDLE_VALUE;
   HANDLE lhsFileHandle = INVALID_HANDLE_VALUE;
   if (lhsStr.isUnicode ())  {
      lhsFileHandle = CreateFileW ((LPCWSTR)lhsStr.utf16 ().elements,
                                   0, 0, NULL,
                                   OPEN_EXISTING,
                                   ((lhs.isDir ()) ?
                                   FILE_FLAG_BACKUP_SEMANTICS :
                                   FILE_ATTRIBUTE_NORMAL),
                                   NULL);
   } else {
      lhsFileHandle = CreateFileA (lhsStr.ascii (),
                                   0, 0, NULL,
                                   OPEN_EXISTING,
                                   ((lhs.isDir ()) ?
                                   FILE_FLAG_BACKUP_SEMANTICS :
                                   FILE_ATTRIBUTE_NORMAL),
                                   NULL);
   }

   if (lhsFileHandle != INVALID_HANDLE_VALUE)  {
      // --- Trying to fetch meta information on the file
      BY_HANDLE_FILE_INFORMATION lhsBuf;
      if (GetFileInformationByHandle (lhsFileHandle, &lhsBuf))  {
         DWORD lhsDev = lhsBuf.dwVolumeSerialNumber;
         DWORD lhsLowerId = lhsBuf.nFileIndexLow;
         DWORD lhsHigherId = lhsBuf.nFileIndexHigh;

         // Opening the file passed as left hand parameter to retrieve pieces of
         // information on it
         rhsFileHandle = INVALID_HANDLE_VALUE;

         if (rhsStr.isUnicode ())  {
            rhsFileHandle = CreateFileW ((LPCWSTR)rhsStr.utf16 ().elements,
                                         0, 0, NULL,
                                         OPEN_EXISTING,
                                         ((rhs.isDir ()) ?
                                         FILE_FLAG_BACKUP_SEMANTICS :
                                         FILE_ATTRIBUTE_NORMAL),
                                         NULL);
         } else {
            rhsFileHandle = CreateFileA (rhsStr.ascii (),
                                         0, 0, NULL,
                                         OPEN_EXISTING,
                                         ((rhs.isDir ()) ?
                                         FILE_FLAG_BACKUP_SEMANTICS :
                                         FILE_ATTRIBUTE_NORMAL),
                                         NULL);
         }

         if (rhsFileHandle != INVALID_HANDLE_VALUE)  {
            // --- Trying to fetch meta information on the file
            BY_HANDLE_FILE_INFORMATION rhsBuf;
            if (GetFileInformationByHandle (rhsFileHandle, &rhsBuf))  {
               DWORD rhsDev = rhsBuf.dwVolumeSerialNumber;
               DWORD rhsLowerId = rhsBuf.nFileIndexLow;
               DWORD rhsHigherId = rhsBuf.nFileIndexHigh;
               // comparison of the unique identifying device number and the parted ID
               result = (lhsDev == rhsDev && lhsHigherId == rhsHigherId && lhsLowerId == rhsLowerId);
            } else  {
               errorMsg = "Cannot get information on '" + rhsStr + "'.";
            }
         } else  {
            errorMsg = "Cannot get the file handle of '" + rhsStr + "'.";
         }
      } else  {
         errorMsg = "Cannot get information on '" + lhsStr + "'.";
      }
   } else  {
      errorMsg = "Cannot get the file handle of '" + lhsStr + "'.";
   }
   if (lhsFileHandle != INVALID_HANDLE_VALUE)  CloseHandle (lhsFileHandle);
   if (rhsFileHandle != INVALID_HANDLE_VALUE)  CloseHandle (rhsFileHandle);
   if (errorMsg != "")  {
      SX_THROW (errorMsg);
   }
   return result;
#else
   // --- Fetching the device ID and the inode number of the left hand side
   struct stat lhsbuf;
   struct stat rhsbuf;
   if (lstat (lhsStr.getElems (), &lhsbuf))  {
      return false;
   }
   dev_t lhsDev = lhsbuf.st_dev;
   ino_t lhsIno = lhsbuf.st_ino;
   // --- Retrieving the device ID and the inode number of the right hand side
   if (lstat (rhsStr.getElems (), &rhsbuf))  {
      return false;
   }
   // Returning the comparison of the device IDs and the inode numbers
   return lhsDev == rhsbuf.st_dev && lhsIno == rhsbuf.st_ino;
#endif /* WIN32 */
}

bool SxFileInfo::operator== (const SxFileInfo &compare) const
{
   return equals (*this, compare);
}
bool SxFileInfo::operator!= (const SxFileInfo &compare) const
{
   return !equals (*this, compare);
}

int SxFileInfo::compareModificationTime (const SxFileInfo &lhs,
                                         const SxFileInfo &rhs)
{
   SX_DBG_TRACE("");
   unsigned long long lhsTime = lhs.lastModifiedMS ();
   unsigned long long rhsTime = rhs.lastModifiedMS ();

   if (lhsTime < rhsTime)  {
      return -1;
   } else if (lhsTime > rhsTime)  {
      return 1;
   } else {
      return 0;
   }
}

int SxFileInfo::compareSize (const SxFileInfo &lhs,
                             const SxFileInfo &rhs)
{
   SX_DBG_TRACE("");
   int64_t lhsSize = lhs.getSize ();
   int64_t rhsSize = rhs.getSize ();
   if (lhsSize < rhsSize)  {
      return -1;
   } else if (lhsSize > rhsSize)  {
      return 1;
   } else {
      return 0;
   }
}

bool SxFileInfo::isAbs (const SxString &str)
{
#ifdef WIN32
   return str.find (":") == 1; // 'X:'
#else
   return str.find ("/") == 0;
#endif /* WIN32 */
}

SxString SxFileInfo::getHomeStr ()
{
   SX_DBG_TRACE("");
   return SxString (getenv("HOME")).substitute ('\\', '/');
}

SxString SxFileInfo::getTmpStr ()
{
   SX_DBG_TRACE("");
#  ifndef WIN32
      return SxString (getenv ("TMPDIR")).substitute ('\\', '/');
#  else
      return SxString (getenv ("TMP")).substitute ('\\', '/');
#  endif
}

// --- Protected methods
SxString const & SxFileInfo::getOrig () const
{
   SX_DBG_TRACE("");
   // Returning the originally set path
   return orig;
}

SxString SxFileInfo::simplifyPath (const SxString &path)
{
   SxString res = path;
   while (res.contains ("\\\\"))  res = res.substitute ("\\\\", "/");
   while (res.contains ("//"))    res = res.substitute ("//", "/");
   return res;
}

SxString SxFileInfo::resolvePath (const SxString &path)
{
   SxList<SxString> tokens = path.tokenize ("/\\");
   SxList<SxString> tokensPath;
   SxList<SxString>::ConstIterator it;
   for (it = tokens.begin (); it != tokens.end (); ++it)  {
      if (*it == "..")  {
         tokensPath.removeLast ();
      }  else if (*it != ".")  {
         tokensPath << (*it);
      }
   }

   SxString res = SxString::join (tokensPath, '/');
   if (path.getSize() > 0 && path(0) == '/') res.prepend ('/');

   return res;
}

//------------------------------------------------------------------------------
// SxFileInfo::Cache
//------------------------------------------------------------------------------

// --- Constructors
// Constructor
SxFileInfo::Cache::Cache ()
  : dirty (true), type (Nothing), existing (false), status (0),
   userID ((sxuid_t)-1), groupID ((sxgid_t)-1),
   mode (0), size (0),
   accessTime (0),accessTimeMS (0),
   modificationTime (0), modificationTimeMS (0)
{
   /* empty */
}

SxFileInfo::Cache::Cache (const Cache &rhs)
  : dirty (rhs.dirty), type (rhs.type), existing (rhs.existing),
   status (rhs.status),
   userID (rhs.userID), groupID (rhs.groupID),
   mode (rhs.mode), size (rhs.size),
   accessTime (rhs.accessTime),
   accessTimeMS (rhs.accessTimeMS),
   modificationTime (rhs.modificationTime),
   modificationTimeMS (rhs.modificationTimeMS),
   absPath (rhs.absPath)
{
   // empty
}

// Destructor
SxFileInfo::Cache::~Cache ()
{
   /* empty */
}

// --- Methods
const bool & SxFileInfo::Cache::isDirty () const
{
   return dirty;
}

SxString const & SxFileInfo::Cache::getAbsPath () const
{
   SX_DBG_TRACE("");
   // Returning the absolute path
   return absPath;
}

const bool & SxFileInfo::Cache::isExisting ()
{
   if (isDirty ()) {
      update ();
   }
   return existing;
}

const int & SxFileInfo::Cache::getStatus ()
{
   if (isDirty ()) {
      update ();
   }
   return status;
}

const SxFileInfo::Cache::Type & SxFileInfo::Cache::getType ()
{
   if (isDirty ()) {
      update ();
   }
   return type;
}

const sxuid_t & SxFileInfo::Cache::getUID ()
{
   if (isDirty ()) {
      update ();
   }
   return userID;
}

const sxgid_t & SxFileInfo::Cache::getGID ()
{
   if (isDirty ()) {
      update ();
   }
   return groupID;
}

const int & SxFileInfo::Cache::getMode ()
{
   if (isDirty ()) {
      update ();
   }
   return mode;
}

const int64_t &SxFileInfo::Cache::getSize ()
{
   if (isDirty ()) {
      update ();
   }
   return size;
}

const unsigned long & SxFileInfo::Cache::getAccessTime ()
{
   if (isDirty ()) {
      update ();
   }
   return accessTime;
}
const unsigned long long & SxFileInfo::Cache::getAccessTimeMS ()
{
   if (isDirty ()) {
      update ();
   }
   return accessTimeMS;
}

const unsigned long & SxFileInfo::Cache::getModificationTime ()
{
   if (isDirty ()) {
      update ();
   }
   return modificationTime;
}
const unsigned long long & SxFileInfo::Cache::getModificationTimeMS ()
{
   if (isDirty ()) {
      update ();
   }
   return modificationTimeMS;
}

int SxFileInfo::Cache::getStatusInfo (SxStructStat *buffer) const
{
   SX_DBG_TRACE("");
#ifdef WIN32
   // Resolving shortcuts contained in the path before a stat command is
   // performed
   SxString str (SxSymLink::resolveShortcutsInPath (getAbsPath ()));
   if (str.isUnicode ()) return sxwstat ((wchar_t *)str.utf16 ().elements,
                                         buffer);
   else                  return sxstat (str.ascii (), buffer);
#else
   return sxstat (getAbsPath ().getElems (), buffer);
#endif /* WIN32 */
}

int SxFileInfo::Cache::getStatusInfoCompletelyResolved (SxStructStat *buffer)
   const
{
   SX_DBG_TRACE("");
#ifdef WIN32
   // Resolving shortcuts contained in the path before a stat command is
   // performed
   SxString str = SxSymLink::resolveShortcutsInAbsPath (absPath);
   if (str.isUnicode ()) return sxwstat ((wchar_t *)str.utf16 ().elements,
                                         buffer);
   else                  return sxstat (str.getElems (), buffer);
#else
   return sxstat (absPath.getElems (), buffer);
#endif /* WIN32 */
}

void SxFileInfo::Cache::setDirty ()
{
   dirty = true;
}

void SxFileInfo::Cache::setAbsPath (const SxString &absPath_)
{
   absPath = SxFileInfo::resolvePath (absPath_);
}

void SxFileInfo::Cache::update ()
{
   SX_DBG_TRACE("");
   // Resetting the type
   setType (Nothing);
   // --- Fetching the stat information on the file system node
   SxStructStat statBuf;
   int error = getStatusInfoCompletelyResolved (&statBuf);
   if (error)  {
      setGID ((sxgid_t)-1);
      setUID ((sxuid_t)-1);
      setMode (0);
      setSize (0);
      setAccessTime (0);
      setAccessTimeMS (0);
      setModificationTime (0);
      setModificationTimeMS (0);

   } else  {
      setGID (statBuf.st_gid);
      setUID (statBuf.st_uid);
      setMode (statBuf.st_mode);
      setSize (statBuf.st_size);

      setAccessTime ((unsigned long) (statBuf.st_atime));
      setModificationTime ((unsigned long) (statBuf.st_mtime));
#     ifdef LINUX
         setAccessTimeMS ((unsigned long long) (statBuf.st_atime * 1000
                          + statBuf.st_atim.tv_nsec / 1000000));
         setModificationTimeMS ((unsigned long long) (statBuf.st_mtime * 1000
                               + statBuf.st_mtim.tv_nsec / 1000000));
#     else
         setAccessTimeMS ((unsigned long long) (statBuf.st_atime * 1000));
         setModificationTimeMS ((unsigned long long) (statBuf.st_mtime * 1000));
#     endif


   }
   setStatus (error);

   // Fetching the pieces of stat-information that may not be resolved
   // completely (shortcuts are not resolved automatically like symbolic links
   // so this has to be done manually to a certain extent!)
   SxStructStat buffer;
   int err = getStatusInfo (&buffer);
   int modePartiallyResolved;
   if (err)  {
      modePartiallyResolved = 0;
   } else  {
      modePartiallyResolved = buffer.st_mode;
   }
   setExisting (err == 0);
   setTypeAttrFile (modePartiallyResolved);
   setTypeAttrDir (modePartiallyResolved);
   setTypeAttrSymLink ();

#  ifdef SX_FS_CACHING_ENABLED
      dirty = false;
#  endif /* SX_FS_CACHING_ENABLED */

   //std::cout << "absPath: \"" << absPath << "\" dirty: " << ((int)dirty) << " type: " << ((int)type)<< " existing: " << ((int)existing) << " groupID: " << groupID<< " userID: " << userID << " mode: " << mode << " size: " << size << " accessTime: " << accessTime << " modificationTime: " << modificationTime << std::endl;//TEST

}

void SxFileInfo::Cache::setStatus (int status_)
{
   status = status_;
}

void SxFileInfo::Cache::setTypeAttrFile (int modePartiallyResolved)
{
   // Ensure that this function is only called from within update ()
   setType ((Type)(type | ((S_ISREG (modePartiallyResolved))?
                           File : Nothing)));
}

void SxFileInfo::Cache::setTypeAttrDir (int modePartiallyResolved)
{
   // Ensure that this function is only called from within update ()
   setType ((Type)(type | ((S_ISDIR (modePartiallyResolved))?
                          Dir : Nothing)));
}

void SxFileInfo::Cache::setTypeAttrSymLink ()
{
   // Ensure that this function is only called from within update ()
#ifdef WIN32
   //setType ((Type)(type |
   //         (sxnumlibs_isShortcut (SxSymLink::resolveShortcutsInPath (getAbsPath ()).getElems ())?
   //          SymLink: Nothing)));

#else
   struct stat lstatBuf;
   // Here it is crucial to use lstat () instead of stat, because otherwise
   // there will never be found a symbolic link!
   int lstatErr = lstat (getAbsPath ().getElems (), &lstatBuf);
   if (!lstatErr) {
      // Testing if the element defined by the absolute path is a symbolic link
      setType ((Type)(type | ((S_ISLNK (lstatBuf.st_mode))? SymLink: Nothing)));
   }
#endif
}

void SxFileInfo::Cache::setType (Type type_)
{
   type = type_;
}

void SxFileInfo::Cache::setExisting (bool existing_)
{
   existing = existing_;
}

void SxFileInfo::Cache::setUID (sxuid_t userID_)
{
   userID = userID_;
}

void SxFileInfo::Cache::setGID (sxgid_t groupID_)
{
   groupID = groupID_;
}

void SxFileInfo::Cache::setMode (int mode_)
{
   mode = mode_;
}

void SxFileInfo::Cache::setSize (int64_t size_)
{
   size = size_;
}

void SxFileInfo::Cache::setAccessTime (unsigned long accessTime_)
{
   accessTime = accessTime_;
}
void SxFileInfo::Cache::setAccessTimeMS (unsigned long long accessTimeMS_)
{
   accessTimeMS = accessTimeMS_;
}

void SxFileInfo::Cache::setModificationTime (unsigned long modificationTime_)
{
   modificationTime = modificationTime_;
}

void SxFileInfo::Cache::setModificationTimeMS (unsigned long long modificationTimeMS_)
{
   modificationTimeMS = modificationTimeMS_;
}


// --------------------------------------------------------------------------

SxFileInfo operator/ (const SxFileInfo &a, const SxString &b)
{
   SX_DBG_TRACE("");
   // Checking input arguments
   SX_CHECK (b != SxString ());
   SX_CHECK (a.getAbsPath ().getSize() > 0, a.getAbsPath ().getSize());
   SX_CHECK (SxFileInfo::isAbs (a.getAbsPath ())); // only storing absolute
                                                   // paths
   // Returning the concatenated absolute path and string
   return SxFileInfo ( a.getAbsPath()
                     + '/'
                     + b);
}

SxFileInfo operator/ (const SxString &a, const SxFileInfo &b)
{
   SX_DBG_TRACE("");
   // Checking input arguments
   SX_CHECK (a != SxString ());
   SX_CHECK (b.getAbsPath ().getSize() > 0, b.getAbsPath ().getSize());
   SX_CHECK (SxFileInfo::isAbs (b.getAbsPath ())); // only storing absolute
                                                   // paths
   // Returning the concatenated string and absolute path
   return SxFileInfo(a + b.getAbsPath());
}

SxFileInfo operator+ (const SxFileInfo &a, const SxString &b)
{
   SX_DBG_TRACE("");
   // Checking input arguments
   SX_CHECK (b != SxString ());
   SX_CHECK (a.getAbsPath ().getSize() > 0, a.getAbsPath ().getSize());
   SX_CHECK (SxFileInfo::isAbs (a.getAbsPath ())); // only storing absolute
                                                   // paths
   // Returning the concatenated absolute path and string
   return SxFileInfo (a.getAbsPath() + b);
}

SxFileInfo operator+ (const SxString &a, const SxFileInfo &b)
{
   SX_DBG_TRACE("");
   // Checking input arguments
   SX_CHECK (a != SxString ());
   SX_CHECK (b.getAbsPath ().getSize() > 0, b.getAbsPath ().getSize());
   SX_CHECK (SxFileInfo::isAbs (b.getAbsPath ())); // only storing absolute
                                                   // paths
   // Returning the concatenated string and absolute path
   return SxFileInfo(a + b.getAbsPath());
}

// --------------------------------------------------------------------------

std::ostream &operator<< (std::ostream &s, const SxFileInfo &in)
{
   SX_DBG_TRACE("");
   s << in.getAbsPath ();
   return s;
}


#ifdef WIN32
   std::wostream& operator<< (std::wostream &s, const SxFileInfo &in)
   {
      s << in.getAbsPath ();
      return s;
   }
#endif

//------------------------------------------------------------------------------
// SxFISortedByTime
//------------------------------------------------------------------------------
SxFISortedByTime::SxFISortedByTime (const SxString &path_)
  : SxFileInfo (path_)
{
   // empty
}

SxFISortedByTime::SxFISortedByTime (const SxFileInfo &rhs_)
  : SxFileInfo (rhs_)
{
   // empty
}

SxFISortedByTime::~SxFISortedByTime ()
{
   // empty
}

bool SxFISortedByTime::operator== (const SxFISortedByTime &rhs) const
{
   return (SxFileInfo::compareModificationTime (*this, rhs) == 0);
}

bool SxFISortedByTime::operator> (const SxFISortedByTime &rhs) const
{
   return (SxFileInfo::compareModificationTime (*this, rhs) == 1);
}

bool SxFISortedByTime::operator< (const SxFISortedByTime &rhs) const
{
   return (SxFileInfo::compareModificationTime (*this, rhs) == -1);
}

//------------------------------------------------------------------------------
// SxFISortedBySize
//------------------------------------------------------------------------------
SxFISortedBySize::SxFISortedBySize (const SxString &path_)
  : SxFileInfo (path_)
{
   // empty
}

SxFISortedBySize::SxFISortedBySize (const SxFileInfo &rhs_)
  : SxFileInfo (rhs_)
{
   // empty
}

SxFISortedBySize::~SxFISortedBySize ()
{
   // empty
}

bool SxFISortedBySize::operator== (const SxFISortedBySize &rhs) const
{
   return (SxFileInfo::compareSize (*this, rhs) == 0);
}

bool SxFISortedBySize::operator> (const SxFISortedBySize &rhs) const
{
   return (SxFileInfo::compareSize (*this, rhs) == 1);
}

bool SxFISortedBySize::operator< (const SxFISortedBySize &rhs) const
{
   return (SxFileInfo::compareSize (*this, rhs) == -1);
}

int sxstat (const char *path, SxStructStat *buf)
{
#ifdef WIN32
   return _stat64 (path, buf);
#else
   return stat (path, buf);
#endif
}

#ifdef WIN32
   int sxwstat (const wchar_t *path, SxStructStat *buf)
   {
      return _wstat64 (path, buf);
   }
#endif


int sxfstat (int fd, SxStructStat *buf)
{
#ifdef WIN32
   return _fstat64 (fd, buf);
#else
   return fstat (fd, buf);
#endif
}
