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
#include <SxFSExploreAction.h>
#include <SxFSAction.h>
#include <SxFSError.h>
#include <SxSortedList.h>
#include <SxRegex.h>

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

// --- Constructors
// Constructor
SxFSExploreAction::SxFSExploreAction ()
{
   SX_DBG_TRACE("");
   // empty
}
// Destructor
SxFSExploreAction::~SxFSExploreAction ()
{
   SX_DBG_TRACE("");
   // empty
}

SxDir SxFSExploreAction::getHome ()
{
   SX_DBG_TRACE("");
   return SxDir (SxFileInfo::getHomeStr ());
}

SxDir SxFSExploreAction::getTmp ()
{
   SX_DBG_TRACE("");
   return SxDir (SxFileInfo::getTmpStr ());
}

SxDir SxFSExploreAction::pwd ()
{
   SX_DBG_TRACE("");

#  ifdef WIN32
      // Attempting to retrieve the current working directory
      SxArray<uint16_t> buffer(SxFileInfo::MAX_PATH_LEN);
      if (!_wgetcwd ((LPWSTR)buffer.elements, SxFileInfo::MAX_PATH_LEN))  {
         // possible reasons for a failure are a too small path length, a lack
         // of search or read privileges, or a memory allocation problem
         SX_THROW ("An error occurred during an attempt to "
                   "read the current working directory. "
                   + SxFSError::getGetcwdErrMsg ());
      }
      return SxDir (SxString::fromUtf16 (buffer.elements));
#  else
      // Attempting to retrieve the current working directory
      SxArray<char> buffer(SxFileInfo::MAX_PATH_LEN);
      if (!getcwd (buffer.elements, SxFileInfo::MAX_PATH_LEN))  {
         // possible reasons for a failure are a too small path length, a lack
         // of search or read privileges, or a memory allocation problem
         SX_THROW ("An error occurred during an attempt to "
                   "read the current working directory. "
                   + SxFSError::getGetcwdErrMsg ());
      }
      return SxDir (SxString::unicodeFromUtf8 (buffer.elements));
#  endif
}

SxList<SxFISortedByTime> SxFSExploreAction::ls_t (const SxFileInfo &path)
{
   SX_DBG_TRACE("");
   SxList<SxFileInfo> results = SxFSExploreAction::ls (path);
   SxList<SxFileInfo>::Iterator itResults;
   SxSortedList<SxFISortedByTime> sortedList;
   for (itResults = results.begin ();
        itResults != results.end ();
        ++itResults)
   {
      sortedList << SxFISortedByTime (*itResults);
   }
   return ((SxList<SxFISortedByTime>) sortedList);
}

SxList<SxFISortedBySize> SxFSExploreAction::ls_S (const SxFileInfo &path)
{
   SX_DBG_TRACE("");
   SxList<SxFileInfo> results = SxFSExploreAction::ls (path);
   SxList<SxFileInfo>::Iterator itResults;
   SxSortedList<SxFISortedBySize> sortedList;
   for (itResults = results.begin ();
        itResults != results.end ();
        ++itResults)
   {
      sortedList << SxFISortedBySize (*itResults);
   }
   return ((SxList<SxFISortedBySize>) sortedList);
}

bool SxFSExploreAction::contains (const SxList<SxFileInfo> &visited,
                                  const SxFileInfo &current)
{
   bool ret (false);
   SxList<SxFileInfo>::ConstIterator it;
   for (it = visited.begin (); !ret && it != visited.end (); ++it)  {
      if (SxFileInfo::equals ((*it), current))  {
         ret = true;
      }
   }
   return ret;
}

void SxFSExploreAction::findName (const SxDir        &curDir,
                                  const SxString     &pattern,
                                  SxList<SxFileInfo> *visited,
                                  SxList<SxFileInfo> *ret,
                                  bool                followLinks)
{
   //std::cout << "Looking for pattern \"" << pattern << "\" in \"";//TEST
   //std::cout << curDir.getAbsPath () << "\"" << std::endl;//TEST

   SX_CHECK (ret);
   SxList<SxDir> childDirs;
   // Getting the entries of the current directory
   (*ret) << SxFSExploreAction::searchNode (curDir, pattern, NULL, &childDirs);
   //(*ret) << SxFSExploreAction::ls (curDir / target, curDir, &childDirs);
   SX_CHECK (visited);
   // Storing that the current directory was visited
   (*visited) << curDir;
   // Traversing the children of the current directory
   SxList<SxDir>::Iterator itChildDirs;
   for (itChildDirs = childDirs.begin ();
        itChildDirs != childDirs.end ();
        ++itChildDirs)
   {
      // Checking if the current element is a directory that still has to be
      // processed
      SxDir const &childDir = (*itChildDirs);
      if (((!followLinks && !childDir.isSymLink ()) ||
           (followLinks && !SxFSExploreAction::contains ((*visited), childDir)))
          && (childDir.getName () != "." && childDir.getName () != ".."))  {
         SxFSExploreAction::findName (childDir, pattern, visited, ret,
                                      followLinks);
      }
   }
}

SxList<SxFileInfo> SxFSExploreAction::find (const SxString &target,
                                            const SxDir &curDir)
{
   SxList<SxFileInfo> ret, visited;
   SxString pat = target;
   prepareRegExpPattern (&pat);
   finishRegExpPattern (&pat);
   SxFSExploreAction::findName (curDir, pat, &visited, &ret);
   return ret;
}

void SxFSExploreAction::prepareRegExpPattern (SxString *pStr)
{
   SX_CHECK (pStr);
   SxString &str = *pStr;
   //std::cout << "passed string: \"" << str << "\""<< std::endl;//TEST
   // Escaping regular expression characters
   // Applying a trick to avoid that the user can specify the
   // regular expression characters '(' and ')', which have to be deleted
   // in the very first step, without screwing up windows' paths
#ifdef WIN32
   str = str.substitute ("\\","/");
#endif /* WIN32 */
   str = str.substitute ("(","\\(");
   str = str.substitute (")","\\)");
   str = str.substitute ("|","\\|");
   // Removing empty {*}-sequence
   str = str.substitute ("{}","");
   // Searching and substituting {*,*,...}-sequences by (*|*|...)-sequences
   SxString pat = "^.*(\\{(.+)\\}).*$";
   SxList<SxString> elements = SxRegex(pat, "P").match (str);
   //std::cout << elements << std::endl;//TEST
   while (elements.getSize () == 4)  {
      SxString values = elements(2).substitute (",", "|");
      str = str.substitute(elements(1), "(" +values + ")");
      elements = SxRegex(pat, "P").match (str);
      //std::cout << elements << std::endl;//TEST
   }
   //std::cout << "prepared string: \"" << str << "\"" << std::endl;//TEST
}

void SxFSExploreAction::finishRegExpPattern (SxString *pStr)
{
   SX_CHECK (pStr);
   SxString &str = *pStr;
   str = SxString("^")
                    + str.substitute ("[","\\[")
                    .substitute ("]","\\]")
                    .substitute ("^","\\^")
                    .substitute ("$","\\$")
                    .substitute ("+","\\+")
                    .substitute (".", "[.]")
                    .substitute ("*", ".*")
                    .substitute ("?", ".")
                    + "$";
}

SxList<SxFileInfo> SxFSExploreAction::ls (const SxFileInfo &src)
{
   SX_DBG_TRACE("");
   SxString const &srcStr = src.getAbsPath ();
   SxList<SxFileInfo>  fileInfos;

   SxList<SxDir> *subDirs = NULL;
   SxString limitStr;
#ifdef WIN32
   limitStr = srcStr.head (3); // 'C:\\'
#else
   limitStr = "/";
#endif /* WIN32 */
   SxFileInfo limit(limitStr);

   // --- Replacing some special tokens in the passed pattern string by their
   //     regular expression counterparts
   SxList<SxString> preparedNodes;
   SxString str = srcStr.right (limit.getAbsPath ());
   prepareRegExpPattern (&str);
   // Again the trick of separator substitution to circumvent problems with
   // regular expressions

   SxList<SxString> nodes = str.tokenize ('/');
   SxList<SxString>::Iterator nodesIt;
   SxList<SxDir> foundDirs;
   SxList<SxDir> newFoundDirs;
   SxString curPreparedNodeStr;
   SxList<SxDir>::Iterator itFoundDirs;
   SxDir curDir;
   for (nodesIt = nodes.begin (); nodesIt != nodes.end (); ++nodesIt)  {
      fileInfos = SxList<SxFileInfo> ();
      curPreparedNodeStr = *nodesIt;
      finishRegExpPattern (&curPreparedNodeStr);
      preparedNodes << curPreparedNodeStr;

      //std::cout << "curPreparedNodeStr = " << curPreparedNodeStr;//TEST
      //std::cout << std::endl;//TEST

      newFoundDirs = SxList<SxDir> ();
      if (nodesIt == nodes.begin())  {
         curDir = limit;
         fileInfos << searchNode (curDir, curPreparedNodeStr, &foundDirs,
                                  subDirs);
      } else  {
         for (itFoundDirs = foundDirs.begin ();
              itFoundDirs != foundDirs.end ();
              ++itFoundDirs)
         {
            curDir = (*itFoundDirs);
            fileInfos << searchNode (curDir, curPreparedNodeStr, &newFoundDirs);
         }
         foundDirs = newFoundDirs;
      }
   }
   // Returning the found pieces of information on the directory
   return fileInfos;
}

// Collect all filenames in given directory. Recursive.
void SxFSExploreAction::listDir (const SxString   path_,
                                 SxList<SxString> *fileList_)
{
   SX_CHECK (fileList_);

#ifdef WIN32
      SxList<SxFileInfo> items = SxFSAction::ls (path_ + "/*");
      for (SxList<SxFileInfo>::ConstIterator it  = items.begin ();
                                             it != items.end ();
                                           ++it)  {
         if        (SxFSAction::test_f (it->getAbsPath ())) {
            // Append valid file to the list
            fileList_->append (it->getAbsPath ());
         } else if (SxFSAction::test_d (it->getAbsPath ())) {
            // Append the folder to the list
            fileList_->append (it->getAbsPath ());
            listDir (it->getAbsPath (), fileList_);
         }
      }
#else
      struct dirent *entry = NULL;
      DIR *dir = NULL;
      struct stat statBuf;
      SxList<SxString> dirList;
      SxList<SxString>::Iterator it;
      SxString name;

      // Is the current iteration root folder.
      bool rootDir = (fileList_->getSize () == 0);

      dir = opendir (path_.getElems ());
      if (!dir)  {
         SX_THROW ("Can't open directory '" + path_ + "': " + sxstrerror());
      }

      // --- directory name on top
      if (path_ != ".")  {
         fileList_->append (path_);
      }

      // --- current level
      while ((entry = readdir (dir)) != NULL)  {
         if (strcmp (entry->d_name, ".") != 0
          && strcmp (entry->d_name, "..") != 0)
         {
             if (path_ != ".")  {

               if (path_.findLast ('/') != (path_.getSize () - 1))
                  name = path_ + SxString("/");
               else name = path_;

               name += SxString::unicodeFromUtf8 (entry->d_name);

             }  else  {
                name = entry->d_name;
             }

             if (entry->d_type == DT_DIR)  {
                dirList.append (name);
             }  else if (entry->d_type == DT_UNKNOWN) {
                if (lstat (name.getElems (), &statBuf) == 0)  {
                  if (S_ISDIR(statBuf.st_mode))  {
                     dirList.append (name);
                  }  else  {
                     fileList_->append (name);
                  }
                }
             }  else  {
                fileList_->append (name);
             }
         }
      }

      closedir (dir);

      // --- sub-directory level
      for (it = dirList.begin (); it != dirList.end (); ++it)  {
         SxFSExploreAction::listDir (*it, fileList_);
      }

      // After the complete of operation remove the root search folder
      if (rootDir) {
         fileList_->removeFirst ();
      }

#endif /* WIN32 */
}


SxList<SxFileInfo> SxFSExploreAction::searchNode (const SxFileInfo &curDir,
                                                  const SxString &pattern,
                                                  SxList<SxDir> *foundDirs,
                                                  SxList<SxDir> *subDirs)
{
   SxList<SxFileInfo>  fileInfos;
   SxString const & curDirStr = curDir.getAbsPath ();
   //std::cout << "SEARCHING... \"" << pattern << "\" in \"" << curDirStr;//TEST
   //std::cout << "\"" << std::endl;//TEST
#ifdef WIN32
   // --- Trying to open the directory and to fetch a file handle to the first
   //     file of the directory
   WIN32_FIND_DATAW fileData;
   HANDLE dir = FindFirstFileW ((LPCWSTR)(curDirStr + "\\*").utf16 ().elements,
                                &fileData);
   if (dir == INVALID_HANDLE_VALUE)  {
      if (GetLastError () == ERROR_FILE_NOT_FOUND)  return fileInfos;
      SX_THROW ("Can't open directory '" + curDirStr
               + "'. " + SxFSError::getOpendirErrMsg ());
   }

   // --- Iterating over all items in the directory
   bool done = false;
   while (!done)  {
      SxString name
         = SxString::fromUtf16 ((const uint16_t *)fileData.cFileName);
      fileInfos << processNode (name, curDir, pattern, foundDirs, subDirs);
      if (!FindNextFileW (dir, &fileData))
      {
         if (GetLastError () == ERROR_NO_MORE_FILES)  {
            done = true;
         } else  {
            SX_THROW ("Can't open directory '" + curDirStr
                        +"'. "+SxFSError::getOpendirErrMsg());
         }
      }
   }
   FindClose (dir);
#else
   // --- Trying to open the directory
   DIR *dir = opendir (curDirStr.getElems ());
   if (!dir)  {
      SX_THROW ("Can't open directory '" + curDirStr
                +"'. "+ SxFSError::getOpendirErrMsg ());
   }

   // Initializing the directory entry
   struct dirent *entry = NULL;

   // --- Traversing the directory entries
   while ( (entry = readdir (dir)) != NULL )  {
      fileInfos << processNode (entry->d_name, curDir, pattern, foundDirs,
                               subDirs);
   }
   // --- Closing the opened directory
   closedir (dir);
#endif
   return fileInfos;
}


SxList<SxFileInfo> SxFSExploreAction::processNode (const SxString   &name,
                                                   const SxFileInfo &curDir,
                                                   const SxString   &pattern,
                                                   SxList<SxDir>    *foundDirs,
                                                   SxList<SxDir>    *subDirs)
{
   //SX_DBG_MSG (name << " " << curDir);
   SxList<SxFileInfo> fileInfos;
   SxFileInfo info;

   SxString absPath = curDir.getAbsPath ();
#  ifdef WIN32
      if (absPath.getSize () == 3) info = SxFileInfo (curDir + name); // 'C:\\'
      else                         info = SxFileInfo (curDir / name);
#  else
      if (absPath == "/") info = SxFileInfo (curDir + name);
      else                info = SxFileInfo (curDir / name);
#  endif /* WIN32 */

   // Treating "." and ".." in a special way (otherwise */somefile would
   // also search in the current and the parent directory)
   if ((SxString (".") != name  &&
        SxString ("..") != name))
   {
      // Performing pattern matching using regular expressions
      if (SxRegex(pattern, "P").match (info.getName ()).getSize () > 0)  {
         fileInfos << info;
         if (foundDirs)  {
            SxDir infoDir (info);
            if (infoDir.exists ())  {
               //std::cout << "FOUND! " << infoDir.getAbsPath ();//TEST
               //std::cout << std::endl;//TEST
               (*foundDirs) << infoDir;
            }
         }
      }
      // Storing all found sub directories apart from "." and ".." provided
      // this is desired
      if (subDirs)  {
         SxDir subDir (info);
         if (subDir.exists ())  {
            (*subDirs) << subDir;
         }
      }
   } else if ((SxString (".") == name && pattern == "^[.]$") ||
              (SxString ("..") == name && pattern == "^[.][.]$"))
   {
      fileInfos << info;
      if (foundDirs)  {
         SxDir infoDir (info);
         if (infoDir.exists ())  {
            //std::cout << "FOUND! " << infoDir.getAbsPath ();//TEST
            //std::cout << std::endl;//TEST
            (*foundDirs) << infoDir;
         }
      }
   }
   return fileInfos;
}

SxList<SxFileInfo> SxFSExploreAction::getFileInfos (const SxFileInfo &path)
{
   SX_DBG_TRACE("");
   SxList<SxFileInfo> infos;
   SxList<SxFileInfo> children = SxFSExploreAction::ls (path / "*");
   SxList<SxFileInfo>::Iterator itChildren;
   for (itChildren = children.begin ();
        itChildren != children.end ();
        ++itChildren)
   {
      const SxFileInfo & curChild = *itChildren;
      //SX_DBG_MSG ("name " << curChild.getName ());
      if (curChild.getName () != "." && curChild.getName () != "..")  {
         infos << curChild;
      }
   }
   return infos;
}


SxList<SxFile> SxFSExploreAction::getFiles (const SxFileInfo &path)
{
   SX_DBG_TRACE("");
   SxList<SxFile> files;
   SxList<SxFileInfo> children = SxFSExploreAction::getFileInfos (path);
   SxList<SxFileInfo>::Iterator itChildren;
   for (itChildren = children.begin ();
        itChildren != children.end ();
        ++itChildren)
   {
      const SxFileInfo & curChild = *itChildren;
      if (curChild.isFile ())  {
         files << curChild;
      }
   }
   return files;
}

SxList<SxDir> SxFSExploreAction::getDirs (const SxFileInfo &path)
{
   SX_DBG_TRACE("");
   SxList<SxDir> directories;
   SxList<SxFileInfo> children = SxFSExploreAction::getFileInfos (path);
   SxList<SxFileInfo>::Iterator itChildren;
   for (itChildren = children.begin ();
        itChildren != children.end ();
        ++itChildren)
   {
      const SxFileInfo & curChild = *itChildren;
      if (curChild.isDir ())  {
         directories << SxDir (curChild.getAbsPath ());
      }
   }
   return directories;
}

SxList<SxSymLink> SxFSExploreAction::getSymLinks (const SxFileInfo &path)
{
   SX_DBG_TRACE("");
   SxList<SxSymLink> symLinks;
   SxList<SxFileInfo> children = SxFSExploreAction::getFileInfos (path);
   SxList<SxFileInfo>::Iterator itChildren;
   for (itChildren = children.begin ();
        itChildren != children.end ();
        ++itChildren)
   {
      const SxFileInfo & curChild = *itChildren;
      if (curChild.isSymLink ())  {
         symLinks << SxSymLink (curChild.getAbsPath ());
      }
   }
   return symLinks;
}
