// --- Including header-files
#include <SxFSAuthorAction.h>
#include <SxFSError.h>
#include <SxRegex.h>
#include <SxException.h>
#include <fcntl.h>
#include <errno.h>
#ifdef WIN32
#  include <shortcut.h>
#  include <windows.h>
#  include <io.h>
#  include <unknwn.h>
#  include <shlobj.h>
#  include <objbase.h>
#  include <olectl.h>
#  include <initguid.h>
#else
#  include <sys/types.h>
#  include <sys/stat.h>
#  include <unistd.h>
#endif /* WIN32 */

void SxFSAuthorAction::chmod (int perm, const SxFileInfo& target)
{
   if (perm < 0 || perm > 777)  {
      SX_THROW ("Can't change permissions of '" + target.getAbsPath ()
                + "' to mode " + SxString(perm) + ". "
                "The permission value is out of range '000' to '777'.");

   }
   SxFSAuthorAction::chmod (target, perm);
}

void SxFSAuthorAction::chmod (const SxString   &modStr,
                              const SxFileInfo &target)
{
   if (modStr.getSize () < 1)  {
      SX_THROW ("Can't change permissions of '"
                + target.getAbsPath () + "' to mode '" + modStr
                + "'. The permission string is empty. "
                  "Enter a string in the same format as this example: "
                  "rwxr--r--");
   }

   SxList<SxString> ls;
   SxList<SxString>::Iterator it;

   if (   modStr(0) == 'u'
       || modStr(0) == 'g'
       || modStr(0) == 'o'
       || modStr(0) == 'a')
   {
      // --- u+r
      try  {
         const SxString expr = "^([ugoa]+)(\\+|-|=)([rwx]+)$";
         ls  = SxRegex (expr, "P").match (modStr);
      } catch (SxException e)  {
         // --- invalid ls result
      }

      if (ls.getSize () != 5)  {
         SX_THROW ("Can't change permissions of '"
                   + target.getAbsPath () + "' to mode '" + modStr
                   + "'. The permission string has wrong syntax. "
                     "Enter a string in the same format as this example: "
                     "u+x");
      }

      SxString const &dom = ls(1);
      SxString const &op = ls(2);
      SxString const &priv = ls(3);
      int bits = 0;
      for (int i = 0; i < priv.getSize (); ++i)  {
         if (priv(i) == 'r')  {
            bits |= 0x04;
         } else if (priv(i) == 'w')  {
            bits |= 0x02;
         } else if (priv(i) == 'x')  {
            bits |= 0x01;
         }
      }
      int userBits (0);
      int groupBits (0);
      int otherBits (0);
      for (int i = 0; i < dom.getSize (); ++i)  {
         if (dom(i) == 'u')  {
            userBits = (bits << 6);
         } else if (dom(i) == 'g')  {
            groupBits = (bits << 3);
         } else if (dom(i) == 'o')  {
            otherBits = bits;
         } else if (dom(i) == 'a')  {
            userBits = (bits << 6);
            groupBits = (bits << 3);
            otherBits = bits;
         }
      }
      bits = userBits | groupBits | otherBits;
      if (op == SxString ("+"))  {
         chmod (target, (target.getMode () | bits));
      } else if (op == SxString ("-"))  {
         chmod (target, (target.getMode () & (~bits)));
      } else  {
         chmod (target, bits);
      }
   }  else  {
      // --- rwxrwxrwx
      try  {
         const SxString triple = "((r|-)(w|-)(x|-))";
         const SxString expr = "^(" + triple + triple + triple + ")$";
         ls  = SxRegex(expr, "P").match (modStr);
      } catch (SxException e)  {
         // --- invalid ls result
      }

      if (ls.getSize () != 15)  {
         SX_THROW ("Can't change permissions of '"
                   + target.getAbsPath () + "' to mode '" + modStr
                   + "'. The permission string has wrong syntax. "
                     "Enter a string in the same format as this example: "
                     "rwxr--r--");
      }

      SxString const &user = ls(2);
      SxString const &group = ls(6);
      SxString const &other = ls(10);
      SX_CHECK (user.getSize () == group.getSize () &&
                group.getSize () == other.getSize () &&
                other.getSize () == 3);
      int userBits = 0;
      int groupBits = 0;
      int otherBits = 0;
      for (int i = 0; i < user.getSize (); ++i)  {
         if (user(i)!= '-')  {
            userBits |= 1 << (2-i);
         }
      }
      userBits <<= 6;
      for (int i = 0; i < group.getSize (); ++i)  {
         if (group(i)!= '-')  {
            groupBits |= 1 << (2-i);
         }
      }
      groupBits <<= 3;
      for (int i = 0; i < other.getSize (); ++i)  {
         if (other(i)!= '-')  {
            otherBits |= 1 << (2-i);
         }
      }
      SxFSAuthorAction::chmod (target, userBits | groupBits | otherBits);
   }
}

void SxFSAuthorAction::chmod (const SxFileInfo& target,
                              SxFileInfo::AccessGroup group,
                              SxFileInfo::AccessType type)
{
   SxString const & targetStr (target.getAbsPath());
   if (!target.exists ())  {
      SX_THROW ("Can't change permissions of '" + targetStr
                + "'. The filename does not exist.");
   }

   int mask = 0;
   if (group == SxFileInfo::Owner)  {
      mask |= ((type & SxFileInfo::Readable)? S_IRUSR : 0)
            | ((type & SxFileInfo::Writable)? S_IWUSR : 0)
            | ((type & SxFileInfo::Executable)? S_IXUSR : 0);
   }
#ifndef WIN32
   if (group == SxFileInfo::Group)  {
      mask |= ((type & SxFileInfo::Readable)   ? S_IRGRP : 0)
                       | ((type & SxFileInfo::Writable)? S_IWGRP : 0)
                       | ((type & SxFileInfo::Executable)? S_IXGRP : 0);
   }
   if (group == SxFileInfo::Other)  {
      mask |= ((type & SxFileInfo::Readable)   ? S_IROTH : 0)
                       | ((type & SxFileInfo::Writable)? S_IWOTH : 0)
                       | ((type & SxFileInfo::Executable)? S_IXOTH : 0);
   }
#endif /* not WIN32 */
   SxFSAuthorAction::chmod (target, mask);
}

void SxFSAuthorAction::chmod (const SxFileInfo& target, int mask)
{
   SxString const &targetStr (target.getAbsPath ());
#  ifdef WIN32
      if (targetStr.isUnicode ()) {

         if (  ::_wchmod ((const wchar_t *)targetStr.utf16 ().elements, mask)
            == -1)
         {
            SX_THROW ("Changing permissions of '" + targetStr + "' to mode "
                      + SxString(mask) + "failed. "
                      + SxFSError::getChmodErrMsg ());
         }
      } else {
         //TODO Test whether this substitute is .NET dependent?!
         if (::_chmod (targetStr.getElems (), mask) == -1)  {
            SX_THROW ("Changing permissions of '" + targetStr + "' to mode "
                      + SxString(mask) + ". " + SxFSError::getChmodErrMsg ());
         }
      }
#  else
      SxString path (targetStr.getElems ());
      SxString error;

      if (mask < 0 || mask > UINT16_MAX) {
         error = "Invalid mode value " + SxString(mask) + ". "
                 "Enter a file mode that fits to octal mask 0177777.";
      }  else  {
         mode_t mode = static_cast<mode_t>(mask);
         if (::chmod (path.ascii(), mode) != 0) {
            error = "chmod() failed: " + sxstrerror ();
         }
      }

      if (error != "")  {
         SX_THROW("Can't change permissions of '"
                  + path + "' to mode 0"
                  + SxString::sprintf("%06o", mask)
                  + ". " + error);
      }


#  endif /* WIN32 */
   target.setDirty ();
}

void SxFSAuthorAction::chgrp (const SxFileInfo& target, sxgid_t gid)
{
   SxString const & targetStr (target.getAbsPath());
   if (!target.exists ())  {
      SX_THROW ("Can't change group of '" + targetStr
                + "'. The filename does not exist.");
   }

   int error = -1;
#ifdef WIN32
   //std::cerr << "Unfortunately chgrp isn't implemented for windows yet.";//TEST
   //std::cerr << std::endl;//TEST
   //SX_EXIT;//TODO Fix (moved to trouble ticket)
   error = 0;
#else
   // Fetching the user id of the file, directory or symbolic link that should
   // get a new group id and performing the operation
   error = ::chown (targetStr.getElems(), target.getUID (), gid);
#endif /* WIN32 */
   if (error)  {
      SX_THROW ("Can't change group of '" + targetStr + "' to gid "
                + SxString(gid) + ". " + SxFSError::getChownErrMsg ());
   }
   target.setDirty ();
}

void SxFSAuthorAction::chown (const SxFileInfo& target, sxuid_t uid)
{
   SxString const & targetStr (target.getAbsPath());
   if (!target.exists ())  {
      SX_THROW ("Can't change owner of '" + targetStr
                + "'. The filename does not exist.");
   }

   int error = -1;
#ifdef WIN32
   //std::cerr << "Unfortunately chown isn't implemented for windows yet.";//TEST
   //std::cerr << std::endl;//TEST
   //SX_EXIT;//TODO Fix (moved to trouble ticket)
   error = 0;
#else
   // Fetching the group id of the file, directory or symbolic link that should
   // get a new user id and performing the operation
   error = ::chown (targetStr.getElems (), uid, target.getGID ());
#endif /* WIN32 */

   if (error)  {
      SX_THROW ("Can't change owner of '" + targetStr + "' to uid "
                + SxString(uid) + ". " + SxFSError::getChownErrMsg ());
   }
   target.setDirty ();
}
