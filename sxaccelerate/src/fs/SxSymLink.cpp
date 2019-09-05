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
#include <SxSymLink.h>
#include <SxFSError.h>
#ifdef WIN32
#  include <shortcut.h>
#  include <windows.h>
#else
#  include <unistd.h>
#endif /* WIN32 */

// --- Constructors
// Constructor
/** \brief Creates a SxSymLink-object */
SxSymLink::SxSymLink (const SxString &path_) : SxFileInfo (path_)
{
   SX_DBG_TRACE("");
   //empty
}

SxSymLink::SxSymLink (const SxFileInfo &rhs) : SxFileInfo (rhs)
{
   SX_DBG_TRACE("");
   // empty
}

// Destructor
/** \brief Destroys a SxSymLink-object */
SxSymLink::~SxSymLink ()
{
   SX_DBG_TRACE("");
   //empty
}

// --- Methods
SxString const SxSymLink::getTarget () const
{
   SX_DBG_TRACE("");
   SxString ret;
   SxString const &str = getAbsPath ();

#ifdef WIN32
   ret = SxSymLink::getShortcutTarget (str);
#else
   ssize_t error = -1;
   char buf[SxFileInfo::BUFFER_SIZE];

   // Trying to figure out the target of the symbolic link
   error = ::readlink (str.getElems (),
                       buf,
                       SxFileInfo::BUFFER_SIZE);
   // Testing if either the symbolic link could not be read
   if (error == -1)  {
      SX_THROW ("Can't get target of the symbolic link '"
               + str + "'. " + SxFSError::getReadlinkErrMsg ());
   } else  {
      ret = SxString (buf, error);
   }
#endif /* WIN32 */

   return ret;
}

bool SxSymLink::exists () const
{
   SX_DBG_TRACE("");
   // Checking whether the absolute path leads to an existing symbolic
   // link
   return isSymLink ();
}

bool SxSymLink::isValid () const
{
   SX_DBG_TRACE("");
   // Checking whether the target of the symbolic link exists
   return SxFileInfo::exists ();
}

#ifdef WIN32
SxString SxSymLink::resolveShortcutsInPath (const SxString &lnkSrc)
{
   SX_DBG_TRACE("");
   return resolveShortcuts (lnkSrc, true);
}

SxString SxSymLink::resolveShortcutsInAbsPath (const SxString &lnkSrc)
{
   SX_DBG_TRACE("");
   return resolveShortcuts (lnkSrc);
}

SxString SxSymLink::resolveShortcuts (const SxString &lnkSrc, bool keepName)
{
   SX_DBG_TRACE("");
   if (lnkSrc.find (".lnk") != -1)  {
      SxList<SxString> tokens = lnkSrc.tokenize ('/');
      SxString curPath;
      bool cont = true;
      SxList<SxString>::Iterator itTokens;
      int i;
      for (itTokens = tokens.begin (), i = 0;
           cont && itTokens != tokens.end ();
           ++itTokens, ++i)
      {
         if (itTokens == tokens.begin ())  {
            // Initializing the current path element with the root element
            curPath = (*itTokens);
         } else  {
            SxString prevPath = curPath;
            // --- Checking whether the current element is a shortcut or the
            //     last token/name that should be preserved
            SxString curElem ((SxFileInfo (curPath) /
                               (*itTokens)).getAbsPath ());
            if (!(keepName && i == tokens.getSize () - 1) &&
                sxnumlibs_isShortcut (curElem.getElems ()))
                 // Willingly avoiding call of "SxSymLink(curElem).exists ()"
                 // as this would lead to circular dependencies
            {
               // Substituting the current element by the shortcut target
               try  {
                  curPath = SxSymLink::getShortcutTarget (curElem);
               } catch (SxException ex)  {
                  curPath = prevPath + SxString ('/') + (*itTokens);
                  cont = false;
                  ex.print ();
               }
            } else  {
               // Appending the current path element
               curPath += SxString ('/') + (*itTokens);
            }
         }
      }
      return curPath;
   } else  {
      return lnkSrc;
   }
}

SxString SxSymLink::getShortcutTarget (const SxString &lnkSrc)
{
   SX_DBG_TRACE("");
   char lnkDest[MAX_PATH];
   char *elem0 = lnkDest;
   char **pLnkDest = &elem0;
   int err = ::sxnumlibs_getShortcutTarget (lnkSrc.getElems (), pLnkDest);
   SxFSError::throwShortcutException (err, lnkSrc);
   return SxString (lnkDest);
}

#endif /* WIN32 */
