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

#include <SxPluginLoader.h>
#include <SxTypeDefs.h> /* sxFnCast */
#include <stdio.h>
#include <stdlib.h>
#ifndef WIN32
#  include <dlfcn.h>
#endif /* WIN32*/

SxPluginLoader::SxPluginLoader () : soHandle (NULL)
{
   // empty
}


SxPluginLoader::~SxPluginLoader ()
{
   // empty
}


bool SxPluginLoader::openSharedObject (const SxString &soName)
{
#  ifdef WIN32
      soHandle = LoadLibrary (soName.ascii());
      if (!soHandle)  {
         sxprintf ("\nERROR: Can't open plugin %s\n", soName.ascii());
         int err = GetLastError();
         LPVOID lpMsgBuf;//, lpDisplayBuf;
         FormatMessage( FORMAT_MESSAGE_ALLOCATE_BUFFER 
                      | FORMAT_MESSAGE_FROM_SYSTEM
                      | FORMAT_MESSAGE_IGNORE_INSERTS,
                      NULL, err, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                      (LPTSTR)&lpMsgBuf, 0, NULL);
//       lpDisplayBuf = (LPVOID)LocalAlloc(LMEM_ZEROINIT,
//                        (lstrlen((LPCTSTR)lpMsgBuf)
//                        +lstrlen((LPCTSTR)lpszFunction)+40)*sizeof(TCHAR));
         sxprintf ("       Windows DLL: %d: %s", err, lpMsgBuf);
         LocalFree (lpMsgBuf);
         return false;
      }
#  else
      soHandle = dlopen (soName.ascii(), RTLD_NOW);
      if (!soHandle)  {
         sxprintf ("\nERROR while loading %s\n", soName.ascii());
         sxprintf ("%s\n", dlerror ());
         return false;
      }
      dlerror ();
#  endif /* WIN32 */
   return true;
}


void SxPluginLoader::closeSharedObject ()
{
#  ifdef WIN32
      FreeLibrary (soHandle);
#  else
      dlclose (soHandle);
#  endif /* WIN32 */
   soHandle = NULL;
}


SxPluginLoader::SymPtr SxPluginLoader::getSymbol (const SxString &name)
{

#  ifdef WIN32
      const char *(*res)() = NULL;
	  *(void **)(&res) = GetProcAddress (soHandle, name.ascii());
      if (!res)  {
         sxprintf ("ERROR: Can't access plugin symbol %s\n", name.ascii());
         SX_EXIT;
      }

#  else
      dlerror (); // clean previous errors

      SymPtr res = NULL;
      res = sxFnCast<SxPluginLoader::SymPtr> (dlsym (soHandle, name.ascii()));
      //res = dlsym (soHandle, name.ascii()); // warning using -pedantic
      //*(void **)(&res) = dlsym (soHandle, name.ascii());  // avoid warning

      const char *error = NULL;
      if ((error = dlerror()) != NULL)  {
         sxprintf ("ERROR: %s\n", error);
         SX_EXIT;
      }
      dlerror (); // clean previous errors
#  endif /* WIN32 */

   return res;
}
