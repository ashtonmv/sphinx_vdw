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
#ifndef _SX_UTIL_H_
#define _SX_UTIL_H_

#include <math.h>
#include <stdlib.h>
#include <ostream>
#include <SxConfig.h>
#include <SxMacroLib.h>
#include <atomic>

#ifdef WIN32
#  if defined(_EXPORT_sxutil)
#     define SX_EXPORT_UTIL __declspec(dllexport)
#  else
#     define SX_EXPORT_UTIL __declspec(dllimport)
#  endif
#else
#  define SX_EXPORT_UTIL
#endif

// --- 1st half (64bit) of shared object ID ("sxutil")
//     md5 -s "sxutil" | echo -n "sxutil" | md5sum
#define SX_COMPONENT_ID 0x588a48a3cb5464f6#795428de46d03c75

#ifdef WIN32
#  define _WINSOCKAPI_ // prevent winsock.h to be included from windows.h
                       // otherwise later include of winsock2.h would show
                       // redefinition errors
#  include <windows.h> // DWORD, GetLastError()
#else
#   include <pthread.h>
#endif /* WIN32 */


class SX_EXPORT_UTIL SxUtil
{
   public:
      SxUtil ();
     ~SxUtil ();

     static SxUtil &getGlobalObj ();

     int nProcs;

#    ifdef USE_SX_LOG
#       ifdef WIN32
           CRITICAL_SECTION logMutex;
#       else
           pthread_mutexattr_t logMutexAttr;
           pthread_mutex_t     logMutex;
#       endif

         void lockLog ();
         void unlockLog ();
#    endif
};

inline int SX_EXPORT_UTIL toInt (float arg)  { return (int)(floor(arg+.5)); }
inline int SX_EXPORT_UTIL toInt (double arg) { return (int)(floor(arg+.5)); }

extern void SX_EXPORT_UTIL sxGetExecPath (char *res, size_t len); // #exception
extern void SX_EXPORT_UTIL sxchmod (const char *path, int mode); // #exception

extern char __flags_;
extern void SX_EXPORT_UTIL setNProcs (int n);
extern int  SX_EXPORT_UTIL getNProcs ();
extern long SX_EXPORT_UTIL sxChunkSize;

#endif /* _SX_UTIL_H_ */
