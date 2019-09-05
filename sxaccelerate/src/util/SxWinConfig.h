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


#ifndef _SX_WIN_CONFIG_H_
#define _SX_WIN_CONDIG_H_

#ifdef MSVC
#  pragma warning(disable:4065)
#  pragma warning(disable:4244)
#  pragma warning(disable:4251)
#  pragma warning(disable:4290)
#  pragma warning(disable:4305)
#  pragma warning(disable:4355)
#  pragma warning(disable:4996)
#endif /* MSVC */

#  ifdef min
#    undef min
#  endif

#  ifdef max
#    undef max
#  endif

#ifdef MSVC
#  include <io.h>
#  define _Exit  _exit
#  define isatty _isatty

   // Windows doesn't even know what ssize_t is!
#ifdef WIN64
   typedef long long ssize_t;
#  define GDB "\"C:\\Program Files (x86)\\Windows Kits\\10\\Debuggers\\x64\\windbg.exe\" -z"
#else
   typedef long ssize_t;
#  define GDB "\"C:\\Program Files (x86)\\Windows Kits\\10\\Debuggers\\x86\\windbg.exe\" -z"
#endif
#endif /* MSVC */

#  define SXDATE     "Unknown date"
#  define SVNTAG     "Unknown SVN tag"
#  define WHOAMI     "Unknown Windows user"
#  define CXX        "MSC"
#  define CXXVERSION "Unknown version"
#  define CXXFLAGS   "Unknown flags"
#  define LDFLAGS    "Unknown LDFLAGS"
#  define MEMTRACER  "Undefined"
#  define USE_SX_LOG 1

// --- include SxUnixConfig.h in MinGW compilation
#ifdef MSVC
#  include <SxWinEnv.h>
#else
#  include <SxUnixConfig.h>
#endif


#endif /* _SX_WIN_CONFIG_H_ */
