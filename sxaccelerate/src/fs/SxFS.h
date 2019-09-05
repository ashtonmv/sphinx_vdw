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
#ifndef _SX_FS_H_
#define _SX_FS_H_

#ifdef WIN32
#  include <SxConfig.h>
#  if defined(_EXPORT_sxfs)
#     define SX_EXPORT_FS __declspec(dllexport)
#  else
#     define SX_EXPORT_FS __declspec(dllimport)
#  endif

#  ifndef S_IRUSR
#     define S_IRUSR    _S_IREAD
#     define S_IWUSR    _S_IWRITE
#     define S_IXUSR    _S_IEXEC
#  endif
#  ifndef S_ISREG
#     define S_ISDIR(m)  (((m) & S_IFMT) == S_IFDIR)
#     define S_ISFIFO(m) (((m) & S_IFMT) == S_IFIFO)
#     define S_ISREG(m)  (((m) & S_IFMT) == S_IFREG)
#  endif

// --- redefine some POSIX functions to Windows style
#  include <direct.h>
#  include <io.h>
inline int chdir (const char *p)         { return _chdir (p); }
inline char *getcwd (char *b, size_t s)  { return _getcwd (b,(int)s); }
inline int rmdir (const char *d)         { return _rmdir (d); }
inline int umask (int m)                 { return _umask (m); }

#  define mode_t     int

#else
#  define SX_EXPORT_FS
#endif

#ifdef WIN32
typedef struct _stati64 SxStructStat;
#else // WIN32
typedef struct stat SxStructStat;
#endif

extern int sxstat  (const char *path, SxStructStat *buf);
#ifdef WIN32
   extern int sxwstat (const wchar_t *path, SxStructStat *buf);
#endif
extern int sxfstat (int fd, SxStructStat *buf);

#endif /* _SX_FS_H_ */
