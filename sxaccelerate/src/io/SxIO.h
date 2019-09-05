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
#ifndef _SX_IO_H_
#define _SX_IO_H_

#ifdef WIN32
#  if defined(_EXPORT_sxio)
#     define SX_EXPORT_IO __declspec(dllexport)
#  else
#     define SX_EXPORT_IO __declspec(dllimport)
#  endif
#else
#  define SX_EXPORT_IO
#endif

#include <SxString.h>
/** \brief Declare filename as output file
    Some output files may be opened and closed multiple times during
    a single run, appending some data each time. Of course, preexisting
    files must be deleted. Sometimes, it is not clear where the initial
    delete should go (e.g. multiple geometry optimizers).
    This routine deletes the file exactly once. Further calls with the
    same filename have no effect.
    For this, the routine keeps internally a list of filenames.
    @filename

*/
void SX_EXPORT_IO deleteFileOnce (const SxString &filename);
#define SX_OUTPUT_FILE(file) deleteFileOnce(file)

/// \brief Report a failure to open a file (used by #sxfopen)
void SX_EXPORT_IO sxfopenError(const char *name, const char *mode);

/** \brief Open a file, crash with error message upon failure
    @param name file name
    @param mode fopen mode "a", "w", or "r"
    @return a file pointer, non-NULL
  */
inline FILE* sxfopen(const char *name, const char *mode)
{
   SX_CHECK(name);
   SX_CHECK(mode);
   SX_CHECK(mode[0] == 'a' || mode[0] == 'r' || mode[0] == 'w');
   FILE *fp = fopen (name, mode);
   if (!fp) sxfopenError (name, mode);
   return fp;
}

/** \brief Open a file, crash with error message upon failure
    @param name file name
    @param mode fopen mode "a", "w", or "r"
    @return a file pointer, non-NULL
  */
inline FILE* sxfopen(const SxString &name, const char *mode)
{
   return sxfopen(name.ascii (), mode);
}




#endif /* _SX_IO_H_ */
