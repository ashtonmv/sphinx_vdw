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

#include <sxfileno.h>

int sxfileno (FILE *fp)
{
#  ifdef WIN32
      return _fileno (fp);
#  else
      return fileno (fp);
#  endif
}
