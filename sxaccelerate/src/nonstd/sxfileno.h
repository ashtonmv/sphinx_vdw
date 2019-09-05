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

#ifndef _SX_FILENO_H_
#define _SX_FILENO_H_

#include <SxNonStd.h>
#include <stdio.h>

SX_EXPORT_NONSTD int sxfileno (FILE *);

// --- overwrite fileno used by bison
#define fileno sxfileno

#endif /* _SX_FILENO_H_ */
