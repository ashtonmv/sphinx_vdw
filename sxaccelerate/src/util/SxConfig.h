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

#ifndef _SX_CONFIG_H_
#define _SX_CONFIG_H_

#ifdef MSVC
#  include <SxWinConfig.h>
#else
#  include <SxUnixConfig.h>
#endif

#define BRANCH "undefined"


#endif /* _SX_CONFIG_H_ */
