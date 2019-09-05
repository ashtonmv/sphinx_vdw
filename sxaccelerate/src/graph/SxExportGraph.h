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

#ifndef _SX_EXPORT_GRAPH_H_
#define _SX_EXPORT_GRAPH_H_

#include <SxConfig.h>

#ifdef WIN32
#  if defined(_EXPORT_sxgraph)
#     define SX_EXPORT_GRAPH __declspec(dllexport)
#  else
#     define SX_EXPORT_GRAPH __declspec(dllimport)
#  endif
#else
#  define SX_EXPORT_GRAPH
#endif

class SX_EXPORT_GRAPH SxGraphLnkDummy
{
   public:
      static void dummy ();
};

#endif /* _SX_EXPORT_GRAPH_H_ */
