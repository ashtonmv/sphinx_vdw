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
#ifndef _SX_IDX_H_
#define _SX_IDX_H_

#include <SxUtil.h>

/**
  @ingroup Numerics
  */
class SxIdx
{
   public:
      int start;
      int end;

      inline SxIdx ()                    { start = end = 0; }
      inline SxIdx (int idx0, int idx1)  
      { 
         SX_CHECK (0 <= idx0 && idx0 <= idx1, idx0, idx1);
         start = idx0; end = idx1; 
      }
};

#endif /* _SX_IDX_H_ */
