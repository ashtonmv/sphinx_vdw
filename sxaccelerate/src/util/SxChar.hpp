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

inline void SxConstChar::Iterator::inc (ssize_t count)
{
#if USE_CHAR_TIMING
   CharTiming::timing.nIncCalls++;
#endif
   SX_CHECK (count >= 0, count);
   if (charObj->isUnicode)  {
      const char *origPtr = charObj->start + byteIdx, *ptr = origPtr;
      utf8Inc (&ptr, count);
      byteIdx += ptr - origPtr;
   }
   else  byteIdx += count;
   SX_CHECK (byteIdx <= charObj->getNBytes ());
      // (We may reach 'limit', but never pass it.)
   charIdx += count;
}
