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

/*
   1. UTF-8: http://en.wikipedia.org/wiki/UTF-8 and RFC3629
   2. UTF-16: http://en.wikipedia.org/wiki/UTF-16 and RFC2781
   3. 7-bit ASCII characters have the same byte representation in UTF-8 as in
      Latin-1.
   4. On Windows, the data type 'wchar_t' is 16 bits wide and unsigned and
      represents UTF-16 encoded characters.
*/

#include <SxChar.h>
#include <SxUnicodeTbl.h>

#include <stdlib.h>

SxConstChar::SxConstChar (const char *start_, const char *limit_,
   bool isUnicode_)
{
   initialize (start_, limit_, isUnicode_);
}

SxConstChar::SxConstChar (const char *start_, ssize_t nBytes, bool isUnicode_)
{
   initialize (start_, start_ + nBytes, isUnicode_);
}

SxConstChar::~SxConstChar ()
{
   // empty
}

void SxConstChar::initialize (const char *start_, const char *limit_,
   bool isUnicode_)
{
   SX_CHECK (start_);
   SX_CHECK (limit_);
   SX_CHECK (limit_ >= start_, (const void *) start_,
      (const void *) limit_);
   // SX_CHECK (*limit_ == '\0', (const void *) limit_); -- not always
   start = start_;
   limit = limit_;
   isUnicode = isUnicode_;
   if ( (isUnicode) && (start < limit) )
      SX_CHECK (!isUtf8Continuation (*start), (uint8_t) *start);
   isDirty = true;
}

bool SxConstChar::is7bit (const char *str)
{
   if (str == NULL)  return true;
   char ch;
   while ( (ch = *str++) != '\0' )  {
      if (!is7bit (ch))  return false; // found a set 8th bit
   }
   return true;
}

uint8_t SxConstChar::analyzeText (const char *str, ssize_t nBytes)
{
   uint8_t retval = 0;
   char ch;
   if (str == NULL)  return 0;
   if (nBytes == -1)  nBytes = (ssize_t) ::strlen (str);
   for (; (nBytes > 0) && ( (ch = *str++) != '\0' ); nBytes--)  {
      if (is7bit (ch))  { } // 7-bit ASCII; nice
      else if (isUtf8Lead (ch) || isUtf8Continuation (ch))  { // likely UTF-8
         if (retval < 1)  retval = 1;
      }
      else  { // unknown; likely 8-bit ASCII
         retval = 2;
         break;
      }
   }
   return retval;
}

SxConstChar::UnicodePoint SxConstChar::utf8ToChar (const char **ptr_)
{
   SX_CHECK (ptr_);
   const char *ptr = *ptr_;
   SX_CHECK (ptr);
   SX_CHECK (!isUtf8Continuation (*ptr), (uint8_t) *ptr);
   SxConstChar::UnicodePoint res;
   char ch = *ptr++, ch2, ch3, ch4;
   if ((ch & 128) == 0)
      res = static_cast<SxConstChar::UnicodePoint>(ch); // the most likely case: 7-bit ASCII char
   else  {
      res = 0;
      if ((ch & (7 << 5)) == (6 << 5))  {
         if ( ((ch2 = *ptr) & (3 << 6)) == (2 << 6) )  {
            res = static_cast<SxConstChar::UnicodePoint>(((ch & 31) << 6) | (ch2 & 63));
            ptr++;
         }
      }
      else if ((ch & (15 << 4)) == (14 << 4))  {
         if ( ( ((ch2 = ptr[0]) & (3 << 6)) == (2 << 6) ) &&
              ( ((ch3 = ptr[1]) & (3 << 6)) == (2 << 6) ) )  {
            res = static_cast<SxConstChar::UnicodePoint>(((ch & 15) << 12)
		                                         | ((ch2 & 63) << 6)
                                                         | (ch3 & 63));
            ptr += 2;
         }
      }
      else if ((ch & (31 << 3)) == (30 << 3))  {
         if ( ( ((ch2 = ptr[0]) & (3 << 6)) == (2 << 6) ) &&
              ( ((ch3 = ptr[1]) & (3 << 6)) == (2 << 6) ) &&
              ( ((ch4 = ptr[2]) & (3 << 6)) == (2 << 6) ) )  {
            res = static_cast<SxConstChar::UnicodePoint>(((ch & 7) << 18)
                                                         | ((ch2 & 63) << 12)
                                                         | ((ch3 & 63) << 6)
                                                         | (ch4 & 63));
            ptr += 3;
         }
      }
      SX_CHECK (res > 0, res);
   }
   *ptr_ = ptr;
   return res;
}

int SxConstChar::charToUtf8 (UnicodePoint u, char *origDest)
{
   SX_CHECK (origDest);
   char *dest = origDest;
   if (u < 0x80)  { // the most likely case: a 7-bit ASCII character
      *dest++ = (char) u;
   }  else if (u < 0x800)  {
      *dest++ = (char) ((3 << 6) | ((u >> 6) & 31));
      *dest++ = (char) ((1 << 7) | (u & 63));
   }  else if (u < 0x10000)  {
      if ( (u < 0xd800) || (u > 0xdfff) )  { // (avoid surrogate pairs)
         *dest++ = (char) ((7 << 5) | ((u >> 12) & 15));
         *dest++ = (char) ((1 << 7) | ((u >> 6) & 63));
         *dest++ = (char) ((1 << 7) | (u & 63));
      }
   }  else if (u < 0x110000)  {
      *dest++ = (char) ((15 << 4) | ((u >> 18) & 7));
      *dest++ = (char) ((1 << 7) | ((u >> 12) & 63));
      *dest++ = (char) ((1 << 7) | ((u >> 6) & 63));
      *dest++ = (char) ((1 << 7) | (u & 63));
   }  else  {
      cout << "Invalid Unicode character " << u << endl;
      SX_EXIT;
   }
   return (int) (dest - origDest);
}

ssize_t SxConstChar::countChars (const char *start_, const char *limit_,
   bool isUnicode_)
{
   SX_CHECK (start_);
   SX_CHECK (limit_);
   SX_CHECK (start_ <= limit_, (const void *) start_,
      (const void *) limit_);
   return (isUnicode_ ? getNUtf8Chars (start_, limit_) : (limit_ - start_));
}

ssize_t SxConstChar::getNChars () const
{
   if (isDirty)  {
      if (isUnicode)  nChars = getNUtf8Chars (start, limit);
      else  nChars = getNBytes ();
      isDirty = false;
   }
   return nChars;
}

void SxConstChar::utf8Inc (const char **ptr_, ssize_t count)
{
   SX_CHECK (ptr_);
   SX_CHECK (count >= 0, count);
   if (count == 0)  return;
   const char *ptr = *ptr_;
   SX_CHECK (ptr);
   SX_CHECK (!isUtf8Continuation (*ptr), (uint8_t) *ptr);
   while (count-- > 0)  {
      SX_CHECK (*ptr != '\0', (uint8_t) *ptr); // don't leave the string
      do { ptr++; } while (!SxConstChar::isUtf8CharStart (*ptr));
   }
   *ptr_ = ptr;
}

void SxConstChar::inc (const char **ptr_, ssize_t count, bool isUnicode_)
{
   SX_CHECK (ptr_);
   SX_CHECK (count >= 0, count);
   if (count == 0)  return;
   if (!isUnicode_)  {
      const char *ptr = *ptr_;
      SX_CHECK (ptr);
      ptr += count;
      *ptr_ = ptr;
   }
   else  utf8Inc (ptr_, count);
}

void SxConstChar::utf8Dec (const char **ptr_, ssize_t count)
{
   SX_CHECK (ptr_);
   SX_CHECK (count >= 0, count);
   if (count == 0)  return;
   const char *ptr = *ptr_;
   SX_CHECK (ptr);
   SX_CHECK (!isUtf8Continuation (*ptr), (uint8_t) *ptr);
   while (count-- > 0)  {
      do { ptr--; } while (!SxConstChar::isUtf8CharStart (*ptr));
   }
   *ptr_ = ptr;
}

ssize_t SxConstChar::getNUtf8Chars (const char *start_, const char *limit_)
{
   SX_CHECK (start_);
   if (limit_ == NULL)  limit_ = start_ + ::strlen (start_);
   ssize_t res = 0;
   for (const char *ptr = start_; ptr < limit_; res++)  utf8Inc (&ptr);
   return res;
}

ssize_t SxConstChar::getNCharsToNBytes (const char *start_, ssize_t nChars_,
   bool isUnicode_)
{
   SX_CHECK (start_);
   SX_CHECK (nChars_ >= 0, nChars_);
   if (!isUnicode_)  return nChars_;
   const char *ptr = start_;
   utf8Inc (&ptr, nChars_);
   return ptr - start_;
}

SxConstChar::UnicodePoint SxConstChar::combineUtf16 (const uint16_t &lead,
   const uint16_t &trail)
{
   SxConstChar::UnicodePoint u = static_cast<SxConstChar::UnicodePoint>((lead
			         & sxUtf16Mask) << sxUtf16Shift);
   return (u | (trail & sxUtf16Mask)) + 0x10000;
}

int SxConstChar::utf16ToUtf8 (const uint16_t **origSrc, char *buffer)
{
   SX_CHECK (origSrc);
   SX_CHECK (buffer);
   const uint16_t *src = *origSrc;
   SX_CHECK (src);
   uint16_t w = *src++;
   SxConstChar::UnicodePoint u;
   if (SxConstChar::isUtf16LeadVal (w))  {
      uint16_t w2 = *src++;
      SX_CHECK (SxConstChar::isUtf16TrailVal (w2), w2);
      u = SxConstChar::combineUtf16 (w, w2);
   }
   else  u = (SxConstChar::UnicodePoint) w;
   *origSrc = src;
   return SxConstChar::charToUtf8 (u, buffer);
}

// (The algorithm requires that wchar_t is UTF-16, as it is on Windows.)
SxArray<char> SxConstChar::wcharsToUtf8 (const uint16_t *origStr,
   ssize_t nWchars)
{
   SxArray<char> res;
   if (origStr == NULL)  return res;
   const uint16_t *str;
   if (nWchars == -1)  {
      nWchars = 0;
      str = origStr;
      while (*str++ != 0)  nWchars++;
   }
   if (nWchars < 1)  return res;
   char dummyBuf[4];
   ssize_t newNBytes = 1;
   str = origStr;
   while (*str != 0)  newNBytes += utf16ToUtf8 (&str, dummyBuf); // just count
   res.resize (newNBytes);
   char *dest = res.elements;
   str = origStr;
   while (*str != 0)  dest += utf16ToUtf8 (&str, dest); // store byte sequence
   *dest = '\0';
   return res;
}

SxArray<uint16_t> SxConstChar::asciiToWChars (const char *src, ssize_t nBytes_)
{
   SxArray<uint16_t> res;
   if (src == NULL)  return res;
   if (nBytes_ == -1)  nBytes_ = (ssize_t) ::strlen (src);
   if (nBytes_ < 1)  return res;
   ssize_t nWchars = nBytes_ + 1;
   res.resize (nWchars);
   uint16_t *dest = res.elements;
   while (nBytes_-- > 0)  *dest++ = (uint16_t) ((uint8_t) (*src++));
   *dest = 0;
   return res;
}

// --- Convert an UTF-8 encoded string to a wchar_t (UTF-16) sequence
SxArray<uint16_t> SxConstChar::toWChars (const char *src, ssize_t nBytes_,
   bool isUnicode_)
{
   SxArray<uint16_t> res;
   if ( (src == NULL) || (nBytes_ < 1) || (*src == '\0') )  return res;
   if (!isUnicode_)  return SxConstChar::asciiToWChars (src, nBytes_);

   ssize_t nWchars = 1;
   SxConstChar str(src, nBytes_, true);
   SxConstChar::Iterator it = str.begin (), itEnd = str.end ();
   for (; it != itEnd; ++it)  { // just count
      UnicodePoint u = *it;
      if (u <= 0xffff)  nWchars++;
      else  nWchars += 2; // need a surrogate pair
   }
   res.resize (nWchars);
   uint16_t *dest = res.elements;
   for (it.resetPos (); it != itEnd; ++it)  { // convert and store the values
      UnicodePoint u = *it;
      if (u <= 0xffff)  *dest++ = (uint16_t) u;
      else  { // need a surrogate pair
         u -= 0x10000;
         SX_CHECK (u < (1 << 20), u); // value must be encodable
         *dest++ = (uint16_t) (u >> sxUtf16Shift); // leading surrogate
         *dest++ = (uint16_t) (u & sxUtf16Mask); // trailing surrogate
      }
   }
   *dest = 0;
   return res;
}

ssize_t SxConstChar::getNUtf8Bytes (const char *src, ssize_t nAsciiBytes)
{
   SX_CHECK (src);
   ssize_t res = 0;
   while (nAsciiBytes-- > 0)  { // calculate the resulting number of bytes
      char dummyBuf[4];
      uint8_t ch = (uint8_t) (*src++);
      SX_CHECK (ch != '\0', ch);
      res += charToUtf8 ((SxConstChar::UnicodePoint) ch, dummyBuf);
   }
   return res;
}

ssize_t SxConstChar::asciiToUtf8 (char *origDest, const char *src,
                                 ssize_t nAsciiBytes)
{
   SX_CHECK (origDest);
   SX_CHECK (src);
   char *dest = origDest;
   while (nAsciiBytes-- > 0)  { // store the UTF-8 byte sequence
      uint8_t ch = (uint8_t) (*src++);
      SX_CHECK (ch != '\0', ch);
      dest += charToUtf8 ((SxConstChar::UnicodePoint) ch, dest);
   }
   return dest - origDest;
}

int SxConstChar::encode (SxConstChar::UnicodePoint u, bool isUnicode_,
   char *buffer)
{
   if (isUnicode_)  return charToUtf8 (u, buffer);
   SX_CHECK (u > 0 && u <= 255, u);
   *buffer = (char) u;
   return 1; // we copied one byte
}

void SxConstChar::set (char *dest, const char *src, ssize_t srcNBytes,
  ssize_t count)
{
   SX_CHECK (dest);
   SX_CHECK (src);
   SX_CHECK (srcNBytes > 0, srcNBytes);
   SX_CHECK (count >= 0, count);
   while (--count >= 0)  {
      for (ssize_t idx = 0; idx < srcNBytes; idx++)  *dest++ = src[idx];
   }
}

int SxConstChar::caseConvDiffCB (const void *a, const void *b)
{
   const CaseConvItem *A = (const CaseConvItem *) a,
      *B = (const CaseConvItem *) b;
   return ((int) A->from) - ((int) B->from);
}

SxConstChar::UnicodePoint
SxConstChar::caseConvLookup (const CaseConvItem *tbl, size_t nElems,
                             UnicodePoint from)
{
   CaseConvItem key;
   key.from = from;
   const CaseConvItem *entry =
      (const CaseConvItem *) ::bsearch (&key, tbl, nElems,
      sizeof (CaseConvItem), caseConvDiffCB);
   return (entry ? (entry->to) : from );
}

SxConstChar::UnicodePoint SxConstChar::toLowerTblLookup (
   SxConstChar::UnicodePoint from)
{
   return caseConvLookup (sxToLowerTbl, sizeof (sxToLowerTbl) /
      sizeof ((sxToLowerTbl)[0]) - 1, from);
   // ("-1" to ignore the trailing "{ 0, 0 }" entry which just avoids trailing
   // commas)
}

SxConstChar::UnicodePoint SxConstChar::toUpperTblLookup (
   SxConstChar::UnicodePoint from)
{
   return caseConvLookup(sxToUpperTbl, sizeof (sxToUpperTbl) /
      sizeof ((sxToUpperTbl)[0]) - 1, from);
   // ("-1" to ignore the trailing "{ 0, 0 }" entry which just avoids trailing
   // commas)
}

void SxConstChar::makeGap (char *str, ssize_t fromByteIdx, ssize_t nInsBytes,
   ssize_t oldNBytes)
{
   SX_CHECK (str);
   for (ssize_t i = oldNBytes - 1; i >= fromByteIdx; i--)  {
      str[i + nInsBytes] = str[i];
   }
}


// ---------------------------------------------------------------------------

SxConstChar::Iterator::Iterator (const SxConstChar *obj)
{
   SX_CHECK (obj);
   charObj = obj;
   resetPos ();
}

SxConstChar::Iterator::~Iterator ()
{
   // empty
}

void SxConstChar::Iterator::setCharIdx (ssize_t charIdx_)
{
   SX_CHECK (charIdx_ >= 0, charIdx_);
   SX_CHECK (charIdx_ <= charObj->getNChars (), charIdx_,
             charObj->getNChars ());
   if (!charObj->isUnicode)  {
      SX_CHECK (charIdx_ <= charObj->getNBytes (), charIdx_,
         charObj->getNBytes ());
      charIdx = byteIdx = charIdx_;
   }  else  {
      charIdx = charIdx_;
      byteIdx = SxConstChar::getNCharsToNBytes (charObj->start, charIdx_, true);
   }
}

void SxConstChar::Iterator::setByteIdx (ssize_t byteIdx_)
{
   SX_CHECK (byteIdx_ >= 0, byteIdx_);
   ssize_t nBytes = charObj->getNBytes ();
   SX_CHECK (byteIdx_ <= nBytes, byteIdx_);
   byteIdx = byteIdx_;
   if (!charObj->isUnicode)  charIdx = byteIdx;
   else  {
      if (byteIdx_ < nBytes)  SX_CHECK (isUtf8CharStart (*(getPtr ())));
      charIdx = getNUtf8Chars (charObj->start, charObj->start + byteIdx);
   }
}

SxConstChar::UnicodePoint SxConstChar::Iterator::operator* () const
{
   const char *ptr = getPtr ();
   SX_CHECK (ptr < charObj->limit, (const void *) ptr,
      (const void *) charObj->limit);
   if (!charObj->isUnicode)  return (UnicodePoint) ((uint8_t) (*ptr));
   return utf8ToChar (&ptr);
}

bool SxConstChar::Iterator::contains (const SxConstChar::UnicodePoint &u)
{
   SxConstChar::Iterator hayIt = *this;
   for (; hayIt.inRange (); ++hayIt)  {
      if (*hayIt != u)  continue;
      *this = hayIt;
      return true; // found
   }
   return false; // not found
}

bool SxConstChar::Iterator::contains (const SxConstChar &needle)
{
   // The standard C library function ::strstr () shouldn't be used here
   // because it might access the string beyond 'limit', possibly accessing
   // invalid data.
   SX_CHECK (needle.getNBytes () > 0); // can't find "nothing"
   SxConstChar::Iterator needleStartIt = needle.begin ();
   SxConstChar::UnicodePoint needleStartU = *needleStartIt;
   SxConstChar::Iterator hayIt = *this;
   while (1)  {
      if (!hayIt.contains (needleStartU))  return false; // quick anchor scan
      SxConstChar::Iterator hayMatchStartIt = hayIt;
      SxConstChar::Iterator needleIt = needleStartIt;
      do  {
         ++needleIt;
         if (!needleIt.inRange ())  {
            // matched the entire needle; set the iterator position to the
            // start of the match:
            *this = hayMatchStartIt;
            return true;
         }
         ++hayIt;
         if (!hayIt.inRange ())  return false; // not enough hay remaining
      }  while (*hayIt == *needleIt);
      hayIt = hayMatchStartIt + 1;
   }
   return false; // not found (unreachable; just to suppress compiler warnings)
}

bool SxConstChar::Iterator::contains (const char *needle_)
{
   SX_CHECK (needle_);
   SxConstChar needle(needle_, needle_ + ::strlen (needle_), false);
   return contains (needle);
}

int SxConstChar::codePointDiffCB (const void *a, const void *b)
{
   const SxConstChar::UnicodePoint *A = (const SxConstChar::UnicodePoint *) a,
      *B = (const SxConstChar::UnicodePoint *) b;
   return ((int) *A) - ((int) *B);
}

bool SxConstChar::Iterator::isBlank () const
{
   if (!charObj->isUnicode)  return SxConstChar::isBlank (*(getPtr ()));
   UnicodePoint u = **this;
   const UnicodePoint *entry = (const UnicodePoint *) ::bsearch (&u, sxBlankTbl,
      sizeof (sxBlankTbl) / sizeof (sxBlankTbl[0]) - 1, sizeof (UnicodePoint),
      codePointDiffCB);
   return (entry != NULL);
}

bool SxConstChar::Iterator::isLineTerm () const
{
   if (!charObj->isUnicode)  return *(getPtr ()) == '\n';
   UnicodePoint u = **this;
   for (uint8_t idx = 0;
      idx < sizeof (sxLineTermTbl) / sizeof ((sxLineTermTbl)[0]) - 1;
      idx++)  {
      if (sxLineTermTbl[idx] == u)  return true; // found
      // (::bsearch () wouldn't make much sense for such a tiny array)
   }
   return false; // not found
}

bool SxConstChar::Iterator::isNewline (bool skip)
{
   if (*(getPtr ()) == '\r')  {
      if (skip)  {
         inc ();
         if (*(getPtr ()) == '\n')  inc (); // CRLF
      }
      return true;
   }  else if (isLineTerm ())  { // (e.g. a simple '\n')
      if (skip)  inc ();
      return true;
   }
   return false;
}

bool SxConstChar::Iterator::isDigit () const
{
   UnicodePoint u = **this;
   return ( (u >= '0') && (u <= '9') );
}

void SxConstChar::Iterator::encode (char **dest_) const
{
   SX_CHECK (dest_);
   char *dest = *dest_;
   SX_CHECK (dest);
   SX_CHECK (inRange ());
   const char *src = getPtr ();
   if (!charObj->isUnicode)  *dest++ = *src;
   else  {
      do  {
         *dest++ = *src++;
      }  while ( (src < charObj->limit) && (!isUtf8CharStart (*src)) );
   }
   *dest_ = dest;
}

SxConstChar::Iterator SxConstChar::begin () const
{
   return SxConstChar::Iterator(this);
}

SxConstChar::Iterator SxConstChar::end () const
{
   SxConstChar::Iterator res(this);
   res.setByteIdx (getNBytes ());
   return res;
}

void SxConstChar::Iterator::inc (ssize_t count)
{
   SX_CHECK (count >= 0, count);
   SX_CHECK (charIdx + count <= charObj->getNChars (), charIdx, count,
      charObj->getNChars ()); // (We may reach 'limit', but never pass it.)
   if (count == 0)  return;
   if (!charObj->isUnicode)  byteIdx += count;
   else  {
      const char *origPtr = charObj->start + byteIdx, *ptr = origPtr;
      utf8Inc (&ptr, count);
      byteIdx += (ptr - origPtr);
   }
   SX_CHECK (byteIdx <= charObj->getNBytes ());
      // (We may reach 'limit', but never pass it.)
   charIdx += count;
}

void SxConstChar::Iterator::dec (ssize_t count)
{
   SX_CHECK (count >= 0, count);
   SX_CHECK (charIdx >= count, charIdx, count);
   if (!charObj->isUnicode)  byteIdx -= count;
   else  {
      const char *origPtr = charObj->start + byteIdx, *ptr = origPtr;
      utf8Dec (&ptr, count);
      byteIdx -= (origPtr - ptr);
   }
   SX_CHECK (byteIdx >= 0, byteIdx);
   charIdx -= count;
}


// ---------------------------------------------------------------------------

SxConstChar::Iterator operator+ (const SxConstChar::Iterator &from,
                                 ssize_t nChars_)
{
   SxConstChar::Iterator res = from;
   res += nChars_;
   return res;
}

SxConstChar::Iterator operator- (const SxConstChar::Iterator &from,
                                 ssize_t nChars_)
{
   SxConstChar::Iterator res = from;
   res -= nChars_;
   return res;
}

