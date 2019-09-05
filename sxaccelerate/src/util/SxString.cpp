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

#include <SxString.h>
#include <SxLog.h>
#include <SxException.h>

//#include <SxCalc.h> not used?
#include <sys/stat.h>
#include <string.h>
#include <stdarg.h>
#include <cmath>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>

#ifdef WIN32
#  include <io.h>
   static int sxInitUnicodeEnv() {
      SetConsoleOutputCP(CP_UTF8);
      SetConsoleCP(CP_UTF8);
      setvbuf(stdout, nullptr, _IOFBF, 1000);
      setvbuf(stderr, nullptr, _IOFBF, 1000);
      return 0;
   }

   static int initVal = sxInitUnicodeEnv();
#endif /* WIN32 */

// By default, the Android system sends stdout and stderr output to /dev/null.
// http://developer.android.com/tools/debugging/debugging-log.html
// sxprintf() will use __android_log_write() from android/log.h
#ifdef SX_ANDROID
#  include <android/log.h>
#endif /* SX_ANDROID */

#ifndef WIN32
#  include <sys/ioctl.h>



#endif /* WIN32 */

const char *SxString::emptyString = "";
//TODO Substitute 10240 by the actual maximal file size in bytes
const ssize_t SxString::MAX_N_CHARS_FILE = 10240;
// Maximal number of characters needed by one integer variable
#define MAX_N_CHARS_INT 64
// Maximal number of characters needed by one float variable
#define MAX_N_CHARS_FLOAT 512
// Maximal number of characters needed by one double variable
#define MAX_N_CHARS_DOUBLE 512

//------------------------------------------------------------------------------

SxString::SxString ()
  : SxArray<char> (), nChars(0), isUnicode_(false), isDirty(false)
{
   // empty
}

SxString::SxString (const char c)
  : SxArray<char> (), nChars(0), isUnicode_(false), isDirty(false)
{
   resize0 (1);
   if (c == 0) SxArray<char>::operator()(0) = '0';
   else        SxArray<char>::operator()(0) = c;
}

SxString::SxString (const char *str)
  : SxArray<char> (), nChars(0), isUnicode_(false), isDirty(false)
{
   if (str)  {
      replace (str);
   }
}

SxString::SxString (const char *str, ssize_t len)
  : SxArray<char> (), nChars(0), isUnicode_(false), isDirty(false)
{
   replace (str, len);
}

SxString::SxString (const uint16_t *str)
  : SxArray<char> (), nChars(0), isUnicode_(false), isDirty(false)
{
   isUnicode_ = true;
   if (str)  {
      replace (str);
   }
}

SxString::SxString (const Buffer &buf)
  : SxArray<char> (), nChars(0), isUnicode_(false), isDirty(false)
{
   ssize_t nChars_, maxNChars = buf.getMaxNChars (), maxNBytes, actualNBytes;
   const char *src;
   const uint16_t *usrc;
   char *dest;
   Buffer::Mode mode = buf.getMode ();
   switch (mode)  {
      case Buffer::ASCII:
         src = buf.getBuffer ();
         nChars_ = (ssize_t) ::strlen (src);
         if (nChars_ > 0)  {
            resize0 (nChars_, false);
            ::memcpy (elements, src, static_cast<size_t>(nChars_));
         }
         break;
      case Buffer::UTF8:
         *this = unicodeFromUtf8 (buf.getBuffer ());
         break;
      case Buffer::UTF16:
         maxNBytes = 4 * maxNChars;
         resize0 (maxNBytes, false);
         dest = elements;
         usrc = reinterpret_cast<const uint16_t *>(buf.getBuffer ());
         while (*usrc != 0)  dest += SxConstChar::utf16ToUtf8 (&usrc, dest);
         actualNBytes = dest - elements;
         if (actualNBytes < maxNBytes)  resize0 (actualNBytes, true); // shrink
         break;
      default:
         SX_EXIT; // invalid mode
         break;
   }
}

#ifdef SX_IOS
SxString::SxString (CFStringRef str)
  : SxArray<char> (), nChars(0), isUnicode_(false), isDirty(false)
{
   if (str) {
      replace (str);
   }
}
#endif /* SX_IOS */

SxString::SxString (const char *s1, const char *s2)
  : SxArray<char> (), nChars(0), isUnicode_(false), isDirty(false)
{
   if (s1 || s2)  {
      concatenate (s1, s2);
   }
}

SxString::SxString (const SxString &s1, const SxString &s2)
  : SxArray<char> (), nChars(0), isUnicode_(false), isDirty(false)
{
   concatenate (s1, s2);
   if (  s1.isUnicode ()
      || s2.isUnicode ()) {
      isUnicode_ = true;
   }
}

SxString::SxString (int i)
  : SxArray<char> (), nChars(0), isUnicode_(false), isDirty(false)
{
   char buffer[MAX_N_CHARS_INT];
   int len = SxString::snprintf (buffer, MAX_N_CHARS_INT, i);
   replace (buffer, len);
}

SxString::SxString (unsigned int i)
  : SxArray<char> (), nChars(0), isUnicode_(false), isDirty(false)
{
   char buffer[MAX_N_CHARS_INT];
   int len = SxString::snprintf (buffer, MAX_N_CHARS_INT, i);
   replace (buffer, len);
}

SxString::SxString (long i)
  : SxArray<char> (), nChars(0), isUnicode_(false), isDirty(false)
{
   char buffer[MAX_N_CHARS_INT];
   int len = SxString::snprintf (buffer, MAX_N_CHARS_INT, i);
   replace (buffer, len);
}

SxString::SxString (unsigned long i)
  : SxArray<char> (), nChars(0), isUnicode_(false), isDirty(false)
{
   char buffer[MAX_N_CHARS_INT];
   int len = SxString::snprintf (buffer, MAX_N_CHARS_INT, i);
   replace (buffer, len);
}

SxString::SxString (long long i)
  : SxArray<char> (), nChars(0), isUnicode_(false), isDirty(false)
{
   char buffer[MAX_N_CHARS_INT];
   int len = SxString::snprintf (buffer, MAX_N_CHARS_INT, i);
   replace (buffer, len);
}

SxString::SxString (unsigned long long i)
  : SxArray<char> (), nChars(0), isUnicode_(false), isDirty(false)
{
   char buffer[MAX_N_CHARS_INT];
   int len = SxString::snprintf (buffer, MAX_N_CHARS_INT, i);
   replace (buffer, len);
}

SxString::SxString (int i, int width)
  : SxArray<char> (), nChars(0), isUnicode_(false), isDirty(false)
{
   SX_CHECK (width >= 0 && width <= MAX_N_CHARS_INT, width, MAX_N_CHARS_INT);
   // --- Writing the passed value to a temporary string variable
   char buffer[MAX_N_CHARS_INT + 1];
   int len = ::sprintf (buffer, "%0*d", width, i);
   SX_CHECK (len >= 0);
   // Converting the temporary string into a SxString-object and assigning it
   replace (buffer, len);
}

SxString::SxString (float i)
  : SxArray<char> (), nChars(0), isUnicode_(false), isDirty(false)
{
   // --- Writing the passed value to a temporary string variable
   char buffer[MAX_N_CHARS_FLOAT + 1];
   int len = ::sprintf (buffer, "%g", i);
   SX_CHECK (len >= 0);
   // Converting the temporary string into a SxString-object and assigning it
   replace (buffer, len);
}

SxString::SxString (double i)
  : SxArray<char> (), nChars(0), isUnicode_(false), isDirty(false)
{
   // --- Writing the passed value to a temporary string variable
   char buffer[MAX_N_CHARS_DOUBLE + 1];
   int len = ::sprintf (buffer, "%g", i);
   SX_CHECK (len >= 0);
   // Converting the temporary string into a SxString-object and assigning it
   replace (buffer, len);
}

SxString::SxString (double i, const SxString &format)
  : SxArray<char> (), nChars(0), isUnicode_(false), isDirty(false)
{
   // --- Writing the passed value to a temporary string variable
   char buffer[MAX_N_CHARS_DOUBLE + 1];
   int len;
   // Using "%g" as default format
   if (format == "")  {
      len = ::sprintf (buffer, "%g", i);
   } else  {
      len = ::sprintf (buffer, format.getElems (), i);
   }
   SX_CHECK (len >= 0);
   // Converting the temporary string into a SxString-object and assigning it
   replace (buffer, len);
}

SxString::SxString (const SxString &in)
  : SxArray<char> (), nChars(0), isUnicode_(false), isDirty(false)
{
   replace (in);
   if (in.isUnicode ()) isUnicode_ = true;
}

SxString::SxString (SxString &&in)
  : SxArray<char> (), nChars(0), isUnicode_(false), isDirty(false)
{
   replace (std::move(in));
}

SxString::SxString (const SxArray<char> &in)
  : SxArray<char> (), nChars(0), isUnicode_(false), isDirty(false)
{
   if (in.getSize () > 0)  {
      if (in.findPos ('\0') >= 0)  {
         replace (in.elements);
      }  else  {
         replace (in.elements, in.getSize ());
      }
   }
}

SxString::~SxString ()
{
   // empty
}

void SxString::updateNChars (ssize_t newNChars)
{
   SX_CHECK (newNChars >= 0, newNChars);
   nChars = newNChars;
   setDirty (false);
}

bool SxString::isEmpty () const
{
   return getSize () == 0;
}

void SxString::resize (ssize_t newNChars, bool keep)
{
   if (newNChars < 1)  {
      removeAll ();
      return; // done
   }
   ssize_t oldNChars = (keep ? getSize () : 0);
   if (newNChars < oldNChars)  { // merely shorten the string
      ssize_t newNBytes = SxConstChar::getNCharsToNBytes (getElems (),
         newNChars, isUnicode_);
      resize0 (newNBytes, keep);
   }  else if (newNChars > oldNChars)  { // extend the string with spaces
      ssize_t oldNBytes = (keep ? getNBytes () : 0);
      ssize_t count = newNChars - oldNChars;
      ssize_t newNBytes = oldNBytes + count;
      resize0 (newNBytes, keep);
      ssize_t nBlank = count;
      if (nBlank > 0)  {
         memset (elements + oldNBytes, ' ', static_cast<size_t>(nBlank));
      }
   }
   if (isUnicode_)  updateNChars (newNChars);
}

void SxString::resize0 (ssize_t nBytes, bool keep)
{
   if (nBytes > 0)  {
      // --- "abc" -> {'a', 'b', 'c', '\0'}
      SxArray<char>::resize (nBytes + 1, keep);
      if (SxArray<char>::getSize () > nBytes)  { // "should" be true
         SxArray<char>::operator()(nBytes) = 0;
         nChars = nBytes;
         if (isUnicode_)  setDirty (true); // callers must call updateNChars ()
         return; // done
      }
   }
   removeAll ();
}

void SxString::set (int c)
{
   ssize_t nChars_ = getSize ();
   if (nChars_ == 0)  return; // nothing to do
   if (c == 0)  c = '0';
   if (!isUnicode_)  {
      SX_CHECK (c > 0 && c <= 255, c);
      ::memset (elements, c, static_cast<size_t>(nChars_));
   }  else  {
      char buffer[4];
      ssize_t oldNBytes = getNBytes ();
      int nBufBytes = SxConstChar::encode (static_cast<SxConstChar::UnicodePoint>(c), true, buffer);
      ssize_t newNBytes = nBufBytes * nChars_;
      if (newNBytes != oldNBytes)  {
         resize0 (newNBytes, false);
         updateNChars (nChars_);
      }
      SxConstChar::set (elements, buffer, nBufBytes, nChars_);
   }
}

void SxString::replace (const char *in)
{
   isUnicode_ = false;
   if (in)  {
      size_t len = ::strlen (in);
      if (len > 0)  {
         replace (in, (ssize_t) len);
         return;
      }
   }

   removeAll ();
}

void SxString::replace (const char *in, ssize_t len)
{
   SX_CHECK (in);
   isUnicode_ = false;
   resize0 (len, false);
   if (len > 0)  {
      ::memcpy (elements, in, static_cast<size_t>(len));
   }
}

void SxString::replace (const SxString &in)
{
   replace (in.getElems (), in.getNBytes ());
   if (in.isUnicode_)  {
      isUnicode_ = true;
      updateNChars (in.getSize ());
   }
}

void SxString::replace (SxString &&in)
{
   SX_CHECK (in.getElems ());
   ssize_t _nBytes = in.getNBytes ();
   ssize_t _size = in.getSize ();
   //removeAll ();

   isUnicode_ = false;

   // move the internal array including terminating '\0'.
   SxArray<char>::operator= (std::move(in));
   nChars = _nBytes;
   if (in.isUnicode_)  {
      isUnicode_ = true;
      updateNChars (_size);
   }
   in.nChars = 0;
   in.elements = NULL;
}

#ifdef SX_IOS

// https://developer.apple.com/library/mac/documentation/CoreFoundation/Conceptual/CFStrings/Articles/AccessingContents.html
//
//     CFStringGetCString() did not work
//     str is expected in UTF-16 encoding
//
// - <https://developer.apple.com/library/mac/documentation/CoreFoundation/Conceptual/CFStrings/Articles/UnicodeBasis.html>: "a CFString object represents an
//   array of 16-bit Unicode characters (UniChar)" which may contain surrogates
// - <http://www.opensource.apple.com/source/CarbonHeaders/CarbonHeaders-18.1/MacTypes.h>: "typedef UInt16 UniChar;", thus unsigned
// - <https://developer.apple.com/library/mac/documentation/CoreFoundation/Reference/CFStringRef/#//apple_ref/c/func/CFStringGetLength>: "Returns the number
//   (in terms of UTF-16 code pairs) of Unicode characters"; that's bogus
//   because two surrogate UTF-16 values form one Unicode character; the word
//   "pair" can refer to the fact that each UTF-16 value consists of two bytes
//   or to surrogate pairs. We shouldn't use that function, but we seem to have
//   no better option...

void SxString::replace (CFStringRef str)
{
   removeAll ();
   isUnicode_ = false;
   if (!str)  return;

   CFRange range;
   range.location = 0;
   range.length = CFStringGetLength(str);
   ssize_t length = static_cast<ssize_t>(range.length);
   if (length < 1)  return;
   SxArray<UniChar> uniChar(length);
   CFStringGetCharacters(str, range, uniChar.elements);

   // --- check whether Unicode is necessary
   const UniChar *origSrc = uniChar.elements, *src = origSrc;
   const UniChar *srcLimit = &(origSrc[length]);
   for (ssize_t i = 0; i < length; ++i)  {
      UniChar val = src[i];
      if ( (val > 255) && (val != 8217) )  {
         isUnicode_ = true;
         break;
      }
   }

   // --- calculate how many bytes are necessary at most
   ssize_t newNBytes = 0, nChars_ = 0;
   if (!isUnicode_)  newNBytes = length;
   else  {
      while (src < srcLimit)  {
         char dummyBuf[4];
         newNBytes += SxConstChar::utf16ToUtf8 (&src, dummyBuf);
      }
   }

   // --- store the bytes, converting to UTF-8 if necessary
   resize0 (newNBytes, false);
   char *origDest = elements, *dest = origDest;
   if (!isUnicode_)  {
      for (ssize_t i = 0; i < length; ++i)  {
         UniChar val = src[i];
         if (val == 8217)  val = '\'';
         *dest++ = (char) val;
      }
   }  else  {
      for (src = origSrc; src < srcLimit; nChars_++)  {
         dest += SxConstChar::utf16ToUtf8 (&src, dest);
      }
   }
   ssize_t actualNBytes = (ssize_t) (dest - origDest);
   if (actualNBytes < newNBytes)  resize0 (actualNBytes, true); // shrink
   if (isUnicode_)  updateNChars (nChars_);
}

#endif /* SX_IOS */

SxString SxString::substitute (const SxString &what,
                               const SxString &with_,
                               ssize_t        nTimes) const
{
   if (what.getSize() < 1 || nTimes == 0 || what == with_)  {
      return *this;
   }

   ssize_t n = containsWhole (what);
   if (n < 1)  {
      return *this;
   }
   if (nTimes >= 0 && n > nTimes)  {
      n = nTimes;
   }

   SxString res;
   bool needUni = isUnicode_ || with_.isUnicode_;
   if (needUni)  res.isUnicode_ = true;
   const SxString &input = ( (needUni && !isUnicode_) ? toUnicode () : *this );
   const SxString &with = ( (needUni && !with_.isUnicode_) ?
      with_.toUnicode () : with_ );

   ssize_t whatNBytes = what.getNBytes ();
   ssize_t whatNBytesTl = n * whatNBytes;
   ssize_t withNBytes = with.getNBytes ();
   ssize_t withNBytesTl = n * withNBytes;
   ssize_t inputNBytes = input.getNBytes ();
   ssize_t newNBytes = inputNBytes - whatNBytesTl + withNBytesTl;
   if (newNBytes < 1)  {
      return res;
   }
   res.resize0 (newNBytes, false);
   char *dest = res.elements;

   ssize_t whatNChars = what.getSize ();
   if (needUni)  {
      ssize_t withNChars = with.getSize ();
      ssize_t inputNChars = input.getSize ();
      res.updateNChars (inputNChars - n * whatNChars + n * withNChars);
   }

   SxConstChar str(input.elements, inputNBytes, needUni);
   SxConstChar::Iterator it = str.begin ();
   ssize_t prevByteIdx = 0;
   while ( (--n >= 0) && (it.contains (what.elements)) )  {
      ssize_t byteIdx = it.getByteIdx ();
      ssize_t keepNBytes = byteIdx - prevByteIdx;
      if (keepNBytes > 0)  { // copy characters before the matched 'what'
         ::memcpy (dest, input.elements + prevByteIdx, static_cast<size_t>(keepNBytes));
         dest += keepNBytes;
      }
      if (withNBytes > 0)  { // copy 'with'
         ::memcpy (dest, with.elements, static_cast<size_t>(withNBytes));
         dest += withNBytes;
      }
      it += whatNChars; // skip 'what' to prepare for the next iteration
      prevByteIdx = it.getByteIdx ();
   }
   ssize_t restNBytes = inputNBytes - prevByteIdx;
   if (restNBytes > 0)  {
      ::memcpy (dest, input.elements + prevByteIdx,
         static_cast<size_t>(restNBytes));
   }
   return res;
}

void SxString::concatenate (const char *str, const char *append_)
{
   isUnicode_ = false;
   if (str || append_)  {
      ssize_t n1 = (str) ? static_cast<ssize_t>(strlen (str)) : 0;
      ssize_t n2 = (append_) ? static_cast<ssize_t>(strlen (append_)) : 0;
      if (n1 > 0 || n2 > 0)  {
         resize0 (n1 + n2, false);
         if (SxArray<char>::getSize() > n1 + n2)  {
            if (n1 > 0)  {
               memcpy (elements, str, static_cast<size_t>(n1));
            }
            if (n2 > 0)  {
               memcpy (elements + n1, append_, static_cast<size_t>(n2));
            }
         }
      }  else  {
         removeAll ();
      }
   }  else  {
      removeAll ();
   }
}

SxString SxString::asciiToUnicode (const char *src, ssize_t nBytes)
{
   SxString res;
   res.isUnicode_ = true;
   if (src == NULL)  return res;
   if (nBytes == -1)  nBytes = (ssize_t) ::strlen (src);
   SX_CHECK (nBytes >= 0, nBytes);
   if ( (nBytes == 0) || (*src == '\0') )  return res;
   ssize_t nUtf8Bytes = SxConstChar::getNUtf8Bytes (src, nBytes); // just count
   if (nUtf8Bytes < 1)  return res;
   res.resize0 (nUtf8Bytes, false);
   res.updateNChars (nBytes);
   SxConstChar::asciiToUtf8 (res.elements, src, nBytes); // convert and store
   return res;
}

SxString SxString::unicodeFromUtf8 (const char *src,
                                    ssize_t nBytes,
                                    ssize_t nChars_)
{
   SxString res;
   res.isUnicode_ = true;
   if (src == NULL)  return res;
   if (nBytes == -1)  nBytes = (ssize_t) ::strlen (src);
   SX_CHECK (nBytes >= 0, nBytes);
   if ( (nBytes == 0) || (*src == '\0') )  return res;
   res.resize0 (nBytes, false);
   ::memcpy (res.elements, src, static_cast<size_t>(nBytes));
   if (nChars_ == -1)  {
      nChars_ = SxConstChar::getNUtf8Chars (res.elements,
                                            res.elements + nBytes);
   }
   SX_CHECK (nChars_ <= nBytes, nChars_, nBytes); // plausibility check
   res.updateNChars (nChars_);
   // If number of character is equal to number of nBytes
   // then the string is ascii
   if (nChars_ == nBytes) res.isUnicode_ = false;
   return res;
}

SxString SxString::fromUtf8 (const char *src)
{
	return SxString::unicodeFromUtf8 (src);
}

SxString SxString::fromUtf16 (const uint16_t *src)
{
	return SxString(src);
}

SxString SxString::fromAscii(const char *src)
{
	return SxString(src);
}

void SxString::equalize (SxString *s1, SxString *s2)
{
   if (s1->isUnicode_ == s2->isUnicode_)  return; // nothing to do
   if (!s1->isUnicode_)  *s1 = s1->toUnicode ();
   else  *s2 = s2->toUnicode ();
}

void SxString::concatenate (const SxString &str, const SxString &append_)
{
   SxString s1 = str, s2 = append_;
   equalize (&s1, &s2);
   concatenate (s1.elements, s2.elements); // concatenate the raw bytes
   isUnicode_ = s1.isUnicode_;
   if (isUnicode_)  updateNChars (s1.getSize () + s2.getSize ());
}

void SxString::replace (const uint16_t *origStr)
{
   isUnicode_ = true;
   if ( (origStr == NULL) || (*origStr == '\0') )  {
      removeAll ();
      return; // done
   }
   SxArray<char> arr = SxConstChar::wcharsToUtf8 (origStr);
   const char *src = arr.elements;
   ssize_t newNBytes = (src ? ((ssize_t) ::strlen (src)) : 0);
   resize0 (newNBytes, false);
   if (newNBytes > 0)  ::memcpy (elements, src, static_cast<size_t>(newNBytes));
   updateNChars (SxConstChar::countChars (elements, elements + newNBytes,
      true));
}

SxArray<uint16_t> SxString::toWChars () const
{
   return SxConstChar::toWChars (getElems (), getNBytes (), isUnicode_);
}

SxArray<uint16_t> SxString::utf16() const
{
	return toWChars();
}

const char *SxString::utf8() const
{
	if (elements == NULL) return "";
	return elements;
}

SxString SxString::subString (ssize_t from) const
{
   return subString (from, getSize () - 1);
}

SxString SxString::subString (ssize_t from, ssize_t to) const
{
   if (getSize () < 1) // can't construct an SxConstChar object from NULL ptr
      return (isUnicode_ ? unicodeFromUtf8 (NULL) : SxString());
   SxConstChar str(elements, getNBytes (), isUnicode_);
   return SxString::subString (str, from, to);
}

SxString SxString::subString (const SxConstChar &str,
                              ssize_t fromCharIdx,
                              ssize_t toCharIdx)
{
   ssize_t nChars = str.getNChars ();
   if (fromCharIdx < 0)  fromCharIdx = 0;
   if (toCharIdx >= nChars)  toCharIdx = nChars - 1;
   ssize_t n = toCharIdx - fromCharIdx + 1;
   bool isUni = str.isUtf8 ();
   const char *src = str.getStart ();
   if (!isUni)  return SxString(src + fromCharIdx, n);
   SxString res;
   res.isUnicode_ = true;
   if (n < 1)  return res; // done
   SxConstChar::Iterator it_(&str), it = it_ + fromCharIdx, it2 = it + n;
   ssize_t fromByte = it.getByteIdx (), toByte = it2.getByteIdx () - 1;
   ssize_t newNBytes = toByte - fromByte + 1;
   return unicodeFromUtf8 (src + fromByte, newNBytes, n);
}

SxString SxString::subStringByBytes (ssize_t fromByteIdx, ssize_t toByteIdx)
   const
{
   SX_CHECK_3VARS (fromByteIdx <= toByteIdx && fromByteIdx >= 0 &&
      toByteIdx < getNBytes (), fromByteIdx, toByteIdx, getNBytes ());
   SxString res;
   const char *src = getElems ();
   if (isUnicode_)  {
      // must not start/end "inside" a character representation
      SX_CHECK (SxConstChar::isUtf8CharStart (src[fromByteIdx]));
      SX_CHECK (!SxConstChar::isUtf8Lead (src[toByteIdx]));
      SX_CHECK (!SxConstChar::isUtf8Continuation (src[toByteIdx + 1]));
      res.isUnicode_ = true;
   }
   ssize_t nBytes = toByteIdx - fromByteIdx + 1;
   if (nBytes > 0)  {
      res.resize0 (nBytes, false);
      char *dest = res.elements;
      ::memcpy (dest, src + fromByteIdx, static_cast<size_t>(nBytes));
      if (isUnicode_)
         res.updateNChars (SxConstChar::getNUtf8Chars (dest, dest + nBytes));
   }
   return res;
}

SxString SxString::head (ssize_t len) const
{
   return subString (0, len - 1);
}

SxString SxString::tail (ssize_t len) const
{
   return subString (getSize () - len, getSize () - 1);
}

SxString SxString::left (const SxString &in) const
{
   ssize_t pos = find (in);
   if (pos <= 0)  return SxString ();

   return subString (0, pos-1);
}

SxString SxString::right (const SxString &in) const
{
   const ssize_t &len = getSize();
   const ssize_t &lenIn = in.getSize();

   ssize_t pos = find (in);
   if (pos < 0)             return SxString ();
   if (pos + lenIn == len)  return SxString ();

   return subString (pos + lenIn);
}

SxConstChar::UnicodePoint SxString::charAt (ssize_t charIdx) const
{
   SX_CHECK (charIdx >= 0 && charIdx <= getSize (), charIdx, getSize ());
   if (charIdx == getSize ())  return 0; // end-of-string marker
   SxConstChar str(elements, getNBytes (), isUnicode_);
   SxConstChar::Iterator it = str.begin () + charIdx;
   return *it;
}

SxString SxString::operator() (ssize_t start_, ssize_t end_) const
{
   return subString (start_, end_);
}

ssize_t SxString::find (const SxString &needle_, ssize_t fromCharIdx) const
{
   if ( (fromCharIdx < 0) || (fromCharIdx >= getSize ()) )  return -1;
   if (needle_.getSize () == 0)  return -1; // can't find "nothing"
   SxConstChar haystack(elements, getNBytes (), isUnicode_);
   SxConstChar::Iterator it = haystack.begin () + fromCharIdx;
   SxConstChar needle(needle_.getElems (), needle_.getNBytes (),
      needle_.isUnicode ());
   if (it.contains (needle))  return it.getCharIdx ();
   return -1;
}

ssize_t SxString::find (const char *needle_, ssize_t fromCharIdx) const
{
   if ( (fromCharIdx < 0) || (fromCharIdx >= getSize ()) )  return -1;
   if ( (!needle_) || (!needle_[0]) )  return -1; // can't find "nothing"
   SxConstChar haystack(elements, getNBytes (), isUnicode_);
   SxConstChar::Iterator it = haystack.begin () + fromCharIdx;
   SxConstChar needle(needle_, needle_ + ::strlen (needle_), false);
   if (it.contains (needle))  return it.getCharIdx ();
   return -1;
}

ssize_t SxString::findLast (const SxString &needle_, ssize_t toCharIdx) const
{
   if (needle_.getSize () == 0)  return -1; // can't find "nothing"
   ssize_t res = -1;
   
   SxConstChar temp (elements, getNBytes (), isUnicode_);
   SxConstChar::Iterator tempIt = temp.end ();
   if (toCharIdx >= 0 && temp.getNChars () > toCharIdx + 1)
      tempIt = temp.begin () + toCharIdx;

   SxConstChar haystack (elements, tempIt.getByteIdx (), isUnicode_);
   SxConstChar needle (needle_.getElems (),
                       needle_.getNBytes (),
                       needle_.isUnicode ());

   SxConstChar::Iterator it = haystack.end ();  
   if (toCharIdx >= 0 && temp.getNChars () > toCharIdx + 1)
      it = haystack.begin () + toCharIdx;

   for (;; it--) {
      if (it.contains (needle)) {
         res = it.getCharIdx ();
         break;
      }
      if (it.getCharIdx () == 0) break;
   }
   return res;
}

ssize_t SxString::findLast (const char *needle_, ssize_t toCharIdx) const
{
   if ( (!needle_) || (!needle_[0]) )  return -1; // can't find "nothing"
   ssize_t res = -1;

   SxConstChar temp (elements, getNBytes (), isUnicode_);
   SxConstChar::Iterator tempIt = temp.end ();
   if (toCharIdx >= 0 && temp.getNChars () > toCharIdx + 1)
      tempIt = temp.begin () + toCharIdx;

   SxConstChar haystack (elements, tempIt.getByteIdx (), isUnicode_);
   SxConstChar needle (needle_, needle_ + ::strlen (needle_), false);

   SxConstChar::Iterator it = haystack.end ();  
   if (toCharIdx >= 0 && temp.getNChars () > toCharIdx + 1)
      it = haystack.begin () + toCharIdx;
   for (;; it--) {
      if (it.contains (needle)) {
         res = it.getCharIdx ();
         break;
      }
      if (it.getCharIdx () == 0) break;
   }
   return res;
}

SxList<ssize_t> SxString::findAll (const SxString &needle_) const
{
   SxList<ssize_t> hits;
   if (needle_.getSize () == 0)  return hits; // can't find "nothing"
   SxConstChar haystack(elements, getNBytes (), isUnicode_);
   SxConstChar::Iterator it = haystack.begin ();
   SxConstChar needle(needle_.getElems (), needle_.getNBytes (),
      needle_.isUnicode ());
   while (it.contains (needle))  {
      hits << it.getCharIdx ();
      ++it;
   }
   return hits;
}

int SxString::getNOccurrences (const SxConstChar &needle, ssize_t fromCharIdx,
                               bool whole, ssize_t needleNChars) const
{
   if (fromCharIdx >= getSize ())  return 0;
   if (fromCharIdx < 0)  fromCharIdx = 0;
   if ( (whole) && (needleNChars == -1) )  needleNChars = needle.getNChars ();
   if (whole)  SX_CHECK (needleNChars >= 0);
   if ( (needleNChars == 0) || (*(needle.getStart ()) == '\0') )
      return 0; // can't find "nothing"
   int res = 0;
   SxConstChar haystack(getElems (), getNBytes (), isUnicode_);
   SxConstChar::Iterator it = haystack.begin () + fromCharIdx;
   while (it.contains (needle))  {
      res++;
      it += (whole ? needleNChars : 1);
   }
   return res;
}

int SxString::contains (const SxString &needle_) const
{
   SxConstChar needle(needle_.getElems (), needle_.getNBytes (),
      needle_.isUnicode ());
   return getNOccurrences (needle, 0, false);
}

int SxString::contains (const char *needle_) const
{
   SX_CHECK (needle_);
   SxConstChar needle(needle_, needle_ + ::strlen (needle_), false);
   return getNOccurrences (needle, 0, false);
}

int SxString::containsWhole (const SxString &needle_, ssize_t fromCharIdx) const
{
   SxConstChar needle(needle_.getElems (), needle_.getNBytes (),
      needle_.isUnicode ());
   return getNOccurrences (needle, fromCharIdx, true, needle_.getSize ());
}

int SxString::containsWhole (const char *needle_, ssize_t fromCharIdx) const
{
   SX_CHECK (needle_);
   ssize_t needleNBytes = static_cast<ssize_t>(::strlen (needle_)); // equals number of characters
   SxConstChar needle(needle_, needle_ + needleNBytes, false);
   return getNOccurrences (needle, fromCharIdx, true, needleNBytes);
}

SxString SxString::adjustRight () const
{
   const ssize_t &len = getSize ();
   SX_CHECK (len >= 0, len);
   SX_CHECK (len > 0);
   // Fetching a string where the leading and tailing whitespaces were removed
   SxString body = trim ();
   // --- Computing the lengths of the strings
   const ssize_t &lenBody = body.getSize ();
   SX_CHECK (lenBody >= 0, lenBody);
   ssize_t leftSpaces = len - lenBody;
   // --- Creating a string and filling it
   ssize_t nBytes = body.getNBytes ();
   SxString res;
   res.resize0 (leftSpaces + nBytes, false);
   if (body.isUnicode_)  {
      res.isUnicode_ = true;
      res.updateNChars (len);
   }
   char *dest = res.elements;
   while (leftSpaces-- > 0)  *dest++ = ' ';
   if (nBytes > 0)  ::memcpy (dest, body.elements, static_cast<size_t>(nBytes));
   return res;
}

SxString SxString::changeCase (bool toLower_) const
{
   const ssize_t &nBytes = getNBytes ();
   SX_CHECK (nBytes >= 0, nBytes);
   SxString res;
   if (isUnicode_)  res.isUnicode_ = true;
   if (nBytes == 0)  return res; // done
   ssize_t nChars_ = getSize ();
   res.resize0 (nBytes, false);
   if (res.isUnicode_)  res.updateNChars (nChars_);
   char *dest = res.elements;
   if (!isUnicode_)  {
      const char *src = elements;
      char ch;
      while ( (ch = *src++) != '\0' )  {
         if (toLower_)  ch = SxConstChar::asciiToLower (ch);
         else  ch = SxConstChar::asciiToUpper (ch);
         *dest++ = ch;
      }
   }  else  {
      SxConstChar str(elements, nBytes, isUnicode_);
      SxConstChar::Iterator it = str.begin (), itEnd = str.end ();
      for (; it != itEnd; ++it)  {
         SxConstChar::UnicodePoint u;
         if (toLower_)  u = it.toLower ();
         else  u = it.toUpper ();
         dest += SxConstChar::encode (u, true, dest);
      }
   }
   return res;
}

SxString SxString::trim  () const
{
   return simplifyWhiteSpace ();
}

SxString SxString::stripWhiteSpace () const
{
   const ssize_t &nBytes = getNBytes ();
   SX_CHECK (nBytes >= 0, nBytes);
   if (!nBytes)  return (isUnicode_ ? asciiToUnicode (NULL) : SxString());

   // --- Running through the leading white spaces
   SxConstChar str(getElems (), nBytes, isUnicode_);
   SxConstChar::Iterator it1 = str.begin (), itEnd = str.end ();
   while ( (it1 != itEnd) && (it1.isBlank ()) )  ++it1;
   // --- Testing if the string exclusively consists of white spaces
   if (it1 == itEnd)  {
      return SxString ();
   }
   // --- Identifying the trailing white spaces
   SxConstChar::Iterator it2 = itEnd;
   do { --it2; } while (it2.isBlank ());
   return subString (it1.getCharIdx (), it2.getCharIdx ());
}

SxString SxString::simplifyWhiteSpace () const
{
   const ssize_t &nBytes = getNBytes ();
   SX_CHECK (nBytes >= 0, nBytes);
   SxString res;
   if (isUnicode_)  res.isUnicode_ = true;
   if (!nBytes)  return res;
   res.resize0 (nBytes, false);
   char *origDest = res.elements, *dest = origDest;
   SxConstChar str(getElems (), nBytes, isUnicode_);
   SxConstChar::Iterator it = str.begin (), itEnd = str.end ();
   bool writeSpace = false;
   ssize_t newNChars = 0;
   for (; it != itEnd; ++it)  {
      if (it.isBlank ())  {
         if (dest > origDest)  {
            // defer; don't copy a space at beginning or end of string
            writeSpace = true;
         }
         continue;
      }
      // --- copy the current character
      if (writeSpace)  { // must write one deferred space character first
         writeSpace = false;
         *dest++ = ' ';
         newNChars++;
      }
      it.encode (&dest);
      newNChars++;
   }
   ssize_t newNBytes = dest - origDest;
   if (newNBytes < nBytes)  res.resize0 (newNBytes, true); // shrink
   if (res.isUnicode_)  res.updateNChars (newNChars);
   return res;
}

SxString SxString::removeWhiteSpace () const
{
   const ssize_t &nBytes = getNBytes ();
   SX_CHECK (nBytes >= 0, nBytes);
   SxString res;
   if (isUnicode_)  res.isUnicode_ = true;
   if (!nBytes)  return res;
   res.resize0 (nBytes, false);
   char *origDest = res.elements, *dest = origDest;
   ssize_t newNChars = 0;
   SxConstChar str(elements, nBytes, isUnicode_);
   SxConstChar::Iterator it = str.begin (), itEnd = str.end ();
   for (; it != itEnd; ++it)  {
      if (!it.isBlank ())  {
         it.encode (&dest);
         newNChars++;
      }
   }
   ssize_t newNBytes = dest - origDest;
   if (newNBytes < nBytes)  res.resize0 (newNBytes, true); // shrink
   if (res.isUnicode_)  res.updateNChars (newNChars);
   return res;
}

SxString SxString::stripComments () const
{
   const ssize_t &nBytes = getNBytes ();
   SX_CHECK (nBytes >= 0, nBytes);
   SxString res;
   if (isUnicode_)  res.isUnicode_ = true;
   if (!nBytes)  return res;
   res.resize0 (nBytes);
   char *origDest = res.elements, *dest = origDest;
   bool insideComment = false;
   ssize_t newNChars = 0;
   SxConstChar str(elements, nBytes, isUnicode_);
   SxConstChar::Iterator it = str.begin (), itEnd = str.end ();
   for (; it != itEnd; ++it)  {
      if (insideComment)  {
         if (!it.isLineTerm ())  continue;
         insideComment = false;
      }
      else if (*it == '#')  {
         insideComment = true;
         continue;
      }
      it.encode (&dest);
      newNChars++;
   }
   ssize_t newNBytes = dest - origDest;
   if (newNBytes < nBytes)  res.resize0 (newNBytes, true); // shrink
   if (res.isUnicode_)  res.updateNChars (newNChars);
   return res;
}

SxString SxString::wrap  (const SxString &prefix_,
                          ssize_t        indent,
                          ssize_t        length_,
                          bool           firstOnly) const
{
   if (getSize() == 0)  return "";  // nothing to wrap;
   ssize_t length = 80;

#  ifndef WIN32
      if (length_ != 0)  length = length_;
      else  {
         // --- get number of columns in terminal
         struct winsize ws;
         if ( (ioctl (STDOUT_FILENO, TIOCGWINSZ, &ws) >= 0) && (ws.ws_col > 0) )
            length = ws.ws_col;
      }
#  endif

   SX_CHECK (length > 0, length);

   // --- <prefix><indent space>
   SxString prefix = prefix_;
   bool needUni = (isUnicode_ || prefix.isUnicode_);
   const SxString &input = ( (needUni && !isUnicode_) ? toUnicode () : *this );
   if (needUni && !prefix.isUnicode_)  prefix = prefix.toUnicode ();
   prefix.append (' ', indent);
   ssize_t prefixNBytes = prefix.getNBytes ();
   ssize_t prefixNChars = prefix.getSize ();
   SX_CHECK (length > prefixNChars, length, prefixNChars);
   length -= prefixNChars; // Calculating the available line length

   // --- input
   ssize_t inputNChars = input.getSize ();
   SxConstChar inp(input.elements, input.getNBytes (), input.isUnicode_);
   SxConstChar::Iterator it = inp.begin (), itEnd = inp.end ();

   // --- find total length of the output, same as wrap but without write
   ssize_t nBytesTl = prefixNBytes; // arbitrary characters
   if ( (needUni) && (firstOnly) )  prefixNBytes = prefixNChars; // spaces
   while (it != itEnd)  {
      ssize_t maxNChars = it.getCharIdx () + length;
      if (maxNChars > inputNChars) maxNChars = inputNChars;
      ssize_t lastSpace = -1;
      while (it.getCharIdx () < maxNChars)  {
         if (it.isNewline ())  break;
         if (*it == ' ')  lastSpace = it.getCharIdx ();
         nBytesTl += it.getNCharsToNBytes (1);
         ++it;
      }
      if (it != itEnd)  {
         if (it.isNewline (true))  {
            nBytesTl++;
         }  else if (*it == ' ')  {
            nBytesTl++;
            do { ++it; } while (it != itEnd && *it == ' '); // skipping blanks
         }  else  {
            if (lastSpace >= 0)  {
               nBytesTl -= it.getByteIdx ();
               it -= (it.getCharIdx () - (lastSpace + 1));
               nBytesTl += it.getByteIdx ();
            }
            nBytesTl += it.getNCharsToNBytes (1);
         }
      }
      if (it != itEnd) nBytesTl += prefixNBytes;
   }

   SxString res;
   res.resize0 (nBytesTl, false);
   char *dest = res.elements;
   it.resetPos ();

   // --- prefix
   const char *src = prefix.elements;
   if (needUni)  prefixNBytes = prefix.getNBytes (); // may have changed
   for (ssize_t i=0; i < prefixNBytes; i++)  *dest++ = *src++;
   if (firstOnly)  {
      prefix.set (' ');
      if (needUni)  prefixNBytes = prefix.getNBytes (); // may have changed
   }

   while (it != itEnd)  {
      // --- wrap one line to length characters (from pos to max in the input)
      ssize_t maxNChars = it.getCharIdx () + length;
      if (maxNChars > inputNChars) maxNChars = inputNChars;

      // --- copy text and find last break, which can be newline or space
      ssize_t lastSpace = -1;
      while (it.getCharIdx () < maxNChars)  {
         if (it.isNewline ())  break;
         if (*it == ' ')  lastSpace = it.getCharIdx ();
         it.encode (&dest);
         ++it;
      }

      // --- 3 cases of what is the first character after wrap length
      if (it != itEnd)  {
         if (it.isNewline (true))  {
            *dest++ = '\n';
         }  else if (*it == ' ')  {
            *dest++ = '\n';
            do { ++it; } while (it != itEnd && *it == ' '); // skipping blanks
         }  else  {
            if (lastSpace >= 0)  {
               // --- move the last split word to the next line
               dest -= it.getByteIdx ();
               it -= (it.getCharIdx () - (lastSpace + 1));
               dest += it.getByteIdx ();
            }
            *dest++ = '\n';
         }
      }

      // --- prefix
      if (it != itEnd) {
         src = prefix.elements;
         for (ssize_t i=0; i < prefixNBytes; i++)  *dest++ = *src++;
      }
   }
   return res;
}

int SxString::toInt (bool *error) const
{
   if (error)  *error = false;
   int value = 0;
   bool failed = false;
   value = SxString::toNumber (elements, &failed, value, 10);
   if (failed)  {
      if (error)  { *error = true; return 0; }
      SX_THROW ("Can't convert string '"+ head(30) +"' to int.");
   }
   return value;
}

int64_t SxString::toInt64 (bool *error) const
{
   if (error)  *error = false;
   int64_t value = 0;
   bool failed = false;
   value = SxString::toNumber (elements, &failed, value, 10);
   if (failed)  {
      if (error)  { *error = true; return 0; }
      SX_THROW ("Can't convert string '"+ head(30) +"' to long.");
   }
   return value;
}

float SxString::toFloat (bool *error) const
{
   if (error)  *error = false;
   const ssize_t &len = getSize ();
   SX_CHECK (len >= 0, len);
   // Checking whether there is something to convert at all
   if (!len)  {
      if (error)  { *error = true; return 0.; }
      SX_THROW ("Can't convert string empty string to float.");
   }
   float value = 0.;

   // --- Trying to convert the string to a float value
   SxArray<char> buffer(len);
   if (sscanf (getElems (), "%f %s", &value, buffer.elements) != 1)  {
#    ifdef MACOSX
        if (strlen(buffer.elements) != 1 || buffer(0) != '.')
#    endif
        {
           if (error)  { *error = true; return 0.; }
           SX_THROW ("Can't convert string '"+ head(30) +"' to float.");
        }
   }
   return value;
}

double SxString::toDouble (bool *error) const
{
   if (error)  *error = false;
   const ssize_t &len = getSize ();
   SX_CHECK (len >= 0, len);
   // Checking whether there is something to convert at all
   if (!len)  {
      if (error)  { *error = true; return 0; }
      SX_THROW ("Can't convert string empty string to double.");
   }
   double value = 0.;

   // --- Trying to convert the string to a double value
   SxArray<char> buffer(len);
   if (sscanf (getElems (), "%lf %s", &value, buffer.elements) != 1)  {
#    ifdef MACOSX
        if (strlen(buffer.elements) != 1 || buffer(0) != '.')
#    endif
       {
          if (error)  { *error = true; return 0; }
          SX_THROW ("Can't convert string '"+ head(30) +"' to double.");
       }
   }
   return value;
}

bool SxString::isInt () const
{
   bool error = true;
   SxString::toNumber<int>(elements, &error, 0, 10);
   return !error;
}

bool SxString::isInt64 () const
{
   bool error = true;
   SxString::toNumber<int64_t>(elements, &error, 0, 10);
   return !error;
}

bool SxString::isFloat () const
{
   const ssize_t &len = getSize ();
   SX_CHECK (len >= 0, len);
   if (!len)  return false;
   float value = 0.;
   SxArray<char> buffer(len);
   if (sscanf (getElems (), "%f %s", &value, buffer.elements) != 1)  return false;
   else  return true;
}

bool SxString::isDouble () const
{
   const ssize_t &len = getSize ();
   SX_CHECK (len >= 0, len);
   if (!len)  return 0;
   double value = 0.;
   SxArray<char> buffer(len);
   if (sscanf (getElems (), "%lf %s", &value, buffer.elements) != 1)  return 0;
   else                                                            return 1;
}

SxList<SxString> SxString::tokenize (const SxString &delimiters,
                                     bool allowEmpty) const
{
   SxList<SxString> res;
   const ssize_t &nBytes = getNBytes ();
   SX_CHECK (nBytes >= 0, nBytes);
   if (!nBytes)  return res;

   enum delimMode { SingleAscii, SeveralAscii, Unicode };
      // (The delimiters string consists of either one ASCII character or
      // several ASCII characters or some general, complicated Unicode stuff.)
   delimMode mode;
   SxArray<SxConstChar::UnicodePoint> arrDelimU;

   const char *delements = delimiters.elements;
   ssize_t delimNBytes = delimiters.getNBytes ();
   SxConstChar dStr(delements, delimNBytes, delimiters.isUnicode_);
   SxConstChar::Iterator dIt = dStr.begin (), dItEnd = dStr.end ();
   if ( (delimiters.isUnicode_) /* && (!SxConstChar::is7bit (delements)) */ )
      mode = Unicode;
   else if (delimNBytes > 1)  mode = SeveralAscii;
   else  mode = SingleAscii;
   SxConstChar input(elements, getNBytes (), isUnicode_);
   SxConstChar::Iterator it = input.begin (), itEnd = input.end ();
   ssize_t prevCharIdx = 0;

   // --- Looking for delimiters inside string
   for (; it != itEnd; ++it)  {
      switch (mode)  {
         case SingleAscii:
            if (*it != (uint8_t) delements[0])  continue;
            break;
         case SeveralAscii:
            if (::strchr(delements,
                static_cast<int>(*it)) == NULL)  continue;
            break;
         case Unicode:
            bool found = false;
            for (dIt.resetPos (); dIt != dItEnd; ++dIt)  {
               if (*it == *dIt)  {
                  found = true;
                  break;
               }
            }
            if (!found)  continue;
            break;
      }
      ssize_t delimCharIdx = it.getCharIdx ();
      if ( (prevCharIdx < delimCharIdx) || (allowEmpty) )
         res << subString (prevCharIdx, delimCharIdx - 1);
      prevCharIdx = delimCharIdx + 1; // skip the delimiter
   }

   // --- Testing whether there is a substring that still has to be appended
   // and doing so if this is the case
   ssize_t charIdx = it.getCharIdx ();
   if ( (prevCharIdx < charIdx) || (allowEmpty) )
      res << subString (prevCharIdx, charIdx - 1);
   // Returning the list of substrings
   return res;
}

SxString SxString::join (const SxList<SxString> &list_,
                         const SxString &delimiter)
{
   // These join () functions start with the assumption that no involved string
   // uses Unicode (because that's the most likely case and the function shall
   // be fast). If the assumption turns out to be wrong because one of the
   // strings uses Unicode, all strings must be converted to Unicode.
   SxString res;

   if (list_.getSize () > 0)  {
      SxList<SxString> list;
      ssize_t nBytes, nBytesTl = 0, nCharsTl = 0;
      SxList<SxString>::ConstIterator it;
      bool needUnicode = delimiter.isUnicode_;
      for (it = list_.begin (); it != list_.end (); ++it)  { // count the bytes
         if ((*it).isUnicode_)  {
            if (!needUnicode)  { // must convert former strings and fix count
               needUnicode = true;
               ssize_t idx, n = list.getSize ();
               for (idx = 0; idx < n; idx++)  {
                  SxString oldS = list(idx), newS = oldS.toUnicode ();
                  list(idx) = newS;
                  nBytesTl = nBytesTl + newS.getNBytes () - oldS.getNBytes ();
               }
            }
            list << *it;
            nBytes = (*it).getNBytes ();
         }  else if (needUnicode)  { // must convert the current string
            SxString s = (*it).toUnicode ();
            list << s;
            nBytes = s.getNBytes ();
         }  else  {
            list << (*it);
            nBytes = (*it).getNBytes ();
         }
         nBytesTl += nBytes;
      }
      SxString delim = (needUnicode && !delimiter.isUnicode_)
                       ? delimiter.toUnicode () : delimiter;
      nBytesTl += delim.getNBytes () * (list.getSize () - 1);
      if (nBytesTl > 0)  {
         res.resize0 (nBytesTl, false);
         if (needUnicode)  res.isUnicode_ = true;
         ssize_t pos = 0;
         it = list.begin();
         if ((*it).getNBytes () > 0)  { // 0
            ::memcpy (res.elements, (*it).elements,
                      static_cast<size_t>((*it).getNBytes ()));
            pos += (*it).getNBytes ();
            nCharsTl += (*it).getSize ();
         }
         for (++it; it != list.end(); ++it)  { // 1 to n
            if (delim.getNBytes () > 0)  {
               ::memcpy (res.elements+pos,delimiter.elements,
                         static_cast<size_t>(delim.getNBytes ()));
               pos += delim.getNBytes ();
            }
            if ((*it).getNBytes () > 0)  {
               ::memcpy (res.elements + pos, (*it).elements,
                         static_cast<size_t>((*it).getNBytes ()));
               pos += (*it).getNBytes ();
               nCharsTl += (*it).getSize ();
            }
         }
         nCharsTl += delim.getSize () * (list.getSize () - 1);
         if (res.isUnicode_)  res.updateNChars (nCharsTl);
      }
   }

   return res;
}

SxString SxString::join (const SxArray<SxString> &arr,
                         const SxString &delimiter)
{
   SxString res;
   ssize_t n = arr.getSize ();

   if (n > 0)  {
      SxList<SxString> list;
      ssize_t nBytes, nBytesTl = 0, nCharsTl = 0;
      bool needUnicode = delimiter.isUnicode_;
      for (ssize_t idx = 0; idx < n; idx++)  {
         SxString str = arr(idx);
         if (str.isUnicode_)  {
            if (!needUnicode)  { // must convert former strings and fix count
               needUnicode = true;
               ssize_t idx_, n_ = list.getSize ();
               for (idx_ = 0; idx_ < n_; idx_++)  {
                  SxString oldS = list(idx_), newS = oldS.toUnicode ();
                  list(idx_) = newS;
                  nBytesTl = nBytesTl + newS.getNBytes () - oldS.getNBytes ();
               }
            }
            list << str;
            nBytes = str.getNBytes ();
         }  else if (needUnicode)  { // must convert the current string
            SxString s = str.toUnicode ();
            list << s;
            nBytes = s.getNBytes ();
         }  else  {
            list << str;
            nBytes = str.getNBytes ();
         }
         nBytesTl += nBytes;
      }
      SxString delim = (needUnicode && !delimiter.isUnicode_)
                       ? delimiter.toUnicode () : delimiter;
      nBytesTl += delim.getNBytes () * (n - 1);
      if (nBytesTl > 0)  {
         res.resize0 (nBytesTl, false);
         if (needUnicode)  res.isUnicode_ = true;
         ssize_t pos = 0;
         SxList<SxString>::ConstIterator it = list.begin();
         if ((*it).getNBytes () > 0)  { // 0
            ::memcpy (res.elements, (*it).elements,
                      static_cast<size_t>((*it).getNBytes ()));
            pos += (*it).getNBytes ();
            nCharsTl += (*it).getSize ();
         }
         for (++it; it != list.end(); ++it)  { // 1 to n
            if (delim.getNBytes () > 0)  {
               ::memcpy (res.elements+pos,delimiter.elements,
                         static_cast<size_t>(delim.getNBytes ()));
               pos += delim.getNBytes ();
            }
            if ((*it).getNBytes () > 0)  {
               ::memcpy (res.elements + pos, (*it).elements,
                         static_cast<size_t>((*it).getNBytes ()));
               pos += (*it).getNBytes ();
               nCharsTl += (*it).getSize ();
            }
         }
         nCharsTl += delim.getSize () * (list.getSize () - 1);
         if (res.isUnicode_)  res.updateNChars (nCharsTl);
      }
   }

   return res;
}




SxString SxString::readStdin (const SxString &defVal)
{
   int size = 1;
   SxString str, newStr; str.resize (size);

   // --- Reading characters from standard in one-by-one
   int curChar;
   for (int i = 0; (curChar = getchar ()) != EOF && curChar != '\n'; ++i) {
      // Testing whether the currently allocated memory is enough
      if (i == size - 1)  {
         int newSize = size + 10240;
         newStr.resize (newSize, true);
         if (size > 0)  str = newStr;
         size = newSize;
      }
      // Inserting the read in character into the string
      str(i) = (char)curChar;
   }
   if (str.getSize () == 1)  return "";

   SxString res = str.trim ();
   if (res == "")  return defVal;
   return res;
}

ssize_t SxString::prepend (int c, ssize_t count)
{
   insert (0, c, count);
   return getSize ();
}

ssize_t SxString::prepend (const char *str)
{
   insert (0, str);
   return getSize ();
}

ssize_t SxString::prepend (const SxString &in)
{
   insert (0, in);
   return getSize ();
}

void SxString::insert (ssize_t newCharIdx, int c, ssize_t count)
{
   if (count > 0 && newCharIdx >= 0 && newCharIdx < getSize () + 1)  {
      if (c == 0)  c = '0';
      char buffer[4];
      int nBufBytes = SxConstChar::encode (static_cast<SxConstChar::UnicodePoint>(c), isUnicode_, buffer);
      ssize_t nInsBytes = count * nBufBytes, nChars_ = getSize ();
      ssize_t nBytes = getNBytes ();
      ssize_t newNBytes = nBytes + nInsBytes;
      resize0 (newNBytes, true);
      if (isUnicode_)  updateNChars (nChars_ + count);

      SxConstChar str(elements, nBytes, isUnicode_);
      SxConstChar::Iterator it = str.begin () + newCharIdx;
      ssize_t byteIdx = it.getByteIdx ();
      SxConstChar::makeGap (elements, byteIdx, nInsBytes, nBytes);
      SxConstChar::set (elements + byteIdx, buffer, nBufBytes, count); // fill
   }
}

void SxString::insert (ssize_t newCharIdx, const char *in)
{
   if (in && in[0] != '\0' && newCharIdx >= 0 && newCharIdx < getSize () + 1) {
      if (isUnicode_)  return insert (newCharIdx, SxString(in));
      ssize_t nBytes = getNBytes (), nChars_ = getSize ();
      ssize_t nInsBytes = static_cast<ssize_t>(::strlen (in));
      ssize_t newNBytes = nBytes + nInsBytes;
      resize0 (newNBytes, true);
      if (isUnicode_)  updateNChars (nChars_ + nInsBytes);
      if (SxArray<char>::getSize() > newNBytes)  {
         SxConstChar str(elements, nBytes, isUnicode_);
         SxConstChar::Iterator it = str.begin () + newCharIdx;
         ssize_t byteIdx = it.getByteIdx ();
         SxConstChar::makeGap (elements, byteIdx, nInsBytes, nBytes);
         ::memcpy (elements + byteIdx, in, static_cast<size_t>(nInsBytes));
      }
   }
}

void SxString::insert (ssize_t newCharIdx, const SxString &in)
{
   if (  in.getSize () == 0
      || newCharIdx < 0
      || newCharIdx > getSize ())
   {
      return; // nothing to do
   }
   if ( (!isUnicode_) && (!in.isUnicode_) )  {
      insert (newCharIdx, in.elements);
      return;
   }

   // --- Create a temporary copies to work with
   SxString s1 = SxString(*this);
   SxString s2 = SxString(in);
   // Make both of them unicode
   equalize (&s1, &s2);

   // Resize the current object to required byte size
   ssize_t n1 = s1.getNBytes ();
   ssize_t n2 = s2.getNBytes ();
   resize0 (n1 + n2, true);
   // Update the SxString
   updateNChars (s1.getSize () + s2.getSize ());

   if (!isUnicode_)  { // byte representation may have changed
      isUnicode_ = true;
      if (n1 > 0)  ::memcpy (elements, s1.elements, static_cast<size_t>(n1));
   }

   SxConstChar str (elements, n1, isUnicode_);
   SxConstChar::Iterator it = str.begin () + newCharIdx;
   ssize_t newBytePos = it.getByteIdx ();
   for (ssize_t i = n1 - 1; i >= newBytePos; i--)  {
      elements[i + n2] = elements[i];
   }
   ::memcpy (elements + newBytePos, s2.elements, static_cast<size_t>(n2));
}

ssize_t SxString::append (const int c, ssize_t count)
{
   insert (getSize (), c, count);
   return getSize ();
}

ssize_t SxString::append (const char *str)
{
   insert (getSize (), str);
   return getSize ();
}

ssize_t SxString::append (const SxString &in)
{
   insert (getSize (), in);
   return getSize ();
}

void SxString::remove (ssize_t charIdx, ssize_t count)
{
   ssize_t nChars_ = getSize ();
   if (count > 0 && charIdx >= 0 && charIdx < nChars_)  {
      if (charIdx + count > nChars_)  {
         count = nChars_ - charIdx;
      }
      ssize_t newNChars = nChars_ - count;
      if (newNChars > 0)  {
         ssize_t oldNBytes = getNBytes ();
         SxConstChar str(elements, oldNBytes, isUnicode_);
         SxConstChar::Iterator it = str.begin () + charIdx, it2 = it + count;
         ssize_t nRemovedBytes = it2.getByteIdx () - it.getByteIdx ();
         ssize_t newNBytes = oldNBytes - nRemovedBytes;
         for (ssize_t i = it.getByteIdx (); i < newNBytes; i++)
            elements[i] = elements[i + nRemovedBytes];
         resize0 (newNBytes, true);
         if (isUnicode_)  updateNChars (newNChars);
      }  else  {
         removeAll ();
      }
   }
}

void SxString::removeElement (int c)
{
   SxConstChar str(elements, getNBytes (), isUnicode_);
   SxConstChar::Iterator it = str.begin (), itEnd = str.end ();
   for (; it != itEnd; ++it)  {
      if (*it == (SxConstChar::UnicodePoint) c)  {
         remove (it.getCharIdx ());
         return; // done
      }
   }
}

void SxString::removeFirst ()
{
   remove (0);
}

void SxString::removeLast ()
{
   const ssize_t &nChars_ = getSize ();
   if (nChars > 0)  remove (nChars_ - 1, 1);
}

void SxString::removeAll ()
{
   SxArray<char>::resize (0);
   SX_CHECK (elements == NULL);
   nChars = 0;
   if (isUnicode_)  setDirty (false);
}

SxString SxString::sprintf (const char *fmt, ...)
{
   if (!fmt)  {
      return SxString ();
   }
   // --- http://linux.die.net/man/3/snprintf
   char tmp [maxFormatLength];
   int n = 0;

   va_list arg;
   va_start (arg, fmt);
   n = ::vsnprintf (tmp, sizeof(tmp), fmt, arg);
   va_end (arg);

   if (n < 0) { SX_EXIT; }
   if (n < maxFormatLength) return SxString (tmp, n);

   // --- fall-back if maxFormatLength was not sufficient
   char *tmp2 = new char[static_cast<size_t>(n+1)];
   va_start (arg, fmt);
   int n2 = ::vsnprintf (tmp2, static_cast<size_t>(n + 1), fmt, arg);
   va_end (arg);
   if (n2 != n) { SX_EXIT; }
   SxString result(tmp2, n);
   delete [] tmp2;
   return result;
}


SX_EXPORT_UTIL SxList<SxString> __lics;

SxString &SxString::operator= (const SxString &in)
{
   if (this != &in)  {
      replace (in);
   }
   return *this;
}

SxString &SxString::operator= (SxString &&in)
{
   if (this != &in)  {
      replace (std::move(in));
   }
   return *this;
}

SxString &SxString::operator= (const SxArray<char> &in)
{
   replace (in);
   return *this;
}

SxString &SxString::operator= (const char *in)
{
   replace (in);
   return *this;
}

bool SxString::operator== (const SxString &in) const
{
   const ssize_t n1 = SxArray<char>::getSize ();
   const ssize_t n2 = in.SxArray<char>::getSize ();
   if (n1 != n2)  {
      return n1 + n2 == 1;
   }  else if (n1 == 0 && n2 == 0)  {
      return true;
   }  else  {
      if (::strncmp (elements, in.elements, (size_t) n1) != 0)  return false;
      if (isUnicode_ == in.isUnicode_)  return true;
      return SxConstChar::is7bit (elements);
      // (If both strings only contain 7-bit ASCII, the setting of isUnicode_
      // doesn't matter. Otherwise: the byte sequences in the strings are
      // identical, but the resulting characters/meanings/interpretations are
      // different.)
   }
}

bool SxString::operator== (const char *str) const
{
   if (elements && str)  {
      if (::strcmp (elements, str) != 0)  return false;
      if (!isUnicode_)  return true;
      return SxConstChar::is7bit (elements);
        // (cf. comment in "operator== (const SxString &in)")
   }  else if (str)  {
      return str[0] == '\0';
   }  else if (elements)  {
      return elements[0] == '\0';
   }  else  {
      return true;
   }
}

bool SxString::operator! () const
{
   return getSize () == 0;
}

bool SxString::operator!= (const SxString &in) const
{
   return !(*this == in);
}

bool SxString::operator!= (const char *str) const
{
   return !(*this == str);
}

bool SxString::operator< (const SxString &in) const
{
   if (elements && in.elements)  {
      return strcmp (elements, in.elements) < 0;
   }  else if (in.elements) {
      return in.elements[0] != '\0';
   }  else  {
      return false;
   }
}

bool SxString::operator< (const char *str) const
{
   if (elements && str)  {
      return strcmp (elements, str) < 0;
   }  else if (str) {
      return str[0] != '\0';
   }  else  {
      return false;
   }
}

bool SxString::operator> (const SxString &in) const
{
   if (elements && in.elements)  {
      return strcmp (elements, in.elements) > 0;
   }  else if (elements) {
      return elements[0] != '\0';
   }  else  {
      return false;
   }
}

bool SxString::operator> (const char *str) const
{
   if (elements && str)  {
      return strcmp (elements, str) > 0;
   }  else if (elements) {
      return elements[0] != '\0';
   }  else  {
      return false;
   }
}

SxString SxString::operator<< (const SxString &in) const
{
   return *this + in;
}

//------------------------------------------------------------------------------

SxString::Buffer::Buffer(ssize_t maxNChars, Mode mode) : buffer(NULL)
{
   SX_CHECK (maxNChars >= 0, maxNChars);
   SX_CHECK (mode == ASCII || mode == UTF8 || mode == UTF16, mode);
   mode_ = mode;
   allocate (maxNChars);
}

SxString::Buffer::Buffer(const SxString &str, Mode mode) : buffer(NULL)
{
   SX_CHECK (mode == ASCII || mode == UTF8 || mode == UTF16, mode);
   mode_ = mode;
   ssize_t nBytes = str.getNBytes (), size;
   if (nBytes < 1)  return; // nothing to do
   const char *src = str.getElems ();
   char *dest;
#ifdef WIN32
   SxArray<uint16_t> wArr;
   const uint16_t *wSrc;
   uint16_t *wDest;
   ssize_t count;
#endif
   switch (mode_)  {
      case ASCII:
         SX_CHECK ( (!str.isUnicode ()) || (SxConstChar::is7bit (src)) );
         size = nBytes + 1;
         resize (size);
         ::memcpy (buffer, src, static_cast<size_t>(size)); // copy including the trailing '\0'
         break;
      case UTF8:
         size = str.getSize () + 1;
         resize (size);
         if ( (str.isUnicode ()) || (SxConstChar::is7bit (src)) )  {
            // we only must copy the bytes
            ::memcpy (buffer, src, static_cast<size_t>(nBytes + 1)); // copy incl. the trailing '\0'
         }  else  { // we must convert
            dest = getBuffer ();
            dest += SxConstChar::asciiToUtf8 (dest, src, nBytes); // convert
            *dest = '\0';
         }
         break;
#ifdef WIN32
      case UTF16:
         // IMPLEMENTME for POSIX too?
         wArr = str.toWChars ();
         wSrc = wArr.elements;
         count = wArr.getSize ();
         resize (count + 1);
         wDest = (uint16_t *) buffer;
         for (; count > 0; count--)  *wDest++ = *wSrc++;
         *wDest = 0;
         break;
#endif
      default:
         SX_EXIT; // invalid mode
         break;
   }
}

SxString::Buffer::~Buffer()
{
   delete [] buffer;
}

void SxString::Buffer::resize (ssize_t newMaxNChars)
{
   SX_CHECK (newMaxNChars >= 0, newMaxNChars);
   if (maxNChars_ == newMaxNChars)  return; // nothing to do
   delete [] buffer;
   buffer = NULL;
   allocate (newMaxNChars);
}

void SxString::Buffer::allocate (ssize_t newMaxNChars)
{
   SX_CHECK (newMaxNChars >= 0, newMaxNChars);
   maxNChars_ = newMaxNChars;
   ssize_t maxPerChar = maxNBytesPerChar (mode_), nElems;
   if (maxPerChar == 4)  nElems = maxNChars_;
   else  {
      SX_CHECK (maxPerChar == 1, maxPerChar);
      nElems = (maxNChars_ + 3) / 4;
   }
   buffer = new uint32_t [static_cast<size_t>(nElems)];
}

uint8_t SxString::Buffer::maxNBytesPerChar (Mode mode)
{
   return ( (mode == ASCII) ? 1 : 4 );
   // - In UTF-8, one character may need 4 bytes.
   // - In UTF-16, one character may need a surrogate pair, which means two
   //   uint16_t values, which means 4 bytes.
}

//------------------------------------------------------------------------------

//SxString operator+ (const char *a, const SxString &b)
//{
//   return SxString(a) + b;
//}


SX_EXPORT_UTIL std::ostream& operator<< (std::ostream &s, const SxString &in)
{
   if (in.getSize () == 0)  return s;
   else                     return s << in.getElems();
}
#ifdef WIN32
   SX_EXPORT_UTIL std::wostream& operator<< (std::wostream &s, const SxString &in)
   {
      if (in.getSize () == 0)  return s;
      else                     return s << in.getElems ();
   }
#endif

////----------------------------------------------------------------------------
////  The C-like printf function is replaced by the sxprintf function.
////  This wrapper ensures that the C printf function cannot be called
////  anymore.
////
//int printf (const char *fmt, ...)
//{
//// fprintf (stdout, "The C-printf() function must not be used anymore.\n");
//// fprintf (stdout, "Call 'sxprintf()' (defined in SxString.h) instead. \n");
//// SX_EXIT;
//   va_list args;
//   va_start (args, fmt);
//   int res = vprintf (fmt, args);
//   va_end (args);
//   return res;
//}
#ifndef va_copy
#  ifdef WIN32
      // dangerous: this hack works only under windows!!!
#     define va_copy(a,b) ((a)=(b))
#  else
#     define va_copy __va_copy
#  endif
#endif

int sxprintf (const char *fmt, va_list arg)
{
   // FIXME: MPI output on master only
   SX_CHECK (fmt);

   const int maxFormatLength = 1024;
   char tmp [maxFormatLength];

   va_list arg2;
   va_copy (arg2, arg);

   int n = ::vsnprintf (tmp, maxFormatLength, fmt, arg);

   if (n >= maxFormatLength)  {
      // --- fall-back if maxFormatLength was not sufficient
      char *tmp2 = new char[static_cast<size_t>(n+1)];
      int n2 = ::vsnprintf (tmp2, static_cast<size_t>(n + 1), fmt, arg2);
      if (n2 != n) { SX_EXIT; }
#     ifdef SX_ANDROID
         __android_log_write (ANDROID_LOG_INFO, SX_LOG_ID, tmp2);
#     elif defined(SX_IOS)
         sxLogMsg (SX_LOG_ID, tmp2);
#     else
         cout << tmp2;
#     endif /* SX_ANDROID */
      delete []tmp2;
   } else if (n > 0)  {
#     ifdef SX_ANDROID
         __android_log_write(ANDROID_LOG_INFO, SX_LOG_ID, tmp);
#     elif defined(SX_IOS)
         sxLogMsg (SX_LOG_ID, tmp);
#     else
         cout << tmp;
#     endif /* SX_ANDROID */
   }

   va_end (arg2);

   return n;
}

int sxprintf (const char *fmt, ...)
{
   if (!fmt)  {
      return 0;
   }
   va_list arg;
   va_start (arg, fmt);
   int n = sxprintf (fmt, arg);
   va_end (arg);
   return n;
}

int sxfprintf (FILE *fp, const char *fmt, ...)
{
   SX_CHECK (fp);
   if (!fmt) return 0;
   va_list arg;
   va_start(arg, fmt);
   int n;
#  ifdef SX_ANDROID
      static const int InfoLvl = 3;
       n = __android_log_vprint (InfoLvl, "SX", fmt, arg);
#  else
      if (fp != stdout)  {
         n = vfprintf (fp, fmt, arg);
      } else {
         n = sxprintf (fmt, arg);
      }
      va_end (arg);
#  endif
   return n;
}

