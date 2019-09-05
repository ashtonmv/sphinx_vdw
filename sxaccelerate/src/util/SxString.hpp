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

/** Template spezialization of the getNBytes template function.
    Used to evaluate the maximum memory consumption of a variable */
//template<>
//ssize_t getNBytes<SxString> (const SxString &in) {
//   return (in.getSize() * sizeof(char));
//}


/** \brief wrapper for standard printf function

    The C-like printf function is replaced by the sxprintf function.
    This wrapper ensures that the C printf function cannot be called
    anymore.
    \sa sxprintf */
//int printf (const char *fmt, ...);

inline ssize_t SxString::getSize () const
{
   SX_CHECK (!isDirty);
   return nChars;
}

inline const char *SxString::ascii () const
{
   SX_CHECK (!isUnicode_);
   return (elements ? elements : emptyString);
}

inline ssize_t SxString::getNBytes () const
{
   ssize_t res = SxArray<char>::getSize ();
   if (res > 0)  res--; // exclude the trailing '\0' byte
   return res;
}

inline const char *SxString::getElems () const
{
   return (elements ? elements : emptyString);
}

inline void SxString::setDirty (bool d)
{
#ifndef NDEBUG
   isDirty = d;
#endif
}

inline SxString SxString::toUnicode () const
{
   SX_CHECK (!isUnicode_);
   return SxString::asciiToUnicode (elements, getNBytes ());
}

inline SxString SxString::toUpper () const { return changeCase (false); }
inline SxString SxString::toLower () const { return changeCase (true); }

SxString SxString::operator+ (const SxString &in) const
{
   if (in.elements)  {
      return SxString (*this, in);
   }  else  {
      return *this;
   }
}

SxString SxString::operator+ (const char *str) const
{
   if (str)  {
      SxString res (elements, str);
      if (isUnicode ()) res.setUnicode ();
      res.updateNChars (nChars + (ssize_t)strlen (str));
      return res;
   }  else  {
      return *this;
   }
}

SxString SxString::operator+ (const char c) const
{
   return SxString (*this, SxString(c));
}

SxString SxString::operator+ (const int num) const {
   return SxString (*this, SxString(num));
}

SxString SxString::operator+ (const long num) const {
   return SxString (*this, SxString(num));
}

SxString SxString::operator+ (const unsigned long num) const {
   return SxString (*this, SxString(num));
}

SxString SxString::operator+ (const long long num) const {
   return SxString (*this, SxString(num));
}

SxString SxString::operator+ (const unsigned long long num) const {
   return SxString (*this, SxString(num));
}

SxString SxString::operator+ (const double num) const {
   return SxString (*this, SxString(num));
}

void SxString::operator+= (const SxString &in)
{
   append (in);
   if (in.isUnicode ()) isUnicode_ = true;
}

void SxString::operator+= (const char *in)
{
   append (in);
}

void SxString::operator+= (const char c)
{
   append (c);
}

void SxString::operator+= (const int num)
{
   append (SxString(num));
}

void SxString::operator+= (const long num)
{
   append (SxString(num));
}

void SxString::operator+= (const unsigned long num)
{
   append (SxString(num));
}

void SxString::operator+= (const long long num)
{
   append (SxString(num));
}

void SxString::operator+= (const unsigned long long num)
{
   append (SxString(num));
}

void SxString::operator+= (const double num)
{
   append (SxString(num));
}

inline SxString operator+ (const char c, const SxString &b)
{
   return SxString (SxString(c), b);
}

inline SxString operator+ (const int num, const SxString &b)
{
   return SxString (SxString(num), b);
}

inline SxString operator+ (const long num, const SxString &b)
{
   return SxString (SxString(num), b);
}

inline SxString operator+ (const unsigned long num, const SxString &b)
{
   return SxString (SxString(num), b);
}

inline SxString operator+ (const long long num, const SxString &b)
{
   return SxString (SxString(num), b);
}

inline SxString operator+ (const unsigned long long num, const SxString &b)
{
   return SxString (SxString(num), b);
}

inline SxString operator+ (const double num, const SxString &b)
{
   return SxString (SxString(num), b);
}

inline SxString operator+ (const char *a, const SxString &b)
{
   return SxString (a, b.elements);
}

inline bool operator== (const char *a, const SxString &b)
{
   if (a && b.elements)  {
      return strcmp (a, b.elements) == 0;
   }  else if (a) {
      return a[0] == '\0';
   }  else if (b.elements) {
      return b.elements[0] == '\0';
   }  else  {
      return true;
   }
}

inline bool operator!= (const char *a, const SxString &b)
{
   return !(a == b);
}

inline bool operator< (const char *a, const SxString &b)
{
   if (a && b.elements)  {
      if (a[0] == '\0') return a[0] != b.elements[0];
      return strcmp (a, b.elements) < 0;
   }  else if (b.elements) {
      return b.elements[0] != '\0';
   }  else  {
      return false;
   }
}

inline bool operator> (const char *a, const SxString &b)
{
   if (a && b.elements)  {
      if (a[0] == '\0') return false;
      return strcmp (a, b.elements) > 0;
   }  else if (a) {
      return a[0] != '\0';
   }  else  {
      return false;
   }
}


template<class T>
inline T SxString::toNumber (bool *error, int base_) const
{
   SX_CHECK (error != NULL);
   *error = true;
   T default_ = 0;
   T result = SxString::toNumber (elements, error, default_, base_);

   return result;
}

template<class T>
inline T SxString::toNumber (const char *str,
                             bool       *error_,
                             T          default_,
                             int        base_)
{
   SX_CHECK (error_ != NULL);
   *error_ = true;

   if (!str || str[0] == '\0')  {
      return default_;
   }

   ssize_t i = 0;

   while (SxConstChar::isBlank (str[i]))  {
      i++;
   }

   // --- sign
   bool neg = false;
   if (str[i] == '+') i++;
   else if (str[i] == '-')  { neg = true; i++; }

   // --- base
   unsigned long long base = 10;
   if (base_ > 1 && base_ <= 36)  {
      base = (unsigned long long)base_;
   }  else  {
      if (str[i] == '0')  {
         base = 8;
         i++;
      }  else if (str[i] == 'b' || str[i] == 'B')  {
         base = 2;
         i++;
      }
      if (str[i] == 'x' || str[i] == 'X') {
         if (base == 8)  {
            base = 16;
         }  else  {
            return default_;
         }
         i++;
      }
      if (base == 8) i--;
   }

   // --- last valid digit
   unsigned long long cutOff;
   if (neg) {
#     ifdef WIN32
#        pragma warning (disable:4146)
#     endif
      cutOff = (unsigned long long)(-(unsigned long long)sxmin<T>());
   }  else {
      cutOff = (unsigned long long)(sxmax<T>());
   }
   int cutLim = (int)(cutOff % base);
   cutOff /= base;

   // --- read number
   int c = 0;
   ssize_t start = i;
   unsigned long long acc = 0;
   while ((c = str[i]) != 0)  {
      if      (SxConstChar::isDigit ((char) c))  c -= '0';
      else if (SxConstChar::isLower ((char) c))  c = c - 'a' + 10;
      else if (SxConstChar::isUpper ((char) c))  c = c - 'A' + 10;
      else                   break;

      if (c >= (int)base)  break;

      if (acc > cutOff || (acc == cutOff && c > cutLim))  {
         // --- overflow
         break;
      }  else  {
         acc *= base;
         acc += (unsigned long long)c;
      }
      i++;
   }

   if (i > start) {
      while (SxConstChar::isBlank (str[i]))  {
         i++;
      }
      if (str[i] == '\0')  {
         *error_ = false;
         if (neg) return (T) (-(T)acc);
         return (T)acc;
      }
   }

   return default_;
}

template<class T>
inline int SxString::snprintfu (char *tmp, size_t n, T value)
{
   SX_CHECK (tmp);
   if (n < 1)  return 0;
   n--;

   size_t i = 0;
   if (value == 0)  {
      if (i < n) tmp[i++] = '0';
   }  else  {
      while (value > 0 && i < n)  {
         tmp[i++] = static_cast<char>((value % 10) + 48);
         value /= 10;
      }
      char c;
      for (size_t j = 0; j < i/2; j++) {
         c = tmp[j];
         tmp[j] = tmp[i - j - 1];
         tmp[i - j - 1] = c;
      }
   }
   tmp[i] = '\0';
   return (int)i;
}

template<class T>
inline int SxString::snprintfs (char *tmp, size_t n, T value)
{
   SX_CHECK (tmp);
   if (n < 1)  return 0;
   n--;

   size_t i = 0;
   if (value == 0)  {
      if (i < n) tmp[i++] = '0';
   }  else  {
      if (value < 0)  {
         while (value != 0 && i < n)  {
            tmp[i++] = static_cast<char>(-(value % 10) + 48);
            value /= 10;
         }
         if (i < n) tmp[i++] = '-';
      }  else  {
         while (value != 0 && i < n)  {
            tmp[i++] = static_cast<char>((value % 10) + 48);
            value /= 10;
         }
      }
      char c;
      for (size_t j = 0; j < i/2; j++) {
         c = tmp[j];
         tmp[j] = tmp[i - j - 1];
         tmp[i - j - 1] = c;
      }
   }
   tmp[i] = '\0';
   return (int)i;
}

template<class T>
inline int SxString::snprintf (char *tmp, size_t n, T value)
{
   return snprintfs (tmp, n, value);
}
template<> inline int SxString::snprintf<unsigned char> (
   char *tmp, size_t n, unsigned char value)
{
   return snprintfu (tmp, n, value);
}
template<> inline int SxString::snprintf<unsigned int> (
   char *tmp, size_t n, unsigned int value)
{
   return snprintfu (tmp, n, value);
}
template<> inline int SxString::snprintf<unsigned long> (
   char *tmp, size_t n, unsigned long value)
{
   return snprintfu (tmp, n, value);
}
template<> inline int SxString::snprintf<unsigned long long> (
   char *tmp, size_t n, unsigned long long value)
{
   return snprintfu (tmp, n, value);
}


// ---------------------------------------------------------------------------

inline char *SxString::Buffer::getBuffer () const
{
   return (char *) buffer;
}

inline ssize_t SxString::Buffer::getNBytes () const
{
   return maxNChars_ * maxNBytesPerChar (mode_);
}

inline ssize_t SxString::Buffer::getMaxNChars () const
{
   return maxNChars_;
}

inline SxString::Buffer::Mode SxString::Buffer::getMode () const
{
   return mode_;
}

