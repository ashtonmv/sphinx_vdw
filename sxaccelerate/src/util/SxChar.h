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

#ifndef _SX_CHAR_H_
#define _SX_CHAR_H_

#include <SxTypeDefs.h>
#include <SxUtil.h>
#include <SxArray.h>

using namespace std;

/** \brief Character handling class.

    This class provides functions to handle text characters and character
    sequences which are encoded in Latin-1 (8-bit ASCII) or UTF-8 (Unicode).
    It is a low-level helper class for SxString.

    \ingroup Tools */

class SX_EXPORT_UTIL SxConstChar
{
   public:
      /** Unicode codepoint (U+0000..U+10FFFF), representing one character */
      typedef uint32_t UnicodePoint;

      /** \brief Character iterator class.

          This class provides functions to iterate back and forth over the
          individual characters of the referenced SxConstChar objects and to
          analyze and convert the characters.

          Every iterator object has a current character/byte index position
          relative to the start of the string (SxConstChar object).

          \ingroup Tools */
      class SX_EXPORT_UTIL Iterator
      {
         public:
            /** Creates a new iterator for the given string and sets the
                current position to the start of the string */
            Iterator (const SxConstChar *obj);
           ~Iterator ();

            /** Returns the current character index position */
            inline ssize_t getCharIdx () const { return charIdx; }
            /** Returns the current byte index position */
            inline ssize_t getByteIdx () const { return byteIdx; }
            /** Sets the current character index position to the given value
                and adapts the byte index position accordingly */
            void setCharIdx (ssize_t charIdx_);
            /** Sets the current byte index position to the given value and
                adapts the character index position accordingly. The given byte
                index must refer to the start of a character representation. */
            void setByteIdx (ssize_t byteIdx_);
            /** Returns a pointer to the SxConstChar object for which the
                iterator was created */
            inline const SxConstChar *getCharObj () const { return charObj; }

            /** Returns a pointer to the memory address which correlates to the
                current index position in the text */
            inline const char *getPtr () const
            { return charObj->start + byteIdx; }
            /** Resets the iterator to the start of the string (SxConstChar
                object), that is, to index position 0 */
            inline void resetPos () { byteIdx = charIdx = 0; }
            /** Returns how many bytes are used to encode the given number of
                characters, starting at the current index position */
            inline ssize_t getNCharsToNBytes (ssize_t nChars_) const
            {
               SX_CHECK_3VARS (charIdx + nChars_ <= charObj->getNChars (),
                  charIdx, nChars_, charObj->getNChars ());
               return SxConstChar::getNCharsToNBytes (getPtr (), nChars_,
                  charObj->isUnicode);
            }

            /** Returns whether "operator*" may be used */
            inline bool inRange () const
            { return (getPtr () < charObj->limit); }
            /** Returns the character which is at the current position */
            UnicodePoint operator* () const;
            /** @{Increments the character index position, that is, moves to the
                  next character */
            inline void operator++ () { inc (); }
            inline void operator++ (int) { inc (); }
            /** @} */
            /** @{Decrements the character index position, that is, moves to the
                  previous character */
            inline void operator-- () { dec (); }
            inline void operator-- (int) { dec (); }
            /** @} */
            /** Increases the character index position by the given number of
                characters */
            inline void operator+= (const ssize_t &nChars_) { inc (nChars_); }
            /** Decreases the character index position by the given number of
                characters */
            inline void operator-= (const ssize_t &nChars_) { dec (nChars_); }
            /** Returns whether both iterators refer to the same index position
                in the same string (SxConstChar object) */
            inline bool operator== (const Iterator &it) const
            {
               return (charObj == it.charObj) &&
                  (getByteIdx () == it.getByteIdx ());
            }
            /** Returns the opposite of "operator==" */
            inline bool operator!= (const Iterator &it) const
            { return !(operator== (it)); }

            /** Increases the current character index position, that is, moves
                to a subsequent character. If 'count' is 1, moves to the next
                character. */
            void inc (ssize_t count = 1);
            /** Decreases the current character index position, that is, moves
                back to a prior character */
            void dec (ssize_t count = 1);

            /** Tries to find the codepoint 'u'. Updates the iterator position
                if found. Returns whether found. */
            bool contains (const SxConstChar::UnicodePoint &u);
            /** Tries to find 'needle'. Updates the iterator position if found.
                Returns whether found. */
            bool contains (const SxConstChar &needle);
            /** Tries to find 'needle'. Updates the iterator position if found.
                Returns whether found. */
            bool contains (const char *needle_);

            /** Returns whether the current character is a "blank"; in ASCII,
                this means ' ' or '\t'. This function is a locale-independent
                generalization of the standard C function isblank (). */
            bool isBlank () const;
            /** Returns whether the current character is a "line terminator".
                In ASCII, this only means '\n'. In Unicode, it encompasses the
                character categories "Zl" and "Zp". */
            bool isLineTerm () const;
            /** Returns whether a newline marker is at the current position. If
                'skip' is true, the marker is skipped (the current position is
                increased). - A newline marker can be a single CR or LF byte or
                the combination CRLF (Carriage Return '\r', Line Feed '\n'). In
                Unicode strings, it can also be any character for which
                isLineTerm () returns 'true'. - Cf.
                <http://en.wikipedia.org/wiki/Newline> */
            bool isNewline (bool skip = false);
            /** Returns whether the current character is an ASCII digit
                character. This function is a locale-independent version of the
                standard C function isdigit (). */
            bool isDigit () const;

            /** Returns the character at the current index position, converted
                to lowercase if appropriate */
            inline UnicodePoint toLower () const
            { return toLowerTblLookup (**this); }
            /** Returns the character at the current index position, converted
                to uppercase if appropriate */
            inline UnicodePoint toUpper () const
            { return toUpperTblLookup (**this); }

            /** Encodes (writes/copies) the current character into the given
                buffer and increases the buffer pointer accordingly */
            void encode (char **dest_) const;

            /** Returns an iterator whose character index position is 'nChars_'
                characters greater than the position of the given iterator */
            friend Iterator operator+(const Iterator &, ssize_t nChars_);

            /** Returns an iterator whose character index position is 'nChars_'
                characters lesser than the position of the given iterator */
            friend Iterator operator-(const Iterator &, ssize_t nChars_);

         protected:
            /** The object to which the iterator refers */
            const SxConstChar *charObj;
            /** The current index position */
            ssize_t charIdx, byteIdx;
      };

      /** @{Creates a new SxConstChar object and initializes it */
      SxConstChar (const char *start_, const char *limit_, bool isUnicode_);
      SxConstChar (const char *start_, ssize_t nBytes, bool isUnicode_);
      /** @} */
     ~SxConstChar ();

      /** Returns a pointer to the memory address at which the byte
          representation of the string begins */
      inline const char *getStart () const { return start; }
      /** Returns the number of bytes which comprise the object text */
      inline ssize_t getNBytes () const { return limit - start; }
      /** Returns the number of characters which comprise the object text */
      ssize_t getNChars () const;
      /** Returns whether the text is a UTF-8 encoded Unicode text */
      inline bool isUtf8 () const { return isUnicode; }

      /** Returns an iterator whose index position is set to the beginning of
          the string */
      Iterator begin () const;
      /** Returns an iterator whose index position is set to the end of the
          string, that is, to the 'limit' byte directly after the text */
      Iterator end () const;

      /** @{Locale-independent versions of standard functions for
            character-class checks with ASCII characters */
      /** Returns whether the given character is a "blank" character, that is,
          a space ' ' or tabulator '\t' character */
      static inline bool isBlank (char ch) { return ch == ' ' || ch == '\t'; }
      /** Returns whether the given character is a digit */
      static inline bool isDigit (char ch) { return ch >= '0' && ch <= '9'; }
      /** Returns whether the given character is a lowercase letter */
      static inline bool isLower (char ch) { return ch >= 'a' && ch <= 'z'; }
      /** Returns whether the given character is an uppercase letter */
      static inline bool isUpper (char ch) { return ch >= 'A' && ch <= 'Z'; }
      /** Returns whether the given character is a letter */
      static inline bool isAlpha (char ch)
      { return isUpper (ch) || isLower (ch); }
      /** Returns true iff the given character is a letter or digit */
      static inline bool isAlnum (char ch)
      { return isAlpha (ch) || isDigit (ch); }
      /** @} */

      /** Converts uppercase 7-bit ASCII letters to the respective lowercase
          letters. This function is a locale-independent replacement for the
          standard C function tolower (). */
      static inline char asciiToLower (char ch)
      { return SxConstChar::isUpper (ch) ? (char) (ch + ('a' - 'A')) : ch; }

      /** Converts lowercase 7-bit ASCII letters to the respective uppercase
          letters. This function is a locale-independent replacement for the
          standard C function toupper (). */
      static inline char asciiToUpper (char ch)
      { return SxConstChar::isLower (ch) ? (char) (ch - ('a' - 'A')) : ch; }

      /** Returns whether the given byte looks like a 7-bit ASCII character */
      static inline bool is7bit (char x)
      {
         return (!(x & 128));
      }

      /** Returns whether the given byte represents a single-byte UTF-8
          character, that is, a character whose entire representation fits into
          one byte. Such bytes represent exactly the 7-bit ASCII characters. */
      static inline bool isUtf8SingleByteChar (char x)
      {
         return is7bit (x);
      }

      /** Returns whether the given byte looks like a leading byte of a
          multi-byte UTF-8 character representation */
      static inline bool isUtf8Lead (char x)
      {
         return ((x & (128 | 64)) == (128 | 64));
      }

      /** Returns whether the given byte looks like the start of a (single-byte
          or multi-byte) UTF-8 character representation */
      static inline bool isUtf8CharStart (char x)
      {
         return isUtf8SingleByteChar (x) || isUtf8Lead (x);
      }

      /** Returns whether the given byte looks like a continuation (a
          non-leading byte) of a multi-byte UTF-8 character representation */
      static inline bool isUtf8Continuation (char x)
      {
         return ((x & (128 | 64)) == 128);
      }

      /** Returns whether the string only contains 7-bit ASCII characters */
      static bool is7bit (const char *str);
      /** Returns 0 for 7-bit ASCII texts, 1 for texts which only contain UTF-8
          bytes and 2 for general texts which contain non-UTF-8 bytes (likely
          8-bit ASCII, Latin-1) */
      static uint8_t analyzeText (const char *str, ssize_t nBytes = -1);

      /** Converts the UTF-8 byte sequence at the given address to one Unicode
          character and increases the pointer accordingly. */
      static SxConstChar::UnicodePoint utf8ToChar (const char **ptr_);
      /** Converts a single Unicode character to a UTF-8 byte sequence and
          writes the bytes to 'dest'. Returns the number of bytes which
          were written. */
      static int charToUtf8 (UnicodePoint u, char *dest);
      /** Converts a UTF-16 value (or a surrogate pair) to a sequence of UTF-8
          bytes and stores them in the 'buffer'. Increases the source pointer
          accordingly. Returns how many bytes were stored. */
      static int utf16ToUtf8 (const uint16_t **origSrc, char *buffer);

      /** Progresses to the start of a subsequent UTF-8 encoded character */
      static void utf8Inc (const char **ptr_, ssize_t count = 1);
      /** Moves back to the start of a prior UTF-8 encoded character
          representation */
      static void utf8Dec (const char **ptr_, ssize_t count = 1);
      /** Progresses to the start of a subsequent character. If 'count' is 1,
          moves to the start of the next character. */
      static void inc (const char **ptr_, ssize_t count, bool isUnicode_);
      /** Returns the number of UTF-8 encoded characters from the address
          'start_' to immediately before the address 'limit_'. If the given
          'limit_' is NULL, it is set to the address of the first found
          zero-byte. */
      static ssize_t getNUtf8Chars (const char *start_,
         const char *limit_ = NULL);
      /** Returns the number of bytes which are used to represent the 'nChars_'
          characters starting at address 'start_' */
      static ssize_t getNCharsToNBytes (const char *start_, ssize_t nChars_,
         bool isUnicode_);
      /** Returns the number of characters from the address 'start_' to
          immediately before the address 'limit_'. This function is roughly an
          opposite of getNCharsToNBytes (). */
      static ssize_t countChars (const char *start_, const char *limit_,
         bool isUnicode_);
      /** Converts a wchar_t (on Windows, this means UTF-16) string to a UTF-8
          encoded string. If the given string is NULL or empty, the result is
          an empty array. */
      static SxArray<char> wcharsToUtf8 (const uint16_t *origSrc,
         ssize_t nWchars = -1);
      /** Converts the given 8-bit ASCII string to a zero-terminated sequence
          of wchar_t values. If the given string is NULL or empty, the result
          is an empty array. */
      static SxArray<uint16_t> asciiToWChars (const char *src,
         ssize_t nBytes_ = -1);
      /** Converts the given (ASCII or UTF-8 encoded) string to a
          zero-terminated sequence of wchar_t (UTF-16) values. If the given
          string is NULL or empty, the result is an empty array. */
      static SxArray<uint16_t> toWChars (const char *src, ssize_t nBytes_,
         bool isUnicode_);
      /** @{Conversion of 8-bit ASCII strings to UTF-8 encoded strings */
      /** Returns how many bytes the UTF-8 representation of the given number
          of characters from the given 8-bit ASCII string would need */
      static ssize_t getNUtf8Bytes (const char *src, ssize_t nAsciiBytes);
      /** Converts the given 8-bit ASCII characters to UTF-8 encoded
          characters and stores them in the buffer 'dest'. Returns how many
          bytes it stored. */
      static ssize_t asciiToUtf8 (char *dest, const char *src,
         ssize_t nAsciiBytes);
      /** @} */
      /** Encodes the given character as a UTF-8 byte sequence or one ASCII
          byte. Writes the bytes into the provided buffer. Returns the number
          of written bytes. */
      static int encode (SxConstChar::UnicodePoint u, bool isUnicode_,
         char *buffer);
      /** Sets/fills the memory from address 'dest' with the given source
          byte sequence 'count' times. (The source bytes usually represent one
          character.) */
      static void set (char *dest, const char *src, ssize_t srcNBytes,
         ssize_t count);

      /** Makes a gap in a text by moving bytes "to the right". The gap can
          then e.g. be filled in an insert () method. */
      static void makeGap (char *str, ssize_t fromByteIdx, ssize_t nInsBytes,
         ssize_t oldNBytes);

      /** UTF-16 surrogate pair shift/mask values */
      enum { sxUtf16Shift = 10, sxUtf16Mask = ((1 << sxUtf16Shift) - 1) };

      /** Helper data type which contains a mapping (e.g. toLower/toUpper
          conversion) of a Unicode codepoint to another */
      typedef struct
      {
         SxConstChar::UnicodePoint from, to;
      } CaseConvItem;

   protected:
      /** The memory range of the underlying text bytes. We may access the
          bytes at start..(limit-1). */
      const char *start, *limit;
      /** The number of characters which are represented by the string
          (SxConstChar object). This is only calculated if necessary. Don't
          access this variable directly, use getNChars () instead! */
      mutable ssize_t nChars;
      /** Whether the value in the member variable 'nChars' is invalid */
      mutable bool isDirty;
      /** Whether the text is a UTF-8 encoded Unicode text */
      bool isUnicode;

      /** Initializes a newly constructed object */
      void initialize (const char *start_, const char *limit_, bool isUnicode_);

      /** Returns the difference between two Unicode codepoints */
      static int codePointDiffCB (const void *, const void *);
      /** Returns the difference between the 'from' values of two CaseConvItems
      */
      static int caseConvDiffCB (const void *, const void *);
      /** @{Performs a Unicode character lookup in a table which refers to
           case conversions (toLower/toUpper conversions) */
      static SxConstChar::UnicodePoint caseConvLookup (
         const CaseConvItem *tbl, size_t nElems,
         SxConstChar::UnicodePoint from);
      static UnicodePoint toLowerTblLookup (UnicodePoint from);
      static UnicodePoint toUpperTblLookup (UnicodePoint from);
      /** @} */

      /** Combines a leading and a trailing UTF-16 surrogate value and returns
          the resulting Unicode character */
      static SxConstChar::UnicodePoint combineUtf16 (const uint16_t &lead,
         const uint16_t &trail);

      /** Returns whether the given value looks like a leading surrogate for a
          UTF-16 surrogate pair */
      static inline bool isUtf16LeadVal (uint16_t x)
      {
         return (x >= 0xd800) && (x <= 0xdbff);
      }

      /** Returns whether the given value looks like a trailing surrogate for a
          UTF-16 surrogate pair */
      static inline bool isUtf16TrailVal (uint16_t x)
      {
         return (x >= 0xdc00) && (x <= 0xdfff);
      }
};

#endif /* _SX_CHAR_H_ */
