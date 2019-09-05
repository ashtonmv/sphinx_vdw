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

#ifndef _SX_STRING_H_
#define _SX_STRING_H_

// --- tell the compiler to check format string against types
#ifdef __GNUC__
#define __SXCHECK_FORMAT(fa,va) __attribute__((format(printf,fa,va)))
#else
#define __SXCHECK_FORMAT(fa,va)
#endif

#include <SxTypeDefs.h>
#include <SxChar.h>
#include <SxUtil.h>

using namespace std;

/** \brief String manipulation class.

    This class contains functions to handle strings.

    \ingroup Tools
    \author Sixten Boeck, boeck@gemmantics.com
    \author Vaclav Bubnik, bubnik@gemmantics.com */
class SX_EXPORT_UTIL SxString : public SxArray<char>
{
   public:

      /** The buffer object allows to store text temporarily, e.g. for system
          calls which require certain encodings (such as the UTF-16 encoding
          for wchar_t on Windows). The buffer text must be zero-terminated. The
          buffer text can be converted to an SxString object with the
          respective SxString constructor.

          Example usage: on Windows, wchar_t means UTF-16 encoding. To reserve
          and use a buffer which can contain up to 100 characters (including
          the zero-termination):
\code
          SxString::Buffer b(100, Buffer::UTF16);
          windows_wchar_syscall (b.getBuffer (), b.getNBytes ());
\endcode

          If you don't know/care whether the preprocessor symbol UNICODE is
          defined during compilation, you don't have to set an explicit mode in
          the constructor:
\code
          SxString::Buffer b(100);
          windows_syscall (b.getBuffer (), b.getNBytes ());
\endcode

          This latter example is useful in cases where the compiler
          automatically chooses a fooA() or fooW() function (ASCII/wchar_t)
          depending on the preprocessor symbol UNICODE on Windows.

          \ingroup Tools
      */
      class SX_EXPORT_UTIL Buffer
      {
         public:

            /** The mode or encoding of buffer text characters. A text may use
                either 8-bit ASCII or UTF-8 or UTF-16. On Windows, the data
                type wchar_t means UTF-16 encoding. */
            enum Mode { ASCII, UTF8, UTF16 };
            /** The default mode, depending on preprocessor symbols */
#ifdef WIN32 // Windows
#ifdef UNICODE
            static const Mode dfltMode = UTF16; // wchar_t on Windows
#else
            static const Mode dfltMode = ASCII;
#endif
#else // POSIX
            static const Mode dfltMode = ASCII;
#endif

            /** Returns a newly created buffer object. The internal text buffer
                may then be filled with up to 'maxNChars' characters (including
                the zero-termination). */
            Buffer(ssize_t maxNChars, Mode mode = dfltMode);
            /** Converts the given string object to a buffer object which uses
                the given mode. Returns the buffer object. */
            Buffer(const SxString &str, Mode mode = dfltMode);
            /** Destructs the buffer object */
           ~Buffer();

            /** Returns a pointer to the beginning of the internal buffer
                memory */
            inline char *getBuffer () const;
            /** Returns the number of bytes which were allocated for the
                internal buffer */
            inline ssize_t getNBytes () const;
            /** Returns the maximum number of characters (including the
                zero-termination) which may be stored in the internal buffer */
            inline ssize_t getMaxNChars () const;
            /** Returns the mode which was given in the constructor call */
            inline Mode getMode () const;

            /** Returns the maximum number of bytes which may be necessary to
                represent one character internally */
            static uint8_t maxNBytesPerChar (Mode mode = dfltMode);

         protected:

            /** Re-allocates the internal buffer memory to adapt to the given
                new maximum number of chars. Throws away any old contents. */
            void resize (ssize_t newMaxNChars);
            /** Allocates the internal buffer */
            void allocate (ssize_t newMaxNChars);

            /** The maximum number of characters (including the
                zero-termination) which may be stored in the internal buffer */
            ssize_t maxNChars_;
            /** The mode/encoding of the buffer text */
            Mode mode_;
            /** The pointer to the beginning of the internal buffer for text
                characters. It uses a 32-bit type only to ensure proper
                alignment for all modes. */
            uint32_t *buffer;
      };

      /** Creates an empty string */
      SxString ();
      /** Creates a new string that contains a single character */
      SxString (const char);
      /** Creates a new string that is a copy of the provided C-like
          string. The supplied string must be \0-terminated!*/
      SxString (const char *);
      /** Same as the constructor without 'ssize_t'. Only a given number of
          characters will be copied into the new string object */
      SxString (const char * str, ssize_t nChars_);
      /** Creates a Unicode string object from the given zero-terminated wide
          string */
      SxString (const uint16_t *);

      SxString (const SxString &);
      SxString (const SxArray<char> &);
      SxString (SxString &&);

      /** Creates a string object from the contents of the given buffer.
          Handles Unicode automatically. */
      SxString (const Buffer &);
#  ifdef SX_IOS
      /** Creates a new string from Core Foundation string */
      SxString (CFStringRef str);
#  endif /* SX_IOS */
      /** Concatenate two \0-terminated strings */
      SxString (const char *, const char *);
      /** Concatenate two strings */
      SxString (const SxString &, const SxString &);
      /** Converts an integer to a string */
      SxString (int);
      /** Converts an unsigned integer to a string */
      SxString (unsigned int);
      /** Converts an long to a string */
      SxString (long);
     /** Converts an unsigned long to a string */
      SxString (unsigned long);
      /** Converts an integer to a string */
      SxString (long long);
      /** Converts an unsigned integer to a string */
      SxString (unsigned long long);
      /** Converts an integer to a string according to the printf-format
          statement "%0#d". The string is right aligned, the empty left
          space is filled with zeros. Width specifies the field width */
      SxString (int i, int width);
      /** Converts a float to a string acoording to the printf format
          statment "%g". */
      SxString (float);
      /** Converts a double to a string acoording to the printf format
          statment "%g". */
      SxString (double);
      /** Converts a double to a string acoording to the supplied printf format
          statment. */
      SxString (double, const SxString &format);

      /** Destroys a string */
      virtual ~SxString ();

       /** Returns the number of characters (not including the trailing '\0'
          value); this may be different from the number of bytes if the string
          may contain general UTF-8 encoded Unicode characters. If you need the
          number of bytes, use the method getNBytes () instead. */
      inline ssize_t getSize () const;

      /** Returns the number of bytes which are used to represent the string
          internally (not including the trailing '\0' byte) */
      inline ssize_t getNBytes () const;

      /** Converts the given 8-bit ASCII string to a Unicode string object and
          returns that */
      static SxString asciiToUnicode (const char *src = NULL,
                                      ssize_t nBytes = -1);
      /** Converts the current non-Unicode string object to a Unicode string
          object and returns that */
      inline SxString toUnicode () const;

      /** Returns a string object created from a Ascii string */
      static SxString fromAscii (const char *src);

      /** Returns a Unicode string object which is created from the given UTF-8
          encoded bytes */
      static SxString unicodeFromUtf8 (const char *src,
                                       ssize_t nBytes = -1,
                                       ssize_t nChars_ = -1);
      /** Returns a Unicode string object created from a UTF-8 encoded bytes.
          The input buffer is expected to be "0x00" terminated. */
      static SxString fromUtf8 (const char *src);

      /** Returns a Utf16 encoded Unicode string object from a UTF16 buffer.
          The input buffer is expected to be "0x00 0x00" terminated. */
      static SxString fromUtf16 (const uint16_t *src);

      /** Returns wide char string representation of SxString including the 
          "0x00 0x00" terminator. */
      SxArray<uint16_t> utf16 () const;

      /** Returns c-string representation of the UTF-8 encoded SxString. */
      const char *utf8 () const;

      bool isEmpty () const;

      /** Returns a pointer to the internal '\0'-terminated string contents;
          not allowed for Unicode string objects */
      inline const char *ascii() const;

      /** Returns a pointer to the internal '\0'-terminated representation of
          the string contents; this works for _any_ string objects and returns
          either the ASCII or the UTF-8 representation. If you are unsure which
          internal encoding is used, you can call the method isUnicode (). -
          This method is mainly intended for use with operating system calls
          and other low-level uses. */
      inline const char *getElems () const;

      /** Returns whether the string object represents a UTF-8 encoded Unicode
          string */
      inline bool isUnicode () const { return isUnicode_; }

      /** Force string object interpretation as UTF-8 encoded Unicode string */
      inline void setUnicode () { isUnicode_ = true; }

      /** \brief This function modifies the size of the string according to the
        passed parameters. In contrast to the resize ()-function of the base
        class it initializes all normally uninitialized character fields with
        blanks. Here the virtual attribute is omitted willingly as it would
        decrement the speed of the array class. Constructs like the following
        are therefore deprecated:
        \code
        SxString str;
        SxArray<char> *arrPtr;
        arrPtr = &str;
        arrPtr->resize (...,...);// This calls SxArray<T>::resize () instead of
                                 // SxString::resize () and thus leads to un-
                                 // initialized values!!!
        \endcode  */
      void resize (ssize_t newNChars, bool keep = false);

      void set (int c);

      /** Replace the string with a C-like string */
      void replace (const char *);
      void replace (const SxString &);

      // replace with rvalue reference
      void replace (SxString &&);

#  ifdef SX_IOS
      /** Replaces the string with a Core Foundation string. Automatically
          enables Unicode if necessary, otherwise disables Unicode. */
      void replace (CFStringRef str);
#  endif /* SX_IOS */

      /** \brief Returns the string with x replaced by y
          \param what The substring to be replaced
          \param with The replacement for 'what'
          \param nTimes all occurrences are replaced by default */
      SxString substitute (const SxString &what,
                           const SxString &with,
                           ssize_t        nTimes = -1) const;


      /** Returns the substring from character index position 'from' to the
          end. */
      SxString subString (ssize_t from) const;
      /** Returns the substring from character index position 'from' to 'to' */
      SxString subString (ssize_t from, ssize_t to) const;
      /** Returns the substring from character index position 'fromCharIdx' to
          'toCharIdx' */
      static SxString subString (const SxConstChar &str,
                                 ssize_t fromCharIdx,
                                 ssize_t toCharIdx);
      /** This low-level method uses byte index positions instead of character
          index positions to calculate the returned substring. */
      SxString subStringByBytes (ssize_t fromByteIdx, ssize_t toByteIdx) const;
      /** Returns substring with up to first len characters. */
      SxString head (ssize_t len) const;
      /** Returns substring with up to last len characters. */
      SxString tail (ssize_t len) const;
      /** Searches for a substring and returns all characters left
          from the found position.
          \verbatim
             SxString a("hello-world");
             SxString b = a.left ("ll");
          \endverbatim
          b contains "he"; If the substring was not found an empty
          string is returned.  */
      SxString left  (const SxString &) const;
      /** Searches for a substring and returns all characters right
          from the found position.
          \verbatim
             SxString a("hello-world");
             SxString b = a.right ("ll");
          \endverbatim
          b contains "o-world"; If the substring was not found an empty
          string is returned. */
      SxString right (const SxString &) const;

      // Registering the operator () from the base class
      using SxArray<char>::operator ();
      /** Returns a substring */
      SxString   operator() (ssize_t start, ssize_t end) const;
      /** Returns the character which is at the given character index position
          of the string object */
      SxConstChar::UnicodePoint charAt (ssize_t charIdx) const;

      /** Returns the character index position of the first occurrence of the
          needle object in the 'this' object. Handles Unicode automatically. */
      ssize_t find (const SxString &needle_, ssize_t fromCharIdx = 0) const;
      /** Returns the position of the first occurrence of the ASCII needle in
          the 'this' object. Handles Unicode automatically. */
      ssize_t find (const char *needle_, ssize_t fromCharIdx = 0) const;

      /** Returns the character index position of the last occurrence of the
          needle object in the 'this' object. Handles Unicode automatically. */
      ssize_t findLast (const SxString &needle_, ssize_t toCharIdx = -1) const;
      /** Returns the position of the last occurrence of the ASCII needle in
          the 'this' object. Handles Unicode automatically. */
      ssize_t findLast (const char *needle_, ssize_t toCharIdx = -1) const;

      /** Returns positions of all occurrences in string */
      SxList<ssize_t> findAll (const SxString &) const;
      /** Returns number of occurrences in string */
      int  contains (const SxString &) const;
      int  contains (const char *) const;
      int  containsWhole (const SxString &, ssize_t fromCharIdx = 0) const;
      int  containsWhole (const char *, ssize_t fromCharIdx = 0) const;

      /** Adjusts the string to the right */
      SxString adjustRight () const;
      /** Returns a string with all characters converted to
          capital letters */
      inline SxString toUpper () const;
      /** Returns a string with all characters converted to lowercase letters */
      inline SxString toLower () const;
      /** Returns a string with all whitespaces (' ' and '\t') removed. */
      SxString trim () const;
      /** \brief Returns a string with all whitespace removed from begin or end

          Whitespaces are empty space ' ' (ASCII: 32) or tabulator (ASCII: 9).
          \par Example
\code
   SxString string = "  This     text  contains  spaces    ";
   SxString s = string.stripWhiteSpace (); // s == "This   text  contains  spaces"
\endcode
      \sa simplifyWhiteSpace
      \sa stripWhiteSpace
      \sa removeWhiteSpace
      \deprecate Use removeWhiteSpace instead
       */
      SxString stripWhiteSpace () const;
      /** \brief Returns a string with all whitespace being simplified

          Whitespaces are empty space ' ' (ASCII: 32) or tabulator (ASCII: 9).
          This function returns a string with whitespaces removed from the
          beginning and the end of the
          string. Furthermore sequences of white spaces will be replaced with a
          single space.
          \par Example
\code
   SxString string = "  This     text  contains  spaces    ";
   SxString s = string.stripWhiteSpace (); // s == "This text contains spaces"
\endcode
      \sa simplifyWhiteSpace
      \sa removeWhiteSpace
       */
      SxString simplifyWhiteSpace () const;
      /** \brief remove all white space in a string

          Whitespaces are empty space ' ' (ASCII: 32) or tabulator (ASCII: 9).
          This function returns the string without any whitespaces.
          \par Example
\code
   SxString string = "  This     text  contains  spaces    ";
   SxString s = string.removeWhiteSpace (); // s == "Thistextcontainsspaces"
\endcode
      \sa simplifyWhiteSpace
      \sa stripWhiteSpace
       */
      SxString removeWhiteSpace () const;
            /** \brief Removes all comments from a string

          This function removes all comments from a string.
\code
abc
def # ghi
# jkl
mno
 # leading white space
\endcode
          would return
\code
abc
def
mno
    <-- 1 white space!
\endcode
*/
      SxString stripComments () const;
      /** Prepend prefix, insert space (alignment), wrap lines,
          extra-indent following lines with prefix size, so that text is
          aligned.
          Scheme:
\verbatim
<prefix><indent space><text text>  |
<    indent space    ><blabla>     |<- end of line
<    indent space    ><text text>  |
\endverbatim

          @param prefix Text to be written in the first line.
          @param indent Number of extra space between prefix and text.
          @param length Length of line.

          @note To get the indent before the prefix, i.e.
\verbatim
<indent space><prefix><text text>  |
<    indent space    ><blabla>     |<- end of line
<    indent space    ><text text>  |
\endverbatim
          try
\code
wrap ( prefix.wrap("",indent), 0 , length);
\endcode
          @author C. Freysoldt
          */
      SxString wrap (const SxString &prefix = "",
                     ssize_t        indent = 0,
                     ssize_t        length = 0,
                     bool           firstOnly = true) const;
      /** Tokenizes a string and returns a list of strings:
          \verbatim
             SxString a ("this is a test program");
             SxList<SxString> words = a.tokenize (' ');
          \endverbatim
          The list 'words' is {"this", "is", "a", "test", "program"}
          \sa join
       */
      SxList<SxString> tokenize (const SxString &delimiters,
                                 bool allowEmpty=false) const;

      /** \brief Join a list of string

          Use this function to create a single string from a list of
          string.

          \par Example:
\code
   SxList<SxString> words = SxList<SxString> () << "this" << "is" << "a"
                                                << "test" << "program";
   cout << SxString::join (words, ":");
   cout << SxString::join (words);
\endcode
   This example code prints the strings "this:is:a:test:program" and
   "thisisatestprogram".
   \sa tokenize
       */
      static SxString join (const SxList<SxString> &list,
                            const SxString &delimiter="");
      static SxString join (const SxArray<SxString> &arr,
                            const SxString &delimiter="");


      /** \brief Read from stdin
        */
      static SxString readStdin (const SxString &defVal="");

      /** \brief Convert printf-arguments to SxString
          Returns NULL string when it is not enough to store the result. */
      static const int maxFormatLength = 1024;
      static SxString sprintf (const char *fmt, ...) __SXCHECK_FORMAT (1,2);

      template<class T> static int snprintf (char *s, size_t n, T number);
      template<class T> static int snprintfu (char *s, size_t n, T number);
      template<class T> static int snprintfs (char *s, size_t n, T number);

      /** Conversion to integer types.
          \verbatim
             SxString a("18446744073709551615");
             unsigned long val = a.toNumber<unsigned long>();

             SxString b("-200");
             signed char s = b.toNumber<signed char>(); // exception
          \endverbatim
      */
      template<class T> T toNumber (bool *error, int base=0) const;

      /** Silent conversion to integer types.
          A default value is used in any error situation.

          \param error Optional variable to return error flag.

          - base 10     SxString::toNumber ("0010", 0, NULL, 10) == 10
          - binary      SxString::toNumber (" b110 ", 0) == 6
          - octal       SxString::toNumber ("  -060", 0) == -48
          - decimal     SxString::toNumber ("   +10", 0) == 10
          - hexadecimal SxString::toNumber ("0xffff", 0) == 65535
          */
      template<class T> static T toNumber (const char *str,
                                           bool       *error,
                                           T          default_=T(),
                                           int        base = 0);


      /** Converts a integer number to an int value */
      int        toInt   (bool *error = NULL) const;
      /** Converts a integer number to an int64_t value */
      int64_t    toInt64 (bool *error = NULL) const;
      /** Converts a integer number to an float value */
      float  toFloat  (bool *error = NULL) const;
      /** Converts a integer number to an double value */
      double toDouble (bool *error = NULL) const;

      /** @{
          \brief Check if the string contains a certain data type

          Strings can be converted into int, int64_t, float, or double
          values. This function returns true if the object contains a string
          which can be converted. */
      bool isInt () const;
      bool isInt64 () const;
      bool isFloat () const;
      bool isDouble () const;
      /** @} */

      // --- Insertion
      ssize_t prepend (int c, ssize_t count = 1);
      ssize_t prepend (const char *str);
      ssize_t prepend (const SxString &);
      void insert (ssize_t newCharIdx, int c, ssize_t count = 1);
      void insert (ssize_t newCharIdx, const char *in);
      void insert (ssize_t newCharIdx, const SxString &in);
      ssize_t append (int c, ssize_t count = 1);
      ssize_t append (const char *str);
      ssize_t append (const SxString &);

      // --- Deletion
      void remove (ssize_t charIdx, ssize_t count=1);
      void removeElement (int c);
      void removeFirst ();
      void removeLast ();
      void removeAll ();

      /** Assignment operator for strings */
      SxString &operator=  (const SxString &);
      SxString &operator=  (SxString &&);
      SxString &operator=  (const SxArray<char> &);
      /** Assignment operator for C-like strings */
      SxString &operator=  (const char *);
      /** Returns the concatination with a string */
      inline SxString  operator+  (const SxString &) const;
      inline SxString  operator+  (const char *) const;
      inline SxString  operator+  (const char) const;
      inline SxString  operator+  (const int) const;
      inline SxString  operator+  (const long) const;
      inline SxString  operator+  (const unsigned long) const;
      inline SxString  operator+  (const long long) const;
      inline SxString  operator+  (const unsigned long long) const;
      inline SxString  operator+  (const double) const;
      /** Appends a string to the current object */
      inline void      operator+= (const SxString &);
      inline void      operator+= (const char *);
      inline void      operator+= (const char);
      inline void      operator+= (const int);
      inline void      operator+= (const long);
      inline void      operator+= (const unsigned long);
      inline void      operator+= (const long long);
      inline void      operator+= (const unsigned long long);
      inline void      operator+= (const double);
      /** Returns true if strings are equal (case-sensitive test!). */
      bool      operator== (const SxString &) const;
      bool      operator== (const char *) const;
      /** \brief support of 'if (!str)' expressions
        */
      bool operator! () const;
      /** Returns true if strings differ at least in one character
          (case-sensitive test!). */
      bool      operator!= (const SxString &) const;
      bool      operator!= (const char *) const;
      /** Compares string alphanumerically */
      bool      operator<  (const SxString &) const;
      bool      operator<  (const char *) const;
      /** Compares string alphanumerically */
      bool      operator>  (const SxString &) const;
      bool      operator>  (const char *) const;

      // allow streaming in SX_LOG
      SxString  operator<<  (const SxString &in) const;

      /** Concatenates "wide" characters provided in the zero-terminated
          wchar_t array. Requires that the string object uses Unicode. */
      void replace (const uint16_t *);

      /** Transforms the string object to a zero-terminated sequence of "wide"
          characters. For an empty string, the result is an empty array.
          Unicode is handled automatically. */
      SxArray<uint16_t> toWChars () const;

   protected:
      /** The number of characters which are currently stored in the string.
          This member variable shouldn't be accessed directly. Normally use
          getSize () to read its current value and updateNChars () to change
          it. Such explicit changes are only necessary with Unicode strings,
          otherwise the method resize0 () performs the required bookkeeping. */
      ssize_t nChars;

      /** Whether the string is a UTF-8 encoded general Unicode string.
      (Otherwise it is an 8-bit ASCII string, encoded in Latin-1.) */
      bool isUnicode_;

      /** Whether the member variable 'nChars' is invalid. It can only become
          invalid with Unicode strings. If it is, the function which modified
          the string must call the method updateNChars () before returning,
          especially after resize0 () was called. */
      bool isDirty;

      /** @{These helper functions manage the member variables 'nChars' and
           'isDirty'. */
      inline void setDirty (bool);
      void updateNChars (ssize_t);
      /** @} */

      /** Replaces the first 'len' characters of the C-like string 'in' */
      void replace (const char *in, ssize_t len);

      /** \0-terminated resize */
      void resize0 (ssize_t nBytes, bool keep = true);

      /** Ensures that either both string objects use Unicode or both only use
          ASCII: if one uses Unicode and the other doesn't, the other is
          converted. */
      static void equalize (SxString *s1, SxString *s2);
      void concatenate (const char *, const char *);
      void concatenate (const SxString &, const SxString &);

      /** Returns how often the entire 'needle' occurs in the string object,
          starting at the character index position 'fromCharIdx'. The variable
          'whole' specifies whether this function shall work like
          containsWhole () or contains (), that is, whether it should skip a
          whole match or only skip the first character of a match before
          looking for the next match. */
      int getNOccurrences (const SxConstChar &needle, ssize_t fromCharIdx,
         bool whole, ssize_t needleNChars = -1) const;

      /** Returns a string object with all uppercase characters converted to
          lowercase or vice versa */
      SxString changeCase (bool toLower_) const;

      // --- Members
      static const ssize_t MAX_N_CHARS_FILE;
      static const char *emptyString;
};

/** Shifts the string to a stream */
SX_EXPORT_UTIL std::ostream& operator<< (std::ostream &s, const SxString &in);
#ifdef WIN32
   SX_EXPORT_UTIL std::wostream& operator<< (std::wostream &s, const SxString &in);
#endif

extern SX_EXPORT_UTIL SxList<SxString> __lics;

/** \brief replacement for the standard C printf function */
SX_EXPORT_UTIL int sxprintf (const char *fmt, ...) __SXCHECK_FORMAT(1,2);
/** \brief replacement for the standard C fprintf function
           if fp might be stdout.
 */
SX_EXPORT_UTIL
int sxfprintf (FILE *fp, const char *fmt, ...) __SXCHECK_FORMAT(2,3);

#define SX_CLEARLINE "\n"

#include <SxString.hpp>

#endif /* _SX_STRING_H_ */


#ifndef SX_SEPARATOR
#   define SX_SEPARATOR  SxString("+---------------------------------------"\
                               "--------------------------------------\n")
#endif

