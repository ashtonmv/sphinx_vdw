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

#ifndef _SX_REGEX_H_
#define _SX_REGEX_H_

#include <SxString.h>
#include <SxMap.h>

/** \brief Perl-like regular expressions
    
   \par Regular expression compile options:
   - m Multiline. ^ and $ match newlines within data.
   - i insensitive. Case insensitive match.
   - x ignore whitespace. Ignore white space and # comments in the pattern.
   - s Singleline. Dot (.) matches anything including new lines.
   - U Ungreedy. Invert greediness of quantifiers.
   - u Unicode/UTF-8. Interpret strings as UTF-8 encoded Unicode, not ASCII.
   - A Anchor. Match only at the first position.
   - P POSIX. Return the matches in the order POSIX regexec () would return.
   - D $ does not match the new line at the end of string. It will be ignored
       if option 'm' is set.
   - l literal. Prepand /Q and append /E to pattern. 

    \par Examples:
\code
    SxRegex("(\d)+(\d)+").match (string);
\endcode

    \author Sixten Boeck, boeck@gemmantics.com */
class SX_EXPORT_UTIL SxRegex
{
   public:
      /** Objects of this class contain offsets and lengths of string matches */
      class SX_EXPORT_UTIL Match
      {
         public:
            Match ();
            /** Special constructor; only allowed for 8-bit ASCII (Latin-1)
                strings, not Unicode */
            Match (ssize_t charOffset_, ssize_t nChars_);
            /** General constructor; may be used for arbitrary SxString objects
            */
            Match (ssize_t charOffset_, ssize_t byteOffset_,
               ssize_t nChars_, ssize_t nBytes_);
           ~Match ();

            /** Sets up the Match object; only allowed for 8-bit ASCII
                (Latin-1) strings, not Unicode */
            void set (ssize_t charOffset_, ssize_t nChars_);
            /** Sets up the Match object; may be used for arbitrary
                SxString objects */
            void set (ssize_t charOffset_, ssize_t byteOffset_,
               ssize_t nChars_, ssize_t nBytes_);

            /** Returns the character offset of the match within the input
                string */
            ssize_t getOffset () const;
            /** Returns the number of characters of the match */
            ssize_t getLength () const;
            /** Returns the byte offset of the match within the input string */
            inline ssize_t getByteOffset () const { return byteOffset; }
            /** Returns the number of bytes of the match */
            inline ssize_t getNBytes () const { return nBytes; }

         protected:
            /** The character offset of the match start within the input string
            */
            ssize_t charOffset;
            /** The byte offset of the match start within the input string */
            ssize_t byteOffset;
            /** The length (number of characters) of the match */
            ssize_t nChars;
            /** The length (number of bytes) of the match */
            ssize_t nBytes;
      };

      SxRegex ();
      SxRegex (const SxString &pattern, const SxString &options="");
     ~SxRegex ();

      /** Compiles the given regex pattern, that is, transforms it into an
          internal representation. Returns whether that worked. */
      bool compile (const SxString &pattern,
                    const SxString &options="",
                    bool           exception=true);

      /** Returns the "capture count", that is, the number of "capturing
          subpatterns" */
      ssize_t getCaptureCount () const;
      ssize_t namedCaptureGroupToNumber (const SxString &name) const;

      /** Returns a list of matched substrings. The list is empty if the
          'input' string didn't match. Otherwise, the first element of the list
          usually contains the entire matched string (and the subsequent
          elements contain parenthesized components of the matched string).
          If you instead want the string components as they would have been
          returned by the old SxString::regexp () method (with Perl symbols:
          first "$`", then parenthesized subexpressions, finally "$'"), use the
          compile option "P" (for "POSIX") before calling match (). */
      SxList<SxString> match (const SxString &input) const;

      /** \brief Match all.

          \par Examples:
\code
   SxRegex re("(\\d+)\\.(\\d+)");

   SxList<SxList<SxString> > match;
   match = re.matchAll ("55.117,14.681");

   cout << match << endl;

   // output:
   // 0: 0: 55.117
   //    1: 55
   //    2: 117
   //
   // 1: 0: 14.681
   //    1: 14
   //    2: 681
\endcode
       */
      SxList<SxList<SxString> > matchAll (const SxString &input) const;

      /** \brief Match to capture groups with names.

          \par Examples:
\code
   SxRegex re("(?<day>\\d+)\\s+(?<month>\\w+)");

   SxMap<SxString,SxString> map;
   map = re.matchToMap ("18 March 2015");

   // map:
   // key     = value
   // "0"     = "18 March"
   // "1"     = "18"
   // "2"     = "March"
   // "day"   = "18"
   // "month" = "March"
\endcode
       */
      SxMap<SxString,SxString> matchToMap (const SxString &input) const;      

      /** \brief Return offsets to captured substrings.
          \param text The input text.
          \param startOffset Start offset in the input text.
          \param capturedSubstrings The results.
          \return the number of used elements in capturedSubstrings.

          \par Examples:
\code
   SxRegex re;

   re.compile ("(\\d+)\\s+(?<month>\\w+)");

   SxArray<SxRegex::Match> substrings;
   ssize_t n = re.matchToOffsets ("18 March 2015", 0, &substrings);
   for (ssize_t i=0; i < n; ++i)  {
      cout << i << ": at " << substrings(i).getOffset ()
                << ", len " << substrings(i).getLength () << endl;
   }

   // --- output
   //     0: at 0, len 8
   //     1: at 0, len 2
   //     2: at 3, len 5
\endcode
       */
      ssize_t matchToOffsets (const SxString &input,
                              ssize_t        fromCharIdx,
                              SxArray<Match> *capturedSubstrings) const;

      /** \brief Run the test unit.
       */
      static void test ();

      /** Singleton class for proper memory management of the sole needed
          pcre2_compile_context */
      class Ctx
      {
         public:
            Ctx ();
           ~Ctx ();
            /** Returns the PCRE2 compile context */
            inline void *getContext () const { return ctx; }

         private:
            /** The PCRE2 compile context */
            /*pcre2_compile_context*/ void *ctx;
      };

   protected:

      // --- compile options (PCRE2 flags)
      uint32_t compileOptions;

      // --- meta-options (for which no PCRE2 flags exist)
      enum { None = 0x00, POSIX = 0x01 };
      uint32_t metaOptions;

      /** Helper class for the proper management of the PCRE2-related
          compiled regular expression pattern code and match data memory */
      class PcreMgr
      {
         public:
            PcreMgr ();
            /** Copy constructor */
            PcreMgr (const PcreMgr &in);
           ~PcreMgr ();
            /** Returns the internal pointer to the managed memory. Returns
                NULL if setPtr () wasn't yet called. */
            inline void *getPtr () const { return ptr_; }
            /** Kind of the managed memory */
            enum Kind { None, Code, Data, Ovector };
            /** Sets the internal pointer and kind to the given values and
                starts the reference-counting */
            void setPtr (void *ptr, Kind kind);
            /** Performs an unref () call and resets the internal state of the
                object */
            void reset ();
            /** Assignment operator */
            PcreMgr &operator= (const PcreMgr &in);

         protected:
            /** Internal pointer to the managed memory */
            void *ptr_;
            /** Kind of the managed memory */
            Kind kind_;
            /** Reference counter; only allocated if necessary; basically
                indicates how many PcreMgr objects refer to the 'ptr_' */
            int *refcount_;

            /** Cleans up (deallocates) old PCRE2 memory, if any */
            void cleanPtr ();
            /** Increments the reference counter (if appropriate) */
            inline void ref () { if (refcount_ != NULL)  (*refcount_)++; }
            /** Decrements the reference counter and cleans up old PCRE2 memory
                if appropriate */
            void unref ();
      };

      /** Managers for PCRE2 pattern code and match data memory */
      PcreMgr codeMgr, dataMgr;
      /** Manager for an adapted oVector; necessary if "metaOptions & POSIX".
          Use findOvector () instead of reading from this directly! */
      mutable PcreMgr ovectorMgr;

      // --- capture
      ssize_t captureCount;
      SxMap<SxString, ssize_t> nameToNumber;
      SxMap<ssize_t, SxString> numberToName;

      /** Roughly the number of matches, e.g. as returned by pcre2_match () */
      mutable int rc_;

      /** Re-initializes the SxRegex object */
      void clean ();

      /** Parses the given compile options */
      void parseCompileOptions (const SxString &options);
      /*pcre2_compile_context*/ void *findCompileContext () const;
      /** Compiles the given regex pattern, that is, transforms it into an
          internal representation */
      void compilePattern (const SxString &);

      /** Returns the appropriate oVector of the regular expression match */
      /*PCRE2_SIZE*/ void *findOvector () const;

      /** Adapts the oVector if the meta-option 'POSIX' was given. The
          resulting values are as if from a POSIX regexec () call */
      void adaptToPosix (/*PCRE2_SIZE*/ size_t textSize) const;
      /** Tries to match the given input string against the pattern once,
          starting at the given byte index within the input string. This
          low-level method is mainly an encapsulation of pcre2_match_8 () and
          provides the basis for all other matching-related methods. It stores
          the number of matches (e.g. the "whole" match and the capturing
          groups) in the member variable 'rc_'; on failure, the value is -1. */
      void matchOncePcre (const SxString &input, ssize_t fromByteIdx) const;
      /** Tries to match the given input string against the pattern once,
          starting at the given byte index within the input string. This
          is the common code for the methods match () and matchAll (). */
      SxList<SxString> matchOnceFrom (const SxString &input,
         ssize_t fromByteIdx) const;

      // --- errors
      static SxString patternNotCompiledError ();
      static SxString patternCompileError (const SxString &pattern,
                                           const SxString &options);
      static SxString unknownOptionError (const SxString &option,
                                          const SxString &validOptions);
      static SxString patternCompileError (const SxString &pattern,
                                           int            col,
                                           const SxString &errmsg);
      static SxString matchDataCreateError ();
      static SxString patternInfoError (const SxString &infoName);
      static SxString captureGroupOutOfRangeError (ssize_t number,
                                                   ssize_t captureCount);
      static SxString namedCaptureGroupNotFoundError (const SxString &name,
                                         const SxMap<SxString, ssize_t> &);
      static SxString matchingError (int errcode);

      // Is pattern literal?
      bool isLiteral;
};

#endif /* _SX_REGEX_H_ */
