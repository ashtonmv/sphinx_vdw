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

#ifndef _SX_SED_H_
#define _SX_SED_H_

#include <SxString.h>
#include <SxRegex.h>

/** \brief Transform text with Perl-like regular expressions
   
   \par Usage:
   SxSed commands are compiled from strings with the following syntax: 
   \b s/pattern/substitution/options.
   Second way to compile a command is to specify the pattern,
   the substitution and the options with three strings.
   Text substitutions with SxSed are possible only after a command
   has been compiled.
   - pattern: pcre syntax http://www.pcre.org/current/doc/html/pcre2syntax.html
   - options:
      - g Global. All matches.
      - m Multiline. ^ and $ match newlines within data.
      - i insensitive. Case insensitive match.
      - x ignore whitespace. Ignore white space and # comments in the pattern.
      - s Singleline. Dot (.) matches anything including new lines.
      - U Ungreedy. Invert greediness of quantifiers.
      - u Unicode/UTF-8. Interpret strings as UTF-8 encoded Unicode, not ASCII.
      - A Anchor. Match only at the first position.

    \par Examples:
\code
   // --- a simple replacement
   SxSed("s/day/night/").subst ("day time");  // night time

   // --- global option
   SxSed("s/a/b/g").subst ("aa");             // bb

   // --- global option, the second constructor
   SxSed("a", "b", "g").subst ("aa");         // bb

   // --- multi-line replacement
   SxSed("s/a/b/").subst ("aa\naa\n");        // ba\nba\n

   // --- replace numerical backreferences
   SxSed("s#(\\d+)/(\\d+)/(\\d+)#$3/$1/$2#")  // 1999/03/26
      .subst ("03/26/1999");

   // --- define and replace named backreference
   SxSed("s/(?<foo>(\\d+))/${foo}./").subst ("10"); // 10.

   // --- define, use and replace named backreference
   SxSed("s/(?<bar>(\\d+)).(?P=bar)/$+{bar}./")     // 1.
      .subst ("1-1");
\endcode
    \par Named backreference:
    - (?<name>...) or (?'name'...) or (?P<name>...) defines the group.
    - \\k<name> or \\k'name' or \\g{name} or (?P=name) is a back-reference.
    - To insert the capture in the replacement string,
      the following syntax can be used:
       - the group's number $1
       - .NET ${name}
       - Perl $+{name}
       - Python \\g<name>

    \par Numerical backreferences:
    - \\g{10} is a back-reference to group 10.
    - $10 inserts the 10th capture in the replacement string.

    \par Symbolic backreferences:
    Symbolic backreferences in the replacement string.
    - $'
    - $`
    - $+
    - $&

   \par Comparison with pcre2_substitute
   - pcre2_substitute advantages
      - Supports Unicode.
      - Uses more pcre2 options (PCRE2_ANCHORED, PCRE2_NOTEMPTY_ATSTART).
   - pcre2_substitute disadvantages
      - Seems to return PCRE2_ERROR_NOMEMORY if the allocated size for
        the output is not big enough. That can mean N repeated calls.
        SxSed first determines the size of the output and then works directly
        with the result value without additional copy and allocations.
      - Because it does not compile the replacement string, it needs to
        parse the backreferences for each call (for every line).
      - It looks it does not have symbolic backreferences by default.

 */
class SX_EXPORT_UTIL SxSed
{
   public:
      SxSed ();
      /** Creates an SxSed object from the given command of the form
          "s/pattern/substitute/options" */
      SxSed (const SxString &command);
      /** Creates an SxSed object from the given pattern and substitute strings
      */
      SxSed (const SxString &pattern_,
             const SxString &substitute_,
             const SxString &options_="");
     ~SxSed ();

      /** \brief Prepare the substitution.
          \param command perl-style command 's/pattern/substitute/options'
          \param exception with value false: compile will return false instead
                           of throwing an exception in a case of an error
       */
      bool compile (const SxString &command_, bool exception=true);

      /** \brief Prepare the substitution.
       */
      bool compile (const SxString &pattern_,
                    const SxString &substitute_,
                    const SxString &options_,
                    bool           exception=true);

      /** Substitutes the input text, possibly splitting and re-assembling
          several lines of which the input may consist. Lines must be separated
          by the '\n' character.
       */
      SxString subst (const SxString &input, ssize_t *nHits = NULL) const;

      /** Substitutes the input text, executing several commands by calling the
          other subst () method for each command of the given command list
          sequentially */
      static SxString subst (const SxString      &input,
                             const SxList<SxSed> &commands_);

      /** \brief Run the test unit.
       */
      static void test ();

   protected:
      enum {
         FlagGlobal = 1u,
         FlagTrivial = 2u
      };

      enum RefFlags : uint8_t {
         None       = 0,
         UpperFirst = 1, // for backRef flag \u
         UpperAll   = 2, // for backRef flag \U
         LowerFirst = 3, // for backRef flag \l
         LowerAll   = 4  // for backRef flag \L
      };

      /** The replacement command which shall be executed */
      SxString command;

      // --- pattern
      SxRegex re;

      // --- substitute
      SxString substitute;
      ssize_t nReferences;
      SxArray<ssize_t> referenceId;
      SxArray<ssize_t> referenceCount;
      SxArray<RefFlags> referenceFlag;
      /** Byte index values */
      SxArray<ssize_t> blockOffset;
      /** Numbers of bytes */
      SxArray<ssize_t> blockLength;

      // --- options
      uint32_t flags;

      int32_t substituteFlags;

      /** @{Compiling a command */
      /** Transforms the given command parameters into an internal format */
      void compileCommand (const SxString &pattern_,
                           const SxString &substitute_,
                           const SxString &options_);
      /** Splits the command string and calls the other compileCommand ()
          method with the resulting strings parts */
      void compileCommand ();
      SxString parseOptions (const SxString &);
      void compilePattern (const SxString &, const SxString &);
      void compileSubstitute (const SxString &);
      void compileBackreferences (const SxString &);
      /** @} */

      /** Stores the given reference */
      void storeReference (ssize_t byteIdx, ssize_t refId, RefFlags flag = RefFlags::None);
      /** Parses and skips the name of a named back-reference and stores the
          reference */
      void storeNamedBackRef (ssize_t destByteIdx, SxConstChar::Iterator &it);
      /** Parses and skips the numerical back-reference and stores the
          reference */
      void storeNumBackRef (ssize_t destByteIdx, SxConstChar::Iterator &it,
                            RefFlags flag = RefFlags::None);
      void validateCaptureCount (ssize_t captureGroupNumber);

      /** Performs substitutions within the given individual input line.
          Formerly named execute (). */
      SxString substLine (const SxString &input, ssize_t *nHits) const;

      void clean ();

      static const ssize_t MAX_SUBPATTERNS;

      // --- errors
      SxString commandCompileError () const;
      static SxString commandSyntaxError ();
      static SxString unknownOptionError (const SxString &option,
                                          const SxString &validOptions);
      static SxString backreferenceDelimiterError (SxConstChar::UnicodePoint);
      static SxString backreferenceOutOfRangeError (ssize_t id,
                                                    ssize_t captureCount_);
      SxString substituteError () const;
      static SxString matchingError (int errcode);
      static SxString outputOverflow (ssize_t, ssize_t, ssize_t);
      static SxString substituteOverflow (ssize_t, ssize_t, ssize_t);
      static SxString inputOverflow (ssize_t, ssize_t, ssize_t);
      static SxString inputOverflow (ssize_t, ssize_t, ssize_t, ssize_t);
};

#endif /* _SX_SED_H_ */
