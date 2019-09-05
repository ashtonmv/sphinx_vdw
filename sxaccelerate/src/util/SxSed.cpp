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

// SxSed can be used only if pcre2 is available
#include <SxSed.h>
#include <SxError.h>

const ssize_t SxSed::MAX_SUBPATTERNS = 100;

/*

   https://regex101.com/

   substitution
      \0       complete match contents
      \1       capture group 1
      $1
      ${foo}
      \{foo}
      \g{foo}
      \g<1>    capture group 1

   flags
      g Global
      m Multiline
      i case sensitive
      x ignore whitespace
      s Singleline
x     u Unicode
x     X eXtended
      U Ungreedy
      A Anchor            PCRE2_ANCHORED Match only at the first position

   groups
?     (?'name'...)   named capture group
      (?<name>...)   named capture group
      (?P<name>...)  named capture group
      (?P=<name>)    match subpattern name

   meta sequences
      \k<name> match subpattern name
x     \k'name' match subpattern name
      \k{name} match subpattern name
      \n       match nth subpattern
      \gn      match nth subpattern
      \g{n}    match nth subpattern

 */
SxSed::SxSed ()
   : nReferences(0),
     flags(0)
{
   // empty
}

SxSed::SxSed (const SxString &command_)
   : nReferences(0),
     flags(0)
{
   compile (command_);
}

SxSed::SxSed (const SxString &pattern_,
              const SxString &substitute_,
              const SxString &options_)
   : nReferences(0),
     flags(0)
{
   compile (pattern_, substitute_, options_);
}

SxSed::~SxSed ()
{
   clean ();
}

void SxSed::clean ()
{
   re = SxRegex();

   nReferences = 0;
   flags = 0;

   command = "";
   substitute = "";
   // (Beware: a simple removeAll () for a string wouldn't reset isUnicode!)
}

bool SxSed::compile (const SxString &command_, bool exception)
{
   SX_TRACE ();

   clean ();
   command = command_;

   try {

      compileCommand ();

   } catch (SxException e) {

      clean ();
      if (!exception) return false;
      SX_THROW (e, "RegexCompilationFailed", commandCompileError ());
   }

   return true;
}

bool SxSed::compile (const SxString &pattern_,
                     const SxString &substitute_,
                     const SxString &options_,
                     bool           exception)
{
   SX_TRACE (); 

   clean ();
   command = "s/" + pattern_ + "/" + substitute_ + "/" + options_;

   try {

      compileCommand (pattern_, substitute_, options_);

   } catch (SxException e) {

      clean ();

      if (!exception)
         return false;

      SX_THROW (e, "RegexCompilationFailed", commandCompileError ());
   }

   return true;
}

SxString SxSed::subst (const SxString &input, ssize_t *nHits) const
{
   SX_TRACE ();

   if (nHits)
      *nHits = 0;

   try {

      ssize_t newlineCharIdx = input.find ("\n");
      if (newlineCharIdx < 0 || newlineCharIdx == input.getSize () - 1)  {
         // --- single line
         return substLine (input, nHits);
      } else {
         // --- multi line
         //     TODO: 'm': options PCRE2_MULTILINE collision?
         SxList<SxString> in = input.tokenize ("\n", true);
         SxList<SxString> out;
         SxList<SxString>::ConstIterator it = in.begin ();
         for (; it.isValid (); ++it)  {
            ssize_t nLineHits;
            out.append (substLine (*it, &nLineHits));
            if (nHits)
               *nHits += nLineHits;
         }
         return SxString::join (out, "\n");
      }

   } catch (SxException e) {
      SX_THROW (e, "RegexSubstitutionFailed", substituteError ());
   }

   return SxString();
}

SxString SxSed::subst (const SxString      &input,
                       const SxList<SxSed> &commands_)
{
   SX_TRACE ();

   SxString result = input;
   SxList<SxSed>::ConstIterator it = commands_.begin ();
   for (; it.isValid (); ++it)  {
      ssize_t nHits;
      result = it->subst (result, &nHits);
      if (nHits < 0)
         return result;
   }

   return result;
}

void SxSed::compileCommand ()
{
   SX_DBG_MSG ("command='" << command << "'");

   SxConstChar cmdStr(command.getElems (), command.getNBytes (),
      command.isUnicode ());
   SxConstChar::Iterator it = cmdStr.begin ();
   if (command.getSize () < 4 || *it != 's')
      SX_THROW ("RegexSyntaxError", commandSyntaxError ());
   ++it; // skip 's'

   // --- find delimiter positions for pattern/substitute/options
   SxConstChar::UnicodePoint delimiter = *it;
   ssize_t iArg = 0;
   ssize_t delimCharPos[2];
   bool escaped = false;

   for (++it; it.inRange (); ++it)  {
      SxConstChar::UnicodePoint u = *it;
      if (escaped)  escaped = false;
      else if (u == '\\')  escaped = true;
      else if (u == delimiter)  {
         if (iArg >= 2)  SX_THROW ("RegexSyntaxError", commandSyntaxError ());
         delimCharPos[iArg++] = it.getCharIdx ();
      }
   }
   if ( (escaped) || (iArg != 2) )  SX_THROW ("RegexSyntaxError", commandSyntaxError ());

   // --- create substrings for pattern/substitute/options arguments; compile
   SxString patternPart = command.subString (2, delimCharPos[0] - 1);
   SxString substitutePart = command.subString (delimCharPos[0] + 1,
      delimCharPos[1] - 1);
   SxString optionsPart = command.subString (delimCharPos[1] + 1,
      command.getSize () - 1);
   compileCommand (patternPart, substitutePart, optionsPart);
}

void SxSed::compileCommand (const SxString &pattern_,
                            const SxString &substitute_,
                            const SxString &options_)
{
   SxString patternOptions = parseOptions (options_);
   compilePattern (pattern_, patternOptions);
   compileSubstitute (substitute_);
}

SxString SxSed::parseOptions (const SxString &options_)
{
   SX_DBG_MSG ("options='" << options_ << "'");

   SxString patternOptions;

   // --- options for pcre2_compile()
   //     http://www.pcre.org/current/doc/html/pcre2_compile.html

   for (ssize_t i=0; i < options_.getSize (); ++i)  {
      switch (options_(i))  {
         case 'g': flags |= static_cast<uint32_t>(FlagGlobal); break;
         case 'm':
         case 'i':
         case 'x':
         case 's':
         case 'U':
         case 'u':
         case 'A':
         // own 'literal' flag
         case 'l':           
            patternOptions += options_(i);
            break;
         case 'T': flags |= static_cast<uint32_t>(FlagTrivial); break;
         default:
            SX_THROW ("RegexUnknownOption", unknownOptionError (options_(i), "gmixsUuATl"));
            break;
      }
   }

   return patternOptions;
}

void SxSed::compilePattern (const SxString &pattern_,
                            const SxString &options_)
{
   re.compile (pattern_, options_);
}

void SxSed::compileSubstitute (const SxString &substitute_)
{
   SX_DBG_MSG ("substitute='" << substitute_ << "'");

   nReferences = 0;
   referenceId.resize (SxSed::MAX_SUBPATTERNS + 2);
   referenceCount.resize (SxSed::MAX_SUBPATTERNS + 2);
   referenceFlag.resize (SxSed::MAX_SUBPATTERNS + 2);
   blockOffset.resize (SxSed::MAX_SUBPATTERNS + 2);
   blockLength.resize (SxSed::MAX_SUBPATTERNS + 2);
   referenceId.set (0);
   referenceCount.set (0);
   blockOffset.set (0);
   blockLength.set (0);

   if (flags & static_cast<uint32_t>(FlagTrivial))  {
      substitute = substitute_;
      blockLength(0) = substitute.getNBytes ();
   } else {
      compileBackreferences (substitute_);
   }

   SX_DBG_MSG ("compiled substitute='" << substitute << "':"
               << substitute.getSize ());

   for (ssize_t i = 0; i <= nReferences; ++i)  {
      SX_DBG_MSG ("block " << i
                  << ": offset=" << blockOffset(i)
                  << ", length=" << blockLength(i));
   }
}

void SxSed::storeReference (ssize_t byteIdx, ssize_t refId, RefFlags flag)
{
   if ( /* (refId < 0) || */ (refId >= referenceId.size) )  {
      SX_THROW (backreferenceOutOfRangeError (refId, MAX_SUBPATTERNS));
      // (Here MAX_SUBPATTERNS is the limit from an application's point of
      // view. The actual internal limit for this class is slightly higher.)
   }
   referenceId(nReferences) = refId;
   referenceCount(refId)++;
   referenceFlag(nReferences) = flag;
   blockLength(nReferences) = byteIdx - blockOffset(nReferences);
   nReferences++;
   blockOffset(nReferences) = byteIdx;
}

void SxSed::storeNamedBackRef (ssize_t destByteIdx, SxConstChar::Iterator &it)
{
   ssize_t res;
   SxConstChar::UnicodePoint delimiter, u = *it;
   switch (u)  {
      case '{': delimiter = '}'; break;
      case '<': delimiter = '>'; break;
      default: SX_THROW ("RegexBackRefDelimiterError", backreferenceDelimiterError (u)); break;
   }
   ++it; // skip the opening delimiter
   ssize_t fromCharIdx = it.getCharIdx ();
   while (*it != delimiter)  ++it; // progress to the matching closing delimiter
   ssize_t toCharIdx = it.getCharIdx () - 1; // name ends before delimiter
   ssize_t nChars = toCharIdx - fromCharIdx + 1;
   SX_CHECK (nChars > 0, nChars); // empty back-reference names are invalid
   SxString name = SxString::subString (*(it.getCharObj ()), fromCharIdx,
      toCharIdx);
   SX_DBG_MSG ("back-reference name '" << name << "'");
   res = re.namedCaptureGroupToNumber (name);
   storeReference (destByteIdx, res);
   ++it; // skip the closing delimiter
}

void SxSed::storeNumBackRef (ssize_t destByteIdx, SxConstChar::Iterator &it,
                             RefFlags flag)
{
   ssize_t value = 0;
   do  {
      value = 10 * value + ((*it) - '0');
      ++it;
   } while ( (it.inRange ()) && (it.isDigit ()) );
   validateCaptureCount (value);
   storeReference (destByteIdx, value, flag);
}

void SxSed::compileBackreferences (const SxString &substitute_)
{
   substitute = substitute_;
      // (starting with an exact copy to ensure appropriate size and isUnicode_)
   char *origDest = substitute.elements, *dest = origDest;
   ssize_t origNBytes = substitute_.getNBytes ();
   bool isUni = substitute_.isUnicode ();
   SxConstChar srcStr(substitute_.getElems (), origNBytes, isUni);
   SxConstChar::Iterator it = srcStr.begin ();

      RefFlags refFlag = RefFlags::None;

   while (it.inRange ())  {
      SxConstChar::UnicodePoint u = *it;
      // For Debug
      // printf("%c\n",u);
      if (u == '\\')  { // start of an escape sequence, maybe a back-reference
         ++it;
         if (!it.inRange ())  { // '\' was at end of string
            *dest++ = '\\';
            break;
         }
         switch (*it)  {
            case '{': // \{name}
               storeNamedBackRef (dest - origDest, it);
               break;
            case 'g': // expecting \g{name} or \g<name>
               ++it;
               storeNamedBackRef (dest - origDest, it);
               break;
            // simple escape character translations:
            case 'a' : *dest++ = 7; ++it; break;
            case 'e' : *dest++ = 27; ++it; break;
            case 'f' : *dest++ = '\f'; ++it; break;
            case 'n' : *dest++ = '\n'; ++it; break;
            case 'r' : *dest++ = '\r'; ++it; break;
            case 't' : *dest++ = '\t'; ++it; break;
            case '\\': *dest++ = '\\'; ++it; break;
            case '0' : *dest++ = '\0'; ++it; break;
            case 'u' : refFlag = RefFlags::UpperFirst; ++it; break;
            case 'U' : refFlag = RefFlags::UpperAll;   ++it; break;
            case 'l' : refFlag = RefFlags::LowerFirst; ++it; break;
            case 'L' : refFlag = RefFlags::LowerAll;   ++it; break;
            default:
               if (!it.isDigit ())  SX_EXIT
               storeNumBackRef (dest - origDest, it, refFlag); // \1
               refFlag = RefFlags::None;
               break;
         }
      }  else if (u == '$')  { // maybe some other style of back-reference
         ++it;
         if (!it.inRange ())  { // '$' was at end of string
            *dest++ = '$';
            break;
         }
         switch (*it)  {
            case '{': // ${name}
               storeNamedBackRef (dest - origDest, it);
               break;
            case '\'': // $'
               storeReference (dest - origDest, MAX_SUBPATTERNS + 1);
               ++it;
               break;
            case '`': // $`
               storeReference (dest - origDest, MAX_SUBPATTERNS);
               ++it;
               break;
            case '+': // $+ or $+{name}
               ++it;
               if ( (it.inRange ()) && (*it == '{') )  { // $+{name}
                  storeNamedBackRef (dest - origDest, it);
               }
               else  storeReference (dest - origDest, re.getCaptureCount ());
               break;
            case '&': // $&
               ++it;
               storeReference (dest - origDest, 0);
               break;
            default:
               if (!it.isDigit ())  *dest++ = '$';
               else {
                  storeNumBackRef (dest - origDest, it, refFlag); // $1
                  refFlag = RefFlags::None;
               }
               break;
         }
      }  else  { // keep literal character
         it.encode (&dest);
         ++it;
      }
   }
   ssize_t actualNBytes = dest - origDest;
   if (actualNBytes < origNBytes)  {
      ssize_t actualNChars = SxConstChar::countChars (origDest, dest, isUni);
      substitute.resize (actualNChars, true); // shrink
   }
   blockLength(nReferences) = actualNBytes - blockOffset(nReferences);
}

void SxSed::validateCaptureCount (ssize_t captureGroupNumber_)
{
   if (captureGroupNumber_ > re.getCaptureCount ())
      SX_THROW ("RegexCaptureOutOfRange", backreferenceOutOfRangeError (
                captureGroupNumber_, re.getCaptureCount ()));
}

SxString SxSed::substLine (const SxString &input, ssize_t *nHits) const
{
   SX_TRACE ();

   SxString substi = substitute;
   if (input.isUnicode () && !substi.isUnicode ())  {
      substi = SxString::asciiToUnicode (substi.getElems (),
         substi.getNBytes ());
   }

   SxString result;
   if (input.isUnicode ()) // result must use Unicode too
      result = SxString::asciiToUnicode (NULL);

   if (nHits)
      *nHits = 0;

   bool isGlobal = (flags & static_cast<uint32_t>(FlagGlobal));
   ssize_t startOffset = 0; // character offset
   ssize_t n = 0;
   ssize_t inputNBytes = input.getNBytes ();
   ssize_t outputNBytes = inputNBytes;
   ssize_t rc = 0;

   SxList<SxArray<SxRegex::Match> > matches;
   SxList<SxString> res;
   SxArray<SxRegex::Match> submatch;
   //SxList<SxRegex::Match> submatch;

   do {

      rc = re.matchToOffsets (input, startOffset, &submatch);

      if (rc > 0)  {
         submatch.resize (SxSed::MAX_SUBPATTERNS + 2, true);

         for (ssize_t i=0; i < rc; ++i)  {
            n = submatch(i).getNBytes () * referenceCount(i); // submatch reference
            outputNBytes += n;
			SX_DBG_MSG("submatch " << i
				<< ": char offset=" << submatch(i).getOffset());
            SX_DBG_MSG(", nChars=" << submatch(i).getLength ()
                        << ", outputNBytes+=" << n
                        << ", outputNBytes=" << outputNBytes);
         }
         // --- example
         //    compile ("s/^(\\w+)\\s+(\\w+)/A/")
         //    subst ("Mon Mar")
         //
         //    submatch 0: 'Mon Mar', offset=0, length=7
         //    submatch 1: 'Mon', offset=0, length=3
         //    submatch 2: 'Mar', offset=4, length=3

         n = substi.getNBytes () - submatch(0).getNBytes ();
         outputNBytes += n;
         SX_DBG_MSG ("replace outputNBytes+=" << n
                     << ", outputNBytes=" << outputNBytes);

         ssize_t entireFromCharIdx = submatch(0).getOffset ();
         ssize_t entireToCharIdx = entireFromCharIdx + submatch(0).getLength ();
         ssize_t entireFromByte = submatch(0).getByteOffset ();
         ssize_t entireToByte = entireFromByte + submatch(0).getNBytes ();

         // --- before
         submatch(SxSed::MAX_SUBPATTERNS).set (0, entireFromCharIdx);
         n = entireFromByte * referenceCount(SxSed::MAX_SUBPATTERNS);
         outputNBytes += n;
         SX_DBG_MSG ("before outputNBytes+=" << n
                     << ", outputNBytes=" << outputNBytes);

         // --- after
         submatch(SxSed::MAX_SUBPATTERNS + 1).set (entireToCharIdx,
                                                   input.getSize ()
                                                   - entireToCharIdx - 1);
         n = (inputNBytes - entireToByte)
           * referenceCount(SxSed::MAX_SUBPATTERNS + 1);
         outputNBytes += n;
         SX_DBG_MSG ("after outputNBytes+=" << n
                     << ", outputNBytes=" << outputNBytes);

         matches.append (submatch);

         if (!isGlobal)
            break;

         if (entireToCharIdx == startOffset)  {
            // --- empty
            if (startOffset < input.getSize())
               startOffset++;
            else
               break;
         } else {
            // --- continue from the first character
            //     after the end of the entire pattern
            startOffset = entireToCharIdx;
         }
      }

   } while (rc > 0);

   SX_DBG_MSG ("nMatches=" << matches.getSize ());
   SX_DBG_MSG ("outputNBytes=" << outputNBytes);

   if (nHits)
      *nHits = matches.getSize ();

   if (outputNBytes < 1)
      return result;

   result.resize (outputNBytes);

   // --- copy and replace
   const char *srcSubst = substi.getElems ();
   const char *src = input.getElems ();
   char *dst = result.elements;
   ssize_t posSrc = 0; // byte index
   ssize_t posDst = 0; // byte index
   ssize_t substNBytes = substi.getNBytes ();

   SxList<SxArray<SxRegex::Match> >::ConstIterator it;
   it = matches.begin ();
   while (it.isValid ())  {
      const SxArray<SxRegex::Match> &match = *it;
      SxRegex::Match m0 = match(0);
      ssize_t m0ByteOffset = m0.getByteOffset ();

      n = m0ByteOffset - posSrc;
      if (posDst + n > outputNBytes + 1)
         SX_THROW (outputOverflow (posDst, n, outputNBytes + 1));
      if (posSrc + n > inputNBytes + 1)
         SX_THROW (inputOverflow (posSrc, n, inputNBytes + 1));
      memcpy (dst, src + posSrc, static_cast<size_t>(n));
      dst += n;
      posDst += n;

      for (ssize_t i = 0; i <= nReferences; ++i)  {
         n = blockLength(i);
         if (posDst + n > outputNBytes + 1 )
            SX_THROW (outputOverflow (posDst, n, outputNBytes + 1));
         if (blockOffset(i) + n > substNBytes + 1 )
            SX_THROW (substituteOverflow (blockOffset(i), n, substNBytes + 1));
         memcpy (dst, srcSubst + blockOffset(i), static_cast<size_t>(n));



         dst += n;
         posDst += n;

         if (   i != nReferences
             && referenceId(i) < SxSed::MAX_SUBPATTERNS + 2
             && referenceId(i) < match.getSize ()
             && (n = match(referenceId(i)).getNBytes ()) > 0)
         {
            ssize_t pos = match(referenceId(i)).getByteOffset ();
            if (posDst + n > outputNBytes + 1)
               SX_THROW (outputOverflow (posDst, n, outputNBytes + 1));
            if (pos + n > inputNBytes + 1)
               SX_THROW (inputOverflow (referenceId(i), pos, n, inputNBytes + 1));

            if (referenceFlag(i) == UpperFirst) {
               memcpy (dst, src + pos, static_cast<size_t>(n));
               dst[0] = (char)toupper(dst[0]);
            } else if (referenceFlag(i) == UpperAll) {
               for (size_t k = 0; k < static_cast<size_t>(n); ++k) {
                  dst[k] = (char)toupper (*(src+pos+k));
               }
            } else if (referenceFlag(i) == LowerFirst) {
               memcpy (dst, src + pos, static_cast<size_t>(n));
               dst[0] = (char)tolower (dst[0]);
            } else if (referenceFlag(i) == LowerAll) {
               for (size_t k = 0; k < static_cast<size_t>(n); ++k) {
                  dst[k] = (char)tolower (*(src+pos+k));
               }
            } else {
               memcpy (dst, src + pos, static_cast<size_t>(n));
            }
            dst += n;
            posDst += n;
         }
      }

      posSrc = m0ByteOffset + m0.getNBytes ();

      ++it;
   }

   n = inputNBytes - posSrc;
   if (posDst + n > outputNBytes + 1)
      SX_THROW (outputOverflow (posDst, n, outputNBytes + 1));
   if (llabs(posSrc + n) > inputNBytes + 1)
      SX_THROW (inputOverflow (posSrc, n, inputNBytes + 1));
   memcpy (dst, src + posSrc, static_cast<size_t>(n));

   return result;
}

static bool sxSedTest (const SxString &command_,
                       const SxString &input,
                       const SxString &result_,
                       ssize_t        nHits_=-1)
{
   SxString result;
   try {
      ssize_t nHits = 0;
      result = SxSed(command_).subst (input, &nHits);
      //cout << "expected '" << result_ << "':" << result_.getSize () << endl;
      //cout << "got      '" << result << "':" << result.getSize () << endl;
      if (nHits_ >= 0 && nHits_ != nHits)
         SX_THROW ("The number of hits " + SxString(nHits)
                   + " is different from the expected value "
                   + SxString(nHits_));

   } catch (SxException e) {
      e.print ();
      return false;
   }

   return (result == result_);
}

void SxSed::test ()
{
   // --- empty
   SX_CHECK (sxSedTest ("s/a/b/",  "", ""));

   // --- replace
   SX_CHECK (sxSedTest ("s/a/b/",  "aa", "ba", 1));

   // --- options
   SX_CHECK (sxSedTest ("s/a/b/g", "aa", "bb", 2));
   SX_CHECK (sxSedTest ("s/a/b/i", "AB", "bB"));

   // --- no match
   SX_CHECK (sxSedTest ("s/a/b/",  "c", "c", 0));

   // --- remove
   SX_CHECK (sxSedTest ("s/b//", "abc", "ac"));

   // --- add
   SX_CHECK (sxSedTest ("s/b/bc/", "abd", "abcd"));

   // --- newlines
   SX_CHECK (sxSedTest ("s/a/b/",  "a\n", "b\n"));
   SX_CHECK (sxSedTest ("s/a/b/",  "a\na", "b\nb"));
   SX_CHECK (sxSedTest ("s/a/b/",  "a\na\n", "b\nb\n"));
   SX_CHECK (sxSedTest ("s/a/b/",  "\n", "\n"));
   SX_CHECK (sxSedTest ("s/a/b/",  "\n\n", "\n\n"));
   SX_CHECK (sxSedTest ("s/a/b/",  "aa\naa\n", "ba\nba\n"));
   SX_CHECK (sxSedTest ("s/a/b/g", "aa\naa\n", "bb\nbb\n", 4));

   // --- compiled backreferences
   SX_CHECK (sxSedTest ("s/(hello)(\\s*=\\s*)/_$1_ $2/",
                        "cout hello=world endl",
                        "cout _hello_ =world endl"));

   // --- backreferences
   SX_CHECK (sxSedTest ("s#(\\d+)/(\\d+)/(\\d+)#$3/$1/$2#",
                        "03/26/1999", "1999/03/26"));
   SX_CHECK (sxSedTest ("s#(\\d+)/(\\d+)/(\\d+)#\\3/\\1/\\2#",
                        "03/26/1999", "1999/03/26"));

   // --- named backreferences
   SX_CHECK (sxSedTest ("s/(?<foo>(\\d+))/$1./", "10", "10."));
   SX_CHECK (sxSedTest ("s/(?<foo>(\\d+))/${foo}./", "10", "10."));
   SX_CHECK (sxSedTest ("s/(?<foo>(\\d+)).(?P=foo)/$+{foo}./", "1-1", "1."));
}

SxString SxSed::commandCompileError () const
{
   return "Cant't compile sed command '" + command + "'";
}

SxString SxSed::commandSyntaxError ()
{
   return "Command is not in perl-style format. "
          "Enter a command in the same format as this example: "
          "s/pattern/substitution/options";
}

SxString SxSed::unknownOptionError (const SxString &option_,
                                    const SxString &validOptions_)
{
   return "Unknown sed option '" + option_ + "'. "
          "Available options are [" + validOptions_ + "].";
}

SxString SxSed::backreferenceDelimiterError (SxConstChar::UnicodePoint u)
{
   return SxString("Unknown delimiter '") + ((int) u)
          + "' to read a named backreference. " +
          "Available syntax is ${name} $+{name} \\k<name> \\g<name>";
}

SxString SxSed::backreferenceOutOfRangeError (ssize_t id_,
                                              ssize_t max_)
{
   return "Backreference " + SxString(id_) + " is out of range. "
          "The maximal number of references is "
          + SxString(max_) + ".";
}

SxString SxSed::substituteError () const
{
   return "Cant't substitute the input text "
          "with sed command '" + command + "'";
}

SxString SxSed::matchingError (int errcode)
{
   return "Matching error " + SxString(errcode) + ".";
}

SxString SxSed::outputOverflow (ssize_t pos, ssize_t n, ssize_t size)
{
   return "Invalid calculation of the output size. "
          "Attempt to write " + SxString(n) + " bytes "
          " to offset " + SxString(pos) + " in the buffer "
          "with size "+ SxString(size) + " bytes.";
}

SxString SxSed::substituteOverflow (ssize_t pos, ssize_t n, ssize_t size)
{
   return "Invalid compilation of the substitute string. "
          "Attempt to read " + SxString(n) + " bytes "
          " from offset " + SxString(pos) + " in the buffer "
          "with size "+ SxString(size) + " bytes.";
}

SxString SxSed::inputOverflow (ssize_t pos, ssize_t n, ssize_t size)
{
   return "Invalid compilation of the substitute string. "
          "Attempt to read " + SxString(n) + " bytes "
          " from offset " + SxString(pos) + " in the buffer "
          "with size "+ SxString(size) + " bytes.";
}

SxString SxSed::inputOverflow (ssize_t i, ssize_t pos, ssize_t n, ssize_t size)
{
   return "Invalid compilation of the substitute string. "
          "Attempt to read block " + SxString(i) + " with "
          + SxString(n) + " bytes "
          " from offset " + SxString(pos) + " in the buffer "
          "with size "+ SxString(size) + " bytes.";
}


