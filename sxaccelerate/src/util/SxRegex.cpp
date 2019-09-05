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

#include <SxConfig.h>
#include <SxRegex.h>
#include <SxError.h>

// --- include 8-bit PCRE2
#ifdef USE_PCRE2
#define PCRE2_CODE_UNIT_WIDTH  8
#include <pcre2.h>
#else
namespace {
   ATTR_NO_RETURN void needPcre2 ()
   {
      cout << "This feature needs pcre2. Configure with --enable-pcre2" << endl;
      SX_EXIT;
   }
}
#endif

SxRegex::Ctx::Ctx ()
{
#ifdef USE_PCRE2
   ctx = pcre2_compile_context_create (NULL);
   if (ctx == NULL)  sxOutOfMemoryHandler ();
   pcre2_set_bsr ((pcre2_compile_context *) ctx, PCRE2_BSR_UNICODE);
#else
	// silently ignore needPcre2 ();
#endif
}

SxRegex::Ctx::~Ctx ()
{
#ifdef USE_PCRE2
   pcre2_compile_context_free ((pcre2_compile_context *) ctx);
#else
	// silently ignore needPcre2 ();
#endif
}

// ---------------------------------------------------------------------------

SxRegex::PcreMgr::PcreMgr () : ptr_(NULL), kind_(None), refcount_(NULL)
{
   // empty
}

SxRegex::PcreMgr::PcreMgr (const PcreMgr &in)
  : ptr_(in.ptr_), kind_(in.kind_), refcount_(in.refcount_)
{
   ref ();
}

SxRegex::PcreMgr::~PcreMgr ()
{
   unref ();
}

void SxRegex::PcreMgr::setPtr (void *ptr, Kind kind)
{
   SX_CHECK (ptr != NULL);
   SX_CHECK (kind == Code || kind == Data || kind == Ovector, kind);
   SX_CHECK (ptr_ == NULL);
   SX_CHECK (kind_ == None, kind_);
   SX_CHECK (refcount_ == NULL);
   refcount_ = new int; // from now on, we need reference-counting
   *refcount_ = 0;
   ptr_ = ptr;
   kind_ = kind;
   ref ();
}

void SxRegex::PcreMgr::reset ()
{
   unref ();
   ptr_ = NULL;
   kind_ = None;
   refcount_ = NULL;
}

SxRegex::PcreMgr &SxRegex::PcreMgr::operator= (const SxRegex::PcreMgr &in)
{
   unref ();
   ptr_ = in.ptr_;
   kind_ = in.kind_;
   refcount_ = in.refcount_;
   ref ();
   return *this;
}

void SxRegex::PcreMgr::unref ()
{
   if (refcount_ == NULL)  return; // no reference-counting; nothing to do
   (*refcount_)--;
   if (*refcount_ > 0)  return; // don't yet deallocate
   cleanPtr ();
   delete refcount_;
   refcount_ = NULL;
}

void SxRegex::PcreMgr::cleanPtr ()
{
#ifdef USE_PCRE2
   if (ptr_ == NULL)  return; // nothing to do
   switch (kind_)  {
      case Code:
         SX_DBG_MSG ("call pcre2_code_free()");
         pcre2_code_free ((pcre2_code *) ptr_);
         break;
      case Data:
         SX_DBG_MSG ("call pcre2_match_data_free()");
         pcre2_match_data_free ((pcre2_match_data *) ptr_);
         break;
      case Ovector:
         delete [] (PCRE2_SIZE *) ptr_;
         break;
      default: SX_THROW ("Invalid PcreMgr kind " + SxString(kind_)); break;
   }
   ptr_ = NULL;
#else
	needPcre2 ();
#endif
}

// ---------------------------------------------------------------------------

SxRegex::Match::Match ()
{
   set (0, 0);
}

SxRegex::Match::Match (ssize_t charOffset_, ssize_t nChars_)
{
   set (charOffset_, nChars_);
}

SxRegex::Match::~Match ()
{
   // empty
}

void SxRegex::Match::set (ssize_t charOffset_, ssize_t nChars_)
{
   charOffset = byteOffset = charOffset_;
   nChars = nBytes = nChars_;
}

void SxRegex::Match::set (ssize_t charOffset_, ssize_t byteOffset_,
   ssize_t nChars_, ssize_t nBytes_)
{
   charOffset = charOffset_;
   byteOffset = byteOffset_;
   nChars = nChars_;
   nBytes = nBytes_;
}

ssize_t SxRegex::Match::getOffset () const
{
   return charOffset;
}

ssize_t SxRegex::Match::getLength () const
{
   return nChars;
}

// ---------------------------------------------------------------------------

SxRegex::SxRegex ()
   : compileOptions(0),
     metaOptions(None),
     captureCount(0),
     isLiteral(false)
{
   // empty
}

SxRegex::SxRegex (const SxString &pattern_,
                  const SxString &options_)
   : compileOptions(0),
     metaOptions(None),
     captureCount(0),
     isLiteral(false)
{
   compile (pattern_, options_);
}

SxRegex::~SxRegex ()
{
   clean ();
}

void SxRegex::clean ()
{
   codeMgr.reset ();
   dataMgr.reset ();
   ovectorMgr.reset ();

   compileOptions = 0;
   metaOptions = None;
   captureCount = 0;  

   nameToNumber.removeAll ();
   numberToName.removeAll ();
}

ssize_t SxRegex::getCaptureCount () const
{
   return captureCount;
}

ssize_t SxRegex::namedCaptureGroupToNumber (const SxString &name_) const
{
   if (!nameToNumber.containsKey (name_))
      SX_THROW ("RegexNamedCaptureNotFound", namedCaptureGroupNotFoundError (name_, nameToNumber));

   return nameToNumber(name_);
}

bool SxRegex::compile (const SxString &pattern_,
                       const SxString &options_,
                       bool           exception)
{

   clean ();

   try {

      parseCompileOptions (options_);
      compilePattern (pattern_);

   } catch (SxException e) {

      clean ();
      if (!exception) return false;
      SX_THROW (e, "RegexCompilationFailed",
                patternCompileError (pattern_, options_));
   }

   return true;
}

void SxRegex::parseCompileOptions (const SxString &options_)
{
#ifdef USE_PCRE2
   SX_DBG_MSG ("options='" << options_ << "'");

   // --- pcre2_compile (.., uint32_t options)
   //     pcre2_match (.., uint32_t options)
   //
   //     http://www.pcre.org/current/doc/html/pcre2_compile.html
   //     http://www.pcre.org/current/doc/html/pcre2_match.html

   for (ssize_t i=0; i < options_.getSize (); ++i)  {
      switch (options_(i))  {
         case 'm': compileOptions |= PCRE2_MULTILINE; break;
         case 'i': compileOptions |= PCRE2_CASELESS; break;
         case 'x': compileOptions |= PCRE2_EXTENDED; break;
         case 's': compileOptions |= PCRE2_DOTALL; break;
         case 'U': compileOptions |= PCRE2_UNGREEDY; break;
         case 'u': compileOptions |= PCRE2_UTF | PCRE2_NO_UTF_CHECK; break;
         case 'A': compileOptions |= PCRE2_ANCHORED; break;
         case 'D': compileOptions |= PCRE2_DOLLAR_ENDONLY; break;
         case 'P': metaOptions    |= POSIX; break;
         case 'l': isLiteral = true; break;          
         default:
            SX_THROW ("RegexUnknownOption", unknownOptionError (options_(i), "mixsUuAPl"));
            break;

         // --- unused pcre2_compile options
         // PCRE2_ALT_BSUX           Alternative handling of \u, \U, and \x
         // PCRE2_AUTO_CALLOUT       Compile automatic callouts
         // PCRE2_DUPNAMES           Allow duplicate names for subpatterns
         // PCRE2_MATCH_UNSET_BACKREF  Match unset back references
         // PCRE2_NEVER_UCP          Lock out PCRE2_UCP, e.g. via (*UCP)
         // PCRE2_NEVER_UTF          Lock out PCRE2_UTF, e.g. via (*UTF)
         // PCRE2_NO_AUTO_CAPTURE    Disable numbered capturing paren-
         //                           theses (named ones available)
         // PCRE2_NO_AUTO_POSSESS    Disable auto-possessification
         // PCRE2_NO_DOTSTAR_ANCHOR  Disable automatic anchoring for .*
         // PCRE2_NO_START_OPTIMIZE  Disable match-time start optimizations
         // PCRE2_UCP                Use Unicode properties for \d, \w, etc.
      }
   }
#else
   needPcre2 (); 
#endif
}

static const SxRegex::Ctx regexCtx; // initialized when the program starts

void *SxRegex::findCompileContext () const
{
#ifdef USE_PCRE2
   if (!(compileOptions & PCRE2_UTF))  return NULL; // nothing special necessary
   return (void *) regexCtx.getContext ();
#else
   needPcre2 ();
   return NULL;
#endif
}

void SxRegex::compilePattern (const SxString &patternIn)
{
#ifdef USE_PCRE2
   SxString pattern_;
   if (isLiteral) pattern_ = "\\Q" + patternIn + "\\E";
   else pattern_ = patternIn;

   SX_DBG_MSG ("pattern='" << pattern_ << "'");
   int errNo = 0;
   PCRE2_SIZE errOffset = 0;

   // --- compile pattern
   PCRE2_SPTR pattern = (PCRE2_SPTR) pattern_.getElems ();
   pcre2_compile_context *ctx = (pcre2_compile_context *) findCompileContext ();
   pcre2_code *code = pcre2_compile (pattern, PCRE2_ZERO_TERMINATED,
                                   compileOptions, &errNo, &errOffset, ctx);
   if (code == NULL)  {
      PCRE2_UCHAR buffer[256];
      pcre2_get_error_message (errNo, buffer, sizeof(buffer));
      SX_THROW ("RegexCompilationFailed",
                patternCompileError (pattern_, (int) errOffset,
                                     SxString((const char *)buffer)));
   }
   codeMgr.setPtr ((void *) code, PcreMgr::Code);

   // --- match data
   //     http://www.pcre.org/current/doc/html/pcre2api.html
   //     Information about a successful or unsuccessful
   //     match is placed in a match data block
   //     A match data block can be used many times,
   //     with the same or different compiled patterns.
   //
   //     The ovector is created to be exactly the right size
   //     to hold all the substrings a pattern might capture.
   //
   pcre2_match_data *data = pcre2_match_data_create_from_pattern (code, NULL);
   if (data == NULL)  SX_THROW ("RegexCompilationFailed", matchDataCreateError ());
   dataMgr.setPtr ((void *) data, PcreMgr::Data);

   // --- captureCount
   uint32_t count = 0;
   errNo = pcre2_pattern_info (code, PCRE2_INFO_CAPTURECOUNT, &count);
   if (errNo != 0)
      SX_THROW ("RegexCompilationFailed",
                patternInfoError ("PCRE2_INFO_CAPTURECOUNT"));
   captureCount = static_cast<ssize_t>(count);
   SX_DBG_MSG ("captureCount=" << captureCount);

   // --- named substrings
   uint32_t nameCount = 0;
   errNo = pcre2_pattern_info (code, PCRE2_INFO_NAMECOUNT, &nameCount);
   if (errNo != 0)
      SX_THROW ("RegexCompilationFailed",
                patternInfoError ("PCRE2_INFO_NAMECOUNT"));
   SX_DBG_MSG ("nameCount=" << nameCount);

   if (nameCount > 0)  {
      PCRE2_SPTR nameTable;
      uint32_t entrySize = 0;

      errNo = pcre2_pattern_info (code, PCRE2_INFO_NAMETABLE, &nameTable);
      if (errNo != 0)
         SX_THROW ("RegexCompilationFailed",
                   patternInfoError ("PCRE2_INFO_NAMETABLE"));

      errNo = pcre2_pattern_info (code, PCRE2_INFO_NAMEENTRYSIZE, &entrySize);
      if (errNo != 0)
         SX_THROW ("RegexCompilationFailed",
                   patternInfoError ("PCRE2_INFO_NAMEENTRYSIZE"));

      PCRE2_SPTR tabptr = nameTable;
      for (uint32_t i = 0; i < nameCount; i++)  {
         ssize_t n = (tabptr[0] << 8) | tabptr[1];
         if (n > captureCount)
            SX_THROW ("RegexCaptureOutOfRange", captureGroupOutOfRangeError (n, captureCount));
         SxString name = reinterpret_cast<const char*>(tabptr + 2);
         nameToNumber(name) = n;
         numberToName(n) = name;
         SX_DBG_MSG ("capture group name " << i << ": number " << n
                     << ", '" << name << "'");
         tabptr += entrySize;
      }
   }
#else
   needPcre2 ();
#endif
}

/*PCRE2_SIZE*/ void *SxRegex::findOvector () const
{
#ifdef USE_PCRE2
   if (metaOptions & POSIX)  return ovectorMgr.getPtr ();
   pcre2_match_data *matchData = (pcre2_match_data *) dataMgr.getPtr ();
   return pcre2_get_ovector_pointer (matchData);
#else
   needPcre2 ();
   return NULL;
#endif
}

void SxRegex::adaptToPosix (size_t textSize) const
{
#ifdef USE_PCRE2
   int oldRc = rc_, rc2 = oldRc + 1;
   PCRE2_SIZE *customOvector = new PCRE2_SIZE [static_cast<size_t>(2 * rc2)];
   if (customOvector == NULL)  sxOutOfMemoryHandler ();
   ovectorMgr.setPtr ((void *) customOvector, PcreMgr::Ovector);
   pcre2_match_data *matchData = (pcre2_match_data *) dataMgr.getPtr ();
   PCRE2_SIZE *origOvector = pcre2_get_ovector_pointer (matchData);
   PCRE2_SIZE *dest = customOvector;
   PCRE2_SIZE fromByte = *origOvector++, limitByte = *origOvector++;
   *dest++ = 0;
   *dest++ = fromByte; // for the text portion _before_ the entire match
   while (--oldRc > 0)  { // copy the positions of captured sub-expressions
      *dest++ = *origOvector++;
      *dest++ = *origOvector++;
   }
   *dest++ = limitByte;
   *dest = textSize; // for the text portion _after_ the entire match
   rc_ = rc2;
#else
   needPcre2 ();
#endif
}

void SxRegex::matchOncePcre (const SxString &input, ssize_t fromByteIdx)
   const
{
#ifdef USE_PCRE2
   SX_TRACE ();
   pcre2_code *code = (pcre2_code *) codeMgr.getPtr ();
   if (!code)  SX_THROW (patternNotCompiledError ());
   PCRE2_SPTR8 text = (PCRE2_SPTR8) input.getElems ();
   PCRE2_SIZE inputNBytes = (PCRE2_SIZE) input.getNBytes ();
   uint32_t matchOptions = 0;
   pcre2_match_data *matchData = (pcre2_match_data *) dataMgr.getPtr ();
   rc_ = pcre2_match_8 (code, text, inputNBytes, static_cast<size_t>(fromByteIdx), matchOptions,
                        matchData, NULL);
   if ( (rc_ <= 0) && (rc_ != PCRE2_ERROR_NOMATCH) )  {
      // --- rc
      //      0 => success, but ovector is not big enough
      //     -1 => failed to match (PCRE2_ERROR_NOMATCH); no problem
      //     -2 => partial match (PCRE2_ERROR_PARTIAL)
      //   < -2 => some kind of unexpected problem
      if ( (rc_ != PCRE2_ERROR_BADOFFSET) || (inputNBytes > 0) )  {
         SX_THROW (matchingError (rc_));
      }
   }
   if (rc_ < 0)  rc_ = -1; // "hide" the exact PCRE2-internal error value
   else if ( (rc_ > 0) && (metaOptions & POSIX) )  adaptToPosix (inputNBytes);
#else
   needPcre2 ();
#endif
}

SxList<SxString> SxRegex::matchOnceFrom (const SxString &input,
   ssize_t fromByteIdx) const
{
#ifdef USE_PCRE2
   SX_TRACE ();
   SxList<SxString> res;
   matchOncePcre (input, fromByteIdx);
   if (rc_ < 1)  return res; // nothing to do
   const PCRE2_SIZE *oVector = (const PCRE2_SIZE *) findOvector ();
   for (ssize_t nCopy = rc_; nCopy > 0; nCopy--)  {
      ssize_t fromByIdx = (ssize_t) *oVector++;
      ssize_t beyond = (ssize_t) *oVector++;
      ssize_t toByteIdx = beyond - 1;
      if (toByteIdx < fromByIdx)  { // especially can be < 0
         res << (input.isUnicode () ? SxString::unicodeFromUtf8 (NULL) :
            SxString());
      }
      else  res << input.subStringByBytes (fromByIdx, toByteIdx);
   }
   return res;
#else
   needPcre2 ();
   SxList<SxString> empty;
   return empty;
#endif
}

SxList<SxString> SxRegex::match (const SxString &input) const
{
   SX_TRACE ();
   return matchOnceFrom (input, (ssize_t ) 0);
}

SxList<SxList<SxString> > SxRegex::matchAll (const SxString &input) const
{
#ifdef USE_PCRE2
   SX_TRACE ();
   SxList<SxList<SxString> > res;
   ssize_t byteOffset = 0;
   ssize_t inputNBytes = input.getNBytes ();
   bool isPosix = ((metaOptions & POSIX) != 0);
   while (1)  {
      SxList<SxString> r = matchOnceFrom (input, byteOffset);
      if (r.getSize () == 0)  break; // no (further) match found
      res << r;
      const PCRE2_SIZE *oVector =
         (const PCRE2_SIZE *) findOvector ();
      ssize_t posBehind = (ssize_t) oVector[isPosix ? (2 * rc_ - 1) : 1];
      if ( (posBehind == (ssize_t) PCRE2_UNSET) || (posBehind <= byteOffset) ) {
         // match was empty
         if (byteOffset >= inputNBytes)  break; // reached end of input
         byteOffset++;
      }
      else  byteOffset = posBehind; // continue _behind_ the entire match
   }
   return res;
#else
   needPcre2 ();
   SxList<SxList<SxString> > empty;
   return empty;
#endif
}

SxMap<SxString,SxString> SxRegex::matchToMap (const SxString &input) const
{
   SxMap<SxString,SxString> result;

   SxList<SxString> substrings = matchOnceFrom (input, 0);

   ssize_t i = 0;
   SxList<SxString>::ConstIterator it = substrings.begin ();
   while (it.isValid ())  {
      SxString number(i);
      result(number) = (*it);

      // --- named capture group
      if (numberToName.containsKey (i))
         result(numberToName(i)) = (*it);
      
      ++i;
      ++it;
   }

   return result;
}

ssize_t SxRegex::matchToOffsets (const SxString &input,
                                 ssize_t        fromCharIdx,
                                 SxArray<Match> *capturedSubstrings_) const
{
#ifdef USE_PCRE2
   SX_TRACE ();

   ssize_t inputNBytes = input.getNBytes ();
   bool isUni = input.isUnicode ();
   const char *text = input.getElems ();
   SxConstChar inp(text, input.getNBytes (), isUni);
   SxConstChar::Iterator it(&inp);
   it.setCharIdx (fromCharIdx);
   ssize_t fromByteIdx = it.getByteIdx ();
   if ( (fromByteIdx < 0) || (fromByteIdx >= inputNBytes) )  return -1;
   SxList<SxString> subs = matchOnceFrom (input, fromByteIdx);
   ssize_t nSubs = subs.getSize ();
   if ( (nSubs == 0) || (!capturedSubstrings_) )  return nSubs; // done
   ssize_t oldSize = capturedSubstrings_->getSize ();
   if (nSubs > oldSize) // need more space
      capturedSubstrings_->resize (nSubs);

   // --- copy the offsets of captured substrings; must convert byte offsets to
   // character offsets!
   const PCRE2_SIZE *oVector = (const PCRE2_SIZE *) findOvector ();
   for (ssize_t i = 0; i < nSubs; ++i)  {
      Match &m = (*capturedSubstrings_)(i);
      ssize_t fromByIdx = (ssize_t) *oVector++;
      ssize_t beyond = (ssize_t) *oVector++;
      if (beyond <= fromByIdx)  continue; // invalid/empty
      if (isUni)  {
         it.setByteIdx (fromByIdx);
         ssize_t fromChIdx = it.getCharIdx ();
         if (fromChIdx < 0)  return -1;
         ssize_t nChars;
         if (beyond > fromByIdx)  {
            ssize_t pos = beyond;
            while (!SxConstChar::isUtf8CharStart (text[--pos]))  { /* loop */ }
            it.setByteIdx (pos);
            ssize_t toChIdx = it.getCharIdx ();
            if (toChIdx < fromChIdx)  return -1;
            nChars = toChIdx - fromChIdx + 1;
         }  else nChars = 0;
         m.set (fromChIdx, fromByIdx, nChars, beyond - fromByIdx);
      }  else m.set (fromByIdx, beyond - fromByIdx); // simple case
   }

   return nSubs;
#else
   needPcre2 ();
   return 0;
#endif
}

#ifndef NDEBUG
static bool sxRegexTest (const SxString &pattern_,
                         const SxString &input,
                         const SxString &result_)
{
   SxString result;
   try {
      result = SxString::join (SxRegex(pattern_).match (input), ",");
      //cout << "expected '" << result_ << "':" << result_.getSize () << endl;
      //cout << "got      '" << result << "':" << result.getSize () << endl;
   } catch (SxException e) {
      e.print ();
      return false;
   }

   return (result == result_);
}
#endif

void SxRegex::test ()
{
#ifdef USE_PCRE2
   SX_CHECK (PCRE2_CODE_UNIT_WIDTH == 8, PCRE2_CODE_UNIT_WIDTH);
#endif

   // --- simple
   SX_CHECK (sxRegexTest ("([a])", "ab", "a,a"));

   // --- no capture groups, the overall matched string
   SX_CHECK (sxRegexTest ("[a]", "a", "a"));

   // --- numerical capture group
   SX_CHECK (sxRegexTest ("([a]).\\1", "aba", "aba,a"));

   // --- named capture group
   SX_CHECK (sxRegexTest ("(?<foo>[a]).(?P=foo)", "aba", "aba,a"));

   // --- match all
   SxList<SxList<SxString> > matchList = SxRegex("(\\d+)\\.(\\d+)")
                             .matchAll ("55.117,14.681");
   SX_CHECK (matchList.getSize () == 2, matchList);
   SxList<SxString> first = matchList.first ();
   SxList<SxString> last = matchList.last ();
   SX_CHECK (first.getSize () == 3
             && first(0) == "55.117"
             && first(1) == "55"
             && first(2) == "117", first);
   SX_CHECK (last.getSize () == 3
             && last(0) == "14.681"
             && last(1) == "14"
             && last(2) == "681", last);

   // --- match map
   SxMap<SxString,SxString> map;
   SX_CHECK ((map = SxRegex("(?<day>\\d+)\\s+(?<month>\\w+)")
                    .matchToMap ("18 March 2015")).getSize () == 5
             && map.containsKey ("0") && map("0") == "18 March"
             && map.containsKey ("1") && map("1") == "18"
             && map.containsKey ("2") && map("2") == "March"
             && map.containsKey ("day") && map("day") == "18"
             && map.containsKey ("month") && map("month") == "March");

   // --- matchToOffsets
   SxArray<SxRegex::Match> substrings;
   SX_CHECK (SxRegex("(\\d+)\\s+(?<month>\\w+)")
             .matchToOffsets ("18 March 24 June", 0, &substrings) == 3
             && substrings.getSize () > 2
             && substrings(0).getOffset () == 0
             && substrings(0).getLength () == 8
             && substrings(1).getOffset () == 0
             && substrings(1).getLength () == 2
             && substrings(2).getOffset () == 3
             && substrings(2).getLength () == 5);

   SX_CHECK (SxRegex("(\\d+)\\s+(?<month>\\w+)")
             .matchToOffsets ("18 March 24 June", 9, &substrings) == 3
             && substrings.getSize () > 2
             && substrings(0).getOffset () == 9
             && substrings(0).getLength () == 7
             && substrings(1).getOffset () == 9
             && substrings(1).getLength () == 2
             && substrings(2).getOffset () == 12
             && substrings(2).getLength () == 4);
}

SxString SxRegex::patternNotCompiledError ()
{
   return "The regex is not compiled.";
}

SxString SxRegex::patternCompileError (const SxString &pattern_,
                                       const SxString &options_)
{
   if (options_ != "")  {
      return "Cant't compile regex pattern '" + pattern_
             + "' with options '" + options_ + "'";
   } else {
      return "Cant't compile regex pattern '" + pattern_ + "'";
   }
}

SxString SxRegex::unknownOptionError (const SxString &option_,
                                      const SxString &validOptions_)
{
   return "Unknown regex option '" + option_ + "'. "
          "Available options are [" + validOptions_ + "].";
}

SxString SxRegex::patternCompileError (const SxString &pattern_,
                                       int            col_,
                                       const SxString &errmsg_)
{
   return "Pattern compilation error in '" + pattern_ 
          + "', col=" + SxString(col_) + ": " + errmsg_;
}

SxString SxRegex::matchDataCreateError ()
{
   return "Cant't create match data for the pattern.";
}

SxString SxRegex::patternInfoError (const SxString &infoName_)
{
   return "Can't get the pattern information " + infoName_ + ".";
}

SxString SxRegex::captureGroupOutOfRangeError (ssize_t number_,
                                               ssize_t captureCount_)
{
   return "Capture group " + SxString(number_) + " is out of range. "
          "The recognized number of capture groups in the pattern is "
          + SxString(captureCount_) + ".";
}

SxString SxRegex::namedCaptureGroupNotFoundError (const SxString &name_,
                                              const SxMap<SxString,ssize_t> &m)
{
   return "Named backreference '" + name_ + "' not found in the pattern. "
          "Available names are [" + SxString::join (m.getKeys (), ", ") + "].";
}

SxString SxRegex::matchingError (int errcode)
{
   return "pcre2_match error " + SxString(errcode) + ".";
}


