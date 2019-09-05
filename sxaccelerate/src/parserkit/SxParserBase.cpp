#include <fstream>
#include <SxParserBase.h>
#include <SxFSAction.h>

// --------------------------------------------------------------------------

SxParserBase::SxParserBase ()
   : scannerPtr(NULL),
     inStream(NULL),
     lineOffset(0),
     traces(false),
     curLineIndent(-1),
     prevLineIndent(-1),
     maxIncludeDepth(10)
{
#  ifndef NDEBUG
      if (::getenv ("SX_DEBUG_PARSER"))
        enableTraces (true);
#  endif /* NDEBUG */
}

SxParserBase::~SxParserBase ()
{
   // empty
}

bool SxParserBase::resume () const
{
   return (errors.getSize() < maxErrors);
}

void SxParserBase::errorHandler (const SxString &msg_,
                                 ssize_t line0, ssize_t col0,
                                 ssize_t line1, ssize_t col1)
{
   SX_TRACE ();
   // -- display non-printable characters
   SxString err = msg_.substitute ("\n","\\n");
   err = SxSed("s#(\t)#\\\\t#g").subst (err);
   err = SxSed("s#(\n)#\\\\n#g").subst (err);
   err = SxSed("s#(\a)#\\\\a#g").subst (err);
   err = SxSed("s#(\b)#\\\\b#g").subst (err);
   err = SxSed("s#(\r)#\\\\r#g").subst (err);

   // -- convert TK_TOKEN to Token
   err = SxSed("s#(TK_)([A-Za-z]+)#\\L$1\\L$2#g").subst (err);
   err = SxSed("s#(tk_)([a-z])#\\u$2#g").subst (err);
   SxString msg = inFile + ":"
                + line0 + "." + col0 + "-" + line1 + "." + (col1-1)
                + ": ";
   msg += err;
   errors << Error (line0, col0, line1, col1, msg);
}

void SxParserBase::setLineOffset (ssize_t offset)
{
   SX_TRACE ();
   lineOffset = offset;
}


void SxParserBase::setSearchPath (const SxString &path)
{
   SX_TRACE ();
   searchPath << "";  // support absolute filenames
#  ifdef WIN32
      SxString sep = ";";
#  else
      SxString sep = ":";
#  endif /* WIN32 */

   searchPath << path.tokenize(sep);
}

void SxParserBase::setMaxIncludeDepth (ssize_t depth)
{
   SX_TRACE ();
   maxIncludeDepth = depth;
}

void SxParserBase::enableTraces (bool state)
{
   SX_TRACE ();
   traces = state;
}

void SxParserBase::validateIncludes (const SxString &pattern) const
{
   SX_TRACE ();

   if (lexIncludes.getSize() > maxIncludeDepth)  {
      SX_THROW ("Too many include levels:\n- "
               + SxString::join (lexIncludes, "\n- "));
   }

   SxSortedList<SxString> fqfns = resolveFiles (pattern);
   SxString filename = fqfns.first ();

   // --- is file readable?
   FILE *fp = fopen (filename.getElems (), "r");
   if (!fp)  {
      SX_THROW ("file open failed for " + filename);
   }
   fclose (fp);
}

void SxParserBase::setIncLoc (int line0, int col0, int line1, int col1)
{
   SX_TRACE ();
   lexIncLoc = SxArray<int> (4);
   lexIncLoc(0) = line0; lexIncLoc(1) = col0;
   lexIncLoc(2) = line1; lexIncLoc(3) = col1;
}

bool SxParserBase::processNextInclude ()
{
   SX_TRACE ();
   if (lexIncStack.getSize() == 0)  return false;

   SxString nextFile = lexIncStack.first ();
   lexIncStack.removeFirst ();

   int dummy = 0;
   popInclude   (&dummy, &dummy, &dummy, &dummy);

   int line0 = lexIncLoc(0), col0 = lexIncLoc(1);
   int line1 = lexIncLoc(2), col1 = lexIncLoc(3);

   SX_DBG_MSG ("   about to open " << nextFile);

   pushIncludes (nextFile, line0, col0, line1, col1);
   SX_DBG_MSG ("returning " << (lexIncStack.getSize() > 0));
   return true;
}

void SxParserBase::setInFile (const SxString &inFile_)
{
   SX_TRACE ();
   inFile = inFile_;
}

bool SxParserBase::pushIncludes (const SxString &pattern,
                                 int line0, int col0,
                                 int line1, int col1)
{
   SX_TRACE ();
   SX_DBG_MSG ("   opening file " << pattern);

   if (lexIncludes.getSize() > maxIncludeDepth)  {
      SX_THROW ("too many include levels.");
   }

   SxSortedList<SxString> fqfns = resolveFiles (pattern);

   SxString filename = fqfns.first ();
   fqfns.removeFirst(); lexIncStack << fqfns;

   // --- is file readable?
   FILE *fp = fopen (filename.getElems (), "r");
   if (!fp)  {
      SX_THROW ("file open failed for " + filename);
   }
   fclose (fp);

   SxArray<int> loc = { line0, col0, line1, col1 };

   lexStreams  << inStream;
   lexIncludes << inFile;
   lexLocs     << loc;

   inFile   = filename;
   inStream = new ifstream (filename.getElems (), ios::in|ios::binary);
   SX_CHECK (inStream->good ());
   return true;
}

bool SxParserBase::popInclude (int *line0, int *col0, int *line1, int *col1)
{
   SX_TRACE ();

   if (lexIncludes.getSize() == 0)  return false;

   delete inStream;
   inStream = lexStreams.last ();
   inFile   = lexIncludes.last ();

   SxArray<int> loc = lexLocs.last();
   *line0 = loc(0); *col0 = loc(1); *line1 = loc(2); *col1 = loc(3);

   lexStreams.removeLast ();
   lexIncludes.removeLast ();
   lexLocs.removeLast ();

   return true;
}

void SxParserBase::pushLex (int state, int line, int col)
{
   SX_TRACE ();
   lexStates.append (state);
   lexStack.append (SxList<SxString>());
   lexLines << line; lexCols << col;
}

void SxParserBase::repushLex (int line, int col)
{
   SX_TRACE ();
   SX_CHECK (lexStates.getSize() > 0);
   lexStates.append (lexStates.last());
   lexStack.append (SxList<SxString>());
   lexLines << line; lexCols << col;
}

void SxParserBase::popLex (int *line, int *col)
{
   SX_TRACE ();
   lexStack.removeLast (); lexStates.removeLast ();
   *line = lexLines.last (); lexLines.removeLast ();
   *col  = lexCols.last (); lexCols.removeLast ();
}

int SxParserBase::getLexState ()
{
   if (lexStates.getSize() == 0)  return -1;
   return lexStates.last ();
}

void SxParserBase::appendLex (const SxString &in)
{
   SX_TRACE ();
   lexStack.last().append (in);
}

SxList<SxString> SxParserBase::collectLex ()
{
   SX_TRACE ();
   SxList<SxString> res = lexStack.last();
   lexStack.last().removeAll ();
   return res;
}

SxString SxParserBase::getLexTag ()
{
   SX_TRACE ();
   SX_CHECK (lexLines.getSize () > 0);
   SX_CHECK (lexCols.getSize () > 0);
   return inFile + ":" + SxString (lexLines.last()) + "." + lexCols.last();
}

void SxParserBase::setIndent (int i)
{
   SX_TRACE ();
   curLineIndent = i;
}

int SxParserBase::getIndent ()
{
   SX_TRACE ();
   return curLineIndent;
}

void SxParserBase::updateIndent ()
{
   SX_TRACE ();
   prevLineIndent = curLineIndent;
}

int SxParserBase::indentDiff ()
{
   SX_TRACE ();
   return curLineIndent - prevLineIndent;
}

bool SxParserBase::unindent ()
{
   SX_TRACE ();
   return (--prevLineIndent > curLineIndent);
}

SxSortedList<SxString> SxParserBase::resolveFiles (const SxString &pattern) const
{
   SX_TRACE ();
   SxSortedList<SxString> fqfns;
   bool useWildcards = (pattern.contains ("*") || pattern.contains("?"));

   SxString curPath = inFile;
   if (useWildcards && SxFile(curPath).exists())  {  // not STDIN or cin
      curPath = SxFile(curPath).getAbsPath();
   }

   SxList<SxFileInfo> items;
   SxList<SxFileInfo>::ConstIterator infoIt;
   SxList<SxString>::ConstIterator it;
   for (it = searchPath.begin(); it != searchPath.end(); ++it)  {
      SxString fqfn = *it + "/" + pattern;
      fqfn = SxString::fromUtf8 (fqfn.ascii ());
      try {
         if (useWildcards)  {
            items = SxFSAction::ls (fqfn);
            for (infoIt = items.begin(); infoIt != items.end(); ++infoIt) {
               const SxString path = infoIt->getAbsPath();
               SX_DBG_MSG ("   test -f " << path);
               if (path != curPath)
                  fqfns << path;
            }
            if (fqfns.getSize() > 0)  return fqfns;
         }  else  {
            SX_DBG_MSG ("   test -f " << fqfn);
            if (SxFile(fqfn).exists())  {
               fqfns = SxSortedList<SxString>() << fqfn;
               return fqfns;
            }
         }
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
   }
   // --- nothing found
   if (!useWildcards)  {
      SX_THROW ("file not found " + pattern);
   }
   return SxSortedList<SxString> ();
}


int SxParserBase::readString (const SxString &in, const SxString &label,
                              ssize_t maxErrors_)
{
   SX_TRACE ();
   inFile = label;
   maxErrors = maxErrors_;
   // -- in.getElems() for utf-8
   inStream = new istringstream (in.getElems ());
   initScanner (traces);
   int res = parse ();
   destroyScanner ();
   delete inStream;
   handleExceptions ();
   return res;
}

int SxParserBase::readStream (std::istream *in, const SxString &label,
                              ssize_t maxErrors_)
{
   SX_TRACE ();
   inFile = label;
   maxErrors = maxErrors_;
   inStream = in;
   initScanner (traces);
   int res = parse ();
   destroyScanner ();
   inStream = NULL;
   handleExceptions ();
   return res;
}


int SxParserBase::readFile (const SxString &filename, ssize_t maxErrors_)
{
   SX_TRACE ();
   inFile = filename;
   maxErrors = maxErrors_;
   if (!SxFSAction::test_f (filename))  {
      SX_THROW ("file not found " + filename);
   }
   // --- is file readable?
   FILE *fp = fopen (filename.getElems (), "r");
   if (!fp)  {
      SX_THROW ("file open failed " + filename);
   }
   fclose (fp);

   std::ifstream *ifStream = new ifstream (filename.getElems (), ios::in|ios::binary);
   SX_CHECK (ifStream->good ());
   inStream = ifStream;
   initScanner (traces);
   int res = parse ();
   destroyScanner ();
   ifStream->close ();
   delete inStream;
   handleExceptions ();
   return res;
}


int SxParserBase::readFunction (const ReaderFunction &cb,
                                const SxString &label, ssize_t maxErrors_)
{
   SX_TRACE ();
   inFile = label;
   maxErrors = maxErrors_;
   inStream = NULL;
   readerCB = cb;
   initScanner (traces);
   int res = parse ();
   destroyScanner ();
   handleExceptions ();
   return res;
}

void SxParserBase::handleExceptions ()
{
   SX_TRACE ();
   if (errors.getSize() > 0)  {
      SxString msg;
      SxList<Error>::ConstIterator it;
      for (it = errors.begin(); it != errors.end(); ++it)  {
#  ifdef WIN32
         msg += it->msg + "\r\n";
#  else
         msg += it->msg + "\n";
#  endif
      }
      SX_THROW (msg);
   }
}

