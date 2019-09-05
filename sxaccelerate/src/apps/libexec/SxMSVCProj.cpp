#include <SxCLI.h>
#include <SxMap.h>
#include <SxRegex.h>
#include <SxDir.h>
#include <SxFileIO.h>
#include <SxFSAction.h>
#include <SxUniqueList.h>

class SxMarker
{
   public:
      typedef SxMap<SxString,SxString>::Iterator      Iterator;
      typedef SxMap<SxString,SxString>::ConstIterator ConstIterator;

      SxMarker (const SxString &border_="*") {
         // for consistency: should be ***marker*** or ###marker###
         SX_CHECK (border_.getSize() >= 1, border_.getSize());

         borderDef = (border_.getSize() == 1) ? border_ + border_ + border_
                                              : border_;
      }

      SxMarker (const SxMap<SxString,SxString> &markers_,
                const SxString &border_="*") : markers (markers_) {
         // for consistency: should be ***marker*** or ###marker###
         SX_CHECK (border_.getSize() >= 1, border_.getSize());

         borderDef = (border_.getSize() == 1) ? border_ + border_ + border_
                                              : border_;
      }

      virtual ~SxMarker () {
         // empty
      }


      virtual SxString function (const SxString &funcId,
                                 const SxString &) const {
         return "***LOOKUP-ERROR:" + funcId + "() undefined***";
      }

      // set new marker
      SxString &operator() (const SxString &key) {
         return markers(key);
      }

      // read marker
      SxString  operator() (const SxString &key) const {
         return markers(key);
      }

      // apply marker
      SxString apply (const SxString &in, const SxString &border_="") const {
         SxString border = borderDef;
         if (border_ != "")  border = border_;

         SxMap<SxString,SxString> m = markers;
         long tmpIdx = 1;

         SxString res = in;
         SxString marker, value, expr, pattern, markerStr;
         SxList<SxString>::Iterator keyIt;
         SxList<double> numVals;
         bool err = false;
         // --- replace unresolved symbols with default values
         SxList<SxString> matches, terms;
         float markerVal, exprVal;
         bool found = true;
         while (found)  {
            found = false;

            // --- replace markers with values
            SxList<SxString> keys = m.getKeys();
            for (keyIt=keys.begin(); keyIt != keys.end(); ++keyIt)  {
               marker = *keyIt;
               value  = m(marker);
               res    = res.substitute (border + marker + border, value);
               err    = false;
               double d = value.toDouble (&err);
               if (err)  numVals << -1.;
               else      numVals << d;
            }

            // --- strip borders
            marker   = res.right(border).left(border);
            pattern  = border + marker + border;
            if (marker.getSize() == 0)  { found = false; break; }

            // --- do not replace error messages
            if (pattern.contains (border+"ERROR:"))        { found = false; break; }
            if (pattern.contains (border+"LOOKUP-ERROR:")) { found = false; break; }

            // --- ***function(marker)***
            try  {
               matches = SxRegex("^([A-Za-z0-9_]+)\\((.*)\\)", "sP").match (marker);
            }
            catch (SxException e) { e.print (); SX_EXIT; }
            if (matches.getSize() == 4)  { 
               SxString funcId = matches(1);
               SxString val = function (funcId, matches(2));
               if (!val.contains (border+"LOOKUP-ERROR:"))  {
                  SxString tmpKey = "_TMP_MARKER_" + SxString(tmpIdx++);
                  m(tmpKey) = val;
                  res = res.substitute (funcId+"("+matches(2)+")", tmpKey);
                  found = true;
                  continue;
               }  else  {
                  res = res.substitute (pattern, val);
               }
            }


            // --- ***expr?if_true:if_false***
            try { matches = SxRegex("([^\\?]+)\\?([^:]*):(.*)", "sP").match (marker);
               if (matches.getSize() == 5)  { 
                  expr  = matches(1);

                  // expr<=>value
                  terms = SxRegex("([^!<=>]+)([!<=>]+)([^!<=>]+)", "sP").match (expr);
                  if (terms.getSize() == 5)  {
                     if (!m.containsKey(terms(1))) markerStr = "";
                     else                          markerStr = m(terms(1));

                     SxList<SxString> tokens =
                        SxRegex("\"([^\"]*)\"", "sP").match (terms(3));
                     if (tokens.getSize() == 3)  terms(3) = "'" + tokens(1) + "'";
                     tokens = SxRegex("'([^']*)'", "sP").match (terms(3));
                     if (tokens.getSize() == 3)  { // string comparison
                        terms(3) = tokens(1);

                        if (terms(2) == "==")  {
                           value = (markerStr == terms(3))
                              ? matches(2)
                              : matches(3);
                        } else if (terms(2) == "!=")  {
                           value = (markerStr != terms(3))
                              ? matches(2)
                              : matches(3);
                        }
                     } else  {                     // numerical comparison
                        markerVal = markerStr.toFloat();
                        exprVal   = terms(3).toFloat();
                        if (terms(2) == "<")
                           value = (markerVal < exprVal) ? matches(2) : matches(3);
                        else if (terms(2) == "<=")
                           value = (markerVal <= exprVal) ? matches(2)
                              : matches(3);
                        else if (terms(2) == ">=")
                           value = (markerVal >= exprVal) ? matches(2)
                              : matches(3);
                        else if (terms(2) == ">")
                           value = (markerVal > exprVal) ? matches(2) : matches(3);
                        else if (terms(2) == "==")
                           value = (fabs(markerVal - exprVal) < 1e-10) ? matches(2)
                              : matches(3);
                        else
                           SX_EXIT;
                     }
                     value = value.substitute ("[#]", markerStr);
                     res   = res.substitute(pattern, value);
                     found = true;
                     continue;
                  }  else  {
                     // check only if marker exists: ***marker?if_true:if_false***
                     value = (m.containsKey(expr)) 
                        ? matches(2) 
                        : matches(3);
                     value = value.substitute ("[#]", markerStr);
                     res   = res.substitute(pattern, value);
                     found = true;
                     continue;
                  }
               } 
            } catch (SxException e) { e.print (); SX_EXIT; }

            // --- ***marker:pre:default:post***
            try { matches = SxRegex("(.+):(.*):(.*):(.*)", "sP").match (marker); } 
            catch (SxException e) { e.print (); SX_EXIT; }
            if (matches.getSize() == 6)  { 
               if (m.containsKey(matches(1)))
                  value = matches(2) 
                     + m(matches(1)) 
                     + matches(4);
               else
                  value = matches(3);
               res = res.substitute(pattern, value);
               found = true;
               continue;
            }

            // --- ***marker:pre:default***
            try { matches = SxRegex("(.+):(.*):(.*)", "sP").match (marker); }
            catch (SxException e) { e.print (); SX_EXIT; }
            if (matches.getSize() == 5)  { 
               if (m.containsKey(matches(1)))
                  value = matches(2) + m(matches(1));
               else
                  value = matches(3);
               res = res.substitute(pattern, value);
               found = true;
               continue;
            }


            // --- ***marker:default***
            try { matches = SxRegex("(.+):(.+)", "sP").match (marker); } 
            catch (SxException e) { e.print (); SX_EXIT; }
            if (matches.getSize() == 4)  { 
               value = (m.containsKey(matches(1)))
                  ? m(matches(1))
                  : matches(2);
               res = res.substitute(pattern, value);
               found = true;
               continue;
            }

            // --- ***marker[id]***
            try { matches = SxRegex("(.+)\\[(.+)\\]", "sP").match (marker); } 
            catch (SxException e) { e.print (); SX_EXIT; }
            if (matches.getSize() == 4)  { 
               if (!m.containsKey(matches(1)))  {
                  res = res.substitute(pattern, "LOOKUP-ERROR:UNKNOWN_MARKER");
                  found = true;
                  continue;
               }
               value = m(matches(1));
               err = false;
               int i = matches(2).toInt (&err);
               SxList<SxString> tokens = value.tokenize (',');
               if (abs(i) >= tokens.getSize() )  {
                  res = ""; //res.substitute(pattern, "LOOKUP-ERROR:BOUNDARY");
                  found = true;
                  continue;
               }
               ssize_t idx = (i >= 0) ? i : tokens.getSize() + i;
               res = res.substitute(pattern, tokens(idx));
               found = true;
               continue;
            }

         }

         return res;
      }

      Iterator      begin () {
         return markers.begin ();
      }
      ConstIterator begin () const {
         return markers.begin ();
      }
      Iterator      end () {
         return markers.end ();
      }
      ConstIterator end () const {
         return markers.end ();
      }

      SxList<SxString> getKeys () const {
         return markers.getKeys ();
      }
      SxList<SxString> getValues () const {
         return markers.getValues ();
      }

      void append (const SxMarker &in) {
         markers.append (in.markers);
      }

      bool containsKey (const SxString &key) const {
         return markers.containsKey (key);
      }
      static SxList<SxString> getUnresolved (const SxString &tmpl, 
            const SxString &border_="*",
            ssize_t maxSize=20) {
         SxUniqueList<SxString> res;

         SxString delim;
         if (border_ == "*")  {
            delim = "\\*";
         } else {
            SX_EXIT;  // delimiter not supported
         }
         SxString border = delim + delim + delim;
         //               \\*\\*\\*([^  \\*    ]+)(\\*\\*\\*)
         SxString expr   = border+"([A-Z1-9_]{1,"+maxSize+"})("+border+")";

         SxList<SxString> matches;
         bool found = true;
         SxString t = tmpl, id;
         while (found)  {
            // --- ***marker***
            try { matches = SxRegex(expr, "sP").match (t); } 
            catch (SxException e) { e.print (); SX_EXIT; }
            if (matches.getSize() == 4)  {
               SxString m = matches(1);
               t = matches.last();  // $'

               // ***expr?true:false*** || ***=expr*** || ***function(val)***
               if (   m.contains ("?") || m.contains ("=")
                   || m.contains ("(")) 
               {
                  continue;
               }

               if ( (id = m.left(":")) != "" )  {
                  res << id;  // ***marker:value***
               } else if ( (id = m.left("[")) != "" )  {
                  res << id;  // ***marker[id]***
               } else {
                  res << m;   // ***marker***
               }

               continue;
            }
            found = false;
         }
         return res;
      }

      void print () const {
         SxMap<SxString,SxString>::ConstIterator it;
         for (it = markers.begin(); it != markers.end(); ++it)  {
            cout << it.getKey () << " = " << it.getValue () << endl;
         }
      }

      //protected:

      SxString borderDef;
      SxMap<SxString,SxString> markers;

};

SxList<SxArray<ssize_t> >
getBlockIdx (const SxString &str,
             const SxString &begin_,
             const SxString &end_,
             bool inclTokens_)
{
   SX_CHECK (!(str.isUnicode () || begin_.isUnicode () || end_.isUnicode ()));
   ssize_t beginOff=0, endOff=0;
   if (inclTokens_)  endOff = end_.getSize();
   else { beginOff = begin_.getSize(); endOff = -1; }

   SxArray<ssize_t> beginIdx = str.findAll (begin_);
   SxArray<ssize_t> endIdx   = str.findAll (end_);
   SX_CHECK (beginIdx.getSize() == endIdx.getSize(),
             beginIdx.getSize(),   endIdx.getSize());

   ssize_t nTokens = beginIdx.getSize();
   ssize_t b, e, i, j;
   SxArray<SxArray<ssize_t> >  blocks (nTokens);
   SxArray<ssize_t> indicies  (2 * nTokens);
   SxArray<bool>    isEnd (2 * nTokens); // false: begin, true: end

   // --- serialize begin and end indicies
   ssize_t bIdx = 0, eIdx = 0, idx = 0;
   ssize_t threshold = -1;
   while (idx < 2*nTokens)  {
      if (bIdx < nTokens)  b = beginIdx(bIdx);
      else                 b = str.getSize(); // just a large value
      if (eIdx < nTokens)  e = endIdx(eIdx);
      else                 e = str.getSize(); // just a large value
      SX_CHECK_3VARS (b < str.getSize() || e < str.getSize(), b, e, str.getSize());
      if (b > threshold && e > threshold)  {
         SX_CHECK (b != e, b);
         if (b < e)  {
            indicies(idx) = b;
            isEnd(idx) = false;
            idx++; bIdx++;
            threshold = b;
         } else {
            indicies(idx) = e;
            isEnd(idx) = true;
            idx++; eIdx++;
            threshold = e;
         }
      } else {
         if (b < threshold)  {
            indicies(idx) = b;
            isEnd(idx) = false;
            idx++; bIdx++;
            threshold = b;
         }
         if (e < threshold)  {
            indicies(idx) = e;
            isEnd(idx) = true;
            idx++; eIdx++;
            threshold = e;
         }
      }
   }

   SxList<SxArray<ssize_t> > indexList, res;
   int lvl=-1, jLvl=-1;
   for (i=0; i < 2*nTokens; ++i)  {
      // --- get next begin
      b = -1;
      for (j=i; j < 2*nTokens; ++j)  {
         if ( !isEnd(j) )  { b = indicies(j); ++lvl; break; }
      }
      if (b == -1)   break;  // no more blocks found

      // --- find matching end
      e = -1;
      jLvl = lvl;
      for (j=i+1; j < 2*nTokens; ++j)  {
         if (!isEnd(j))  ++jLvl;   // nested begin
         else {
            if (jLvl == lvl)  {    // end
               e = indicies(j);
               break;
            }  else  {
               --jLvl;
               SX_CHECK (jLvl >= lvl, jLvl, lvl);
            }
         }
      }
      SX_CHECK (e != -1);  // no matching END found
      lvl = jLvl;

      indexList << (SxArray<ssize_t> (SxList<ssize_t>() << (b + beginOff)
               << (e + endOff)));
   }


   SxList<SxArray<ssize_t> >::ConstIterator it;
   for (it = indexList.begin(); it != indexList.end(); ++it)  {
      b = (*it)(0);
      e = (*it)(1);
      if (!inclTokens_)  {
         bool whitespace = true;
         if (b < str.getSize()-1 && (str)(b) == '\n')   b++;
         whitespace = true;
         while (whitespace && e >  b)  {
            char c = (str)(e);
            if      (c == ' ' || c == '\t')  { e--; continue; }
            else if (c == '\n')              { e--; whitespace = false; }
            else                             { whitespace = false; }
         }
      }
      SxArray<ssize_t> arr(2); arr(0) = b; arr(1) = e;
      res << arr;
   }
   return res;
}

SxList<SxString>
getBlocks (const SxString &str,
           const SxString &begin_,
           const SxString &end_,
           bool inclTokens_)
{
   SxList<SxString> res;
   SxList<SxArray<ssize_t> > indicies = getBlockIdx (str, begin_, end_, inclTokens_);
   SxList<SxArray<ssize_t> >::ConstIterator it;
   for (it = indicies.begin(); it != indicies.end(); ++it)  {
      ssize_t b = (*it)(0);
      ssize_t e = (*it)(1);
      res << str.subString (b, e);
   }
   return res;
}

SxString substBlocks (const SxString &str,
                      const SxMap<SxString,SxString> &markers,
                      bool inclTokens,
                      const SxString &leftBorder,
                      const SxString &rightBorder)
{
   SxString res = str;

   SxList<SxArray<ssize_t> > ranges;
   SxList<SxArray<ssize_t> >::ConstIterator it;

   const SxList<SxString> &keys = markers.getKeys();
   SxList<SxString>::ConstIterator keyIt;
   // --- replace markers with values
   for (keyIt=keys.begin(); keyIt != keys.end(); ++keyIt)  {
      SxString marker = *keyIt;
      SxString value  = markers(marker);
      SxString begin_ = leftBorder + " BEGIN " + marker + " " + rightBorder;
      SxString end_   = leftBorder + " END "   + marker + " " + rightBorder;
      ranges = getBlockIdx (res, begin_, end_, inclTokens);
      for (it = ranges.fromLast(); it != ranges.toFirst(); --it)  {
         res = res.subString (0, (*it)(0)-1)
            + value
            + res.subString ((*it)(1));
      }
   }
   return res;
}

SxString getBlock (const SxString &tmpl, const SxString &id)
{
   SxList<SxString> blocks = getBlocks (tmpl, "<!-- BEGIN " + id + " -->",
         "<!-- END " + id + " -->", false);
   SX_CHECK_VAR (blocks.getSize() == 1, blocks);
   return blocks.first ();
}

int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   SxString cwd    = cli.option ("--cwd", "Path from where file list should be "
         "generated").toString (".");
   SxString projName = cli.option ("--name", "project name").toString ();
   SxString logHash  = cli.option ("--hash", "project hash id").toString ();
   SxString depth  = cli.option ("--depth", "Difference to top level tree").toString ();
   SxString topsrc = cli.option ("--srctop", "Top level of source tree").toString ();
   SxString accltop= cli.option ("--sxaccltop", "Top level of SxAccelerate").toString ("");
   SxString deps   = cli.option ("--deps", 
         "List of dependencies"
         ).toString ("").substitute(" ",";").substitute(",", ";");
   SxString mode   = cli.option ("--mode", "Compilation mode").toString ();
   cli.finalize ();

   if (topsrc.contains ("/src"))  topsrc = depth + topsrc.right("/src");

   if (accltop == "")  accltop = "../" + topsrc;
   if (accltop.contains ("/src"))  accltop = depth + accltop.right("/src");

   SxList<SxString> modes;
   SxList<SxString>::ConstIterator mIt;

   SxString configType;
   if (mode == "withApps")  modes << "libOnly" << "withBin" << "withSBin" << "withLibExec";
   else                     modes << mode;

   for (mIt = modes.begin(); mIt != modes.end(); ++mIt)  {

      // --- dependencies
      SxList<SxString> tokens, depLibs, depIncs;
      tokens = deps.tokenize(';');
      SxList<SxString>::ConstIterator tIt;
      for (tIt = tokens.begin(); tIt != tokens.end(); ++tIt)  {
         SxString entry = tIt->substitute("_SolutionDir_", "$(SolutionDir)")
            .substitute("_IntDir_",      "$(IntDir)")
            .substitute("_OutDir_",      "$(OutDir)")
            .substitute("/", "\\");

         if (entry.contains ("!"))  {
            // explicit CINCLUDE / LIBPATH specified
            depIncs << entry.left ("!");
            depLibs << entry.right ("!");
         }  else  {
            // external library?
            if (entry.contains (".lib") || entry.contains (".a"))  {
               depLibs << *tIt;
            }  else  {
               cout << entry << endl;
               depIncs << ( entry );
               depLibs << ( entry.tokenize('\\').last() + ".lib" );
            }
         }
      }
      depIncs << ( "$(SolutionDir)\\src\\windows" );

      SxString macros = "";
      if      (*mIt == "libOnly")      configType = "DynamicLibrary";
      else if (*mIt == "withBin")      configType = "Application";
      else if (*mIt == "withSBin")     configType = "Application";
      else if (*mIt == "withLibExec")  configType = "Application";
      else {
         cout << "ERROR: mode " << mode << " is not supported.\n";
         SX_QUIT;
      }

      if (mode == "withApps" && (*mIt == "withBin" || *mIt == "withSBin" || *mIt == "withLibExec"))  {
         macros = "SX_STANDALONE";
         depLibs.prepend ( projName + ".lib" );
      }


      try { SxFSAction::cd (cwd); } catch (SxException e) { e.print (); SX_QUIT; }

      SxString tmplPath = SxDir::getExecPath() + "/../system/sxtemplate.vcxproj";
      // support --with-bintarget
      if (!SxFSAction::test_f (tmplPath))
         tmplPath = SxDir::getExecPath() + "/../../system/sxtemplate.vcxproj";
      // fallback to source tree
      if (!SxFSAction::test_f (tmplPath))
         tmplPath = SxString(SXACCELERATE_SRC) + "/system/sxtemplate.vcxproj";

      if (!SxFSAction::test_f (tmplPath))
         SX_EXIT;


      SxString tmpl = SxFileIO::readLines (tmplPath);

      SxString headerTmpl = getBlock (tmpl, "Header");
      SxString sourceTmpl = getBlock (tmpl, "Source");
      SxString lexTmpl    = getBlock (tmpl, "Lex");
      SxString yaccTmpl   = getBlock (tmpl, "Yacc");
      SxString winEnvTmpl = getBlock (tmpl, "MSVC_Env_win64d");

      // --- inspect file structure
      SxList<SxString> headers, sources, lexes, yaccs;
      SxList<SxFileInfo> headerItems, sourceItems, lexItems, yaccItems;
      try { 
         headerItems = SxFSAction::ls ("*.h"); 
         sourceItems = SxFSAction::ls ("*.cpp"); 
         lexItems    = SxFSAction::ls ("*.lpp"); 
         yaccItems   = SxFSAction::ls ("*.ypp"); 
      }
      catch (SxException e)  { e.print (); SX_EXIT; }

      SxList<SxFileInfo>::ConstIterator fIt;

      if (configType == "DynamicLibrary")  {

         for (fIt = headerItems.begin(); fIt != headerItems.end(); ++fIt)
            headers << headerTmpl.substitute ("***Header.filename***", fIt->getName ());
         for (fIt = sourceItems.begin(); fIt != sourceItems.end(); ++fIt)
            sources << sourceTmpl.substitute ("***Source.filename***", fIt->getName ());
         for (fIt = lexItems.begin(); fIt != lexItems.end(); ++fIt)  {
            lexes << lexTmpl.substitute ("***Lex.filename***", fIt->getName ());
            SxString tabCpp = sourceTmpl.substitute ("***Source.filename***",
                  "$(IntDir)\\"+fIt->getName().substitute (".lpp", ".tab.cpp"));
            if (!sources.contains (tabCpp))  sources << tabCpp;
         }
         for (fIt = yaccItems.begin(); fIt != yaccItems.end(); ++fIt)  {
            yaccs << yaccTmpl.substitute ("***Yacc.filename***", fIt->getName ());
            SxString yyCpp = sourceTmpl.substitute ("***Source.filename***",
                  "$(IntDir)\\"+fIt->getName().substitute (".ypp", ".yy.cpp"));
            if (!sources.contains (yyCpp))   sources << yyCpp;
         }

         // --- build blocks
         SxMap<SxString,SxString> blocks;
         blocks("Header")   = SxString::join (headers, "\r\n");
         blocks("Source")   = SxString::join (sources, "\r\n");
         blocks("Lex")      = SxString::join (lexes,   "\r\n");
         blocks("Yacc")     = SxString::join (yaccs,   "\r\n");
         SxString msvcEnv = projName == "sxutil" ? winEnvTmpl : "";
         blocks("MSVC_Env_win32d") = msvcEnv;
         blocks("MSVC_Env_win32r") = msvcEnv;
         blocks("MSVC_Env_win32p") = msvcEnv;
         blocks("MSVC_Env_win64d") = msvcEnv;
         blocks("MSVC_Env_win64r") = msvcEnv;
         blocks("MSVC_Env_win64p") = msvcEnv;
         tmpl = substBlocks (tmpl, blocks, true, "<!--", "-->");

         SxMarker marker;
         marker("NAME")        = projName;
         marker("SX_LOG_ID")   = projName;
         marker("SX_LOG_HASH") = logHash;
         marker("DEPTH")       = depth.substitute("/", "\\");
         marker("TOPSRC")      = depth.substitute("/", "\\");
         marker("SXACCLTOP")   = accltop.substitute("/", "\\");
         marker("MACROS")      = macros;
         marker("DEPINCS")     = SxString::join (depIncs, ";").substitute("/", "\\");
         marker("DEPLIBS")     = SxString::join (depLibs, ";");
         marker("CONFIG_TYPE") = configType;
         marker("CONFIG_TYPE") = configType;

         SxString out = marker.apply (tmpl);

         // --- remove \r\n at the end
         ssize_t nChars = 0;
         if (out(out.getSize()-1)        == '\n')  nChars++;
         if (out(out.getSize()-nChars-1) == '\r')  nChars++;
         out = out(0,out.getSize()-nChars-1);
         SxFileIO::write (out, cwd + "/" + projName + ".vcxproj", 0644);

      }  else  {

         for (fIt = sourceItems.begin(); fIt != sourceItems.end(); ++fIt)  {
            SxString project = fIt->getName ().substitute (".cpp", "");
            SxString srcName = sourceTmpl.substitute ("***Source.filename***", fIt->getName ());

            // --- build blocks
            SxMap<SxString,SxString> blocks;
            blocks("Header") = "";
            blocks("Source") = srcName;
            blocks("Lex")    = "";
            blocks("Yacc")   = "";

            SxMarker marker;
            marker("NAME")    = projName;
            marker("DEPTH")   = depth.substitute("/", "\\");
            marker("TOPSRC")  = depth.substitute("/", "\\");
            marker("SXACCLTOP") = accltop.substitute ("/", "\\");
            marker("MACROS")  = macros;
            marker("DEPINCS") = SxString::join (depIncs, ";").substitute("/", "\\");
            marker("DEPLIBS") = SxString::join (depLibs, ";");
            marker("CONFIG_TYPE")  = configType;

            SxString out = marker.apply (
                  substBlocks (tmpl, blocks, true, "<!--", "-->")
                  );

            // --- remove \r\n at the end
            ssize_t nChars = 0;
            if (out(out.getSize()-1)        == '\n')  nChars++;
            if (out(out.getSize()-nChars-1) == '\r')  nChars++;
            out = out(0,out.getSize()-nChars-1);

            SxFileIO::write (out, cwd + "/" + project + ".vcxproj", 0644);
         }
      }
   }
   return 0;
}
