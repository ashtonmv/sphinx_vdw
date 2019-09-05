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

#include <SxCLI.h>
#include <sys/types.h>
#include <unistd.h>
#ifdef WIN32
#  include <process.h>
#endif
#include <math.h>
#include <SxTime.h>
#include <SxUniqueList.h>
#include <SxException.h>
#include <SxLog.h>

/// List of functions to execute when exiting
SxList<void (*)()> SxCLI::atExit;

void SxCLI::Log::disable (const SxString &why)
{
   Log &logger = get ();
   logger.enabled = false;
   if (logger.reason.getSize () > 0)
      logger.reason += "\n" + why;
   else
      logger.reason = why;
   if (logger.sxTee)  {
#ifndef NDEBUG
      cout << "In " << __FILE__ << ":" << __LINE__ << endl;
      // If you get here, your startup procedure is not
      // well designed. Log file redirection is started by SxCLI::finalize
      // (or the first CLI argument).
      // Any automatic disabling of redirection (e.g. by MPI slaves)
      // must be done before.
#endif
      cout << "WARNING: this process disables file logging because "
           << why << endl;
      cout << "However, log file redirection has already started." << endl;
      logger.sxTee->disable ();
      logger.sxTee = SxPtr<SxRedirect> ();
   }
}

// ------------------------------  CliArg -----------------------------------
SxCLI::CliArg::CliArg (SxCLI*         parentIn,
                       const SxString m,
                       const SxString s,
                       const SxString l,
                       const bool     optIn,
                       const bool     eoIn)
   : shortDescription (eoIn ? "" : s),longDescription(l),
     optional(optIn), emptyOption(eoIn), used(false), parent(parentIn)
{
   SX_CHECK (parentIn);
   stickyOption = parent->stickyDefault;
   group  = parent->currentGroup;
   SX_CHECK (group);
   if (m.getSize () == 0) return;
   // --- identify and check option marks
   if (m.contains ("|"))
      mark << m.tokenize ('|');
   else
      mark << m;

   if (group == &parent->autoOptionGroup)  {
      if (mark.getSize () > 1)  {
         for (int i = 0; i < mark.getSize (); i++)  {
            const SxString &thisMark = mark(i);
            // This is an internal option
            // --- use only option names not used before
            if (!parent->autoOptionNames.contains (thisMark))  {
               mark.remove (i);
               i--;
            }
         }
         // check for unrecoverable name clash with internal option
         SX_CHECK (mark.getSize () > 0);
      }
      // limit number of marks variants for internal options to 3
      while (mark.getSize () > 3) mark.removeLast ();
   } else {
      for (int i = 0; i < mark.getSize (); i++)  {
         const SxString &thisMark = mark(i);
         // remove used marks from available automatic option names
         if (parent->autoOptionNames.contains (thisMark))
            parent->autoOptionNames.removeElement (thisMark);
      }
   }
#ifndef NDEBUG
   // --- check for name collisions with previous options
   for (int im = 0; im < mark.getSize (); ++im)  {
      const SxString &myMark = mark(im);
      for (SxList<CliArg>::ConstIterator optIt = parent->optionList.begin ();
           optIt != parent->optionList.end (); ++optIt)
      {
         if (optIt->fits(myMark))  {
            cout << "Name collision in SxCLI. '"
                 << myMark << "' cannot be used here\nbecause ";
            optIt->printMark (false);
            cout << " might consume this option." << endl;
            cout << "Rename or reorder your options." << endl;
            SX_EXIT;
         }
      }
   }
#endif
}

int SxCLI::CliArg::fits (const SxString& s) const
{
   ssize_t length;
   if (mark.getSize () == 0) return 0;
   for (int i = 0; i < mark.getSize (); i++)  {
      length = mark(i).getSize ();
#ifdef WIN32
      if (mark(i)(0) == '-' && s(0) == '/')  {
         const SxString &m = mark(i);
         if (m(1) == '-')  {
            // --XXX => /XXX
            if (m.subString(2) == s.subString(1)) return int(length);
         } else if (s.getSize () == length)  {
            // -XXX => /XXX
            if (m.subString(1) == s.subString(1)) return int(length);
         }
      } else 
#endif
      {
         if (length > s.getSize ()) continue;
         if (s(0, length - 1) == mark(i) ) return int(length);
      }
   }
   return 0;
}

bool SxCLI::CliGroup::checkExclusions ()
{
   SX_CHECK (initOption);
   if (excludedGroups.getSize () == 0) return true;
   for (SxList<int>::Iterator it = excludedGroups.begin ();
         it != excludedGroups.end ();
         it++)
   {
      if (initOption->parent->groups(*it).status == Required)  {
         initOption->parent->exclusionError (*this, initOption->parent->groups(*it));
         return false;
      }
   }
   return true;
}

SxString SxCLI::CliArg::getArgument ()
{
   SxList<SxString> &args_ = parent->arguments(parent->iLoop);
   SX_CHECK (group);
   // no arguments at all?
   if (args_.getSize () == 0)  {
      if (!optional)  {
        if (group->status == Required)  {
           if (parent->mode == Normal)  {
              // --- error: required option not given
              printMark ();
              cout << " is required!" << endl;
           }
           parent->setError ();
        } else if (group->status == Optional)  {
           group->status = MissingRequired;
           SX_CHECK (group->initOption == NULL);
           group->initOption = this;
        }
      }
      return "";
   }

   // --- no mark? take first arg.
   if (mark.getSize () == 0)  {
      if (group->status == Optional || group->status == Required)
      {
         SxString s = args_(0);
         args_.remove (0);
         if (stickyOption && parent->iLoop + 1 < parent->arguments.getSize ())
            parent->arguments(parent->iLoop + 1).append (s);
         group->status = Required;
         return s;
      } else {
         // allow arguments in groups: silently ignored if group is forbidden
         // we can't issue an error, because we don't know if this is "our"
         // argument. So we just leave it. If it was ours, finalize () will
         // issue an "unknown option" error.
         return "";
      }
   }

   int iArg;
   int length = 0;
   SxString arg = "";
   for (iArg = 0; iArg < args_.getSize (); iArg++)  {
      
      // does mark fit?
      length = fits (args_(iArg));
      if (length == 0) continue;
      // --- mark fits arg
      arg = args_(iArg);
      args_.remove(iArg);
      if (stickyOption && !exists(parent->iLoop+1))
         parent->arguments(parent->iLoop + 1).append (arg);
      break;
   }
   if (arg == "")  {
      if (!optional && group->status == Required)  {
         // --- error: required option not given
         if (parent->mode == Normal)  {
            printMark ();
            cout << " is required!" << endl;
         }
         parent->setError ();
      }
      if (!optional && group->status == Optional)  {
         group->status = MissingRequired; // not yet an error
         SX_CHECK (group->initOption == NULL);
         group->initOption = this;
      }
      return arg; // == ""
   }
   used = true;

   // --- Option is given.
   if (!group->setRequired (this))
      return ""; // error occured
   
   if (emptyOption)  {
      if (length < arg.getSize () && arg(length) == '=')  {
         optValue = arg.subString (length); // get the =value part
         arg.resize (length, true); // return mark without =value stuff
      }
   } else {
      // --- get the argument
      if (length == arg.getSize () )  {
         // --- argument is not attached to mark
         if (iArg < args_.getSize ())  {
            arg = args_(iArg);
            args_.remove(iArg);
         } else {
            if (parent->mode == Completion)  {
               if (shortDescription.getSize () > 0)
                  cout << '<' << shortDescription << '>' << endl;
               else
                  cout << "<??\?>" << endl;
               runAtExitFuncs ();
               exit (0);
            }
            // argument is missing
            cout << "Option '" << arg << "': ";
            if (shortDescription.getSize () > 0)
               cout << shortDescription;
            else 
               cout << "Argument";
            cout << " is missing." << endl;
            parent->setError (); 
            return "";
         }
      } else { 
         // --- argument is attached to mark
         if (arg(length) == '=') // allow for --option=value style
            arg = arg.subString (length + 1);
         else
            arg = arg.subString (length);
      }
   }
   return arg;
}

namespace {
   void printString (const SxString &var, const SxString &txt)
   {
      cout << var << " = \"";
      for (int i = 0; i < txt.getSize (); ++i)  {
         char c = txt(i);
         if (c == '\n')  {
            cout << "\\n\"" << endl;
            for (int j = 0; j < var.getSize (); j++) cout << ' ';
            cout << " + \"";
         } else if (c == '"')
            cout << "''";
         else
            cout << c;
      }
      cout << '"' << ';' << endl;
   }
}

void SxCLI::printout () const
{
   cout << "format cli;" << endl;
   cout << "progName = \"" << programName << "\";" << endl;
   printString ("description", preUsageMessage);
   cout << "authors = \"" << authors << "\";" << endl;
   cout << "version = \"" << versionText << "\";" << endl;

   SxList<CliGroup>::ConstIterator gIt = groups.begin ();
   cout << "groups {" << endl;
   for ( ; gIt != groups.end (); gIt++)
      gIt->printout ();
   cout << "}" << endl << endl;

   cout << "options {" << endl;
   SxList<CliArg>::ConstIterator opt = optionList.begin ();
   for ( ; opt != optionList.end (); opt++)
      opt->printout ();
   cout << "}" << endl;

}

void SxCLI::CliArg::printout () const
{
   const char *shift="   ";
   cout << shift << "option {" << endl;
   if (mark.getSize () > 0)  {
      cout << shift << "   marks = [\"" << mark(0) << '"';
      for (int i = 1; i < mark.getSize (); ++i)
         cout << ", \"" << mark(i) << '"';
      cout << "];" << endl;
   }
   if (emptyOption)
      cout << shift << shift << "flag;" << endl;
   if (optional)
      cout << shift << shift << "optional;" << endl;
   // looping is not yet covered
   // if (stickyOption)
   //   cout << shift << shift << "keepsValueForLoops;" << endl;
   printString ("      help", longDescription);
   if (shortDescription.getSize () > 0)  {
      cout << shift << shift;
      if (emptyOption)  {
         cout << "optionalValue;" << endl << shift << shift
              << "shortValue";
      } else {
         cout << "short";
      }
      cout << " = \"" << shortDescription.stripComments() << "\";" << endl;
      if (shortDescription.contains ("#"))  {
         cout << shift << shift
              << "tag = \"" << shortDescription.right("#") << "\";" << endl;
      }
   } else if (!emptyOption)  {
      cout << shift << shift << "short = \"???\";" << endl;
   }
   if (defaultValue.getSize () > 0)
      cout << shift << shift << "defaultInfo = \"" << defaultValue << "\";"
           << endl;
   SX_CHECK (group);
   cout << shift << shift << "groupId = " << group->id << ';' << endl;
   cout << shift << shift << "groupName = \"" << group->name << "\";" << endl;
   cout << shift << '}' << endl;
}

void SxCLI::CliGroup::printout () const
{
   const char *shift="   ";
   cout << shift << "group {" << endl;
   cout << shift << shift << "id = "   << id << ";" << endl;
   cout << shift << shift << "name = \"" << name << "\";" << endl;
   if (descr.getSize () > 0)
      printString ("      descr", descr);
   if (excludedGroups.getSize () > 0)  {
      cout << shift << shift << "excludes = [" << excludedGroups(0);
      for (int i = 1; i < excludedGroups.getSize (); ++i)
         cout << ", " << excludedGroups(i);
      cout << "]; // group ids" << endl;
   }
   if (requiredGroups.getSize () > 0)  {
      cout << shift << shift << "requires = [" << requiredGroups(0);
      for (int i = 1; i < requiredGroups.getSize (); ++i)
         cout << ", " << requiredGroups(i);
      cout << "]; // group ids" << endl;
   }
   cout << shift << '}' << endl;
}

void SxCLI::CliArg::printMark (bool capital) const
{
   if (mark.getSize () > 0)  {
      cout << (capital ? 'O' : 'o') << "ption ";
#ifdef WIN32
      if (mark(0)(0) == '-')  {
         const SxString &m = mark(0);
         SX_CHECK (m.getSize () > 1);
         cout << '/' << m.subString ((m(1) == '-') ? 2 : 1);
      } else // print mark as is:
#endif
      {
         cout << mark(0);
      }
   } else {
      cout << "<" << shortDescription.stripComments() << ">";
   }
   if (mark.getSize () > 1) cout << " or alternatives";
}

void SxCLI::CliArg::tooSmallError (const SxString& val, const SxString& min)
{
   if (parent->mode == Completion) return;
   printMark ();
   cout << ": Value is " << val << ", but must at least " << min << "!" << endl;
   parent->setError ();
}

void SxCLI::CliArg::tooLargeError (const SxString& val, const SxString& max)
{
   if (parent->mode == Completion) return;
   printMark ();
   cout << ": Value is " << val << ", but must not be larger than ";
   cout << max << "!" << endl;
   parent->setError ();
}

void SxCLI::CliArg::setDefault(const SxString& def,
                               const SxString& min,
                               const SxString& max)
{
   if (def.getSize () + min.getSize () + max.getSize () == 0) return;
   if (def.getSize () > 0)
      defaultValue = SxString("default: ") + def;
   if (min.getSize () > 0)
      defaultValue += SxString ("; >= ") + min; 
   if (max.getSize () > 0)
      defaultValue += SxString ("; <= ") + max; 
   if (defaultValue(0) == ';') // remove leading "; "
      defaultValue = defaultValue.subString (2);   
}

int SxCLI::CliArg::toInt (Given<int> def, Given<int> min, Given<int> max)
{
   // --- check def, min and max are consistent
   SX_CHECK (!(min.given && max.given && (min.value > max.value)),
               min.value, max.value);
   SX_CHECK (!(min.given && def.given && (min.value > def.value)),
               min.value, def.value);
   SX_CHECK (!(def.given && max.given && (def.value > max.value)),
               def.value, max.value);

   if (!hasValue ())  {
      // --- set option behaviour
      optional      = optional && def.given; // can't be optional without default
      emptyOption   = false;
   }

   // --- set up default and range description
   setDefault(def.given ? SxString(def.value) : "",
              min.given ? SxString(min.value) : "",
              max.given ? SxString(max.value) : "");

   int i = (def.given ? def.value : 0);
   // stop if error occured
   //if (parent->error) return i;

   // --- get argument and integer from that (or default)
   SxString arg = hasValue () ? getValue () : getArgument ();
   if (arg != "")  {
      try {
         i = arg.toInt ();
      } catch (SxException e) {
         cout << "Can't convert '" << arg << "' to integer number.\n"; 
         parent->setError ();
         return i;
      }
      // --- check range
      if (min.given && (min.value > i))
         tooSmallError (SxString(i), SxString(min.value));
      if (max.given && (max.value < i))
         tooLargeError (SxString(i), SxString(max.value));
   }
   return i;
}
   
long int SxCLI::CliArg::toLong (Given<long int> def,
                                Given<long int> min,
                                Given<long int> max)
{
   // --- check def, min and max are consistent
   SX_CHECK (!(min.given && max.given && (min.value > max.value)),
               min.value, max.value);
   SX_CHECK (!(min.given && def.given && (min.value > def.value)),
               min.value, def.value);
   SX_CHECK (!(def.given && max.given && (def.value > max.value)),
               def.value, max.value);

   // --- set option behaviour
   if (!hasValue ()) {
      optional      = optional && def.given; // can't be optional without default
      emptyOption   = false;
   }

   // --- set up default and range description
   setDefault(def.given ? SxString(def.value) : "",
              min.given ? SxString(min.value) : "",
              max.given ? SxString(max.value) : "");

   long int i = (def.given ? def.value : 0);
   // stop if error occured
   //if (parent->error) return i;

   // --- get argument and integer from that (or default)
   SxString arg = hasValue () ? getValue () : getArgument ();
   if (arg != "")  {
      try  {
         i = arg.toInt64 ();
      } catch (SxException e)  {
         cout << "Can't convert '" << arg << "' to integer number.\n";
         parent->setError ();
         return i;
      }
      // --- check range
      if (min.given && (min.value > i))
         tooSmallError (SxString(i), SxString(min.value));
      if (max.given && (max.value < i))
         tooLargeError (SxString(i), SxString(max.value));
   }
   return i;
}

int64_t SxCLI::CliArg::toInt64 (Given<int64_t> def,
                                Given<int64_t> min,
                                Given<int64_t> max)
{
   // --- check def, min and max are consistent
   SX_CHECK (!(min.given && max.given && (min.value > max.value)),
               min.value, max.value);
   SX_CHECK (!(min.given && def.given && (min.value > def.value)),
               min.value, def.value);
   SX_CHECK (!(def.given && max.given && (def.value > max.value)),
               def.value, max.value);

   // --- set option behaviour
   if (!hasValue ()) {
      optional      = optional && def.given; // can't be optional without default
      emptyOption   = false;
   }

   // --- set up default and range description
   setDefault(def.given ? SxString(def.value) : "",
              min.given ? SxString(min.value) : "",
              max.given ? SxString(max.value) : "");

   int64_t i = (def.given ? def.value : 0);
   // stop if error occured
   //if (parent->error) return i;

   // --- get argument and integer from that (or default)
   SxString arg = hasValue () ? getValue () : getArgument ();
   if (arg != "")  {
      try  {
         i = arg.toInt64 ();
      } catch (SxException e)  {
         cout << "Can't convert '" << arg << "' to 64-bit integer number.\n";
         parent->setError ();
         return i;
      }
      // --- check range
      if (min.given && (min.value > i))
         tooSmallError (SxString(i), SxString(min.value));
      if (max.given && (max.value < i))
         tooLargeError (SxString(i), SxString(max.value));
   }
   return i;
}

   
float SxCLI::CliArg::toFloat (Given<float> def,
                              Given<float> min,
                              Given<float> max)
{
   // --- check def, min and max are consistent
   SX_CHECK (!(min.given && max.given && (min.value > max.value)),
               min.value, max.value);
   SX_CHECK (!(min.given && def.given && (min.value > def.value)),
               min.value, def.value);
   SX_CHECK (!(def.given && max.given && (def.value > max.value)),
               def.value, max.value);

   // --- set option behaviour
   if (!hasValue ())  {
      optional      = optional && def.given; // can't be optional without default
      emptyOption   = false;
   }

   // --- set up default and range description
   setDefault(def.given ? SxString(def.value) : "",
              min.given ? SxString(min.value) : "",
              max.given ? SxString(max.value) : "");

   float f = (def.given ? def.value : 0.f);
   // stop if error occured
   //if (parent->error) return f;

   // --- get argument and float from that (or default)
   SxString arg = hasValue () ? getValue () : getArgument ();
   if (arg != "")  {
      try  {
         f = arg.toFloat ();
      } catch (SxException e)  {
         cout << "Can't convert '" << arg << "' to number.\n";
         parent->setError ();
         return f;
      }
      // --- check range
      if (min.given && (min.value > f))
         tooSmallError (SxString(f), SxString(min.value));
      if (max.given && (max.value < f))
         tooLargeError (SxString(f), SxString(max.value));
   }
   return f;
}

double SxCLI::CliArg::toDouble (Given<double> def, 
                                Given<double> min, 
                                Given<double> max)
{
   // --- check def, min and max are consistent
   SX_CHECK (!(min.given && max.given && (min.value > max.value)),
               min.value, max.value);
   SX_CHECK (!(min.given && def.given && (min.value > def.value)),
               min.value, def.value);
   SX_CHECK (!(def.given && max.given && (def.value > max.value)),
               def.value, max.value);

   // --- set option behaviour
   if (!hasValue ()) {
      optional      = optional && def.given; // can't be optional without default
      emptyOption   = false;
   }

   // --- set up default and range description
   setDefault(def.given ? SxString(def.value) : "",
              min.given ? SxString(min.value) : "",
              max.given ? SxString(max.value) : "");

   double f = (def.given ? def.value : 0.);
   // stop if error occured
   //if (parent->error) return f;

   // --- get argument and double from that (or default)
   SxString arg = hasValue () ? getValue () : getArgument ();
   if (arg != "")  {
      try  {
         f = arg.toDouble ();
      } catch (SxException e)  {
         cout << "Can't convert '" << arg << "' to number.\n";
         parent->setError ();
         return f;
      }
      // --- check range
      if (min.given && (min.value > f))
         tooSmallError (SxString(f), SxString(min.value));
      if (max.given && (max.value < f))
         tooLargeError (SxString(f), SxString(max.value));
   }
   return f;
}

SxString SxCLI::CliArg::toString (Given<SxString> def)
{
   // --- set option behaviour
   if (!hasValue ())  {
      optional      = optional && def.given; // can't be optional without default
      emptyOption   = false;
   }

   // --- set up default
   if (def.given && (def.value.getSize () > 0))
      defaultValue = SxString ("default: ") + SxString (def.value);

   // stop if error occured
   //if (parent->error) return "";

   // --- get argument (or default)
   if (hasValue ()) return getValue ();
   SxString arg = getArgument ();
   if ((arg.getSize () == 0) && def.given)
      arg = def.value;
   return arg;
}

int SxCLI::CliArg::toChoice (Given<int> def)
{
   if (!hasValue ()) emptyOption = true;
   SxString arg = getArgument ();
   if (arg.getSize () == 0)  {
      if (def.given) return def.value;
      return (optional ? -1 : 0);
   }
   int pos = int(mark.findPos (arg));
#ifdef WIN32
   if (pos < 0 && arg(0) == '/')  {
      arg(0) = '-';
      pos = (int)mark.findPos (arg);
      if (pos < 0)
         pos = (int)mark.findPos ("-" + arg);
   }
#endif
   SX_CHECK (pos >= 0);
   return pos;
}

SxList<double> SxCLI::CliArg::toList3 (const SxString &separators,
                                       const SxString &braces)
{
   // --- set option behaviour
   if (!hasValue ())  {
      optional = false;
      emptyOption = false;
   }

   SxList<double> res;
   res << 0. << 0. << 0.;

   // stop if error occured
   //if (parent->error) return res;

   // --- get argument
   SxString arg = hasValue () ? getValue () : getArgument ();
   if (arg.getSize () == 0) return res;

   SxString orig = arg; // for error messages
   
   // --- take away braces
   for (int i = 0; i < braces.getSize (); i+=2)
      if (arg.contains (braces(i)))  {
         if (!arg.contains (braces(i+1)))  {
            cout << "Parsing error in vector '" << orig << "':" << endl;
            cout << "Opening is '" << braces(i) << "', but closing '";
            cout << braces(i+1) << "' not found." << endl;
            parent->setError ();
            return res;
         }
         arg = arg.right(braces(i)).left(braces(i+1));
         break;
      }

   // substitute separator
   char sep;
   for (int i = 0; i < separators.getSize (); i++)  {
      sep = separators(i);
      int nsep = arg.contains (sep);
      if (nsep == 0) continue;
      if (nsep != 2)  {
         cout << "Error while parsing vector '" << orig;
         cout << "'.\nFound separator '" << sep << "' ";
         cout << nsep << " time";
         if (nsep > 1) cout << 's';
         cout << ", instead of twice." << endl;
         parent->setError ();
         return res;
      }
      arg = arg.substitute (sep," ", 2);
      break;
   }

   // perform conversion
   int nconv = sscanf (arg.ascii (), "%lf %lf %lf", &res(0), &res(1), &res(2));
   if (nconv != 3)  {
      cout << "Error while parsing vector '" << orig << "'." << endl;
      cout << "After removing braces and separators we have '" << arg;
      cout << "'." << endl << "Only " << nconv;
      cout << " element(s) could be properly converted." << endl;
      parent->setError ();
   }
   return res;
}

SxList<int> SxCLI::CliArg::toIntList3 (const SxString &separators,
                                       const SxString &braces)
{
   SxList<double> resDouble = toList3 (separators, braces);
   SxList<int> res;
   for (int i = 0; i < 3; ++i)
      res << int(lround(resDouble(i)));
   return res;
}

SxList<SxString> SxCLI::CliArg::toList ()
{
   SX_CHECK (group);
   emptyOption = false;

   SxList<SxString> list;
   if (parent->error) return list;

   if (mark.getSize () == 0)  {
      list = parent->arguments(parent->iLoop);
      parent->arguments(parent->iLoop).resize(0);
      return list;
   }

   while (!parent->error && exists (true))
      list << getArgument ();
   if (list.getSize () == 0 && !optional && group->status == Required)  {
      printMark ();
      cout << " is required at least once!" << endl;
      parent->setError ();
   }
   
   return (parent->error ? SxList<SxString> () : list);
}

SxList<double> SxCLI::CliArg::toDoubleList (const SxString &separators)
{
   SX_CHECK (group);
   if (!hasValue ()) emptyOption = false;

   SxList<double> list;
   //if (parent->error) return list;
   
   // get argument
   SxString listString = hasValue () ? getValue () : getArgument ();
   if (parent->error || listString.getSize () == 0) return list;

   // find used separator
   char sep = 0;
   for (int i = 0; i < separators.getSize () ; ++i)
      if (listString.contains(sep = separators(i))) break; 
   
   // split argument on sep
   SxList<SxString> tokens = listString.tokenize (sep);

   // to Double
   SxList<SxString>::Iterator it;
   for (it = tokens.begin (); it != tokens.end (); ++it)  {
      try  {
         list << (*it).toDouble ();
      } catch (SxException e)  {
         cout << "Can't convert '" << (*it) << "' to number.\n";
         // append is not performed
      }
   }

   return list;

}

SxList<int> SxCLI::CliArg::toIntList (const SxString &separators)
{
   SX_CHECK (group);
   if (!hasValue ()) emptyOption = false;

   SxList<int> list;
   //if (parent->error) return list;
   
   // get argument
   SxString listString = hasValue () ? getValue () : getArgument ();
   if (parent->error || listString.getSize () == 0) return list;

   // find used separator
   char sep = 0;
   for (int i = 0; i < separators.getSize () ; ++i)
      if (listString.contains(sep = separators(i))) break; 
   
   // split argument on sep
   SxList<SxString> tokens = listString.tokenize (sep);

   // to integer
   SxList<SxString>::Iterator it;
   for (it = tokens.begin (); it != tokens.end (); ++it)  {
      try {
         list << (*it).toInt ();
      } catch (SxException e) {
         cout << "Can't convert '" << (*it) << "' to integer.\n";
         // append is not performed
      }
   }

   return list;

}

SxList<int> SxCLI::CliArg::toIdxList ()
{
   if (!hasValue ()) emptyOption = false;
   SxList<int> idxList;
   //if (parent->error) return idxList;
   SxString listString = hasValue () ? getValue () : getArgument ();
   if (listString.getSize () == 0) return idxList;
   // cut listString in ","-separated items
   SxList<SxString> items = listString.tokenize (',');
   // alternative separator ":" compatible to toDoubleList
   if (items.getSize () == 1 && listString.contains (":"))
      items = listString.tokenize (':');
   SxList<SxString>::Iterator item;
   try  {
      for (item = items.begin (); item != items.end (); ++item)
      {
         // --- parse element
         ssize_t minusPos = (*item).find("-");
         if ((minusPos == 0) || (minusPos == (*item).getSize () - 1) )  {
            printMark ();
            cout << ": error in list " << listString << endl;
            cout << "In " << (*item) << ": no '-' in ";
            cout << (minusPos == 0 ? "first" : "last");
            cout << " position allowed. Use explicit limits." << endl;
            parent->setError ();
            return SxList<int> ();
         }
         if ( minusPos > 0 )  {
            // --- "<n>-<m>" gives numbers from n to m
            int from = (*item).left("-").toInt();
            int to = (*item).right("-").toInt();

            if ( from > to)  {
               printMark ();
               cout << ": error in list " << listString << endl;
               cout << "In " << (*item) << ": " << from << " > " << to << endl;
               parent->setError ();
               return SxList<int> ();
            }
            // add numbers to list (start from 0 instead of from 1)
            for (int i = from; i <= to; i++)  idxList << (i - 1);
         } else {
            // --- just a number (start from 0 instead of from 1)
            idxList << ((*item).toInt()-1);
         }
      }
   } catch (SxException e)  {
      e.print ();
      parent->setError ();
      idxList.resize (0);
   }
   return idxList;
}

SxList<int64_t> SxCLI::CliArg::toIdxList64 ()
{
   if (!hasValue ()) emptyOption = false;
   SxList<int64_t> idxList;
   //if (parent->error) return idxList;
   SxString listString = hasValue () ? getValue () : getArgument ();
   if (listString.getSize () == 0) return idxList;
   // cut listString in ","-separated items
   SxList<SxString> items = listString.tokenize (',');
   // alternative separator ":" compatible to toDoubleList
   if (items.getSize () == 1 && listString.contains (":"))
      items = listString.tokenize (':');
   SxList<SxString>::Iterator item;
   try  {
      for (item = items.begin (); item != items.end (); ++item)
      {
         // --- parse element
         ssize_t minusPos = (*item).find("-");
         if ((minusPos == 0) || (minusPos == (*item).getSize () - 1) )  {
            printMark ();
            cout << ": error in list " << listString << endl;
            cout << "In " << (*item) << ": no '-' in ";
            cout << (minusPos == 0 ? "first" : "last");
            cout << " position allowed. Use explicit limits." << endl;
            parent->setError ();
            return SxList<int64_t> ();
         }
         if ( minusPos > 0 )  {
            // --- "<n>-<m>" gives numbers from n to m
            int64_t from = (*item).left ("-").toInt64 ();
            int64_t to   = (*item).right("-").toInt64 ();

            if ( from > to )  {
               printMark ();
               cout << ": error in list " << listString << endl;
               cout << "In " << (*item) << ": " << from << " > " << to << endl;
               parent->setError ();
               return SxList<int64_t> ();
            }
            // add numbers to list (start from 0 instead of from 1)
            for (int64_t i = from; i <= to; i++)  idxList << (i - 1);
         } else  {
            // --- just a number (start from 0 instead of from 1)
            idxList << ((*item).toInt64 ()-1);
         }
      }
   } catch (SxException e)  {
      e.print ();
      parent->setError ();
      idxList.resize (0);
   }
   return idxList;
}

bool SxCLI::CliArg::exists (int iLoop_, bool keep)
{
   SX_CHECK (parent);
   if (iLoop_ == -1) iLoop_ = parent->iLoop;
   // return true if iLoop exceeds number of loops to prevent sticky option
   // propagation
   if (iLoop_ >= parent->arguments.getSize ()) return true;
   
   // no existence checks for unnamed options! 
   // look at arguments(getILoop ()).getSize () if you need to find an argument
   SX_CHECK (mark.getSize () != 0);
   SX_CHECK (group);
   SxList<SxString> &args_ = parent->arguments(parent->iLoop);
   
   if (!keep)  {
      emptyOption = true;
   }

   int length;
   for (int iArg = 0; iArg < args_.getSize (); iArg++)  {
      // does mark fit ?
      length = fits (args_(iArg));
      if (length == 0) continue;
      if (emptyOption && (length < args_(iArg).getSize ()))  {
         if (args_(iArg)(length) == '=')  {
            optValue = args_(iArg).subString(length);
         } else {
            continue;
         }
      }
      // remove arg unless it should be kept
      if (!keep)  {
         // propagate sticky flag to next loop
         if (stickyOption && parent->iLoop + 1 < parent->arguments.getSize ())
           parent->arguments(parent->iLoop + 1) << args_(iArg);
         if (group->status == MissingRequired)  {
            // --- issue error if a required option of this group is missing
            group->initOption->printMark ();
            cout << " is required if '" << args_(iArg);
            cout << "' is given." << endl;
            parent->setError ();
         }
         args_.remove(iArg);
         group->setRequired (this);
         used = true;
      }
      return true;
   }
   if (!keep && !optional)  {
      if (group->status == Required)  {
         // --- issue error if flag is missing, although it is required
         // for the current group
         printMark ();
         cout << " is required if ";
         group->initOption->printMark ();
         cout << " is given." << endl;
         parent->setError ();
      } else {
         group->status = MissingRequired;
         group->initOption = this;
      }
   }
   return false;
}

// ------------------------- SxCLI::CliGroup --------------------------------
SxCLI::CliGroup::CliGroup (int newId)
{
   id = newId;
   status = Optional;
   initOption = NULL;
}

const int SxCLI::generalGroup = 0;

void SxCLI::CliGroup::printName ()
{
   if (name.getSize () > 0)
      cout << name;
   else
      cout << "Group " << id;
}

// ------------------------- SxCLI ------------------------------------------

void SxCLI::init (const SxList<SxString> &args_)
{
   arguments.resize(1);  // single loop
   arguments(0) = args = args_;

   // no errors so far
   error = false;
   iLoop = 0;

   // --- get program name
   const ssize_t BUFLEN = 10240;
   char cPath[BUFLEN];
   sxGetExecPath (cPath, BUFLEN);
   execPath = SxString (cPath);
   
   programName = execPath.substitute("\\", "/");

#  ifdef WIN32
      SxList<SxString> tokens = execPath.tokenize ('\\');
      tokens.removeLast ();
      execDir = SxString::join (tokens, "\\");
#  else
      SxList<SxString> tokens = execPath.tokenize ('/');
      tokens.removeLast ();
      execDir = "/" + SxString::join (tokens, "/");
#  endif

   if (programName.contains ("/")) 
      programName = programName.tokenize ('/').last ();
   // cut off final ".x"
   SX_CHECK (programName.getSize () > 2, programName.getSize ());
   if (programName.subString (programName.getSize () - 2) == ".x")  {
      programName.resize (programName.getSize () - 2, true);
      programName = programName.toLower ();
   }
   // cut off leading lt-
   SX_CHECK (programName.getSize () >= 3, programName.getSize ());
   if (programName.subString (0,2) == "lt-")  {
      SX_CHECK (programName.getSize () > 3, programName.getSize ());
      programName = programName.subString(3);
   }

   autoOptionGroup.name = "Standard options";

   init ();
   
}

void SxCLI::init ()
{
   autoOptionNames = SxString("-h|--help|--usage|help|HELP|usage"
                              "|--opts|--about|-v|--version|--log|--quiet"
                              "|--debug|--memcheck|--no-exceptions"
                              "|--completion|--sxprintcli"
                     ).tokenize ('|');
   // hidden option: completion mode
   currentGroup = &autoOptionGroup;
   mode = Normal;
   if (iLoop == 0)  {
      if (option ("--completion", "complete line").toBool ())
         mode = Completion;
      else if (option ("--sxprintcli", "file", "print CLI in sx-style").toBool ())
      {
         if (last ().hasValue ())
            Log::get ().sxTee = SxPtr<SxRedirect>::create 
                               (std::cout, last ().getValue ().ascii (), false);
         mode = PrintSxFormat;
      }
   }

   // --- restart setting up the options
   optionList.resize (0);
   helpOptionPtr = NULL;


   // --- initialize general group
   groups.resize(0);
   groups << CliGroup (generalGroup);
   currentGroup = &groups.last ();
   currentGroup->status = Required; // general is always required
   currentGroup->name   = "General options";
}


SxCLI::~SxCLI ()
{
   // empty
}

void SxCLI::setLoopSeparator (const SxString &mark)
{
   // setLoopSeparator may be called only once
   SX_CHECK (arguments.getSize () == 1);

   SxList<SxString>::Iterator allArgIt = arguments(0).begin (),
                              endArg = arguments(0).end ();
   CliArg separator;
   separator.mark = mark.tokenize ('|');
   
   // we keep arguments(0) as original until we finished loop separation
   arguments.append (SxList<SxString> ());

   for (; allArgIt != endArg; ++allArgIt)  {
      if (separator.fits (*allArgIt) && arguments.last ().getSize () > 0)  {
         // new loop: append next loop list
         arguments.append (SxList<SxString> ());
      }
      // append argument to current loop
      arguments.last ().append (*allArgIt);
   }
   // remove original argument list
   arguments.remove(0);
} 

bool SxCLI::looping ()
{
   if (optionList.getSize () == 0) return true; // first round
   init ();
   // increase loop counter and check it is still valid
   return (++iLoop < arguments.getSize ());
}

void SxCLI::printUsage (int errorCode, bool printLong)
{
   SX_CHECK (optionList.getSize () != 0);

   // --- make sure that internal options are listed last
   ssize_t insertPoint = -1;
   while (optionList.last ().group != &autoOptionGroup)  {
      if (insertPoint == -1)  {
         ssize_t i = optionList.getSize () - 1;
         for (; i; --i)
            if (optionList(i).group == &autoOptionGroup) break;
         for (; i; --i)
            if (optionList(i).group != &autoOptionGroup) break;
         insertPoint = i + 1;
      }
      // move last element into list before the internal options
      optionList.insert (insertPoint, optionList.last ());
      optionList.removeLast ();
   } 

   if ((printLong || !helpOptionPtr) && preUsageMessage.getSize () > 0)  {
      cout << endl << preUsageMessage.wrap () << endl;
   }
   cout << endl << "Usage:" << endl;
   SxString line;
   int iMark;
   ssize_t maxMarkSize = 0; // longest mark
   ssize_t markSize; // length of current mark
   int currentId = generalGroup;
   for (int i = 0; i < optionList.getSize (); i++)  {
      if (optionList(i).group->id != currentId)  {
         // --- new group
         CliGroup *group = optionList(i).group;
         if (currentId < 1)  {
            // starts with '{' unless general group
            line += (group->id > 0) ? " {" : " ";
         } else  {
            if (group->excludedGroups.contains (currentId))  {
               // consecutive excluded groups are "{group1} | {group2}"
               line += "} | {";
            } else  {
               // group ends with ']'
               line += "} ";
               if (group->id > 0) line += "{"; // start group
            }
         }
         currentId = group->id;
      } else  {
         line += " "; // put space between options
      }
      if (optionList(i).optional) line += "[";
      // --- loop over alternatives
      for (iMark = 0; iMark < optionList(i).mark.getSize (); iMark++)  {
         if (iMark != 0) line += "|";
         SxString m = optionList(i).mark(iMark);
#ifdef WIN32
         // substitute - and -- with /
         if (m(0) == '-')  {
            if (m(1) == '-') m = m.subString(1);
            m(0) = '/';
         }
#endif
         // --- find maxMarkSize for longDescription alignment below
         markSize = m.getSize ();
         if (markSize > maxMarkSize) maxMarkSize = markSize;
         // ---  add mark alternative to line
         line += m;
      }
      if (!optionList(i).emptyOption)  {
         if (optionList(i).mark.getSize () > 0) line += " ";
         if (optionList(i).shortDescription.getSize () > 0)
            line += "<" + optionList(i).shortDescription.stripComments() + ">";
         else
            line += "<??\?>";
      } else if (optionList(i).shortDescription.getSize () > 0)  {
         line += "[=" + optionList(i).shortDescription.stripComments() + "]";
      }

      if (optionList(i).optional)
         line += "]";
   }
   if (currentId > 0) 
      line += "}"; // close last group
   cout << line.subString(1).wrap(programName, 1) << endl; 
   if (!printLong && helpOptionPtr)  {
      helpOptionPtr->printMark ();
      cout << " lists the detailed help." << endl;
      cout.flush ();
      if (errorCode != -1)  {
   	runAtExitFuncs ();
         _Exit (errorCode);
      } else {
         return;
      }
   }

   currentId = generalGroup;
   for (int i = 0; i < optionList.getSize (); i++)  {
      SX_CHECK (optionList(i).group);
      if (optionList(i).group->id != currentId)  {
         CliGroup *group = optionList(i).group;
         if (group->name == "...")  continue;
         currentId = group->id;
         cout << endl;
         group->printName ();
         if (group->excludedGroups.getSize () > 0)  {
            cout << " (excludes ";
            for (SxList<int>::Iterator it = group->excludedGroups.begin ();
                 it != group->excludedGroups.end ();
                 it++)
            {
               if (it != group->excludedGroups.begin ()) cout << ", ";
               groups(*it).printName ();
            }
            cout << ")";
         }
         if (group->requiredGroups.getSize () > 0)  {
            cout << " (requires ";
            for (SxList<int>::Iterator it = group->requiredGroups.begin ();
                 it != group->requiredGroups.end ();
                 it++)
            {
               if (it != group->requiredGroups.begin ()) cout << ", ";
               groups(*it).printName ();
            }
            cout << ")";
         }
         cout << endl;
         if (group->descr.getSize () > 0)
            cout << group->descr.wrap () << endl;
      }
      if (optionList(i).mark.getSize () > 0)  {
         // --- loop over all but last alternative ... 
         for (iMark = 0; iMark < optionList(i).mark.getSize () - 1; iMark++)  {
            cout << "  ";
#ifdef WIN32
            const SxString &m = optionList(i).mark(iMark);
            if (m(0) == '-')  {
               cout << '/' 
                    << m.subString ((m(1) == '-') ? 2 : 1)
                    << endl;
            } else // just print mark
#endif
            cout << optionList(i).mark(iMark) << endl;
         }
         // ... because the last alternative is needed to align default
         line = optionList(i).mark.last ();
#ifdef WIN32
         if (line(0) == '-')  {
            if (line(1) == '-') line = line.subString(1);
            line(0) = '/';
         }
#endif
      } else {
         // no mark: use short description
         line = optionList(i).shortDescription.stripComments ();
      }
      cout << "  ";
      if (optionList(i).defaultValue.getSize () > 0)  {
         line.resize(maxMarkSize, true);
         cout << line << "    (" << optionList(i).defaultValue << ")" << endl;
      } else {
         cout << line << endl;
      }
      cout << optionList(i).longDescription.wrap("", maxMarkSize + 6);
      cout << endl;
   }
   cout.flush ();
   if (errorCode != -1) { runAtExitFuncs (); _Exit (errorCode); }
}

void SxCLI::printCompOptions () const
{
   cout << endl;
   cout << SX_SEPARATOR;
   if (versionText.getSize() > 0)
      cout << "| Version:             " << versionText << endl;
   if (release.getSize() > 0)
      cout << "| Release number:      " << release << endl;
   cout << "| Compiled by:         " << WHOAMI  << endl;
   cout << "| Date:                " << dateStr << endl;
   long int timeNumber;
   if (sscanf(dateStr.ascii (), "%ld", &timeNumber) == 1)  {
      // translate number to human-readable form
      cout << "|                      " << SxTime::ctime(timeNumber) << endl;
   }
//   cout << "| C++ compiler:        " << CXX << " [" << CXXVERSION << "]\n";
   cout << "| C++ compiler:        " << CXX << " [" << CXXVERSION << "]" << endl;
   cout << SxString(CXXFLAGS).wrap ("| Compiler options:    ")
           .substitute("\n ","\n|") << endl;
#  ifdef FC
      cout << "| Fortran compiler:    " << FC << endl;
      cout << "| Fortran options:     " << FFLAGS << endl;
#  endif
   cout << "| Linker options:      " << LDFLAGS << endl;
   cout << SX_SEPARATOR;
   cout << endl;
   cout.flush ();
   runAtExitFuncs ();
   _Exit (0);
}

void SxCLI::addAutoOptions ()
{
   if (autoOptionNames.getSize () == 0) return;
   // --- switch to automatic option mode
   CliGroup *backupCurrent = currentGroup;
   currentGroup = &autoOptionGroup;

   // --- handle options
   // note: all automatic options with more than one mark alternative
   //       must allow for overriding some of these marks. For this,
   //       add all option marks to autoOptionNames. From the list of
   //       remaining alternatives, only the first three will be used.

   // "--log-components"
   if (option ("--log-components", "options",
               "Define which components should enable logging").toBool ())
   {
      SxString arg = last().getValue ();
      try {
         last().getValue().tokenize(',').foreach ([&](auto it) {
            SxLog::enable (it->ascii());
         });
      } catch (...) { }
   }


   // "--log"
   SxString logName = programName + ".log"; 
   bool createLog = option ("--log", "filename",
      "Create a log file. The creation of log files can also be enabled "
      "by setting the environment variable SX_LOG_STDOUT").toBool();
   if (last ().hasValue ()) logName = last ().getValue ();
   if (::getenv ("SX_LOG_STDOUT"))  createLog = true;
   if (logName.getSize () == 0) createLog = false;

   // "--quiet"
   bool quiet = option ("--quiet", "Do not print to stdout").toBool ();

#  ifndef NDEBUG
   // "--debug"
#  ifndef SX_MOBILE
      if (option ("--debug", "Invoke debugger, e.g., --debug core").toBool ())  {
         SxString cmd;
#     ifdef WIN32
            cmd  = SxString (GDB) + " ";
#     else
            cmd  = SxString(libtool) + " --mode=execute ";
            cmd += SxString(GDB) + " ";
            cmd += execPath;
#     endif
         while (arguments(iLoop).getSize () > 0)  {
            cmd += " " + arguments(iLoop).first ();
            arguments(iLoop).removeFirst ();
         }
         cout << "Executing " << cmd.ascii() << endl; cout.flush ();
         int r = system (cmd.ascii());
         runAtExitFuncs ();
         _Exit (r);
      }

      // "--memcheck"
      if (option ("--memcheck", "options", "Invoke memory tracer").toBool ())  {
         if (SxString(MEMTRACER) == "")  {
            cout << "no memory tracer was configured" << endl;
            SX_EXIT;
         }
         SxString cmd;
#     ifdef WIN32
         cmd = SxString("drmemory.exe -batch") + " ";
#     else
         cmd  = SxString(libtool) + " --mode=execute ";
         cmd += SxString(MEMTRACER) + " ";
#     endif /* WIN32 */
         if (last ().hasValue ()) cmd += last ().getValue () + " ";
         cmd += SxString(execPath);
         for (int i = 1; i < args.getSize(); ++i)  {
            SxString arg = args(i);
            if (! last ().fits (arg))
               cmd += " " + arg;
         }
         cout << "Executing " << cmd.ascii() << endl; cout.flush ();
         int r = system (cmd.ascii());
         runAtExitFuncs ();
         _Exit(r);
      }
#  endif /* SX_MOBILE */

   // "---no-exceptions"
   bool segFaultExcpt = option ("--no-exceptions", 
                                "Throw segfaults rather than handle exceptions"
                        ).toBool ();
   if (segFaultExcpt)  SxException::causeSegFault ();

#  endif /* NDEBUG */

   // "--opts"
   if (option ("-v|--version|--opts", "Print compiler options").toBool())
      printCompOptions ();

   // "--help"
   helpOptionPtr = &option ("-h|--help|--usage|help|HELP|usage",
                            "this help message");
   if (helpOptionPtr->toBool () && mode == Normal)  {
      mode = PrintHelp;
      setError (); // stop other options from emitting errors
   }

   if (mode != Normal) {
      createLog = false;
      quiet = false;
   }

   // "--about"
   if (option ("--about", "print information about package").toBool ()
       && mode == Normal)
   {
      cout << SX_SEPARATOR;
      if (copyright.getSize() > 0)  {
         cout << SxString::sprintf ("| %-20s %s",
                    (programName + " " + versionText).ascii(), 
                    authors.ascii()) << endl;
         cout << "|" << endl;
         cout << copyright;
      } else {
         if (release.getSize () > 0)
            cout << "| " << package << " release " << release << endl;
         cout << "| based on S/PHI/nX package by S. Boeck, J. Neugebauer et al."
              << endl;
         cout << "| " << programName << " ";
         if (authors.getSize () > 0)
            cout << "by " << authors << endl << "| Version ";
         cout << versionText << endl;
      }
      cout << SX_SEPARATOR;
      cout.flush ();
      runAtExitFuncs ();
      _Exit (0);
   }
   if (createLog && ! Log::get ().enabled)  {
      cout << "Disabling logging: " << Log::get ().reason << endl;
      createLog = false;
   }

   // --- redirection
   if (createLog)
      Log::get ().sxTee = SxPtr<SxRedirect>::create (
                             std::cout, logName.ascii(), !quiet
                          );
   if (!createLog && quiet)
      Log::get ().sxTee = SxPtr<SxRedirect>::create (
                             std::cout, SxRedirect::DevNull
                          );

   // --- switch back to user options
   currentGroup = backupCurrent;
   autoOptionNames.resize (0);
}

void SxCLI::printCompletion ()
{
   if (arguments(iLoop).getSize () > 0)  {
      for (int ia = 0; ia < arguments(iLoop).getSize () - 1; ++ia)  {
         cerr << "Unknown option " << arguments(iLoop)(ia) << endl;
      }
      // try completion
      SxString item = arguments(iLoop).last ();
      for (int io = 0; io < optionList.getSize (); ++io)  {
         const CliArg &arg = optionList(io);
         if (arg.used  || arg.group->status == Forbidden) continue;
         for (int im = 0; im < arg.mark.getSize (); ++im)  {
            if (arg.mark(im).getSize () > item.getSize ())  {
               if (arg.mark(im)(0, item.getSize () - 1) == item)  {
                  cout << arg.mark(im) << endl;
               }
            }
         }
      }
      runAtExitFuncs ();
      exit (0);
   }
   SxUniqueList<int> incompleteGroups;
   for (int io = 0; io < optionList.getSize (); ++io)  {
      const CliArg &arg = optionList(io);
      /* Check if a required option is missing from a group that is in use
         (status=Required). If so, we will only print the required options
         from these groups below, in order to provide guidance through complex
         CLIs.
         However, exclude the general group from this behavior. If no special
         option groups are being used, it feels unnatural to see all options only
         after the generally required ones have been provided.
      */
      if (!arg.used && !arg.optional && arg.group->status == Required
          && arg.group->id != generalGroup)
         incompleteGroups << arg.group->id;
   }
   bool showAll = incompleteGroups.getSize () == 0;

   for (int io = 0; io < optionList.getSize (); ++io)  {
      const CliArg &arg = optionList(io);
      bool show = showAll 
                || (incompleteGroups.contains(arg.group->id) && !arg.optional);
      if (!arg.used  && arg.group->status != Forbidden && show) {
         for (int im = 0; im < arg.mark.getSize (); ++im)
            cout << ' ' << arg.mark(im) << endl;
      }
   }
   runAtExitFuncs ();
   exit (0);
}

void SxCLI::finalize (enum ExtraArgs extra)
{
   setGroup (generalGroup);
   addAutoOptions ();

   if (mode == Completion) printCompletion ();
   if (mode == PrintHelp) printUsage (0);
   if (mode == PrintSxFormat) { printout (); runAtExitFuncs (); exit (0); }

   bool ignoreUnknown = (extra == IgnoreExtraArgs);
   bool unknownArgs = (arguments(iLoop).getSize () > 0);
   if (!error && (!unknownArgs || ignoreUnknown)) return; // everything is fine

   // Tell about unknown options only if no other error occured
   if (unknownArgs && !error)  {
      cout << "Unknown option";
      if (arguments(iLoop).getSize () > 1) cout << "s";
      cout << ":" << endl; 
      for (int iArg = 0; iArg < arguments(iLoop).getSize (); iArg++)
         cout << arguments(iLoop)(iArg) << " ";
      cout << endl;
   }
   printUsage (1, false); // short option list and exit!

}

SxCLI::CliArg& SxCLI::option (const SxString& mark,
                              const SxString& shortDescription,
                              const SxString& longDescription)
{
   SX_CHECK (mark.getSize () != 0);
   SX_CHECK (shortDescription.getSize () != 0);
   SX_CHECK (longDescription.getSize () != 0);
   
   optionList << CliArg (this, mark, shortDescription, longDescription);

   return optionList.last ();
}

SxCLI::CliArg& SxCLI::option (const SxString& mark,
                              const SxString& longDescription)
{
   SX_CHECK (mark.getSize () != 0);
   SX_CHECK (  (mark != "..." && longDescription.getSize () != 0)
            || (mark == "..."));
   
   optionList << CliArg (this, mark, "", longDescription, true, true);

   return optionList.last ();
}


SxCLI::CliArg& SxCLI::argument (const SxString& shortDescription,
                                const SxString& longDescription)
{
   SX_CHECK (shortDescription.getSize () != 0);
   SX_CHECK (longDescription.getSize () != 0);
   
   addAutoOptions ();
   optionList << CliArg (this, "", shortDescription, longDescription, false);
   return optionList.last ();
}

// ------- grouping

int SxCLI::newGroup (const SxString &groupName, const SxString &descr)
{
   SX_CHECK (currentGroup);
   int newId = currentGroup->id + 1;
   CliGroup group (newId);
   // get free number
   while (groups.contains (group)) 
      group.id++;
   group.name = groupName;
   group.descr = descr;
   groups << group;
   currentGroup = &groups.last ();
   return currentGroup->id;
}

void SxCLI::excludeGroup (int id)
{
   SX_CHECK (currentGroup);
   // get position of group id
   int pos = int (groups.findPos (CliGroup (id)));
   if (pos == -1)   {
      cout << "BUG in command line parser: group " << id << " is not defined.";
      cout << endl;
      cout << "Only previously defined options can be excluded." << endl;
      cout << "Fix command line parsing in 'main'." << endl;
      SX_EXIT;
   }
   CliGroup &group1 = *currentGroup;
   int currentGroupPos = int(groups.findPos (group1));
   CliGroup &group2 = groups(pos);
   group1.excludedGroups.append (pos);
   group2.excludedGroups.append (currentGroupPos);
   if (group1.status == Required && group2.status == Required)  {
      exclusionError (group1, group2);
   } else {
      if (group1.status == Required) group2.status = Forbidden;
      if (group2.status == Required) group1.status = Forbidden;
   }
}

void SxCLI::requireGroup (int id)
{
   SX_CHECK (currentGroup);
   // get position of group id
   int pos = int (groups.findPos (CliGroup (id)));
   if (pos == -1)   {
      cout << "BUG in command line parser: group " << id << " is not defined.";
      cout << endl;
      cout << "Only previously defined options can be required." << endl;
      cout << "Fix command line parsing in 'main'." << endl;
      SX_EXIT;
   }
   CliGroup &group1 = *currentGroup;
   CliGroup &group2 = groups(pos);
   group1.requiredGroups.append (pos);
   if (group1.status == Required)  {
      group2.setRequired (group1.initOption);
   }
}

bool SxCLI::CliGroup::setRequired (CliArg* currentOption)
{
   // if first option of this group, set group.initOption
   if (initOption == NULL) initOption = currentOption;
   if (status == MissingRequired)  {
      // --- issue error if a required option of this group has not been given
      if (currentOption->parent->mode == Normal)  {
         initOption->printMark ();
         cout << " is required if ";
         currentOption->printMark (false);
         cout << " is given." << endl;
      }
      currentOption->parent->setError ();
      return false;
   }
   if (status == Forbidden)  {
      currentOption->parent->exclusionError (*currentOption->group,*this);
      return false;
   }
   status = Required;
   // --- check exclusions
   if (!checkExclusions ()) return false;
   // recurse into required groups that are not yet in status Required
   if (requiredGroups.getSize () == 0) return true;
   for (SxList<int>::Iterator it = requiredGroups.begin ();
         it != requiredGroups.end ();
         it++)
   {
      if (currentOption->parent->groups(*it).status != Required)  {
         if (!currentOption->parent->groups(*it).setRequired (currentOption))
            return false;
      }
   }
   return true;
}

void SxCLI::exclusionError (const CliGroup &g1, const CliGroup &g2)
{
   SX_CHECK (g1.initOption);
   SX_CHECK (g2.initOption);
   if (g1 == g2)  {
      // --- Attempt to use an option of a forbidden group
      // => print out all conflicts
      for (SxList<int>::ConstIterator it  = g1.excludedGroups.begin ();
                                      it != g1.excludedGroups.end (); it++)
      {
         const CliGroup &gEx = groups (*it);
         if (gEx.status == Required) exclusionError (g1, gEx);
      }
      return;
   }

   g1.initOption->printMark ();
   cout << " can't be used together with ";
   g2.initOption->printMark (false); // "option..." instead of "Option..."
   cout << ".\n";
   setError ();
}

void SxCLI::setGroup (int id)
{
   // get position of group id
   ssize_t pos = groups.findPos (CliGroup (id));
   if (pos == -1)  {
      newGroup ();
      return;
   }
   currentGroup = &groups(pos);
}

bool SxCLI::groupAvailable (int id)
{
   ssize_t pos = groups.findPos (CliGroup (id));
   if (pos == -1)   {
      cout << "BUG in command line parser: group " << id << " is not defined.";
      cout << endl;
      SX_EXIT;
   }
   return (groups(pos).status == Required);
}


void SxCLI::version (const SxString &versionText_)
{
   versionText = versionText_;
}

void SxCLI::runAtExitFuncs ()
{
   SxList<void (*)()>::Iterator it = atExit.begin ();
   // execute all functions in list
   for ( ; it != atExit.end (); it++)  (*(*it))();
}
