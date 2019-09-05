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


#ifndef _SX_CLI_H_
#define _SX_CLI_H_

#ifndef SXLIBTOOL
#   define SXLIBTOOL ""
#endif

#include <SxString.h>
#include <SxConfig.h>
#include <SxUtil.h>
#include <SxList.h>
#include <SxPtr.h>
#include <SxRedirect.h>

//#include <SxVector3.h>
/// Parameter class: knows whether it is initialized or not
// Does not work for bool!
template <class T>
class Given
{
   public:
      T value;
      bool given;
      Given (bool b = false)  {given = b;}
      // inline Given (T val)  {given = true; value = val;}
      // integer constants are of type int
      // floating point number consts like 0. are of type double
      // string consts like "" are of type const char*
      // we must be able to convert these to T 
      // The ideal constructor would be
      // template <class U>
      // inline Given (U val) {given = true; value = T(val);}
      Given (int val) {given = true; value = T(val);}
      Given (double val) {given = true; value = T(val);}
      Given (const char* val) {given = true; value = T(val);}
      ~Given () {}
};

/** @brief Command Line Interface
    @author Christoph Freysoldt freyso@fhi-berlin.mpg.de

    This class provides a command line option parser. The idea is to declare
    the available options with the necessary help messages, and to use this
    information for parsing the command line options and producing a
    help message (usage).

  Each option has the following features
  - one or several alternative marks to indicate which option is supplied
    (e.g. "--waves")
  - or no mark at all, so that the option (here called argument) is defined by
    its position
  - a short description of what is expected (e.g. "file")
  - a long description of the option (e.g. "the wave function used for...")
  - a default value (e.g. "waves.sxb")
  - for numbers: lower and upper bounds
  - an option can be set to be optional, or required
  - an option can be set to have no value, i.e. boolean options
  - an option can be a multiple-choice option

  The whole parsing consists of three phases: first, the option has to be
  declared with an #option or #argument command. It returns the option
  CliArg which is used in phase 2 to extract and convert the
  given argument to the desired
  type. In this phase, the option behaviour (optional/empty) is determined
  from the type. Options without a default value are considered to be
  required.
  In phase 3, initiated by the #finalize command, the parsing is completed:
  if unknown options are left or any errors in the parsing occured,
  the help message is printed and the program exits.

  \par Examples:

-  A simple one to show how simple it is to define options.
\code
#include <SxCLI.h>
#include <SxString.h>

int main (int argc, char** argv)
{
   SxCLI cli (argc, argv); // Command Line Parsing

   // init vars from options
   int i      = cli.option ("-i|--integer",
                            "number",
                            "an integer number").toInt(0,0,100);
   float f    = cli.option ("-f|--float",
                            "number",
                            "a real number").toFloat (0.);
   SxString s = cli.option ("-s",
                            "string",
                            "some string").toString("string1");
   bool sayHi = cli.option ("--hi", "prints hi").toBool ();
   double d   = cli.option ("--double",
                            "number",
                            "a very exact number").toDouble ();
   // d is required, because no default is given

   // arguments are always last!
   SxString command = cli.argument ("command",
                                    "some imaginary action").toString ();

   cli.finalize ();

   cout << "Integer: " << i << endl;
   cout << "Float: "   << f << endl;
   cout << "Double: "  << d << endl;
   cout << "String: "  << s << endl;
   if (sayHi) cout << "Hi!" << endl;
   cout << "Command: " << command << endl;
   cout.flush ();
}
\endcode

- Multiple-choice options
\code
#include <SxCLI.h>
#include <SxString.h>

int main (int argc, char** argv)
{
   SxCLI cli (argc, argv);

   cli.authors = "C. Freysoldt";
   cli.preUsageMessage = "A stupid, yet illustrative CLI example.";

   enum Number { Two, One, Four, Three, None };

   SxList<enum Number> translateNumber;
   translateNumber << Two << One << Four << Three;

   enum Number n;
   int pos = cli.option("--one|--two|--three|--four", "a number description")
             .toChoice ();

   if (pos < 0)
      n = None
   else
      n = (SxList<enum Number> () << One << Two << Three << Four)(pos);

   cli.finalize ();

   switch (n)  {
      case One   : cout << "1"; break;
      case Two   : cout << "2"; break;
      case Three : cout << "3"; break;
      case Four  : cout << "4"; break;
      case None  : cout << "no number"; break;
   }
   cout << endl;
}
\endcode

- Manual parsing of options.
Non-standard type initialization requires usually a manual conversion
of the option into the necessary information. Do NEVER exit during the
parsing, rather issue a parser error:
\code
int main (int argc, char** argv)
{
   SxCLI cli (argc, argv);
   SxString code = cli.option ("-x",
                               "code",
                               "a code: a letter followed by a number, e.g. l2")
                   .toString (""); // the "" suppresses the requiredness

   // set the default description (can be anything)
   cli.last ().defaultValue = "(default: determine code from weather)";

   // --- manual parsing
   char codeLetter; int codeNumber = -1;
   bool useCode = code.getSize () != 0;
   if (useCode)  {
      if (myString.getSize () != 2)  {
         cli.last ().printMark (); // -> "Option -x"
         cout << ": code '" << code << "' has wrong length; must be 2" << endl;
         // issue error
         cli.setError ();
      } else {
         codeLetter = code(0);
         codeNumber = code(1,1).toInt ();
         // + check that codeLetter is a letter,
         // and codeNumber is a valid number
      }
   } else {
      // code from weather
   }

   cli.finalize ();
   if (useCode) cout << "Code letter: '" << codeLetter << "'" << endl;
}
\endcode

- Two ways of running a program, with different options
\code
int main (int argc, char** argv)
{
   SxCLI cli (argc, argv);

   cli.preUsageMessage =
      "This program can be run in two modes: the INTEGER mode requires the"
      "options -i and -j, the (default) REAL mode requires the option -f.";

   int intGroup = cli.newGroup("integer mode");

   // If any of the following two options is given, the other one
   // is considered to be required. If none is given, the complete
   // group is missing and none of the options is expected. Optional
   // options would be declared in a standard way (by a default)
   int i = cli.option ("-i", "integer", "a number").toInt ();
   int j = cli.option ("-j", "integer", "a second number").toInt ();

   int realGroup = cli.newGroup("real mode");
   cli.excludeGroup (intGroup);

   double f = cli.option ("-f", "real number", "the real number")
              .toDouble (0.);

   cli.setGroup(cli.generalGroup);

   SxString format = cli.option ("--name","name","a name").toString ("");

   cli.finalize ();

   bool intMode = cli.groupAvailable (intGroup);
   if (intMode)
      cout << "i = " << i << "; j = " << j << endl;
   else
      cout << "f = " << f;

   if (name.getSize () > 0)
      cout << "name = " << name << endl;

}
\endcode

  More examples can be found in the existing add-ons.

  Hidden automatic options:
  - --completion   suggests possible completions of the current command line
  - --sxprintcli     print out the command line interface in sx-format
  */
class SX_EXPORT_UTIL SxCLI
{
   public:
      class SX_EXPORT_UTIL CliGroup;
      /// Option class
      /**  The SxString and number conversions take the optional default
           value as parameter. If it is not given (and the Given class is used
           to find out about that) the option is required.
           The number conversions take analogous min and max values. If you
           want to set a range without a default, you can use false as value
           for the default.
           @sa Given
        */
      class SX_EXPORT_UTIL CliArg  {
         public:
            /** \brief Constructor
              @param parentIn pointer to SxCLI object
              @param m marks (separated by |)
              @param s short description
              @param l long description
              @param optIn true if optional
              @param eoIn  true if option is a flag
              */
            CliArg (SxCLI* parentIn,
                    SxString m = "",
                    SxString s = "",
                    SxString l = "",
                    bool     optIn = true,
                    bool     eoIn = false);
            /// Constructor
            CliArg () : optional(false) { }
            /// Destructor
            ~CliArg () { /* empty */ }

            /// The option symbols (e.g. "--help")
            SxList<SxString> mark;
            /// Description for command line and errors
            SxString shortDescription;
            /// Description for printUsage ()
            SxString longDescription;
            /// Contains the description of default and range 
            SxString defaultValue;
         protected:
            /// Contains optional value for flag-like options
            SxString optValue;
         public:
            /// Does the flag have a value
            bool hasValue () const { return optValue.getSize () > 0; }
            /// Get optional value
            SxString getValue () const {
               // getValue must only be used for flags with optional values
               SX_CHECK (emptyOption);
               // first character is always '='
               return optValue.subString (1);
            }

            /// \name Behaviour of option
            //@{
            /// True if option is not required
            bool optional;
            /// True if parameter after mark is needed (default)
            bool emptyOption;
            /// For completion mode:
            bool used;
            /// Set option non-optional
            /** @param req If true, the option is non-optional, i.e. must be 
              given. The general default for options is optional, it can be 
              changed with this option.
              \example
\code
   SxList<int> idxList
      = cli.option ("--list", "list",
                    "comma-separated list of indices or index ranges a-b")
        .required ().toIdxList ();
\endcode
             */
            inline CliArg& required (bool req = true)  {
               optional = !req;
               return *this;
            }

            /// Sticky options keep their values when a new loop is started
            bool stickyOption;

            /// Sticky options keep their values when a new loop is started
            inline CliArg &sticky(bool stickyIn = true)  {
               // no sticky arguments in loops -> would always be propagated
               SX_CHECK(!stickyIn || mark.getSize () > 0);

               stickyOption = stickyIn;
               return *this;
            }
            //@}


            /// \name Auxiliary functions
            //@{
            /// Get the argument (or used option for empty options)
            SxString getArgument ();
            /** \brief Print mark (or short description)
              @param capital if true, print "Option", else "option"

              This is the way the option is identified to the user
              in printouts. If several marks are possible, it prints
              the first and adds " or alternatives". If there is no
              mark (argument) the short description is used instead.
              */
            void printMark (bool capital = true) const;
            /// Print out option in sx-style
            void printout () const;
         protected:
            /// Print 'val lower than min' error message, set error
            void tooSmallError (const SxString&, const SxString&);
            /// Print 'val larger than max' error message, set error
            void tooLargeError (const SxString&, const SxString&);
            /// Set defaultValue to default and range description
            void setDefault (const SxString&, const SxString&, const SxString&);
            //@}
         public:
            /// @name Argument type conversion
            //@{
            /// Convert to int
            int toInt (Given<int> def = false,
                       Given<int> min = false,
                       Given<int> max = false);
            /// Convert to long int
            long int toLong (Given<long int> def = false,
                             Given<long int> min = false,
                             Given<long int> max = false);
            /// Convert to int64
            int64_t toInt64 (Given<int64_t> def = false,
                             Given<int64_t> min = false,
                             Given<int64_t> max = false);
            /// Convert to float
            float toFloat (Given<float> def = false,
                           Given<float> min = false,
                           Given<float> max = false);
            /// Convert to double
            double toDouble (Given<double> def = false,
                             Given<double> min = false,
                             Given<double> max = false);
            /// Convert to SxString
            SxString toString (Given<SxString> def = false);
            /** \brief Returns id for multiple-flag options
              \param def default choice.
              \return the id of the matching flag, starting from 0.
                      If no flag is on the command line, the return value is
                      -# the default if any OR
                      -# -1 if the option is optional (see #required) OR
                      -# 0 if the option is required (stops at finalize)

              \par Example
              \code
int idFunctional = cli.option ("--lda|--pbe|--pbe_lda")
                   .required ()
                   .toChoice ();
SxXC::XCFunctional functional = (SxList<SxXC::XCFunctional> ()
                                << SxXC::LDA     // --lda     => 0
                                << SxXC::PBE     // --pbe     => 1
                                << SxXC::PBE_LDA // --pbe_lda => 2
                                )(idFunctional);
              \endcode

              The example shows also how the id can be translated to any
              type. The return value of 0 for missing required option supports
              this feature.
              */
            int toChoice (Given<int> def = false);
            /// Convert to SxList<SxString> (multiple options)
            SxList<SxString> toList ();
            /** This converts a comma-separated (or colon-separated) list
                of indices or index-ranges a-b into a list of int. As
                human-readable indices start from 1, while C++ starts from 0,
                there is a human-to-C++ conversion, i.e. 1 is subtracted from
                each index.
                \example
                "1,4-7" is expanded into internal 0,3,4,5,6
                @note The default is the empty list, so this option is optional
                by default. Use #required to override this behaviour.
              \brief Convert to SxList<int>
              */
            SxList<int> toIdxList ();
            /** This converts a comma-separated (or colon-separated) list
            of indices or index-ranges a-b into a list of 64-bit integers. As
            human-readable indices start from 1, while C++ starts from 0,
            there is a human-to-C++ conversion, i.e. 1 is subtracted from
            each index.
            \example
            "1,4-7" is expanded into internal 0,3,4,5,6
            @note The default is the empty list, so this option is optional
            by default. Use #required to override this behaviour.
            \brief Convert to SxList<int>
            */
            SxList<int64_t> toIdxList64 ();
            /** \brief Convert to 3-list of doubles

              @param separators each character is a possible separator
              @param braces     possible pairs of braces
                This is typically used for reading SxVector3<Double>.
                \code
   SxVector3<Double> vec(cli.opt("--vec","vector","a vector").toList3 ());
                \endcode

                The input format of a vector (v1,v2,v3) is
                <open-brace> v1 <sep> v2 <sep> v3 <close-brace>

                The braces can be any of the pairs () {} [] <>, or
                can be omitted. However, they must match.

                The separator can be ',' or '|' or or ':' or just a space.
                Again, the separators must be all the same.

                These defaults can be changed by the arguments.

                This conversion cannot have a default, and is required.
                If you want to implement a default, do it this way
                \code
SxCLI::CliArg *opt;

// --- declare option
opt = cli.option ("--vector","vector","a 3D vector");
opt->defaultValue = "(default: {1, 2, 3})";  // set printed default

// --- parse option
SxVector<Double> vec(1., 2., 3.);            // set internal default
// override default by command line option if any
if (opt->exists ()) vec = SxVector3<Double> (opt->toList3 ());

// make it optional (AFTER toList3 () !)
opt->required (false);
                \endcode

                @return The return value is an SxList<double> with three
                        valid elements, which may be
                        converted into a SxVector3<Double>.
              */
            SxList<double> toList3 (const SxString &separators = ",|:",
                                    const SxString &braces = "(){}[]<>");
            /** \brief Convert to 3-list of integers

                @note This is a quick hack with an intermediate
                      transformation to doubles. SxCLI therefore doesn't
                      complain about non-integer numbers here and rounds them
                      to integers.
                @note For detailed description: see #toList3
                @return The return value is an SxList with three
                        valid elements, which may be 
                        converted into a SxVector3<Int>.
            */
            SxList<int> toIntList3 (const SxString &separators = ",|:x",
                                    const SxString &braces = "(){}[]<>");
            /** \brief Convert to list of doubles
              @param separators string of possible separator characters.
                                For example ":," would split on ':' if
                                any ':' present, or split on ',' if not

              Converts argument to a list of doubles. The numbers have to
              be separated by one of the characters in separators argument.


              You may use this for reading in ranges. The recommended way
              for this is
              \code
double min, max;
{  // this brace makes minMax a local object

   SxList<double> minMax
      = cli.option("-r|--range","min:max","a colon-separated range")
     // .required ()   // uncomment if there's no default range
        .toDoubleList ();

   // Set defaults
   // cli.last ().defaultValue = "default: ...";
   // min = ...;
   // max = ...;

   if (!cli.error)  { // if there's an error, suppress all further parsing
      if (minMax.getSize () == 2)  {
         min = minMax(0);
         max = minMax(1);

         // check range
         if (min > max)  {
            cout << "Illegal range: min exceeds max" << endl;
            cli.setError ();
         }
      } else {
         cout << "Illegal range format: must be <min>:<max>" << endl;
         cli.setError ();
      }
   }
}
             \endcode

              */
            SxList<double> toDoubleList (const SxString &separators = ":,;/");
            /** \brief Convert to list of integers

              @param separators string of possible separator characters.
                                For example ":," would split on ':' if
                                any ':' present, or split on ',' if not

              Converts argument to a list of integers. The numbers have to
              be separated by one of the characters in separators argument.

              \note For indices, toIdxList may be more convenient.
              \sa toIdxList

              */
            SxList<int> toIntList (const SxString &separators = ":,;/");

         protected:
            /** \brief Convert to bool or look for option
               @param iLoop loop id (-1 for current loop of parent)
               @param keep set true if you just want to know about the
              existence of that option in loop iLoop but convert it later.
              @note This is the internal routine.
             */
            bool exists (int iLoop, bool keep = true);
         public:
            /** Convert to bool or look for option
                @param keep set true if you just want to know about the
              existence of that option but convert it later.
             */ 
            bool exists (bool keep = true)  { return exists(-1, keep); }
            /// Convert to bool
            bool toBool ()  { return exists (-1, false); }
            //@}
            /// The group the option belongs to
            CliGroup *group;
            /// Does mark fit? if so, give length of mark, if not 0
            int fits (const SxString&) const;
         protected:

            /// The SxCLI object to which this CliArg belongs.
            SxCLI *parent;

            friend class CliGroup;
      }; // class CliArg
      // Give CliArg access to our protected members
      friend class CliArg;
      enum GroupStatus {Optional, Required, MissingRequired, Forbidden};
      class SX_EXPORT_UTIL CliGroup  {
         public:
            /// Constructor
            explicit CliGroup (int idNumber = -1);
            /// Destructor
            ~CliGroup () {}

            /// Current status
            GroupStatus status;

            /// Identification number
            int id;
            /// Description
            SxString name;
            /// Long description
            SxString descr;
            /// Prints group description
            void printName ();
            /// Print out group in sx-style
            void printout () const;
            /// List of incompatible groups
            SxList<int> excludedGroups;
            /// List of required groups
            SxList<int> requiredGroups;

            /// Option that determined the status (needed for error messages)
            CliArg* initOption;

            /// Comparison (compare id)
            inline bool operator== (const CliGroup& x) const {
               return (id == x.id);
            }
            /** @return true if compatible, false if not. 
              \brief Issue error if incompatible with other options
              */
            bool checkExclusions ();

            /** \brief Set this group to required
                @return false if there is a problem
            */
            bool setRequired (CliArg* currentOption);
      };

      /// Message to be printed as first line(s) in printUsage
      SxString preUsageMessage;
      /// Program name (determined automatically, but can be overwritten)
      SxString programName;
      /// Authors shown in version
      SxString authors;

   protected:
      /// List of options
      SxList<CliArg> optionList;
      /// Mode
      enum Mode {Normal, PrintHelp, Completion, PrintSxFormat} mode;
      /// Pointer to help option (so short-style printUsage can refer to it)
      CliArg *helpOptionPtr;
      /// original command line arguments
      SxList<SxString> args;

   public:
      /// Constructor: pass argc and argv (T = char and/or wchar_t)
      template<class T>
      inline SxCLI (int argc_, T **argv_, const SxString & = SXLIBTOOL);
      /// Destructor
      virtual ~SxCLI ();

      /**
        First, the preUsageMessage is printed. Then the usage help is printed
        for all options given so far. You should never call this before you
        haven't declared all options.

        \param errorCode
        Usually, printUsage calls exit with errorCode.
        Only if errorCode is -1, the program continues.
        \param printLong
        If set to false, only the abbreviated option list is printed.
        @brief Prints usage
       */
      void printUsage(int errorCode = -1, bool printLong = true);

      /** \brief Print possible continuations of command line
        */
      void printCompletion ();

      /// Print out CLI in sx-style
      void printout () const;

      /** \brief print compiler and link options

          Whenever a SFHIngX program is called with '--opts' detailed
          information about the compilation and linkage is printed and
          the program stops.

       */
      virtual void printCompOptions() const;

      /** What to do with unknown arguments */
      enum ExtraArgs { IgnoreExtraArgs, ComplainExtraArgs };
      /**
        \param ExtraArgs What to do with unknown arguments.
        Use "IgnoreExtraArgs" if they are parsed by other library routines.
        This is incompatible with arguments (options without mark).

        \brief
        ALWAYS last: Check if all arguments are parsed and no error occured
        */
      void finalize (enum ExtraArgs = ComplainExtraArgs);
      /// True if an error occured
      bool error;
      /// The arguments that are not yet parsed (iLoop:iArg)
      SxList<SxList<SxString> > arguments;

   protected:
      /// Internal loop counter
      int iLoop;
   public:
      /// Get loop id
      int getILoop () const { return iLoop; }

      /// Check for more loops
      bool looping ();

      /// Define argument splitting the argument list into loop blocks
      void setLoopSeparator (const SxString &mark);

      bool stickyDefault;

      /// Set error true
      inline void setError () { error = true; }

   /// \name Declaring options   
   //@{
      /// Declare an option
      CliArg& option (const SxString& mark,
                      const SxString& shortDescription,
                      const SxString& longDescription);
      /// Declare an (boolean) option without argument
      CliArg& option (const SxString& mark,
                      const SxString& longDescription);
      /// Declare an argument (no mark)
/**
  Arguments are options without marks. They have a slightly different behaviour
  than "marked" options in that they take the first argument left. They must 
  be given after all the marked options. They are non-optional by default, but 
  this can be changed by required (false).
  @sa #option
  @sa CliArg::required
*/
      CliArg& argument (const SxString& shortDescription,
                        const SxString& longDescription);
/**
  \brief Get last option again

  This function is used to get the last option again to finetune its
  behaviour. Its main use is for explaining default behaviour.
  \example
  \code
  SxString storeFile = cli.option ("-s","file","store result").toString ("");
  // now the defaultValue is "", which isn't printed
  cli.last ().defaultValue = "don't store results";
  \endcode
  */
      CliArg& last ()  { return optionList.last (); }
/**
     SxCLI tries to generate some version information automatically.
     By default, it tells the SFHIngX release it comes from, e.g.

     sxaddon from SPHInX 1.1.0

     Note that the release is that from SxCLI (assuming that you
     recompile the complete code frequently).

     If you don't like that or have your own version numbering,
     you have to trigger the version manually

     \code
cli.version ("1.3a");
     \endcode

     @brief Automatic version number
*/
      void version (const SxString &versionText_);
   ///@}
   public:

      /// Return number of arguments
      int getArgc () const;
      /// Return command line arguments
      const char **getArgv () const;

   protected:
      /// The list of groups
      SxList<CliGroup> groups;

      /// The SxCLI auto-option group
      CliGroup autoOptionGroup;

      /** \brief Option names reserved for SxCLI

        SxCLI automatically adds default options. This list
        contains the namespace reserved for this purpose. You
        may employ some options for your own purposes; SxCLI
        will try to employ alternative names for its auto options.

        */
      SxList<SxString> autoOptionNames;

      /// Handle the automatic options
      void addAutoOptions ();

      /// the actual constructor.
      virtual void init (const SxList<SxString> &args_);

      /// Initialize options/groups
      void init ();

      /// Issue error because options are incompatible
      void exclusionError (const CliGroup &g1, const CliGroup &g2);


      /// Version text
      SxString versionText;
      /// Contains the file compilation string given from the make environment
      SxString dateStr;
      /// Path to libtool
      SxString libtool;
      /// Path to executable (or libtool wrapper)
      SxString execPath;
      /// Directory of the executable (or libtool wrapper)
      SxString execDir;
      /// SPHInX release
      SxString release;
      /// Package name (if other than SPHInX)
      SxString package;
      /// Copyright information
      SxString copyright;

   public:
      /// Current group
      CliGroup *currentGroup;
      /// The standard group
      static const int generalGroup;
   /// \name Grouping options
   //@{
      /**
        @param name group name (or description), i.e. "energy mode"
        @return the new group id
        \brief Start new group
        */
      int newGroup (const SxString &groupName = "",
                    const SxString &descr = "");
      /// Declare current group as exclusive alternative to some other group
      void excludeGroup (int id);
      /// Declare some other group required if current group is used
      void requireGroup (int id);
      /// Set current group
      void setGroup (int id);
      /// Is any option of this group available?
      bool groupAvailable (int id);

      /// Path to executable (or libtool wrapper)
      SxString getExecPath () const { return execPath; }
      SxString getExecDir ()  const { return execDir; }
      SxString getAppName ()  const { return programName; }
   //@}

   /// Singleton class to control log file redirection
      class Log {
         private:
            friend class SxCLI;

            /// Whether logging is enabled
            bool enabled;
            /// Why logging is not enabled
            SxString reason;
            Log () : enabled(true) {}

            /// private singleton
            static Log &get ()
            {
               static Log theOne;
               return theOne;
            }
            /// Control output behavior
            SxPtr<SxRedirect> sxTee;

         public:
            /// Disable logging
            static void disable (const SxString &why);
      };

      /// Functions to execute at CLI shutdown
      static SxList<void (*)()> atExit;
   protected:
      /// Execute functions at exit
      static void runAtExitFuncs ();

};


// --- make this inline because we want all the macros
//     to be expanded from within the add-on
//     UNIX:    <T> = {char}
//     Windows: <T> = {char, wchar_t}
template<class T>
SxCLI::SxCLI (int argc_, T **argv_, const SxString &libtool_)
   : stickyDefault(true), libtool(libtool_)
{
   /// initialize the dateStr variable.
#ifdef DATE
   dateStr = SxString(DATE);
#elif defined SXDATE
   dateStr = SxString(SXDATE);
#else
   dateStr = SxString(__DATE__);
#endif
#ifdef SXPACKAGE
   package = SXPACKAGE;
   if (package.toLower () == "sphinx") package = "S/PHI/nX";
#else
   package = "S/PHI/nX";
#endif
#ifdef SXRELEASE
   release = SXRELEASE;
   versionText = "from " + package + " " SXRELEASE;
#endif
#ifdef SXCOPYRIGHT
   copyright = SxString(SXCOPYRIGHT);
#endif
   SxList<SxString>  myArgs;
   for (int iArg = 1; iArg < argc_; iArg++)
      if (sizeof(T) == sizeof(char))
         myArgs.append(SxString::fromUtf8 ((const char *)argv_[iArg]));
      else if (sizeof(T) == sizeof(uint16_t))
         myArgs.append(SxString::fromUtf16 (reinterpret_cast<const uint16_t *>(argv_[iArg])));
      else SX_EXIT; // not implemented

   init (myArgs);
}

#endif // _SX_CLI_H_
