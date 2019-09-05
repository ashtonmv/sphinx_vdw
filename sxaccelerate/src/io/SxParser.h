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

#ifndef _SX_PARSER_H_
#define _SX_PARSER_H_

#include <SxSymbolTable.h>
#include <SxIO.h>
#include <SxPtr.h>
#include <SxMap.h>
#include <SxString.h>
#include <SxFile.h>
#include <SxUniqueList.h>
#ifdef WIN32
#  include <sxfileno.h>
#endif

/**  \brief SPHInX input file parser for hierarchical ASCII files.

     The SPHInX input file format is similar to XML, the extensible
     markup language (http://www.xml.org). It suppports:
     - hierarchical representation of data
     - nesting files using a 'include' statement
     - including files from remote host using SSL network connections
     - file validation with type definition files (grammar definition)
     - simple arithmetics '+', '-', '*', '/', etc.
     - definition of vectors and matrices
     - event and hierarchy oriented programming interface
       (in contrast to XML's separation of SAX and DOM)
     - C-like input search path

     \par Example:
\code
   SxParser parser;
   SxParser::Table tree = parser.read ("input.sx");
\endcode

  \author  Sixten Boeck, boeck@mpie.de
  */
class SX_EXPORT_IO SxParser
{
   public:

      /** \brief Shortcut for result of the SxParser::read function */
      typedef SxConstPtr<SxSymbolTable>  Table;
      
      /** \brief Create a new parser

          The default constructor opens a new parser. The object tree can
          get retreived with the SxParser::read function.  */
      SxParser ();

      /** \brief Destructor

          This function cleans up memory used during parsing. */
      ~SxParser ();

      /** \brief Read an input file.

          This is the main function of the SPHInX input file parser for
          ASCII files. The only important input paramerter is the file 
          name which ought to be parsed. It returns an object containing
          the read input data hierarchically as a SxSymbolTable object.

          \par Example:
\code
   SxParser parser;
   SxParser::Table tree = parser.read ("input.sx");
   try  {
      const SxSymbolTable *myGroup = tree->getGroup("myGroup");
      cout << "my value = " << myGroup->get ("myValue")->toDouble();
   } catch (SxException e)  {
      e.print ();
   }

   // --- both tree and myGroup do not need to be cleaned up!
\endcode 
      More information about reading files can be found in SxSymbolTable.

      It is possible to provide a default validator which is taken when
      the user has not provided one using the "format" statement. This 
      feature is useful if subgroups are to be read out. 

      \b Attention: This function is \b not reantrant! That means that
      SxParser::read cannot be used from different threads simultaniously!
      This lack is due to the underlying lexical analyzer (lex) and
      compiler compiler (yacc).

      \param file        filename of the input file, e.g. "input.sx"
      \param validator_  default validator if user has none provided,
                         e.g. "std/myGrammar.std"
      \return            The hierarchically parsed data of the input file.
      \sa SxSymbolTable */
      SxParser::Table read (const SxString &file,
                      const SxString &validator_="");
      /** \brief return a list of included files

          This function returns a list of all files which have been included
          in the input file. This list is evaluated during the execution of
          SxParser::read.

          \par Example
\verbatim
   SxParser parser;
   parser.read ("input.sx");
   SxList<SxFile> files = parser.getNestedFiles ();
   SxList<SxFile>::Iterator it;
   for (it = files.begin(); it != files.end(); it++)  {
      cout << (*it)->getName() << endl;
   }
\endverbatim
      */
      SxList<SxFile> getNestedFiles () const;
      
      SxFile locateFile (const SxSymbol *) const;
      SxFile locateFile (const SxSymbolTable *) const;
      int locateLine (const SxSymbol *) const;
      int locateLine (const SxSymbolTable *) const;

      /** \brief Enable/disable verbose mode of the parser

          In the verbose mode every line is being printed after reading
          and before parsing. This feature is useful for debugging only.
          Note, that the verbose mode of the parser can also be controlled
          directly from the input file by adding a line
\verbatim
   verboseMode on;
\endverbatim
      */
      void setVerbose (bool verbose=true);

      /** \brief Return SPHInX search path

          The SPHInX search path is a set of paths where the parser looks
          during parser whenever a nested file is should be included via
          'include \<somefile\>;'. It is preset at compile time and can be
          changed at runtime with the SX_INCLUDE_PATH parameter
       */
        static SxString getDefaultSearchPath ();
        SxString getSearchPath () const;

        void setSearchPath (const SxString &);


         /** \brief Switch validation mode

             This function can switch on/off the validation of the input
             file. Note, this function should be used for debugging purposes
             only.*/
         void setValidation (bool state=true);
         
         /** \brief Exclude group from parsing.
         
             The name can be comma separated list of group names. */
         void excludeGroup (const SxString &name);


   protected:

      /** \brief if this flag is set to true no validation is performed. */
      bool freeFormat;

      SxString searchPath;
      
      SxUniqueList<SxString> excludeGroups;
};



/** This buffer contains the complete string to parse */
extern SxString SxParser_buffer;
/** current line number of parser */
extern int SxParser_lineNumber;
/** The filename of the file which is being parsed. 
    Used for output in SxParser_error */
extern SxString SxParser_filename;
/** Verbose mode for debugging. Prints every line while reading. 
    The verbose mode can easily be switched on by adding the line
\verbatim
   verboseMode on;
\endverbatim
    in the input file.
 */
extern bool SxParser_verbose;
/** The dimension of the tensor recently read in. */
extern int dimension;
/** The complete symbol stack */
extern SxList<SxSymbol> SxParser_symbolStack;
/** group level counter */
extern int SxParser_groupLevel;
/** List of all included files */
extern SxList<SxFile> SxParser_nestedFiles;
/** Prototype of the lexical analyzation routine */
extern int SxParser_lex ();
/** set of global variables */
extern SxList<SxString> SxParser_globalVars;
/** Prints out the IO error message and stops the program. */
void SxParser_error (const SxString &);
/** The main entrance point of the parser */
void SxParser_parseFile (const SxString &filename, 
                         const SxString &path="",
                         const SxUniqueList<SxString> &
                               excludeGroups=SxUniqueList<SxString>());
/** @{
    \brief finds a string somewere in the search path */
SX_EXPORT_IO SxList<SxString> SxParser_findInPath (const SxString &);
SX_EXPORT_IO SxList<SxString> SxParser_findInPath (const SxString &, 
                                                   const SxString &);
SX_EXPORT_IO SxList<SxString> SxParser_findInPath (const SxString &, 
                                                   const SxList<SxString> &);
/**
  @} 
*/


void SxParser_init ();




void initParser ();


#endif /* _SX_PARSER_H_ */
