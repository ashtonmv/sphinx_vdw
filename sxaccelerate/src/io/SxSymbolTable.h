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
#ifndef _SX_SYM_TABLE_H_
#define _SX_SYM_TABLE_H_

#include <stdio.h>
#include <math.h>
#include <SxIO.h>
#include <SxList.h>
#include <SxMap.h>
#include <SxMultiMap.h>
#include <SxArray.h>
#include <SxString.h>
#include <SxException.h>
#include <iostream>


/**
  @ingroup Tools
  */
class SX_EXPORT_IO SxSymbol
{
   public:

      enum SymbolType { Unknown, Number, String, Vector, Variable, 
                        Function, StringFunction, List };

      /** \brief Define enum types for defined and undefined symbols*/
      enum Validity { Undefined, Defined };

      enum AutoCast { NoAutoCast, DoAutoCast };
      
      SxString name;
      int      type;
      /** \brief Keep track whether symbol has been defined already with a 
                 valid value */
      Validity status;
      bool     initialized;
      bool     binary;  // included via HERE documents
      
      SxString parserFilename;
      int      parserLineNumber;

      double               val;
      SxString             str;
      double             (*func)(double);
      SxString           (*strfunc)(const SxString &);
      SxList<SxSymbol>    *valList;  // [1, 2, 3] or [[..], [..]]
      
      SxSymbol ();
      SxSymbol (const SxString &name, double val);
      SxSymbol (const SxString &name, const SxString &str, 
                enum AutoCast autoCast = NoAutoCast);
      SxSymbol (const SxString &name, double (*func)(double));
      SxSymbol (const SxString &name,
                SxString (*strfunc)(const SxString &));
      SxSymbol (const SxSymbol &);
      ~SxSymbol ();

      SxSymbol &operator= (const SxSymbol &);
      
      SxSymbol  operator+ (const SxSymbol &);
      SxSymbol  operator- (const SxSymbol &);
      SxSymbol  operator* (const SxSymbol &);
      SxSymbol  operator/ (const SxSymbol &);
      SxSymbol  operator^ (const SxSymbol &);

      void append (double);
      void append (const SxString &);
      void append (SxSymbol &);
      void prepend (SxSymbol &);
      void print () const;

      SymbolType      getType () const;
      const SxString &getName () const;
      /** \brief Returns true if the symbol has been assigned with a value */
      bool            isDefined () const;
      /** \brief Controls the variable SxSymbol::status.

          This function can be used when a symbol should explicitly be set
          as valid (defined) or invalid (undefined). */
      void setDefined (bool defined=true);

      SxSymbol flatten (); 
      SxList<double>   toList () const;
      SxList<int>      toIntList () const;
      SxList<SxString> toStringList () const;
      const SxString  &toString () const;
      bool             toBool () const;
      bool             toAttribute () const;
      int              toInt () const;
      double           toReal () const;
      int      getRank () const;
      SxList<int> getDimensions () const;

      /** \brief cause an error message for using undefined variables

          This function prints an error message and exits. It is used
          to be called from operators acting on undefined symbols. */
      void undefVarError (const SxSymbol &) const;

};

std::ostream &operator<< (std::ostream &s, const SxSymbol &in);
SxSymbol operator- (const SxSymbol &in);  // unary minus

/**
  After parsing an SX-input file an object of this class contains
  the complete node tree (similar to a DOM tree, see http://www.xml.org).
  Each node is either a SxSymbolTable object again or a SxSymbol.

  \parUsage:
  This sections contains examples to
     -# parsing a node tree
     -# reading optional parameters
     
  (1) Let's assume, you would like to read an input file that contains 
  an atomic structure:
  \verbatim
    system  {
       species  {
          name = "A";
          atom  {  coord = [0.0, 0.0, 0.0];
          atom  {  coord = [0.5, 0.5, 0.5];
       }
       species  {
          name = "B";
          atom  {  coord = [0.5, 0.0, 0.0];
          atom  {  coord = [0.0, 0.0, 0.5];
       }
    }
  \endverbatim
  To read in such a tree your source code should look like
  \verbatim
     SxSymbolTable *system, *species;
     SxList<double> coords;
     try  {
        system = symbolTable->getGroup ("system");
        for (species  = system->getGroup("species");
             species != NULL;
             species  = species->nextSibling ("species"))
        {
           cout << "Reading species " << species->get ("name")->toString() << endl;
           for (atom  = species->getGroup("atom");
                atom != NULL;
                atom  = atom->nextSibling ("atom"))
             {
                coord = atom->get("coords")->toList();
                coord.print();
             }  
        }
     }  catch (SxException e)  {
        e.print ();
     }
  \endverbatim
  \par
  (2) Sometimes it is useful to have optional parameters. In that case every
  exception would stop evaluating any further statement in the 'try' block.
  To overcome that problem catch EACH optional parameter:
  \verbatim
     double optParam1 = 3.1415;  // default value for optParam1
     try  { optParam1 = symbol-get("PI")->toReal(); }
     catch (SxException) { }

     double optParam1 = 2.71828; // default value for optParam2
     try  { optParam1 = symbol-get("E")->toReal(); }
     catch (SxException) { }
  \endverbatim

  The input file parser you can even use as a simplified interpreted
  language. Put all your commands into a group (e.g. "main"). Each
  "function" is then symbolized as a SxSymbolTable group. Lets inspect
  the following input file
  \verbatim
     main {
        steepestDescent { param1=1; param2=2; }
        conjGradient { paramA=3; paramB = 4; }
     }
  \endverbatim
  This lines of code will cause the executed program first to perfrom
  some steepestDescent algorithm followed by an conjugate gradient.
  The corresponding code fragment looks like
  \verbatim
     SxSymbolTable *main = NULL;
     try  { main = symbolTable->getGroup("main"); }
     catch (SxException e) { e.print (); exit(1); }

     // --- parse commands in sequential order
     SxSymbolTable *cmd;
     for (cmd  = main->begin();
          cmd != NULL; 
          cmd  = cmd->nextSibling ());
     { 
        cout << "About to evaluate command " << cmd->getName() << endl;
        evaluate (cmd);
     }
  \endverbatim

  Here "evaluate" is a function with a input parameter (SxSymbolTable *).

  @author Sixten Boeck, boeck@fhi-berlin.mpg.de
  */
class SX_EXPORT_IO SxSymbolTable
{
   public:
      SxMap<SxString, SxSymbol> table;
      SxSymbolTable    *parent;
      const SxString    name;
      const int         level;
      SxList<SxSymbolTable *> children;
      SxList<SxSymbol *> tmpList;  // for yyParser
      bool              deprecateFormat;
      
      SxString parserFilename;
      int      parserLineNumber;

      static SxSymbolTable *&getGlobalPtr ();

      SxSymbolTable ();
      SxSymbolTable (SxSymbolTable *parentPtr, const SxString &lvlName);
      ~SxSymbolTable ();

      /** Returns the pointer to a symbol or NULL. It must not
          be deleted after usage (it is managed by this class)! 
          @exception SxException::MISSING_SYMBOL the symbol was not found
       */
      const SxSymbol *get (const SxString &, // #exception
                           bool localOnly=false, 
                           bool throwException=true, 
                           const SxString &parentPath="") const;

            SxSymbol *get (const SxString &, // #exception
                           bool localOnly=false, 
                           bool throwException=true, 
                           const SxString &parentPath="");

      /** Returns the name of the current group. */
      SxString getName () const;

      bool hasAttribute (const SxString &) const;
      /** 
        @exception SxException::MISSING_GROUP  the group was not found 
       */
      SxSymbolTable *getGroup (const SxString &, // #exception
                               bool recursive=false) const;

      const SxList<SxSymbol> &getChildrenList ()const{return table.getValues();}
      SxSymbolTable *begin () const;
      SxSymbolTable *nextSibling (const SxString &) const;
      SxSymbolTable *nextSibling () const;
      int       getNItems     (const SxString &, bool localOnly=false) const;
      bool      contains      (const SxString &, bool localOnly=false) const;
      bool      containsGroup (const SxString &) const;
      SxSymbol *append        (const SxString &);
      SxSymbol *append        (const SxString &var,  // #exception
                               const SxString &val);
      void      append        (const SxMap<SxString,SxString> &vars); // #exception
      void      remove        (SxSymbol *);
      void      addChild      (SxSymbolTable *);

      SxString  getPathName () const;
      bool useDeprecateFormat () const { return topLevel()->deprecateFormat; }

      SxSymbolTable *topLevel ();
      const SxSymbolTable *topLevel () const;


      // push temporary list;  [[..],[..]]
      SxSymbol *pushList ();
      SxSymbol *popList ();

      // Top level wraper for profiling
      void validateTable (const SxSymbolTable &);
      
      ssize_t validate (const SxSymbolTable &);

      void printList () const;
      void printHash (const SxString &indent1=" ", 
                      const SxString &indent2="   ",
                      int idx=0) const;
      void print () const;

      void validationError (const SxString &, const SxSymbol *);
      void validationError (const SxString &, const SxSymbolTable *);
};

// defined in SxParser.lpp
extern "C" SxString sxReadFile (const SxString &);

#endif /* _SX_SYM_TABLE_H_ */
