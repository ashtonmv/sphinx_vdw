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

#ifndef _SX_SIMPLE_PARSER_H_
#define _SX_SIMPLE_PARSER_H_

#include <SxIO.h>
#include <SxParser.h>
#include <SxSymbolTable.h>
#include <SxStack.h>

// ------------------------------------------------------------------
/** \brief Utility class for macro-based input file parsing

    \b SxClass = S/PHI/nX simplified parser: group hierarchy

    SxSimpleParser is the top-level container for keeping track of the
    group hierarchy in nested parsing, and provides the "top-level"
    hidden variable.

    \author C. Freysoldt */
class SX_EXPORT_IO SxSimpleParser
{
   protected:
      /// Stack of outer groups (pushed/popped via SxParserPush)
      SxStack<const SxSymbolTable*> outerGroups;
      /// The current symbol table entry
      const SxSymbolTable *currentGroup;
      friend class SxParserPush;
      /// Location of most recent macro call
      const char *sourceFile;
      /// Location of most recent macro call
      int sourceLine;
   public:
      const SxSymbolTable* current () const { return currentGroup; }
      /// Set current position (for error message)
      void where (const char* file, int line)
      {
         sourceFile = file;
         sourceLine = line;
      }
      /// Print source location of problematic conversion
      void error () const
      {
         cout << "Error while evaluating symbol tree in "
              << sourceFile << ':' << sourceLine << endl;
      }
      /// Constructor
      SxSimpleParser (const SxSymbolTable *table, const char* file, int line)
         : currentGroup(table),
           sourceFile (file),
           sourceLine (line)
      { /* empty */ }
      /// Return true if there is something to parse
      operator bool () { return currentGroup; }
      /// get next group with same name as current one
      void next ()  {
         if (currentGroup)
            currentGroup = currentGroup->nextSibling (currentGroup->getName ());
      }
};

/** \brief Utility class for macro-based input file parsing

    \b SxClass = S/PHI/nX simplified parser: symbols

    SxParserGet is the interface for evaluating symbols -> variables

    \author C. Freysoldt */
class SX_EXPORT_IO SxParserGet {
   protected:
      SxSimpleParser &theParser;
   public:
      const SxSymbol *symbol;
      const SxString name;
      /// get double from symbol or default value
      double operator|| (double defaultValue) const;
      /// get string from symbol or default value
      SxString operator|| (const char *defaultValue) const;
      /// get string from symbol or default value
      SxString operator|| (const SxString &defaultValue) const;
      /// get int from symbol or default value
      int operator|| (int defaultValue) const;
      /// get symbol for name alternatives
      const SxParserGet & operator|| (const SxParserGet &in) const {
         return (symbol || !in.symbol) ? *this : in;
      }
      /// get int from symbol
      operator int () const;
      /// get double from symbol
      operator double () const;
      /// get SxString from symbol
      operator SxString () const;
      /// get integer array from symbol
      operator SxArray<int> () const;
      /// Direct access to the SxSymbol
      const SxSymbol* operator-> () const;
      /// Flag conversion
      bool toBool () const
      {
         return symbol ? symbol->toAttribute () : false;
      }
      /// Constructor
      SxParserGet (const char *name_,
                   SxSimpleParser &parser,
                   const char* sourceFile,
                   int sourceLine)
         : theParser(parser), name(name_)
      {
         theParser.where (sourceFile, sourceLine);
         symbol = theParser.current()->contains (name)
                ? theParser.current()->get (name) : NULL;
      }
};

/// \brief Conditional overwrite value from symbol tree
/// (only if symbol is available in tree)
inline void SX_EXPORT_IO operator<< (bool& var, const SxParserGet &me)  {
   if (me.symbol) var = me.symbol->toBool ();
}

/// \brief Conditional overwrite value from symbol tree
/// (only if symbol is available in tree)
inline void SX_EXPORT_IO operator<< (int& var, const SxParserGet &me)  {
   if (me.symbol) var = me;
}

/// \brief Conditional overwrite value from symbol tree
/// (only if symbol is available in tree)
inline void SX_EXPORT_IO operator<< (double& var, const SxParserGet &me)  {
   if (me.symbol) var = me;
}

/// \brief Conditional overwrite value from symbol tree
/// (only if symbol is available in tree)
inline void SX_EXPORT_IO operator<< (SxArray<int>& var, const SxParserGet &me)
{
   if (me.symbol) var = me->toIntList ();
}

/** \brief Utility class for macro-based input file parsing

    \b SxClass = S/PHI/nX simplified parser: symbol groups

    SxParserPush is the interface for nesting groups

    \author C. Freysoldt */
class SX_EXPORT_IO SxParserPush
{
   protected:
      SxSimpleParser &theParser;
   public:
      SxParserPush (const char *groupName,
                    SxSimpleParser &parser,
                    const char* sourceFile,
                    int sourceLine);
      operator bool () const { return theParser.currentGroup; }
      ~SxParserPush ();
};

// three macros to get a unique variable name from the line-number
#define MERGE2(x,y) x ## y
#define MERGE(x,y) MERGE2(x,y)
#define UNIQUENAME(x) MERGE(x,__LINE__)

#define SYMBOLPARSE(x)       if (SxSimpleParser mySimpleParser = SxSimpleParser(x,__FILE__,__LINE__))

#define SYMBOLGROUP(x)       if (SxParserPush UNIQUENAME(level) = SxParserPush(x, mySimpleParser, __FILE__,__LINE__))
#define FOREACH_SYMBOLGROUP(x) if (SxParserPush UNIQUENAME(level) = SxParserPush(x, mySimpleParser, __FILE__,__LINE__)) \
   for ( ; mySimpleParser.current () ; mySimpleParser.next ())
#define SYMBOLGROUP_TABLE mySimpleParser.current ()
#define SYMBOL_COUNT(x) mySimpleParser.current ()->containsGroup (x) \
   ? mySimpleParser.current ()->getGroup (x)->getNItems (x) : 0;
#define HAVE_SYMBOL(x)  (mySimpleParser.current ()->contains (x))

#define SYMBOLGET(x) SxParserGet (x,mySimpleParser,__FILE__,__LINE__)

#endif /* _SX_SIMPLE_PARSER_H_ */
