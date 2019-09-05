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
#include <SxSimpleParser.h>


double SxParserGet::operator|| (double defaultValue) const
{
   double res = defaultValue;
   if (symbol) {
      try {
         res = symbol->toReal ();
      } catch (SxException e)  {
         theParser.error ();
         e.print ();
         SX_EXIT;
      }
   }
   return res;
}

/// get string from symbol or default value
SxString SxParserGet::operator|| (const char *defaultValue) const
{
   if (symbol) {
      SxString res;
      try {
         res = symbol->toString ();
      } catch (SxException e)  {
         theParser.error ();
         e.print ();
         SX_EXIT;
      }
      return res;
   }
   return defaultValue;
}

/// get string from symbol or default value
SxString SxParserGet::operator|| (const SxString &defaultValue) const
{
   if (symbol) {
      SxString res;
      try {
         res = symbol->toString ();
      } catch (SxException e)  {
         theParser.error ();
         e.print ();
         SX_EXIT;
      }
      return res;
   }
   return defaultValue;
}

/// get int from symbol or default value
int SxParserGet::operator|| (int defaultValue) const
{
   int res = defaultValue;
   if (symbol) {
      try {
         res = symbol->toInt ();
      } catch (SxException e)  {
         theParser.error ();
         e.print ();
         SX_EXIT;
      }
   }
   return res;
}
/// get int from symbol
SxParserGet::operator int () const
{
   try {
      if (symbol) return symbol->toInt ();
      theParser.current()->get (name); // raise error
   } catch (SxException e)  {
      theParser.error ();
      e.print ();
      SX_EXIT;
   }
   return -1;
}
/// get string from symbol
SxParserGet::operator SxString () const
{
   try {
      if (symbol) return symbol->toString ();
      theParser.current()->get (name); // raise error
   } catch (SxException e)  {
      theParser.error ();
      e.print ();
      SX_EXIT;
   }
   return -1;
}
/// get double from symbol
SxParserGet::operator double () const
{
   try {
      if (symbol) return symbol->toReal ();
      theParser.current()->get (name); // raise error
   } catch (SxException e)  {
      theParser.error ();
      e.print ();
      SX_EXIT;
   }
   return -1.;
}

/// get int array from symbol
SxParserGet::operator SxArray<int> () const
{
   try {
      if (symbol) return symbol->toIntList ();
      theParser.current()->get (name); // raise error
   } catch (SxException e)  {
      theParser.error ();
      e.print ();
      SX_EXIT;
   }
   return SxArray<int> ();
}

/// Direct access to the SxSymbol
const SxSymbol* SxParserGet::operator-> () const
{
   if (symbol) return symbol;
   try {
      return theParser.current()->get (name);
   } catch (SxException e)  {
      theParser.error ();
      e.print ();
      SX_EXIT;
   }
}

SxParserPush::SxParserPush (const char *groupName,
                            SxSimpleParser &parser,
                            const char* sourceFile,
                            int sourceLine)
   : theParser(parser)
{
   SX_CHECK (theParser.currentGroup);
   // set location for error messages
   theParser.where (sourceFile, sourceLine);
   // put current group on stack
   theParser.outerGroups.push (theParser.currentGroup);
   // --- now try to get the requested group
   SxString name = groupName;
   theParser.currentGroup = theParser.currentGroup->containsGroup (name)
                          ? theParser.currentGroup->getGroup (name)
                          : NULL;
   // fallback: the first group may match the table we started with
   if (!theParser.currentGroup && theParser.outerGroups.getSize () == 1
       && theParser.outerGroups.top ()->name == name)
   {
      theParser.currentGroup = theParser.outerGroups.top ();
   }
}

SxParserPush::~SxParserPush ()  {
   theParser.currentGroup = theParser.outerGroups.pop ();
}
