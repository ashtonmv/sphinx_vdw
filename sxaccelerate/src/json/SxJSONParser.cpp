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

#include <SxJSONParser.h>

SxJSONParser::SxJSONParser ()
   : SxParserBase (),
     SxParserAst ()
{
   SX_TRACE ();
   ssize_t rootId = 0;
   push (rootId);
}

SxJSONParser::~SxJSONParser ()
{
   // empty
}

void SxJSONParser::push (ssize_t id)
{
   SX_TRACE ();
   stack.append (&(*(ast->begin(sx::Forward, id))));
}

SxGProps *SxJSONParser::push (const SxVariantType::DataType &type_)
{
   SX_TRACE ();
   SxGProps *res = addNode (type_);
   stack.append (res);
   return res;
}

SxGProps *SxJSONParser::pop ()
{
   SX_TRACE ();
   SX_CHECK (errors.getSize() > 0 || stack.getSize() > 0);
   SxGProps *res = stack.last ();
   stack.removeLast ();
   return res;
}

SxGProps *SxJSONParser::peek ()
{
   SX_TRACE ();
   return stack.last ();
}


SxGProps *SxJSONParser::getParent ()
{
   SX_TRACE ();
   ssize_t i = stack.getSize()-2;
   if (i >= 0)  return stack(i);
   else         return stack(0);  // top level
}

SxGProps *SxJSONParser::getRoot ()
{
   SX_TRACE ();
   SX_CHECK (stack.getSize ()>0, stack.getSize ());
   return stack.first ();
}

int SxJSONParser_parse (SxJSONParser *); // defined by LEX
int SxJSONParser::parse ()
{
   return SxJSONParser_parse (this);     // enter LEX
}
