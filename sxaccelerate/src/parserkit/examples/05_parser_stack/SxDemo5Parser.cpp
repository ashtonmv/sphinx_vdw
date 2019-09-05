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

#include <SxDemo5Parser.h>

SxDemo5Parser::SxDemo5Parser ()
   : SxParserBase (),
     SxDemo5Ast ()
{
   SX_TRACE ();
   ssize_t rootId = 0;
   push (rootId);
}

SxDemo5Parser::~SxDemo5Parser ()
{
   // empty
}

void SxDemo5Parser::push (ssize_t id)
{
   SX_TRACE ();
   stack.append (id);
}


SxDemo5AstNode &SxDemo5Parser::pop ()
{
   SX_TRACE ();
   SX_CHECK (errors.getSize() > 0 || stack.getSize() > 0);
   ssize_t id = stack.last ();
   stack.removeLast ();
   return *(ast.begin(sx::Forward, id)); 
}

SxDemo5AstNode &SxDemo5Parser::getCurrent ()
{
   SX_TRACE ();
   return *(ast.begin(sx::Forward, stack.last())); 
}


SxDemo5AstNode &SxDemo5Parser::getParent ()
{
   SX_TRACE ();
   ssize_t i = stack.getSize()-2;
   if (i >= 0)  return *(ast.begin(sx::Forward, stack(i))); 
   else         return *(ast.begin(sx::Forward, 0));  // top level
}


int SxDemo5Parser_parse (SxDemo5Parser *); // defined by LEX
int SxDemo5Parser::parse ()
{
   return SxDemo5Parser_parse (this);     // enter LEX
}
