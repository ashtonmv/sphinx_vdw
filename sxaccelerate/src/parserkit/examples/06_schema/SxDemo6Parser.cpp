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

#include <SxDemo6Parser.h>

SxDemo6Parser::SxDemo6Parser ()
   : SxParserBase (),
     SxDemo6Ast ()
{
   SX_TRACE ();
   ssize_t rootId = 0;
   push (rootId);
}

SxDemo6Parser::~SxDemo6Parser ()
{
   // empty
}

void SxDemo6Parser::push (ssize_t id)
{
   SX_TRACE ();
   stack.append (id);
}


SxDemo6AstNode &SxDemo6Parser::pop ()
{
   SX_TRACE ();
   SX_CHECK (errors.getSize() > 0 || stack.getSize() > 0);
   ssize_t id = stack.last ();
   stack.removeLast ();
   return *(ast.begin(sx::Forward, id)); 
}

SxDemo6AstNode &SxDemo6Parser::getCurrent ()
{
   SX_TRACE ();
   return *(ast.begin(sx::Forward, stack.last())); 
}


SxDemo6AstNode &SxDemo6Parser::getParent ()
{
   SX_TRACE ();
   ssize_t i = stack.getSize()-2;
   if (i >= 0)  return *(ast.begin(sx::Forward, stack(i))); 
   else         return *(ast.begin(sx::Forward, 0));  // top level
}


int SxDemo6Parser_parse (SxDemo6Parser *); // defined by LEX
int SxDemo6Parser::parse ()
{
   return SxDemo6Parser_parse (this);     // enter LEX
}
