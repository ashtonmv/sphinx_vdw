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

#include <SxDemo4Parser.h>

SxDemo4Parser::SxDemo4Parser ()
   : SxParserBase (),
     SxDemo4Ast ()
{
   SX_TRACE ();
   ssize_t rootId = 0;
   push (rootId);
}

SxDemo4Parser::~SxDemo4Parser ()
{
   // empty
}

void SxDemo4Parser::push (ssize_t id)
{
   SX_TRACE ();
   stack.append (id);
}


SxDemo4AstNode &SxDemo4Parser::pop ()
{
   SX_TRACE ();
   SX_CHECK (errors.getSize() > 0 || stack.getSize() > 0);
   ssize_t id = stack.last ();
   stack.removeLast ();
   return *(ast.begin(sx::Forward, id)); 
}

SxDemo4AstNode &SxDemo4Parser::getCurrent ()
{
   SX_TRACE ();
   return *(ast.begin(sx::Forward, stack.last())); 
}


SxDemo4AstNode &SxDemo4Parser::getParent ()
{
   SX_TRACE ();
   ssize_t i = stack.getSize()-2;
   if (i >= 0)  return *(ast.begin(sx::Forward, stack(i))); 
   else         return *(ast.begin(sx::Forward, 0));  // top level
}


int SxDemo4Parser_parse (SxDemo4Parser *); // defined by LEX
int SxDemo4Parser::parse ()
{
   return SxDemo4Parser_parse (this);     // enter LEX
}
