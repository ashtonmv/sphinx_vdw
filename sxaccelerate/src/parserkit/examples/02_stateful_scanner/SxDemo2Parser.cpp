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

#include <SxDemo2Parser.h>

SxDemo2Parser::SxDemo2Parser ()
   : SxParserBase (),
     SxParserAst ()
{
   SX_TRACE ();
}

SxDemo2Parser::~SxDemo2Parser ()
{
   // empty
}

int SxDemo2Parser_parse (SxDemo2Parser *); // defined by LEX
int SxDemo2Parser::parse ()
{
   return SxDemo2Parser_parse (this);     // enter LEX
}
