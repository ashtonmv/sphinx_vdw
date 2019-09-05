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

#include <SxDemo1Parser.h>

SxDemo1Parser::SxDemo1Parser ()
   : SxParserBase (),
     SxParserAst ()
{
   SX_TRACE ();
}

SxDemo1Parser::~SxDemo1Parser ()
{
   // empty
}

int SxDemo1Parser_parse (SxDemo1Parser *); // defined by LEX
int SxDemo1Parser::parse ()
{
   return SxDemo1Parser_parse (this);     // enter LEX
}
