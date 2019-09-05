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

#include <SxDemo3Parser.h>

SxDemo3Parser::SxDemo3Parser ()
   : SxParserBase (),
     SxParserAst ()
{
   SX_TRACE ();
}

SxDemo3Parser::~SxDemo3Parser ()
{
   // empty
}

int SxDemo3Parser_parse (SxDemo3Parser *); // defined by LEX
int SxDemo3Parser::parse ()
{
   return SxDemo3Parser_parse (this);     // enter LEX
}
