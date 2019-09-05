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

#ifndef _SX_Demo1_PARSER_H_
#define _SX_Demo1_PARSER_H_

#include <SxDemo1.h>
#include <SxParserBase.h>
#include <SxParserAst.h>

/** \brief ...

    \b SxClass = S/PHI/nX ...

    ....
    */
class SX_EXPORT_DEMO_1 SxDemo1Parser : public SxParserBase,
                                       public SxParserAst
{
   public:

      SxDemo1Parser ();
     ~SxDemo1Parser ();

      // --- defined in SXPARSER_FOOTER
      virtual void initScanner (bool);
      virtual void destroyScanner ();


   protected:
      virtual int parse ();

};

#endif /* _SX_Demo1_PARSER_H_ */
