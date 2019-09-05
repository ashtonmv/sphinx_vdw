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

#ifndef _SX_JSON_PARSER_H_
#define _SX_JSON_PARSER_H_

#include <SxJSON.h>
#include <SxParserBase.h>
#include <SxParserAst.h>

/** \brief ...

    \b SxClass = S/PHI/nX ...

    ....
    */
class SX_EXPORT_JSON SxJSONParser : public SxParserBase,
                                    public SxParserAst
{
   public:

      SxJSONParser ();
     ~SxJSONParser ();

      // --- yacc stack interface
      void push (ssize_t);
      SxGProps *push (const SxVariantType::DataType &);
      SxGProps *pop ();
      SxGProps *peek ();
      SxGProps *getParent ();

      SxGProps *getRoot ();

      // --- defined in SXPARSER_FOOTER
      virtual void initScanner (bool);
      virtual void destroyScanner ();


   protected:

      // --- yacc stack
      SxList<SxGProps *> stack;

      virtual int parse ();

};

#endif /* _SX_JSON_PARSER_H_ */
