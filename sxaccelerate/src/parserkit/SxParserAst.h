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

#ifndef _SX_PARSER_AST_H_
#define _SX_PARSER_AST_H_

#include <SxParserKit.h>
#include <SxString.h>
#include <SxGraph.h>
#include <SxGProps.h>


/** \brief ...

    \b SxClass = S/PHI/nX ...

    ....

    \author John Doe, johndoe@example.org */
class SX_EXPORT_PARSER_KIT SxParserAst
{
   public:

      SxParserAst ();
     ~SxParserAst ();

      bool validateKey (const SxString &key_) const;

      SxGProps *addNode (const SxVariantType::DataType &);

      void addEdge (ssize_t, ssize_t);
      void addEdge (const SxGProps *, const SxGProps *);
      void prependEdge (ssize_t, ssize_t);

      SxGProps *getNode (ssize_t);

      const SxPtr<SxGraph<SxGProps> > &getAst () const;

      void printAst () const;

   protected:

      ssize_t idx;
      SxPtr<SxGraph<SxGProps> > ast;

};

#endif /* _SX_PARSER_AST_H_ */
