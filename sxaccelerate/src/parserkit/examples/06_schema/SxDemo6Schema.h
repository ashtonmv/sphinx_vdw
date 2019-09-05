#ifndef _SX_DEMO6_SCHEMA_H_
#define _SX_DEMO6_SCHEMA_H_

#include <SxDemo6Ast.h>
#include <SxGraph.h>

typedef typename SxDemo6AstNode::Type Type;

class SX_EXPORT_DEMO_6 SxDemo6Schema
{
   public:

      SxDemo6Schema ();
     ~SxDemo6Schema ();
      ssize_t validate (const SxGraph<SxDemo6AstNode> &dataG) const;
};

#endif /* _SX_DEMO6_SCHEMA_H_ */
