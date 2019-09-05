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

#include <SxDemo6Ast.h>

std::ostream &operator<< (std::ostream &s, const SxDemo6AstNode &in)
{
   SxString op, type;
   switch (in.type)  {
      case SxDemo6AstNode::Group        : op = "group"; break;
      case SxDemo6AstNode::KeyVal       : op = ":";  break;
      default                           : /* empty: no operator */;
   }
   s << in.id << ":";
   if (op != "")  s << "[" << op << "] ";
   else           s << in.data;
   return s;
}


// ----------------------------------------------------------------------------

SxDemo6AstNode::SxDemo6AstNode () : id(-1)
{
   SX_TRACE ();
}

SxDemo6AstNode::SxDemo6AstNode (ssize_t id_, const Type &type_)
   : id(id_), type(type_)
{
   SX_TRACE ();
}

SxDemo6AstNode::SxDemo6AstNode (ssize_t id_, const Type &type_, const SxVariant &data_)
   : id(id_), type(type_), data(data_)
{
   SX_TRACE ();
}

SxDemo6AstNode::~SxDemo6AstNode ()
{
	// empty
}

bool SxDemo6AstNode::operator== (const SxDemo6AstNode &in) const
{
   SX_TRACE ();
   return in.id == id;
}

// ----------------------------------------------------------------------------

SxDemo6Ast::SxDemo6Ast () : idx(-1)
{
   SX_TRACE ();
   addNode (SxDemo6AstNode::Group, "root");
}

SxDemo6Ast::~SxDemo6Ast ()
{
   SX_TRACE ();
}

ssize_t SxDemo6Ast::addNode (const SxDemo6AstNode::Type &type_)
{
   SX_TRACE ();
   ast.createNode (SxDemo6AstNode (++idx, type_));
   return idx;
}

ssize_t SxDemo6Ast::addNode (const SxDemo6AstNode::Type &type_,
                             const SxVariant &data_)
{
   SX_TRACE ();
   ast.createNode (SxDemo6AstNode(++idx, type_, data_));
   return idx;
}

void SxDemo6Ast::addEdge (ssize_t a, ssize_t b)
{
   SX_TRACE ();
   ast.createEdge (SxDemo6AstNode(a), SxDemo6AstNode(b));
}


void SxDemo6Ast::prependEdge (ssize_t a, ssize_t b)
{
   SX_TRACE ();
   ast.createEdge (SxDemo6AstNode(a), SxDemo6AstNode(b));
}


SxDemo6AstNode &SxDemo6Ast::getNode (ssize_t i)
{
   SX_TRACE ();
   return *(ast.begin (sx::Forward, SxDemo6AstNode(i)));
}

const SxGraph<SxDemo6AstNode> &SxDemo6Ast::getAst () const
{
   return ast;
}

void SxDemo6Ast::printAst () const
{
   SX_TRACE ();
   ast.print (ast.begin(), cout);
}
