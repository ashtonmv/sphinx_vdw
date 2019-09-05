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

#include <SxDemo5Ast.h>

std::ostream &operator<< (std::ostream &s, const SxDemo5AstNode &in)
{
   SxString op, type;
   switch (in.type)  {
      case SxDemo5AstNode::Group        : op = "group"; break;
      case SxDemo5AstNode::KeyVal       : op = ":";  break;
      default                           : /* empty: no operator */;
   }
   s << in.id << ":";
   if (op != "")  s << "[" << op << "] ";
   else           s << in.data;
   return s;
}


// ----------------------------------------------------------------------------

SxDemo5AstNode::SxDemo5AstNode () : id(-1)
{
   SX_TRACE ();
}

SxDemo5AstNode::SxDemo5AstNode (ssize_t id_, const Type &type_)
   : id(id_), type(type_)
{
   SX_TRACE ();
}

SxDemo5AstNode::SxDemo5AstNode (ssize_t id_, const Type &type_, const SxVariant &data_)
   : id(id_), type(type_), data(data_)
{
   SX_TRACE ();
}

SxDemo5AstNode::~SxDemo5AstNode ()
{
	// empty
}

bool SxDemo5AstNode::operator== (const SxDemo5AstNode &in) const
{
   SX_TRACE ();
   return in.id == id;
}

// ----------------------------------------------------------------------------

SxDemo5Ast::SxDemo5Ast () : idx(-1)
{
   SX_TRACE ();
   addNode (SxDemo5AstNode::Group, "root");
}

SxDemo5Ast::~SxDemo5Ast ()
{
   SX_TRACE ();
}

ssize_t SxDemo5Ast::addNode (const SxDemo5AstNode::Type &type_)
{
   SX_TRACE ();
   ast.createNode (SxDemo5AstNode (++idx, type_));
   return idx;
}

ssize_t SxDemo5Ast::addNode (const SxDemo5AstNode::Type &type_,
                             const SxVariant &data_)
{
   SX_TRACE ();
   ast.createNode (SxDemo5AstNode(++idx, type_, data_));
   return idx;
}

void SxDemo5Ast::addEdge (ssize_t a, ssize_t b)
{
   SX_TRACE ();
   ast.createEdge (SxDemo5AstNode(a), SxDemo5AstNode(b));
}


void SxDemo5Ast::prependEdge (ssize_t a, ssize_t b)
{
   SX_TRACE ();
   ast.createEdge (SxDemo5AstNode(a), SxDemo5AstNode(b));
}


SxDemo5AstNode &SxDemo5Ast::getNode (ssize_t i)
{
   SX_TRACE ();
   return *(ast.begin (sx::Forward, SxDemo5AstNode(i)));
}

const SxGraph<SxDemo5AstNode> &SxDemo5Ast::getAst () const
{
   return ast;
}

void SxDemo5Ast::printAst () const
{
   SX_TRACE ();
   ast.print (ast.begin(), cout);
}
