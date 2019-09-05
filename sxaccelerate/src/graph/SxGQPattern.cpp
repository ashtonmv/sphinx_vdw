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

#include <SxGQPattern.h>
#include <SxGQExprBase.h>

size_t sxHash (const SxGQPattern &in) {
   return SxHashFunction::hash (in.getId ());
}


SxGQPattern::SxGQPattern () : id(-1), relType(RelationType::OrderedDirect) {}
SxGQPattern::SxGQPattern (ssize_t id_) : id(id_), relType(RelationType::OrderedDirect) {}
SxGQPattern::SxGQPattern (ssize_t id_, const SxPtr<SxGQExprBase> &expr_) :
                          id(id_), expr(expr_), relType(RelationType::OrderedDirect) {}

SxGQPattern::SxGQPattern (const SxGQPattern &in)
{
   id      = in.id;
   expr    = in.expr;
   relType = in.relType;
}

SxGQPattern::~SxGQPattern () {}

SxGQPattern &SxGQPattern::operator= (const SxGQPattern &in)
{
   if (this == &in)  return *this;
   id      = in.id;
   expr    = in.expr;
   relType = in.relType;
   return *this;
}

bool SxGQPattern::operator== (const SxGQPattern &in) const
{
   return id == in.id;
}

void SxGQPattern::setRelType (int rType)
{
   relType = static_cast<RelationType>(rType);
}

SxGQPattern::RelationType SxGQPattern::getRelType () const
{
   return relType;
}

ssize_t SxGQPattern::getId () const
{
   return id;
}

void SxGQPattern::setExpr (const SxPtr<SxGQExprBase> &expr_)
{
   SX_CHECK (expr_.getPtr ());
   expr = expr_;
}

bool SxGQPattern::eval (const SxGraph<SxGProps>::ConstIterator &it,
                      const Selection &sel) const
{
   SX_CHECK (expr.getPtr ());
   return expr->eval (it, sel);
}

bool SxGQPattern::matchAll (const SxGraph<SxGProps>::ConstIterator &it,
                          const SelSet &sels) const
{
   SX_CHECK (expr.getPtr ());
   return expr->matchAll (it, sels);
}

