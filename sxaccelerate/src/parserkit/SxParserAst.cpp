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

#include <SxParserAst.h>


SxParserAst::SxParserAst () : idx(-1)
{
   SX_TRACE ();
   ast = SxPtr<SxGraph<SxGProps> >::create ();
   SxGProps *n = addNode (SxVariantType::Group);
   n->setProperty ("__sx_Key", "root");
}

SxParserAst::~SxParserAst ()
{
   SX_TRACE ();
}

bool SxParserAst::validateKey (const SxString &key_) const
{
   return ( key_.find ("__sx_") != 0 );
}

SxGProps *SxParserAst::addNode (const SxVariantType::DataType &type_)
{
   SX_TRACE ();
   SxGProps n(++idx);
   SxVariant v;
   v.setType (type_);
   n.setProperty ("__sx_Value", v);
   return const_cast<SxGProps *>(&(*(ast->createNode (n))));
}

void SxParserAst::addEdge (ssize_t a, ssize_t b)
{
   SX_TRACE ();
   ast->createEdge (SxGProps(a), SxGProps(b));
}

void SxParserAst::addEdge (const SxGProps *a, const SxGProps *b)
{
   SX_TRACE ();
   ast->createEdge (*a, *b);
}

void SxParserAst::prependEdge (ssize_t a, ssize_t b)
{
   SX_TRACE ();
   ast->createEdge (SxGProps(a), SxGProps(b));
}

SxGProps *SxParserAst::getNode (ssize_t i)
{
   SX_TRACE ();
   return &(*(ast->begin (sx::Forward, SxGProps(i))));
}

const SxPtr<SxGraph<SxGProps> > &SxParserAst::getAst () const
{
   return ast;
}

void SxParserAst::printAst () const
{
   SX_TRACE ();
   ast->print (ast->begin(), cout);
}
