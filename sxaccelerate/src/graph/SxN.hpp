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

sx::N::N () : propName(""),
              op(SxGQExprBase::OpType::None),
              captureEnabled(false) { }

sx::N::N (SxString prop_)
{
   SX_CHECK (prop_ != "");
   propName = prop_;
   op = SxGQExprBase::OpType::None;
   captureEnabled = false;
}

SxString sx::N::getPropName () const
{
   return propName;
}

SxString sx::N::getCName () const
{
   return captureName;
}

sx::N &sx::N::c (SxString name)
{
   captureEnabled = true;
   captureName = name;
   return *this;
}

bool sx::N::isCaptured () const
{
   return captureEnabled;
}

sx::N &sx::N::any ()
{
   op = SxGQExprBase::OpType::Any;
   return *this;
}

sx::N::operator SxPtr<SxGQExprBase>() const
{
   SxPtr<SxGQExpr> e = SxPtr<SxGQExpr>::create (propName, val,
                                                captureName,
                                                op,
                                                captureEnabled);
   return e;
}

template<class T>
sx::N &sx::N::eq(const T &propVal)
{
   op  = SxGQExprBase::OpType::Equal;
   val = SxVariant(propVal);
   return *this;
}

template<class T>
N &N::neq(const T &propVal)
{
   op  = SxGQExprBase::OpType::NotEqual;
   val = SxVariant(propVal);
   return *this;
}

