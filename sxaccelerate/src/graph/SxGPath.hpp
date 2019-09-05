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

SxGPath::SxGPath () { }
SxGPath::SxGPath (const SxPtr<SxGraph<SxGProps> > &gPtr_,
                  const SxGPath::Selection &sel_)
{
   SX_CHECK (gPtr_.getPtr ());
   SX_CHECK (sel_.getPtr ());

   gPtr = gPtr_;
   sel  = sel_;
}

SxGPath::~SxGPath () { }

ssize_t SxGPath::getSize () const
{
   return sel->getSize ();
}

SxGProps &SxGPath::operator() (ssize_t idx)
{
   ssize_t nodeIdx = -1;
   SX_CHECK (idx < sel->getSize (), idx, sel->getSize ());
   nodeIdx = (*sel)( idx );
   auto it = gPtr->begin (nodeIdx);
   SX_CHECK (it.isValid ());
   return *it;
}

template<class Fn>
void SxGPath::foreach (Fn fn)
{
   SX_CHECK (gPtr.getPtr ());
   SxGraph<SxGProps>::Iterator it = gPtr->begin (sel);
   for (; it != gPtr->end (); ++it) {
      fn (it);
   }
}

template<class Fn>
void SxGPath::foreach (Fn fn) const
{
   SX_CHECK (gPtr.getPtr ());
   SxGraph<SxGProps>::ConstIterator it = gPtr->begin (sel);
   for (; it != gPtr->end (); ++it) {
      fn (it);
   }
}
