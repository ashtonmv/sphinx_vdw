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

SxGPathList::SxGPathList (const SxPtr<SxGraph<SxGProps> > &gPtr_,
                          const SelSet &sels_)
{
   SX_CHECK (gPtr_.getPtr ());
   SX_CHECK (sels_.getPtr ());

   gPtr     = gPtr_;
   sels     = sels_;
}

SxGPathList::SxGPathList (const SxPtr<SxGraph<SxGProps> > &gPtr_,
                          const Selection &sel_)
{
   SX_CHECK (gPtr_.getPtr ());
   SX_CHECK (sel_.getPtr ());

   gPtr     = gPtr_;
   sels = SelSet::create ();
   sels->append (sel_);
}

SxGPathList::~SxGPathList () { }

ssize_t SxGPathList::getSize () const
{
   return sels->getSize ();
}

ssize_t SxGPathList::getPathSize () const
{
   if (sels->getSize () > 0)
      return sels->first()->getSize ();
   else
      return 0;
}

SxGPath SxGPathList::operator() (ssize_t selIdx)
{
   SX_CHECK (selIdx < sels->getSize (), selIdx, sels->getSize ());
   return SxGPath(gPtr, (*sels)(selIdx));
}

SxGPathList::Iterator SxGPathList::begin ()
{
   return typename SxGPathList::Iterator (this, sels->begin ());
}
SxGPathList::ConstIterator SxGPathList::begin () const
{
   return typename SxGPathList::ConstIterator (this, sels->begin ());
}

SxGPathList::Iterator SxGPathList::begin (ssize_t idx)
{
   SX_CHECK (idx < sels->getSize (), idx, sels->getSize ());
   return typename SxGPathList::Iterator (this, sels->begin (idx));
}

SxGPathList::ConstIterator SxGPathList::begin (ssize_t idx) const
{
   SX_CHECK (idx < sels->getSize (), idx, sels->getSize ());
   return typename SxGPathList::ConstIterator (this, sels->begin (idx));
}

SxGPathList::Iterator SxGPathList::end ()
{
   return typename SxGPathList::Iterator (this, sels->end ());
}

SxGPathList::ConstIterator SxGPathList::end () const
{
   return typename SxGPathList::ConstIterator (this, sels->end ());
}

template<class Fn>
void SxGPathList::foreach (Fn fn)
{
   sx::foreach (begin (), end (), fn);
}
template<class Fn>
void SxGPathList::foreach (Fn fn) const
{
   sx::foreach (begin (), end (), fn);
}


std::ostream &operator<< (std::ostream &s,
                          const SxGPathList::ConstIterator &in)
{
   ssize_t i = 0;
   in->foreach ([&](auto it) {
                  if (i != 0)  s << " - ";
                  s << it->getId ();
                  i++;
               });
   s << "\n";
   return s;
}

std::ostream &operator<< (std::ostream &s,
                          const SxGPathList::Iterator &in)
{
   ssize_t i = 0;
   in->foreach ([&](auto it) {
                  if (i != 0)  s << " - ";
                  s << it->getId ();
                  i++;
               });
   s << "\n";
   return s;
}

