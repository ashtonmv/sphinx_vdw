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

template <class T>
SxPtr<SxOperatorBase<T> > SxOperatorBase<T>::getCopy (T *) const
{
   using namespace std;
   cout << "Missing virtual getCopy for current SxOperatorBase." << endl;
   SX_EXIT;
   return  SxPtr<SxOperatorBase<T> > ();
}

template<class T>
T SxOperator<T>::operator*(const T &in) const
{
   SX_CHECK (first);
   if (next)
      return next->operator*( first->operator*(in) );
   else
      return first->operator*(in);
}

template<class T>
void SxOperator<T>::applyInPlace(T &in) const
{
   SX_CHECK (first);
   first->applyInPlace (in);
   if (next)
      next->applyInPlace (in);
}

template <class T>
static inline bool isIdentity (const SxConstPtr<SxOperatorBase<T> > &ptr)
{
   return dynamic_cast<const typename SxOperator<T>::Identity*>(&*ptr);
}

template <class T>
void SxOperator<T>::append (const SxConstPtr<SxOperatorBase<T> > &ap)
{
   SX_CHECK (ap);
   if (!first)  {
      // first in pipe
      SX_CHECK (!next);
      first = ap;
      return;
   }
   // Identity is a no-op
   if (isIdentity(ap)) return;
   if (!next)  {
      if (isIdentity(first))  {
         first = ap; // first was a no-op
      } else {
         // next slot is still free
         next = ap;
      }
      return;
   }

   // --- append after next

   // try: maybe next is SxOperator?
   const SxOperator<T> *nextIsPipe
      = dynamic_cast<const SxOperator<T> *>(next.getPtr ());
   if (nextIsPipe && next.refCounter == 0)  {
      // next is SxOperator with no further references: append to that
      // for this, we cast away constness
      const_cast<SxOperator<T> *>(nextIsPipe)->append (ap);
      return;
   }

   // create new pipe item 
   SxPtr<SxOperator<T> > newItem
     = SxPtr<SxOperator<T> >::create (next);
   *newItem << ap ;
   next = newItem;
}

// ----------------------------------------------------------------------

template <class T>
SxOperator<T> & SxOperator<T>::operator= (const SxOperatorBase<T> &op)
{
   const SxOperator<T> *isPipe 
      = dynamic_cast<const SxOperator<T>*>(&op);
   if (isPipe)  {
      // avoid superfluous indirections and copy first
      first = isPipe->first;
      next = isPipe->next;
   } else {
      first = op.getCopy ((T*)NULL);
      next = SxConstPtr<SxOperatorBase<T> > ();
   }
   return *this;
}

template <class T>
SxOperator<T> & 
SxOperator<T>::operator= (const SxConstPtr<SxOperatorBase<T> > &ap)
{
   SX_CHECK (ap);
   first = ap;
   next = SxConstPtr<SxOperatorBase<T> > ();
   return *this;
}
   
template <class T>
void SxOperator<T>::prepend (const SxConstPtr<SxOperatorBase<T> > &pre)
{
   SX_CHECK (pre);
   if (first && !isIdentity(first))
      next = getCopy ((T*)NULL);
   first = pre;
}

