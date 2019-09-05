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

#ifndef _SX_LIST_H_
#define _SX_LIST_H_

#include <SxError.h>
#include <SxMemConsumer.h>
#include <SxSelection.h>
#include <SxContainer.h>
#include <SxIterator.h>
#include <SxUtil.h>
#include <SxIterator.h>
#include <SxRange.h>
#include <stdio.h>
#include <iostream>
#include <SxAlg.h>
#include <SxPtr.h>

template<class T>
class SxListNode {
   public:
      T elem;
      SxListNode<T> *prev, *next;

      SxListNode () : prev(NULL), next(NULL) { }
      SxListNode (const SxListNode<T> &in)
         : elem(in.elem), prev(in.prev), next(in.next) { }

      SxListNode &operator= (const SxListNode<T> &in) {
         if (this == &in)  return *this;
         elem = in.elem;
         prev = in.prev;
         next = in.next;
         return *this;
      }
};

template<class T,class Node,class Container>
class SxListIterators
{
   public:
      template<class SValue,class SNode, class SContainer,class IT>
      class State
      {
         SX_ITERATOR_STATE
         // friend to State<const ...>
         //template<class CV,class CN,class CC,class CI> friend class State;
	 friend SxListIterators;

         public:

            State (sx::Direction dir_=sx::Forward,
                   SContainer *container_=NULL, SNode *node_ = NULL)
               : dir(dir_), container(container_), node(node_) { }

            State (const State &in, sx::ItCopyMode cmode_ = sx::CopyAll)
               : dir(in.dir), container(in.container), node(in.node)
            {
               SX_UNUSED (cmode_);
            }

            State (State &&in, sx::ItCopyMode cmode_ = sx::CopyAll)
               : dir(in.dir), container(in.container), node(in.node)
            {
               SX_UNUSED (cmode_);
               in.dir = sx::Undefined; in.container = NULL; in.node = NULL;
            }

            typename Container::SelIdx getSelIdx ()  {
               SX_CHECK (node);
               return node;
            }
            // non-const to const cast
            template<class CV,class CN, class CC,class CI>
            State (const State<CV,CN,CC,CI> &in)
                  : dir(in.dir), container(in.container), node(in.node) { }

         protected:
            SContainer *container;
            SNode *node;
            // --- SX_ITERATOR / SX_CONST_ITERATOR callbacks
            void copy (const IT &in_) {
               container = in_.container; node = in_.node; dir = in_.dir;
            }
            void move (IT &&in_) {
               dir = in_.dir; container = in_.container; node = in_.node;
               in_.dir = sx::Undefined; in_.container = NULL; in_.node = NULL;
            }
            IT insertElem (ssize_t newPos, const SValue &elem) {
               SX_CHECK (container);
               return container->insert (newPos, elem);
            }
            IT prependElem (const SValue &elem) {
               SX_CHECK (container);
               return container->prepend (elem);
            }
            IT appendElem (const SValue &elem) {
               SX_CHECK (container);
               return container->append (elem);
            }
            void next () {
               if (node) {
                  node = node->next;
               }
            }
            void prev () {
               if (node) {
                  node = node->prev;
               }
            }
            bool valid () const { return node != NULL; }
            SValue &getRef () { SX_CHECK (node); return  node->elem; }
            SValue *getPtr () { SX_CHECK (node); return &node->elem; }
            bool equal (const IT &in) const { return node == in.node; }

      };

      class Iterator : public State<T,Node,Container,Iterator>
      {
         SX_ITERATOR_NO_LAMBDAS(T,Container,Iterator)
         public:
            Iterator (sx::Direction dir_=sx::Forward, Container *c=NULL,
                      Node *n=NULL)
               : State<T,Node,Container,Iterator> (dir_, c, n) { }

            Iterator (const Iterator &in, sx::ItCopyMode cmode_ = sx::CopyAll)
               : State<T,Node,Container,Iterator> (in, cmode_) { }

            Iterator (Iterator &&in, sx::ItCopyMode cmode_ = sx::CopyAll)
               : State<T,Node,Container,Iterator> (std::move(in), cmode_) { }

      };

      class ConstIterator
         : public State<const T,const Node,const Container,ConstIterator>
      {
         SX_CONST_ITERATOR_NO_LAMBDAS(T,Container,ConstIterator)
         public:
            ConstIterator (sx::Direction dir_=sx::Forward,
                           const Container *c=NULL, const Node *n=NULL)
               : State<const T,const Node,const Container,ConstIterator>
                 (dir_, c, n) { }

            ConstIterator (const ConstIterator &in,
                           sx::ItCopyMode cmode_ = sx::CopyAll)
               : State<const T,const Node,const Container,ConstIterator>
                 (in, cmode_) { }

            ConstIterator (ConstIterator &&in,
                           sx::ItCopyMode cmode_ = sx::CopyAll)
               : State<const T,const Node,const Container,ConstIterator>
                 (std::move(in), cmode_) { }

            // non-const to const cast
            ConstIterator (const Iterator &in)
               : State<const T,const Node,const Container,ConstIterator> (in)
            { }
      };
};



/** This class has been implemented in case the standard template library
    (STL) is not available on some platforms.
    It just provides a doubly linked list.
    @ingroup Tools
    @author Sixten Boeck
    */
template<class T,template<class,class,class> class ItPair=SxListIterators>
class SxList : public SxMemConsumer,
               public SxThis<SxList<T,ItPair> >
{
   public:

      typedef SxListNode<T> Node;

      // --- Iterators
      typedef SxList<T,ItPair> Container;
      typedef typename ItPair<T,Node,Container>::Iterator      Iterator;
      typedef typename ItPair<T,Node,Container>::ConstIterator ConstIterator;

      // --- SxSelection:
      typedef T               TElem;
      typedef Node            *SelIdx;
      typedef SxCList<Node*>  SelStorage;
      typedef Container       SelContainer;

      SX_CONTAINER (Container)

      Node *firstElement;
      Node *lastElement;
      ssize_t size;

      // Constructors
      SxList ();
      SxList (const SxList<T,ItPair> &);
      SxList (const SxCList<T> &);
      SxList (SxList<T,ItPair> &&);
      SxList (const std::initializer_list<T> &);

      // conversion from sx::find
      SxList (const SxList<Iterator> &);
      SxList (const SxList<ConstIterator> &);

      template<long Begin,long End,long Step=1>
         SxList (const SxRange<Begin,End,Step> &range);
      ~SxList ();

      // Operations
      // -- Randomaccess
      inline T &operator() (ssize_t);
      inline const T &operator() (ssize_t) const;

      // -- Assignment
      SxList<T,ItPair> &operator=  (const SxList<T,ItPair> &);
      SxList<T,ItPair> &operator=  (SxList<T,ItPair> &&);
      SxList<T,ItPair> &operator=  (const std::initializer_list<T> &);

      // -- Equality
      bool       operator== (const SxList<T,ItPair> &) const;

      bool       operator!= (const SxList<T,ItPair> &) const;

      // -- Insertion
      inline SxList<T,ItPair> &operator<< (const T &);

      // move insertion
      inline SxList<T,ItPair> &operator<< (T &&);

      inline SxList<T,ItPair> &operator<< (const SxList<T,ItPair> &);
      inline SxList<T,ItPair> &operator<< (SxList<T,ItPair> &&);

      // initializer_list append
      inline SxList<T,ItPair> &operator<< (const std::initializer_list<T> &);

      template<long Begin,long End,long Step=1>
         SxList<T,ItPair> &operator<< (const SxRange<Begin,End,Step> &range);


      // Methods
      // -- Insertion
      Iterator append (const T &);
      Iterator append (T &&);
      Iterator append (const SxList<T,ItPair> &);
      Iterator append (SxList<T,ItPair> &&);
      Iterator append (const std::initializer_list<T> &);
      Iterator append (const Iterator &it, const T &);
      Iterator append (const Iterator &it, T &&);
      Iterator insert (ssize_t newPos, const T &);
      Iterator insert (ssize_t newPos, T &&);
      Iterator prepend (const T &);
      Iterator prepend (T &&);
      Iterator prepend (const Iterator &it, const T &);
      Iterator prepend (const Iterator &it, T &&);

      // -- Search
      ssize_t findPos (const T &) const;
      bool contains (const T &) const;

      // --Size
      inline ssize_t getSize () const { return size; }
      void resize (ssize_t);

      // -- Deletion
      inline void remove (ssize_t);
      inline void removeItem (Node *ptr);
      inline void removeItem (Iterator *it);
      inline void removeElement (const T &);
      inline void removeFirst ();
      inline void removeLast ();
      inline void removeAll ();

      // Accessfunctions
      inline T &first ();
      inline const T &first () const;
      inline T &last ();
      inline const T &last () const;

      inline Iterator      getIterator (const SelIdx &);
      inline ConstIterator getIterator (const SelIdx &) const;
      inline ConstIterator getConstIterator (const Iterator &) const;
      SxPtr<SxList<T,ItPair> > getContainer () const;

      // Iteratorfunctions
      inline typename SxList<T,ItPair>::Iterator begin ();
      inline typename SxList<T,ItPair>::ConstIterator begin () const;
      inline typename SxList<T,ItPair>::Iterator begin
                      (ssize_t idx, sx::Direction dir = sx::Forward);
      inline typename SxList<T,ItPair>::ConstIterator begin
                (ssize_t idx, sx::Direction dir = sx::Forward) const;
      inline typename SxList<T,ItPair>::Iterator end ();
      inline typename SxList<T,ItPair>::ConstIterator end () const;
      inline typename SxList<T,ItPair>::Iterator fromLast ();
      inline typename SxList<T,ItPair>::ConstIterator fromLast () const;
      inline typename SxList<T,ItPair>::Iterator toFirst ();
      inline typename SxList<T,ItPair>::ConstIterator toFirst () const;

};

//------------------------------------------------------------------------------
// SxList-class
//------------------------------------------------------------------------------
template<class T,template<class,class,class> class ItPair>
SxList<T,ItPair>::SxList ()
   : SxMemConsumer (),
     SxThis<SxList<T,ItPair> > ()
{
   size        = 0;
   lastElement = firstElement = NULL;
}

template<class T,template<class,class,class> class ItPair>
SxList<T,ItPair>::SxList (const SxList<T,ItPair> &in)
   : SxMemConsumer (in),
     SxThis<SxList<T,ItPair> > ()
{
   SX_CHECK (in.size >= 0, in.size);

   firstElement = lastElement = NULL;
   size = in.size;
   if (in.size > 0)  {
      Node *srcPtr   = in.firstElement;
      Node *prevDest = NULL;
      firstElement   = new Node;
      Node *destPtr  = firstElement;
      while ( srcPtr )  {
         destPtr->elem = srcPtr->elem;
         destPtr->prev = prevDest;
         destPtr->next = NULL;
         prevDest      = destPtr;

         srcPtr = srcPtr->next;
         if (srcPtr)  destPtr = destPtr->next = new Node;
         lastElement = destPtr;
      }
   }
   TRACK_MALLOC (*this, 1);
}

template<class T,template<class,class,class> class ItPair>
SxList<T,ItPair>::SxList (SxList<T,ItPair> &&in)
   : SxMemConsumer (in),
     SxThis<SxList<T,ItPair> > ()
{
   SX_CHECK (in.size >= 0, in.size);

   firstElement = lastElement = NULL;
   size = in.size;
   if (in.size > 0)  {
      firstElement    = in.firstElement;
      lastElement     = in.lastElement;

      in.firstElement = NULL;
      in.lastElement  = NULL;
      in.size         = 0;
   }
   TRACK_MALLOC (*this, 1);
}

template<class T,template<class,class,class> class ItPair>
SxList<T,ItPair>::SxList (const std::initializer_list<T> &in)
   : SxThis<SxList<T,ItPair> > ()
{
   SX_CHECK (in.size() >= 0, in.size());

   firstElement = lastElement = NULL;
   size = in.size();
   if (size > 0)  {
      ssize_t count = 0;
      Node *prevDest = NULL;
      firstElement   = new Node;
      Node *destPtr  = firstElement;
      for ( auto &elem : in )  {
         destPtr->elem = elem;
         destPtr->prev = prevDest;
         destPtr->next = NULL;
         prevDest      = destPtr;

         ++count;
         if (count < size)  destPtr = destPtr->next = new Node;
         lastElement = destPtr;
      }
   }
   TRACK_MALLOC (*this, 1);
}

template<class T,template<class,class,class> class ItPair>
SxList<T,ItPair>::SxList (const SxCList<T> &in)
   : SxThis<SxList<T,ItPair> > ()
{
   SX_CHECK (in.size >= 0, in.size);

   firstElement = lastElement = NULL;
   size = in.size;
   if (in.size > 0)  {
      typename SxCList<T>::Node *srcPtr = in.firstElement;
      Node *prevDest = NULL;
      firstElement   = new Node;
      Node *destPtr  = firstElement;
      while ( srcPtr )  {
         destPtr->elem = srcPtr->elem;
         destPtr->prev = prevDest;
         destPtr->next = NULL;
         prevDest      = destPtr;

         srcPtr = srcPtr->next;
         if (srcPtr)  destPtr = destPtr->next = new Node;
         lastElement = destPtr;
      }
   }
   TRACK_MALLOC (*this, 1);
}

template<class T,template<class,class,class> class ItPair>
SxList<T,ItPair>::SxList (const SxList<Iterator> &in)
   : SxThis<SxList<T,ItPair> > ()
{
   firstElement = lastElement = NULL;
   size = 0;
   typename SxList<Iterator>::ConstIterator it;
   for (it = in.begin(); it != in.end(); ++it)  append (**it);
}

template<class T,template<class,class,class> class ItPair>
SxList<T,ItPair>::SxList (const SxList<ConstIterator> &in)
   : SxThis<SxList<T,ItPair> > ()
{
   firstElement = lastElement = NULL;
   size = 0;
   typename SxList<ConstIterator>::ConstIterator it;
   for (it = in.begin(); it != in.end(); ++it) append (**it);
}

template<class T,template<class,class,class> class ItPair>
template<long Begin,long End,long Step>
SxList<T,ItPair>::SxList (const SxRange<Begin,End,Step> &range)
   : SxMemConsumer (),
     SxThis<SxList<T,ItPair> > ()
{
   size        = 0;
   lastElement = firstElement = NULL;
   SxList<T,ItPair>::operator<< (range);
}


template<class T,template<class,class,class> class ItPair>
SxList<T,ItPair>::~SxList ()
{
   removeAll ();
}

// Operations
// -- Assignment
template<class T,template<class,class,class> class ItPair>
SxList<T,ItPair> &SxList<T,ItPair>::operator= (const SxList<T,ItPair> &in)
{
   SX_CHECK (in.size >= 0, in.size);
   if ( this == &in ) return *this;

   removeAll ();

   firstElement = lastElement = NULL;
   size = in.size;
   if (in.size > 0)  {
      Node *srcPtr   = in.firstElement;
      Node *prevDest = NULL;
      firstElement   = new Node;
      lastElement    = firstElement;
      Node *destPtr  = firstElement;
      while ( srcPtr )  {
         destPtr->elem = srcPtr->elem;
         destPtr->prev = prevDest;
         destPtr->next = NULL;
         prevDest      = destPtr;

         srcPtr        = srcPtr->next;
         if (srcPtr)  destPtr = destPtr->next = new Node;
         lastElement = destPtr;
      }
   }

   TRACK_MALLOC (*this, 1);
   return *this;
}

// -- Move assignment
template<class T,template<class,class,class> class ItPair>
SxList<T,ItPair> &SxList<T,ItPair>::operator= (SxList<T,ItPair> &&in)
{
   SX_CHECK (in.size >= 0, in.size);
   if ( this == &in ) return *this;

   removeAll ();

   firstElement = lastElement = NULL;
   size = in.size;
   if (in.size > 0)  {
      firstElement    = in.firstElement;
      lastElement     = in.lastElement;

      in.firstElement = NULL;
      in.lastElement  = NULL;
      in.size         = 0;
   }

   TRACK_MALLOC (*this, 1);
   return *this;
}

template<class T,template<class,class,class> class ItPair>
SxList<T,ItPair> &SxList<T,ItPair>::operator= (const std::initializer_list<T> &in)
{
   SX_CHECK (in.size () >= 0, in.size ());

   removeAll ();

   firstElement = lastElement = NULL;
   size = in.size();
   if (size > 0)  {
      ssize_t count = 0;
      Node *prevDest = NULL;
      firstElement   = new Node;
      Node *destPtr  = firstElement;
      for ( auto &elem : in )  {
         destPtr->elem = elem;
         destPtr->prev = prevDest;
         destPtr->next = NULL;
         prevDest      = destPtr;

         ++count;
         if (count < size)  destPtr = destPtr->next = new Node;
         lastElement = destPtr;
      }
   }

   TRACK_MALLOC (*this, 1);
   return *this;
}


// -- Equality
template<class T,template<class,class,class> class ItPair>
bool SxList<T,ItPair>::operator== (const SxList<T,ItPair> &rhs) const
{
   if (getSize () != rhs.getSize ())  {
      return false;
   } else  {
      ConstIterator it;
      ConstIterator itRhs;
      for (itRhs = rhs.begin (), it = begin (); it != end (); ++itRhs, ++it)  {
         if ( (*it) != (*itRhs))  return false;
      }
      return true;
   }
}

template<class T,template<class,class,class> class ItPair>
bool SxList<T,ItPair>::operator!= (const SxList<T,ItPair> &rhs) const
{
   return !operator== (rhs);
}

// -- Randomaccess
template<class T,template<class,class,class> class ItPair>
T &SxList<T,ItPair>::operator() (ssize_t i)
{
   SX_CHECK (i >= 0 && i < size, i, size);
   Node *ptr;
   ssize_t j = i;
// if ( i+1 <= size /2 )  {
   if ( i+1 <= size>>1 )  {   // forward
      ptr = firstElement;
      for (j=0; j<i; j++)  ptr = ptr->next;
   }  else  {             // reverse
      ptr = lastElement;
      for (j=size; j>i+1; j--)  ptr = ptr->prev;
   }
   return ptr->elem;
}

template<class T,template<class,class,class> class ItPair>
const T &SxList<T,ItPair>::operator() (ssize_t i) const
{
   SX_CHECK (i >= 0 && i < size, i, size);
   Node *ptr;
   ssize_t j = i;
// if ( i+1 <= size /2 )  {
   if ( i+1 <= size>>1 )  {   // forward
      ptr = firstElement;
      for (j=0; j<i; j++)  ptr = ptr->next;
   }  else  {             // reverse
      ptr = lastElement;
      for (j=size; j>i+1; j--)  ptr = ptr->prev;
   }
   return ptr->elem;
}

// -- Insertion
template<class T,template<class,class,class> class ItPair>
SxList<T,ItPair> &SxList<T,ItPair>::operator<< (const T &t)
{
   append (t);
   return *this;
}

template<class T,template<class,class,class> class ItPair>
SxList<T,ItPair> &SxList<T,ItPair>::operator<< (T &&t)
{
   append (std::move(t));
   return *this;
}


template<class T,template<class,class,class> class ItPair>
SxList<T,ItPair> &SxList<T,ItPair>::operator<< (const SxList<T,ItPair> &l)
{
   append (l);
   return *this;
}

template<class T,template<class,class,class> class ItPair>
SxList<T,ItPair> &SxList<T,ItPair>::operator<< (SxList<T,ItPair> &&l)
{
   append (std::move(l));
   return *this;
}

template<class T,template<class,class,class> class ItPair>
SxList<T,ItPair> &SxList<T,ItPair>::operator<< (const std::initializer_list<T> &l)
{
   append (l);
   return *this;
}

template<class T,template<class,class,class> class ItPair>
template<long Begin,long End,long Step>
SxList<T,ItPair> &SxList<T,ItPair>::operator<< (const SxRange<Begin,End,Step> &range)
{
   typedef typename SxRange<Begin,End,Step>::ConstIterator ConstIt;
   range.foreach ([this](ConstIt it) { this->append ((T)(*it)); });
   return *this;
}

// Methods
// -- Insertion
template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::Iterator
SxList<T,ItPair>::append (const T &e)
{
   Node *entry = new Node;
   entry->elem = e;

   entry->prev = lastElement;
   entry->next = NULL;

   if ( lastElement )  lastElement->next = entry;
   lastElement = entry;

   if ( !firstElement ) firstElement = entry;
   size++;
   TRACK_MALLOC (*this, 1);
   return typename SxList<T,ItPair>::Iterator(sx::Forward, this, lastElement);
}

template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::Iterator
SxList<T,ItPair>::append (T &&e)
{
   Node *entry = new Node;
   entry->elem = std::move(e);

   entry->prev = lastElement;
   entry->next = NULL;

   if ( lastElement )  lastElement->next = entry;
   lastElement = entry;

   if ( !firstElement ) firstElement = entry;
   size++;
   TRACK_MALLOC (*this, 1);
   return typename SxList<T,ItPair>::Iterator(sx::Forward, this, lastElement);
}

template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::Iterator
SxList<T,ItPair>::append (const SxList<T,ItPair> &in)
{
   typename SxList<T,ItPair>::ConstIterator it;
   for (it = in.begin(); it != in.end(); it++)   append (*it);

   return typename SxList<T,ItPair>::Iterator(sx::Forward, this, lastElement);
}

template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::Iterator
SxList<T,ItPair>::append (SxList<T,ItPair> &&in)
{
   if (in.firstElement) {
      in.firstElement->prev = lastElement;
      if (lastElement) lastElement->next = in.firstElement;
      lastElement = in.lastElement;
   }

   if (!firstElement) firstElement = in.firstElement;
   size+= in.size;

   in.firstElement = NULL;
   in.lastElement  = NULL;
   in.size         = 0;

   return typename SxList<T,ItPair>::Iterator(sx::Forward, this, lastElement);
}

template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::Iterator
SxList<T,ItPair>::append (const std::initializer_list<T> &in)
{
   for (auto &elem : in)   append (elem);

   return typename SxList<T,ItPair>::Iterator(sx::Forward, this, lastElement);
}

template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::Iterator
SxList<T,ItPair>::append (const Iterator &it_, const T &in)
{
   if (getSize() == 0)  return append (in);
   if (*it_ == first()) return insert (1, in);
   if (*it_ == last())  return append (in);

   Iterator it = it_;
   it++;

   Node *ptr = it.ptr;

   Node *entry     = new Node;
   entry->elem     = in;
   entry->prev     = ptr->prev;
   entry->next     = ptr;
   ptr->prev->next = entry;
   ptr->prev       = entry;
   size++;

   return typename SxList<T,ItPair>::Iterator(sx::Forward, this, lastElement);
}

template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::Iterator
SxList<T,ItPair>::append (const Iterator &it_, T &&in)
{
   if (getSize() == 0)  return append (std::move(in));
   if (*it_ == first()) return insert (1, std::move(in));
   if (*it_ == last())  return append (std::move(in));

   Iterator it = it_;
   it++;

   Node *ptr = it.ptr;

   Node *entry     = new Node;
   entry->elem     = std::move(in);
   entry->prev     = ptr->prev;
   entry->next     = ptr;
   ptr->prev->next = entry;
   ptr->prev       = entry;
   size++;

   return typename SxList<T,ItPair>::Iterator(sx::Forward, this, lastElement);
}


template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::Iterator
SxList<T,ItPair>::insert (ssize_t newPos, const T &in)
{
   SX_CHECK (newPos >= 0 && newPos <= getSize(), newPos, getSize());
   if (newPos == 0)  {
      return prepend (in);
   }
   if (newPos == getSize())  {
      return append (in);
   }

   Node *ptr = firstElement;
   for (ssize_t i=0; i < newPos; i++)  {
      ptr = ptr->next;
   }

   Node *entry     = new Node;
   entry->elem     = in;
   entry->prev     = ptr->prev;
   entry->next     = ptr;
   ptr->prev->next = entry;
   ptr->prev       = entry;
   size++;

   return typename SxList<T,ItPair>::Iterator(sx::Forward, this, entry);
}

template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::Iterator
SxList<T,ItPair>::insert (ssize_t newPos, T &&in)
{
   SX_CHECK (newPos >= 0 && newPos <= getSize(), newPos, getSize());
   if (newPos == 0)  {
      return prepend (std::move(in));
   }
   if (newPos == getSize())  {
      return append (std::move(in));
   }

   Node *ptr = firstElement;
   for (ssize_t i=0; i < newPos; i++)  {
      ptr = ptr->next;
   }

   Node *entry     = new Node;
   entry->elem     = std::move(in);
   entry->prev     = ptr->prev;
   entry->next     = ptr;
   ptr->prev->next = entry;
   ptr->prev       = entry;
   size++;

   return typename SxList<T,ItPair>::Iterator(sx::Forward, this, entry);
}

template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::Iterator
SxList<T,ItPair>::prepend (const Iterator &it, const T &in)
{
   if (getSize() == 0)  return append (in);

   if (*it == first())  return prepend (in);

   Node *ptr = it.ptr;

   Node *entry     = new Node;
   entry->elem     = in;
   entry->prev     = ptr->prev;
   entry->next     = ptr;
   ptr->prev->next = entry;
   ptr->prev       = entry;
   size++;

   return typename SxList<T,ItPair>::Iterator(sx::Forward, this, entry);
}

template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::Iterator
SxList<T,ItPair>::prepend (const Iterator &it, T &&in)
{
   if (getSize() == 0)  return append (std::move(in));

   if (*it == first())  return prepend (std::move(in));

   Node *ptr = it.ptr;

   Node *entry     = new Node;
   entry->elem     = std::move(in);
   entry->prev     = ptr->prev;
   entry->next     = ptr;
   ptr->prev->next = entry;
   ptr->prev       = entry;
   size++;

   return typename SxList<T,ItPair>::Iterator(sx::Forward, this, entry);
}


template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::Iterator
SxList<T,ItPair>::prepend (const T &e)
{
   Node *entry    = new Node;
   entry->elem    = e;
   entry->prev    = NULL;
   entry->next    = firstElement;

   if ( firstElement )  firstElement->prev = entry;
   firstElement = entry;

   if ( !lastElement ) lastElement = entry;

   size++;
   TRACK_MALLOC (*this, 1);

   return typename SxList<T,ItPair>::Iterator(sx::Forward, this, firstElement);
}

template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::Iterator
SxList<T,ItPair>::prepend (T &&e)
{
   Node *entry    = new Node;
   entry->elem    = std::move(e);
   entry->prev    = NULL;
   entry->next    = firstElement;

   if ( firstElement )  firstElement->prev = entry;
   firstElement = entry;

   if ( !lastElement ) lastElement = entry;

   size++;
   TRACK_MALLOC (*this, 1);

   return typename SxList<T,ItPair>::Iterator(sx::Forward, this, firstElement);
}

// -- Search
template<class T,template<class,class,class> class ItPair>
ssize_t SxList<T,ItPair>::findPos (const T &e) const
{
   ssize_t i=0;
   Node *ptr = firstElement;
   while ( ptr )  {
      //if ( ptr->elem ==  e )  return i;
//      if ( isEqual (ptr->elem, e) )  return i;
      if (ptr->elem == e)  return i;
      ptr = ptr->next;
      i++;
   }

   return -1;
}

template<class T,template<class,class,class> class ItPair>
bool SxList<T,ItPair>::contains (const T &in) const
{
   return findPos (in) == -1 ? false : true;
}

//template<class T>
//bool SxList<T,ItPair>::isEqual (const T &a, const T &b) const
//{
//   ssize_t nBytes = sizeof (a);
//   char *aPtr = (char *)&a;
//   char *bPtr = (char *)&b;
//   for (ssize_t i=0; i<nBytes; i++)  {
//      if ( *aPtr != *bPtr ) return false;
//      aPtr++;
//      bPtr++;
//   }
//
//   return true;
//}

// -- Size
template<class T,template<class,class,class> class ItPair>
void SxList<T,ItPair>::resize (ssize_t newSize)
{
   SX_CHECK ( newSize >= 0, newSize );

   Node *entry = firstElement;

   if (size == newSize)  return;

   if ( newSize == 0)     {  // clear entire Node

      removeAll ();

   } else if ( newSize < size )  {  // truncate end of Node

      for (ssize_t i=0; i < newSize-1; i++)  {
         entry = entry->next;
      }
      Node *ptr = entry->next, *tmp = NULL;
      while (ptr)  {
         tmp = ptr;
         ptr = ptr->next;
         delete tmp;
      }
      lastElement = entry;
      lastElement->next = NULL;

   } else {                         // add elements to the end

      ssize_t diffSize = newSize - size;
      for (ssize_t i=0; i < diffSize; i++)  {
         entry       = new Node;
         entry->prev = lastElement;
         entry->next = NULL;

         if ( lastElement )  lastElement->next = entry;
         lastElement = entry;

         if ( !firstElement ) firstElement = entry;
      }

   }
   size = newSize;
   TRACK_MALLOC (*this, 1);
}

// -- Deletion
template<class T,template<class,class,class> class ItPair>
void SxList<T,ItPair>::removeItem (Node *ptr)
{
   if (ptr)  {
      Node *prevPtr = ptr->prev;
      Node *nextPtr = ptr->next;

      if (ptr == firstElement) firstElement = nextPtr;
      if (ptr == lastElement)  lastElement = prevPtr;
      if (prevPtr) prevPtr->next = nextPtr;
      if (nextPtr) nextPtr->prev = prevPtr;

      delete ptr;
      size--;
   }
}

template<class T,template<class,class,class> class ItPair>
void SxList<T,ItPair>::removeItem (Iterator *it)
{
   SX_CHECK (it);
   if ((*it).isForward ()) removeElement (*(*it)++);
   else removeElement (*(*it)--);
   //*it = Iterator();
}

template<class T,template<class,class,class> class ItPair>
void SxList<T,ItPair>::remove  (ssize_t idx)
{
   SX_CHECK (idx >= 0 && idx < size, idx, size);

   Node *ptr = firstElement;
   for (ssize_t i=0; i < idx; i++)  {
      ptr = ptr->next;
   }
   removeItem (ptr);
}

template<class T,template<class,class,class> class ItPair>
void SxList<T,ItPair>::removeElement  (const T &e)
{
   Node *ptr = firstElement;
   while ( ptr )  {
      if (ptr->elem == e)  {
         removeItem (ptr);
         return;
      }
      ptr = ptr->next;
   }
}

template<class T,template<class,class,class> class ItPair>
void SxList<T,ItPair>::removeFirst ()
{
   removeItem (firstElement);
}

template<class T,template<class,class,class> class ItPair>
void SxList<T,ItPair>::removeLast ()
{
   removeItem (lastElement);
}

template<class T,template<class,class,class> class ItPair>
void SxList<T,ItPair>::removeAll ()
{
   if (size > 0)  {
      Node *ptr = lastElement;
      Node *tmp;
      // --- delete entire Node
//      while ( ptr )  {
//         tmp = ptr;
//         ptr = ptr->prev;
//         delete tmp;
//      }
      while (size > 1)
      {
         tmp = ptr;
         ptr = ptr->prev;
         delete tmp;
         size--;
      }
      delete ptr;
   }
   size = 0;
   firstElement = lastElement = NULL;
}

// Accessfunctions
template<class T,template<class,class,class> class ItPair>
T &SxList<T,ItPair>::first ()
{
   SX_CHECK (size > 0, size);
   SX_CHECK (firstElement);
   Node *ptr = firstElement;
   return ptr->elem;
}

template<class T,template<class,class,class> class ItPair>
const T &SxList<T,ItPair>::first () const
{
   SX_CHECK (size > 0, size);
   SX_CHECK (firstElement);
   Node *ptr = firstElement;
   return ptr->elem;
}


template<class T,template<class,class,class> class ItPair>
T &SxList<T,ItPair>::last ()
{
   SX_CHECK (size > 0, size);
   SX_CHECK (lastElement);
   Node *ptr = lastElement;
   return ptr->elem;
}

template<class T,template<class,class,class> class ItPair>
const T &SxList<T,ItPair>::last () const
{
   SX_CHECK (size > 0, size);
   SX_CHECK (lastElement);
   Node *ptr = lastElement;
   return ptr->elem;
}

template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::Iterator
SxList<T,ItPair>::getIterator (const typename
                               SxList<T,ItPair>::SelIdx &idx)
{
   SX_CHECK (size > 0, size);
   SX_CHECK (idx);
   return typename SxList<T,ItPair>::Iterator (sx::Forward, this, idx);
}

template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::ConstIterator
SxList<T,ItPair>::getIterator (const typename
                               SxList<T,ItPair>::SelIdx &idx) const
{
   SX_CHECK (size > 0, size);
   SX_CHECK (idx);
   return typename SxList<T,ItPair>::ConstIterator (sx::Forward, this, idx);
}

template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::ConstIterator
SxList<T,ItPair>::getConstIterator (const typename
                                    SxList<T,ItPair>::Iterator &it) const
{
   return it;
}

template<class T,template<class,class,class> class ItPair>
SxPtr<SxList<T,ItPair> > SxList<T,ItPair>::getContainer () const
{
   return this->getThis ();
}

template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::Iterator SxList<T,ItPair>::begin ()
{
   return typename SxList<T,ItPair>::Iterator (sx::Forward, this, firstElement);
}

template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::ConstIterator SxList<T,ItPair>::begin () const
{
   return typename SxList<T,ItPair>::ConstIterator (sx::Forward, this, firstElement);
}

template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::Iterator
SxList<T,ItPair>::begin (ssize_t idx, sx::Direction dir )
{
   SX_CHECK (idx >= 0 && idx < size, idx, size);
   ssize_t i = 0;
   Node *ptr = firstElement;
   while (ptr) {
      if (i == idx) {
         return typename SxList<T,ItPair>::Iterator (dir, this, ptr);
      }
      ptr = ptr->next;
      i += 1;
   }
   return end ();
}

template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::ConstIterator
SxList<T,ItPair>::begin (ssize_t idx, sx::Direction dir ) const
{
   SX_CHECK (idx >= 0 && idx < size, idx, size);
   ssize_t i = 0;
   Node *ptr = firstElement;
   while (ptr) {
      if (i == idx) {
         return typename SxList<T,ItPair>::ConstIterator (dir, this, ptr);
      }
      ptr = ptr->next;
      i += 1;
   }
   return end ();
}

template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::Iterator SxList<T,ItPair>::end ()
{
   return typename SxList<T,ItPair>::Iterator ();
}

template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::ConstIterator SxList<T,ItPair>::end () const
{
   return typename SxList<T,ItPair>::ConstIterator ();
}

template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::Iterator SxList<T,ItPair>::fromLast ()
{
   if ( lastElement )
      return typename SxList<T,ItPair>::Iterator (sx::Backward, this,
                                                  lastElement);
   else
      return typename SxList<T,ItPair>::Iterator ();
}

template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::ConstIterator SxList<T,ItPair>::fromLast () const
{
   if ( lastElement )
      return typename SxList<T,ItPair>::ConstIterator (sx::Backward, this,
                                                       lastElement);
   else
      return typename SxList<T,ItPair>::ConstIterator ();
}

template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::Iterator SxList<T,ItPair>::toFirst ()
{
   return typename SxList<T,ItPair>::Iterator ();
}


template<class T,template<class,class,class> class ItPair>
typename SxList<T,ItPair>::ConstIterator SxList<T,ItPair>::toFirst () const
{
   return typename SxList<T,ItPair>::ConstIterator ();
}

//------------------------------------------------------------------------------
// global functions
//------------------------------------------------------------------------------
template<class T,template<class,class,class> class ItPair>
std::ostream &operator<< (std::ostream &s, const SxList<T,ItPair> &in)
{
   ssize_t i = 0;
   typename SxList<T,ItPair>::ConstIterator it;
   for (it = in.begin (); it != in.end (); it++, i++)
      s << i << ": " << *it << std::endl;
   return s;
}

#ifdef WIN32
template<class T,template<class,class,class> class ItPair>
std::wostream& operator<< (std::wostream &s, const SxList<T,ItPair> &in)
{
   ssize_t i = 0;
   typename SxList<T,ItPair>::ConstIterator it;
   for (it = in.begin (); it != in.end (); it++, i++)
      s << i << ": " << *it << std::endl;
   return s;
}
#endif


template<class T,template<class,class,class> class ItPair>
size_t getNBytes (const SxList<T,ItPair> &in)
{
   typename SxList<T,ItPair>::ConstIterator it;
   size_t nBytes = 0;
   for (it = in.begin(); it != in.end(); ++it)  {
      nBytes += getNBytes (*it)
              + sizeof (typename SxList<T,ItPair>::Node *)   // SxList<T,ItPair>::Node::prev
              + sizeof (typename SxList<T,ItPair>::Node *);  // SxList<T,ItPair>::Node::next
   }
   nBytes += sizeof (in);  // = sizeof (SxList<T,ItPair>::firstElement)
                           // + sizeof (SxList<T,ItPair>::lastElement)
                           // + sizeof (SxList<T,ItPair>::size)
                           // + sizeof (SxMemConsumer::className)
                           // + sizeof (SxMemConsumer::parent)
                           // + sizeof (SxMemConsumer::memCounterId)

   return nBytes;
}

namespace sx
{
   template<class Iterator, class Cond>
   SxList<Iterator> findAll (Iterator it, Iterator end, Cond cond)
   {
      SxList<Iterator> res;
      for (/* empty */; it != end; ++it)  {
         if ( cond(it) )  res << it;
      }
      return res;
   }
}



#endif /* _SX_LIST_H_ */
