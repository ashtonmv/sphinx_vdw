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

#ifndef _SX_CLIST_H_
#define _SX_CLIST_H_

#include <SxError.h>
#include <SxMemConsumer.h>
#include <SxUtil.h>
#include <stdio.h>
#include <iostream>
#include <SxAlg.h>
#include <SxIterator.h>
//#include <SxContainer.h>

template<class T>
class SxCListNode {
   public:
      T elem;
      SxCListNode<T> *prev, *next;

      SxCListNode () : prev(NULL), next(NULL) { }
      SxCListNode (const SxCListNode<T> &in)
         : elem(in.elem), prev(in.prev), next(in.next) { }

      SxCListNode &operator= (const SxCListNode<T> &in) {
         if (this == &in)  return *this;
         elem = in.elem;
         prev = in.prev;
         next = in.next;
         return *this;
      }
};

template<class>
class SxCList;

template<class T,class Node,class Container>
class SxCListIterators
{
   public:
      template<class SValue,class SNode, class SContainer,class IT>
      class State
      {
         SX_ITERATOR_STATE
         // friend to State<const ...>
         template<class CV,class CN,class CC,class CI> friend class State;

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
template<class T>
class SxCList : public SxMemConsumer
{
   public:

      typedef SxCListNode<T> Node;

      // --- Iterators
      typedef SxCList<T> Container;
      typedef typename SxCListIterators<T,Node,Container>::Iterator      Iterator;
      typedef typename SxCListIterators<T,Node,Container>::ConstIterator ConstIterator;


      Node *firstElement;
      Node *lastElement;
      ssize_t size;

      // Constructors
      SxCList ();
      SxCList (const SxCList<T> &);
      SxCList (SxCList<T> &&);

      ~SxCList ();

      // Operations

      // -- Assignment
      SxCList<T> &operator=  (const SxCList<T> &);
      SxCList<T> &operator=  (SxCList<T> &&);

      // -- Equality
      bool       operator== (const SxCList<T> &) const;

      bool       operator!= (const SxCList<T> &) const;

      // -- Insertion
      inline SxCList<T> &operator<< (const T &);

      // move insertion
      inline SxCList<T> &operator<< (T &&);


      // Methods
      // -- Insertion
      SxCList<T> &append (const T &);
      SxCList<T> &append (T &&);

      // -- Search
      ssize_t findPos (const T &) const;
      bool contains (const T &) const;

      // --Size
      inline ssize_t getSize () const { return size; }

      // -- Deletion
      inline void removeItem (Node *ptr);
      inline void removeElement (const T &);
      inline void removeAll ();

      // Accessfunctions
      inline T &first ();
      inline const T &first () const;
      inline T &last ();
      inline const T &last () const;

      // Iteratorfunctions
      inline typename SxCList<T>::Iterator begin ();
      inline typename SxCList<T>::ConstIterator begin () const;
      inline typename SxCList<T>::Iterator begin
                      (ssize_t idx, sx::Direction dir = sx::Forward);
      inline typename SxCList<T>::ConstIterator begin
                (ssize_t idx, sx::Direction dir = sx::Forward) const;
      inline typename SxCList<T>::Iterator end ();
      inline typename SxCList<T>::ConstIterator end () const;
      inline typename SxCList<T>::Iterator fromLast ();
      inline typename SxCList<T>::ConstIterator fromLast () const;
      inline typename SxCList<T>::Iterator toFirst ();
      inline typename SxCList<T>::ConstIterator toFirst () const;

};

//------------------------------------------------------------------------------
// SxCList-class
//------------------------------------------------------------------------------
template<class T>
SxCList<T>::SxCList () : SxMemConsumer ()
{
   size        = 0;
   lastElement = firstElement = NULL;
}

template<class T>
SxCList<T>::SxCList (const SxCList<T> &in) : SxMemConsumer (in)
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

template<class T>
SxCList<T>::SxCList (SxCList<T> &&in)
   : SxMemConsumer (in)
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


template<class T>
SxCList<T>::~SxCList ()
{
   removeAll ();
}

// Operations
// -- Assignment
template<class T>
SxCList<T> &SxCList<T>::operator= (const SxCList<T> &in)
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
template<class T>
SxCList<T> &SxCList<T>::operator= (SxCList<T> &&in)
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


// -- Equality
template<class T>
bool SxCList<T>::operator== (const SxCList<T> &rhs) const
{
   if (getSize () != rhs.getSize ())  {
      return false;
   } else  {

      Node *inPtr = rhs.firstElement;
      Node *ptr   = firstElement;
      while (ptr && inPtr) {
         if (*ptr != *inPtr)  return false;
         ptr   = ptr->next;
         inPtr = inPtr->next;
      }
      return true;
   }
}

template<class T>
bool SxCList<T>::operator!= (const SxCList<T> &rhs) const
{
   return !operator== (rhs);
}

// -- Insertion
template<class T>
SxCList<T> &SxCList<T>::operator<< (const T &t)
{
   return append (t);
}

template<class T>
SxCList<T> &SxCList<T>::operator<< (T &&t)
{
   return append (std::move(t));
}


// Methods
// -- Insertion
template<class T>
SxCList<T> &
SxCList<T>::append (const T &e)
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
   return *this;
}

template<class T>
SxCList<T> &
SxCList<T>::append (T &&e)
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
   return *this;
}

// -- Search
template<class T>
ssize_t SxCList<T>::findPos (const T &e) const
{
   ssize_t i=0;
   Node *ptr = firstElement;
   while ( ptr )  {
      if (ptr->elem == e)  return i;
      ptr = ptr->next;
      i++;
   }

   return -1;
}

template<class T>
bool SxCList<T>::contains (const T &in) const
{
   return findPos (in) == -1 ? false : true;
}

// -- Deletion

template<class T>
void SxCList<T>::removeItem (Node *ptr)
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

template<class T>
void SxCList<T>::removeElement  (const T &e)
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

template<class T>
void SxCList<T>::removeAll ()
{
   if (size > 0)  {
      Node *ptr = lastElement;
      Node *tmp;
      // --- delete entire Node
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

template<class T>
typename SxCList<T>::Iterator SxCList<T>::begin ()
{
   return typename SxCList<T>::Iterator (sx::Forward, this, firstElement);
}

template<class T>
typename SxCList<T>::ConstIterator SxCList<T>::begin () const
{
   return typename SxCList<T>::ConstIterator (sx::Forward, this, firstElement);
}

template<class T>
typename SxCList<T>::Iterator
SxCList<T>::begin (ssize_t idx, sx::Direction dir )
{
   SX_CHECK (idx >= 0 && idx < size, idx, size);
   ssize_t i = 0;
   Node *ptr = firstElement;
   while (ptr) {
      if (i == idx) {
         return typename SxCList<T>::Iterator (dir, this, ptr);
      }
      ptr = ptr->next;
      i += 1;
   }
   return end ();
}

template<class T>
typename SxCList<T>::ConstIterator
SxCList<T>::begin (ssize_t idx, sx::Direction dir ) const
{
   SX_CHECK (idx >= 0 && idx < size, idx, size);
   ssize_t i = 0;
   Node *ptr = firstElement;
   while (ptr) {
      if (i == idx) {
         return typename SxCList<T>::ConstIterator (dir, this, ptr);
      }
      ptr = ptr->next;
      i += 1;
   }
   return end ();
}

template<class T>
typename SxCList<T>::Iterator SxCList<T>::end ()
{
   return typename SxCList<T>::Iterator ();
}

template<class T>
typename SxCList<T>::ConstIterator SxCList<T>::end () const
{
   return typename SxCList<T>::ConstIterator ();
}

template<class T>
typename SxCList<T>::Iterator SxCList<T>::fromLast ()
{
   if ( lastElement )
      return typename SxCList<T>::Iterator (sx::Backward, this,
                                            lastElement);
   else
      return typename SxCList<T>::Iterator ();
}

template<class T>
typename SxCList<T>::ConstIterator SxCList<T>::fromLast () const
{
   if ( lastElement )
      return typename SxCList<T>::ConstIterator (sx::Backward, this,
                                                lastElement);
   else
      return typename SxCList<T>::ConstIterator ();
}

template<class T>
typename SxCList<T>::Iterator SxCList<T>::toFirst ()
{
   return typename SxCList<T>::Iterator ();
}


template<class T>
typename SxCList<T>::ConstIterator SxCList<T>::toFirst () const
{
   return typename SxCList<T>::ConstIterator ();
}


#endif /* _SX_LIST_H_ */
