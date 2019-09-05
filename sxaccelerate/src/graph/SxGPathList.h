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

#ifndef _SX_G_PATH_LIST_H_
#define _SX_G_PATH_LIST_H_

#include <iostream>
#include <SxGraph.h>
#include <SxGProps.h>
#include <SxGPath.h>

/** \brief Path List Class

    \b SxGPathList = SPHInX Path List Class

    SxGPathList allows to store a graph query result(selections).
    The individual selection/path is selected by providing it's
    index in () operator, which returns the individual SxGPath
    object. A second () operator can be used to access individual
    node out of the selected path object.

 */
class SxGPathList
{
   public:

      typedef SxPtr<SxUniqueList<ssize_t> > Selection;
      typedef SxPtr<SxList<Selection> >  SelSet;

      typedef SxGPath               T;
      typedef SxGPathList   Container;

      // ------- Iterators -------

      template<class SValue,class Container,class IT>
      class State
      {
         SX_ITERATOR_STATE
         // frient to State<const...>
         template<class CV,class CC,class CIT> friend class State;

         public:
            State () : dir(sx::Forward), container(NULL) { }
            State (Container *c,
                   SxList<Selection>::Iterator it)
                   : dir(sx::Forward), container(c), selIt(it) {
               if (valid ())
                  currObj = SValue(container->gPtr, *selIt);
            }

            State (const State &in, sx::ItCopyMode cmode_ = sx::CopyAll)
               : dir(in.dir), container(in.container), selIt(in.selIt),
                 currObj(in.currObj)
            {
               SX_UNUSED (cmode_);
            }

            State (State &&in, sx::ItCopyMode cmode_ = sx::CopyAll)
               : dir(in.dir), container(in.container), selIt(in.selIt),
                 currObj(in.currObj)
            {
               SX_UNUSED (cmode_);
               in.dir = sx::Undefined; in.container = NULL;
            }
            // non-const to const cast
            template<class CV,class CC,class CIT>
            State (const State<CV, CC,CIT> &in)
               : dir(in.dir),container(in.container), selIt(in.selIt),
                 currObj(in.currObj) { }

         protected:
            Container *container;
            ssize_t idx;
            SxList<Selection>::Iterator selIt;
            SxGPath currObj;

            void copy (const IT &in) {
               dir = in.dir; container = in.container; selIt = in.selIt;
               currObj = in.currObj;
            }
            void move (IT &&in) {
               dir = in.dir; container = in.container; selIt = in.selIt;
               currObj = in.currObj;
               in.dir = sx::Undefined; in.container = NULL;
            }

            IT insertElem (ssize_t newPos, const SValue &elem) {
               SX_CHECK (container);
               SX_UNUSED (newPos, elem);
               SX_EXIT;
            }
            IT prependElem (const SValue &elem) {
               SX_CHECK (container);
               SX_UNUSED (elem);
               SX_EXIT;
            }
            IT appendElem (const SValue &elem) {
               SX_CHECK (container);
               SX_UNUSED (elem);
               SX_EXIT;
            }

            // --- SX_ITERATOR / SX_CONST_ITERATOR callbacks
            void next () {
               SX_CHECK (container);
               ++selIt;
               if (valid ())
                  currObj = SxGPath(container->gPtr, *selIt);
            }
            void prev () {
               SX_CHECK (container);
               --selIt;
               if (valid ())
                  currObj = SxGPath(container->gPtr, *selIt);
            }
            bool valid () const { return (container != NULL && selIt.isValid ()); }
            SValue &getRef () {
               SX_CHECK (container);
               SX_CHECK (valid ());
               return currObj;
            }
            SValue *getPtr () {
               SX_CHECK (container);
               SX_CHECK (valid ());
               return &currObj;
            }
            bool equal (const IT &in) const {
               return (selIt == in.selIt && container == in.container);
            }

      };


      class Iterator : public State<T,Container,Iterator>
      {
         friend std::ostream &operator<< (std::ostream &s, const Iterator &in);
         SX_ITERATOR_NO_LAMBDAS (T,Container,Iterator)
         public:
            Iterator () : State<T,Container,Iterator> () { }
            Iterator (Container *c,
                      SxList<Selection>::Iterator it)
                      : State<T,Container,Iterator> (c, it)
                      { }

            Iterator (const Iterator &in, sx::ItCopyMode cmode_ = sx::CopyAll)
               : State<T,Container,Iterator> (in, cmode_) { }

            Iterator (Iterator &&in, sx::ItCopyMode cmode_ = sx::CopyAll)
               : State<T,Container,Iterator> (std::move(in), cmode_) { }
           ~Iterator () { }
      };

      class ConstIterator : public State<const T,const Container,ConstIterator>
      {
         friend std::ostream &operator<< (std::ostream &s, const ConstIterator &in);
         SX_CONST_ITERATOR_NO_LAMBDAS (T,Container,ConstIterator)
         public:
            ConstIterator ()
               : State<const T,const Container,ConstIterator> ()
                 { }
            ConstIterator (const Container *c,
                           SxList<Selection>::Iterator it)
               : State<const T,const Container,ConstIterator> (c, it)
                 { }

            ConstIterator (const ConstIterator &in,
                           sx::ItCopyMode cmode_ = sx::CopyAll)
               : State<const T,const Container,ConstIterator> (in, cmode_)
                 { }

            ConstIterator (ConstIterator &&in,
                           sx::ItCopyMode cmode_ = sx::CopyAll)
               : State<const T,const Container,ConstIterator>
                 (std::move(in), cmode_) { }

            // non-const to const cast
            ConstIterator (const Iterator &in)
               : State<const T,const Container,ConstIterator> (in)
                 { }

           ~ConstIterator () { }
      };


      SxGPathList (const SxPtr<SxGraph<SxGProps> > &gPtr_,
                   const SelSet &sels_);
      SxGPathList (const SxPtr<SxGraph<SxGProps> > &gPtr_,
                   const Selection &sel_);
     ~SxGPathList ();

      ssize_t getSize () const;
      ssize_t getPathSize () const;
      SxGPath operator() (ssize_t selIdx);

      // --------------------

      // Iteratorfunctions
      inline Iterator begin ();
      inline ConstIterator begin () const;
      inline Iterator begin (ssize_t idx);
      inline ConstIterator begin (ssize_t idx) const;
      inline Iterator end ();
      inline ConstIterator end () const;

      template<class Fn>
      void foreach (Fn fn);

      template<class Fn>
      void foreach (Fn fn) const;

   protected:
      SxPtr<SxGraph<SxGProps> > gPtr;
      SelSet sels;
};

inline std::ostream &operator<< (std::ostream &s,
                                 const SxGPathList::ConstIterator &in);

inline std::ostream &operator<< (std::ostream &s,
                                 const SxGPathList::Iterator &in);

#include <SxGPathList.hpp>

#endif /* _SX_G_PATH_LIST_H_ */
