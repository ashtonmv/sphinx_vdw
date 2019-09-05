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

#ifndef _SX_SELECTION_H_
#define _SX_SELECTION_H_

#include <SxAlg.h>
#include <SxIterator.h>
#include <SxCList.h>
#include <SxContainer.h>
#include <SxPtr.h>

/** \brief Selection type for all containers

    \b SxSelection = SPHInX Container Selection Class

    This class can be used to apply selection filters on
    container objects by using the where(..) function that accepts
    a lambda function. The lambda function decides for each
    item in the container if it satisfies the filter criteria.
    The resulting SxSelection object can be used
    to apply further filters on the selection.


\code
#include <SxArray.h>
#include <SxSelection.h>
...

SxArray<int> a(5);

for (int i = 0; i < 5; ++i) {
   a(i) = i+1;
}

SxSelection<SxArray<T> > sel;
sel = a.where ([](auto &it) { return *it % 2 == 0; });

sel.foreach ([](auto it) { cout << *it << " ";});
\endcode

or simpler

\code
a.where (...).foreach ([](auto it) { cout << *it...})
\endcode
 */
template<class C>
class SxSelection
{
   public:
      typedef typename C::TElem           TElem;
      typedef typename C::SelIdx          SelIdx;
      typedef SxSelection<C>              Container;
      typedef typename C::SelStorage      SelStorage;
      typedef typename C::SelContainer    SelContainer;

      template<class> friend class SxArray;
      template<class, template<class,class,class> class> friend class SxList;
      template<class,class,class> friend class SxMap;
      template<class,class,template<class> class,template<class> class,
               template<class,class,class> class> friend class SxGraph;


      template<class SValue,class Container,class SIT,class IT>
      class State
      {
         SX_ITERATOR_STATE
         // frient to State<const...>
         template<class CV,class CC,class CSIT,class CIT> friend class State;

         public:
            State (sx::Direction dir_=sx::Forward, Container *c=NULL,
                   const SIT &sIt_= SIT())
                   : dir(dir_), container(c), it(sIt_) { }

            State (const State &in, sx::ItCopyMode cmode_ = sx::CopyAll)
               : dir(in.dir), container(in.container), it(in.it) {
               SX_UNUSED (cmode_);
            }

            State (State &&in, sx::ItCopyMode cmode_ = sx::CopyAll)
               : dir(in.dir), container(in.container), it(in.it)
            {
               SX_UNUSED (cmode_);
               in.dir = sx::Undefined; in.container = NULL;
            }
            // non-const to const cast
            template<class CV,class CC,class CSIT,class CIT>
            State (const State<CV,CC,CSIT,CIT> &in)
               : dir(in.dir),container(in.container), it(in.it) { }

            typename Container::SelIdx getSelIdx () const {
               SX_CHECK (container);
               SX_CHECK (it.isValid ());
               return *it;
            }

         protected:
            Container *container;
            SIT it;

            void copy (const IT &in) {
               dir = in.dir; container = in.container; it = in.it;
            }
            void move (IT &&in) {
               dir = in.dir; container = in.container; it = in.it;
               in.dir = sx::Undefined; in.container = NULL;
            }
            IT insertElem (ssize_t newPos, const SValue &elem) {
               SX_UNUSED (newPos, elem);
               SX_EXIT;
            }
            IT prependElem (const SValue &elem) {
               SX_UNUSED (elem);
               SX_EXIT;
            }
            IT appendElem (const SValue &elem) {
               SX_UNUSED (elem);
               SX_EXIT;
            }
            // --- SX_ITERATOR / SX_CONST_ITERATOR callbacks
            void next () {
               SX_CHECK (container);
               ++it;
            }
            void prev () {
               SX_CHECK (container);
               --it;
            }
            bool valid () const { return it.isValid (); }
            SValue &getRef () {
               SX_CHECK (container);
               SX_CHECK (it.isValid ());
               return *(container->getIterator(it));
            }
            SValue *getPtr () {
               SX_CHECK (container);
               SX_CHECK (it.isValid ());
               return ((container->getIterator(it)).operator->());
            }
            bool equal (const IT &in) const {
               return (it == in.it && container == in.container);
            }
      };

      class Iterator : public State<TElem,Container,typename SelStorage::Iterator,Iterator>
      {
         SX_ITERATOR (TElem,Container,Iterator)
         typedef typename SelStorage::Iterator SIT;
         public:
            Iterator (sx::Direction dir_=sx::Forward, Container *c=NULL,
                      const SIT &sIt_=SIT())
                      : State<TElem,Container,SIT,Iterator> (dir_, c, sIt_)
                      { }

            Iterator (const Iterator &in, sx::ItCopyMode cmode_ = sx::CopyAll)
               : State<TElem,Container,SIT,Iterator> (in, cmode_) { }

            Iterator (Iterator &&in, sx::ItCopyMode cmode_ = sx::CopyAll)
               : State<TElem,Container,SIT,Iterator> (std::move(in), cmode_) { }
      };

      class ConstIterator : public State<const TElem,const Container,
                                         typename SelStorage::ConstIterator,
                                         ConstIterator>
      {
         SX_CONST_ITERATOR (TElem,Container,ConstIterator)
         typedef typename SelStorage::ConstIterator CSIT;
         public:
            ConstIterator (sx::Direction dir_=sx::Forward,
                           const Container *c=NULL, const CSIT &sIt_=CSIT())
               : State<const TElem,const Container,CSIT,
                       ConstIterator> (dir_, c, sIt_)
                 { }

            ConstIterator (const ConstIterator &in,
                           sx::ItCopyMode cmode_ = sx::CopyAll)
               : State<const TElem,const Container,CSIT,
                       ConstIterator> (in, cmode_)
                 { }

            ConstIterator (ConstIterator &&in,
                           sx::ItCopyMode cmode_ = sx::CopyAll)
               : State<const TElem,const Container,CSIT,ConstIterator>
                 (std::move(in), cmode_) { }

            // non-const to const cast
            ConstIterator (const Iterator &in)
               : State<const TElem,const Container,CSIT,ConstIterator> (in)
                 { }
      };

      SX_CONTAINER(Container)

      SxSelection ();
      SxSelection (const SxCList<SelIdx> &);
      SxSelection (const SxPtr<SelStorage> &);
     ~SxSelection ();

      void setContainer (SxPtr<C> ptr);

      ssize_t getSize () const;
      SxPtr<C> getContainer () const;


      // sort two lists(selected elems and their indices) at once
      template<class Iterator1, class Iterator2>
      void selSort (const Iterator1 &it1, const Iterator1 &lastIt1,
                    const Iterator2 &it2, const Iterator2 &lastIt2,
                    SxCBoundPtr<int, const decltype(Iterator1())&,
                    const decltype(Iterator1())&> comp);


      // sort two lists(selected elems and their indices) at once
      template<class Iterator1, class Iterator2>
      void selSort (const Iterator1 &it1, const Iterator1 &lastIt1,
                    const Iterator2 &it2, const Iterator2 &lastIt2);


      void removeDuplicates (SxCList<TElem> *lst, SxCList<SelIdx> *lstIdx);

      // compute union
      SxSelection<C> operator|  (const SxSelection<C> &in);
      // compute intersection
      SxSelection<C> operator&  (const SxSelection<C> &in);
      // compute complement
      SxSelection<C> operator-  (const SxSelection<C> &in);
      // compute symmetric complement
      SxSelection<C> operator!= (const SxSelection<C> &in);

      /** --- fetch container iterators */
      inline typename C::Iterator getIterator (const typename SelStorage::Iterator &it);
      inline typename C::ConstIterator getIterator (const typename SelStorage::ConstIterator &it) const;
      inline typename C::ConstIterator getConstIterator (const Iterator &it) const;

      /** --- Iterators */
      ConstIterator begin () const;
      Iterator begin ();

      ConstIterator begin (ssize_t idx_,
                           sx::Direction dir_ = sx::Forward) const;
      Iterator begin (ssize_t idx_,
                      sx::Direction dir_ = sx::Forward);


      ConstIterator end () const;
      Iterator end ();

      ConstIterator fromLast () const;
      Iterator fromLast ();

      ConstIterator toFirst () const;
      Iterator toFirst ();

   protected:
      /** storage for selected items */
      SxPtr<SelStorage> selected;

      /** pointer to original container */
      SxPtr<C> container;

};

#include <SxSelection.hpp>

#endif /*_SX_SELECTION_H_*/

