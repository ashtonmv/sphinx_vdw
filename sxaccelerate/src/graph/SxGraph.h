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

#ifndef _SX_GRAPH_H_
#define _SX_GRAPH_H_

#include <SxUniqueList.h>
#include <SxArrayList.h>
#include <SxHashedArray.h>
#include <SxIterator.h>
#include <SxContainer.h>
#include <SxGraphStorage.h>
#include <SxPair.h>
#include <SxSet.h>
#include <SxMap.h>
#include <SxTime.h>
#include <SxPtr.h>
#include <SxAlg.h>
#include <SxSignals.h>


namespace sx
{
   enum GraphEvent { BeforeEdgeInEvent  = 0x01,
                     BeforeEdgeOutEvent = 0x02,
                     AfterEdgeInEvent   = 0x04,
                     AfterEdgeOutEvent  = 0x08 };
}

template<class T>
class SxGraphNode
{
   public:
      T element;
      SxArray<SxPair<ssize_t, ssize_t> > in;
      SxArray<SxPair<ssize_t, ssize_t> > out;
      ssize_t prev;
      ssize_t next;

      SxGraphNode () : in((ssize_t)0,(ssize_t)16),out((ssize_t)0,(ssize_t)16),
                       prev(-1),next(-1) { }
};

template<class T>
class SxGraphEdge
{
   public:
      T element;
      ssize_t in;
      ssize_t out;
      ssize_t prev;
      ssize_t next;
      SxGraphEdge (T val = T(), ssize_t in_ = -1, ssize_t out_ = -1,
                   ssize_t prev_ = -1, ssize_t next_ = -1)
                  : element (val), in (in_), out (out_),
                    prev(prev_), next(next_) { }

      ssize_t getOther (ssize_t v) {
         if (in == v)        return out;
         else if (out == v)  return in;
         else                SX_EXIT; // invalid vertex
      }
};


// Provides a default data type for edges.
// Contains auto incremented numeric value
// which results in unique hash value for each object.
class SxDefaultType {

   ssize_t val;
   public:
      SxDefaultType (): val(getVal()) {}
      ~SxDefaultType (){};
      operator ssize_t() const{ return val; }
      bool operator== (const SxDefaultType &in) { return val == in.val; }
      inline static ssize_t getVal() { static ssize_t vv(0); return vv++; }
};


template<class SValue,class SNode,class SContainer,class IT>
class SxGraphItState
{
   SX_ITERATOR_STATE
   // friend to State<const...>
   template<class CV,class CN,class CC,class CI> friend class SxGraphItState;

   public:

      using Selection = SxPtr<SxUniqueList<ssize_t> >;

      SxGraphItState ();
      SxGraphItState (sx::Direction dir_, SContainer *c, ssize_t idx_,
                      const Selection &sel);
      SxGraphItState (const SxGraphItState &in,
                      sx::ItCopyMode cMode_ = sx::CopyAll);
      SxGraphItState (SxGraphItState &&in, sx::ItCopyMode cMode_ = sx::CopyAll);
      // non-const to const cast
      template<class CV,class CN,class CC,class CI>
         SxGraphItState (const SxGraphItState<CV,CN,CC,CI> &in_)
         : container(in_.container), node(in_.node) { copy(in_); }

      SxGraphItState &operator= (const SxGraphItState &);
      SxGraphItState &operator= (SxGraphItState &&);

      // --- SxGraph specifics
      Selection getUnvisited () const;


      IT in  (ssize_t) const;
      IT out (ssize_t) const;
      SValue &getIn (ssize_t);                 // deprecate, use in(ssize_t)
      SValue &getOut (ssize_t);                // deprecate, use out(ssize_t)
      ssize_t getSizeIn () const;
      ssize_t getSizeOut () const;
      ssize_t findEdgeOut (ssize_t idx);
      ssize_t findEdgeIn (ssize_t idx);

      inline ssize_t getDepth () const { return pathLen; }
      inline SNode *getNode () { return node; }
      inline ssize_t getIdx () const { return idx; }
      inline typename SContainer::SelIdx getSelIdx () const {
         SX_CHECK (container);
         return idx;
      }

      IT neighbors (ssize_t maxDepth=1) const;   // excludes parent
      IT hops (ssize_t maxDepth) const;   // includes parent
      IT reverse () const;


   protected:

      SContainer *container;   // graph
      SNode *node;
      Selection selection;
      ssize_t idx;
      ssize_t distMapIdx;
      ssize_t maxDepth;  // 0 if unlimited
      ssize_t pathLen;
      typename SxList<ssize_t>::ConstIterator selectionIt;
      //SxList<SxPair<ssize_t,ssize_t> > distMap;
      SxHashedArray<SxPair<ssize_t,ssize_t> > distMap;
      SxSet<ssize_t> visited;

      // --- selected neighbors
      mutable ssize_t selectedSizeIn;
      mutable ssize_t selectedSizeOut;
      mutable SxArray<ssize_t> selectedIn;
      mutable SxArray<ssize_t> selectedOut;

      void procNext ();
      void procPrev ();

      // --- SX_ITERATOR / SX_CONST_ITERATOR callbacks
      template<class CIT> void copy (const CIT &,
                                     sx::ItCopyMode cMode_ = sx::CopyAll);
      void move (SxGraphItState &&, sx::ItCopyMode cMode_ = sx::CopyAll);

      IT insertElem (ssize_t newPos, const SValue &elem);
      IT prependElem (const SValue &elem);
      IT appendElem (const SValue &elem);
      void next ();
      void prev ();
      bool valid () const { return (container != NULL && node != NULL); }
      SValue &getRef ();
      SValue *getPtr ();
      bool equal (const IT &in) const;

};

template<class T,class Node,class Container>
class SxGraphIterators
{
   public:
      using Selection = SxPtr<SxUniqueList<ssize_t> >;

      template<class V,class N,class C,class IT>
      using State = SxGraphItState<V,N,C,IT>;

      class Iterator : public State<T,Node,Container,Iterator>
      {
         SX_ITERATOR(T,Container,Iterator)
         public:

            using EdgeEvent = SxSignal<const Iterator&,
                                       const Iterator&,
                                       const char *>;


            Iterator ()
               : State<T,Node,Container,Iterator> () { }
            Iterator (sx::Direction dir_, Container *c, ssize_t idx_,
                      const Selection &sel)
               : State<T,Node,Container,Iterator> (dir_, c, idx_, sel)
            { }
            Iterator (const Iterator &in, sx::ItCopyMode cMode_ = sx::CopyAll)
               : State<T,Node,Container,Iterator> (in, cMode_)
                 { }
            Iterator (Iterator &&in, sx::ItCopyMode cMode_ = sx::CopyAll)
               : State<T,Node,Container,Iterator>
                 (std::move(in), cMode_) { }

            // conversion from State
            Iterator (const State<T,Node,Container,Iterator> &in_)
               : State<T,Node,Container,Iterator> (in_) { }
      };

      class ConstIterator
         : public State<const T,const Node,const Container,ConstIterator>
      {
         SX_CONST_ITERATOR(T,Container,ConstIterator)
         public:

            ConstIterator ()
               : State<const T,const Node,const Container,ConstIterator> () { }
            ConstIterator (sx::Direction dir_, const Container *c, ssize_t idx_,
                           const Selection &sel)
               : State<const T,const Node,const Container,ConstIterator> (
                    dir_, c, idx_, sel) { }
            ConstIterator (const ConstIterator &in,
                           sx::ItCopyMode cMode_ = sx::CopyAll)
               : State<const T,const Node,const Container,ConstIterator>
                 (in, cMode_) { }
            ConstIterator (ConstIterator &&in,
                           sx::ItCopyMode cMode_ = sx::CopyAll)
               : State<const T,const Node,const Container,ConstIterator>
                 (std::move(in), cMode_) { }

            // conversion from State
            ConstIterator (
               const State<const T,const Node,const Container,ConstIterator> &in_)
               : State<const T,const Node,const Container,ConstIterator> (in_)
                 { }

            // non-const to const cast
            ConstIterator (const Iterator &in_)
               : State<const T,const Node,const Container,ConstIterator>(in_)
                 { }

      };
};


template<class T, template<class> class V>
class SxGraphStorage
{
   template<class GT,class E,template<class> class NT,template<class> class ET,
            template<class,class,class> class ItPair> friend class SxGraph;
   template<class CV,class CN,class CC,class CI> friend class SxGraphItState;

   protected:

      using ElemType = V<T>;

      //SxMap<ssize_t, ElemType> elems;
      SxHashedArray<ElemType>  elems;
      SxArray<ssize_t>         hashTable;

      ssize_t firstElement;
      ssize_t lastElement;
      ssize_t firstFree;

      ssize_t nElems;

      ssize_t hashSize;
      ssize_t hashUsed;
      ssize_t hashBound;

   public:

      SX_GRAPH_STORAGE(T,V);

      SxGraphStorage ();
      SxGraphStorage (const SxGraphStorage<T,V> &);
      SxGraphStorage (SxGraphStorage<T,V> &&);
      SxGraphStorage<T,V> &operator= (const SxGraphStorage<T,V> &);
      SxGraphStorage<T,V> &operator= (SxGraphStorage<T,V> &&);
      ~SxGraphStorage ();

      ssize_t     addElem (const T &);
      ssize_t     addElem (T &&);
      void        removeElem (ssize_t idx);

      bool    containsElem (ssize_t idx) const;
      size_t  getNBytes () const;
      void    removeElems ();

            V<T> &getElem (ssize_t idx);
      const V<T> &getElem (ssize_t idx) const;
      ssize_t     getIdx (const T &);

      // --- hash table
      ssize_t hashFindPos (const T &) const;
      bool    hashContains (const T &, ssize_t *);
      void    hashResize (ssize_t newSize_);
};

template<class T,class E = SxDefaultType,template<class> class NT = SxGraphNode,
         template<class> class ET = SxGraphEdge,
         template<class,class,class> class ItPair=SxGraphIterators>
class SxGraph : public SxThis<SxGraph<T,E,NT,ET,ItPair> >
{
   public:
      template<class CV,class CN,class CC,class CI> friend class SxGraphItState;

      using Node = NT<T>;
      using Edge = ET<E>;

      // --- Iterator
      using Container     = SxGraph<T,E,NT,ET,ItPair>;
      using Iterator      = typename ItPair<T,Node,Container>::Iterator;
      using ConstIterator = typename ItPair<T,Node,Container>::ConstIterator;



      using Selection = SxPtr<SxUniqueList<ssize_t> >;


      // --- Event API
      using EdgeEvent = SxSignal<const Iterator&,
                                 const Iterator&,
                                 const char *>;

      // --- SxSelection:
      typedef T                     TElem;
      typedef ssize_t               SelIdx;
      typedef SxUniqueList<SelIdx>  SelStorage;
      typedef Container             SelContainer;

      SX_CONTAINER (Container)

      SxGraph ();
      SxGraph (SxPtr<SxGraphStorage<T, NT> >, SxPtr<SxGraphStorage<E, ET> >);
      SxGraph (const SxGraph<T,E,NT,ET,ItPair> &);
      SxGraph (SxGraph<T,E,NT,ET,ItPair> &&);
      ~SxGraph ();

      SxGraph<T,E,NT,ET,ItPair> &operator=(const SxGraph<T,E,NT,ET,ItPair> &);
      SxGraph<T,E,NT,ET,ItPair> &operator=(SxGraph<T,E,NT,ET,ItPair> &&);

      // -- Insertion
      ConstIterator createNode (const T &);
      ConstIterator createNode (T &&);
      void createEdge (const T &from, const T &to, const E &elem=E());
      void createEdge (ConstIterator &, ConstIterator &, const E &elem=E());

      // -- Deletion
      void removeElement (const T &);
      void removeEdge (const T &, const T &);
      void removeSelection (const Selection &selection);
      void removeAll ();

      // -- Search
      bool containsNode (const T &) const;
      ssize_t findPos (const T &) const;

      ssize_t findPath      (const T         &from,
                             const T         &to,
                             ssize_t         maxPathLen=0,
                             const Selection &selection=Selection()) const;

      SxList<ConstIterator>
         getPath (const T         &from,
                  const T         &to,
                  ssize_t         maxPathLen=0,
                  const Selection &selection=Selection()) const;

      SxList<ConstIterator>
         getPath (const ConstIterator &from,
                  const ConstIterator &to,
                  ssize_t             maxPathLen=0,
                  const Selection     &selection=Selection()) const;


      SxList<T> topsort (const Selection &selection=Selection()) const;

      SxArray<SxArray<T> > getCycles (const Selection &selection=Selection()) const;
      SxList<SxList<ssize_t> > getCyclesIdx(const Selection &selection=Selection()) const;


      // -- Size
      inline ssize_t getSize () const { return nodes->nElems; }
      inline ssize_t getNEdges () const { return edges->nElems; }
      size_t getNBytes () const;

      // Iteratorfunctions
      Iterator      begin (ssize_t nodeIdx, sx::Direction dir_=sx::Forward);
      ConstIterator begin (ssize_t nodeIdx, sx::Direction dir_=sx::Forward) const;

      Iterator      begin (const Selection &selection=Selection());
      ConstIterator begin (const Selection &selection=Selection()) const;
      Iterator      begin (sx::Direction dir_, const T      &from,
                           const Selection &selection=Selection());
      ConstIterator begin (sx::Direction dir_, const T      &from,
                           const Selection &selection=Selection()) const;
      Iterator      begin (const T &from,
                          const Selection &selection=Selection());
      ConstIterator begin (const T &from,
                           const Selection &selection=Selection()) const;
      Iterator beginIn   (const T &from,
                          const Selection &selection=Selection());
      Iterator beginBoth (const T &from,
                          const Selection &selection=Selection());
      Iterator end ();

      ConstIterator beginIn   (const T         &from,
                               const Selection &selection=Selection()) const;
      ConstIterator beginBoth (const T         &from,
                               const Selection &selection=Selection()) const;
      ConstIterator end () const;

      Iterator fromLast ()  {SX_EXIT;}
      ConstIterator fromLast () const {SX_EXIT;}

      Iterator toFirst ()  {SX_EXIT;}
      ConstIterator toFirst () const {SX_EXIT;}


      // --- SxSelection functions
      inline Iterator      getIterator (const SelIdx &);
      inline ConstIterator getIterator (const SelIdx &) const;
      inline ConstIterator getConstIterator (const Iterator &) const;
      SxPtr<Container> getContainer () const;

      // --- event handling
      EdgeEvent &getSignal (const Iterator &,sx::GraphEvent);

      void print (Iterator, std::ostream &, const SxString &options="") const;
      void print (ConstIterator, std::ostream &, const SxString &options="") const;
      void printDebug () const;

      void print () const;

      void print (const SxUniqueList<ssize_t> &visited,
                  std::ostream                &,
                  const SxString              &options="") const;

      // --- findPath
      void nextLevel (ssize_t                          ptr,
                      ssize_t                          pathLen,
                      ssize_t                          maxPathLen,
                      const Selection                  &selection,
                      SxList<SxPair<ssize_t,ssize_t> > *queue,
                      SxSet<ssize_t>                   *visited) const;

      // --- getCycles
      bool circuit (ssize_t                        idx,
                    ssize_t                        start,
                    SxList<ssize_t>                &stack,
                    SxArray<bool>                  &visited,
                    SxArray<bool>                  &blocked,
                    SxArray<SxArray<ssize_t> > &B,
                    SxList<SxList<ssize_t> >       *circuits) const;
      void unblock (ssize_t, SxArray<bool> &, SxArray<SxArray<ssize_t> >&) const;

      void removeNode (ssize_t nodeIdx);

      bool hasOutEdge (ssize_t, ssize_t) const;
      ssize_t findEdgeOut (ssize_t, ssize_t) const;
      ssize_t findEdgeIn (ssize_t, ssize_t) const;
      void unlinkOut (ssize_t, ssize_t);
      void unlinkIn (ssize_t, ssize_t);



   protected:

      // --- events
      SxMap<ssize_t,EdgeEvent> sigBeforeEdgeIn;
      SxMap<ssize_t,EdgeEvent> sigBeforeEdgeOut;
      SxMap<ssize_t,EdgeEvent> sigAfterEdgeIn;
      SxMap<ssize_t,EdgeEvent> sigAfterEdgeOut;

      SxPtr<SxGraphStorage<T, NT> > nodes;
      SxPtr<SxGraphStorage<E, ET> > edges;

      void handleEdgeEvent (unsigned int e, ssize_t from, ssize_t to);


};


#include <SxGraph.hpp>

#endif /* _SX_GRAPH_H_ */
