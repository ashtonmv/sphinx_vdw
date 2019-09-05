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

#ifndef _SX_GRAPH_STORAGE_H_
#define _SX_GRAPH_STORAGE_H_

#define SX_GRAPH_STORAGE(TT,VV)                                               \
   public:                                                                    \
      ssize_t add (const TT &elem)                                            \
      {                                                                       \
         return this->addElem (elem);                                         \
      }                                                                       \
      ssize_t add (TT &&elem)                                                 \
      {                                                                       \
         return this->addElem (std::move(elem));                              \
      }                                                                       \
      void remove (ssize_t idx)                                               \
      {                                                                       \
         this->removeElem (idx);                                              \
      }                                                                       \
      VV<TT> &get (ssize_t idx)                                               \
      {                                                                       \
         return this->getElem (idx);                                          \
      }                                                                       \
      const VV<TT> &get (ssize_t idx) const                                   \
      {                                                                       \
         return this->getElem (idx);                                          \
      }                                                                       \
      ssize_t findPos (const TT &elem)                                        \
      {                                                                       \
         return this->getIdx (elem);                                          \
      }                                                                       \
      bool contains (ssize_t idx)                                             \
      {                                                                       \
         return this->containsElem (idx);                                     \
      }                                                                       \
      void removeAll () {                                                     \
         this->removeElems ();                                                \
      }

#endif /*_SX_GRAPH_STORAGE_H_*/
