#ifndef _SX_CONTAINER_H_
#define _SX_CONTAINER_H_

#define SX_CONTAINER_NO_LIST(Container)                                      \
   public:                                                                   \
      template<class Fn> void foreach (Fn fn)                                \
      {                                                                      \
         sx::foreach (begin(), end(), fn);                                   \
      }                                                                      \
      template<class Fn> void foreach (Fn fn) const                          \
      {                                                                      \
         sx::foreach (begin(), end(), fn);                                   \
      }                                                                      \
      template<class Elem> typename Container::Iterator find (Elem elem)     \
      {                                                                      \
         return sx::find (begin(), end(), elem);                             \
      }                                                                      \
      template<class Elem>                                                   \
      typename Container::ConstIterator find (Elem elem) const               \
      {                                                                      \
         return sx::find (begin(), end(), elem);                             \
      }                                                                      \
      template<class Cond>                                                   \
      typename Container::Iterator findCond (Cond cond)                      \
      {                                                                      \
         return sx::findCond (begin(), end(), cond);                         \
      }                                                                      \
      template<class Cond>                                                   \
      typename Container::ConstIterator findCond (Cond cond) const           \
      {                                                                      \
         return sx::findCond (begin(), end(), cond);                         \
      }                                                                      \
      void sort (SxCBoundPtr<int, const typename Container::Iterator&,       \
                 const typename Container::Iterator&> comp)                  \
      {                                                                      \
         if (getSize () == 0)  return;                                       \
         sx::sort (begin(), fromLast(), comp);                               \
      }                                                                      \
      void sort ()                                                           \
      {                                                                      \
         if (getSize () == 0)  return;                                       \
         sx::sort (begin(), fromLast());                                     \
      }                                                                      \
      void qsort (SxCBoundPtr<int, const typename Container::Iterator&,      \
                  const typename Container::Iterator&> comp)                 \
      {                                                                      \
         if (getSize () == 0) return;                                        \
         sx::qsort (begin (), fromLast (), comp);                            \
      }                                                                      \
      void qsort ()                                                          \
      {                                                                      \
         if (getSize () == 0)  return;                                       \
         sx::qsort (begin(), fromLast());                                    \
      }                                                                      \
      void sortDesc (SxCBoundPtr<int, const typename Container::Iterator&,   \
                     const typename Container::Iterator&> comp)              \
      {                                                                      \
         if (getSize () == 0)  return;                                       \
         sx::sort (fromLast(), begin(), comp);                               \
      }                                                                      \
      void sortDesc ()                                                       \
      {                                                                      \
         if (getSize () == 0)  return;                                       \
         sx::sort (fromLast(), begin());                                     \
      }                                                                      \
      void qsortDesc (SxCBoundPtr<int, const typename Container::Iterator&,  \
                      const typename Container::Iterator&> comp)             \
      {                                                                      \
         if (getSize () == 0)  return;                                       \
         sx::qsort (fromLast(), begin(), comp);                              \
      }                                                                      \
      void qsortDesc ()                                                      \
      {                                                                      \
         if (getSize () == 0)  return;                                       \
         sx::qsort (fromLast(), begin());                                    \
      }                                                                      \

#define SX_CONTAINER(Container)                                              \
      SX_CONTAINER_NO_LIST(Container)                                        \
      SxSelection<typename Container::SelContainer>                          \
      where (SxCBoundPtr<bool,                                               \
             typename Container::SelContainer::ConstIterator> f)             \
      {                                                                      \
         SxCList<typename Container::SelIdx> lst;                            \
         auto it = begin ();                                                 \
         for (;it != end (); ++it) {                                         \
            if (f (getConstIterator(it))) {                                  \
               lst.append (it.getSelIdx ());                                 \
            }                                                                \
         }                                                                   \
         SxSelection<typename Container::SelContainer> sel(lst);             \
         sel.container = this->getContainer ();                              \
         return sel;                                                         \
      }                                                                      \
      template<class Cond>                                                   \
      SxSelection<typename Container::SelContainer> findAll (Cond cond)      \
      {                                                                      \
         SxCList<typename Container::SelIdx> lst;                            \
         auto it = begin ();                                                 \
         for (;it != end (); ++it) {                                         \
            if (cond (getConstIterator (it))) {                              \
               lst.append (it.getSelIdx ());                                 \
            }                                                                \
         }                                                                   \
         SxSelection<typename Container::SelContainer> sel(lst);             \
         sel.container = this->getContainer ();                              \
         return sel;                                                         \
      }                                                                      \




#endif /* _SX_CONTAINER_H_ */
