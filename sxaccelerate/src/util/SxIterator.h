#ifndef _SX_ITERATOR_H_
#define _SX_ITERATOR_H_

#include <SxConfig.h>
#include <cstddef>  // NULL

#define SX_ITERATOR_STATE                                                    \
   public:                                                                   \
      sx::Direction getDirection () const {return this->dir; }               \
      bool isForward () const { return (this->dir == sx::Forward); }         \
      void setForward ()      { this->dir = sx::Forward; }                   \
      void setBackward ()     { this->dir = sx::Backward; }                  \
      void setReverse () {                                                   \
         this->dir = (this->dir == sx::Forward)? sx::Backward : sx::Forward; \
      }                                                                      \
   protected:                                                                \
      sx::Direction dir;                                                     \

#ifndef NDEBUG
#   define SX_ITERATOR_DBG(Val,Container,Iterator)                           \
      template<class Fn>                                                     \
      Iterator &sxDbgMsg (const char *file, long line,                       \
                          const char *func,const char *logId, Fn f)          \
      {                                                                      \
         f (file, line, func, logId, *this);                                 \
         return this->getIterator ();                                        \
      }                                                                      \
      Iterator &sxBreak (const char *file, long line,                        \
                         const char *func)                                   \
      {                                                                      \
         ::sxBreak (file, line, func);                                       \
         return this->getIterator ();                                        \
      }
#else
#   define SX_ITERATOR_DBG(Val,Container,Iterator)                           \
      template<class Fn>                                                     \
      Iterator &sxDbgMsg (const char *, long, const char *,                  \
                          const char *, Fn)                                  \
      {                                                                      \
         return this->getIterator ();                                        \
      }                                                                      \
      Iterator &sxBreak (const char *, long, const char *)                   \
      {                                                                      \
         return this->getIterator ();                                        \
      }
#endif


#define SX_ITERATOR_LAMBDAS(Val,Container,Iterator)                          \
      template<class Fn_> void foreach (Fn_ fn_) {                           \
         SX_CHECK (this->container);                                         \
         sx::foreach (*this, this->container->end(),fn_);                    \
      }                                                                      \
      template<class Fn_> Iterator find (Fn_ fn_) {                          \
         SX_CHECK (this->container);                                         \
         return sx::find (*this, this->container->end(), fn_);               \
      }                                                                      \
      template<class Fn_> Iterator findCond (Fn_ fn_) {                      \
         SX_CHECK (this->container);                                         \
         return sx::findCond (*this, this->container->end(), fn_);           \
      }                                                                      \
      template<class Fn_>                                                    \
      SxSelection<typename Container::SelContainer> findAll (Fn_ fn_) const {\
         SX_CHECK (this->container);                                         \
         SxCList<typename Container::SelIdx> lst;                            \
         auto it = this->container->begin ();                                \
         for (;it != this->container->end (); ++it) {                        \
            if (fn_ (this->container->getConstIterator(it))) {               \
               lst.append (it.getSelIdx ());                                 \
            }                                                                \
         }                                                                   \
         SxSelection<typename Container::SelContainer> sel(lst);             \
         sel.container = this->container->getContainer ();                   \
         return sel;                                                         \
      }



#define SX_CONST_ITERATOR_LAMBDAS(Val,Container,ConstIterator)               \
      template<class Fn_> void foreach (Fn_ fn_) const {                     \
         SX_CHECK (this->container);                                         \
         sx::foreach (*this, this->container->end(), fn_);                   \
      }                                                                      \
      template<class Fn_> ConstIterator find (Fn_ fn_) const {               \
         SX_CHECK (this->container);                                         \
         return sx::find (*this, this->container->end(), fn_);               \
      }                                                                      \
      template<class Fn_> ConstIterator findCond (Fn_ fn_) const {           \
         SX_CHECK (this->container);                                         \
         return sx::findCond (*this, this->container->end(), fn_);           \
      }                                                                      \
      template<class Fn_>                                                    \
      SxSelection<typename Container::SelContainer> findAll (Fn_ fn_) const {\
         SX_CHECK (this->container);                                         \
         SxCList<typename Container::SelIdx> lst;                            \
         auto it = this->container->begin ();                                \
         for (;it != this->container->end (); ++it) {                        \
            if (fn_ (this->container->getConstIterator(it))) {               \
               lst.append (it.getSelIdx ());                                 \
            }                                                                \
         }                                                                   \
         SxSelection<typename Container::SelContainer> sel(lst);             \
         sel.container = this->container->getContainer ();                   \
         return sel;                                                         \
      }


// Using macros to simplyfy defining individual SxIterator types

/** \brief Generic non-const iterator
 *
 *  This version is w/o lambda function support
 *
 *  \author Sixten Boeck */
#define SX_ITERATOR_NO_LAMBDAS(Val,Container,Iterator)                       \
   protected:                                                                \
      Iterator &getIterator () { return *(Iterator *)this; }                 \
                                                                             \
   public:                                                                   \
                                                                             \
      typedef Val Value;                                                     \
                                                                             \
      bool operator== (const Iterator &in_) const {                          \
         return this->equal (in_);                                           \
      }                                                                      \
      bool operator!= (const Iterator &in_) const {                          \
         return !this->equal (in_);                                          \
      }                                                                      \
      Iterator &operator= (const Iterator &in_) {                            \
         if (this == &in_) return  this->getIterator ();                     \
         this->copy (in_);                                                   \
         return this->getIterator ();                                        \
      }                                                                      \
      Iterator &operator= (Iterator &&in_) {                                 \
         if (this == &in_) return this->getIterator ();                      \
         this->move (std::move(in_));                                        \
         return this->getIterator ();                                        \
      }                                                                      \
      bool operator! () const { return !this->valid(); }                     \
      operator bool () const { return this->valid(); }                       \
      bool isValid () const { return this->valid(); }                        \
                                                                             \
      /* postfix++ / postfix--  */                                           \
      Iterator operator++ (int)  {                                           \
         Iterator clone (this->getIterator());                               \
         this->next ();                                                      \
         return clone;                                                       \
      }                                                                      \
      Iterator operator-- (int)  {                                           \
         Iterator clone (this->getIterator());                               \
         this->prev ();                                                      \
         return clone;                                                       \
      }                                                                      \
                                                                             \
      /* ++prefix / --prefix */                                              \
      Iterator &operator++ () { this->next (); return this->getIterator(); } \
      Iterator &operator-- () { this->prev (); return this->getIterator(); } \
                                                                             \
      const Value &operator*  () const {                                     \
         return const_cast<Iterator *>(this)->getRef();                      \
      }                                                                      \
      Value &operator*  () {                                                 \
         return this->getRef();                                              \
      }                                                                      \
      const Value *operator-> () const {                                     \
         return const_cast<Iterator *>(this)->getPtr();                      \
      }                                                                      \
      Value *operator-> () {                                                 \
         return this->getPtr();                                              \
      }                                                                      \
      Iterator &print (std::ostream &s=std::cout) {                          \
         s << *this;                                                         \
         return this->getIterator();                                         \
      }                                                                      \
      Iterator &print (std::ostream &s = std::cout) const {                  \
         s << *this;                                                         \
         return const_cast<Iterator *>(this)->getIterator();                 \
      }                                                                      \
      template<class Elem>                                                   \
      Iterator insert (ssize_t newPos, const Elem &elem) {                   \
         return this->insertElem (newPos, elem);                             \
      }                                                                      \
      template<class Elem>                                                   \
      Iterator prepend (const Elem &elem) {                                  \
         return this->prependElem (elem);                                    \
      }                                                                      \
      template<class Elem>                                                   \
      Iterator append (const Elem &elem) {                                   \
         return this->appendElem (elem);                                     \
      }                                                                      \
      Iterator forward () {                                                  \
         Iterator res(*this, sx::CopyItData);                                \
         res.setForward ();                                                  \
         return res;                                                         \
      }                                                                      \
      Iterator backward () {                                                 \
         Iterator res(*this, sx::CopyItData);                                \
         res.setBackward ();                                                 \
         return res;                                                         \
      }                                                                      \
      Iterator reverse () {                                                  \
         Iterator res(*this, sx::CopyItData);                                \
         res.setReverse ();                                                  \
         return res;                                                         \
      }                                                                      \
      template<class Fn>                                                     \
      Iterator &sxLog (const char *logId, uint32_t logHash,                  \
                       const char *file, long line,                          \
                       const char *func, Fn f)                               \
      {                                                                      \
         f (logId, logHash, file, line, func, *this);                        \
         return this->getIterator ();                                        \
      }                                                                      \
      SX_ITERATOR_DBG(Val,Container,Iterator)

#define SX_ITERATOR(Val,Container,Iterator)                                  \
      SX_ITERATOR_NO_LAMBDAS(Val,Container,Iterator)                         \
      SX_ITERATOR_LAMBDAS(Val,Container,Iterator)


/** \brief Generic const iterator
 *
 *  This version is w/o lambda function support
 *
 *  \author Sixten Boeck */
#define SX_CONST_ITERATOR_NO_LAMBDAS(Val,Container,ConstIterator)            \
   protected:                                                                \
      ConstIterator &getIterator () { return *(ConstIterator *)this; }       \
                                                                             \
   public:                                                                   \
                                                                             \
      typedef Val Value;                                                     \
                                                                             \
      bool operator== (const ConstIterator &in_) const {                     \
         return this->equal (in_);                                           \
      }                                                                      \
                                                                             \
      bool operator!= (const ConstIterator &in_) const {                     \
         return !this->equal (in_);                                          \
      }                                                                      \
      ConstIterator &operator= (const ConstIterator &in_) {                  \
         if (this == &in_) return  this->getIterator ();                     \
         this->copy (in_);                                                   \
         return this->getIterator ();                                        \
      }                                                                      \
      ConstIterator &operator= (ConstIterator &&in_) {                       \
         if (this == &in_) return this->getIterator ();                      \
         this->move (std::move(in_));                                        \
         return this->getIterator ();                                        \
      }                                                                      \
      bool operator! () const { return !this->valid(); }                     \
      operator bool () const { return this->valid(); }                       \
      bool isValid () const { return this->valid(); }                        \
                                                                             \
      /* postfix++ / postfix-- */                                            \
      ConstIterator operator++ (int)  {                                      \
         ConstIterator clone (this->getIterator());                          \
         this->next ();                                                      \
         return clone;                                                       \
      }                                                                      \
      ConstIterator operator-- (int)  {                                      \
         ConstIterator clone (this->getIterator());                          \
         this->prev ();                                                      \
         return clone;                                                       \
      }                                                                      \
                                                                             \
      /* ++prefix / --prefix */                                              \
      ConstIterator &operator++ () {                                         \
         this->next ();                                                      \
         return this->getIterator();                                         \
      }                                                                      \
      ConstIterator &operator-- () {                                         \
         this->prev ();                                                      \
         return this->getIterator();                                         \
      }                                                                      \
                                                                             \
      const Value &operator*  () const {                                     \
         return const_cast<ConstIterator *>(this)->getRef();                 \
      }                                                                      \
      const Value *operator-> () const {                                     \
         return const_cast<ConstIterator *>(this)->getPtr();                 \
      }                                                                      \
      ConstIterator &print (std::ostream &s=std::cout) {                     \
         s << *this;                                                         \
         return this->getIterator();                                         \
      }                                                                      \
      ConstIterator &print (std::ostream &s = std::cout) const {             \
         s << *this;                                                         \
         return const_cast<ConstIterator *>(this)->getIterator();            \
      }                                                                      \
      ConstIterator begin () const {                                         \
         SX_CHECK (this->container);                                         \
         return const_cast<ConstIterator *>(this)->container->begin ();      \
      }                                                                      \
      ConstIterator end () const {                                           \
         SX_CHECK (this->container);                                         \
         return const_cast<ConstIterator *>(this)->container->end ();        \
      }                                                                      \
      ConstIterator forward () {                                             \
         ConstIterator res(*this, sx::CopyItData);                           \
         res.setForward ();                                                  \
         return res;                                                         \
      }                                                                      \
      ConstIterator backward () {                                            \
         ConstIterator res(*this, sx::CopyItData);                           \
         res.setBackward ();                                                 \
         return res;                                                         \
      }                                                                      \
      ConstIterator reverse () {                                             \
         ConstIterator res(*this, sx::CopyItData);                           \
         res.setReverse ();                                                  \
         return res;                                                         \
      }                                                                      \



#define SX_CONST_ITERATOR(Val,Container,Iterator)                            \
      SX_CONST_ITERATOR_NO_LAMBDAS(Val,Container,Iterator)                   \
      SX_CONST_ITERATOR_LAMBDAS(Val,Container,Iterator)                      \
      SX_ITERATOR_DBG(Val,Container,ConstIterator)




#endif /* _SX_ITERATOR_H_ */
