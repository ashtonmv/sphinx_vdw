#ifndef _SX_RANGE_H_
#define _SX_RANGE_H_

#include <SxAlg.h>
#include <SxIterator.h>
#include <SxError.h>
#include <limits.h>

/** \brief SxRange class
 *
 *  This class is mainly a demonstrator of a MINIMAL data type used to
 *  define an SxIterator / SxConstIterator type.
 *
 * \author Sixten Boeck */
template<long Begin,long End,long Step=1>
class SxRange
{
   public:

      typedef SxRange<Begin,End,Step> Container;

      template<class Container,class IT>
      class State
      {
         // friend to State<const...>
         template<class CC,class CIT> friend class State;

         public:
            State (Container *c=NULL, long v=LONG_MAX)
               : container(c), val(v) { }

            // non-const to const cast
            template<class CC,class CIT>
            State (const State<CC,CIT> &in)
               : container(in.container), val(in.val) { }
			bool isForward() const { return true; }
         protected:
            Container *container;
            long val;

            // --- SX_ITERATOR / SX_CONST_OTERATOR callbacks
            void copy (const IT &in) { container=in.container; val=in.val; }
            void next () {
               val += ( (Begin<=End && Step>0) || (Begin>End && Step<0) )
                    ?  Step
                    : -Step;
               if (  (Begin <= End && (val > End || val < Begin))
                  || (Begin >  End && (val < End || val > Begin)))
                  val = LONG_MAX;
            }
            void prev () {
               val -= ( (Begin<=End && Step>0) || (Begin>End && Step<0) )
                    ?  Step
                    : -Step;
               if (  (Begin <= End && (val > End || val < Begin))
                  || (Begin >  End && (val < End || val > Begin)))
                  val = LONG_MAX;
            }
            bool valid () const { return val != LONG_MAX; }
            long &getRef () { return  val; }
            long *getPtr () { return &val; }
            bool equal (const IT &in) const { return val == in.val; }

      };

      class Iterator : public State<Container,Iterator>
      {
         SX_ITERATOR_NO_LAMBDAS(long,Container,Iterator)
         public:
            Iterator (Container *c=NULL, long v=LONG_MAX)
               : State<Container,Iterator>(c,v) { }
      };

      class ConstIterator : public State<const Container,ConstIterator>
      {
         SX_CONST_ITERATOR_NO_LAMBDAS(long,Container,ConstIterator)
         public:
         ConstIterator (const Container *c=NULL, long v=LONG_MAX)
            : State<const Container,ConstIterator> (c, v) { }
         // non-const to const cast
         ConstIterator (const Iterator &in)
            : State<const Container,ConstIterator> (in) { }
      };

      SX_CONTAINER_NO_LIST (Container)

      ConstIterator begin () const { return ConstIterator (this,Begin); }
           Iterator begin ()       { return Iterator (this,Begin); }
      ConstIterator end () const   { return ConstIterator (); }
           Iterator end ()         { return Iterator (); }
      ConstIterator fromLast () const { return begin (); }
      Iterator fromLast ()         { return begin (); }
      ConstIterator toFirst () const { return end (); }
      Iterator toFirst ()          { return end (); }
      ssize_t getSize ()           { return 0; }
};


#endif /* _SX_RANGE_H_ */
