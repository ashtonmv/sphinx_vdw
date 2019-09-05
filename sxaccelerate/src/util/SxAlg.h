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

#ifndef _SX_ALG_H_
#define _SX_ALG_H_
#include <iostream>
#include <SxCBoundPtr.h>
/** \brief Simple foreach implementation

    Usage:

    SxList<SxString> list; list << ...;
    typedef SxList<SxString>::Iterator Iterator;
    sx::foreach (list.begin(), list.end(), [](Iterator &it) {
       cout << *it << endl;
    });

    or simpler

    list.foreach ([](Iterator &it) { cout << *it << endl; });


    \author Sixten Boeck, boeck@dacs-labs.com */

// see also
//    SxIterator.h:  SX_ITERATOR_LAMBDAS
//                   SX_CONST_ITERATOR_LAMBDAS
//    SxContainer.h  SX_CONTAINER

namespace sx
{
   // Enum for iterator direction
   enum Direction {Undefined, Forward, Backward, Both};

   enum ItCopyMode {NoCopy=0x01, CopyItMeta=0x02, CopyItData=0x04, CopyAll=0x08};

   // iterate over all elements
   template<class Iterator, class Function>
   void foreach (Iterator it, Iterator end, Function fn)
   {
      bool forward = it.isForward ();
      if (forward) {
         for (/* empty */; it != end; ++it)  {
            fn (it);
         }
      } else {
         for (/* empty */; it != end; --it)  {
            fn (it);
         }
      }
   }

   // find first element match and exit, return it.end() if not found
   template<class Iterator, class Elem>
   Iterator find (Iterator it, Iterator end, Elem elem)
   {
      Iterator res = end;
      bool forward = it.isForward ();
      if (forward) {
         for (/* empty */; it != end; ++it)  {
            if ( *it == elem )  {  res = it; break; }
         }
      } else {
         for (/* empty */; it != end; --it)  {
            if ( *it == elem )  {  res = it; break; }
         }
      }
      return res;
   }

   // find first condition match and exit, return it.end() if not found
   template<class Iterator, class Cond>
   Iterator findCond (Iterator it, Iterator end, Cond cond)
   {
      Iterator res = end;
      bool forward = it.isForward ();
      if (forward) {
         for (/* empty */; it != end; ++it)  {
            if ( cond(it) )  {  res = it; break; }
         }
      } else {
         for (/* empty */; it != end; --it)  {
            if ( cond(it) )  {  res = it; break; }
         }
      }
      return res;
   }

   // findAll : declared in SxList.h

   // sort using bubble sort [it, lastIt]
   template<class Iterator>
   void sort (Iterator it, Iterator lastIt,
              SxCBoundPtr<int, const decltype(Iterator())&,
                               const decltype(Iterator())&> comp)
   {
      SX_CHECK(it.isValid ());
      SX_CHECK(lastIt.isValid ());// should point to a valid element

      if (it == lastIt) return;

      Iterator tmpIt = it, tmpIt2 = it, itPrev = it;
      bool forward = it.isForward ();

      if (forward)  ++lastIt;
      else          --lastIt;

      while (tmpIt2 != lastIt) {
         SX_CHECK(tmpIt2.isValid ());
         tmpIt = it;
         itPrev = it;
         if (forward)  ++tmpIt;
         else          --tmpIt;
         SX_CHECK(tmpIt.isValid ());
         while (tmpIt != lastIt) {
            SX_CHECK(tmpIt.isValid ());
            SX_CHECK(itPrev.isValid ());

            if (comp(tmpIt,itPrev) == -1) {
               auto val = *tmpIt;
               *tmpIt = *itPrev;
               *itPrev = val;
            }
            if (forward) {
               ++tmpIt; ++itPrev;
            } else {
               --tmpIt; --itPrev;
            } 
         }
         if (forward)  ++tmpIt2;
         else          --tmpIt2;
      }
   }

   template<class Iterator>
   void sort (Iterator it, Iterator lastIt)
   {
      sx::sort (it, lastIt, [](const Iterator &itA, const Iterator &itB)->int
               { if (*itA < *itB)        return -1;
                 else if (*itA == *itB)  return  0;
                 else                    return  1;
               });

   }

   // sort using quick sort [it, lastIt]
   template<class Iterator>
   void qsort (Iterator it, Iterator lastIt,
               SxCBoundPtr<int, const decltype(Iterator())&,
                                const decltype(Iterator())&> comp)
   {
      SX_CHECK(it.isValid ());
      SX_CHECK(lastIt.isValid ());// should point to a valid element

      if (it == lastIt) return;

      bool forward = it.isForward ();
      Iterator pivotIt;
      Iterator leftIt = it;
      if (forward)  ++leftIt; 
      else          --leftIt;
      SX_CHECK(leftIt.isValid ());
      Iterator rightIt = lastIt;
      auto pivotVal = *it;

      while (leftIt != rightIt) {
         while ((leftIt != rightIt) && (comp(leftIt,it) != 1)) {
            if (forward)  ++leftIt;
            else          --leftIt;
            SX_CHECK(leftIt.isValid ());
         }
         while ((leftIt != rightIt) && (comp(rightIt,it) == 1)) {
            if (forward)  --rightIt;
            else          ++rightIt;
            SX_CHECK(rightIt.isValid ());
         } 
         if (leftIt == rightIt) break;
         auto val = *leftIt;
         *leftIt = *rightIt;
         *rightIt = val;
      }

      if (comp(leftIt,it) == 1) {
         if (forward)  --leftIt;
         else          ++leftIt;
      } 
      *it = *leftIt;
      *leftIt = pivotVal;
      pivotIt = leftIt;

      Iterator pivPrevIt = pivotIt;
      Iterator pivNextIt = pivotIt;
      if (forward)  {--pivPrevIt; ++pivNextIt;}
      else          {++pivPrevIt; --pivNextIt;}

      if (it != pivotIt)  qsort (it, pivPrevIt, comp);
      if (lastIt != pivotIt)  qsort (pivNextIt, lastIt, comp);

   }

   template<class Iterator>
   void qsort (Iterator it, Iterator lastIt)
   {
      sx::qsort<Iterator> (it, lastIt, [](const Iterator &itA, const Iterator &itB)->int
                { if (*itA < *itB)        return -1;
                  else if (*itA == *itB)  return  0;
                  else                    return  1;
                });
   }

}

#endif /* _SX_ALG_H_ */
