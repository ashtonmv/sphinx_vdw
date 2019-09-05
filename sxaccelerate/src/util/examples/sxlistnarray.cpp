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

#include <SxUtil.h>
#include <SxList.h>
#include <SxArray.h>
#include <stdio.h>
#include <iostream>
#include <SxTimer.h>


template<class T> void listPrint (const SxList<T> &list);  


int main ()
{
   const size_t n = 10000;
   SxList<int> myList,listCopy;
   SxTimer timerBoard;

   timerBoard.init (1,true);

   timerBoard.start (0);
   
   // --- Initialisation of a list by adressing the elements directly.
   //     Acces of list elemets with "()"
   //     Very unperformant. Do not do this.
   myList.resize(n);
   for (size_t i = 0; i < n; i++)  {
      myList(i) = i;
   }
   timerBoard.stop (0);

   printf("for-loop initialization adress element directly of %.0e elements",\
         1.0*n);
   printf(" took %.2g seconds.\n", timerBoard.getTime (0));
   
   
   timerBoard.reset (0);
   myList.removeAll ();
   timerBoard.start (0);


   // --- Initialisation of a list with the append operator "<<". 
   //     Quite performant, but not as fast as the iterator
   
   for (size_t i = 0; i < 1000*n; i++)  {
      myList << i;
   }
   timerBoard.stop(0);
 
   printf("for-loop initialization using \"<<\" of %.0e elements took ",\
         1000.0*n);
   printf("%.2g seconds.\n", timerBoard.getTime (0)); 

   
   timerBoard.reset (0);
   myList.removeAll ();
   timerBoard.start (0);

   // --- The way how to setup a list: USE THE ITERATOR
   myList.resize(1000*n);
   SxList<int>::Iterator it;
   size_t data = 0;
   for (it = myList.begin (); it != myList.end (); it++)  {
      *it = data;
      data++;
   }
   

   timerBoard.stop (0);

   printf("iterator initialization of %.0e elements took %.2g seconds.\n",\
         1000.0*n, timerBoard.getTime (0)); 


   // --- Working with lists, SxList is designed for stacks
   //     so you can delete everywhere but just append 
  
   cout << endl << endl;
   myList.resize (10);
  
   cout << "This is the start list:\n";
   listPrint (myList);

   listCopy = myList;

   listCopy.removeFirst ();
   cout << "Deleting first element:\n";
   listPrint(listCopy);

   cout << "Deleting last element:\n";
   listCopy.removeLast();
   listPrint(listCopy);

   cout << "Deleting element \"6\":\n";
   listCopy.removeElement (6);
   listPrint(listCopy);
   
   cout << "Deleting third element:\n";
   // Remeber element counting starts with 0
   listCopy.remove (3-1);
   listPrint(listCopy);

   cout << "Appending \"666\":\n";
   listCopy.append (666); 
   listPrint(listCopy);

   myList.resize(1000*n);
   data = 0;
   for (it = myList.begin (); it != myList.end (); it++)  {
      *it = data;
      data++;
   }

   cout << endl << endl;
   // --- Now we copy the list in an array and compare the acces time.
   //     Acces of array elements with "()"

   SxArray<int> myArray (myList);


   timerBoard.reset(0);
   timerBoard.start(0);

   cout << "center element of list: " << myList(500*n-1) << endl;

   timerBoard.stop(0);

   cout << "printout took "<<  timerBoard.getTime(0) << " seconds." << endl;

   timerBoard.reset(0);
   timerBoard.start(0);

   cout << "center element of array: " << myArray(500*n-1) << endl;

   timerBoard.stop(0);

   cout << "printout took "<<  timerBoard.getTime(0) << " seconds."<< endl;

   return 0;
}

template <class T> void listPrint (const SxList<T> &list)  
{
   typename SxList<T>::ConstIterator it;
   for (it = list.begin (); it != list.end (); it++)  {
      cout << *it << " ";
   }
   cout << endl << endl;
}

