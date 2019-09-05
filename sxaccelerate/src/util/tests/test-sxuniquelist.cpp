
//  SxUniqueList<> test

// --- SPHInX
#include <SxString.h>
#include <SxCLI.h>
#include <SxUniqueList.h>

// --- errors in release and debug mode
#ifdef NDEBUG
#   define SX_TEST(expr)                                                     \
           if ( !(expr) )  {                                                 \
              std::cout << std::endl << "ASSERTATION FAILED in "             \
                        << __FILE__ << ", line " << __LINE__ << "!\n";       \
              std::cout << #expr "\n\n";                                     \
              SX_EXIT;                                                       \
           }
#else
#   define SX_TEST(p) SX_CHECK (p)
#endif /* NDEBUG */

template<class H>
void testInsert ()
{
   SxUniqueList<int,H> list;
   
   SX_TEST (list.getSize() == 0);
   SX_TEST (list.findPos (1) < 0);

   list << 1;
   SX_TEST (list.getSize() == 1);
   SX_TEST (list.findPos (1) == 0);
   SX_TEST (list.contains (1));
   SX_TEST (list.contains (100) == false);
   SX_TEST (list == SxList<int>() << 1);
   
   list << 1 << 1 << 2;
   SX_TEST (list.getSize() == 2);
   SX_TEST (list.findPos (1) == 0);
   SX_TEST (list.findPos (2) == 1);
   SX_TEST (list.contains (1));
   SX_TEST (list.contains (2));
   SX_TEST (list.contains (100) == false);
   SX_TEST (list == SxList<int>() << 1 << 2);
   
   SxUniqueList<int,H> list2 = list;
   SX_TEST (list2.getSize() == 2);
   SX_TEST (list2.contains (1));
   SX_TEST (list2.contains (2));
   SX_TEST (list2 == SxList<int>() << 1 << 2);

   SxUniqueList<int,H> list3(list);
   SX_TEST (list3.getSize() == 2);
   SX_TEST (list3.contains (1));
   SX_TEST (list3.contains (2));
   SX_TEST (list3 == SxList<int>() << 1 << 2);
   
   SxUniqueList<int,H> list4;
   list4 << list;
   SX_TEST (list4.getSize() == 2);
   SX_TEST (list4.contains (1));
   SX_TEST (list4.contains (2));
   SX_TEST (list4 == SxList<int>() << 1 << 2);
   list4 << list;
   SX_TEST (list4 == SxList<int>() << 1 << 2);
   list4 = (SxList<int>() << 2 << 3);
   SX_TEST (list4.getSize() == 2);
   SX_TEST (list4.contains (1) == false);
   SX_TEST (list4.contains (2));
   SX_TEST (list4.contains (3));
   SX_TEST (list4 == SxList<int>() << 2 << 3);

   list.prepend (0);
   SX_TEST (list.getSize() == 3);
   SX_TEST (list.findPos (0) == 0);
   SX_TEST (list.findPos (1) == 1);
   SX_TEST (list.findPos (2) == 2);
   SX_TEST (list.contains (0));
   SX_TEST (list.contains (1));
   SX_TEST (list.contains (2));
   SX_TEST (list.contains (100) == false);
   SX_TEST (list == SxList<int>() << 0 << 1 << 2);
   list.prepend (0);
   SX_TEST (list == SxList<int>() << 0 << 1 << 2);
   
   list.append (5);
   SX_TEST (list.getSize() == 4);
   SX_TEST (list.contains (0));
   SX_TEST (list.contains (1));
   SX_TEST (list.contains (2));
   SX_TEST (list.contains (5));
   SX_TEST (list.contains (100) == false);
   SX_TEST (list == SxList<int>() << 0 << 1 << 2 << 5);
   list.append (5);
   SX_TEST (list == SxList<int>() << 0 << 1 << 2 << 5);

   list.insert (3, 4);
   SX_TEST (list.getSize() == 5);
   SX_TEST (list.contains (0));
   SX_TEST (list.contains (1));
   SX_TEST (list.contains (2));
   SX_TEST (list.contains (4));
   SX_TEST (list.contains (5));
   SX_TEST (list.contains (100) == false);
   SX_TEST (list == SxList<int>() << 0 << 1 << 2 << 4 << 5);
   list.insert (3, 4);
   SX_TEST (list == SxList<int>() << 0 << 1 << 2 << 4 << 5);

   // --- 100, resize hash table (keep=true)
   SxUniqueList<int,H> list5;
   SxUniqueList<SxString,H> list5str;
   for (int i=0; i < 100; ++i) {
      list5 << i;
      list5str << SxString(i);
   }
   for (int i=0; i < 100; ++i) {
      SX_TEST (list5.contains (i));
      SX_TEST (list5.findPos (i) == i);
      SX_TEST (list5str.contains (i));
      SX_TEST (list5str.findPos (i) == i);
   }

   //SX_TEST ((list |= 1) == 0);
   //SX_TEST ((list |= 2) == 1);
   //SX_TEST ((list |= 5) == 4);
}

template<class H>
void testRemove ()
{
   SxUniqueList<int,H> list;
   list << 0 << 1 << 2 << 3;

   list.removeFirst ();
   SX_TEST (list.getSize() == 3);
   SX_TEST (list.contains (0) == false);
   SX_TEST (list.contains (1));
   SX_TEST (list.contains (2));
   SX_TEST (list.contains (3));
   SX_TEST (list.contains (100) == false);
   SX_TEST (list == SxList<int>() << 1 << 2 << 3);
   
   list.removeLast ();
   SX_TEST (list.getSize() == 2);
   SX_TEST (list.contains (1));
   SX_TEST (list.contains (2));
   SX_TEST (list.contains (3) == false);
   SX_TEST (list.contains (100) == false);
   SX_TEST (list == SxList<int>() << 1 << 2);

   list.removeElement (2);
   SX_TEST (list.getSize() == 1);
   SX_TEST (list.contains (1));
   SX_TEST (list.contains (2) == false);
   SX_TEST (list.contains (100) == false);
   SX_TEST (list == SxList<int>() << 1);

   list.removeAll ();
   SX_TEST (list.getSize() == 0);
   SX_TEST (list.contains (1) == false);
   SX_TEST (list.contains (100) == false);
   SX_TEST (list == SxList<int>());
   
   list << (SxList<int>() << 1 << 2);
   SX_TEST (list == SxList<int>() << 1 << 2);

   list << (SxList<int>() << 2 << 3 << 4);
   SX_TEST (list == SxList<int>() << 1 << 2 << 3 << 4);
}

template<class H>
SxUniqueList<SxString,H> foo (const SxUniqueList<SxString,H> &list)
{
   SxUniqueList<SxString,H> res;
   res << list;
   return res;
}

template<class H>
void testCopy ()
{
   SxUniqueList<SxString,H> list;
   
   list << "a" << "b";
   list = foo<H> (list);
   list.removeElement ("a");
   
   SX_TEST (list == SxList<SxString>() << "b");
}

template<class H>
void test ()
{
   testInsert<H> ();
   testRemove<H> ();
   testCopy<H> ();
}

int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   cli.finalize ();
   
   test<SxHashFunction> ();
   test<SxNull> ();
               
   printf ("ok\n");
   
   return 0;
}
