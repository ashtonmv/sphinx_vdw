
//  SxMap<> and SxMultiMap<> test
//  
//  This test is automatic. It will terminate with SX_EXIT
//  if the results do not match with expected values.
//


// --- SPHInX
#include <SxString.h>
#include <SxCLI.h>
#include <SxTime.h>
#include <SxMap.h>
#include <SxMultiMap.h>

#ifdef HAVE_TR1
#  include <tr1/unordered_map>
#endif

#ifdef NDEBUG
#   define SX_TEST(p)                       \
              if (!(p))  {                  \
                 printf ("fail: %s\n", #p); \
                 return 1;                  \
              }
#else
#   define SX_TEST(p) SX_CHECK(p);
#endif /* NDEBUG */

template<class T>
bool testN (ssize_t n)
{
   T map;
   
   ssize_t i;
   ssize_t count;
   ssize_t idx;
   typename T::Key tmp;
   
   SxArray<int> dataInt(n);
   SxArray<int> sequence(n);
   
   for (i = 0; i < n; i++)  {
      dataInt(i) = (typename T::Key)i;
      sequence(i) = (typename T::Key)i;
   }
   for (i = 0; i < n; i++)  {
      idx = rand() % n;
      tmp = dataInt(i);
      dataInt(i) = dataInt(idx);
      dataInt(idx) = tmp;
   }
   
   SX_DBG_MSG ("insert all keys");
   for (i = 0; i < n; i++)  {
      map(dataInt(i)) = (typename T::Key)i;
   }
   SX_TEST (map.getSize() == n);
   SX_TEST (SxArray<typename T::Key>(map.getKeys ()) == dataInt);
   SX_TEST (SxArray<typename T::Value>(map.getValues ()) == sequence);

   SX_DBG_MSG ("find all keys");
   count = 0;
   for (i = 0; i < n; i++)  {
      if (map.containsKey (dataInt(i)))  {
         count++;
      }
   }
   SX_TEST (count == n);

   count = 0;
   for (i = 0; i < n; i++)  {
      if ((idx = map.findKey (dataInt(i))) >= 0)  {
         //SX_TEST (map.getValue(idx) == map(dataInt(i)));
         count++;
      }
   }
   SX_TEST (count == n);
      
   SX_DBG_MSG ("remove all keys");
   SxList<typename T::Key> keys = map.getKeys ();
   //SxList<typename T::Value> values = map.getValues ();
   count = 0;
   for (i = 0; i < n; i++)  {
      keys.removeElement (dataInt(i));
      //values.removeElement (map(dataInt(i)));
      if (map.removeKey (dataInt(i)))  {
         count++;
      }
      if (map.containsKey (dataInt(i)))  {
         cout << "ERROR map still contains the key, remove at " << i << endl;
         return false;
      }
      SX_TEST (map.getKeys () == keys);
      //SX_TEST (map.getValues () == values);
      ssize_t j;
      for (j = i + 1; j < n; j++)  {
         if (!map.containsKey (dataInt(j)))  {
            cout << "ERROR lost key, remove at " << i << endl;
            return false;
         }
      }
   }
   SX_TEST (count == n);
   SX_TEST (map.getSize () == 0);
   SX_TEST (map.getKeys () == SxList<int> ());
   SX_TEST (map.getValues () == SxList<int> ());
      
   return true;
}


template<class T>
bool testBase ()
{
   T map;
   
   SX_DBG_MSG ("test empty map");
   SX_TEST (map.getSize() == 0);
   SX_TEST (map.containsKey (1) == false);
   SX_TEST (map.containsValue (1) == false);
   SX_TEST (map.getKeys () == SxList<int> ());
   SX_TEST (map.getValues () == SxList<int> ());
   SX_TEST (map.begin () == map.end ());
   SX_TEST (map.removeKey (1) == false);
   SX_TEST (map.removeValue (1) == false);
   
   SX_DBG_MSG ("insert one key");
   map(0) = 5;
   SX_TEST (map.getSize() == 1);
   SX_TEST (map.containsKey (0));
   SX_TEST (map.containsKey (1) == false);
   SX_TEST (map.containsValue (5));
   SX_TEST (map.containsValue (4) == false);
   SX_TEST (map.findKey (0) >= 0);
   SX_TEST (map.findKey (1) < 0);
   SX_TEST (map.getKeys () == SxList<int>() << 0);
   SX_TEST (map.getValues () == SxList<int>() << 5);
   SX_TEST (map.begin () != map.end ());
   SX_TEST (map.removeKey (1) == false);
   SX_TEST (map.removeValue (1) == false);
   
   SX_DBG_MSG ("remove one key");
   SX_TEST (map.removeKey (0));
   SX_TEST (map.getSize() == 0);
   SX_TEST (map.containsKey (1) == false);
   SX_TEST (map.containsValue (1) == false);
   SX_TEST (map.getKeys () == SxList<int> ());
   SX_TEST (map.getValues () == SxList<int> ());
   SX_TEST (map.begin () == map.end ());
   SX_TEST (map.removeKey (1) == false);
   SX_TEST (map.removeValue (1) == false);
   
   SX_DBG_MSG ("remove one value");
   map(0) = 5;
   SX_TEST (map.removeValue (5));
   SX_TEST (map.getSize() == 0);
   SX_TEST (map.containsKey (1) == false);
   SX_TEST (map.containsValue (1) == false);
   SX_TEST (map.getKeys () == SxList<int> ());
   SX_TEST (map.getValues () == SxList<int> ());
   SX_TEST (map.begin () == map.end ());
   SX_TEST (map.removeKey (1) == false);
   SX_TEST (map.removeValue (1) == false);
   
   SX_DBG_MSG ("remove and insert");
   map(0) = 10;
   map(1) = 11;
   SX_TEST (map.getSize() == 2);
   SX_TEST (map.containsKey (0));
   SX_TEST (map.containsKey (1));
   SX_TEST (map.containsKey (2) == false);
   SX_TEST (map.containsValue (10));
   SX_TEST (map.containsValue (11));
   SX_TEST (map.containsValue (12) == false);
   SX_TEST (map.findKey (0) >= 0);
   SX_TEST (map.findKey (1) >= 0);
   SX_TEST (map.findKey (2) < 0);
   SX_TEST (map.getKeys () == SxList<int>() << 0 << 1);
   SX_TEST (map.getValues () == SxList<int>() << 10 << 11);
   SX_TEST (map.begin () != map.end ());
   SX_TEST (map.removeKey (2) == false);
   SX_TEST (map.removeValue (12) == false);
   SX_TEST (map.removeKey (0));
   map(2) = 12;
   SX_TEST (map.getSize() == 2);
   SX_TEST (map.containsKey (0) == false);
   SX_TEST (map.containsKey (1));
   SX_TEST (map.containsKey (2));
   SX_TEST (map.containsValue (10) == false);
   SX_TEST (map.containsValue (11));
   SX_TEST (map.containsValue (12));
   SX_TEST (map.getKeys () == SxList<int>() << 1 << 2);
   SX_TEST (map.getValues () == SxList<int>() << 11 << 12);
   SX_TEST (map.begin () != map.end ());
   SX_TEST (map.removeKey (0) == false);
   SX_TEST (map.removeValue (10) == false);
   map.removeAll ();
   SX_TEST (map.getSize() == 0);
   SX_TEST (map.containsKey (0) == false);
   SX_TEST (map.containsKey (1) == false);
   SX_TEST (map.containsKey (2) == false);
   SX_TEST (map.containsValue (10) == false);
   SX_TEST (map.containsValue (11) == false);
   SX_TEST (map.containsValue (12) == false);
   SX_TEST (map.getKeys () == SxList<int>());
   SX_TEST (map.getValues () == SxList<int>());
   
   return true;
}

template<class T>
bool testCopy ()
{
   T map;
   
   map(0) = 10;
   map(1) = 11;
   
   T map2 = map;
   SX_TEST (map2.getSize() == 2);
   SX_TEST (map2.containsKey(0));
   SX_TEST (map2.containsKey(1));
   SX_TEST (map2(0) == 10);
   SX_TEST (map2(1) == 11);
   SX_TEST (map2.getKeys () == SxList<int>() << 0 << 1);
   SX_TEST (map2.getValues () == SxList<int>() << 10 << 11);
   
   T map3(map);
   SX_TEST (map3.getSize() == 2);
   SX_TEST (map3.containsKey(0));
   SX_TEST (map3.containsKey(1));
   SX_TEST (map3(0) == 10);
   SX_TEST (map3(1) == 11);
   SX_TEST (map3.getKeys () == SxList<int>() << 0 << 1);
   SX_TEST (map3.getValues () == SxList<int>() << 10 << 11);
      
   return true;
}

template<class T>
bool test ()
{
   if (!testBase<T>()) return false;
   if (!testCopy<T>()) return false;
   if (!testN<T>(16))  return false;
         
   return true;
}

bool testMultiMap ()
{
   SxMultiMap<int,int> map;
   
   map.append (0, 1);
   SX_TEST (map.first () == 1);
   SX_TEST (map.last () == 1);
   map.append (0, 2);
   SX_TEST (map.first () == 1);
   SX_TEST (map.last () == 2);
   
   SX_TEST (map.getSize() == 2);
   SX_TEST (map.containsValue (3) == false);
   SX_TEST (map.containsValue (1));
   SX_TEST (map.containsValue (2));
   SX_TEST (map.containsKey (1) == false);
   SX_TEST (map.containsKey (0));
   SX_TEST (map.getValues(0) == SxList<int>() << 1 << 2);
   //SX_TEST (map.getKeys ()   == SxList<int>() << 0 << 0);
   SX_TEST (map.getKeys ()   == SxList<int>() << 0);
   SX_TEST (map.getValues () == SxList<int>() << 1 << 2);
   
   SX_TEST (map.removeValue (1));
   SX_TEST (map.getSize() == 1);
   SX_TEST (map.containsValue (1) == false);
   SX_TEST (map.containsValue (2));
   SX_TEST (map.containsKey (1) == false);
   SX_TEST (map.containsKey (0));
   SX_TEST (map.getValues(0) == SxList<int>() << 2);
   SX_TEST (map.getKeys ()   == SxList<int>() << 0);
   SX_TEST (map.getValues () == SxList<int>() << 2);
   
   return true;
}

template<class T>
void benchInt (ssize_t n, const int *data, double *t1, double *t2, double *t3)
{
   T map;
   double t;
   ssize_t i;
   
   t = SxTime::getRealTime ();
   for (i = 0; i < n; i++)  {
      map(data[i]) = (typename T::Key)i;
   }
   *t1 = SxTime::getRealTime () - t;
   
   ssize_t count = 0;
   t = SxTime::getRealTime ();
   for (i = 0; i < n; i++)  {
      if (map.containsKey (data[i]))  {
         count++;
      }
   }
   *t2 = SxTime::getRealTime () - t;
   std::cout << count << std::endl;

   t = SxTime::getRealTime ();
   for (i = 0; i < n; i++)  {
      map.removeKey (data[i]);
   }
   *t3 = SxTime::getRealTime () - t;
}

template<class T>
void benchIntStd (ssize_t n, const int *data, double *t1,double *t2,double *t3)
{
   T map;
   double t;
   ssize_t i;
   
   t = SxTime::getRealTime ();
   for (i = 0; i < n; i++)  {
      map[data[i]] = (int)i;
   }
   *t1 = SxTime::getRealTime () - t;
   
   ssize_t count = 0;
   t = SxTime::getRealTime ();
   for (i = 0; i < n; i++)  {
      if (map.find (data[i]) != map.end ())  {
         count++;
      }
   }
   *t2 = SxTime::getRealTime () - t;
   std::cout << count << std::endl;

   t = SxTime::getRealTime ();
   for (i = 0; i < n; i++)  {
      map.erase (data[i]);
   }
   *t3 = SxTime::getRealTime () - t;
}


void benchN (ssize_t n)
{
   ssize_t i;
   ssize_t idx;
   int tmp;
   
   SxArray<int> dataInt(n);
   
   for (i = 0; i < n; i++)  {
      dataInt(i) = (int)i;
   }
   for (i = 0; i < n; i++)  {
      idx = rand() % n;
      tmp = dataInt(i);
      dataInt(i) = dataInt(idx);
      dataInt(idx) = tmp;
   }
   
   SxArray<SxString> id;
   id.append ("Map");
   id.append ("MultiMap");
#ifdef HAVE_TR1
   id.append ("std::unordered_map");
   //id(3) = "Map (nohash)";
#endif
   
   SxArray<double> t1(id.getSize ());
   SxArray<double> t2(id.getSize ());
   SxArray<double> t3(id.getSize ());
   
   benchInt<SxMap<int,int> > (n, dataInt.elements, &t1(0), &t2(0), &t3(0));
   benchInt<SxMultiMap<int,int> > (n, dataInt.elements, &t1(1), &t2(1), &t3(1));

#ifdef HAVE_TR1
   benchIntStd<std::tr1::unordered_map<int,int> > (
      n, dataInt.elements, &t1(2), &t2(2), &t3(2));
#endif
//benchInt<SxMap<int,int,SxNull> > (n, dataInt.elements, &t1(3),&t2(3),&t3(3));
   
   SxArray<ssize_t> sortIdx;
   SxArray<SxString> sortId;
   
   sortIdx = t1.getSortIdx ();
   sortId = id;
   sortId.sortByIdx (sortIdx);
   t1.sortByIdx (sortIdx);
   
   std::cout << "\ninsert [ms]" << std::endl;
   for (i = 0; i < id.getSize (); i++)  {
      std::cout << "   " << t1(i)*1.e3 << " : " << sortId(i) << std::endl;
   }

   sortIdx = t2.getSortIdx ();
   sortId = id;
   sortId.sortByIdx (sortIdx);
   t2.sortByIdx (sortIdx);
   
   std::cout << "\nfind [ms]" << std::endl;
   for (i = 0; i < id.getSize (); i++)  {
      std::cout << "   " << t2(i)*1.e3 << " : " << sortId(i) << std::endl;
   }

   sortIdx = t3.getSortIdx ();
   sortId = id;
   sortId.sortByIdx (sortIdx);
   t3.sortByIdx (sortIdx);
   
   std::cout << "\nremove [ms]" << std::endl;
   for (i = 0; i < id.getSize (); i++)  {
      std::cout << "   " << t3(i)*1.e3 << " : " << sortId(i) << std::endl;
   }
}


void bench (ssize_t n)
{
   benchN (n);
}


int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   ssize_t size  = cli.option ("-n", "int", "n elements").toInt (100000);
   cli.finalize ();
   
   if (!test<SxMap<int,int> > ())  {
      printf ("failed\n");
      return 1;
   }
   
   if (!test<SxMap<int,int,SxNull> > ())  {
      printf ("failed\n");
      return 1;
   }
   
   if (!test<SxMultiMap<int,int> > ())  {
      printf ("failed\n");
      return 1;
   }
   
   if (!testMultiMap ())  {
      printf ("failed\n");
      return 1;
   }
   
   printf ("ok\n");
   
   bench (size);
   
   return 0;
}
