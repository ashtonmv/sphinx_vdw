
//  SxHashFunction test
//  

// --- SPHInX
#include <SxString.h>
#include <SxCLI.h>
#include <SxHashFunction.h>

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


bool test ()
{
   cout << "hash(10)=" << SxHashFunction::hash (10) << endl;
   
   // --- 3.141592653589793
   cout << "hash(pi)=" << SxHashFunction::hash (3.14159265) << endl;
   cout << "hash(pi)=" << SxHashFunction::hash (3.141592653) << endl;
   
   SxList<double> list = SxList<double>() << 3.14159265 << 3.141592653;
   cout << "hash(SxList<double>)=" << SxHashFunction::hash (list) << endl;

   SxArray<double> a = list;
   cout << "hash(SxArray<double>)=" << SxHashFunction::hash (a) << endl;
   
   SX_TEST (SxHashFunction::hash (list) == SxHashFunction::hash (a));
   
   size_t h = SxHashFunction::hash (3.14159265);
   cout << "hash(pi)=" << h << endl;
   h = SxHashFunction::combine (h, SxHashFunction::hash (3.14159265));
   cout << "combine(pi)=" << h << endl;

   SxString null;
   SxString empty = "";
   SX_TEST (SxHashFunction::hash (null) == SxHashFunction::hash (empty));
   
   return true;
}


int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   cli.finalize ();
   
   if (!test ())  {
      printf ("failed\n");
      return 1;
   }
         
   printf ("ok\n");
   
   return 0;
}
