#include <SxCLI.h>
#include <SxJSONParser.h>
#include <SxSchema.h>
#include <SxGQuery.h>
#include <SxTimer.h>

typedef typename SxGQExprBase::SelSet SelSet;
using namespace sx;


void readAtoms () {

   SxJSONParser demo;
   SxJSONParser schemaParser;
   try {
      demo.setSearchPath (".");
      schemaParser.setSearchPath (".");
      SxTimer t(1);
      t.start ();

      demo.readFile ("atoms.dat");

      t.stop ();

      cout << "readfile time:" << t.getTime () << endl;

      t.reset ();
      t.start ();
      schemaParser.readFile ("atomSchema.dat");


      SxSchema schema (schemaParser.getAst ());

      schema.convert ();
      t.stop ();
      std::cout << "schema conversion time: " << t.getTime () << std::endl;

      t.reset ();
      t.start ();

      std::cout << "schema res: "
                << schema.validate (demo.getAst ())
                << endl;
      t.stop ();
      std::cout << "validation time: " << t.getTime () << std::endl;

      SxPtr<SxGraph<SxGProps> > ast = demo.getAst ();
      auto it = ast->begin (0);
      ++it; // skip dummy root

      // 1.6 secs, 100k
      // it is slow because we cannot further limit the search space
      // less than depth 2. query will search for specie in all atoms

      /*size_t c = 0;
      SxGQuery q = (N("sxKey") == "species") ^ (N("sxKey") == "atom");
      t.reset ();
      t.start ();
      SelSet sels = q.matchAll (ast, it.neighbors (2));
      cout << "sel size: " << sels->getSize () << endl;
      for (auto selIt = sels->begin (); selIt != sels->end (); ++selIt)
      {
         auto it = ast->begin ((**selIt)(1)).out (0);
         SxList<SxVariant>::Iterator elemIt = it->getProperty ("sxValue").begin ();
         for (; elemIt.isValid (); ++elemIt) {
            elemIt->toDouble ();
         }
         ++c;
      }
      t.stop ();
      std::cout << "number of atoms: " << c++ << std::endl;
      std::cout << "atoms query time: " << t.getTime () << std::endl;*/
   

      // 1.8 secs, 100k
      // Searching for specie is fast if search is limited to single
      // level. But then again we search for atoms and coords, where
      // large search space causes a lot of failed matches.

      /*SxGQuery q = (N("sxKey") == "species");

      t.reset ();
      t.start ();
      SelSet sels = q.matchAll (ast, it.neighbors ());
      std::cout << "number of species: " << sels->getSize () << std::endl;

      size_t c = 0;
      SxList<double> lst;
      SxGQuery q2 = (N("sxKey") == "atom") ^ (N("sxKey") == "coords");
      for (auto selIt = sels->begin (); selIt != sels->end (); ++selIt)
      {
         SelSet aSels = q2.matchAll (ast, ast->begin ((**selIt)(0)).neighbors ());

         for (auto aselIt = aSels->begin (); aselIt != aSels->end (); ++aselIt)
         {
            auto it = ast->begin ((**aselIt)(1));
            SxList<SxVariant>::Iterator elemIt = it->getProperty ("sxValue").begin ();
            for (; elemIt.isValid (); ++elemIt) {
               elemIt->toDouble ();
            }
            c++;
         }

      }

      t.stop ();
      std::cout << "number of atoms: " << c << std::endl;
      std::cout << "atom coords: " << lst << std::endl;
      std::cout << "atoms query time: " << t.getTime () << std::endl;*/

      
      // 0.9 secs, 100k
      // faster then others because we first search all the
      // species in targeted area and collect the atoms only
      // from those individual species. hence minimizing the 
      // number of nodes to be checked
      SxGQuery q = (N("sxKey") == "species");

      t.reset ();
      t.start ();
      SelSet sels = q.matchAll (ast, it.neighbors ());
      std::cout << "number of species: " << sels->getSize () << std::endl;
      t.stop ();
      std::cout << "query time: " << t.getTime () << std::endl;

      size_t c = 0;
      t.reset ();
      t.start ();
      SxList<double> lst;
      SxGQuery q2 = (N("sxKey") == "atom");
      for (auto selIt = sels->begin (); selIt != sels->end (); ++selIt)
      {
         SelSet aSels = q2.matchAll (ast, ast->begin ((**selIt)(0)).neighbors ());

         for (auto aselIt = aSels->begin (); aselIt != aSels->end (); ++aselIt)
         {
            auto coordIt = ast->begin ((**aselIt)(0)).out (0);
            SxList<SxVariant>::Iterator elemIt = coordIt->getProperty ("sxValue").begin ();
            for (; elemIt.isValid (); ++elemIt) {
               elemIt->toDouble ();
            }
            c++;
         }

      }
      t.stop ();
      std::cout << "number of atoms: " << c << std::endl;
      std::cout << "atom coords: " << lst << std::endl;
      std::cout << "atoms query time: " << t.getTime () << std::endl;

   } catch (SxException e) {
      e.print ();
   }

}


int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   cli.finalize ();


   SxJSONParser demo;
   SxJSONParser demo2;
   try {
      demo.setSearchPath ("."); 
      demo2.setSearchPath ("."); 

      // -- parse invalid comment
      //demo.readString ("*/abc");

      // -- parse from file
      //SxString s = SxString::fromUtf8 (SxString::read ("abc.dat").getElems());

      // -- parse UTF-8 string and test file include
      //SxString s = SxString::fromUtf8 (u8"\n#include<./abc.dat> \"παράδειγμα\":\"aad\\n\\nds\xADfe\"");

      // -- parse multilevel comments
      //SxString s = SxString::fromUtf8 (u8"/*co/*mme*/nt*/\"παράδειγμα\":\"aad\\n\\nds\xADfe\"");

      // -- parse UTF-8 string with comments and new lines

      //SxString s = SxString::unicodeFromUtf8 (u8"/*comment*/\"παράδειγμα\":\"ตัว\\n\\nอย่าง\"");
      //SxString s = SxString::unicodeFromUtf8 (u8"/*comment*/\"παράδειγμα\":\"ตัวอย่าง\"");
      //SxString s2 = SxString::unicodeFromUtf8 (u8"ตัวอย่าง");
      //demo.readString (s);
      

      // -- example json file
      demo.readFile ("ary.dat");

      // -- example json schema file
      demo2.readFile ("arySch.dat");
      SxSchema schema (demo2.getAst ());

      // -- optimize schema tree/graph
      schema.convert ();

      std::cout << "Schema res: " << schema.validate (demo.getAst ()) << std::endl;

   } catch (SxException e) {
      e.print ();
   }

   return 0;
}
