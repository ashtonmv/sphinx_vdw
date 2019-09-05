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

#include <SxCLI.h>
#include <SxMetadata.h>
#include <SxTimer.h>

class Foo {};

ostream &operator<< (ostream& out, const class Foo&)
{
   out << "foo";
   return out;
}

enum MyTime { GetMeta } ;

SX_REGISTER_TIMERS(MyTime)
{
   regTimer (GetMeta, "get meta");
}

enum Set1 { Lala = 0, Huhu = 1 };
enum Set2 { Mama = 0, Meme = 1 };

int main (int argc, char **argv)
{
   // --- command line parsing
   SxCLI cli(argc, argv);
   int Nj = cli.option ("-N", "int", "number of repetitions").toInt (100);
   cli.finalize ();

   SxMetadata m;


   m.attach ("A", int(3));
   m.attach ("B", 4.3);
   m.attach ("C", SxString ("hallo") );
   m.attach (Foo ());

   for (int i = 0; i < 1000; ++i)  {
      m.attach (i, sqrt(double(i)));
   }

   m.attach (Lala, SxString ("lala"));
   m.attach (Mama, SxString ("mama"));
   m.attach (Meme, 1.);

   cout << "A=" << m.get<int> ("A") << endl;
   cout << "B=" << m.get<double> ("B") << endl;
   cout << "C=" << m.get<SxString> ("C") << endl;
   cout << "foo: " << m.get<Foo> () << endl;
   cout << "lala=" << Lala << ": " << m.get<SxString> (Lala) << endl;
   cout << "Mama=" << Mama << ": " << m.get<SxString> (Mama) << endl;

   SxString &refC = m.get ("C");
   refC = "tschoe";
   cout << "C=" << m.get<SxString> ("C") << endl;

   m.remove<Foo> ();
   m.remove (Mama);

   SX_START_TIMER (GetMeta);
   SxString key = "B";
   for (int j = 0; j < Nj; ++j)  {
      double sum = 0.;
      for (int i = 0; i < 1000; ++i)  {
         //sum += m.get<double> (i);
         sum += m.get<double> (key);
         //sum += m.get<double> (Meme);
      }
   }
   SX_STOP_TIMER (GetMeta);
   printTiming ();

}
