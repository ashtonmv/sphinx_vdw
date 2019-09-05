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
#include <SxString.h>
#include <SxCLI.h>
#include <stdio.h>
#include <SxTimer.h>
#include <SxError.h>
#include <SxMarker.h>
#include <SxMemMonitor.h>
#include <SxCalc.h>

double mySqrt (double val)
{
   SX_CHECK (val >= 0, val);
   return sqrt(val);
}


class FuncLib : public SxMarker
{
   public:
      FuncLib () : SxMarker () { }
      virtual ~FuncLib () { /* empty */ }

      virtual SxString function (const SxString &funcId, 
                                 const SxString &in_) const
      {
         if (funcId == "sqrt")  {
            bool err = false;
            SxString in;
            if (markers.containsKey (in_))  in = markers(in_);
            else                            in = in_;

            double x = in.toDouble (&err);
            if (err)  return "sqrt(NaN)";
            return SxString ( ::sqrt(x) );
         } else if (funcId == "upper")  {
            if (markers.containsKey(in_))  
               return markers(in_).toUpper();
            else
               return in_.toUpper ();
         } else if (funcId == "lower") {
            if (markers.containsKey(in_))  
               return markers(in_).toLower();
            else
               return in_.toLower ();
         } else
            return "";
      }
};



/** \example strings.cpp

  In this example demonstrates the usage of the SFHIngX string manipulation 
  class.

  \brief  Basic usage of strings
  \author Sixten Boeck
  */
int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   cli.finalize ();

   SxString s = "123";
   cout << "LONG = " << s.toInt64 () << endl;
   return 0;

   cout << "Sqrt(-5) = " << mySqrt(-5) << endl;



   SxList<int> list;
   list << 10 << 20 << 30;
   cout << list << endl;


   return 0;

/*


//   SxRedirect tee (std::cout, "abc.log");
//   SxLog log (&tee);
////   log.printf ("%s %12.6f\n", "hi", 1./3.);
// cout << SxString::sprintf ("%s, %12.6f\n", "ab", 1./3.);
// sxprintf ("%s, %12.6f\n", "abc", 1./3.);
// log.printf ("%s, %12.6f\n", "def", 1./3.);
//
//     cout << "1: this is a C++ stream buffer" << endl;
//     cout << "2: this is a C++ stream buffer\n";
//     cout << "3: this is a C++ stream buffer\n";
//
//
//   return 0;

   SxString text = "  This    text 		contains   some whitespaces.    ";
   cout << "text:              '" << text                      << "'" << endl;
   cout << "trim:              '" << text.trim()               << "'" << endl;
   cout << "stripWhiteSpace:   '" << text.stripWhiteSpace()    << "'" << endl;
   cout << "simplifyWhiteSpace: '" << text.simplifyWhiteSpace() << "'" << endl;
   cout << "removeWhiteSpace:  '" << text.removeWhiteSpace()   << "'" << endl;

   text = text.simplifyWhiteSpace ();
   cout << "left:              '" << text.left ("contains")    << "'" << endl;
   cout << "right:             '" << text.right ("contains")    << "'" << endl;
   cout << "substitute:        '" << text.substitute("text", "string") 
                                  << "'" << endl;

   // --- regular expressions
   SxString date = "Mon Mar 28 11:08:31 CEST 2006";
   cout << "date:              '" << date << "'" << endl;
   SxList<SxString> res;
   try  { res = date.regexp ("^(\\w+)\\s+(\\w+)\\s+(\\d+)\\s+([0-9]+):([0-9]+):([0-9]+)\\s+\\w+\\s+([0-9]+)"); }
   catch (SxException e)  { e.print (); SX_EXIT; }
   SX_CHECK_VAR (res.getSize() == 9, res.getSize());
   cout << "   Day:            '" << res(1) << "'" << endl;
   cout << "   Month:          '" << res(2) << "'" << endl;
   cout << "   Day of Month:   '" << res(3) << "'" << endl;
   cout << "   Hours:          '" << res(4) << "'" << endl;
   cout << "   Min:            '" << res(5) << "'" << endl;
   cout << "   Sec:            '" << res(6) << "'" << endl;
   cout << "   Year:           '" << res(7) << "'" << endl;

   // --- text markers
   FuncLib markers;
   markers("ALAT") = "10.68";
   markers("VOL") = "1000";
   markers("LIST") = "abc,def,ghi,jkl";
   cout << markers.apply("a = ***ALAT*** Bohr\n");
   cout << markers.apply("a = ###ALAT### Bohr\n");
   cout << markers.apply("b = ***BLAT:unknown*** Bohr\n");
   cout << markers.apply("ALAT ***ALAT?does:does not*** exist.\n");
   cout << markers.apply("BLAT ***BLAT?does:does not*** exist.\n");
   cout << markers.apply("volume is ***VOL==1000?one thousand ([#]):a number([#])*** Bohr^3\n");
   cout << markers.apply("sqrt(***VOL*** Bohr) is ***sqrt(VOL)*** Bohr\n");
   cout << markers.apply("cbrt(***VOL*** Bohr) is ***=$VOL^(1/3)*** Bohr\n");
   cout << markers.apply("case conversion: upper ***upper(aBcDe)***\n");
   cout << markers.apply("                 lower ***lower(aBcDe)***\n");
   cout << markers.apply("volume is ***VOL<1000?less:greater/equal*** than 1000 Bohr^3\n");
   cout << markers.apply("volume is ***VOL>=1000?greater/equal:less*** than 1000 Bohr^3 and b=***BLAT:unknown***\n");
   cout << markers.apply("list: ***LIST[-6]***\n");

   // --- text blocks
   SxString tmpl = SxString::read ("Makefile.am");
   SxMap<SxString,SxString> blocks;
   blocks ("SxStringBlock_Example") = "This text block has been replaced";
   cout << tmpl.substBlocks (blocks, true, "###", "###");
   SxList<SxString> tmpls = tmpl.getBlocks ("### BEGIN SxStringBlock_Example ###",
                                            "### END SxStringBlock_Example ###",
                                            false);
   cout << ">" << tmpls(0) << "<\n";

   // --- number conversion
   float f1, f2;
   try                    { f1 = SxString("123.456").toFloat (); } 
   catch (SxException e)  { e.print (); f1 = 0.;}
   try                    { f2 = SxString("123x456").toFloat (); } 
   catch (SxException e)  { e.print (); f2 = 0.;}
   cout << "f1 = " << f1 << endl;
   cout << "f2 = " << f2 << endl;

   return 0; */
}

