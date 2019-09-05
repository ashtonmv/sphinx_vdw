#include <SxCLI.h>
#include <SxDemo1Parser.h>
#include <SxDemo2Parser.h>
#include <SxDemo3Parser.h>
#include <SxDemo4Parser.h>
#include <SxDemo5Parser.h>
#include <SxDemo6Parser.h>
#include <SxDemo6Schema.h>
#include <SxGQuery.h>


typedef typename SxGQExprBase::SelSet SelSet;
using namespace sx;


void exampleDemo1 () {
   cout << "\nexample demo1 begin\n";

   SxDemo1Parser parser;

   try {
      // -- read/parse an alpha/numeric unqouted string
      parser.readString ("dfasdlkfasdfl");
   } catch (SxException e) {
      e.print ();
   }

   cout << "\nexample demo1 end\n";
}

/*
Demo2Parser is a stateful parser that supports
multiline comments and quoted strings. 
*/
void exampleDemo2 () {
   cout << "\nexample demo2 begin\n";

   SxDemo2Parser parser;

   try {
      // -- read/parse an ascii unqouted string
      parser.readString ("/**/\"abc\":\"def\"");
   } catch (SxException e) {
      e.print ();
   }

   cout << "\nexample demo2 end\n";
}

/*
Demo3Parser is a stateful parser that supports
multiline comments and quoted UTF-8 strings.
The parser expects key value pairs of strings. 
*/
void exampleDemo3 () {
   cout << "\nexample demo3 begin\n";

   SxDemo3Parser parser;

   try {
      SxString s = SxString::fromUtf8 (u8"/*co/*mme*/nt*/\"παράδειγμα\":\"aad\\n\\nds\"");
      parser.readString (s);
   } catch (SxException e) {
      e.print ();
   }

   cout << "\nexample demo3 end\n";
}

/*
Demo4Parser is a stateful parser that supports
multiline comments, quoted UTF-8 strings and
C like include statements that allows to
include input data from other files.
The parser expects key value pairs of strings. 
*/
void exampleDemo4 () {
   cout << "\nexample demo4 begin\n";

   SxDemo4Parser parser;

   try {
      SxString s = SxString::unicodeFromUtf8 (u8"/*comment*/\"παράδειγμα\":\"ตัว\\n\\nอย่าง\"");
      parser.readString (s);
   } catch (SxException e) {
      e.print ();
   }

   cout << "\nexample demo4 end\n";
}

/*
Demo5Parser additionally provides a simple 
JSON parser that supports key/value pairs.
The values can be string/int/float or an
object. The parser uses stack to store the
intermediate results.
*/
void exampleDemo5 () {
   cout << "\nexample demo5 begin\n";

   SxDemo5Parser parser;

   try {
      SxString s = SxString::fromUtf8 (u8"{\"k1\":\"v1\", \
                                           \"k2\":22, \
                                           \"k3\":33.4 \
                                          }");
      parser.readString (s);
   } catch (SxException e) {
      e.print ();
   }

   cout << "\nexample demo5 end\n";
}

/*
Demo5Parser provides a simple JSON
parser that supports key/value pairs.
The values can be string/int/float or an
object. The parser uses stack to store the
intermediate results. Additionally, it
provides simple schema validation. 
*/
void exampleDemo6 () {
   cout << "\nexample demo6 begin\n";

   SxDemo6Parser parser;

   try {
      SxString s = SxString::fromUtf8 (u8"{\"k1\":\"v1\", \
                                           \"k2\":22, \
                                           \"k3\":33.4 \
                                          }");
      parser.readString (s);
      SxDemo6Schema schema;
      // -- validate data based on given schema
      // -- just verifies that each string is below 10 in length
      if (schema.validate (parser.getAst ())) {
         cout << "Schema validation successfull\n";
      } else {
         cout << "Schema validation failed\n";
      }
   } catch (SxException e) {
      e.print ();
   }

   cout << "\nexample demo6 end\n";
}

int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   SxString demoId = cli.option ("--id", "string",
                                 "name of the demo").toString ("");
   cli.finalize ();

   if (demoId == "demo1") {
      exampleDemo1 ();
   } else if (demoId == "demo2") {
      exampleDemo2 ();
   } else if (demoId == "demo3") {
      exampleDemo3 ();
   } else if (demoId == "demo4") {
      exampleDemo4 ();
   } else if (demoId == "demo5") {
      exampleDemo5 ();
   } else if (demoId == "demo6") {
      exampleDemo6 ();
   } else if (demoId == "") {
      exampleDemo1 ();
      exampleDemo2 ();
      exampleDemo3 ();
      exampleDemo4 ();
      exampleDemo5 ();
      exampleDemo6 ();
   } else {
      cerr << "invalid demo name\n";
      return 1;
   }

   return 0;
}
