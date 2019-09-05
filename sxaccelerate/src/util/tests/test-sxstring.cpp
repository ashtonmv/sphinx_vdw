// Before you execute this test program, make sure that your terminal was set
// to the text encoding UTF-8. Otherwise, some weird glyphs appear on the
// screen instead of the nice non-ASCII characters.

// After executing, you might want to check manually whether a file was created
// whose name begins with "sxstring-testfile-". Both the filename and the file
// must contain correct Asian characters.

#include <SxString.h>
#include <SxRegex.h>
#include <SxSed.h>
#include <SxCLI.h>

#define INCOMPLETE (1)

// --- errors in release and debug mode
#ifdef NDEBUG
#   define SX_TEST(expr)                                                     \
           if ( !(expr) )  {                                                 \
              std::cout << std::endl << "ASSERTATION FAILED in "             \
                        << __FILE__ << ", line " << __LINE__ << "!\n"        \
                        << "Function: " << SX_PRETTY_FUNCTION << "\n";       \
              std::cout << #expr "\n\n";                                     \
              SX_EXIT;                                                       \
           }
#else
#   define SX_TEST(p) SX_CHECK (p)
#endif /* NDEBUG */

#define SX_TEST_EX(p,b)            \
   {                               \
      bool catchEx = false;        \
      try {                        \
         SX_TEST (p);              \
      }  catch (...)  {            \
         catchEx = true;           \
      }                            \
      SX_TEST (catchEx == (b));    \
   }

int invalidBist ()
// built-in self tests
{
#if INCOMPLETE || (defined(NDEBUG))
   puts ("NOT running SxStringBist ()!");
#else
   puts ("Running SxStringBist ()");
   if (SxString::SxStringBist ())  return 1;
   if (SxUnicode::SxUnicodeBist ())  return 1;
#endif

   return 0;
}

int invalidResize ()
{
   SX_TRACE ();
   SxString str;
   
   str.resize (100);
   SX_TEST (str.ascii() && strlen (str.ascii ()) == 100);
   
   str.resize (50, false);
   SX_TEST (str.ascii() && strlen (str.ascii ()) == 50);
   
   SxString zero('\0');
   SX_TEST (zero.getSize() == 1 && zero.ascii() && strlen (zero.ascii ()) == 1);
   
   zero.append ('\0');
   SX_TEST (zero.getSize() == 2 && zero.ascii() && strlen (zero.ascii ()) == 2);
   
   str.resize (5, false);
   str(2) = '0';
   SX_TEST (str.getSize() == 5);
   SX_TEST (str.ascii() && strlen (str.ascii ()) == 5);
   SX_TEST (str == "  0  ");
   
   str.resize (6, true);
   SX_TEST (str == "  0   ");
   str.resize (3, true);
   SX_TEST (str == "  0");
   
   return 0;
}

int invalidConstructor ()
{
   SX_TRACE ();
   SxString empty;
   SX_TEST (empty.getSize () == 0);
   SX_TEST (empty.ascii () != NULL);
   SX_TEST (strlen (empty.ascii()) == 0);
   SX_TEST (!empty);
   SX_TEST (empty.isEmpty ());
   
   SxString a ('a');
   SX_TEST (a.getSize () == 1);
   SX_TEST (!a.isEmpty ());
   SX_TEST (a.ascii ());
   SX_TEST (strlen (a.ascii()) == 1);
   SX_TEST (strcmp ("a", a.ascii ()) == 0);

   SxString abc2("abc", 2);
   SX_TEST (!abc2 == false);
   SX_TEST (abc2.getSize() == 2);
   SX_TEST (abc2.isEmpty() == false);
   SX_TEST (strlen (abc2.ascii()) == 2);
   SX_TEST (strcmp ("ab", abc2.ascii()) == 0);
   
   SxString abc ("abc");
   SX_TEST (abc.getSize () == 3);
   SX_TEST (!abc.isEmpty ());
   SX_TEST (abc.ascii ());
   SX_TEST (strlen (abc.ascii()) == 3);
   SX_TEST (strcmp (abc.ascii (), "abc") == 0);

   SxString abcd ("abc", "d");
   SX_TEST (abcd.getSize () == 4 && strcmp (abcd.ascii (), "abcd") == 0);

   SxString ab ("ab", (const char*)NULL);
   SX_TEST (ab.getSize () == 2 && strcmp (ab.ascii (), "ab") == 0);

   SxString cd ((const char*)NULL, "cd");
   SX_TEST (cd.getSize () == 2 && strcmp (cd.ascii (), "cd") == 0);

   SxString nulla ((const char*)NULL, (const char*)NULL);
   SX_TEST (nulla.getSize () == 0 && nulla.ascii ());

   SxString nullb ((const char *)NULL);
   SX_TEST (nullb.getSize () == 0 && nullb.ascii ());
   
   SxString abcCopy (abc);
   SX_TEST (abcCopy.getSize () == 3);
   SX_TEST (!abcCopy.isEmpty ());
   SX_TEST (abcCopy.ascii ());
   SX_TEST (strlen (abcCopy.ascii()) == 3);
   SX_TEST (strcmp (abcCopy.ascii (), "abc") == 0);

   SxString empty2 (empty);
   SX_TEST (empty2.getSize () == 0);
   SX_TEST (empty2.isEmpty ());
   SX_TEST (empty2.ascii () != NULL);
   SX_TEST (strlen (empty2.ascii()) == 0);

   return 0;
}

int invalidSet ()
{
   SX_TRACE ();
   SxString stars;
   
   stars.set ('*');
   SX_TEST (!stars);
   
   stars.resize (3);
   stars.set ('*');
   SX_TEST (stars.getSize () == 3 && stars == "***");
   
   return 0;
}

int invalidFormat ()
{
   SX_TRACE ();
   SxString str;
   int a = 5;
   int b = 3;
   str = SxString::sprintf ("%d plus %d is %d", a, b, a + b);
   SX_TEST (str.getSize () == 13 && str == "5 plus 3 is 8");

   // --- should produce gcc warning message
   //str = SxString::sprintf ("%f plus %d is %d", a, b);
   //SX_TEST (str.ascii() && "invalid format string");

   SxArray<char> big (SxString::maxFormatLength);
   big.set (32);
   big(SxString::maxFormatLength - 1) = 0;
   str = SxString::sprintf ("%s", big.elements);
   SX_TEST (str.ascii() && strcmp (str.ascii(), big.elements) == 0
            && "fully used");
   
   // --- fall-back if maxFormatLength was not sufficient
   big(0) = 'b';
   str = SxString::sprintf ("a%s", big.elements);
   SX_TEST (str.getSize() == SxString::maxFormatLength && str.head(2)=="ab");

   return 0;
}

int invalidNumber ()
{
   SX_TRACE ();
   SX_TEST (SxString::toNumber ("0", -1) == 0);
   SX_TEST (SxString::toNumber ("10", 0) == 10);
   SX_TEST (SxString::toNumber ("   +10", 0) == 10);
   SX_TEST (SxString::toNumber ("   -28", 0) == -28);
   SX_TEST (SxString::toNumber ("  b110", 0) == 6);
   SX_TEST (SxString::toNumber ("  -060", 0) == -48);
   SX_TEST (SxString::toNumber ("0xffff", 0) == 65535);
   
   SX_TEST (SxString("1").toNumber<int>() == 1);
   SX_TEST (SxString("001").toNumber<int>() == 1);
   SX_TEST (SxString("10").toNumber<int>() == 10);
   SX_TEST (SxString("+28").toNumber<int>() == 28);
   SX_TEST (SxString("-28").toNumber<int>() == -28);
   SX_TEST (SxString("110").toNumber<int>(2) == 6);
   SX_TEST (SxString("-60").toNumber<int>(8) == -48);
   SX_TEST (SxString("10").toNumber<int>(10) == 10);
   SX_TEST (SxString("ffff").toNumber<unsigned int>(16) == 65535);
   SX_TEST (SxString("-ffff").toNumber<int>(16) == -65535);
      
   bool failed = false;
   SX_TEST (SxString::toNumber ("    ", 0, &failed) == 0 && failed);
   SX_TEST (SxString::toNumber ("    ", 0) == 0);
   SX_TEST (SxString::toNumber<int> ("    ") == 0);
   SX_TEST (SxString::toNumber (" - 1", 0, &failed) == 0 && failed);
   SX_TEST (SxString::toNumber ("  ff", 0, &failed) == 0 && failed);
   SX_TEST (SxString::toNumber (" 008", 0, &failed) == 0 && failed);
   SX_TEST (SxString::toNumber (" 0xg", 0, &failed) == 0 && failed);
  
   SX_TEST_EX (SxString("1").toNumber<int>() == 1, false); 
   SX_TEST_EX (SxString(" 10km").toNumber<int>() == 0, true);
   SX_TEST_EX (SxString("32-bit").toNumber<int>() == 0, true);
   SX_TEST_EX (SxString("6/2").toNumber<int>() == 0, true);
   SX_TEST_EX (SxString("1+1").toNumber<int>() == 0, true);
   SX_TEST_EX (SxString("e").toNumber<int>() == 0, true);

   const SxString null;
   SX_TEST_EX (null.toNumber<int>() == 1, true);

   const SxString empty (" ");
   SX_TEST_EX (empty.toNumber<int>() == 1, true);
   
   try { SxString("invalid number").toNumber<int>(); }
   catch (SxException e) { } //{ e.print (true); }
   
   return 0;
}

int invalidNumber2 ()
{
   // --- numbers, SxString(int), isInt, toInt
   SxString sInt((int)28);
   SX_TEST (strcmp (sInt.ascii(), "28") == 0);
   SX_TEST (sInt.isInt() == true);
   SX_TEST_EX (sInt.toInt() == 28, false);
   
   SxString sIntErr("a");
   SX_TEST (sIntErr.isInt() == false);
   SX_TEST_EX (sIntErr.toInt() == 10, true);

   SxString sIntForm((int)7, 3);
   SX_TEST (strcmp (sIntForm.ascii(), "007") == 0);
   
   SxString sUInt((unsigned int)27);
   SX_TEST (strcmp (sUInt.ascii(), "27") == 0);
   
   SxString sLong((long)16);
   SX_TEST (strcmp (sLong.ascii(), "16") == 0);
   SX_TEST (sLong.isInt64 () == true);
   SX_TEST_EX (sLong.toInt64 () == 16, false);

   SxString sFloat((float)3.1234f);
   SX_TEST (strcmp (sFloat.ascii(), "3.1234") == 0);
   SX_TEST (sFloat.isFloat() == true);
   SX_TEST_EX (sFloat.toFloat() == 3.1234f, false);
   
   SX_TEST_EX (SxString("123.456").toFloat() == 123.456f, false);
   SX_TEST_EX (SxString("123x456").toFloat() == 0.f, true);
   
   SxString sDouble((double)3.1234);
   SX_TEST (strcmp (sDouble.ascii(), "3.1234") == 0);
   SX_TEST (sDouble.isDouble() == true);
   SX_TEST_EX (sDouble.toDouble() == 3.1234, false);

   SxString sDoubleForm(3.123456789, "%.9f");
   SX_TEST (strcmp (sDoubleForm.ascii(), "3.123456789") == 0);

   return 0;
}

int test (const SxString &str, int base)
{
   return SxString::toNumber (str.ascii(), -1, NULL, base);
}

int invalidNumber3 ()
{
   SX_TRACE ();

   // --- auto base
   SX_TEST(test ("110",   0) == 110);
   SX_TEST(test ("0110",  0) == 72);
   SX_TEST(test ("b110",  0) == 6);
   SX_TEST(test ("B110",  0) == 6);
   SX_TEST(test ("0x110", 0) == 272);
   SX_TEST(test ("0X110", 0) == 272);

   SX_TEST(test ("0 110",  0) == -1);
   SX_TEST(test ("b 110",  0) == -1);
   SX_TEST(test ("0x 110", 0) == -1);
   SX_TEST(test ("b",  0) == -1);
   SX_TEST(test ("0x", 0) == -1);

   // --- binary
   SX_TEST(test ("110",   2) == 6);
   SX_TEST(test ("0110",  2) == 6);
   SX_TEST(test ("b110",  2) == -1);
   SX_TEST(test ("0x110", 2) == -1);
   SX_TEST(test ("-110",  2) == -6);
   // --- octal
   SX_TEST(test ("110",   8) == 72);
   SX_TEST(test ("0110",  8) == 72);
   SX_TEST(test ("b110",  8) == -1);
   SX_TEST(test ("0x110", 8) == -1);
   SX_TEST(test ("-110",  8) == -72);
   // --- decimal
   SX_TEST(test ("110",   10) == 110);
   SX_TEST(test ("0110",  10) == 110);
   SX_TEST(test ("0010",  10) == 10);
   SX_TEST(test ("b110",  10) == -1);
   SX_TEST(test ("0x110", 10) == -1);
   SX_TEST(test ("-0110",  10) == -110);
   // --- hex
   SX_TEST(test ("110",   16) == 272);
   SX_TEST(test ("0110",  16) == 272);
   SX_TEST(test ("b110",  16) == 45328);
   SX_TEST(test ("0x110", 16) == -1);
   SX_TEST(test ("-110",   16) == -272);

   // --- space
   SX_TEST(test (" 110",   0) == 110);
   SX_TEST(test ("110 ",   0) == 110);
   SX_TEST(test (" 110 ",  0) == 110);
   SX_TEST(test (" 110,",  0) == -1);
   SX_TEST(test (",110 ",  0) == -1);
   SX_TEST(test ("",       0) == -1);
   SX_TEST(test (" ",      0) == -1);
   SX_TEST(test (",",      0) == -1);

   // --- sign
   SX_TEST(test ("+110",   0) == 110);
   SX_TEST(test ("+0110",  0) == 72);
   SX_TEST(test ("+b110",  0) == 6);
   SX_TEST(test ("+0x110", 0) == 272);
   SX_TEST(test ("-110",   0) == -110);
   SX_TEST(test ("-0110",  0) == -72);
   SX_TEST(test ("-b110",  0) == -6);
   SX_TEST(test ("-0x110", 0) == -272);

   SX_TEST(test (" +110",   0) == 110);
   SX_TEST(test (" -110",   0) == -110);
   SX_TEST(test (" + 110",   0) == -1);
   SX_TEST(test (" - 110",   0) == -1);

   SX_TEST(SxString("0").isInt());
   SX_TEST(SxString("1").isInt());
   SX_TEST(SxString(" 1").isInt());
   SX_TEST(SxString("\t1").isInt());
   SX_TEST(SxString("1 ").isInt());
   SX_TEST(SxString(" 1 ").isInt());
   SX_TEST(SxString("").isInt() == false);
   SX_TEST(SxString("a").isInt() == false);
   SX_TEST(SxString("1,").isInt() == false);
   SX_TEST(SxString(",1").isInt() == false);
   SX_TEST(SxString("1 ,").isInt() == false);
   SX_TEST(SxString(", 1").isInt() == false);

   SX_TEST(SxString("0").isInt64 ());

   return 0;
}

// test with the previous version
/*
int toInt (const SxString &str, bool *error)
{
   SX_CHECK (error);
   *error = false;
   int value = 0;
   SxArray<char> buffer(str.getSize ());
   if (sscanf (str.ascii (), "%d %s", &value, buffer.elements) != 1)  {
      *error = true;
      return 0;
   }
   return value;
}
bool isInt (const SxString &str)
{
   int value = 0;
   SxArray<char> buffer(str.getSize ());
   if (sscanf (str.ascii (), "%d %s", &value, buffer.elements) != 1)
      return false;
   else
      return true;
}
void testInt (const SxString &str)
{
   bool error = false;
   cout << "string " << str << endl;
   cout << "   old isInt " << isInt (str) << endl;
   cout << "   old toInt " << toInt (str, &error) << endl;
   cout << "   old err   " << error << endl;
   cout << "   new isInt " << str.isInt () << endl;
   cout << "   new toInt " << str.toInt (&error) << endl;
   cout << "   new err   " << error << endl;
   cout << "   new isInt64 " << str.isInt64 () << endl;
   cout << "   new toInt64 " << str.toInt64 (&error) << endl;
   cout << "   new err    " << error << endl;
}
*/

template<class T>
inline bool testNum (const SxString &str_, T default_)
{
   bool error = true;
   T v = SxString::toNumber (str_.ascii(), default_, &error, 10);
   SxString str = SxString(v);
   return (str == str_);
}

#define SX_TEST_U(str,res)                                \
   SX_TEST(testNum<unsigned int>(str, 0) == res);         \
   SX_TEST(testNum<unsigned long>(str, 0) == res);        \
   SX_TEST(testNum<unsigned long long>(str, 0) == res);

#define SX_TEST_S(str,res)                                \
   SX_TEST(testNum<signed char>(str, 0) == res);          \
   SX_TEST(testNum<int>(str, 0) == res);                  \
   SX_TEST(testNum<long>(str, 0) == res);                 \
   SX_TEST(testNum<long long>(str, 0) == res);            

#define SX_TEST_NUM(str,res)                              \
   SX_TEST_U(str,res);                                    \
   SX_TEST_S(str,res);

int invalidNumber4 ()
{
   // 8                    -128 .. 127
   //                         0 .. 255
   // 32            -2147483648 .. 2147483647
   //                         0 .. 4294967295
   // 64   -9223372036854775808 .. 9223372036854775807
   //                         0 .. 18446744073709551615
   // 128  -170141183460469231731687303715884105728
   //                         0 .. 340282366920938463463374607431768211455

   // testInt ("2147483647"); // ok
   // testInt ("2147483648");
   // testInt ("18446744073709551616");
   // testInt ("340282366920938463463374607431768211455");
   // testInt ("340282366920938463463374607431768211456");

   // --- common numbers
   SX_TEST_NUM("0",true);
   SX_TEST_S("127",true);
   SX_TEST_S("-128",true);
   SX_TEST_U("255",true);
   SX_TEST_U("-1",false);
   SX_TEST_NUM("",false);
   SX_TEST_NUM("str",false);

   // --- used to choose template T
   signed char c = 1;
   int i = 1;
   long l = 1;
   unsigned char uc = 1;
   unsigned int ui = 1;
   unsigned long ul = 1;

   // --- signed char
   SX_TEST( testNum ("0", c));
   SX_TEST( testNum ("-128", c));
   SX_TEST( testNum ("127", c));
   SX_TEST(!testNum ("-129", c));
   SX_TEST(!testNum ("128", c));

   // --- unsigned char
   SX_TEST( testNum ("0", uc));
   SX_TEST( testNum ("255", uc));
   SX_TEST(!testNum ("-1", uc));
   SX_TEST(!testNum ("256", uc));

   // --- int
   SX_TEST(testNum ("0", i));
   SX_TEST(testNum ("-2147483648", i));
   SX_TEST(testNum ("2147483647", i));
   SX_TEST(SxString("2147483647").isInt ());
   SX_TEST_EX(SxString("2147483647").toInt() == 2147483647, false);
   if (sizeof(int) == 4) {
      SX_TEST(!testNum ("-2147483649", i));
      SX_TEST(!testNum ("2147483648", i));
      SX_TEST(!SxString("2147483648").isInt ());
      SX_TEST_EX (SxString("2147483648").toInt() == 0, true);
   }
   if (sizeof(int) == 8) {
      SX_TEST( testNum ("-9223372036854775808", i));
      SX_TEST( testNum ( "9223372036854775807", i));
      SX_TEST( SxString( "9223372036854775807").isInt ());
      SX_TEST(!testNum ("-9223372036854775809", i));
      SX_TEST(!testNum ( "9223372036854775808", i));
      SX_TEST(!SxString( "9223372036854775808").isInt ());

      int maxInt64 = sxmax<int>();
      SX_TEST_EX(SxString("9223372036854775807").toInt() == maxInt64, false);
      SX_TEST_EX(SxString("9223372036854775808").toInt() == 0, true);
   }
   if (sizeof(int) == 16) {
      SX_TEST( SxString("-170141183460469231731687303715884105728").isInt());
      SX_TEST( SxString( "170141183460469231731687303715884105727").isInt());
      SX_TEST( testNum ("-170141183460469231731687303715884105728", i));
      SX_TEST( testNum ( "170141183460469231731687303715884105727", i));
      SX_TEST(!SxString("-170141183460469231731687303715884105729").isInt());
      SX_TEST(!SxString( "170141183460469231731687303715884105728").isInt());
      SX_TEST(!testNum ("-170141183460469231731687303715884105729", i));
      SX_TEST(!testNum ( "170141183460469231731687303715884105728", i));
   }

   // --- long
   SX_TEST(testNum ("0", l));
   SX_TEST(testNum ("-2147483648", l));
   SX_TEST(testNum ( "2147483647", l));
   SX_TEST(SxString( "2147483647").isInt64 ());
   SX_TEST_EX(SxString("2147483647").toInt64 () == 2147483647, false);
   if (sizeof(long) == 8) {
      SX_TEST( testNum ("-9223372036854775808", l));
      SX_TEST( testNum ( "9223372036854775807", l));
      SX_TEST( SxString( "9223372036854775807").isInt64 ());
      SX_TEST(!testNum ("-9223372036854775809", l));
      SX_TEST(!testNum ( "9223372036854775808", l));
      SX_TEST(!SxString( "9223372036854775808").isInt64 ());

      long maxLong64 = sxmax<long>();
      SX_TEST_EX(SxString("9223372036854775807").toInt64 () == maxLong64,false);
      SX_TEST_EX(SxString("9223372036854775808").toInt64 () == 0, true);
   }
   if (sizeof(long) == 16) {
      SX_TEST( SxString("-170141183460469231731687303715884105728").isInt64 ());
      SX_TEST( SxString( "170141183460469231731687303715884105727").isInt64 ());
      SX_TEST( testNum ("-170141183460469231731687303715884105728", l));
      SX_TEST( testNum ( "170141183460469231731687303715884105727", l));
      SX_TEST(!SxString("-170141183460469231731687303715884105729").isInt64 ());
      SX_TEST(!SxString( "170141183460469231731687303715884105728").isInt64 ());
      SX_TEST(!testNum ("-170141183460469231731687303715884105729", l));
      SX_TEST(!testNum ( "170141183460469231731687303715884105728", l));
   }

   // --- unsigned int
   SX_TEST( testNum ("0", ui));
   SX_TEST( testNum ("4294967295", ui));
   SX_TEST(!testNum ("-1", ui));
   if (sizeof(unsigned int) == 4) {
      SX_TEST(!testNum ("4294967296", ui));
   }
   if (sizeof(unsigned int) == 8) {
      SX_TEST( testNum ("18446744073709551615", ui));
      SX_TEST(!testNum ("18446744073709551616", ui));
   }
   if (sizeof(unsigned int) == 16) {
      SX_TEST( testNum ("340282366920938463463374607431768211455", ui));
      SX_TEST(!testNum ("340282366920938463463374607431768211456", ui));
   }

   // --- unsigned long
   SX_TEST( testNum ("0", ul));
   SX_TEST( testNum ("4294967295", ul));
   SX_TEST(!testNum ("-1", ul));
   if (sizeof(unsigned long) == 8) {
      SX_TEST( testNum ("18446744073709551615", ul));
      SX_TEST(!testNum ("18446744073709551616", ul));
   }
   if (sizeof(unsigned long) == 16) {
      SX_TEST( testNum ("340282366920938463463374607431768211455", ul));
      SX_TEST(!testNum ("340282366920938463463374607431768211456", ul));
   }

   return 0;
}

int invalidIndex ()
{
   SX_TRACE ();
   SxString abc ("abc");

//   SX_TEST_EX (abc.at(-1) == 'a' && "out of range Exception", true);
//   SX_TEST_EX (abc.at(3) == 'a' && "out of range Exception", true);
//   SX_TEST_EX (SxString().at(0) == 'a', true);
   
   SX_TEST (abc(0) == 'a' && abc(1) == 'b' && abc(2) == 'c');
   SX_TEST (abc.at(0) == 'a' && abc.at(1) == 'b' && abc.at(2) == 'c');
   
//   abc.set (0, 'A');
//   abc.set (1, 'B');
//   abc.set (2, 'C');
//   SX_TEST (abc(0) == 'A' && abc(1) == 'B' && abc(2) == 'C');

   return 0;
}

int invalidReplace ()
{
   SxString xyz = SxString("xyz");
   SxString s = SxString("abc");
   s.replace ("xyz");
   SX_TEST (s == xyz);
   xyz.replace ("");
   SX_TEST (xyz == "");
   SX_TEST (xyz.getNBytes () == 0);

   return 0;
}

int invalidSubString ()
{
   SX_TRACE ();
   SX_TEST (!SxString().subString (0, 0));
   SX_TEST (!SxString().subString (-1, 1));
   SX_TEST (!SxString().head (1));
   SX_TEST (!SxString().tail (1));

   SX_TEST (!SxString("").subString (0, 0));
   SX_TEST (!SxString("").head (1));
   SX_TEST (!SxString("").tail (1));

   SX_TEST (SxString("a").subString (0, 0) == "a");
   SX_TEST (SxString("a").head (1) == "a");
   SX_TEST (SxString("a").tail (1) == "a");

   SX_TEST (SxString("abc").subString (1, 1) == "b");
   SX_TEST (SxString("abc").head (1) == "a");
   SX_TEST (SxString("abc").tail (1) == "c");

   SX_TEST (SxString("abcde").subString (2, 3) == "cd");
   SX_TEST (SxString("abcde").head (2) == "ab");
   SX_TEST (SxString("abcde").tail (2) == "de");

   SX_TEST (SxString("abc").subString (1, 3) == "bc");
   SX_TEST (SxString("abc").subString (-1, 3) == "abc");
   SX_TEST (SxString("abc").head (5) == "abc");
   SX_TEST (SxString("abc").tail (5) == "abc");
   SX_TEST (!SxString("abc").subString (2, 1));
   SX_TEST (!SxString("abc").head (0));
   SX_TEST (!SxString("abc").tail (0));
   SX_TEST (!SxString("abc").head (-4));
   SX_TEST (!SxString("abc").tail (-4));
   
   SxString text = "This text contains some whitespaces.";
   SX_TEST (text.left ("contains") == "This text ");
   SX_TEST (text.right ("contains") == " some whitespaces.");
   SX_TEST (text.substitute ("text", "string") ==
            "This string contains some whitespaces.");

   return 0;
}

int invalidSubstitute ()
{
   SX_TRACE ();
   SxString null;
   
   SX_TEST (!SxString().substitute ("a", "b"));
   
   SX_TEST (SxString("").substitute ("a", "b") == "");
   SX_TEST (SxString("").substitute ("", "b") == "");
   SX_TEST (SxString("").substitute ("a", "") == "");
   SX_TEST (SxString("").substitute ("", "") == "");
   
   SX_TEST (SxString("abc").substitute (null, null) == "abc");
   SX_TEST (SxString("abc").substitute (null, "b") == "abc");
   SX_TEST (SxString("abc").substitute ("a", null) == "bc");
   
   SX_TEST (SxString("abc").substitute ("", "") == "abc");
   SX_TEST (SxString("abc").substitute ("", "b") == "abc");
   SX_TEST (SxString("abc").substitute ("a", "") == "bc");
   SX_TEST (SxString("abc").substitute ("a", "b") == "bbc");
   
   SX_TEST (SxString("a").substitute ("a", "abc") == "abc");
   SX_TEST (SxString("abc").substitute ("abc", "a") == "a");
   SX_TEST (SxString("aaa").substitute ("aa", "b") == "ba");
   SX_TEST (SxString("ab").substitute ("abc", ".") == "ab");
   SX_TEST (SxString("abd").substitute ("abc", ".") == "abd");

   SX_TEST (SxString("abcabd").substitute ("ab", "xyz") == "xyzcxyzd");
   SX_TEST (SxString("abcabd").substitute ("ab", "xyz", 1) == "xyzcabd");
   SX_TEST (SxString("abcabd").substitute ("ab", "") == "cd");
   SX_TEST (SxString(" s tri  ng").substitute (" ", "") == "string");

   return 0;
}

int invalidFind ()
{
   SX_TRACE ();
   SX_TEST (SxString().find (SxString()) < 0);
   SX_TEST (SxString().find ("") < 0);
   SX_TEST (SxString().find ("x") < 0);

   SX_TEST (SxString("").find (SxString()) < 0);
   SX_TEST (SxString("").find ("") < 0);
   SX_TEST (SxString("").find ("x") < 0);

   SX_TEST (SxString("x").find (SxString()) < 0);
   SX_TEST (SxString("x").find ("") < 0);
   SX_TEST (SxString("x").find ("x") == 0);

   SX_TEST (SxString("...").find (".") == 0);
   SX_TEST (SxString("...").find ("x") < 0);
   SX_TEST (SxString(".x.").find ("x") == 1);
   SX_TEST (SxString(".xx").find ("xx") == 1);
   SX_TEST (SxString(".xx.").find ("xx") == 1);
   SX_TEST (SxString(".xxx").find ("xxx") == 1);

   SX_TEST (SxString().findLast (SxString()) < 0);
   SX_TEST (SxString().findLast ("") < 0);
   SX_TEST (SxString().findLast ("x") < 0);

   SX_TEST (SxString("").findLast (SxString()) < 0);
   SX_TEST (SxString("").findLast ("") < 0);
   SX_TEST (SxString("").findLast ("x") < 0);

   SX_TEST (SxString("x").findLast (SxString()) < 0);
   SX_TEST (SxString("x").findLast ("") < 0);
   SX_TEST (SxString("x").findLast ("x") == 0);
   
   SX_TEST (SxString(".").findLast (".") == 0);
   SX_TEST (SxString("...").findLast (".") == 2);
   SX_TEST (SxString("...").findLast ("x") < 0);
   SX_TEST (SxString(".x.").findLast ("x") == 1);
   
   SX_TEST (SxString("img.pcx.tga").findLast (".") == 7);
   SX_TEST (SxString("img.pcx.tga").findLast (".tga") == 7);
   SX_TEST (SxString("img.pcx").findLast (".", 4) < 0);

   return 0;
}

int invalidContains ()
{
   SX_TRACE ();
   SX_TEST (SxString().contains (SxString()) == 0);
   SX_TEST (SxString().contains ("") == 0);
   SX_TEST (SxString().contains ("a") == 0);

   SX_TEST (SxString("").contains (SxString()) == 0);
   SX_TEST (SxString("").contains ("") == 0);
   SX_TEST (SxString("").contains ("a") == 0);

   SX_TEST (SxString("a").contains (SxString()) == 0);
   SX_TEST (SxString("a").contains ("") == 0);
   SX_TEST (SxString("a").contains ("a") == 1);
   SX_TEST (SxString("a").contains ("ab") == 0);
   
   SX_TEST (SxString("aaa").contains ("a") == 3);
   SX_TEST (SxString("aaa").contains ("aa") == 2 && "overlapped");
   SX_TEST (SxString("aaa").containsWhole ("aa") == 1 && "non-overlapped");
   SX_TEST (SxString("strstr, strcasestr").contains ("str") == 4);

   return 0;
}

int invalidInsert ()
{
   SX_TRACE ();
   SxString c;
   c.insert (1, 'y');
   SX_TEST (c.getSize() == 0);
   c.insert (0, 'y');
   SX_TEST (c.getSize() == 1 && c(0) == 'y');
   c.insert (-1, 'x');
   SX_TEST (c.getSize() == 1 && c(0) == 'y');
   c.insert (2, 'x');
   SX_TEST (c.getSize() == 1 && c(0) == 'y');
   c.insert (0, 'x', 0);
   SX_TEST (c.getSize() == 1 && c(0) == 'y');
   c.insert (0, 'x', -1);
   SX_TEST (c.getSize() == 1 && c(0) == 'y');
   c.insert (1, 'x', 2);
   SX_TEST (c.getSize() == 3 && c.ascii() && strcmp (c.ascii (), "yxx") == 0);

   SxString s;
   s.insert (0, "xy");
   SX_TEST (s.getSize() == 2 && s(0) == 'x' && s(1) == 'y');
   s.insert (1, "");
   s.insert (1, (const char*)NULL);
   SX_TEST (s.getSize() == 2 && s(0) == 'x' && s(1) == 'y');
   s.insert (1, "z");
   SX_TEST (s.getSize() == 3 && s.ascii() && strcmp (s.ascii (), "xzy") == 0);
   
   return 0;
}

int invalidPrepend ()
{
   SX_TRACE ();
   SxString c;
   c.prepend ('y');
   SX_TEST (c.getSize() == 1 && c(0) == 'y');
   c.prepend ('x', 0);
   SX_TEST (c.getSize() == 1 && c(0) == 'y');
   c.prepend ('x', -1);
   SX_TEST (c.getSize() == 1 && c(0) == 'y');
   c.prepend ('x', 2);
   SX_TEST (c.getSize() == 3 && c.ascii() && strcmp (c.ascii (), "xxy") == 0);

   SxString s;
   s.prepend ("xy");
   SX_TEST (s.getSize() == 2 && s(0) == 'x' && s(1) == 'y');
   s.prepend ("");
   s.prepend ((const char*)NULL);
   SX_TEST (s.getSize() == 2 && s(0) == 'x' && s(1) == 'y');
   s.prepend ("z");
   SX_TEST (s.getSize() == 3 && s.ascii() && strcmp (s.ascii (), "zxy") == 0);

   return 0;
}

int invalidAppend ()
{
   SX_TRACE ();
   SxString c;
   c.append ('y');
   SX_TEST (c.getSize() == 1 && c(0) == 'y');
   c.append ('x', 0);
   SX_TEST (c.getSize() == 1 && c(0) == 'y');
   c.append ('x', -1);
   SX_TEST (c.getSize() == 1 && c(0) == 'y');
   c.append ('x', 2);
   SX_TEST (c.getSize() == 3 && c.ascii() && strcmp (c.ascii (), "yxx") == 0);
   
   SxString s;
   s.append ("xy");
   SX_TEST (s.getSize() == 2 && s(0) == 'x' && s(1) == 'y');
   s.append ("");
   s.append ((const char*)NULL);
   SX_TEST (s.getSize() == 2 && s(0) == 'x' && s(1) == 'y');
   s.append ("z");
   SX_TEST (s.getSize() == 3 && s.ascii() && strcmp (s.ascii (), "xyz") == 0);
   
   SxString abc;
   SX_TEST (abc.getSize() == 0 && abc.ascii () && strlen (abc.ascii()) == 0);
   abc.append ("");
   SX_TEST (abc.getSize() == 0 && abc.ascii () && strlen (abc.ascii()) == 0);
   abc.append (SxString());
   SX_TEST (abc.getSize() == 0 && abc.ascii () && strlen (abc.ascii()) == 0);
   abc.append ("a");
   SX_TEST (abc.getSize() == 1 && abc(0) == 'a');
   abc.append ("");
   SX_TEST (abc.getSize() == 1 && abc(0) == 'a');
   abc.append (SxString());
   SX_TEST (abc.getSize() == 1 && abc(0) == 'a');
   abc.append (SxString("bc"));
   SX_TEST (abc.getSize() == 3 && strcmp (abc.ascii (), "abc") == 0);
   
   return 0;
}

int invalidRemove ()
{
   SX_TRACE ();
   SxString abc("abcde");
   abc.remove (-1);
   SX_TEST (abc.getSize() == 5 && strcmp (abc.ascii(), "abcde") == 0);
   abc.remove (1);
   SX_TEST (abc.getSize() == 4 && strcmp (abc.ascii(), "acde") == 0);
   abc.remove (0, 2);
   SX_TEST (abc.getSize() == 2 && strcmp (abc.ascii(), "de") == 0);
   abc.remove (1, 3);
   SX_TEST (abc.getSize() == 1 && strcmp (abc.ascii(), "d") == 0);
   abc.remove (0);
   SX_TEST (abc.getSize() == 0 && abc.ascii() && strlen (abc.ascii()) == 0);
   abc.remove (1);
   SX_TEST (abc.getSize() == 0 && abc.ascii() && strlen (abc.ascii()) == 0);

   SxString xyz("xyzw");
   xyz.removeFirst ();
   SX_TEST (xyz.getSize() == 3 && strcmp (xyz.ascii(), "yzw") == 0);
   xyz.removeLast ();
   SX_TEST (xyz.getSize() == 2 && strcmp (xyz.ascii(), "yz") == 0);
   xyz.removeAll ();
   SX_TEST (xyz.getSize() == 0 && xyz.ascii() && strlen (xyz.ascii()) == 0);
   
   return 0;
}

int invalidCase ()
{
   SX_TRACE ();
   SX_TEST (SxString().toUpper () == SxString());
   SX_TEST (SxString().toLower () == SxString());
   SX_TEST (SxString("").toUpper () == "");
   SX_TEST (SxString("").toLower () == "");

   SX_TEST (SxString("a").toUpper () == "A");
   SX_TEST (SxString("A").toLower () == "a");
   
   SX_TEST (SxString("a").toLower () == "a");
   SX_TEST (SxString("A").toUpper () == "A");

   SX_TEST (SxString("aBc1").toUpper () == "ABC1");
   SX_TEST (SxString("AbC1").toLower () == "abc1");
   
   return 0;
}

int invalidPlus ()
{
   SX_TRACE ();
   SxString null;
   SxString empty ("");
   SxString a ("a");
   SxString b ("b");

   SX_TEST (null + null  == null);
   SX_TEST (null + empty == empty);
   SX_TEST (null + "a"   == "a");
   SX_TEST (null + b     == "b");
   
   SX_TEST (empty + null  == empty);
   SX_TEST (empty + empty == empty);
   SX_TEST (empty + "a"   == "a");
   SX_TEST (empty + b     == "b");
   
   SX_TEST ("a" + null  == "a");
   SX_TEST ("a" + empty == "a");
   SX_TEST ("a" + b     == "ab");
   
   SX_TEST (a + null  == "a");
   SX_TEST (a + empty == "a");
   SX_TEST (a + "b"   == "ab");
   SX_TEST (a + b     == "ab");
   
   SxString str;
   str = null;
   SX_TEST (str == null);
   str = empty;
   SX_TEST (str == empty);
   str = "a";
   SX_TEST (str == "a");
   str = b;
   SX_TEST (str == "b");
   str = "abc";
   SX_TEST (str == "abc");
   
   str = null;
   str += null;
   SX_TEST (str == null);

   str = null;
   str += empty;
   SX_TEST (str == empty);
   str += empty;
   SX_TEST (str == empty);

   str = null;
   str += a;
   SX_TEST (str == a);
   str += b;
   SX_TEST (str == "ab");
   str += null;
   SX_TEST (str == "ab");
   str += empty;
   SX_TEST (str == "ab");
   
   // --- disable operator- which otherwise returns invalid memory address
   //     the current operator- is implemented with SX_EXIT
   //
   // str = null;
   // const char *ptr = str - 1;
   //
   // ASSERTATION FAILED in SxString.h, line 717!
   // !"no match for 'operator-' in SxString."
   
   // --- overload operator+ to fix:
   //     - const char *ptr = SxString("abc") + 5;
   //     - const char *ptr = 5 + SxString("abc");
   //     which otherwise return invalid memory address
   SX_TEST (SxString() + 1 == "1");
   SX_TEST (SxString() + '\0' == "0");
   SX_TEST (SxString() + '\t' == "\t");
   SX_TEST (SxString() + 'a' == "a");
   SX_TEST (SxString() + 255 == "255");
   SX_TEST (SxString("pi ") + 3.14 == "pi 3.14");
   
   SX_TEST (SxString("1 + 1 is ") + (1 + 1) == "1 + 1 is 2");
   
   SX_TEST ('\0' + SxString(" km") == "0 km");
   SX_TEST ('5' + SxString(" km") == "5 km");
   SX_TEST (5 + SxString(" km") == "5 km");
   SX_TEST (1e-6 + SxString(" eps") == "1e-06 eps");

   str = null;
   str += '\0';
   SX_TEST (str == "0");
   str += '1';
   SX_TEST (str == "01");
   str += 5;
   SX_TEST (str == "015");
   str += 0.2;
   SX_TEST (str == "0150.2");
   
   return 0;
}

int invalidCompare ()
{
   SX_TRACE ();
   SxString null;
   SxString empty ("");
   SxString a ("a");
   SxString b ("b");

   // ==
   SX_TEST (null == null);
   SX_TEST (null == empty);
   SX_TEST (null == "");
   SX_TEST (!(null == "a"));

   SX_TEST (empty == null);
   SX_TEST (empty == empty);
   SX_TEST (empty == "");
   SX_TEST (!(empty == "a"));

   SX_TEST (!(a == null));
   SX_TEST (!(a == empty));
   SX_TEST (a == "a");
   SX_TEST (!(a == "b"));
   SX_TEST (!(a == "ab"));
   SX_TEST (!(a == "A"));
   
   SX_TEST ("" == null);
   SX_TEST ("" == empty);
   SX_TEST (!("" == a));

   SX_TEST (!("a" == null));
   SX_TEST (!("a" == empty));
   SX_TEST ("a" == a);
   SX_TEST (!("a" == b));

   // <
   SX_TEST (!(null < null));
   SX_TEST (!(null < empty));
   SX_TEST (!(null < ""));
   SX_TEST (null < "a");
   SX_TEST (null < a);

   SX_TEST (!(empty < null));
   SX_TEST (!(empty < empty));
   SX_TEST (!(empty < ""));
   SX_TEST (empty < "a");
   SX_TEST (empty < a);
   
   SX_TEST (!("" < null));
   SX_TEST (!("" < empty));
   SX_TEST ("" < a);

   SX_TEST (!(a < null));
   SX_TEST (!(a < empty));
   SX_TEST (!(a < "a"));
   SX_TEST (a < "b");
   SX_TEST (a < "ab");
   SX_TEST (!(a < "A"));

   SX_TEST (!("a" < null));
   SX_TEST (!("a" < empty));
   SX_TEST (!("a" < a));
   SX_TEST ("a" < b);

   // >
   SX_TEST (!(null > null));
   SX_TEST (!(null > empty));
   SX_TEST (!(null > ""));
   SX_TEST (!(null > "a"));

   SX_TEST (!(empty > null));
   SX_TEST (!(empty > empty));
   SX_TEST (!(empty > ""));
   SX_TEST (!(empty > "a"));

   SX_TEST (!("" > null));
   SX_TEST (!("" > empty));
   SX_TEST (!("" > a));
   
   SX_TEST (a > null);
   SX_TEST (a > empty);
   SX_TEST (!(a > "a"));
   SX_TEST (!(a > "b"));
   SX_TEST (!(a > "ab"));
   SX_TEST (a > "A");

   SX_TEST ("a" > null);
   SX_TEST ("a" > empty);
   SX_TEST (!("a" > a));
   SX_TEST (!("a" > b));

   return 0;
}

int invalidRegexp ()
{
   SxString date = "Mon Mar 28 11:08:31 CEST 2006";
   SxList<SxString> res;
   
   try  {
      SxRegex::test ();
      res = SxRegex("^(\\w+)\\s+(\\w+)\\s+(\\d+)\\s+([0-9]+):"
                    "([0-9]+):([0-9]+)\\s+\\w+\\s+([0-9]+)",
                    "P").match (date);
   }  catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   
   SX_TEST (res.getSize() == 9);
   SX_TEST (res(0) == "");
   SX_TEST (res(1) == "Mon");
   SX_TEST (res(2) == "Mar");
   SX_TEST (res(3) == "28");
   SX_TEST (res(4) == "11");
   SX_TEST (res(5) == "08");
   SX_TEST (res(6) == "31");
   SX_TEST (res(7) == "2006");
   SX_TEST (res(8) == "");
   
   return 0;
}

int invalidSed ()
{
   SxSed::test ();

   SxSed cmd("s/\\\\/\\\\\\\\/g");
   SxString path = cmd.subst ("C:\\foo\\bar\\");
   SX_TEST (path == "C:\\\\foo\\\\bar\\\\");
   
   path = SxSed("\\Q{{DXDRIVE}}\\E","E:\\","g").subst("{{DXDRIVE}}/foo/bar");
   SX_TEST (path == "E:\\/foo/bar");

   path = SxSed("{{DXDRIVE}}","E:\\","gl").subst("{{DXDRIVE}}/foo/bar");
   SX_TEST (path == "E:\\/foo/bar");

   path = SxSed("\\{\\{DXDRIVE\\}\\}","E:\\","g").subst("{{DXDRIVE}}/foo/bar");
   SX_TEST (path == "E:\\/foo/bar");

   path = path.substitute("/","\\");
   SX_TEST (path == "E:\\\\foo\\bar");



   return 0;
}

int invalidWhiteSpace ()
{
   SxString text = "  This    text 		contains   some whitespaces.    ";

   SX_TEST (text.stripWhiteSpace()
            == "This    text 		contains   some whitespaces.");   
   SX_TEST (text.trim() == "This text contains some whitespaces.");
   SX_TEST (text.simplifyWhiteSpace() =="This text contains some whitespaces.");
   SX_TEST (text.removeWhiteSpace() == "Thistextcontainssomewhitespaces.");
   
   return 0;
}

int invalidEmptyString ()
{
   SxArray<char> emptyArray(1); emptyArray(0) = '\0';
   
   // --- create
   SX_TEST (SxString("").size == 0);
   SX_TEST (SxString("", "").size == 0);
   
   SX_TEST (SxString(emptyArray).size == 0);

   // --- resize
   SxString alpha = "alpha";
   alpha.resize (0);
   SX_TEST (alpha.size == 0);
   
   // --- replace
   alpha = "alpha";
   alpha.replace ("");
   SX_TEST (alpha.size == 0);
   
   SX_TEST (SxString("a").substitute("a","").size == 0);
   SX_TEST (SxString(" ").removeWhiteSpace().size == 0);
   SX_TEST (SxString("#").stripComments().size == 0);
   
   // --- insert
   SxString beta;
   beta.insert (0, "");
   SX_TEST (beta.size == 0);
   
   beta.append ("");
   SX_TEST (beta.size == 0);
   
   // --- remove
   alpha = "a";
   alpha.remove (0);
   SX_TEST (alpha.size == 0);

   alpha = "a";
   alpha.removeElement ('a');
   SX_TEST (alpha.size == 0);

   alpha = "a";
   alpha.removeLast ();
   SX_TEST (alpha.size == 0);
   
   // --- operator =
   alpha = "";
   SX_TEST (alpha.size == 0);
   
   alpha = "alpha";
   beta = "";
   alpha = beta;
   SX_TEST (alpha.size == 0);
   
   alpha = "alpha";
   alpha = emptyArray;
   SX_TEST (alpha.size == 0);

   // --- operator +
   beta = "";
   alpha = "" + beta;
   SX_TEST (alpha.size == 0);

   alpha = beta + "";
   SX_TEST (alpha.size == 0);

   return 0;
}

template<class T>
int invalidJoinT ()
{
   T data;
   SX_TEST (SxString::join(data).size == 0);
   SX_TEST (SxString::join(data, ":").size == 0);
   
   data.append ("");
   SX_TEST (SxString::join(data).size == 0);
   SX_TEST (SxString::join(data, ":").size == 0);

   data.append ("");
   SX_TEST (SxString::join(data).size == 0);
   SX_TEST (SxString::join(data, ":") == ":");

   data.append ("a");
   SX_TEST (SxString::join(data) == "a");
   SX_TEST (SxString::join(data, ":") == "::a");
   
   T dataA;
   dataA.append ("a");
   SX_TEST (SxString::join(dataA) == "a");
   SX_TEST (SxString::join(dataA, ":") == "a");

   T dataAB;
   dataAB.append ("a");
   dataAB.append ("b");
   SX_TEST (SxString::join(dataAB) == "ab");
   SX_TEST (SxString::join(dataAB, ":") == "a:b");
   
   // --- words
   T words;
   words.append ("this");
   words.append ("is");
   words.append ("a");
   words.append ("test");
   words.append ("program");
   SX_TEST (SxString::join(words) == "thisisatestprogram");
   SX_TEST (SxString::join(words, ":") == "this:is:a:test:program");
   
   return 0;
}

int invalidJoin ()
{
   invalidJoinT<SxArray<SxString> > ();
   invalidJoinT<SxList<SxString> > ();

   return 0;
}

int invalidGetBlocks ()
{
   SxString orig = "BEGIN\nfoo\nEND";
   SxList<SxArray<ssize_t> > blocks = orig.getBlockIdx ("BEGIN", "END", false);
   SX_TEST (blocks.getSize () == 1);
   SxArray<ssize_t> arr = blocks(0);
   SxString substr = orig.subString (arr(0), arr(1));
   SxList<SxString> direct = orig.getBlocks ("BEGIN", "END", false);
   SX_TEST (direct.getSize () == 1);
   SX_TEST (arr(0) == 6);
   SX_TEST (arr(1) == 8);
   SX_TEST (substr == "foo");

   direct = orig.getBlocks ("BEGIN", "END", true);
   SX_TEST (direct.getSize () == 1);
   SX_TEST (direct(0) == orig);

   SxString bar = "BEGINbarEND"; // without newlines
   blocks = bar.getBlockIdx ("BEGIN", "END", false);
   SX_TEST (blocks.getSize () == 1);
   SxArray<ssize_t> bArr = blocks(0);
   SX_TEST (bArr(0) == 5);
   SX_TEST (bArr(1) == 7);
   SxList<SxString> bDirect = bar.getBlocks ("BEGIN", "END", true);
   SX_TEST (bDirect.getSize () == 1);
   SX_TEST (bDirect(0) == bar);

   // IMPLEMENTME: more complicated tests of getBlockIdx (), getBlocks (),
   // substBlocks ()!

   return 0;
}

static void inspectRawBytes (const char *label, const char *str)
{
   if (str == NULL)  str = "";
   cout << "inspectRawBytes (): " << label << " - *";
   char ch;
   const char *s = str;
   while ( (ch = *s++) != '\0' )  {
      char vis = ( (ch != ' ') ? ch : '_' ); // ((char) 183) );
      cout << vis;
   }
   cout << "*, strlen=" << ::strlen (str);
   while ( (ch = *str++) != '\0' )  {
      char vis = ( (ch != ' ') ? ch : '_' ); // ((char) 183) );
         // (visible replacement for space charater; "middot")
      cout << "; *" << vis << "*/";
      printf("%u", (uint8_t) ch);
   }
   cout << endl;
}

static void inspectAny (const char *label, const SxString &obj)
{
   printf ("inspectAny (): %s; *%s*; isUnicode=%u; nBytes=%ld; nChars=%ld\n",
      label, obj.getElems (), obj.isUnicode (), obj.getNBytes (),
      obj.getSize ());
   inspectRawBytes (label, obj.elements);
#if !INCOMPLETE
   obj.inspect (label, SxString::ifAll);
#endif
}

static inline void inspectUnicode (const char *label, const SxString &obj)
{
   inspectAny (label, obj); // (same implementation, different semantics)
}

// special German characters; we must define strings with 8-bit ASCII
// characters carefully because compilers mangle such characters and ruin tests
#define ASC8(value) ((char) (value))
#define AUML (ASC8(228))
#define OUML (ASC8(246))
#define UUML (ASC8(252))
#define SZLIG (ASC8(223))

// Japanese (Hiragana/Katakana) characters
#define YA (0x3084) // "や"
#define ZO (0x305e) // "ぞ"
#define DE (0x30c7) // "デ"

// Chinese characters
#define FORCE (0x529b) // "力"
#define PLAN  (0x4f01) // "企"; plan a project
#define BALL_BMP (0x7403) // "球"
#define BALL (0x2f801) // "丸"; another character for "ball" is the above U+7403

// various
#define EURO (0x20ac)

#define SPSTR " " // one space
#define SPSTR3 SPSTR SPSTR SPSTR

static SxString unicodeFromAscii (const char *str)
{
   return SxString::asciiToUnicode (str);
}

static SxString unicodeFromUtf8 (const char *str)
{
   return SxString::unicodeFromUtf8 (str);
}

int invalidUnicode ()
{
   SxConstChar foo("abcdef", 6, false);
   SxConstChar::Iterator it = foo.begin ();

   /* char a = 'ä';
   SxString aobj = SxString (a);
   inspectAny ("aobj", aobj); */

   static const char strUmlaut[] = { 'a', 'b', 'c', ' ', AUML, OUML, UUML, ' ',
      SZLIG, 0 };
   // static const char strUmlaut[] = "abc äöü ß";
   SxString s1 = unicodeFromAscii (NULL);
   inspectRawBytes ("strUmlaut", strUmlaut);
   SxString s2 = unicodeFromAscii (strUmlaut);
   inspectUnicode ("s2 A", s2);
   SxString s3 = unicodeFromAscii (NULL);
   SxString s4 = s2.toUpper ();
   SxString s5 = s4.toLower ();
   SxString s6 = unicodeFromUtf8 (s2.getElems ());
   s3.append ("xyz ");
   s3.append ((int) ((uint8_t) AUML), 5);
   inspectUnicode ("s3", s3);
   s3.append (" #");
   // inspectUnicode ("s1", s1);
   // inspectUnicode ("s2", s2);
   // inspectUnicode ("s3", s3);
   inspectUnicode ("s4", s4);
   // inspectUnicode ("s5", s5);
   // inspectUnicode ("s6", s6);
   SX_TEST (s1.getNBytes () == 0);
   SX_TEST (s2.getNBytes () == 9 + 4);
   inspectUnicode ("s2 B", s2);
   SX_TEST (s2.getSize () == 9);
   SX_TEST (s3.getSize () == 4 + 5 * 1 + 2);
   SX_TEST (s3.getNBytes () == 4 + 5 * 2 + 2);
   SX_TEST (s4.getNBytes () == s2.getNBytes ());
   SX_TEST (s5 == s2);
   SX_TEST (s6.getNBytes () == s2.getNBytes ());
   SX_TEST (s6.getSize () == s2.getSize ());
   SX_TEST (s6 == s2);

   static const char strH1[] = { AUML, 0 }, strH2[] = { OUML, 0 };
   SxString h1 = unicodeFromAscii (strH1);
   SxString h1double = SxString (h1, h1);
   SxString h2 = unicodeFromAscii (strH2);

   cout << std::endl;
   SxString h3 = s2.substitute (h1, h2);
   inspectUnicode ("h1", h1);
   inspectUnicode ("h2", h2);
   inspectUnicode ("s2", s2);
   inspectUnicode ("h3", h3);

   cout << std::endl;
   inspectUnicode ("s3", s3);
   SxString h4 = s3.substitute (h1, h2);
   inspectUnicode ("h4", h4);

   SxString r0 = "to be replaced";
   // inspectAny ("r0_1", r0);
   r0.replace (s2);
   // inspectAny ("r0_2", r0);
   SX_TEST (r0 == s2);
   SX_TEST (r0.getNBytes () == s2.getNBytes ());
   SX_TEST (r0.getSize () == s2.getSize ());

   SxString b1 = s3.subString (0);
   SxString b2 = s3.subString (5, 6);
   SxString b3 = s3.subString (10);
   SxString b4 = s3.subString (1, 4);
   SxString b5 = s3.left (' ');
   SxString b6 = s3.right (' ').right (' ');
   SxString b7 = s3.left (h1);
   // inspectUnicode ("b1", b1);
   // inspectUnicode ("b2", b2);
   // inspectUnicode ("b3", b3);
   // inspectUnicode ("b4", b4);
   // inspectUnicode ("b5", b5);
   // inspectUnicode ("b6", b6);
   // inspectUnicode ("b7", b7);
   SX_TEST (b1.getSize () == s3.getSize ());
   inspectUnicode ("b2", b2);
   SX_TEST (b2.getSize () == 2);
   SX_TEST (b2.getNBytes () == 4);
   SX_TEST (b3.getSize () == 1);
   SX_TEST (b3.getNBytes () == 1);
   SX_TEST (b4.head (3) == "yz ");
   SX_TEST (b4.getSize () == 4);
   SX_TEST (b4.getNBytes () == 5);
   SX_TEST (b4.tail (1).getNBytes () == 2);
   SX_TEST (b5.getNBytes () == 3);
   SX_TEST (b6.getNBytes () == 1);

   SX_TEST (s3.find (h1) == 4);
   SX_TEST (s3.find (h1, 6) == 6);
   SX_TEST (s3.findLast (h1) == 8);

   SxString a1 = unicodeFromAscii (NULL);
   inspectUnicode ("a1-null", a1);
   a1.append (SPSTR3 "foo" SPSTR SPSTR);
   inspectUnicode ("a1-foo", a1);
   inspectUnicode ("h1", h1);
   a1.append (h1);
   a1.append ("bar" SPSTR3);
   a1.append (h2);
   a1.append (SPSTR3);
   a1.append ("#" SPSTR3);
   // inspectUnicode ("a1", a1);
   SxString a2 = a1.adjustRight ();
   inspectUnicode ("a1", a1);
   inspectUnicode ("a2", a2);
   SX_TEST (a1.getNBytes () == a2.getNBytes ());
   SX_TEST (a1.getNBytes () == 21 * 1 + 2 * 2);
   SX_TEST (a1.getSize () == a2.getSize ());
   SxString a3 = a2.subString (20, 20);
   inspectUnicode ("a3", a3);
   SX_TEST (a3 == h2);

   cout << std::endl;
   SxString a4 = unicodeFromAscii (NULL);
   inspectUnicode ("a4-null", a4);
   a4.append (YA);
   inspectUnicode ("a4-YA", a4);
   a4.append (' ');
   inspectUnicode ("a4-YA2", a4);
   a4.append (ZO);
   inspectUnicode ("a4-ZO", a4);
   a4.append (' ');
   a4.append (DE);
   inspectUnicode ("a4-DE", a4);
   a4.append (' ');
   inspectUnicode ("a4", a4);
   a4.append (FORCE);
   inspectUnicode ("a4-FO", a4);
   a4.append (' ');
   a4.append (PLAN);
   inspectUnicode ("a4-PL", a4);
   a4.append (' ');
   a4.append (BALL);
   inspectUnicode ("a4", a4);
   SX_TEST (a4.getNBytes () == 24);

   SX_TEST (s3.contains (h1) == 5);
   SX_TEST (s3.containsWhole (h1) == 5);
   SX_TEST (s3.contains (h1double) == 4);
   SX_TEST (s3.containsWhole (h1double) == 2);
   SX_TEST (s3.contains ("z") == 1);
   SX_TEST (s3.contains ("#") == 1);

   SxList<SxString> tok = s3.tokenize (h1, true);
   /* cout << "tokenized s3 at auml:" << std::endl;
   for (ssize_t idx = 0; idx < tok.getSize (); idx++)  {
      char buf[100];
      SxString::sprintf (buf, "tok%ld", idx);
      inspectAny (buf, tok(idx));
   } */
   SX_TEST (tok.getSize () == 6);

   tok = s3.tokenize (SPSTR, true);
   /* cout << "tokenized s3 at space:" << std::endl;
   for (ssize_t idx = 0; idx < tok.getSize (); idx++)  {
      char buf[100];
      SxString::sprintf (buf, "tok%ld", idx);
      inspectAny (buf, tok(idx));
   } */
   SX_TEST (tok.getSize () == 3);

   tok = s2.tokenize (SPSTR, true);
   /* cout << "tokenized s2 at space:" << std::endl;
   for (ssize_t idx = 0; idx < tok.getSize (); idx++)  {
      char buf[100];
      SxString::sprintf (buf, "tok%ld", idx);
      inspectAny (buf, tok(idx));
   } */
   SX_TEST (tok.getSize () == 3);

   SX_TEST (unicodeFromAscii ("#").stripComments ().size ==
      0);
   SxString c1 = unicodeFromAscii ("foo#bar\nboo\ncoo\ndoo#x");
   // inspectUnicode ("c1", c1);
   SxString c2 = c1.stripComments ();
   // inspectUnicode ("c2", c2);
   SX_TEST (c2.getNBytes () == 15);
   SX_TEST (c2.getSize () == 15);

#ifdef WIN32
   // Use only characters in the BMP (U+0000..U+ffff) in these tests! Other
   // characters don't fit into one 16-bit wchar_t on Windows.
   cout << endl << "Windows-specific tests (wchar_t):" << endl;
   static const wchar_t wstrHello[] = { 'H', 'e', 'l', 'l', 'o', '\0' };
   SxString wHello(wstrHello);
   cout << "wstrHello[]: *" << wstrHello << "*; *" << wHello << "*" << endl;

   SxString w1 = "foo bar";
   w1.replace (wstrHello);
   inspectUnicode ("w1", w1);
   SX_TEST (w1.getNBytes () == 5);
   SX_TEST (w1.getSize () == 5);

   SxArray<wchar_t> w2Obj = w1.toWChars ();
   wchar_t *w2 = w2Obj.elements;
   cout << "w2";
   wchar_t wch;
   ssize_t wCount = 0;
   if (w2 != NULL)  {
      while ( (wch = *w2++) != 0 )  {
         cout << "; " << wch;
         wCount++;
      }
   }
   cout << endl;
   SX_TEST (wCount == 5);

   static const wchar_t wstrUmlaut[] = { 'a', 'b', 'c', ' ', (uint8_t) AUML,
      (uint8_t) OUML, (uint8_t) UUML, ' ', (uint8_t) SZLIG, 0 };
   SxString w3(wstrUmlaut);
   inspectUnicode ("w3", w3);
   inspectUnicode ("s2", s2);
   SX_TEST (w3.getNBytes () == 9 + 4);
   SX_TEST (w3.getSize () == 9);

   static const wchar_t wstrAsian[] = { YA, ' ', ZO, ' ', DE, ' ', BALL_BMP,
      '\0' };
   SxString w4(wstrAsian);
   inspectUnicode ("w4", w4);
   SX_TEST (w4.getNBytes () == (4 * 3) + (3 * 1));
      // four "complicated" characters, three spaces
   SX_TEST (w4.getSize () == 7);
#endif

   SxString filename = SxString::asciiToUnicode ("sxstring-testfile-");
   filename.append (YA);
   filename.append ('-');
   filename.append (FORCE);
   SxString filetext = a4;
   filetext.append (" foo\n");
   filetext.write (filename);
   SxString filetext2 = SxString::readBinary (filename, -1);
   SxString filetext3 = SxString::unicodeFromUtf8 (filetext2.getElems (),
      filetext2.getNBytes ());
   inspectUnicode ("filetext", filetext);
   inspectUnicode ("filetext2", filetext2);
   inspectUnicode ("filetext3", filetext3);
   SX_TEST (filetext == filetext3);

   return 0;
}

#define SX_CALL_TESTFUNC(id) puts ("Calling " #id " ()"); if (id ())  return 1;

int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   cli.finalize ();

   try {
   SX_CALL_TESTFUNC (invalidBist)
   SX_CALL_TESTFUNC (invalidResize)
   SX_CALL_TESTFUNC (invalidConstructor)
   SX_CALL_TESTFUNC (invalidSet)
   SX_CALL_TESTFUNC (invalidFormat)
   SX_CALL_TESTFUNC (invalidNumber)
   SX_CALL_TESTFUNC (invalidNumber2)
   SX_CALL_TESTFUNC (invalidNumber3)
   SX_CALL_TESTFUNC (invalidNumber4)

   SX_CALL_TESTFUNC (invalidIndex)
   SX_CALL_TESTFUNC (invalidSubString)
   SX_CALL_TESTFUNC (invalidReplace)
   SX_CALL_TESTFUNC (invalidSubstitute)

   SX_CALL_TESTFUNC (invalidFind)
   SX_CALL_TESTFUNC (invalidContains)

   SX_CALL_TESTFUNC (invalidInsert)
   SX_CALL_TESTFUNC (invalidPrepend)
   SX_CALL_TESTFUNC (invalidAppend)
   SX_CALL_TESTFUNC (invalidRemove)

   SX_CALL_TESTFUNC (invalidCase)
   SX_CALL_TESTFUNC (invalidPlus)
   SX_CALL_TESTFUNC (invalidCompare)

   SX_CALL_TESTFUNC (invalidRegexp)
   SX_CALL_TESTFUNC (invalidSed)
   SX_CALL_TESTFUNC (invalidWhiteSpace)
   SX_CALL_TESTFUNC (invalidEmptyString)
   SX_CALL_TESTFUNC (invalidJoin)
   SX_CALL_TESTFUNC (invalidGetBlocks)

   SX_CALL_TESTFUNC (invalidUnicode)

   sxprintf ("test ok\n");

   } catch (SxException e) {
      e.print (true);
      SX_EXIT;
   }

   return 0;
}
