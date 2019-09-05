#include <stdio.h>
#include <string.h>
#include <SxError.h>
#include <SxList.h>
#include <SxString.h>
#include <SxFunction.h>
#include <SxCLI.h>

void foo (const SxString &i) { printf ("foo(%s)\n", i.ascii()); }
void fooIso (double val) { printf ("global iso: %g\n", val); }

class A
{
   public:
      A (int i_) : i(i_) { }

      void fooV () { printf ("A::foo()\n"); }
      void foo (const SxString &s) { printf ("A::foo (%s)\nA::i=%d\n", s.ascii(), i); }
      void foo1 (const SxString &s) { printf ("A::foo (%s)\nA::i=%d\n", s.ascii(), i); }
      void fooD (double d) { printf ("A::fooD (%g)\n", d); }
      double fooDD (double d1, double d2) { printf ("A::fooDD (%g,%g)\n", d1, d2); return d1 * d2; }

   protected:
      int i;
};

/* Complicated stuff: multiply derived classes, virtual functions, virtual base classes
class B1 {
   public: int dummyB1;
   virtual void bar (const SxString &s) {
      printf ("B2::bar (%s)\n", s.ascii ());
   }
   virtual ~B1 () {}
};

class B2 {
   public:
   int dummy;
   virtual void foo (const SxString &s) {
      printf ("B2::foo (%s)\n", s.ascii ());
   }
   virtual ~B2 () {}
};

class C1 : virtual public B1
{
   public:
   void foo (const SxString &s) {
      printf ("C1::foo (%s) dummyB1=%i\n", s.ascii (), B1::dummyB1);
   }
   virtual void bar (const SxString &s) {
      printf ("C1::bar (%s)\n", s.ascii ());
   }
   void setC1 (int i) { B1::dummyB1 = i; }
   virtual ~C1 () {}
};

class C2 : virtual public B1, public B2, public C1 {
   public:
   virtual void foo (const SxString &s) {
      printf ("C2::foo (%s) dummyB1=%i\n", s.ascii (), B1::dummyB1);
   }
   virtual void bar (const SxString &s) {
      printf ("C2::bar (%s)\n", s.ascii ());
   }
   void setC2 (int i) { B1::dummyB1 = i; }
   virtual ~C2 () {}
};
*/


int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   cli.finalize ();

   SxPtr<A> a = SxPtr<A>::create (10);
   SxPtr<A> b = SxPtr<A>::create (22);

//   // --- create callbacks with signature "void (*)(const SxString &)"
   SxFunction<void> stringCB0;
   SxFunction<void,const SxString &> stringCB1, stringCB2, stringCB3, stringCB4;

//   stringCB0 = SxFunction<void,void>::create (a, &A::fooV);
//   stringCB0();
   stringCB1 = SxFunction<void,const SxString&>::create (&foo);
   stringCB2 = SxFunction<void,const SxString&>::create (a, &A::foo);
   stringCB3 = SxFunction<void,const SxString&>::create (a, &A::foo1);
   stringCB4 = stringCB2;

   if (stringCB1 != stringCB1)  {
      SX_EXIT;
   }
   if (stringCB1 == stringCB2)  {
      SX_EXIT;
   }

   stringCB2 = stringCB1;

   // --- create callbacks with signature "double (*)(double,double)"
//   SxFunction<double,double,double> doubleCB1;
//   doubleCB1 = SxFunction<double,double,double>::create (a, &A::fooDD);

//   printf ("doubleCB = %g\n", doubleCB1(10.,20.));



// a = SxPtr<A>::create (100);
 printf ("222: %p\n", (void*)stringCB2.receiver.getPtr());

   // ------------------

   stringCB1("hi global");
   stringCB2("hi A");
   stringCB3("hi B");
   stringCB4("hi A'");

   // Call foo(const SxString&) through bound ptr with signature (const char*)
   SxFunction<void,const char *> charCB1
      = SxFunction<void, const char*>::create(a, &A::foo);

   charCB1 ("test");

// SX_EXIT; // TODO: deregister SxPtr <--> SxFunction + valgrind


   // --- execute callbacks

   /*
   SxPtr<C2> c2 = SxPtr<C2>::create ();
   SxFunction<void,const SxString &> bp = SxFunction<void,const SxString&>::create (c2, &B2::foo);
   c2->setC2 (4);
   c2->setC1 (3);
   bp ("expect 3:");
   c2->setC2 (4);
   bp ("expect 4:");
   bp = SxFunction<void,const SxString&>::create (SxPtr<B2> (c2), &B2::foo);
   c2->setC1 (3);
   bp ("expect 3:");
   c2->setC2 (4);
   bp ("expect 4:");
   bp = SxFunction<void,const SxString&>::create (c2, &C1::foo);
   c2->setC1 (3);
   bp ("expect 3:");
   c2->setC2 (4);
   bp ("expect 4:");
   bp = SxFunction<void,const SxString&>::create (SxPtr<B1>(c2), &B1::bar);
   c2->setC1 (3);
   bp ("expect 3:");
   c2->setC2 (4);
   bp ("expect 4:");
   bp = SxFunction<void,const SxString&>::create (c2, &C1::bar);
   c2->setC1 (3);
   bp ("expect 3:");
   c2->setC2 (4);
   bp ("expect 4:");
   */
   
   return 0;
}
