#include <stdio.h>
#include <SxPtr.h>
#include <SxList.h>
#include <SxCLI.h>

class A
{
   public:
      A () { }
      virtual ~A ()  { }
      virtual void foo () const { printf ("a\n"); }
};


class B : public A
{
   public:
      B () : A () { }
      virtual ~B ()  { }
      virtual void foo () const { printf ("b\n"); }
      
      
};


SxList<SxPtr<A> >  abc ()
{
   SxList<SxPtr<A> >  list;
   list << SxPtr<A>::create ();
   list << SxPtr<B>::create ();
   return list;
}

SxList<SxPtr<A> >  def ()
{
   SxList<SxPtr<A> >  list;
   list << SxPtr<B>::create ();
   list << SxPtr<A>::create ();
   return list;
}

int main (int argc, char **argv)
{
   // --- for '--memcheck'
   SxCLI cli(argc, argv);
   cli.finalize ();

   SxPtr<A> a1 = SxPtr<A>::create ();
   SxPtr<B> b1;
   b1 = SxPtr<B>::create ();

   a1->foo ();  // prints "a"
   b1->foo ();  // prints "b"

   // --- overwrite a1
   a1 = b1;
   a1->foo ();  // prints "b"

   // --- using derived classes
   SxPtr<A> a2 = SxPtr<B>::create ();
   a2->foo ();  // prints "b"
   
   // --- using constness
   SxConstPtr<A> a3 = SxConstPtr<A>::create ();
   SxConstPtr<A> a4 = SxConstPtr<B>::create ();
   a3->foo ();  // prints "a"
   a4->foo ();  // prints "b"

   // --- typecasts non-const to const
   a3 = a1;
   a4 = b1;
   a3->foo ();  // prints "a"
   a4->foo ();  // prints "b"
   
   // --- using lists and auto pointers
   SxList<SxPtr<A> >::Iterator it;

   SxList<SxPtr<A> >  list = abc ();
   for (it = list.begin(); it != list.end(); ++it)
      (*it)->foo ();  // prints "a b"

   list = def ();
   for (it = list.begin(); it != list.end(); ++it)
      (*it)->foo ();  // printd "b a"
  
   return 0;
}


