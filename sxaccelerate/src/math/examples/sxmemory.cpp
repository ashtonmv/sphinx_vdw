#include <SxError.h>
#include <SxStack.h>
#include <SxVector.h>
#include <SxPtr.h>

int main ()
{
   int i;
   // -------------- STACKS
   // simply push things on the stack
   // you can retrieve the pushed elements by pop (last in - first out)
   SxStack<int> stack;
   // push 0,1,2,3,...
   for (i = 0; i < 100000; ++i)
      stack.push(i);

   cout << "stack size: " << stack.getSize () << endl;

   // pop 99999,99998,99997,...
   for (i = 0; i < 100000; ++i)
      if (stack.pop() + i != 99999) { cout << "Error" << endl; }

   cout << "stack size: " << stack.getSize () << endl;

   // But stacks can do more ...
   // Stacks are used to build up arrays when the size is not yet known
   // stacks are very efficient for that (much faster than lists)
   // once the final size is reached, export the data to SxArray or SxVector

   for (i = 0; i < 100000; ++i) stack.push(i);
   cout << "stack size: " << stack.getSize () << endl;
   SxArray<int> array(stack);
   cout << "array size: " << array.getSize () << endl;
   for (i = 0; i < 5; ++i) cout << array(i);
   cout << "..." << array(array.getSize () - 1) << endl;
   
   // Stacks can be used for all kinds of types, e.g. strings
   SxStack<SxString> strStack;

   strStack.push ("-one");
   strStack << "-two"; // same as strStack.push("-two")
   strStack << "-three";
   cout << strStack.pop ();
   cout << strStack.pop ();
   cout << strStack.pop ();
   cout << endl;

   strStack.push ("completely");
   strStack.push ("useless");
   strStack.push ("data");
   // If you really want to discard all the data:
   strStack.removeAll ();

//   // Objects on the stack really die when popped or exported, which makes
//   // SxStack<SxPtr<?> > a useful thing for very big objects.
//   // make an example for int (which of course is not big)
//   SxStack<SxPtr<int> > ptrStack;
//   SxPtr<int> ptr = SxPtr<int>::create (1);
//   
//   // NEVER DO THIS IN REAL CODE! This is just to demonstrate the
//   // life-time of pointers on the stack
//   SxRefCounter ref = ptr.refCounter; // reference counter of recently created object
//
//   // now 1 SxPtr points to object: ptr
//   cout << "nRef = " << *ref << " should be 0" << endl;
//   ptrStack << ptr;
//   // now 2 SxPtrs point to object: ptr and a stack element
//   cout << "nRef = " << *ref << " should be 1" << endl;
//   ptr = SxPtr<int> ();
//   // now 1 SxPtr points to object: stack
//   cout << "nRef = " << *ref << " should be 0" << endl;
//   ptr = ptrStack.pop ();
//   // now 1 SxPtr points to object: ptr
//   cout << "nRef = " << *ref << " should be 0" << endl;
//   ptr = SxPtr<int> ();
//   // now object is dead and ref is invalid. You can get a nice
//   // segmentation fault when using ref now. That's why you never do what I
//   // did here.
   
}
