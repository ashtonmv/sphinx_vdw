#include <SxCLI.h>
#include <SxRedirect.h>
#include <SxLogBufThread.h>

class A : public SxSystemThread
{
   public:
      virtual void main ()
      {
         cout << "He" << "llo";
      }
};

class B : public SxSystemThread
{
   public:
      virtual void main ()
      {
         cout << "World" << "!";
      }
};
                                     
int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   cli.finalize ();

   SxLogBufThread logBuf;
   logBuf.setFile ("logfile.log");
   SxRedirect tee (std::cout, &logBuf);
   logBuf.start ();
      
   A a;
   B b;
   
   a.start ();
   b.start ();
   
   a.wait ();
   b.wait ();

   cout << endl;
   
   // --- Output:
   // 11/15/12 11:49:53: World!Hello
   //
   // HeWorldllo! would be also possible
      
   return 0;
}
