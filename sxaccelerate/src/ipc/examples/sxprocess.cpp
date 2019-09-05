#include <SxProcess.h>
#include <SxUtil.h>

// --- define an own process class (used for the 2nd demonstration step)
class MyLsProc : public SxProcess
{
   public:
      MyLsProc (const SxString &cmdLine) : SxProcess (cmdLine) { }
      virtual void readFromStdout (const SxString &line)
      {
         printf ("process line %s from STDOUT\n", line.ascii());
         fflush (stdout);
      }
};

int main ()
{
   // --- execute an external program
   int status = 0;
#  ifdef WIN32
   SxString cmd = "c:/msys/1.0/bin/ls.exe";
#  else
   SxString cmd = "/bin/ls";
#  endif
   SxProcess ls (cmd + " -la");
   try                    { status = ls.start ();  }
   catch (SxException e)  { e.print (); }

   printf ("-----\n");

   // --- execute external program and read its STDOUT
   MyLsProc lsProc (cmd + " -la");
   lsProc.setCommunication (SxProcess::Stdout);
   try                    { status = lsProc.start ();  }
   catch (SxException e)  { e.print (); }

   MyLsProc lsProc1 = MyLsProc (cmd);
   lsProc1.setCommunication (SxProcess::Stdout);
   try                    { status = lsProc1.start ();  }
   catch (SxException e)  { e.print (); }

   return 0;
}
