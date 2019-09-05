#include <SxProcess.h>
#include <SxUtil.h>
#include <SxTime.h>

int main ()
{
   // --- execute an external program
   int status = 0;
#  ifdef WIN32
   SxString cmd = "cmd.exe";
#  else
   SxString cmd = "xterm";
#  endif
   
   // First we test the environment by openening a new terminal
   // where possible child-processes can be started
   // If sucessfull the program should only go on when all child processes have been closed
   // In addition, a termination of the host process should also terminate all child processes
   SxProcess proc (cmd);
   cout << "Start process" << endl;
   try                    { status = proc.start ();  }
   catch (SxException e)  { e.print (); }
   cout << "Process ended with exit status " << status << endl;

   // The second part creates a new terminal in the no hangup mode.
   // The host process will terminate directly, but the child processes 
   // should remain open
   cmd += " &";
   proc = SxProcess (cmd);
   cout << "Start process in nohangup mode" << endl;
   try { status = proc.start(); SxTime::msleep (10000);}
   catch (SxException e) { e.print(); }
   cout << "Process ended with exit status " << status << endl;
   cout << "Check that child processes are still open." << endl;

   // The third part is a Windows only feature.
   // We open an application "notepad", suspend it and resume it after 
   // a specified amount of time

#  ifdef WIN32
   cmd = "notepad.exe";
   proc = SxProcess (cmd);
   try {
      proc.run ();
      SxTime::msleep (10000);
      cout << "Suspend application" << endl;
      proc.suspend ();
      SxTime::msleep (20000);
      cout << "Resume application" << endl;
      proc.resume ();
      proc.wait ();
   } catch (SxException e) {e.print ();}
#  endif

   // The fourth part checks if the process is automatically 
   // terminated after specific amount of time.
#  ifdef WIN32
   cmd = "notepad.exe";
   proc = SxProcess(cmd);
   try {
      // set process lifetime to 5 seconds
      proc.setTimeout (5);
      proc.run ();
      cout << "Application started. Waiting for timeout kill." << endl;
      proc.wait ();
   } catch (SxException e) {e.print ();}
#  endif

   return 0;
}
