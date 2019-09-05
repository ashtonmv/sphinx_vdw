#include <SxProcess.h>

class MyProc : public SxProcess
{
   public:
      MyProc (const SxString         &cmdLine,
              const SxList<SxString> &cmdArgs);
      virtual SxString writeToStdin ()
      {
         return "hello";
      }
};

MyProc::MyProc (const SxString         &cmdLine,
                const SxList<SxString> &cmdArgs)
   : SxProcess (cmdLine, cmdArgs)
{
   setCommunication (SxProcess::Stdin);
}
 
// write "hello world!" to sxprocstdin.txt
int main ()
{
   SxList<SxString> cmdArgs;
   cmdArgs << "-c" << "cat > sxprocstdin.txt";
   
   MyProc proc ("/bin/bash", cmdArgs);

   proc.run ();
   
   // --- continue writing to STDIN
   proc.write (" world");
   proc.write ("!");
   
   int exitStatus = proc.wait ();
   
   return exitStatus;
}
