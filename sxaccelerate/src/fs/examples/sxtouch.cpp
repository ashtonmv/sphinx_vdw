// --- Including header-files
#include <iostream>
#include <SxCLI.h>
#include <SxString.h>
#include <SxFile.h>
#include <SxFSAction.h>

int main (int argc, char **argv)
{
   int ret = 0;
   // --- Using the command line interface to retrieve the program arguments
   SxCLI cli (argc, argv);
   SxString defaultTarget ("defaultTargetFile");
   SxString target = cli.option ("-t|--target",
                                 "target of the touch-operation"
                                ).toString (defaultTarget.ascii ());
   cli.finalize ();
   
   SxFile targetFile (target);
   SX_CHECK (!targetFile.exists ());
   try  {
      SxFSAction::touch (targetFile);
   } catch (SxException ex)  {
      ex.print ();
      ret = 1;
   }
   return ret;
}
