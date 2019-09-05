// --- Including header-files
#include <iostream>
#include <SxCLI.h>
#include <SxString.h>
#include <SxDir.h>
#include <SxFSAction.h>

int main (int argc, char **argv)
{
   int ret = 0;
   // --- Using the command line interface to retrieve the program arguments
   SxCLI cli (argc, argv);
   SxString defaultTarget ("defaultTargetDir");
   SxString target = cli.option ("-t|--target",
                                 "target of the mkdir-operation"
                                ).toString (defaultTarget.ascii ());
   bool enableMkdir_p = cli.option ("-p",
                                    "Causes that target is used with \"mkdir -p\""
                                   ).toBool ();
   cli.finalize ();
   SxDir targetDir (target);
   
   SX_CHECK (!targetDir.exists ());
   try  {
      if (enableMkdir_p)  {
         SxFSAction::mkdir_p (targetDir);
      } else  {
         SxFSAction::mkdir (targetDir);
      }
   } catch (SxException ex)  {
      ex.print ();
      ret = 1;
   }
   return ret;
}
