// --- Including header-files
#include <iostream>
#include <SxCLI.h>
#include <SxString.h>
#include <SxFileInfo.h>
#include <SxSymLink.h>
#include <SxFSAction.h>

int main (int argc, char **argv)
{
   int ret = 0;
   // --- Using the command line interface to retrieve the program arguments
   SxCLI cli (argc, argv);
   SxString permStr = cli.option ("-p|--perm",
                                 "permission"
                                ).toString ("700");
   SxString targetStr = cli.option ("-t|--target",
                                 "dest of the mv-operation"
                                ).toString ("");
   cli.finalize ();

   if (permStr.isInt ())  {
      try  {
         int perm (permStr.toInt ());
         SxFileInfo target (targetStr);
         SxFSAction::chmod (perm, target);
      } catch (SxException ex)  {
         ex.print ();
         ret = 1;
      }
   } else  {
      try  {
         SxFileInfo target (targetStr);
         SxFSAction::chmod (permStr, target);
      } catch (SxException ex)  {
         ex.print ();
         ret = 1;
      }
/*      std::cerr << "Error: The passed permission has to be an integer value.";
      std::cerr << std::endl;
      ret = 1;*/
   }
   return ret;
}
