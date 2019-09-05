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
   SxString targetStr = cli.option ("-t|--target",
                                 "target of the rm-operation"
                                ).toString ("");
   bool recursive = cli.option ("-r",
                                "Enables recursive removal"
                               ).toBool ();
   cli.finalize ();

   if (targetStr != SxString ())  {
      SxFileInfo target (targetStr);
      try  {
         if (recursive) {
            SxFSAction::rm_r (target);
         } else  {
            SxFSAction::rm (target);
         }
      } catch (SxException ex)  {
         ex.print ();
         ret = 1;
      }
   } else  {
      std::cerr << "No target was specified." << std::endl;
      ret = 1;
   }
   return ret;
}
