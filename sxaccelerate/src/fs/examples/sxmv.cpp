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
   SxString defaultSrc ("defaultSrc");
   SxString defaultDest ("defaultDest");
   SxString srcStr = cli.option ("-s|--src",
                                 "source of the mv-operation"
                                ).toString (defaultSrc.ascii ());
   SxString destStr = cli.option ("-d|--dest",
                                 "dest of the mv-operation"
                                ).toString (defaultDest.ascii ());
   cli.finalize ();

   SxFileInfo src (srcStr);
   SxFileInfo dest (destStr);
   try  {
      SxFSAction::mv (src, dest);
   } catch (SxException ex)  {
      ex.print ();
      ret = 1;
   }
   return ret;
}
