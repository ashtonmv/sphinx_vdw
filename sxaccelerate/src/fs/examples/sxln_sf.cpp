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
   SxString defaultTarget ("defaultTargetSymLink.lnk");
   SxString defaultPath ((SxFSAction::pwd ()/"defaultPathSymLink").getAbsPath ());
   SxString target = cli.option ("-t|--target",
                                 "target of the ln-sf-operation"
                                ).toString (defaultTarget.ascii ());
   SxString path = cli.option ("-p|--path",
                                 "path of the ln-sf-operation"
                                ).toString (defaultPath.ascii ());
   cli.finalize ();

   SxSymLink targetSymLink (target);
   SxFileInfo pathSymLink (path);
   SX_CHECK (!targetSymLink.exists ());
   try  {
      SxFSAction::ln_sf (pathSymLink, targetSymLink);
   } catch (SxException ex)  {
      ex.print ();
      ret = 1;
   }
   return ret;
}
