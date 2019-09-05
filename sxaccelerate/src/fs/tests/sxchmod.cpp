// ---------------------------------------------------------------------------
//
//           The general purpose cross platform C/C++ framework
//
//                       S x A c c e l e r a t e
//
//           Home:       https://www.sxlib.de
//           License:    Apache 2
//           Authors:    see src/AUTHORS
//
// ---------------------------------------------------------------------------

#include <SxFSAction.h>
#include <SxCLI.h>

int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   SxString target = cli.option ("--input|-i", "Filename").toString ();
   SxString perms  = cli.option ("--perms|-p", "Permissions").toString ();
   cli.finalize (); 
   
   try  {
      int mode = 0;
      bool error = false;
      mode = perms.toInt (&error);
      if (!error)  {
         SxFSAction::chmod (mode, target);
      }  else  {
         SxFSAction::chmod (perms, target);
      }
   } catch (SxException e)  {
      cout << "ERROR: " << e.getMessage () << endl;
   }
   
   return 0;
}
