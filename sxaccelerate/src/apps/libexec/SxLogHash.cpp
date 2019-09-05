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


#include <SxCLI.h>
#include <SxString.h>
#include <SxLog.h>

int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   SxString id = cli.option ("--id", "string",
                             "Identifier of the shared object").toString ();
   cli.finalize ();

   SX_LOG ("aaabbbccc");
   SX_DBG_MSG ("**** DEBUG MSG ****");
   //SxLog::enable ("sxutil");
   SxString a;

   cout << SX_FILE << std::endl;
   cout << SxLog::djb2 (id.ascii()) << endl;
   return 0;
}
