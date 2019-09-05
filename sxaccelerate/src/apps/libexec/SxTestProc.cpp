#include <SxCLI.h>
#include <SxTime.h>
#include <iostream>
#include <string>

int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   SxString coutStr = cli.option ("--stdout", "string",
                        "Print provided message string to STDOUT "
                        "or redirect input from STDIN when '-' is "
                        "provided."
   ).toString ("");
   SxString cerrStr = cli.option ("--stderr", "string",
                        "Print provided message string to STDERR "
                        "or redirect input from STDIN when '-' is "
                        "provided."
   ).toString ("");
   int runtimeSec   = cli.option ("--runtime", "int",
                        "Run (sleep) for the provided number of seconds "
                        "before generating any output."
   ).toInt (0);
   int exitCode     = cli.option ("--exit", "int",
                        "Return with the specified exit code."
   ).toInt (0);

   cli.finalize ();

   cout << "sxtestproc started." << endl;

   if (runtimeSec > 0)  SxTime::msleep (static_cast<unsigned int>(runtimeSec) * 1000U);

   if (coutStr == "-" || cerrStr == "-")  {
      SxString input;
      SxList<SxString> lines;
      for (std::string stdStr; std::getline (std::cin, stdStr); )  {
         lines << stdStr.c_str ();
      }
      input = SxString::join (lines, "\n");

      if (coutStr == "-")  coutStr = input;
      if (cerrStr == "-")  cerrStr = input;
   }

   if (coutStr != "")  cout << coutStr << endl;
   if (cerrStr != "")  cerr << cerrStr << endl;

   return exitCode;
}
