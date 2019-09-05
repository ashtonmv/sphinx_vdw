// --- Including header-files
#include <iostream>
#include <SxCLI.h>
#include <SxString.h>
#include <SxFile.h>
#include <SxFSAction.h>
#include <SxTimer.h>
#include <ctime>

int main (int argc, char **argv)
{
   int ret = 0;
   // --- Using the command line interface to retrieve the program arguments
   SxCLI cli (argc, argv);
   SxString targetStr = cli.option ("--target",
                                    "target of the ls-operation"
                                   ).toString ("");
   cli.finalize ();
   try  {
      SxFSAction::cd (SxDir ().getAbsPath ());
   } catch (SxException ex)  {
	   ex.print ();
   }
   SxTimer timer = SxTimer (1, true);
   double elapsed = 0.;
   if (targetStr != SxString ())  {
	   std::cout << "Target: \""<< targetStr << "\"" << std::endl;//TEST
	   std::cout << "Current directory: \""<< SxFSAction::pwd ().getAbsPath () << "\"" << std::endl;//TEST
      try  {
         SxList<SxFileInfo> results;
         timer.start ();
         results = SxFSAction::find (targetStr);
		 std::cout << "Number of found results: "<< results.getSize () << std::endl;//TEST
         timer.stop ();
         elapsed = timer.getTime ();
         SxList<SxFileInfo>::Iterator itResults;
         for (itResults = results.begin ();
              itResults != results.end ();
              ++itResults)
         {
            SxFileInfo const &curRes = (*itResults);
            if (curRes.isSymLink ())  {
               std::cout <<  "l";
            } else if (curRes.isFile ())  {
               std::cout <<  "-";
            } else if (curRes.isDir ())  {
               std::cout <<  "d";
            } else  {
               SX_EXIT;
            }
            std::cout << curRes.getPerms () << " ";
            std::cout <<  curRes.getUID () << " ";
            std::cout <<  curRes.getGID () << " ";
            std::cout <<  curRes.getSize () << " ";
            tm *modTime;
            time_t lastMod = curRes.lastModified ();
            modTime = localtime (&lastMod);
            std::cout << modTime->tm_year+1900 << "-";
            std::cout << modTime->tm_mon+1 << "-";
            std::cout << modTime->tm_mday  << " ";
            std::cout << modTime->tm_hour <<":" << modTime->tm_min << ":";
            std::cout << modTime->tm_sec << " ";
            std::cout <<  curRes.getAbsPath ();
            if (curRes.isSymLink ())  {
               std::cout << " -> " <<  SxSymLink (curRes).getTarget ();
            }
            std::cout << std::endl;
         }
      } catch (SxException ex)  {
         ex.print ();
         ret = 1;
      }
   } else  {
      std::cerr << "No target was specified." << std::endl;
      ret = 1;
   }
   std::cout << "----------------------" << std::endl;//TEST
   std::cout << "Time elapsed: " << elapsed << std::endl;//TEST
   std::cout << "----------------------" << std::endl;//TEST
   std::cout << "Press any key..." << std::endl;//TEST
   std::flush (std::cout);
   char c;
   std::cin >> c;
   return ret;
}
