// --- Including header-files
#include <iostream>
#include <SxCLI.h>
#include <SxString.h>
#include <SxFile.h>
#include <SxFSAction.h>
#include <ctime>

void outp (const SxFileInfo &curRes)
{
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

int main (int argc, char **argv)
{
   int ret = 0;
   // --- Using the command line interface to retrieve the program arguments
   SxCLI cli (argc, argv);
   SxString targetStr = cli.option ("--target",
                                 "target of the ls-operation"
                                ).toString ("");
   bool sortByTime = cli.option ("-t",
                                 "target of the ls-operation"
                                ).toBool ();
   bool sortBySize = cli.option ("-s",
                                 "target of the ls-operation"
                                ).toBool ();
/*   bool recursive = cli.option ("-r",
                                "Enables recursive removal"
                               ).toBool ();*/
   cli.finalize ();

   if (targetStr != SxString ())  {
      SxFileInfo target (targetStr);
      try  {
         //if (recursive) {
         //   SxFSAction::ls_r (target);
         //} else  {
         SxList<SxFISortedByTime> resultsSortedByTime;
         SxList<SxFISortedBySize> resultsSortedBySize;
         SxList<SxFileInfo> resultsUnsorted;
         if (sortByTime && sortBySize)  {
            cerr << "You can either sort by time or by size but ";
            cerr << "not both at once.";
            return 2;
         } else if (sortByTime)  {
            std::cout << "sorted by time..."<< std::endl;//TEST
            resultsSortedByTime = SxFSAction::ls_t (target);
            SxList<SxFISortedByTime>::Iterator itResults;
            for (itResults = resultsSortedByTime.begin (); 
                 itResults != resultsSortedByTime.end (); 
                 ++itResults)
            {
               SxFileInfo const &curRes = (*itResults);
               outp (curRes);
            }        
         } else if (sortBySize)  {
            std::cout << "sorted by size..."<< std::endl;//TEST
            resultsSortedBySize = SxFSAction::ls_S (target);
            SxList<SxFISortedBySize>::Iterator itResults;
            for (itResults = resultsSortedBySize.begin (); 
                 itResults != resultsSortedBySize.end (); 
                 ++itResults)
            {
               SxFileInfo const &curRes = (*itResults);
               outp (curRes);
            }
         } else  {
            resultsUnsorted = SxFSAction::ls (target);
            SxList<SxFileInfo>::Iterator itResults;
            for (itResults = resultsUnsorted.begin (); 
                 itResults != resultsUnsorted.end (); 
                 ++itResults)
            {
               SxFileInfo const &curRes = (*itResults);
               outp (curRes);
            }
         }
         
         //}
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
