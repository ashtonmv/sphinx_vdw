#include <SxCLI.h>
#include <SxProcess.h>
#include <SxLogBuf.h>
#include <SxRedirect.h>

int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   cli.finalize ();

   // --- daemonize this process   
   SxProcess::daemon ();

   // --- use logBuf for stdout to preppend current time to each output  
   SxLogBuf logBuf;
   logBuf.setFile ("logfile.log");
   SxRedirect teeOut (std::cout, &logBuf);
   std::streambuf *cerrbuf = std::cerr.rdbuf ();
   std::cerr.rdbuf (std::cout.rdbuf ());
   
   cout << "Daemon started" << endl;
   cerr << "error" << endl;
   
   std::cerr.flush ();
   std::cerr.rdbuf (cerrbuf);
   
   // cat logfile.log
   //    10/31/12 13:47:58: Daemon started
   //    10/31/12 13:47:58: error
   
   return 0;
}
