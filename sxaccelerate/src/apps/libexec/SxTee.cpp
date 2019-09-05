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
#include <SxProcess.h>
#include <SxEmergency.h>
#include <SxUtil.h>
#include <SxCLI.h>
#ifndef WIN32
#  include <sys/wait.h>
#endif

class SxTee : public SxProcess, public SxEmergency
{
   public:
      SxTee (const SxString &logfile, const SxString &cmdIn, 
             const SxList<SxString> &args);

      virtual void readFromStdout (const SxString &);
      virtual void readFromStderr (const SxString &);

      virtual void dumpData (int);

      FILE *fp;

      ~SxTee ();

};

SxTee::SxTee (const SxString &logfile, const SxString &cmdIn, 
              const SxList<SxString> &args)
   : SxProcess (), fp(NULL)
{
   // open output
   fp = fopen( logfile.ascii (), "w");
   if (!fp)  {
      cout << "Cannot open file '" << logfile << "' for writing.\n";
      SX_QUIT;
   }
   setCommand (cmdIn);
   setArguments (args);
   setCommunication (Stdout);
   lineWise = false;
}

SxTee::~SxTee ()
{
   if (fp) fclose(fp);
}

void SxTee::readFromStderr (const SxString &line)
{
   fprintf(fp, "%s", line.ascii ());
   fprintf(stderr, "%s", line.ascii ());
   fflush(fp);
   fflush(stderr);
}

void SxTee::readFromStdout (const SxString &line)
{
   fprintf(fp, "%s", line.ascii ());
   printf("%s", line.ascii ());
   fflush(fp);
   fflush(stdout);
}

void SxTee::dumpData(int signum)
{
   if (!processRunning) return;
   // SIGNAL handling: kill child process with signum
   softTerminate (signum);
   // TODO: read out child's output buffers and dump them
}

int main (int argc, char **argv)
{
// initSPHInXCOM ();
   SxCLI cli(argc, argv);
   cli.setLoopSeparator ("-e"); // abuse looping feature
   cli.stickyDefault = false;
   SxString logFile = cli.option("-o", "file", "log file name").toString ();
   bool join = cli.option ("-j|--join","join output and error channels")
               .toBool ();
   // declared, but not parsed here, parsed manually below
   cli.option ("-e", "command", "command to execute").required ();
   if (!cli.error && cli.last ().exists ())  {
      cout << "-e option must be given after -o option!" << endl;
      cli.setError ();
   }
   // now we finish CLI parsing, the rest is manual
   cli.finalize ();
   
   if (cli.arguments.getSize () < 1)  {
      cout << "SxTee: Missing -e option" << endl;
      cli.printUsage ();
      SX_QUIT;
   }
   
   cli.looping ();
   // parse SxCLI option
   SxString cmd = cli.option("-e","command","command to execute").toString ();
   if (cli.error) { SX_QUIT; }
   
   if (cmd == argv[0])  {
      cerr << "Self-tee detected. Stopping to prevent infinite loop." << endl;
      SX_QUIT;
   }
   
   // join higher loops
   for (int iLoop = 2; iLoop < cli.arguments.getSize (); ++iLoop)
      cli.arguments(1).append (cli.arguments(iLoop));

   SxTee tee(logFile, cmd, cli.arguments(1));
   if (join) 
      tee.setCommunication (SxProcess::Channels(  SxProcess::Stdout 
                                                | SxProcess::Stderr));
   
   // sxtee's internal error exit code 
   int status = -126;
   try {
      status = tee.start ();
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
   return status;
}

