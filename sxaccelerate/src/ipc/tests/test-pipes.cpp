// A test program for the standard stream pipe implementation in the class
// SxProcess.

#include <SxUtil.h>
#include <SxProcess.h>
#include <unistd.h>
#ifndef WIN32
#include <sys/select.h>
#endif

// The following symbols are missing in the Windows header files:
#ifndef STDIN_FILENO
#define STDIN_FILENO  0
#endif
#ifndef STDOUT_FILENO
#define STDOUT_FILENO 1
#endif
#ifndef STDERR_FILENO
#define STDERR_FILENO 2
#endif

SxString progName; // name of executable program file

void testLsSimple ()
{
   cout << "TESTING: simple 'ls/ipconfig' command" << endl;

   int status = 0;
#ifdef WIN32
   SxProcess process("ipconfig");
#else
   SxProcess process("ls -l -s");
#endif
   process.setKillTree (false);
   try {
      status = process.start ();
   } catch (SxException e)  {
      e.print ();
   }
   printf ("process exit status: %d\n\n", status);
}

void testLsStdout ()
{
   cout << "TESTING: 'ls/ipconfig' command with readFromStdout ()" << endl;

   class MyProcess : public SxProcess
   {
      public:
         inline MyProcess(const SxString &cmdLine) : SxProcess(cmdLine) { }
         virtual void readFromStdout (const SxString &line)
         {
            cout << "parent received line from child's STDOUT: " << line << endl;
         }
   };

#ifdef WIN32
   MyProcess process("ipconfig");
#else
   MyProcess process("ls -l -s");
#endif
   process.setKillTree (false);
   process.enableBuffer (SxProcess::Stdout);
   try                    { process.start (); }
   catch (SxException e)  { e.print (); }
   cout << endl;
}

void testBuiltinStderr ()
{
   cout << "TESTING: built-in with readFromStderr ()" << endl;

   class MyProcess : public SxProcess
   {
      public:
         inline MyProcess(const SxString &cmdLine) : SxProcess(cmdLine) { }
         virtual void readFromStdout (const SxString &line)
         {
            // must not happen in this test; the child writes only to STDERR
            cout << "BUG: parent received line from child's STDOUT: " << line << endl;
         }
         virtual void readFromStderr (const SxString &line)
         {
            cout << "parent received line from child's STDERR: " << line << endl;
         }
   };

   MyProcess process(progName + " --builtin-child write-stderr");
   process.setKillTree (false);
   process.enableBuffer (SxProcess::Stderr);
   try                    { process.start (); }
   catch (SxException e)  { e.print (); }
   cout << endl;
}

void testBuiltinSynched ()
{
   cout << "TESTING: built-in with readFromStdout (), synchronizing STDOUT and "
           "STDERR" << endl;

   class MyProcess : public SxProcess
   {
      public:
         inline MyProcess(const SxString &cmdLine) : SxProcess(cmdLine) { }
         virtual void readFromStdout (const SxString &line)
         {
            cout << "parent received line from child's STDOUT: " << line << endl;
         }
         virtual void readFromStderr (const SxString &line)
         {
            // must not happen in this test because the child's STDERR output
            // is "synchronized" into the STDOUT stream
            cout << "BUG: parent received line from child's STDERR: " << line << endl;
         }
   };

   MyProcess process(progName + " --builtin-child write-both");
   process.setKillTree (false);
   process.enableBuffer (SxProcess::StdOutErr);
   process.synchronize (SxProcess::Stdout | SxProcess::Stderr);
   try                    { process.start (); }
   catch (SxException e)  { e.print (); }
   cout << endl;
}

enum { bufferSize = 10240 };

typedef struct
{
   SxString label;
   size_t used; // number of bytes
   char data[bufferSize];
} Buffer;

/** Writes a text line if available in the buffer. Returns whether the caller
    should call again because more lines may already be available in the
    buffer. */
bool showLine (Buffer &buf)
{
   SxString str = SxString (buf.data, buf.used);
   ssize_t pos = str.find ("\n");
   if (pos < 0)  {
      if (buf.used == bufferSize)  SX_THROW ("Line too long");
      return false; // must read some more first in order to get a '\n'
   }
   SxString line = str.head (pos + 1);
   cout << "parent received line from child's " << buf.label << ": " << line;
   str = str.tail (str.getSize () - pos - 1);
   ssize_t nBytes = str.getNBytes ();
   ::memcpy (buf.data, str.getElems (), nBytes);
   buf.used = nBytes;
   return true;
}

/** Reads from the file descriptor into the buffer. Returns whether the caller
   should continue select()ing/calling with this file descriptor. */
bool readBuffer (Buffer &buf, int fd)
{
   ssize_t n = ::read (fd, buf.data + buf.used, bufferSize - buf.used);
   if (n < 0)  {
      int err = errno;
      if (err == EINTR)  return true; // try again
      SX_THROW ("Can't read into buffer: " + sxstrerror (err));
   }
   if (n == 0)  return false; // EOF
   buf.used += n;
   while (showLine (buf))  { /* loop */ }
   return true;
}

#ifndef WIN32

/** Reads text lines from both STDOUT and STDERR */
void readStdOutErrLineWise (SxProcess &process)
{
   bool tryStdout = true, tryStderr = true;
   int fdOut = process.getPipeFd (STDOUT_FILENO),
       fdErr = process.getPipeFd (STDERR_FILENO),
       fdLimit = ::max (fdOut, fdErr) + 1;
   Buffer bufOut, bufErr;
   bufOut.used = bufErr.used = 0;
   bufOut.label = "STDOUT";
   bufErr.label = "STDERR";
   while (tryStdout || tryStderr)  {
      // --- prepare file descriptor sets
      fd_set fdset;
      FD_ZERO (&fdset);
      if (tryStdout)  FD_SET (fdOut, &fdset);
      if (tryStderr)  FD_SET (fdErr, &fdset);
      fd_set readfds = fdset, errorfds = fdset;

      // --- possibly sleep on the file descriptors until content becomes
      //     available
      int status = ::select (fdLimit, &readfds, NULL, &errorfds, NULL);
      if (status < 0)  {
         int err = errno;
         if (err == EINTR)  continue; // try again
         SX_THROW ("Can't read; " + sxstrerror (err));
      }

      // --- handle weird errors which shouldn't happen
      if ( (tryStdout) && (FD_ISSET (fdOut, &errorfds)) )  {
         cout << "Unknown exceptional condition for STDOUT" << endl;
         tryStdout = false;
      }
      if ( (tryStderr) && (FD_ISSET (fdErr, &errorfds)) )  {
         cout << "Unknown exceptional condition for STDERR" << endl;
         tryStderr = false;
      }

      // --- handle the content
      if ( (tryStdout) && (FD_ISSET (fdOut, &readfds)) &&
         (!readBuffer (bufOut, fdOut)) )  {
         tryStdout = false;
      }
      if ( (tryStderr) && (FD_ISSET (fdErr, &readfds)) &&
         (!readBuffer (bufErr, fdErr)) )  {
         tryStderr = false;
      }
   }
}

void testBuiltinStdOutErrSequential ()
{
   cout << "TESTING: built-in, reading from STDOUT and STDERR, sequential"
        << endl;

   class MyProcess : public SxProcess
   {
      public:
         inline MyProcess(const SxString &cmdLine) : SxProcess(cmdLine) { }
         virtual void readFromStdout (const SxString &line)
         {
            cout << "parent received line from child's STDOUT: " << line << endl;
         }
         virtual void readFromStderr (const SxString &line)
         {
            cout << "parent received line from child's STDERR: " << line << endl;
         }
   };

   // Here we must read "manually" because SxProcess::read () can't know from
   // which stream we want to read.
   MyProcess process(progName + " --builtin-child write-both");
   process.setKillTree (false);
   process.enableBuffer (SxProcess::StdOutErr);
   try                    { process.run ();
                            readStdOutErrLineWise (process);
                            process.wait (); }
   catch (SxException e)  { e.print (); }
   cout << endl;
}

void testBuiltinStdOutErrMixed ()
{
   cout << "TESTING: built-in, reading from STDOUT and STDERR, likely mixed"
        << endl;

   class MyProcess : public SxProcess
   {
      public:
         inline MyProcess(const SxString &cmdLine) : SxProcess(cmdLine) { }
         virtual void readFromStdout (const SxString &line)
         {
            cout << "parent received line from child's STDOUT: " << line << endl;
         }
         virtual void readFromStderr (const SxString &line)
         {
            cout << "parent received line from child's STDERR: " << line << endl;
         }
   };

   // Here we must read "manually" because SxProcess::read () can't know from
   // which stream we want to read.
   MyProcess process(progName + " --builtin-child write-mixed");
   process.setKillTree (false);
   process.enableBuffer (SxProcess::StdOutErr);
   try                    { process.run ();
                            readStdOutErrLineWise (process);
                            process.wait (); }
   catch (SxException e)  { e.print (); }
   cout << endl;
}

#endif // #ifndef WIN32

void testLongerStdoutPipe ()
{
   cout << "TESTING: 'ls | grep | cat' command with readFromStdout ()" << endl;

   class MyProcess : public SxProcess
   {
      public:
         inline MyProcess(const SxString &cmdLine) : SxProcess(cmdLine) { }
         virtual void readFromStdout (const SxString &line)
         {
            cout << "parent received line from child's STDOUT: " << line << endl;
         }
         virtual void readFromStderr (const SxString &line)
         {
            // should not happen in this test; the child processes should write
            // only to STDOUT
            cout << "PROBLEM: parent received line from child's STDERR: " << line << endl;
         }
   };

   // MyProcess process("ls -l -s | grep e | cat -"); would be nice
   MyProcess process(progName + " --builtin-child longer-pipe");
   process.setKillTree (false);
   process.enableBuffer (SxProcess::Stdout);
   try                    { process.start (); }
   catch (SxException e)  { e.print (); }
   cout << endl;
}

/** Lets the parent process write some text to the child's STDIN stream */
void testStdin ()
{
   cout << "TESTING: built-in ping pong with writeToStdin ()" << endl;

   class MyProcess : public SxProcess
   {
      public:
         inline MyProcess(const SxString &cmdLine) : SxProcess(cmdLine) { }
         virtual SxString writeToStdin ()
         {
            SxString res;
            res  = "First line passed from parent to child\n";
            res += "Second line passed from parent to child\n";
            res += "Third line passed from parent to child\n";
            return res;
         }
         virtual void readFromStdout (const SxString &line)
         {
            cout << "parent received line from child's STDOUT: " << line << endl;
         }
   };

   MyProcess process(progName + " --builtin-child pingpong");
   process.setKillTree (false);
   process.enableBuffer (
      (SxProcess::Channels) (SxProcess::Stdin | SxProcess::Stdout));
   try                    { process.run ();
                            process.closeStdStream (STDIN_FILENO); // EOF
                            process.read ();
                            process.wait (); }
   catch (SxException e)  { e.print (); }
}

/** The child process writes some text to its STDOUT stream */
void writeStdout ()
{
   cout << "First line written by child program to stdout" << endl;
   cout << "Second line written by child program to stdout" << endl;
   cout << "Third line written by child program to stdout" << endl;
}

/** The child process writes some text to its STDERR stream */
void writeStderr ()
{
   cerr << "First line written by child program to stderr" << endl;
   cerr << "Second line written by child program to stderr" << endl;
   cerr << "Third line written by child program to stderr" << endl;
}

/** The child process writes long texts to its STDOUT and STDERR streams */
void writeMixed ()
{
   // Here we want to write so many lines that they require more than one
   // ::read () operation per buffer. We write lines alternatingly ("mixed")
   // to STDOUT and STDERR.
   for (int i = 1; i < bufferSize / 30; i++)  {
      cout << "Line #" << i << " written by child program to stdout" << endl;
      cerr << "Line #" << i << " written by child program to stderr" << endl;
   }
}

/** The child process reads what the parent wrote to the child's STDIN stream
    and sends this text back via STDOUT to the parent ("ping pong"), with
    slight modifications */
void readAndSendBack ()
{
   while (!::feof (stdin))  {
      char buffer[1024];
      if (::fgets (buffer, sizeof (buffer), stdin) != buffer)  {
         int err = errno;
         if (::feof (stdin))  break;
         SX_THROW ("fgets () failed: " + sxstrerror (err));
      }
      size_t len = ::strlen (buffer);
      if (len > 500) // (to avoid buffer overflows)
         SX_THROW ("Line too long"); // should never happen with this program
      if ( (len > 0) && (buffer[len - 1] == '\n') )
         buffer[--len] = '\0'; // remove trailing newline
      ::sprintf (buffer + len, "; passed back by the child, old strlen %lu", len);
      if (::puts (buffer) == EOF)  {
         int err = errno;
         SX_THROW ("puts () failed: " + sxstrerror (err));
      }
   }
}

int main (int argc, char **argv)
{
#ifndef WIN32
   signal (SIGPIPE, SIG_IGN);
#endif
   if (argc == 0)  SX_THROW ("No command-line");
   argc--;
   progName = SxString(*argv++);
   if (argc == 0)  { // start the basic tests in the parent process
      testLsSimple ();
      testLsStdout ();
      testBuiltinStderr ();
      testBuiltinSynched ();
#ifndef WIN32 // (we don't have separate STDERR pipes on Windows)
      testBuiltinStdOutErrSequential ();
      testBuiltinStdOutErrMixed ();
      testLongerStdoutPipe ();
#endif
      testStdin ();
   }  else if ( (argc == 2) && (!::strcmp (argv[0], "--builtin-child")) )  {
      // executing as a child program; write some fantasy texts to streams
      if (!::strcmp (argv[1], "write-stdout"))  writeStdout ();
      else if (!::strcmp (argv[1], "write-stderr"))  writeStderr ();
      else if (!::strcmp (argv[1], "write-both"))  {
         writeStdout ();
         writeStderr ();
      }
      else if (!::strcmp (argv[1], "write-mixed"))  writeMixed ();
#ifndef WIN32 // (Windows usually doesn't know "grep" etc.)
      else if (!::strcmp (argv[1], "longer-pipe"))  {
        ::system ("ls -l -s | grep e | cat -");
      }
#endif
      else if (!::strcmp (argv[1], "pingpong"))  readAndSendBack ();
      else  SX_THROW ("unknown built-in process \"" + SxString(argv[1]) + "\"");
   }  else  SX_THROW ("Bad command-line arguments");

   return 0;
}
