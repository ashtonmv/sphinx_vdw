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

#ifndef _SX_PROCESS_H_
#define _SX_PROCESS_H_

#include <SxString.h>
#include <SxTimerThread.h>

#ifdef WIN32
   typedef DWORD pid_t;
#endif /* WIN32 */

/** \brief Start external programs and communicate with them.

    \b SxProcess = SPHInX Process Launcher

    This class allows to execute external programs by creating new subprocesses.
    SPHInX can communicate with a subprocess by reading from or writing to the
    streams stdin, stdout, and stderr of the subprocess. Additionally, a new
    process can be suspended or terminated.

    On POSIX operating systems, the SIGPIPE signal handler should be set to
    SIG_IGN (or to a "reasonable" function). By default, a problem with one of
    the inter-process pipes causes the operating system to terminate the
    program entirely, which looks like a crash. Cf. the manual pages signal(2)
    and signal(7).

    \par Example 1
    Execute the UNIX command "ls -l -s" without any further inter-process
    communication.
\code
   int status=0;
   SxProcess process ("ls -l -s");
   try {
      status = process.start ();
   } catch (SxException e)  {
      e.print ();
   }
   printf ("process exit status = %d\n", status);
\endcode

    \par Example 2
      The following code snippet demonstrates how to read from STDOUT of
      the subprocess. Here the class MyProcess has been derived from SxProcess
      in order to overload the SxProcess::readFromStdout () function. 
\code
   class MyProcess : public SxProcess
   {
      public:
         MyProcess (const SxString &cmdLine) : SxProcess (cmdLine) { }
         virtual void readFromStdout (const SxString &line)
         {
            printf ("process line %s from STDOUT\n", line.getElems ());
         }
   };

   signal (SIGPIPE, SIG_IGN);
   MyProcess process ("ls -l -s");
   process.setCommunication (SxProcess::StdOut);
   try                    { process.start (); }
   catch (SxException e)  { e.print (); }
\endcode

    \par Example 3
       If the external program writes data to STDOUT \b and STDERR both
       channels can be synchronized and combined to a new STDOUT.
\code
   class MyProcess : public SxProcess
   {
      public:
         MyProcess (SxString &cmdLine) : SxProcess (cmdLine) { }
         virtual void readFromStdout (const SxString &line)
         {
            printf ("process line %s from STDOUT/STDERR\n", line.getElems ());
         }
   };

   signal (SIGPIPE, SIG_IGN);
   MyProcess process ("ls -l -s");
   process.setCommunication (SxProcess::StdOut | SxProcess::Stderr);
   process.synchronize (SxProcess::StdOut | SxProcess::Stderr);
   try                    { process.start (); }
   catch (SxException e)  { e.print (); }
\endcode

    \ingroup group_os
    \author  Sixten Boeck, boeck@mpie.de */
class SX_EXPORT_IPC SxProcess
{
   public:

      /** \brief Communication channels

          This enum type defines the communication channels between
          the parent's and the child's process.  The communication
          channels can be changed with SxProcess::setCommunication (). */
      enum Channels {
         /** \brief No pipe communication

            Child process writes to parent's STDOUT and parent's STDERR and
            reads from parent's STDIN. */
         None = 0x00,
         /** \brief Child can read from STDIN pipe.

            The parent can write data to the child process's STDIN channel.
            To do so, the derived class has to overload the
            SxProcess::writeToStdin () function. */
         StdIn = 0x01,
         Stdin = 0x01, // for backward compatibility
         /** \brief Child can write to STDOUT pipe.

            The parent can read data from the child process's STDOUT channel.
            In order to do so, the derived class has to overload the
            SxProcess::readFromStdout () function. */
         StdOut = 0x02,
         Stdout = 0x02, // for backward compatibility
         /** \brief Child can write to STDERR pipe.

            The parent can read data from the child process's STDERR channel.
            In order to do so, the derived class has to overload the
            SxProcess::readFromStderr () function. */
         StdErr = 0x04,
         Stderr = 0x04, // for backward compatibility

         StdOutErr = StdOut | StdErr,

         StdInOut = StdIn | StdOut,

         StdInErr = StdIn | StdErr,

         StdInOutErr = StdOutErr | StdIn
      };

      /** \brief Channel synchronization model

          This enumeration type describes the synchronization relationships
          between output channels of the child process. The synchronization
          mode can be modified using the SxProcess::synchronize () function. */
      enum Sync {
         /** \brief  No synchronization of channels

             STDOUT and STDERR are not synchronized. Both are handled
             separately. This is the default. */
         NoSync,
         /** \brief Combine STDOUT and STDERR

             The STDOUT and STDERR channels of the child process are
             "synchronized": STDERR is the same as STDOUT. From the perspective
             of the parent process, both channels are combined into STDOUT. */
         SyncOutErr
      };

      enum Mode {
         /** \brief run process as usual child with session ID of parent */
         Child,
         /** \brief run child process with session ID of its own; if you use
             this mode, you must set 'killTree' to 'false' */
         Daemon
      };

      /** \brief Default constructor

          Constructs an SxProcess. Neither the command nor the argument list
          is initialized. Before calling SxProcess::start (), they have to be
          set up using SxProcess::setCommand () and SxProcess::setArguments ()
          resp. SxProcess::addArgument (). */
      SxProcess (Mode mode_ = Child);

      /** \brief Construct a SxProcess with a given command line

          This constructor initializes the command and arguments which will
          be executed when calling SxProcess::start */
      SxProcess (const SxString &cmdLine, Mode = Child);

      /** \brief Construct a SxProcess with a command and an argument list

          This constructor initializes the command as well as an argument
          list. The command will be started when calling SxProcess::start */
      SxProcess (const SxString &cmd, const SxList<SxString> &argList_,
                 Mode = Child);

      /** \brief Destructor */
      virtual ~SxProcess ();

      /** \brief get value of an environement variable

          This function returns the value of an environment variable.
          Usually environment variables are set with shell commands
          (bash: export, csh: setenv). If the variable does not exist,
          an empty SxString is returned.
          \sa SxProcess::setenv
       */
      static SxString getenv (const SxString &varName);
      /** \brief set environment variable

          This function is equivalent to the bash shell's export or to
          the (t)csh's setenv.
          \sa SxProcess::getenv */
      static void setenv (const SxString &varName, const SxString &value);

      /** \brief daemonize the current process

           \par Example:
\code
   SxProcess::daemon ("/", "/var/log/sxdaemon.log");
   cout << "Daemon started" << endl;
\endcode

           \par Example:
\code
   // --- deamon in the current working directory (should not be a mount point)
   SxProcess::daemon ();

   SxLogBuf logBuf;
   logBuf.setFile ("/var/log/sxdaemon.log");
   SxRedirect tee (std::cout, &logBuf);

   cout << "Daemon started" << endl;
\endcode
      \sa examples/sxdaemon.cpp
      */
      static void daemon (int umask_=022,
                          const SxString &chdir="", const SxString &logfile="");

      /** \brief Return process id

          Returns process identification number
      */
      pid_t getPid ();

      pid_t getChildProcPID ();
      pid_t getChildProcThread0ID ();

      /** \brief Return parent pid

          Returns parent or calling process identification

          \par Windows:
          It seems there is no process tree in windows.
          Parent pid is valid if parent is still alive otherwise it can be
          already another process which took the free pid number.
      */
      static pid_t getPPid (); // #exception

      /** \brief Set inter-process communication model.

          This function can be used to enable/disable the read-from
          or write-to channels. The default is SxProcess::None. That means
          that the process writes to the parent's stdout and stderr channels.

          \par Example:
          Let's assume the class SxMyProcess is derived from SxProcess and
          overloads SxProcess::readFromStdout
\code
SxMyProcess process;
process.setCommunication (SxProcess::StdOut);
\endcode
          \sa writeToStdin
          \sa readFromStdout
          \sa readFromStderr
       */
      void setCommunication (Channels);
      void enableBuffer (Channels);

      /** \brief Change the current synchronization model.

          The two channels STDOUT and STDERR and be treated separately
          or a single combined and synchronized channel. In the latter case
          the combined channel is STDOUT again. Then STDERR is not used
          anylonger.
          The default is not to synchronize STDOUT and STDERR.

          \par Example: Do not synchronize STDOUT and STDERR (default)
\code
   SxProcess process;
   process.synchronize (SxProcess::None);
\endcode
          \par Example: Synchronize STDOUT and STDERR
\code
   SxProcess process;
   process.synchronize (SxProcess::StdOut | SxProcess::Stderr);
\endcode

      \par mask Bitwise OR'ed channels which should be synchronize.
                Note, that only STDOUT and STDERR can be synchronized. */
      void synchronize (int mask);

      /** \brief Sets or overwrites the command.

          The given command is executed when SxProcess::start is invoked. */
      void setCommand (const SxString &);
      /** \brief Sets the entire command line

          The given string is tokenized by whitespaces in order to separate
          the command and the arguments. */
      void setCommandLine (const SxString &);

      /** \brief Clean up internal argument list

          This function cleans the argument list of a command. It is useful in
          combination with SxProcess::addArgument. */
      void removeAllArguments ();

      /** \brief Sets or overwrites the argument list.

          The given argument list is taken when SxProcess::start is invoked. */
      void setArguments (const SxList<SxString> &);

      /** \brief Add an argument to the list of arguments.

          Append argument to the end of the process's argument list. The
          argument list is used when invoking SxProcess::start. */
      void addArgument (const SxString &);

      /** \brief Returns the current list of arguments.

          This function returns a copy of the current argument list.  */
      SxList<SxString> getArguments () const;

      /** \brief Set process's working directory.

          On default the process is executed in the current working
          directory (of the parent's process). By calling this function
          a different working directory can be chosen for the child process.

          \b Note, the provided working directory \b must be an absolute
          path!  */
      void setWorkingDirectory (const SxString &);

      /** \brief Setting the input buffer size.

          By default, the size of the Stdin, Stdout, and Stderr communication
          buffers is set to 10k. However, if the child process writes data
          without explicit flushing of the buffer no intermediate data are
          sent to one of the readFromXXX functions. To overcome this problem
          the buffer size can be reduced. So flushing can be enforced. */
      void setBufferSize (ssize_t);

      /** \brief Kills a process if wallclock runtime exceeds the timeout.

\code
   // --- kill the process after 20 seconds
   SxProcess proc ("cat -");
   proc.setTimeout (20.);
   proc.start ();
\endcode
      */
      void setTimeout (double);

      /** \brief Kill the process tree internally.

          Enabled by default.

          Internal kill calls in this SxProcess will call kill() on
          the process together with its subtree. */
      void setKillTree (bool enabled);

      /** \brief Execute a command.

          This function creates a child process, sets it up and lets it start
          the external program. The command to be started must have been
          specified e.g. with SxProcess::setCommand () and the
          argument list with SxProcess::setArguments () resp.
          SxProcess::addArgument ().
          The system environment variables are equal to those of the calling
          process. However, the environment of the child process can be
          provided by the optional string list argument. Every list item
          must have a syntax like "key=value".
          The process state can be inquired with ::isRunning. It can be
          terminated with ::softTerminate or ::kill.
          \throw SxException::ForkFailed  process did not have enough
                                          memory to copy process stack
          \throw SxException::ExecFailed  process could be be executed.
          \param env                      list of environment variables used
                                          for the child process
          \setenv PATH
          \setenv LD_LIBRARY_PATH
          \return exit status of the child process

          \verbatim
            try {
               SxProcess proc("/bin/ls");
               proc.enableBuffer (SxProcess::StdOut);
               int status = proc.start ();
               cout << proc.getBuffer (SxProcess::StdOut) << endl;
            } catch (SxException e) {
               e.print ();
            }
          \endverbatim

          \verbatim
            try {
               SxProcess proc("/bin/ls");
               proc.enableBuffer (SxProcess::StdOut);
               proc.run ();
               proc.read ();
               int status = proc.wait ();
               cout << proc.getBuffer (SxProcess::StdOut) << endl;
            } catch (SxException e) {
               e.print ();
            }
          \endverbatim
       */
      virtual int start (const SxList<SxString> &env=SxList<SxString>()); // #exception

      virtual void run (const SxList<SxString> &env=SxList<SxString>()); // #exception

      /** Reads from the STDOUT resp. STDERR stream of the child process and
          provides the resulting texts to readFromStdout () resp.
          readFromStdErr () as appropriate. If you configured separate channels
          for STDOUT and STDERR, this method can't know which of the streams is
          meant, so you must e.g. overload this method with a custom method
          which reads from both streams, e.g. as indicated by a select () or
          poll () system call (possibly with a zero timeout parameter). */
      virtual void read (); // #exception
      virtual void readBuffer (int FD);

      /** Writes the given string to the STDIN stream of the child process */
      virtual void write (const SxString &); // #exception
      /** Writes the given number of bytes to the STDIN stream of the child
          process */
      virtual void write (const char *data, ssize_t origNBytes); // #exception
      virtual int wait (); // #exception

      /** \brief Inform derived class that execution of child process has
                 finished.

          This function is called when the child process has finished. That 
          can be due to normal exit of the external program or due to
          explicit termination (::softTerminate, ::kill). If you need to be
          informed about program termination overload this function in the
          derived class. */
      virtual void finishExecution ();

      /**  \brief Return exit status of the process.

           This function returns the exit status from the child process
           after it was launched with SxProcess::start. Note, that the
           exit status is also returned from SxProcess::start directly.

           \par Example:
\code
   SxProcess process ("ls -l");
   try  { process.start (); }
   catch (SxException e)  { e.print (); }

   printf ("Return status of 'ls': %d\n", process.exitStatus ());
\endcode
      */
      int exitStatus () const;

      /**  \brief Return the signal number which caused the process to exit.

           This function returns the signal number from the exit status
           of the child process after it was launched with SxProcess::start.

           Succesful process execution returns value 0.

           \par Windows:
           Timeouts and errors from waiting until the process
           terminates will return some value other than 0. SIGKILL and SIGTERM
           values are returned if the process was terminated using SxProcess
           interface. Ending the process from the outside, for example from
           taskmgr, will return exit signal 0. In such case at least the exit
           status will be non zero. */
      int exitSignal () const;

      SxString getBuffer (Channels);


      /** \brief Returns true if the child process is running currently.

          As long as the child process is running this function returns
          true.  */
      bool isRunning () const; 

      /** \brief Suspend the whole Process Tree under Windows */
      void suspend ();

      /** \brief Resumes the whole process tree under Windows */
      void resume ();

      /** \brief Suspend a single thread under Windows */
      void suspendThread (int threadId);

      /** \brief Resumes a single thread under Windows */
      void resumeThread (int threadId);


      /** \brief Try to terminate a running child process gently.

           This function tries to terminate a child process softly by
           sending a SIGTERM signal. The external program (child process) can
           respond by trapping the SIGTERM signal. So, an emergency function
           of the external program can be invoked, which, for example, writes
           data before exiting.
           An immediate termination can be done with the SxProcess::kill. */
      void softTerminate (int signalId = SIGTERM); // #exception

      /** \brief Kill child process immediately.

           This function terminates a child process immediately by
           sending a SIGKILL signal. The external program (child process) cannot
           respond by trapping the SIGKILL signal. So, no emergency function
           of the external program can be invoked.
           A soft termination can be invoked by SxProcess::softTerminate. */
      void kill (); // #exception

      /** \brief Send the specified signal to the specified process.
          \return exit status if waitForTermination is enabled */
      static int kill     (pid_t pid,                 // #exception
                           int   signalId = SIGTERM,
                           bool  tree = false,
                           bool  waitForTermination = false);

      /** \brief Writing to the child's STDIN channel

           When an inter-process communication using the child's STDIN channel
           has been enabled (::setCommunication(::Stdin)) the ::writeToStdin
           function has to be overloaded. The returned string then is passed to
           the child's STDIN channel.

           \par Example:
\code
class MyProcess : public SxProcess
{
   public:
      MyProcess ();
      virtual ~MyProcess ();

      virtual SxString writeToStdin ()
      {
         SxString res;
         res  = "This is the first line to be passed to the CHILD's STDIN\n";
         res += "This is the second line to be passed to the CHILD's STDIN\n";
         res += "This is the third line to be passed to the CHILD's STDIN\n";
         return res;
      }
};
\endcode
      \sa      readFromStdout
      \sa      readFromStderr
      \return  String which is passed to the STDIN of the child process.
       */
      virtual SxString writeToStdin ();

      /** \brief Enable/Disable lineWise reading of output channels.
          Enabled by default. */
      void setLineWise (bool enable);

      /** \brief Read from child's STDOUT channel.

           When an inter-process communication using the child's STDOUT channel
           has been enabled (::setCommunication(::StdOut)) or the STDERR
           channel of the child process has been combined and synchronized
           with STDOUT (::synchronize(SxProcess::SyncOutErr))
           the ::readFromStdout function has to be overloaded.
           The provided argument contains a line which was received from
           STDOUT or from the combined STDOUT/STDERR.

           \par Example:
\code
class MyProcess : public SxProcess
{
   public:
      MyProcess ();
      virtual ~MyProcess ();

      virtual void readFromStdout (const SxString &line) const
      {
         cout << "Received line: " << line << " from external program.\n";
      }
};
\endcode
      \sa writeToStdin
      \sa readFromStderr
      \sa synchronize
       */
      virtual void readFromStdout (const SxString &);

      /** \brief Read from child's STDERR channel.

           When an inter-process communication using the child's STDERR channel
           has been enabled (::setCommunication(::Stderr)) the
           ::readFromStderr function has to be overloaded.
           The provided argument contains a line which was received from
           STDERR.

           Note, if the STDERR has been combined and synchronized with the
           STDOUT channel using ::synchronize the output is redirected to
           ::readFromStdout!

           \par Example:
\code
class MyProcess : public SxProcess
{
   public:
      MyProcess ();
      virtual ~MyProcess ();

      virtual void readFromStderr (const SxString &line) const
      {
         cout << "Received error line: " << line
              << " from external program.\n";
      }
};
\endcode
      \sa writeToStdin
      \sa readFromStdout
      \sa synchronize
       */
      virtual void readFromStderr (const SxString &line);

      /** Closes the parent process end of a standard stream pipe. The parameter
          must be 0 or 1 on Windows resp. 0..2 on POSIX (STDIN_FILENO,
          STDOUT_FILENO or STDERR_FILENO, referring to the desired standard
          stream) */
      void closeStdStream (int stdStreamFd);

#     ifdef WIN32
         /** \brief If enabled, will start the process with Windows
               SW_HIDE hidden window flag*/
         void setHiddenWindow (bool value);
#     endif

#ifndef WIN32
      /** Returns a pipe file descriptor for the parent process. The parameter
          must be in the range 0..2 (STDIN_FILENO, STDOUT_FILENO or
          STDERR_FILENO, referring to the desired standard stream). */
      int getPipeFd (int stdStreamFd) const;
#endif

   protected:

      /** \brief File descriptor numbers of the default channels */
      enum ChannelNumbers { StdinId = 0, StdoutId = 1, StderrId = 2 };

      /** \brief Command which is executed by SxProcess::start */
      SxString          cmd;
      /** \brief Command's arguments used by SxProcess::start */
      SxList<SxString>  argList;

      /** \brief current communication model */
      Channels channels;

      /** \brief buffers for stdout and stderr */
      Channels useBuffers;
      SxString outBuffer, errBuffer;

      /** \brief run process as daemon, i.e., create new session id */
      Mode mode;

      /** \brief current synchronization model */
      Sync sync;

      /** \brief Working directory for the child process. */
      SxString          workDir;

      /** \brief The id of the current (or last invoked) process. */
      pid_t pid;

      /** \brief Use tree flag in internal kill calls. */
      bool killTree;

      /** \brief Indicates whether the child process is running or not.*/
      mutable SxMutex statusMutex;
      bool  processRunning;

      /// Whether to read output linewise
      bool lineWise;

      /// Do not suspend execution of the calling thread for wait status
      bool noHang;

      /** \brief Exit status of previous run of external program. */
      int processExitStatus;

      /** \brief Signal number which caused the previous run to exit. */
      int processExitSignal;

      /** \brief Length of the input buffer */
      ssize_t bufLen;

#     ifdef WIN32
         /** \brief Flag whenever should the Windows process be started with
                    SW_HIDE hidden window flag*/
         bool hiddenWindow;
#     endif

      /** \brief How many seconds to wait before the child process exits.
          The default value (x < 0) will wait forever. */
      double processTimeout;

      /** \brief check to see if child process has closed it's StdErr pipe*/
      bool eofStdErr;

      /** \brief check to see if child process has closed it's StdOut pipe*/
      bool eofStdOut;

      void readLineWise (bool useStderr = false);
      void readBufferLineWise (int FD, SxPtr<SxArray<char>> *bufPtr,
                               SxPtr<ssize_t> *curIdxPtr);

      void initPid ();
      /* Marks all pipe-related file descriptors/handles as invalid */
      void initPipes ();

      /** The ends of parent-child standard stream pipes */
      enum PipeEnds {
         childEnd = 0x01, parentEnd = 0x02, bothEnds = childEnd | parentEnd
      };

      /** This low-level function closes one end of a parent-child standard
          stream pipe. */
      void closePipeIdx (int stdStreamFd, uint8_t idx);

      /** Closes one or both ends of a parent-child standard stream pipe. This
          function may only be called in the parent process. */
      void closePipe (int stdStreamFd, PipeEnds ends);

      /** Ensures that both ends of all existing parent-child pipes are closed
      */
      void closePipes ();

      /** Returns whether the channels STDOUT and STDERR are (resp. shall be)
          combined into STDOUT */
      inline bool combinedStdout () const
      { return (sync == SyncOutErr) && ((channels & StdOutErr) == StdOutErr); }

#ifdef WIN32
      /** Windows-specific data type for pipes. Handle array index 0 contains
          the reading end of the pipe; array index 1 contains the writing end.
      */

      typedef struct { HANDLE handles[2]; OVERLAPPED oOverlap; } WindowsPipe;
      WindowsPipe pipes[3]; // pipe #0 for STDIN, pipe #1 for STDOUT

      HANDLE hEvents[3];
      class SX_EXPORT_IPC SxProcExecuter : public SxSystemThread
      {
         public:
            SxProcess *obj;
            SxString cmd;
            HANDLE mainProc;

            int processExitStatus; // After thread terminates
            int processExitSignal; // After thread terminates

            // Thread started the process and obtained pid
            pid_t pid;
            pid_t thread0ID;
            bool gotPid;
#           ifdef WIN32
               bool hiddenWindow;
#           endif
            /** Handle for the Windows JobObject which can capture all
                (grand-)child processes */
            HANDLE jobObjectHndl;
            HANDLE jobSigPort;
            SxMutex mutex;
            SxThreadCond condition;
            pid_t waitForPid ();

            SxProcExecuter (SxProcess *, const SxString &);
           ~SxProcExecuter ();
            virtual void main ();
            /** Creates a process object */
            virtual void createProcess (const SxString &cmd);
            /** Creates a JobObject which can capture all (grand-)child
                processes */
            void createJobObject ();
            /** Removes the JobObject and terminates all (grand-)child
                processes */
            void removeJobObject ();
            void removeProcess ();
            void waitForProcess (HANDLE hProcess);
            void waitForJob ();
            ssize_t getNProcesses ();
            SxArray<HANDLE> getProcessHandles (DWORD access);
            void closeProcessHandles (SxArray<HANDLE> &procs);
      };
      SxPtr<SxProcExecuter> childProc;
      static void winkill (pid_t pid, UINT exitCode, bool wait);
      unsigned char getRndByte () const;
      SxString getUUIDv4 ();
#else /* WIN32 */
      /** POSIX-specific data type for pipes. File descriptor array index 0
          contains the reading end of the pipe; array index 1 contains the
          writing end. */
      typedef struct { int fds[2]; } PosixPipe;
      PosixPipe pipes[3]; // for STDIN, STDOUT and STDERR
      class SxProcessTimer : public SxTimerThread
      {
         public:
            SxProcessTimer (pid_t pid_, bool killTree);
            virtual void interval ();
         protected:
            pid_t pid;
            int n;
            bool killTree;
      };
      SxPtr<SxProcessTimer> processTimer;
#endif

};
//SX_REGISTER_CLASS (SxProcess);
#endif /* _SX_PROCESS_H_ */
