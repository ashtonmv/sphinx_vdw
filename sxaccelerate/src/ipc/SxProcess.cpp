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
#include <SxConfig.h> 
#include <SxTime.h>
#include <stdlib.h> 
#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h> 
#endif
#ifdef HAVE_SYS_WAIT_H
#  include <sys/wait.h> 
#endif
#ifdef HAVE_SYS_RESOURCE_H
#  include <sys/resource.h> 
#endif
#ifdef HAVE_UNISTD_H
#  include <unistd.h> 
#endif
#ifdef HAVE_SIGNAL_H
#  include <signal.h> 
#endif
#include <sys/stat.h>
#ifdef HAVE_SYS_RESOURCE_H
#  include <sys/resource.h>
#endif
#include <fcntl.h>
#include <errno.h>

#ifdef WIN32
#  include <windows.h>
//#  include <process.h>
#  include <sxputenv.h>
#  include <direct.h>
#  include <tlhelp32.h> // CreateToolhelp32Snapshot
#  include <aclapi.h>
#  include <VersionHelpers.h>
#  include <cstdlib>
#  include <ctime>
#else
#  include <dirent.h>
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

void SxProcess::suspendThread (int threadID)
{
   SX_TRACE ();
#  ifdef WIN32
   HANDLE hThread = OpenThread (THREAD_ALL_ACCESS, FALSE, threadID);
   SuspendThread (hThread);
   CloseHandle (hThread);
#  else
   //TODO: Implement me
   (void)threadID; // avoid warning
   SX_EXIT;
#  endif
}

void SxProcess::resumeThread (int threadID)
{
   SX_TRACE ();
#  ifdef WIN32
   HANDLE hThread = OpenThread (THREAD_ALL_ACCESS, FALSE, threadID);
   ResumeThread (hThread);
   CloseHandle (hThread);
#  else
   //TODO: Implement me
   (void)threadID; // avoid warning
   SX_EXIT;
#  endif
}

void SxProcess::suspend () 
{
   SX_TRACE ();

#  ifdef WIN32
      typedef LONG (NTAPI *NtSuspendProcess)(IN HANDLE ProcessHandle);
      NtSuspendProcess pfnNtSuspendProcess 
         = (NtSuspendProcess)GetProcAddress(GetModuleHandle("ntdll"), "NtSuspendProcess");

      if (  childProc.getPtr ()
         && IsWindows8OrGreater ())  {
         SxArray<HANDLE> processHandles = childProc->getProcessHandles (PROCESS_ALL_ACCESS);
         ssize_t nProcs = processHandles.getSize ();
         for (ssize_t iProc = 0; iProc < nProcs; ++iProc)  {
            pfnNtSuspendProcess(processHandles(iProc));
            CloseHandle(processHandles(iProc));
         }
      } else  {
         HANDLE processHandle = OpenProcess(PROCESS_ALL_ACCESS, FALSE, pid);
         pfnNtSuspendProcess(processHandle);
         CloseHandle(processHandle);
      }
#  else
      SX_EXIT; // to be implemented
#  endif
}

void SxProcess::resume () 
{
   SX_TRACE ();

#  ifdef WIN32
      typedef LONG (NTAPI *NtResumeProcess)(IN HANDLE ProcessHandle);
      NtResumeProcess pfnNtResumeProcess 
         = (NtResumeProcess)GetProcAddress(GetModuleHandle("ntdll"), "NtResumeProcess");

      if (  childProc.getPtr ()
         && IsWindows8OrGreater ())  {
         SxArray<HANDLE> processHandles = childProc->getProcessHandles (PROCESS_ALL_ACCESS);
         ssize_t nProcs = processHandles.getSize ();
         for (ssize_t iProc = 0; iProc < nProcs; ++iProc)  {
            pfnNtResumeProcess(processHandles(iProc));
            CloseHandle(processHandles(iProc));
         }
      } else  {
         HANDLE processHandle = OpenProcess(PROCESS_ALL_ACCESS, FALSE, pid);
         pfnNtResumeProcess(processHandle);
         CloseHandle(processHandle);
      }
#  else
      SX_EXIT; // to be implemented
#  endif
}


bool SxProcess::isRunning () const
{
   SX_TRACE ();
   bool result;
   SX_MUTEX (statusMutex)  {
      result = processRunning;
   }

   return result;
}

#ifdef WIN32
SxProcess::SxProcExecuter::SxProcExecuter (SxProcess *obj_, const SxString &cmd_)
   : SxSystemThread (),
     obj(obj_), cmd(cmd_),
     processExitStatus(0),
     processExitSignal(0),
     pid(0),
     gotPid(false),
     hiddenWindow(false),
     jobObjectHndl(NULL)
{
   // empty
}

SxProcess::SxProcExecuter::~SxProcExecuter ()
{
   if ( (obj->killTree) && (jobObjectHndl != NULL) )  removeJobObject ();
   if (mainProc)  CloseHandle (mainProc);
}

void SxProcess::SxProcExecuter::createProcess (const SxString &cmd)
{
   SX_TRACE ();
   SX_CHECK(jobObjectHndl);

   // --- Initialize process info
   PROCESS_INFORMATION procInfo;
   ::ZeroMemory(&procInfo, sizeof (PROCESS_INFORMATION));

   // --- Initialize startup information 
   STARTUPINFO si;
   ::ZeroMemory(&si, sizeof (STARTUPINFO));
   si.cb = sizeof (STARTUPINFO);

   // ---Adjust startup information for the usage of pipes
   bool usePipes = (obj->channels & (StdIn | StdOut | StdErr)) != 0;
   if (obj->sync != NoSync)  {
      si.hStdError = obj->pipes[STDOUT_FILENO].handles[1];
      si.hStdOutput = obj->pipes[STDOUT_FILENO].handles[1];
   } else  {
      if (obj->channels & StdErr)
         si.hStdError = obj->pipes[STDERR_FILENO].handles[1];
      else
         si.hStdError = GetStdHandle(STD_ERROR_HANDLE);
      if (obj->channels & StdOut)
         si.hStdOutput = obj->pipes[STDOUT_FILENO].handles[1];
      else
         si.hStdOutput = GetStdHandle(STD_OUTPUT_HANDLE);
   }
   if (obj->channels & StdIn)
      si.hStdInput = obj->pipes[STDIN_FILENO].handles[0];
   else
      si.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
   if (usePipes)  si.dwFlags |= STARTF_USESTDHANDLES;
   // --- hide window
   //     cmd window should not show up when starting processes
   if (hiddenWindow) {
      si.wShowWindow = SW_HIDE;
      si.dwFlags |= STARTF_USESHOWWINDOW;
   }

   // --- prepare Security Identification
   //PSID sidPtr = NULL;
   //SID_IDENTIFIER_AUTHORITY sidWorld = SECURITY_WORLD_SID_AUTHORITY;
   //AllocateAndInitializeSid(&sidWorld, 1, SECURITY_WORLD_RID, 0, 0, 0, 0, 0, 0, 0, &sidPtr);

   // --- setup an access control list entry
   //EXPLICIT_ACCESS ea;
   //ZeroMemory(&ea, sizeof(EXPLICIT_ACCESS));
   //ea.grfAccessPermissions = STANDARD_RIGHTS_ALL;
   //ea.grfAccessMode = GRANT_ACCESS;
   //ea.grfInheritance = CONTAINER_INHERIT_ACE;
   //ea.Trustee.TrusteeForm = TRUSTEE_IS_SID;
   //ea.Trustee.TrusteeType = TRUSTEE_IS_WELL_KNOWN_GROUP;
   //ea.Trustee.ptstrName = (LPSTR)sidPtr;

   // --- setup an access control list (pointer)
   //ACL *aclPtr = NULL;
   //DWORD err = SetEntriesInAcl(1, &ea, NULL, &aclPtr);
   //if (err != ERROR_SUCCESS) {
   //   SX_THROW("Can't setup access control: " + sxstrerror(err));
   //}

   // --- prepare security descriptor
   //SECURITY_DESCRIPTOR *sdPtr = NULL;
   //sdPtr = (SECURITY_DESCRIPTOR*)LocalAlloc(LPTR, SECURITY_DESCRIPTOR_MIN_LENGTH);
   //InitializeSecurityDescriptor(sdPtr, SECURITY_DESCRIPTOR_REVISION);
   //SetSecurityDescriptorDacl(sdPtr, TRUE, aclPtr, FALSE);

   // --- prepare security attributes
   //SECURITY_ATTRIBUTES sa;
   //sa.nLength = sizeof(SECURITY_ATTRIBUTES);
   //sa.lpSecurityDescriptor = sdPtr;
   //sa.bInheritHandle = TRUE;

   // -- create suspended process
   int success;
   if (cmd.isUnicode ()) {
      LPSTARTUPINFOW siWPtr = (STARTUPINFOW *)&si;
      success = CreateProcessW (NULL, (LPWSTR)cmd.utf16 ().elements,        
         NULL, NULL, FALSE, CREATE_SUSPENDED,
         NULL, NULL, siWPtr, &procInfo);
      SX_DBG_MSG ("Process " << cmd << " created: " 
                  << (success ? "successful" : "failed"));
   } else  {
      success = CreateProcessA (NULL, (LPSTR)cmd.ascii (),
         NULL, NULL, FALSE, CREATE_SUSPENDED,
         NULL, NULL, &si, &procInfo);
      SX_DBG_MSG ("Process " << cmd << " created: "
                  << (success ? "successful" : "failed"));
   }

   if (success == 0) {
      DWORD err = ::GetLastError ();
      SX_THROW("Can't create process: " + sxstrerror(err));
   }

   // assign processes in job object for windows 8 and higher. Windows 7 does not allow for
   // multiple job objects.
   if (jobObjectHndl != NULL && !obj->noHang && IsWindows8OrGreater ())  {
      success = ::AssignProcessToJobObject(jobObjectHndl, procInfo.hProcess);
      mainProc = procInfo.hProcess;

      if (success == 0) {
         DWORD err = ::GetLastError ();
         CloseHandle (procInfo.hProcess);
         CloseHandle (procInfo.hThread);
         SX_THROW("Can't add child process to job object: " + sxstrerror(err) +
            ((err == ERROR_ACCESS_DENIED) ?
               "\nOn Windows Vista, add a proper UAC manifest to the application!"
               : ""));
      }
   }
   // Close pipe handles that are only needed by child process
   if (obj->sync == NoSync) {
      if (obj->channels & StdErr)
         CloseHandle(obj->pipes[STDERR_FILENO].handles[1]);
   }
   if (obj->channels & StdOut)
      CloseHandle(obj->pipes[STDOUT_FILENO].handles[1]);
   if (obj->channels & StdIn)
      CloseHandle(obj->pipes[STDIN_FILENO].handles[0]);
   // --- start process
   ResumeThread (procInfo.hThread);

   // --- signal that process is now running with following pid
   SX_MUTEX(mutex) {
      pid = procInfo.dwProcessId;
      thread0ID = procInfo.dwThreadId;
      gotPid = true;
      condition.wakeOne();
   }

   SX_MUTEX (obj->statusMutex)  {
      obj->processRunning = true;
   }

   // --- Clean up
   //if (aclPtr != NULL)  LocalFree(aclPtr);
   //if (sidPtr != NULL)  LocalFree(sidPtr);
   //if (sdPtr != NULL)   LocalFree(sdPtr);
   //CloseHandle (procInfo.hProcess);
   CloseHandle (procInfo.hThread);

   SX_DBG_MSG ("Process " << cmd << " running");
}

void SxProcess::SxProcExecuter::createJobObject ()
{
   SX_TRACE ();
   // --- prepare Security Identification
   //PSID sidPtr = NULL;
   //SID_IDENTIFIER_AUTHORITY sidWorld = SECURITY_WORLD_SID_AUTHORITY;
   //AllocateAndInitializeSid(&sidWorld, 1, SECURITY_WORLD_RID, 0, 0, 0, 0, 0, 0, 0, &sidPtr);
   
   // --- setup an access control list entry
   //EXPLICIT_ACCESS ea;
   //ZeroMemory (&ea, sizeof (EXPLICIT_ACCESS));
   //ea.grfAccessPermissions = STANDARD_RIGHTS_ALL;
   //ea.grfAccessMode = GRANT_ACCESS;
   //ea.grfInheritance = CONTAINER_INHERIT_ACE;
   //ea.Trustee.TrusteeForm = TRUSTEE_IS_SID;
   //ea.Trustee.TrusteeType = TRUSTEE_IS_WELL_KNOWN_GROUP;
   //ea.Trustee.ptstrName = (LPSTR) sidPtr;

   // --- setup an access control list (pointer)
   //ACL *aclPtr = NULL;
   //DWORD err = SetEntriesInAcl (1, &ea, NULL, &aclPtr);
   //if (err != ERROR_SUCCESS)  {
   //   SX_THROW ("Can't setup access control: " + sxstrerror (err));
   //}

   // --- prepare security descriptor
   //SECURITY_DESCRIPTOR *sdPtr = NULL;
   //sdPtr = (SECURITY_DESCRIPTOR*) LocalAlloc(LPTR,SECURITY_DESCRIPTOR_MIN_LENGTH);
   //InitializeSecurityDescriptor(sdPtr,SECURITY_DESCRIPTOR_REVISION);
   //SetSecurityDescriptorDacl(sdPtr,TRUE,aclPtr,FALSE);

   // --- prepare security attributes
   //SECURITY_ATTRIBUTES sa;
   //sa.nLength = sizeof(SECURITY_ATTRIBUTES);
   //sa.lpSecurityDescriptor = NULL;
   //sa.bInheritHandle = TRUE;

   // --- create the job object
   HANDLE hJobObj = jobObjectHndl = ::CreateJobObject(NULL, NULL);
   //HANDLE hJobObj = jobObjectHndl = ::CreateJobObject(&sa, NULL);
   if (hJobObj == NULL)  {
      DWORD err = ::GetLastError ();
      SX_THROW ("Can't create JobObject for child processes: " +
                sxstrerror (err));
   }
   
   // --- specify job object limits
   JOBOBJECT_EXTENDED_LIMIT_INFORMATION jobInfo;
   ::ZeroMemory (&jobInfo, sizeof (jobInfo));
   jobInfo.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE | JOB_OBJECT_LIMIT_BREAKAWAY_OK;
   // limit maximal application time if required
   if (obj->processTimeout > 0)  {
      jobInfo.BasicLimitInformation.PerJobUserTimeLimit.QuadPart = obj->processTimeout * 1e7; // time in 100 ns
      jobInfo.BasicLimitInformation.LimitFlags = jobInfo.BasicLimitInformation.LimitFlags | JOB_OBJECT_LIMIT_JOB_TIME;
   }
   if (!::SetInformationJobObject(hJobObj, JobObjectExtendedLimitInformation,
      &jobInfo, sizeof (jobInfo)))  {
      DWORD err = ::GetLastError ();
      SX_THROW ("Cannot set extended limit information for job object: "
                + sxstrerror (err));
   }

   jobSigPort = CreateIoCompletionPort (INVALID_HANDLE_VALUE, NULL, 0, 1);
   SX_CHECK (jobSigPort != NULL);

   JOBOBJECT_ASSOCIATE_COMPLETION_PORT port;
   port.CompletionKey = jobObjectHndl;
   port.CompletionPort = jobSigPort;
   SetInformationJobObject (jobObjectHndl, JobObjectAssociateCompletionPortInformation, &port, sizeof(port));

   // --- Clean up
   //if (aclPtr != NULL)  LocalFree(aclPtr);
   //if (sidPtr != NULL)  LocalFree(sidPtr);
   //if (sdPtr  != NULL)  LocalFree(sdPtr);
 }

void SxProcess::SxProcExecuter::removeJobObject ()
{
   SX_TRACE ();
   SX_CHECK (jobObjectHndl);
   UINT exitCode = ERROR_PROCESS_ABORTED;
   BOOL success = ::TerminateJobObject (jobObjectHndl, 0);
   if (!success) {
      DWORD err = ::GetLastError();
      SX_THROW("Can't terminate JobObject: " +
                sxstrerror(err));
   }
   CloseHandle (jobObjectHndl);
   jobObjectHndl = NULL;

   CloseHandle (jobSigPort);
   jobSigPort = NULL;

   SX_MUTEX (obj->statusMutex)  {
      obj->processRunning = false;
   }
   processExitStatus = static_cast<int>(exitCode);
}

void SxProcess::SxProcExecuter::removeProcess ()
{
   SX_TRACE ();
   UINT exitCode = ERROR_PROCESS_ABORTED;
   HANDLE proc = OpenProcess (SYNCHRONIZE | PROCESS_TERMINATE | PROCESS_QUERY_INFORMATION, false, pid);
   if (proc == NULL)  {
      DWORD err = GetLastError ();
      SX_DBG_MSG ("Get Process for pid " + SxString(pid)
                  + " failed: " + sxstrerror (err));
   }
   int success = ::TerminateProcess(proc, exitCode);
   if (!success)  {
      DWORD err = ::GetLastError();
      SX_DBG_MSG ("ERROR: Termination of child process after timeout failed: "
                  << sxstrerror(err));
      CloseHandle (proc);
      SX_THROW ("Failed to terminate process [" + SxString(pid) + "]: " + sxstrerror(err));
   }
   // --- TerminateProcess is async
   // WaitForSingleObject (hObject, INFINITE);
   DWORD error = ::WaitForSingleObject(proc, INFINITE);
   if (error != WAIT_OBJECT_0) {
      DWORD err = ::GetLastError();
      SX_DBG_MSG("ERROR: Waiting for child process to exit failed: "
                 << sxstrerror(err));
      processExitSignal = 1;
      processExitStatus = 1;
      SX_THROW ("Failed to wait for process [" + SxString(pid) + "] termination: " + sxstrerror(err));
   }
   CloseHandle (proc);
   processExitStatus = static_cast<int>(exitCode);
}

void SxProcess::SxProcExecuter::main ()
{
   SX_TRACE ();

   try {
      // --- Create job object for process collection and clean shutdown   
      createJobObject ();
 
      // --- create child process and wait until we recieve the pid
      createProcess (cmd);
   
      waitForPid ();
      SX_MUTEX (obj->statusMutex)  {
         obj->processRunning = true;
      }

      // --- waitForJob does not work for Windows 7:
      //     Windows 7 does not support nested jobs, so attaching jobs to
      //     job objects is not possible if the job has been assigned to a
      //     job object already. This, however hinders the kill tree approach.
      //     So we go back to the most easy approach and just wait for the 
      //     main process to end: NO ATOMATIC CLEANUP OF CHILD PROCESSES IN CASE OF CRASH OR 
      //     EXTERNAL SHUTDOWN!!!
      if (!obj->noHang)  {
         if (obj->killTree && IsWindows8OrGreater ())  {
            SX_DBG_MSG ("SxProcess waits for job object");
            waitForJob ();
         } else  {
            // wait for main process
            HANDLE proc = OpenProcess (SYNCHRONIZE | PROCESS_TERMINATE |
                                       PROCESS_QUERY_INFORMATION, false, pid);
            if (proc == NULL)  {
               DWORD err = GetLastError ();
               SX_THROW ("Get Process for pid " + SxString(pid)
                         + " failed: " + sxstrerror (err));
            }
            SX_DBG_MSG ("SxProcess waits for main process " 
                        << (long)proc << " with pid " << pid);
            waitForProcess (proc);
            CloseHandle (proc);
         } 
      }

      SX_DBG_MSG ("SxProcess finished");
      SX_MUTEX (obj->statusMutex)  {
         obj->processRunning = false;
      }

      // close pipe to make reader stop reading
      // This is why we need SxProcExecuter thread.
      obj->closePipeIdx (STDOUT_FILENO, 1);
   
   } catch (SxException e) {
      e.print ();
      SX_THROW (e, "Process interrupted unexpectedly");
   } catch (...)  {
      DWORD err = GetLastError ();
      SX_THROW ("Process interrupted unexpectedly: " + sxstrerror (err) 
               + "(" + sxprintf ("0x%x", err) + ")");
   }
}

void SxProcess::SxProcExecuter::closeProcessHandles (SxArray<HANDLE> &procs)
{
   SX_TRACE ();
   ssize_t nProcs = procs.getSize ();
   for (ssize_t iProc = 0; iProc < nProcs; ++iProc) {
      if (procs(iProc) != NULL) CloseHandle(procs(iProc));
   }
}

SxArray<HANDLE> SxProcess::SxProcExecuter::getProcessHandles (DWORD access)
{
   SX_TRACE ();
   SX_CHECK (jobObjectHndl);
   SxList<HANDLE> procs;
   // Get process IDs
   JOBOBJECT_BASIC_PROCESS_ID_LIST test;
   DWORD bufferSize;
   QueryInformationJobObject(jobObjectHndl, JobObjectBasicProcessIdList,
      &test, sizeof(test), &bufferSize);
   PJOBOBJECT_BASIC_PROCESS_ID_LIST procIDs = (PJOBOBJECT_BASIC_PROCESS_ID_LIST)malloc (bufferSize);
   QueryInformationJobObject(jobObjectHndl, JobObjectBasicProcessIdList,
      procIDs, bufferSize, NULL);

   DWORD nProcs = procIDs->NumberOfAssignedProcesses;
   DWORD nProcsInList = procIDs->NumberOfProcessIdsInList;
   if (nProcsInList < nProcs) {
      cout << "ERROR: Buffer for process list to small. We exit here!";
      SX_EXIT;
   }

   for (DWORD iProc = 0; iProc < nProcs; ++iProc) {
      DWORD pid = procIDs->ProcessIdList[iProc];
      HANDLE proc = OpenProcess (access, false, pid);
      if (proc == NULL) {
         DWORD err = GetLastError ();
         SX_DBG_MSG ("WARNING: Open process for pid " << SxString(pid)
                     <<" failed: " << err << ": " << sxstrerror (err));
      } else  {
         procs << proc;
      }
   }


   free(procIDs);

   return procs;
}

ssize_t SxProcess::SxProcExecuter::getNProcesses ()
{
   SX_TRACE ();
   // jobObjectHndl can be NULL
   JOBOBJECT_BASIC_PROCESS_ID_LIST procIDs;
   QueryInformationJobObject(jobObjectHndl, JobObjectBasicProcessIdList,
      &procIDs, sizeof(procIDs), NULL);

   return procIDs.NumberOfAssignedProcesses;
 
}

void SxProcess::SxProcExecuter::waitForProcess (HANDLE hProcess)
{
   SX_TRACE ();
   SX_CHECK(hProcess);

   // --- exit values like crash
   //processExitSignal = 1;
   processExitStatus = 1;

   // --- timeout
   DWORD timeout = INFINITE;
   if (obj->processTimeout > 0.) {
      timeout = static_cast<DWORD>(obj->processTimeout * 1e3);
   }

   // --- wait
   DWORD error = ::WaitForSingleObject(hProcess, timeout);
   if (error == WAIT_OBJECT_0) {
      // --- OK, get exit Code
      DWORD exitCode;
      BOOL status = ::GetExitCodeProcess(hProcess, &exitCode);
      if (!status) {
         DWORD err = ::GetLastError();
         SX_DBG_MSG ("ERROR: Can't get exit code of child process: "
                     << sxstrerror(err));
      }
      processExitStatus = static_cast<int>(exitCode);
   } else if (error == WAIT_FAILED) {
      DWORD err = ::GetLastError();
      SX_THROW("Waiting for child process failed: "
         + sxstrerror(err));
   } else {
      // WAIT_TIMEOUT
      UINT exitCode = ERROR_PROCESS_ABORTED;
      processExitStatus = static_cast<int>(exitCode);
      int success = ::TerminateProcess(hProcess, exitCode);
      if (!success)  {
         DWORD err = ::GetLastError();
         SX_THROW ("Termination of child process after timeout failed: "
            + sxstrerror(err));
      }
      // --- TerminateProcess is async
      // WaitForSingleObject (hObject, INFINITE);
     DWORD error = ::WaitForSingleObject(hProcess, INFINITE);
     if (error != WAIT_OBJECT_0) {
        DWORD err = ::GetLastError();
        SX_THROW ("Waiting for child process to exit after termination failed: "
                   + sxstrerror(err));
        
     }
   }
}

void SxProcess::SxProcExecuter::waitForJob ()
{
   SX_TRACE ();
   SX_CHECK (jobObjectHndl);
   SX_CHECK (jobSigPort);
   
   BOOL status = TRUE;
   DWORD sigCode = 0;

   ULONG_PTR completionKey;
   LPOVERLAPPED overlapped;

   // --- assume job has crashed
   processExitStatus = 1;

   // --- the following is based on a code snippet taken from
   //     https://flylib.com/books/en/4.419.1.42/1/
   //      which demonstrates how to wait for a job object

   // Wait as long as processes are running
   // Job Objects have there own time limit management
   // Runtime limits are set during job object creation-> see createJobObject
   do  {
      status = GetQueuedCompletionStatus (jobSigPort, &sigCode, &completionKey, &overlapped, INFINITE);
   } while ( status && sigCode != JOB_OBJECT_MSG_ACTIVE_PROCESS_ZERO);

   if (!status)  {
      DWORD err = ::GetLastError ();
      SX_THROW ("Unable to query job status: " + sxstrerror(err));
   }

   if (sigCode == JOB_OBJECT_MSG_ACTIVE_PROCESS_ZERO)  { // process exited normally
      // --- OK, get exit Code of first process
      DWORD exitCode;
      BOOL status = ::GetExitCodeProcess (mainProc, &exitCode);
      if (!status) {
         DWORD err = ::GetLastError ();
         CloseHandle (mainProc);
         mainProc = NULL;
         SX_THROW ("Failed to get exit code of child process: "
                     + sxstrerror(err));
      }
      processExitStatus = static_cast<int>(exitCode);

   }

   CloseHandle (mainProc);
   mainProc = NULL;
     
}

pid_t SxProcess::SxProcExecuter::waitForPid ()
{
   SX_TRACE ();
   // --- wait until CreateProcess starts new process and obtains pid
   pid_t result = 0;
   SX_MUTEX (mutex)  {
      while (!gotPid)  {
         condition.wait (&mutex);
      }
      result = pid;
   }
   return result;
}

#endif /* WIN32 */

#ifndef WIN32
SxProcess::SxProcessTimer::SxProcessTimer (pid_t pid_, bool killTree_)
   : pid(pid_),
     n(0),
     killTree(killTree_)
{
   // empty
}

void SxProcess::SxProcessTimer::interval ()
{
   if (n == 0)  {
      try { SxProcess::kill (pid, SIGTERM, killTree); }
      catch (...) { }
   }  else  {
      try { SxProcess::kill (pid, SIGKILL, killTree); }
      catch (...) { }
      stop ();
   }
   
   n++;
}
#endif /* not WIN32 */

// --------------------------------------------------------------------------

SxProcess::SxProcess (Mode mode_)
   : channels (SxProcess::None),
     useBuffers (SxProcess::None),
     mode (mode_),
     sync (SxProcess::NoSync),
     pid(0),
     killTree (true),
     processRunning (false),
     lineWise(true),
     noHang(false),
     processExitStatus (0),
     processExitSignal (0),
     bufLen (10240),
#    ifdef WIN32
        hiddenWindow (false),
#    endif
     processTimeout(-1.),
     eofStdErr(false),
     eofStdOut(false)
{
   srand((unsigned int)time(0));
   initPid ();
   initPipes ();
}
 
 
SxProcess::SxProcess (const SxString &cmdLine, Mode mode_)
   : channels (SxProcess::None),
     useBuffers (SxProcess::None),
     mode (mode_),
     sync (SxProcess::NoSync),
     pid(0),
     killTree (true),
     processRunning (false),
     lineWise(true),
     noHang(false),
     processExitStatus (0),
     processExitSignal (0),
     bufLen (10240),
#    ifdef WIN32
        hiddenWindow (false),
#    endif
     processTimeout(-1.),
     eofStdErr(false),
     eofStdOut(false)
{ 
   srand((unsigned int)time(0));
   initPid ();
   initPipes ();
   setCommandLine (cmdLine); 
} 
 
 
SxProcess::SxProcess (const SxString &cmd_, const SxList<SxString> &argList_,
                      Mode mode_)
   : channels (SxProcess::None),
     useBuffers (SxProcess::None),
     mode (mode_),
     sync (SxProcess::NoSync),
     pid(0),
     killTree (true),
     processRunning (false),
     lineWise(true),
     noHang(false),
     processExitStatus (0),
     processExitSignal (0),
     bufLen (10240),
#    ifdef WIN32
        hiddenWindow (false),
#    endif
     processTimeout(-1.),
     eofStdErr(false),
     eofStdOut(false)
{ 
   srand((unsigned int)time(0));
   initPid ();
   initPipes ();
   setCommand (cmd_);
   setArguments (argList_); 
} 
 
 
SxProcess::~SxProcess ()
{ 
   closePipes ();
} 
 
 
SxString SxProcess::getenv (const SxString &varName)
{ 
   return SxString (::getenv (varName.getElems ()) ); 
} 
 
 
void SxProcess::setenv (const SxString &varName, const SxString &value)
{
#ifdef WIN32
   sxputenv (const_cast<char *>((varName.trim() + "=" + value).getElems ()));
#else
   const int OVERWRITE = 1; 
 ::setenv (varName.getElems (), value.getElems (), OVERWRITE); // BSD 4.3 
// putenv (static_cast<char *>(varName.trim() + "=" + value).getElems ()); //POSIX 
#endif
}

void SxProcess::daemon (int umask_, const SxString &chdir_,
                        const SxString &logfile_)
{
#ifdef WIN32
   SX_EXIT; // daemons are extra windows services
#else
   if (umask_ < 0)  {
      SX_DBG_MSG("ERROR: Invalid umask value " << umask_);
      SX_QUIT; 
   }

   struct rlimit rl;
   if (getrlimit (RLIMIT_NOFILE, &rl) < 0)  {
      int err = errno;
      SX_DBG_MSG("ERROR: getrlimit(RLIMIT_NOFILE) failed: "
                 << sxstrerror (err));
      SX_QUIT;
   }

   // --- create daemon 
   pid_t pid = fork ();
   if (pid < 0)  {
      int err = errno;
      SX_DBG_MSG("ERROR: Can't fork daemon: fork() failed: "
                 << sxstrerror (err));
      SX_QUIT;
   }

   if (pid > 0)  {
      // --- parent
      exit (EXIT_SUCCESS);
   }

   pid_t sid = setsid ();
   if (sid < 0)  {
      int err = errno;
      SX_DBG_MSG("ERROR: Can't create process group: setsid() failed: "
                 << sxstrerror (err));
      SX_QUIT;
   }

   // --- file creation mask
   umask (static_cast<mode_t>(umask_));
   
   // --- daemon working directory
   if (chdir_ != "")  {
      if (::chdir(chdir_.getElems()) != 0) {
         int err = errno;
         SX_THROW("Can't change the current directory to '"
                  + chdir_ + "': " + sxstrerror (err));
      }
   }

   // --- close all file descriptors
   int nFile = 1024;
   if (rl.rlim_max != RLIM_INFINITY)  {
      // rlim_t rl.rlim_max: unsigned long long or unsigned long
      nFile = static_cast<int>(rl.rlim_max);
   }
   for (int i=0; i < nFile; i++)  {
      close (i);
   }

   // --- open new file descriptors
   SxString output = (logfile_ != "") ? logfile_ : "/dev/null";

   int fd0 = open ("/dev/null", O_RDONLY);
   int fd1 = open (output.getElems (), O_WRONLY|O_CREAT|O_APPEND,
      S_IRUSR|S_IWUSR);
   int fd2 = (fd1 == STDOUT_FILENO) ? dup(fd1) : -1;

   if (   fd0 != STDIN_FILENO
       || fd1 != STDOUT_FILENO
       || fd2 != STDERR_FILENO)
   {
      SX_DBG_MSG("ERROR: Can't create daemon stdin|stdout|stderr");
      SX_QUIT;
   }
#endif /* WIN32 */
}


void SxProcess::initPid ()
{
#  ifdef WIN32
      pid = GetCurrentProcessId ();
#  else
      pid = getpid ();
#  endif
}

void SxProcess::initPipes ()
{
#ifdef WIN32
   for (uint8_t i = STDIN_FILENO; i <= STDERR_FILENO; i++)  {
      pipes[i].handles[0] = pipes[i].handles[1] = INVALID_HANDLE_VALUE;
   }
#else
   for (uint8_t i = STDIN_FILENO; i <= STDERR_FILENO; i++)  {
      pipes[i].fds[0] = pipes[i].fds[1] = -1;
   }
#endif
}

pid_t SxProcess::getPid ()
{
   return pid;
}

pid_t SxProcess::getPPid ()
{
#  ifdef WIN32
      //msdn.microsoft.com/en-us/library/windows/desktop/ms686701(v=vs.85).aspx
      //     Take a snapshot of all processes in the system
      //     "32" part in names can/cannot be ignored on 64
      pid_t result = 0;
      HANDLE hProcessSnap;
      hProcessSnap = CreateToolhelp32Snapshot (TH32CS_SNAPPROCESS, 0);
      if (hProcessSnap == INVALID_HANDLE_VALUE)  {
         DWORD err = ::GetLastError ();
         SX_THROW ("CreateToolhelp32Snapshot(TH32CS_SNAPPROCESS,0) failed: "
                   + sxstrerror (err));
      }

      PROCESSENTRY32 pe32;
      pe32.dwSize = sizeof(PROCESSENTRY32);
      if (!Process32First (hProcessSnap, &pe32))  {
         DWORD err = ::GetLastError ();
         CloseHandle (hProcessSnap);
         SX_THROW ("Process32First() failed: " + sxstrerror (err));
      }

      // --- walk the snapshot of processes and find myPid process
      DWORD myPid = GetCurrentProcessId ();
      do {
         if (myPid == pe32.th32ProcessID)  {
            result = pe32.th32ParentProcessID;
            break;
         }
      } while (Process32Next (hProcessSnap, &pe32));
      CloseHandle (hProcessSnap);

      return result;
#  else
      return getppid ();
#  endif /* WIN32 */
}


#ifdef WIN32
pid_t SxProcess::getChildProcPID ()
{
   SX_CHECK (childProc);
   return childProc->pid;
}

pid_t SxProcess::getChildProcThread0ID ()
{
   SX_CHECK (childProc);
   return childProc->thread0ID;
}
#endif

void SxProcess::setKillTree (bool enabled_)
{
   killTree = enabled_;
}

void SxProcess::setWorkingDirectory (const SxString &workDir_)
{
   workDir = workDir_;
}

void SxProcess::setCommunication (Channels channels_)
{
   channels = channels_;
}

void SxProcess::enableBuffer (Channels channels_)
{
   setCommunication (channels_);
   useBuffers = channels_;
}


void SxProcess::synchronize (int mask)
{
   if      (mask ==  None)              {  sync = NoSync;  }
   else if (mask == (StdOut | StdErr))  {  sync = SyncOutErr; }
   else  {
      SX_EXIT; // it makes only sense to combine STDOUT and STDERR
   }
}


void SxProcess::setCommand (const SxString &cmd_)
{
   cmd = cmd_;
}

void SxProcess::setCommandLine (const SxString &cmdLine)
{
   removeAllArguments ();

   SxList<SxString>::Iterator it;
   SxList<SxString> args = cmdLine.tokenize(' ');
   SX_CHECK (args.getSize() > 0);
   it = args.begin();
   cmd = *it++;
   while (it != args.end())  {
      argList << *it++;
   }
}


void SxProcess::setBufferSize (ssize_t newSize)
{
   bufLen = newSize;
}


void SxProcess::removeAllArguments ()
{
   argList.removeAll ();
}


void SxProcess::setArguments (const SxList<SxString> &argList_)
{
   argList = argList_;
}


void SxProcess::addArgument (const SxString &arg)
{
   argList << arg;
}

SxList<SxString> SxProcess::getArguments () const
{
   return argList;
}

void SxProcess::setTimeout (double processTimeout_)
{
   processTimeout = processTimeout_;
}

int SxProcess::start (const SxList<SxString> &env)
{
   run (env);
   read ();
   return wait ();
}

#ifdef WIN32
unsigned char SxProcess::getRndByte () const
{
   unsigned int rndVal = 0;
   unsigned int limit = RAND_MAX - (RAND_MAX % 256);

   do {
      rndVal = rand ();

   } while (rndVal >= limit);

   return rndVal % 256;
}

SxString SxProcess::getUUIDv4 ()
{
   ssize_t i;
   const ssize_t nBytes = 16;
   SxArray<unsigned char> data(nBytes);
   for (i=0; i < 16; ++i)  {  // generate 128 random bits
      data(i) = getRndByte ();
   }
   // --- RFC 4122, 4.4
   data(6) = 0x40 | (data(6) & 0x0f);
   data(8) = 0x80 | (data(8) & 0x3f);
   SxString res; // "XXXXXXXX-XXXX-4XXX-XXXX-XXXXXXXXXXXX";
   for (i = 0; i < 16; ++i) {
      res += SxString::sprintf ("%02X", data(i));
      if (i == 3 || i == 5 || i == 7 || i == 9) {
         res += "-";
      }
   }
   return res;
}
#endif
void SxProcess::run (const SxList<SxString> &env)
{
   if (killTree && (mode == SxProcess::Daemon))  {
      SX_THROW ("The combination 'killTree and Mode::Daemon' makes no sense");
   }
   pid = 0;
   processExitStatus = 0;
   processExitSignal = 0;
   outBuffer = errBuffer = "";
   eofStdOut=false;
   eofStdErr=false;

   // --- noHang
   //     cmd_ is cmd without optional '&' no hang symbol
   //     argList_ is argList without optional '&'
   noHang = false;
   SxString cmd_ = cmd;
   SxList<SxString> argList_ = argList;
   if (argList_.getSize () > 0) {
      if (argList_.last().tail(1) == "&")  {
         noHang = true;
         argList_.last().removeLast ();
         if (argList_.last() == "")  {
            argList_.removeLast ();
         }
      }
   }  else if (cmd_.tail(1) == "&") {
      noHang = true;
      cmd_.removeLast ();
   }
   if (noHang)  {
      // --- allow only Stdin, NOHANG would otherwise block and wait on read
      if (channels & Stdin) channels = Stdin;
      else                  channels = None;
   }
   if (cmd_ == "")  {
      SX_THROW ("SxProcess::run() command can't be empty string");
   }

   bool usePipes = (channels & (StdIn | StdOut | StdErr)) != 0;

   SxString          key, value;
   SxList<SxString>  keyValuePair;
 
   // --- save environment variables of current process
   SxString binPath (getenv("PATH"));
   SxString libPath (getenv("LD_LIBRARY_PATH"));

   // --- setup process's environment variables
   SxList<SxString>::ConstIterator it;
   for (it = env.begin(); it != env.end(); ++it)  {
      keyValuePair = (*it).tokenize('=');

      // --- check 'key=value'
      SX_CHECK (keyValuePair.getSize() == 2, keyValuePair.getSize());

      setenv (key, value);
   }
   if (binPath.getSize() > 0)
      setenv ("PATH", binPath);
#  if defined WIN32
      // libPath = binPath
#  elif defined __APPLE__
      if (libPath.getSize() > 0)
         setenv ("DYLD_LIBRARY_PATH", libPath);
#  else
      if (libPath.getSize() > 0)
         setenv ("LD_LIBRARY_PATH", libPath);
#  endif

   // --- change working directory for child
   if (workDir != "")  {
#  ifdef WIN32
      if (workDir.isUnicode()) {
         if (::_wchdir ((LPCWSTR)workDir.utf16 ().elements) != 0) {
            int err = errno;
            SX_THROW("Can't change the current directory to '"
                     + workDir + "'. ");
         }
      } else {
         if (::chdir (workDir.ascii ()) != 0) {
            int err = errno;
            SX_THROW("Can't change the current directory to '"
                     + workDir + "'. ");
         }
      }
#  else
      if (::chdir(workDir.getElems()) != 0) {
         int err = errno;
         SX_THROW("Can't change the current directory to '"
                  + workDir + "': " + sxstrerror (err));
      }
#  endif
   }

#  ifdef WIN32
      // --- Windows allows spaces: if spaces are used we must quote
      //     the command or the elements of the argument list
      if (cmd_.contains (" ")) cmd_ = "\"" + cmd_ + "\"";
      for (ssize_t i = 0; i < argList_.getSize (); ++i)  {
         if (argList_(i).contains (" ")) argList_(i) = "\"" + argList_(i) + "\"";
      }
      // --- only cmd_ can be passed to CreateProcess
      //     argList could be set with setArguments()
      if (argList_.getSize () > 0) {
         cmd_ += " " + SxString::join (argList_, " ");
      }

      SECURITY_ATTRIBUTES saAttr;
      ::ZeroMemory(&saAttr, sizeof(saAttr));

      PSID pAdminSID = NULL;
      PACL pACL = NULL;
      PSECURITY_DESCRIPTOR pSD = NULL;
      EXPLICIT_ACCESS ea[1];
      SID_IDENTIFIER_AUTHORITY SIDAuthNT = SECURITY_NT_AUTHORITY;
      ZeroMemory (&ea, sizeof(EXPLICIT_ACCESS));
      if (!AllocateAndInitializeSid(&SIDAuthNT, 2,
                                    SECURITY_BUILTIN_DOMAIN_RID,
                                    DOMAIN_ALIAS_RID_ADMINS,
                                    0, 0, 0, 0, 0, 0, &pAdminSID)) {

         SX_THROW("PipeFailed", "Can't Allocate SID");
      }
      // Allow admins all access
      ea[0].grfAccessPermissions = KEY_ALL_ACCESS;
      ea[0].grfAccessMode = SET_ACCESS;
      ea[0].grfInheritance = NO_INHERITANCE;
      ea[0].Trustee.TrusteeForm = TRUSTEE_IS_SID;
      ea[0].Trustee.TrusteeType = TRUSTEE_IS_GROUP;
      ea[0].Trustee.ptstrName = (LPSTR)pAdminSID;

      DWORD dRes = SetEntriesInAcl(1, ea, NULL, &pACL);
      if (dRes != ERROR_SUCCESS) {
         SX_THROW("PipeFailed", "Can't set ACL entries for security descriptor");
      }
      pSD = (PSECURITY_DESCRIPTOR) LocalAlloc(LPTR,
                                              SECURITY_DESCRIPTOR_MIN_LENGTH);
      if (pSD == NULL)
         SX_THROW("PipeFailed", "Can't allocate security descriptor");
      if (!InitializeSecurityDescriptor(pSD,
                                        SECURITY_DESCRIPTOR_REVISION)) {
         SX_THROW("PipeFailed", "Can't initialize security descriptor");
      }

      saAttr.nLength = sizeof (SECURITY_ATTRIBUTES);
      saAttr.bInheritHandle = TRUE;
      saAttr.lpSecurityDescriptor = pSD;
      SxString name;

      if (usePipes)  {
         if ((channels & StdErr) && !combinedStdout()) {
            ::ZeroMemory(&pipes[STDERR_FILENO].oOverlap, sizeof(OVERLAPPED));
            pipes[STDERR_FILENO].oOverlap.hEvent = CreateEvent(NULL, FALSE,
                                                               TRUE, NULL);
            if (pipes[STDERR_FILENO].oOverlap.hEvent == NULL) {
               SX_THROW("PipeFailed", "Create event failed for STDERR");
            }
            name = "\\\\.\\pipe\\";
            name += getUUIDv4 ();
            LPTSTR pipeName = TEXT(name.elements);
            HANDLE tmpErrHandle = CreateNamedPipe(pipeName,
                                                  PIPE_ACCESS_INBOUND
                                                  | FILE_FLAG_OVERLAPPED,
                                                  PIPE_TYPE_BYTE |
                                                  PIPE_READMODE_BYTE|PIPE_WAIT,
                                                  1, BUFSIZ*sizeof(TCHAR),
                                                  BUFSIZ*sizeof(TCHAR),
                                                  0, &saAttr);
            if (tmpErrHandle == INVALID_HANDLE_VALUE) {
               SX_THROW("PipeFailed", "Can't create STDERR pipe.");
            }
            pipes[STDERR_FILENO].handles[1] = CreateFile(pipeName,
                                                         FILE_WRITE_DATA
                                                         |SYNCHRONIZE,
                                                         0,
                                                         &saAttr,
                                                         OPEN_EXISTING,
                                                         FILE_ATTRIBUTE_NORMAL,
                                                         0);
            ::ZeroMemory(&pipes[STDERR_FILENO].handles[0], sizeof(HANDLE));
            if (!DuplicateHandle(GetCurrentProcess(), tmpErrHandle,
                                 GetCurrentProcess(),
                                 &(pipes[STDERR_FILENO].handles[0]),
                                 0, FALSE, DUPLICATE_SAME_ACCESS)) {
               SX_THROW("PipeFailed", "Can't duplicate STDERR pipe handle");
            }
            if (!CloseHandle(tmpErrHandle))
               SX_THROW("PipeFailed", "Close handle failed for STDERR");
         
         }
         if (channels & StdOut) {
            // --- create pipe for child's STDOUT
            ::ZeroMemory(&pipes[STDOUT_FILENO].oOverlap, sizeof(OVERLAPPED));
            pipes[STDOUT_FILENO].oOverlap.hEvent = CreateEvent(NULL, FALSE,
                                                               TRUE, NULL);
            if (pipes[STDOUT_FILENO].oOverlap.hEvent == NULL) {
               SX_THROW("PipeFailed", "Create event failed for STDOUT");
            }
            name = "\\\\.\\pipe\\";
            name += getUUIDv4 ();
            LPTSTR pipeName = TEXT(name.elements);
            HANDLE tmpOutHandle = CreateNamedPipe(pipeName,
                                                  PIPE_ACCESS_INBOUND |
                                                  FILE_FLAG_OVERLAPPED,
                                                  PIPE_TYPE_BYTE |
                                                  PIPE_READMODE_BYTE | PIPE_WAIT,
                                                  1, BUFSIZ * sizeof(TCHAR),
                                                  BUFSIZ * sizeof(TCHAR),
                                                  0, &saAttr);
            if (tmpOutHandle == INVALID_HANDLE_VALUE) {
               SX_THROW("PipeFailed", "Can't create STDOUT pipe.");
            }
            pipes[STDOUT_FILENO].handles[1] = CreateFile(pipeName,
                                                         FILE_WRITE_DATA |
                                                         SYNCHRONIZE,
                                                         0,
                                                         &saAttr,
                                                         OPEN_EXISTING,
                                                         FILE_ATTRIBUTE_NORMAL,
                                                         0);
            if (!DuplicateHandle(GetCurrentProcess(), tmpOutHandle,
                                 GetCurrentProcess(),
                                 &(pipes[STDOUT_FILENO].handles[0]),
                                 0, FALSE, DUPLICATE_SAME_ACCESS)) {
               SX_THROW("PipeFailed", "Can't duplicate STDOUT pipe handle");
            }
            if (!CloseHandle(tmpOutHandle))
               SX_THROW("PipeFailed", "Close handle failed for STDOUT");
         }
         if (channels & StdIn) {
            // --- create pipe for child's STDIN
            ::ZeroMemory(&pipes[STDIN_FILENO].oOverlap, sizeof(OVERLAPPED));
            pipes[STDIN_FILENO].oOverlap.hEvent = CreateEvent(NULL, FALSE,
                                                              TRUE, NULL);
            if (pipes[STDIN_FILENO].oOverlap.hEvent == NULL) {
               SX_THROW("PipeFailed", "Create event failed for STDIN");
            }
            name = "\\\\.\\pipe\\";
            name += getUUIDv4 ();
            LPTSTR pipeName = TEXT(name.elements);
            HANDLE tmpInHandle = CreateNamedPipe(pipeName,
                                                 PIPE_ACCESS_OUTBOUND |
                                                 FILE_FLAG_OVERLAPPED,
                                                 PIPE_TYPE_BYTE |
                                                 PIPE_READMODE_BYTE | PIPE_WAIT,
                                                 1, BUFSIZ * sizeof(TCHAR),
                                                 BUFSIZ * sizeof(TCHAR),
                                                 0, &saAttr);
            if (tmpInHandle == INVALID_HANDLE_VALUE) {
               SX_THROW("PipeFailed", "Can't create STDIN pipe.");
            }

            pipes[STDIN_FILENO].handles[0] = CreateFile(pipeName,
                                                        FILE_READ_DATA |
                                                        SYNCHRONIZE,
                                                        0,
                                                        &saAttr,
                                                        OPEN_EXISTING,
                                                        FILE_ATTRIBUTE_NORMAL,
                                                        0);
            if (!DuplicateHandle(GetCurrentProcess(), tmpInHandle,
                                 GetCurrentProcess(),
                                 &(pipes[STDIN_FILENO].handles[1]), 0,
                                 FALSE, DUPLICATE_SAME_ACCESS)) {
               SX_THROW("PipeFailed", "Can't duplicate STDIN pipe handle");
            }
            if (!CloseHandle(tmpInHandle))
               SX_THROW("PipeFailed", "Close handle failed for STDIN");
         }

      }

      childProc = SxPtr<SxProcExecuter>::create (this, cmd_);
      childProc->hiddenWindow = hiddenWindow;
      childProc->start ();
      pid = childProc->waitForPid ();

      SX_MUTEX (statusMutex)  {
         processRunning = true; 
      }

      if (channels & Stdin)  {  // --- write to child's STDIN channel
         write (writeToStdin ());
      }

#else /* WIN32 */

   // --- Unix

   // --- process setup
   //     (C-like setup requires conventional memory handling)
   char **argv = NULL;
   argv = new char * [static_cast<size_t>(argList_.getSize()) + 2];
   int i=0;
   argv[i] = new char [static_cast<size_t>(cmd_.getSize()) + 1];
   sprintf (argv[i++], "%s", cmd_.getElems ());

   for (it = argList_.begin(); it != argList_.end(); ++it)  {
      argv[i] = new char [static_cast<size_t>((*it).getSize()) + 1];
      sprintf (argv[i++], "%s", (*it).getElems ());
   }
   argv[i] = 0;

   // --- Combine C-like with C++-like error handling
   //     (see 'goto error' statements)
   SxString exceptionMsg;
   SxString exceptionFile;
   int exceptionLine = 0;

   // --- create pipes for desired standard channels
   if (usePipes)  {
      for (i = STDIN_FILENO; i <= STDERR_FILENO; i++)  {
         if ( (i == STDERR_FILENO) && (combinedStdout ()) )  break; // done
         if ( (channels & (1 << i)) && (::pipe (pipes[i].fds) != 0) )  {
            int err = errno;
            exceptionMsg  = "Can't create POSIX pipe: pipe () failed: "
                          + sxstrerror (err);
            exceptionFile = __FILE__;
            exceptionLine = __LINE__;
            goto error;  // C-like error handling
         }
      }
   }

   processRunning = true;

   // --- fork process into parent and child process
   pid = fork ();
   if (pid == -1)  {
      int err = errno;
      exceptionMsg  = "Can't fork new process: fork() failed: "
                    + sxstrerror (err);
      exceptionFile = __FILE__;
      exceptionLine = __LINE__;
      goto error;  // C-like error handling
   }

   // --- start process timeout
   if (pid > 0 && processTimeout > 0.)  {
      processTimer = SxPtr<SxProcessTimer>::create (pid, killTree);
      processTimer->setInterval (processTimeout);
      processTimer->start ();
   }

   // --- daemonize process?
   if (mode == Daemon)  {
      pid_t sid = setsid ();
      if (sid == -1)  {
         int err = errno;
         exceptionMsg  = "Can't set session id of new process: "
                         "setsid() failed: " + sxstrerror (err);
         exceptionFile = __FILE__;
         exceptionLine = __LINE__;
         goto error;  // C-like error handling
      }
   }

   // --- child process: execute external program
   if (pid == 0)  {
      // --- close the pipe ends which are only used by the parent process
      if (channels & StdIn)  {
         closePipeIdx (STDIN_FILENO, 1); // child doesn't write to its stdin
         dup2 (pipes[STDIN_FILENO].fds[0], STDIN_FILENO); // the reading end
      }
      if (channels & StdOut)  {
         closePipeIdx (STDOUT_FILENO, 0); // child doesn't read from its stdout
         dup2 (pipes[STDOUT_FILENO].fds[1], STDOUT_FILENO); // the writing end
      }
      if (channels & StdErr)  {
         if (combinedStdout ())  { // just dup the stdout channel
            dup2 (pipes[STDOUT_FILENO].fds[1], STDERR_FILENO);
         }  else  { // we want a separate channel
            closePipeIdx (STDERR_FILENO, 0); // child doesn't read from its stderr
            dup2 (pipes[STDERR_FILENO].fds[1], STDERR_FILENO); // the writing end
         }
      }

      // --- execute the program
      execvp (cmd_.getElems (), argv);

      // --- we never should arrive here
      int err = errno;
      exceptionMsg  = "Can't execute program '" + cmd_
                    + "': execvp() failed: " + sxstrerror (err);
      exceptionFile = __FILE__;
      exceptionLine = __LINE__;
      goto error;  // C-like error handling
   }

   // --- parent process: read from / write to channels of child process;
   //     first, close all pipe ends which are only used by the child process
   if (channels & StdIn) // parent doesn't read from child stdin
      closePipeIdx (STDIN_FILENO, 0);
   if (channels & StdOut) // parent doesn't write to child stdout
      closePipeIdx (STDOUT_FILENO, 1);
   if ( (channels & StdErr) && (!combinedStdout ()) )  {
      closePipeIdx (STDERR_FILENO, 1); // parent doesn't write to child stderr
   }

   if (channels & Stdin)  {  // --- write to child's STDIN channel
      write (writeToStdin ());
   }
#endif /* WIN32 */

#  ifndef WIN32
error:
   // --- cleanup memory
   for (i=0; i < argList_.getSize()+2; ++i)  delete [] argv[i];
   delete [] argv;

   // --- did an error occur or do we just clean up?
   if (exceptionMsg != "")  {
       if (pid == 0)  exit (128); // the OS cleans up the child process nicely
       else           throw SxException (exceptionMsg.getElems (),
                                         exceptionFile.getElems (),
                                         exceptionLine);
      // The parent process pipes are closed in the SxProcess destructor.
   }
   if (pid == 0)  exit(0);
#  endif
}

void SxProcess::write (const SxString &data)
{
   write (data.elements, data.getNBytes ());
}

void SxProcess::write (const char *data, ssize_t origNBytes)
{
   SX_CHECK (channels & Stdin);
   if (origNBytes < 1)  return;
#  ifdef WIN32
      DWORD nBytes;
      BOOL status = WriteFile (pipes[STDIN_FILENO].handles[1], data,
                               static_cast<DWORD>(origNBytes), &nBytes,
                               NULL);
     if (!status) {
       SX_THROW("PipeFailed", "Can't write to pipe: " + sxstrerror(GetLastError()));
     }
     // close the pipe to cause EOF
     ::CloseHandle (pipes[STDIN_FILENO].handles[1]);

      //if (!status)  {
      //   char winmsg[256];
      //   sxgetWinError (GetLastError (), winmsg, 256);
      //   sxprintf ("Problem writing to pipe: %s\n", winmsg);
      //}
#  else /* WIN32 */
      int fd = pipes[STDIN_FILENO].fds[1];
      ::write (fd, data, static_cast<size_t>(origNBytes));
      // close the pipe to cause EOF
      closePipe (STDIN_FILENO, parentEnd);

      //ssize_t nBytes = ::write (fd, data, (size_t)origNBytes);
      //if (nBytes != origNBytes)  {
      //   printf ("Problem writing to pipe\n");
      //}
#  endif /* WIN32 */
}

void SxProcess::readBufferLineWise (int FD, SxPtr<SxArray<char>> *bufferPtr,
                                    SxPtr<ssize_t> *curIdxPtr)
{
   ssize_t rc = 0;
   char c = 0;
   const ssize_t MAX_N_CHARS_FILE = 2147483647; // see SxString
#ifndef WIN32
   int fd = pipes[FD].fds[0];
   // make read non-blocking
   fcntl (fd, F_SETFL, fcntl (fd, F_GETFL, 0) | O_NONBLOCK);
#endif

   do  {
#     ifdef WIN32
         rc = 0;
         DWORD nr;
         // Try to ReadFile without actually reading to see 
         // if there is something to read
         BOOL status;
         status = ReadFile(pipes[FD].handles[0], &c, 0, &nr, &pipes[FD].oOverlap);
         if (!status) {
            DWORD err = ::GetLastError();
            if (err == ERROR_IO_PENDING) {
               // IO Pending
               return;
            } else if (err == ERROR_BROKEN_PIPE) {
               if (FD == STDERR_FILENO)  { eofStdErr = true; }
               else                      { eofStdOut = true; }
            } else {
               SX_THROW("PipeFailed", "Can't read from pipe: " + sxstrerror(err));
            }
         } else {
            status = ReadFile(pipes[FD].handles[0], &c, 1, &nr, &pipes[FD].oOverlap);
            // --- possible (!status) is 'The pipe has been ended'
            if (status) {
               rc = static_cast<ssize_t>(nr);
            } else {
               DWORD err = ::GetLastError();
               if (err == ERROR_IO_PENDING) {
                  // IO Pending
                  return;
               } else if (err == ERROR_BROKEN_PIPE) {
                  if (FD == STDERR_FILENO)  { eofStdErr = true; }
                  else                      { eofStdOut = true; }
               } else {
                  SX_THROW("PipeFailed", "Can't read from pipe: " + sxstrerror(err));
               }
            }
         }
#     else /* WIN32 */
         rc = ::read (fd, &c, 1);
         if (rc == -1) {
            // IO is pending
            return;
         } else if (rc == 0) {
            if (FD == STDERR_FILENO) { eofStdErr = true; }
            else                     { eofStdOut = true; }
         }
#     endif /* WIN32 */

      if (rc == 1 && c != '\n') {
         if (*(*curIdxPtr) >= MAX_N_CHARS_FILE) {
            SX_THROW ("PipeFailed", "Line from stdout/stderr is too long "
                      "(max " + SxString(MAX_N_CHARS_FILE) + ").");
         }
         if ((*bufferPtr)->getSize () <= *(*curIdxPtr)) {
            (*bufferPtr)->resize (*(*curIdxPtr) + bufLen, true);
         }
         (*(*bufferPtr))(*(*curIdxPtr)) = c;
         (*(*curIdxPtr))++;
      } else if (rc == 1 || *(*curIdxPtr) > 0) {
         if ( (*(*curIdxPtr) > 0) &&
              ((*bufferPtr)->elements[(*(*curIdxPtr)) - 1] == '\r') ) {
            // --- simplify CR+LF to '\n'
            (*(*curIdxPtr))--;
         }
         SxString str((*bufferPtr)->elements, (*(*curIdxPtr)));
         if (FD == STDERR_FILENO)  readFromStderr (str);
         else                      readFromStdout (str);
         (*(*curIdxPtr)) = 0;
      }
   } while (rc == 1);

}

void SxProcess::readLineWise (bool useStderr)
{
   SX_CHECK (bufLen > 0, bufLen);

   SxPtr<SxArray<char>> outBufferPtr, errBufferPtr;
   if (channels & StdOut)
      outBufferPtr = SxPtr<SxArray<char>>::create (bufLen);
   if (useStderr)
      errBufferPtr = SxPtr<SxArray<char>>::create (bufLen);
   SxPtr<ssize_t> outIdxPtr, errIdxPtr;
   outIdxPtr = SxPtr<ssize_t>::create (0);
   errIdxPtr = SxPtr<ssize_t>::create (0);

#  ifdef WIN32
      DWORD nCount = 0;
      if (channels & StdOut) {
         hEvents[nCount] = pipes[STDOUT_FILENO].oOverlap.hEvent;
         nCount++;
      }
      if (useStderr) {
         hEvents[nCount] = pipes[STDERR_FILENO].oOverlap.hEvent;
         nCount++;
      }
      while (1) {
         if (useStderr) {
            if ((channels & StdOut)) {
               if (eofStdOut && eofStdErr)
                  break;
            } else if (eofStdErr) {
               break;
            }
         } else {
            if ((channels & StdOut) && eofStdOut)  break;
         }
         DWORD dwWait = WaitForMultipleObjects(nCount,
                                               hEvents,
                                               FALSE,
                                               INFINITE);
         DWORD i = dwWait - WAIT_OBJECT_0;
         if (i > (nCount - 1))
            SX_THROW("PipeFailed", "Signaled index out of range");
         if ((channels & StdOut) && i == 0 && !eofStdOut)
            readBufferLineWise(STDOUT_FILENO, &outBufferPtr, &outIdxPtr);
         if (useStderr && i == (nCount-1) && !eofStdErr)
            readBufferLineWise(STDERR_FILENO, &errBufferPtr, &errIdxPtr);
      }
#  else
      struct timeval tv;
      fd_set set;
      int res = 0;
      while (1) {
         if (useStderr) {
            if ((channels & StdOut)) {
               if (eofStdOut && eofStdErr)
                  break;
            } else if (eofStdErr) {
               break;
            }
         } else {
            if ((channels & StdOut) && eofStdOut)  break;
         }
         tv.tv_sec = 1;
         tv.tv_usec = 0;
         FD_ZERO(&set);
         if (channels & StdOut)
            FD_SET(pipes[STDOUT_FILENO].fds[0], &set);
         if (useStderr)
            FD_SET(pipes[STDERR_FILENO].fds[0], &set);
         res = select (FD_SETSIZE, &set, NULL, NULL, &tv);
         if (res >= 1) {
            if ((channels & StdOut) && !eofStdOut)
               readBufferLineWise (STDOUT_FILENO,
                                   &outBufferPtr,
                                   &outIdxPtr);

            if (useStderr && !eofStdErr)
               readBufferLineWise (STDERR_FILENO,
                                   &errBufferPtr,
                                   &errIdxPtr);
         }
     }
#  endif
}

void SxProcess::readBuffer (int FD)
{
   SX_CHECK (bufLen > 0, bufLen);
   SxList<SxPtr<SxArray<char>>> bufferList;
   SxList<ssize_t> bufferSizeList;
   // --- read from child's STDOUT/STDERR
   ssize_t nBytes   = 0;
#ifndef WIN32
   int fd = pipes[FD].fds[0];
   // make read non-blocking
   fcntl (fd, F_SETFL, fcntl (fd, F_GETFL,0) | O_NONBLOCK);
#endif

   do {
      SxPtr<SxArray<char> > buffer = SxPtr<SxArray<char> >::create (bufLen);
      size_t nRead = static_cast<size_t>(bufLen);
#     ifdef WIN32
         DWORD nr;
         DWORD availBytes;
         // Try to ReadFile without actually reading to see
         // if there is something to read
         BOOL status;
         status = ReadFile(pipes[FD].handles[0], buffer->elements,
                           0, &nr, &pipes[FD].oOverlap);
         if (!status) {
            DWORD err = ::GetLastError();
            // if IO is pending then break
            if (err == ERROR_IO_PENDING) { break; }
            if (err == ERROR_BROKEN_PIPE) // child process closed pipe; done
            {
               if (FD == STDERR_FILENO)  { eofStdErr = true; }
               else                      { eofStdOut = true; }
               break;
            }
            SX_THROW("PipeFailed", "Can't read from pipe: " + sxstrerror(err));
         }
         else {
            // Find out total number of bytes available to be read
            status = PeekNamedPipe(pipes[FD].handles[0], NULL, 0, 
                                   NULL, &availBytes, NULL);
            if (!status) {
               SX_THROW("PipeFailed", "Can't read from pipe: " + sxstrerror(GetLastError()));
            }
            status = ReadFile(pipes[FD].handles[0], buffer->elements,
                              availBytes, &nr, &pipes[FD].oOverlap);
            // --- possible (!status) is 'The pipe has been ended'
            if (status) {
               nBytes = static_cast<ssize_t>(nr);
            }
            else {
               DWORD err = ::GetLastError();
               if (err == ERROR_IO_PENDING) { break; }
               if (err == ERROR_BROKEN_PIPE)  // child process closed pipe; done
               {
                  if (FD == STDERR_FILENO) { eofStdErr = true; }
                  else                     { eofStdOut = true; }
                  break;
               }
               SX_THROW("PipeFailed", "Can't read from pipe: " + sxstrerror(err));
            }
         }
#     else /* WIN32 */
         SX_DBG_MSG ("read nBytes=" << nBytes << " bufLen=" << nRead);
         nBytes = ::read (fd, buffer->elements, nRead);
         if (nBytes == -1){
            // IO pending
            break;
         } else if (nBytes <= 0) {
            if (FD == STDERR_FILENO)  eofStdErr = true;
            else                      eofStdOut = true;
            break;
         }

#     endif /* WIN32 */
      bufferList << buffer;
      bufferSizeList << nBytes;
   } while (nBytes > 0);

   // --- total size
   ssize_t i = 0;
   ssize_t n = 0;
   SxArray<ssize_t> bufferSizeArr (bufferSizeList);
   for (i = 0; i < bufferSizeArr.getSize (); i++)  {
      n += bufferSizeArr(i);
   }

   // if total bytes read is zero then return
   if (n == 0) return;

   // --- concatenate all read buffers
   SxArray<char> buffer(n + 1);

   size_t offset = 0;
   size_t len = 0;
   i = 0;
   SxList<SxPtr<SxArray<char> > >::ConstIterator it;
   for (it = bufferList.begin (); it != bufferList.end(); ++it)  {
      len = (size_t)bufferSizeArr(i);
      ::memcpy (buffer.elements + offset, (*it)->elements, len);
      offset += len;
      i++;
   }
   buffer(n) = '\0';  // explicit NULL termination
   if (FD == STDERR_FILENO)  readFromStderr (buffer.elements);
   else                      readFromStdout (buffer.elements);
}

void SxProcess::read ()
{
   if (!(channels & (StdOut | StdErr)))  {
      return;
   }
   bool useStderr = (channels & StdErr) && (!combinedStdout ());
   if (lineWise) {
      readLineWise (useStderr);
      return;
   }

#  ifdef WIN32
      DWORD nCount = 0;
      if (channels & Stdout) {
         hEvents[nCount] = pipes[STDOUT_FILENO].oOverlap.hEvent;
         nCount++;
      }
      if (useStderr) {
         hEvents[nCount] = pipes[STDERR_FILENO].oOverlap.hEvent;
         nCount++;
      }
      while (1) {
         if (useStderr) {
            if ((channels & StdOut)) {
               if (eofStdOut && eofStdErr)
                  break;
            } else if (eofStdErr) {
               break;
            }
         } else {
            if ((channels & StdOut) && eofStdOut)  break;
         }
         // wait for signal from pipes
         DWORD dwWait = WaitForMultipleObjects(nCount,
                                               hEvents,
                                               FALSE,
                                               INFINITE);
         DWORD i = dwWait - WAIT_OBJECT_0;
         if (i > (nCount - 1))
            SX_THROW("PipeFailed", "Signaled index out of range");
         if ((channels & StdOut) && i == 0 && !eofStdOut)
            readBuffer(STDOUT_FILENO);
         if (useStderr && i == (nCount-1) && !eofStdErr)
            readBuffer(STDERR_FILENO);
      }
#  else
      struct timeval tv;
      fd_set set;
      int res = 0;
      while (1) {
         if (useStderr) {
            if ((channels & StdOut)) {
               if (eofStdOut && eofStdErr)
                  break;
            } else if (eofStdErr) {
               break;
            }
         } else {
            if ((channels & StdOut) && eofStdOut)  break;
         }
         tv.tv_sec = 1;
         tv.tv_usec = 0;
         FD_ZERO(&set);
         if (channels & StdOut)
            FD_SET(pipes[STDOUT_FILENO].fds[0], &set);
         if (useStderr)
            FD_SET(pipes[STDERR_FILENO].fds[0], &set);
         res = select (FD_SETSIZE, &set, NULL, NULL, &tv);
         if (res >= 1) {
            if ((channels & StdOut) && !eofStdOut)  readBuffer (STDOUT_FILENO);
            if (useStderr && !eofStdErr)            readBuffer (STDERR_FILENO);
         }
     }
#  endif
}

int SxProcess::wait ()
{
   // --- close the writing end of the child STDIN pipe to cause EOF
   closePipe (STDIN_FILENO, parentEnd);

#ifdef WIN32
   // --- wait for end of process
   childProc->wait ();
   //if (childProc->getError () != "") {
   //   // ...
   //}

   closePipes (); // pipes are no longer needed

   processExitStatus = childProc->processExitStatus;
   processExitSignal = childProc->processExitSignal;

   if (processExitStatus == static_cast<int>(ERROR_PROCESS_ABORTED))  {
      processExitSignal = 9; // SIGKILL
   }  else if (processExitStatus == static_cast<int>(ERROR_CONTROL_C_EXIT))  {
      processExitSignal = 15; // SIGTERM
   }

#else /* WIN32 */
   // --- wait for termination of child process
   int options = 0;
   if (noHang)  {
      options = WNOHANG;
   }
   pid_t ret = ::waitpid (pid, &processExitStatus, options);
   if (ret < 0)  {
      if (errno == ECHILD)  {
         // --- process pid does not exist
         processExitSignal = SIGKILL;
      }
      // errno EINTR, EINVAL
   } else if (WIFEXITED(processExitStatus))  {
      processExitStatus = WEXITSTATUS(processExitStatus);
   } else if (WIFSIGNALED(processExitStatus))  {
      // --- get and translate signal
      processExitSignal = WTERMSIG(processExitStatus);
      SxString name;
      switch (processExitSignal)  {
         case SIGKILL : name="SIGKILL"; break;
         case SIGFPE  : name="SIGFPE"; break;
         case SIGSEGV : name="SIGSEGV"; break;
         case SIGILL  : name="SIGILL"; break;
         default      : name=SxString(processExitSignal);
      }

      SxString msg = "Received signal " + name;

//#ifdef WCOREDUMP
//      if (WCOREDUMP())  {
//         msg += " (core dumped)";
//      }
//#endif /* WCOREDUMP */
      msg += ".\n";
      // write msg to error channel
      if (channels & Stderr)  readFromStderr (msg);
   }

   // --- finish process timeout
   if (processTimer.getPtr ())  {
      processTimer->stop ();
      processTimer->wait ();
      processTimer = SxPtr<SxProcessTimer>();
   }
#endif /* WIN32 */

   if (processExitStatus == 128)  {
      SX_THROW ("Can't execute program '" + cmd
                + "': Exit status is 128");
   }

   finishExecution ();

   processRunning = false;

   return processExitStatus;
}

int SxProcess::exitStatus () const
{
   if (processRunning)  return 0;
   return processExitStatus;
}
 

int SxProcess::exitSignal () const
{
   if (processRunning)  return 0;
   return processExitSignal;
}

#ifdef WIN32
// List of System Error Codes (exitCode)
// msdn.microsoft.com/en-us/library/windows/desktop/ms681381(v=vs.85).aspx
void SxProcess::winkill (pid_t pid_, UINT exitCode_, bool wait_)
{
   DWORD access = PROCESS_TERMINATE;
   HANDLE hProcess = OpenProcess (access, FALSE, pid_);
   if (hProcess)  {
      BOOL status = TerminateProcess (hProcess, exitCode_);
      if (!status)  {
         DWORD err = ::GetLastError ();
         SX_THROW ("TerminateProcess() pid " + SxString(pid_)
                   + " failed: " + sxstrerror (err));
      }
      if (wait_)  {
         DWORD result = WaitForSingleObject (hProcess, INFINITE);
         if (result != WAIT_OBJECT_0)   {
            sxprintf ("Wait for single object unexpectedly quit: %x\n", result);
            if (result == WAIT_FAILED)  {
               DWORD err = GetLastError ();
               SX_THROW ("Wait for pid " + SxString(pid_)
                         + " to terminate failed: " + sxstrerror (err));
            }
         }
      }
   } else  {
      DWORD err = GetLastError ();
      SX_THROW ("Failed to obtain Process handle. Winkill for pid " + SxString(pid_)
         + " failed: " + sxstrerror (err));
   }
}
#endif /* WIN32 */

int SxProcess::kill (pid_t pid_, int signalId_, bool tree_, bool wait_)
{
   int exitStatus = 0;

   if (!tree_)  {
#     ifdef WIN32
         SX_UNUSED (signalId_);
         if (pid_ < 0)  {
            pid_ = -pid_;
         }
         // exitCode ERROR_CONTROL_C_EXIT as SIGTERM
         winkill (pid_, ERROR_PROCESS_ABORTED, wait_);
#     else
         ::kill (pid_, signalId_);
         if (wait_)  {
            ::waitpid (pid_, &exitStatus, 0);
            if (WIFEXITED(exitStatus))  {
               exitStatus = WEXITSTATUS(exitStatus);
            }
         }
#     endif
   }  else  {
      SxList<pid_t> killList;

#     ifdef WIN32
         killList << pid_;
#     else
         // --- if both /proc/*/status and command ps are not available
         //     kill tree will send signal only the root process pid_

         SxList<pid_t> stack;
         stack << pid_;

         while (stack.getSize () > 0)  {
            pid_t p = stack.last ();
            stack.removeLast ();
            killList.prepend (p);

            if (signalId_ == SIGKILL)  {
               // --- do not spawn new processes in p
               //     other signalId_ might leave some procs running
               kill (p, SIGSTOP, false);
            }

            // --- get pid of all children of p
            SxList<SxString> lines;

            // --- 1) try /proc/<pid>/status to read ppid
            SxString path = "/proc";
            DIR *dir = opendir (path.getElems ());
            if (dir)  {
               struct dirent *entry = NULL;
               while ((entry = readdir (dir)) != NULL)  {
                  if (entry->d_type == 4)  {
                     SxString name(entry->d_name);
                     if (name.isInt ())  {
                        SxString filename = path + "/" + name + "/stat";
                        FILE *f = fopen (filename.getElems (), "r");
                        if (f) {
                           int ppid = 0;
                           fscanf (f, "%*d %*s %*s %d ", &ppid);
                           if (ppid == static_cast<int>(p))  {
                              lines << (name + " " + SxString(ppid));
                           }
                           fclose (f);
                        }
                     }
                  }
               }
               closedir (dir);
            }  else  {
               // --- 2) fallback to try ps command
               try {
                  SxProcess psProc("ps -e -o pid= -o ppid=");
                  psProc.enableBuffer (SxProcess::StdOut);
                  if (psProc.start () == 0)  {
                     SxString psOutput = psProc.getBuffer (SxProcess::StdOut);
                     lines = psOutput.tokenize ('\n');
                  }
               }  catch (SxException e)  {
                  // Can't execute program ps
               }
            }

            // --- put children of p on the stack
            SxList<SxString>::Iterator it;
            for (it = lines.begin(); it != lines.end(); ++it)  {
               SxArray<SxString> tokens = it->tokenize (' ');
               if (   tokens.getSize () == 2
                   && tokens(0).isInt ()
                   && tokens(1).isInt ())
               {
                  pid_t tpid  = tokens(0).toInt();
                  pid_t tppid = tokens(1).toInt();
                  if (tppid == p)  {
                     stack << tpid;
                  }
               }
            }
         }
#     endif /* WIN32 */

      SxList<pid_t>::ConstIterator it = killList.begin ();
      while (it.isValid ())  {
         exitStatus = kill (*it, signalId_, false, wait_);
         ++it;
      }
   }

   return exitStatus;
}

void SxProcess::softTerminate (int signalId_)
{
   if (processRunning)  {
      processExitStatus = kill (pid, signalId_, killTree, true);
      SX_MUTEX (statusMutex)  {
         processRunning = false;
      }
   }
}

void SxProcess::kill ()
{
   if (!processRunning)  return; // nothing to do

#ifdef WIN32
   if (killTree && IsWindows8OrGreater ())  {
      childProc->removeJobObject ();
      // (Can't get a processExitStatus here because the job processes are
      // killed asynchronously.)
      return; // done
   } else {
      childProc->removeProcess ();
   }
   int SIGKILL = -1; // not used
#else
   processExitStatus = kill (pid, SIGKILL, killTree, true);
#endif
   SX_MUTEX(statusMutex) {
      processRunning = false;
   }
}

void SxProcess::finishExecution ()
{
   // empty
}

void SxProcess::closePipeIdx (int stdStreamFd, uint8_t idx)
{
   SX_CHECK (idx == 0 || idx == 1, idx);
#ifdef WIN32
   SX_CHECK (stdStreamFd == STDIN_FILENO || stdStreamFd == STDOUT_FILENO
             || stdStreamFd == STDERR_FILENO, stdStreamFd);
   if (pipes[stdStreamFd].handles[idx] != INVALID_HANDLE_VALUE)  {
      ::CloseHandle (pipes[stdStreamFd].handles[idx]);
      pipes[stdStreamFd].handles[idx] = INVALID_HANDLE_VALUE;
   }
#else
   SX_CHECK (stdStreamFd >= STDIN_FILENO && stdStreamFd <= STDERR_FILENO,
             stdStreamFd);
   if (pipes[stdStreamFd].fds[idx] >= 0)  {
      ::close (pipes[stdStreamFd].fds[idx]);
      pipes[stdStreamFd].fds[idx] = -1;
   }
#endif
}

void SxProcess::closePipe (int stdStreamFd, PipeEnds ends)
{
   if (ends & childEnd)
      closePipeIdx (stdStreamFd, (stdStreamFd == STDIN_FILENO) ? 0 : 1);
   if (ends & parentEnd)
      closePipeIdx (stdStreamFd, (stdStreamFd == STDIN_FILENO) ? 1 : 0);
}

void SxProcess::closePipes ()
{
#ifdef WIN32
   for (uint8_t i = STDIN_FILENO; i <= STDERR_FILENO; i++)
      closePipe (i, bothEnds);
#else
   for (uint8_t i = STDIN_FILENO; i <= STDERR_FILENO; i++)
      closePipe (i, bothEnds);
#endif
}

void SxProcess::closeStdStream (int stdStreamFd)
{
   closePipe (stdStreamFd, parentEnd);
}

#ifdef WIN32
void SxProcess::setHiddenWindow (bool value)
{
   SX_TRACE ();
   SX_CHECK (!isRunning ());
   SX_CHECK (!childProc.getPtr ());
   hiddenWindow = value;
}
#endif

#ifndef WIN32
int SxProcess::getPipeFd (int stdStreamFd) const
{
   SX_CHECK (stdStreamFd >= STDIN_FILENO && stdStreamFd <= STDERR_FILENO,
             stdStreamFd);
   int res = pipes[stdStreamFd].fds[(stdStreamFd == STDIN_FILENO) ? 1 : 0];
   SX_CHECK (res >= 0, res);
   return res;
}
#endif

SxString SxProcess::getBuffer (Channels c)
{
   switch (c)  {
      default     :  SX_EXIT; break;
      case StdOut :  SX_CHECK (useBuffers & StdOut);
                     return outBuffer;
      case StdErr :  SX_CHECK (useBuffers & StdErr);
                     return errBuffer;
   }
   return "";
}


SxString SxProcess::writeToStdin ()
{
   // --- in order to communicate with child's STDIN, the derived class
   //     has to overload writeToStdin()
   //     Otherwise do not activate the STDIN inter-process communication.
   SX_EXIT;
   return SxString();
}


void SxProcess::setLineWise (bool enable)
{
   lineWise = enable;
}


void SxProcess::readFromStdout (const SxString &data)
{
   // --- in order to communicate with child's STDOUT, the derived class
   //     has to overload readFromStdout ()
   //     Otherwise do not activate the STDOUT inter-process communication.
   SX_CHECK (useBuffers != None);
   outBuffer += data + "\n";
}


void SxProcess::readFromStderr (const SxString &data)
{
   // --- in order to communicatw with child's STDERR, the derived class
   //     has to overload readFromStderr ()
   //     Otherwise do not activate the STDERR inter-process communication.
   SX_CHECK (useBuffers != None);
   errBuffer += data + "\n";
}

