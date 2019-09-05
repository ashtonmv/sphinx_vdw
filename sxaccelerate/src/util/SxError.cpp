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

#include <SxError.h>
#include <SxUtil.h>
#include <SxTime.h>
#include <SxRedirect.h>
#include <stdio.h>
#include <iostream>
#include <errno.h>
#include <SxConfig.h>

#include <stdarg.h>

#ifndef WIN32

#  include <sys/types.h>
#  include <sys/wait.h>
#  include <unistd.h>


#endif

#ifdef SX_ANDROID
#  include <android/log.h>
#endif /* SX_ANDROID */


// --------------------------------------------------------------------------
void SxDebug::log (int level_, const char *tag_, const char *func_,
                   const char *fmt_, ...)
{
   if (fmt_)  {
#     ifdef SX_ANDROID
         int level = 0;
         switch (level_) {
            case 2:
               level = ANDROID_LOG_VERBOSE;
               break;
            case 3:
               level = ANDROID_LOG_DEBUG;
               break;
            case 4:
               level = ANDROID_LOG_INFO;
               break;
            case 5:
               level = ANDROID_LOG_WARN;
               break;
            case 6:
               level = ANDROID_LOG_ERROR;
               break;
         }
#     else
         SX_UNUSED (level_);
#     endif /* SX_ANDROID */

      va_list arg;
      va_start (arg, fmt_);
#     ifdef SX_ANDROID
         __android_log_vprint (level, tag_, fmt_, arg);
#     else
         sxprintf ("___%s: %s: ", tag_, func_);
         vfprintf (stderr, fmt_, arg);
         sxprintf ("\n");
#     endif /* SX_ANDROID */
      va_end (arg);
   }
}


int SxDebug::isEnabled (const char *file, const char *func, const char *feat)
{
   int res = 0;
   if (getenv ("SX_DEBUG_SIMPLE"))  {
      res = 1;
   }
   
   if (getenv ("SX_DEBUG"))  {
      res = 2;
   }
   
   if (file)  {
      char *filename = getenv ("SX_DEBUG_FILE");
      if (filename &&
          strlen (filename) > 0 && 

#  ifndef WIN32
                 strstr (filename, file))
#  else
                 _stricmp (filename, file) == 0)
#  endif
      {

         // --- missmatch
         //     filename_ = "SxArray.h"
         //     filename  = "SxError.h"
         return res > 0 ? res : 2;
      }
   }

   if (func)  {
      char *strFunctions = getenv ("SX_DEBUG_FUNC");
      if (strFunctions
         &&  strlen (strFunctions) > 0
         &&  strstr (strFunctions, func))
      {
         return res > 0 ? res : 2;
      }
   }

   if (feat)  {
      char *strFeature = getenv ("SX_DEBUG_FEAT");
      if (  strFeature
         && strlen (strFeature) > 0
         && strstr (strFeature, feat))
      {
         return res > 0 ? res : 2;
      }
   }
   return res;
}

int SxDebug::isEnabled (const char *file, const char *func, const char *sxLogId,
                        const char *feat)
{
   int res = isEnabled (file, func, feat);
   if (res > 0)  return res;

   if (sxLogId) {
      char *strComponents = getenv ("SX_DEBUG_COMP");
      if (strComponents &&
          strlen (strComponents) > 0 &&
          strstr (strComponents, sxLogId))
      {
         return res > 0 ? res : 2;
      }
   }

   return res;
}


std::ostream &SxDebug::tag (std::ostream &s,
                            const char *file, long line, const char *func,
                            int tron)
{
   if      (tron == 1)  s << "___";
   else if (tron == 2)  s << "___" <<file<< ":" <<line<< ": " <<func<< ": ";
   return s;
}
#  ifdef WIN32

std::wostream &SxDebug::tag(std::wostream &s,
	const char *file, long line, const char *func,
	int tron)
{
	if (tron == 1)
      s << "___";
	else if (tron == 2)
      s << "___" << file << ":" << line
        << ": " << func << ": ";
	return s;
}
#  endif

void SxDebug::stackDump (int skipLvls)
{
#  ifndef WIN32

      char cmdBuf[10240];
#     if defined(LINUX)
         char nameBuf[512];
         char pidBuf[30];
         sprintf (pidBuf, "--pid=%d", getpid());
         nameBuf[readlink("/proc/self/exe", nameBuf, 511)]=0;
         sprintf (cmdBuf, "gdb --batch -n -ex thread -ex bt %s %s 2>&1"
                          "| grep '^#' | tail -n+%d"
                          "| sed 's/^#\\([0-9]\\+\\)//g;"
                                 "s/ in / /g;"
                                 "s/ at / /g;"
                                 "s/\\s\\+0x[0-9A-Fa-f]\\+//g;"
                                 "s/main (/main(/g;"
                                 "s/, /,/g'"
                          "| awk '{print NR\":\",$NF\":\",$1}'",
                    nameBuf, pidBuf, skipLvls+1);

#     elif defined(MACOSX)
         sprintf (cmdBuf, "lldb attach -p %d "
                               "--batch -o 'bt' -o 'detach' -o 'quit' 2>&1"
                          "| grep frame "
                          "| sed 's/*/ /g;"
                                 "s/^\\(.*\\)`\\(.*\\)$/\\2/g;"
                                 "s/ at /:/g;s/ \\[[^]]*\\]//g'"
                          "| tail -n+%d",
                  getpid(), skipLvls+1);

#     elif defined(SX_MOBILE)
         SX_EXIT;
#     else
#        error "Unsupported OS"
#     endif

      int childPid = fork();
      if (!childPid) {           
          printf ("\nCall stack:\n");
          execlp ("/bin/sh", "sh", "-c", cmdBuf, NULL);
          printf ("\n\n");
          return;
      } else {
          waitpid (childPid,NULL,0);
      }
#  endif


#  ifdef WIN32
#     define MAX_STACK_SIZE 1024
#     define MAX_FNNAME_SIZE 1024

      void *stack[MAX_STACK_SIZE];
      HANDLE proc = GetCurrentProcess ();
      SymInitialize (proc, NULL, TRUE);
      WORD nFrames = CaptureStackBackTrace (skipLvls, MAX_STACK_SIZE, stack, NULL);
      SYMBOL_INFO *sym = (SYMBOL_INFO *)malloc (sizeof (SYMBOL_INFO) + (MAX_FNNAME_SIZE - 1) * sizeof(TCHAR));
      sym->MaxNameLen = MAX_FNNAME_SIZE;
      sym->SizeOfStruct = sizeof (SYMBOL_INFO);
      DWORD offset = 0;
      IMAGEHLP_LINE64 *line = (IMAGEHLP_LINE64 *)malloc (sizeof(IMAGEHLP_LINE64));
      line->SizeOfStruct = sizeof (IMAGEHLP_LINE64);
      if (nFrames > 0)  printf ("\nCall stack:\n");
      for (int i = 0; i < nFrames; i++)  {
         DWORD64 addr = (DWORD64)(stack[i]);
         SymFromAddr (proc, addr, NULL, sym);
         if (SymGetLineFromAddr64 (proc, addr, &offset, line))  {
            printf ("%d: %s:%lu: %s (0x%llX)\n", i, line->FileName, line->LineNumber, sym->Name, sym->Address);
         } else {
//          printf ("%d: SymGetLineFromAddr64 failed (error %lu)\n", iLvl, GetLastError());
            printf ("%d: ---:--- %s (0x%llX).\n", i, sym->Name, sym->Address);
        }
    }
#   endif
}



bool SxDebug::dump ()
{
   int skipLvls = 2;  // skip SxDebug::dump,stackDump
   SxDebug::stackDump (skipLvls);

#ifdef WIN32
   // --- create helper thread in suspended mode
   DWORD threadId = 0;
   HANDLE hThread = CreateThread (NULL, 0, writeCB, NULL, CREATE_SUSPENDED,
                                  &threadId);

   if (hThread) {
      // --- suspend all other threads to create memory snapshot
      enumerateThreads (SuspendThread, threadId);

      // --- resume helper thread and write dump
      ResumeThread (hThread) ;
      WaitForSingleObject (hThread, INFINITE);

      DWORD dumpRes = 0;
      GetExitCodeThread (hThread, &dumpRes);
      CloseHandle (hThread);

      // --- resume all other threads
      enumerateThreads (ResumeThread, threadId);

      printf ("Segmentation fault. Core dumped. (core.dmp)\n");

      return dumpRes == 0;
   }
#  endif /* WIN32 */
   return true;
}

#ifdef WIN32
BOOL CALLBACK SxDebug::filterCB (PVOID inParam,
                                 const PMINIDUMP_CALLBACK_INPUT cbInput,
                                 PMINIDUMP_CALLBACK_OUTPUT)
{
   if (cbInput->CallbackType == CancelCallback)  return FALSE;

   switch (cbInput->CallbackType) {
      case ModuleCallback:        return TRUE;
      case ThreadCallback:        return TRUE;
      case ThreadExCallback:      return TRUE;
      case IncludeThreadCallback: {
         // --- exclude minidump's helper thread from dump
         if (cbInput->IncludeThread.ThreadId == ::GetCurrentThreadId() )       
            return FALSE;
         else
            return TRUE;
      }
      case IncludeModuleCallback: return TRUE;
      case MemoryCallback:        return TRUE;
   }

   return FALSE;
}

DWORD CALLBACK SxDebug::writeCB (LPVOID)
{
   HANDLE fp = CreateFileA ("core.dmp", GENERIC_WRITE, 0, NULL,
                            CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
   if (fp == INVALID_HANDLE_VALUE)  return 0;
   int type = MiniDumpNormal
            | MiniDumpWithDataSegs
            | MiniDumpWithPrivateReadWriteMemory
            | MiniDumpWithHandleData
            | MiniDumpWithFullMemory
            | MiniDumpWithFullMemoryInfo
            | MiniDumpWithThreadInfo
            | MiniDumpWithUnloadedModules
            | MiniDumpWithProcessThreadData;

   MINIDUMP_CALLBACK_INFORMATION callback = { 0 };
   callback.CallbackParam = NULL;
   callback.CallbackRoutine = filterCB;

   BOOL bRet = MiniDumpWriteDump (GetCurrentProcess(), GetCurrentProcessId(),
                                  fp, (MINIDUMP_TYPE)type, NULL, NULL,
                                  &callback);
   CloseHandle (fp);

   return bRet;
}

void SxDebug::enumerateThreads (DWORD (WINAPI *cb)(HANDLE),
                                DWORD skipThreadId)
{
   HANDLE hSnapshot = CreateToolhelp32Snapshot (TH32CS_SNAPTHREAD, 0);
   if (hSnapshot != INVALID_HANDLE_VALUE)  {
      THREADENTRY32 thread;
      thread.dwSize = sizeof (thread);
      if (Thread32First (hSnapshot, &thread))  {
         do {
            if (  thread.th32OwnerProcessID == GetCurrentProcessId()
               && thread.th32ThreadID != skipThreadId 
               && thread.th32ThreadID != GetCurrentThreadId())
            {
               HANDLE hThread = OpenThread (THREAD_SUSPEND_RESUME, FALSE,
                                            thread.th32ThreadID);
               if (hThread) {
                  cb (hThread);
                  CloseHandle (hThread);
               }
            }
         } while (Thread32Next (hSnapshot, &thread));
      }

      CloseHandle( hSnapshot );
   }
}
#endif /* WIN32 */



// --------------------------------------------------------------------------


void sxExit (int exitCode)
{
#  ifndef SX_ANDROID
      fflush (stdout);
      fflush (stderr);
      std::cout.flush ();
      std::cerr.flush ();
#  endif

#  ifdef NDEBUG
      exit (exitCode);
#  else
      SX_UNUSED (exitCode);
      SxDebug::dump ();
      abort ();
#  endif

#  ifdef SX_IOS
      __builtin_unreachable();
#   endif /* SX_IOS */
}

void sxQuit ()
{
   exit (1);

#  ifdef SX_IOS
      __builtin_unreachable();
#  endif /* SX_IOS */
}

void sxOutOfMemoryHandler ()
{
   sxprintf ("\nSxAccelerate: memory allocation failed.\n");
   sxExit (2);
}

void sxUncaughtExceptionHandler ()
{
   sxprintf ("\nSxAccelerate: unexpected exception\n");
   sxExit (3);
}


SxString sxstrerror ()
{
#  ifdef WIN32
      return sxstrerror (GetLastError ());
#  else
      return sxstrerror (errno);
#  endif /* WIN32 */
}

SxString sxstrerror (int errnum)
{
   SxString msg;
#  ifdef SX_ANDROID
      SX_EXIT; // not yet implemented
#  else
#     ifdef WIN32
         char buf[1024];
         if (::strerror_s (buf, sizeof (buf), errnum) != 0)  msg = "Unknown error";
         else  msg = SxString(buf);
#     else
#        ifdef HAVE_STRERROR_R
            char buf[1024]; // message is truncated if buf is too small
#           if defined(STRERROR_R_CHAR_P) || defined(__USE_GNU)
               char const *str = strerror_r (errnum, buf, sizeof (buf));
               if (str)  msg = str;
               else      msg = "Unknown error";
#           else
               strerror_r (errnum, buf, sizeof (buf));
               msg = buf;
#           endif /* STRERROR_R_CHAR_P */
#        else
            msg = strerror (errnum);
#        endif /* HAS_STRERROR_R */
#     endif /* WIN32 */
#  endif /* SX_ANDROID */

   return msg;
}

#ifdef WIN32
SxString sxstrerror (DWORD sysErrCode)
{
   LPVOID lpMsgBuf;
   DWORD ret = FormatMessage (FORMAT_MESSAGE_ALLOCATE_BUFFER |
                              FORMAT_MESSAGE_FROM_SYSTEM |
                              FORMAT_MESSAGE_IGNORE_INSERTS,
                              NULL,
                              sysErrCode,
                              MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                              (LPTSTR) &lpMsgBuf,
                              0, NULL);
   if (ret)  {
      SxString msg = SxString((LPCSTR)lpMsgBuf, static_cast<int>(ret));
      LocalFree (lpMsgBuf);
      return msg;
   }  else  {
      return "FormatMessage () error";
   }
}
#endif /* WIN32 */

#ifdef SX_IOS
SxString sxstrerror (CFErrorRef error)
{
   SxString result;
   if (error != NULL) {
      CFIndex errorCode = CFErrorGetCode(error);
      if (errorCode != 0) {
         result = "error " + SxString(errorCode);
         CFStringRef desc = CFErrorCopyDescription(error);
         if (desc) {
            result += ": " + SxString(desc);
            CFRelease(desc);
         }
      }
   }
   return result;
}
#endif /* SX_IOS */

// --------------------------------------------------------------------------

// This construct installs the above handlers for
// - new_handler: new failed to allocate memory
// - unexpected : an exception was not caught
// before the main program starts
static class SxInstallHandlers {
   public:
      SxInstallHandlers ()  {
         std::set_new_handler (sxOutOfMemoryHandler);
         std::set_unexpected (sxUncaughtExceptionHandler);
      }
} theHandlerInstaller;

SxTrace::SxTrace (int enabled_, const char *filename_, long line_,
                  const char *function_, const char *functionShort_,
                  const char *dbgMsg_)
   : enabled (enabled_),
     line(line_)
{
   if (enabled > 1)  {
      strncpy (filename, filename_, 79);
      filename[79] = '\0';
      strncpy (functionName, function_, 79);
      functionName[79] = '\0';
      strncpy (functionNameShort, functionShort_, 79);
      functionNameShort[79] = '\0';

      dbgMsg[0] = '\0';
      if (dbgMsg_)  {
         strncpy (dbgMsg, dbgMsg_, 79);
         dbgMsg[79] = '\0';
      }


      SxUtil::getGlobalObj().lockLog ();
      std::cerr << "___" << filename <<":" << line << " "
                         << functionNameShort <<": BEGIN";
      if (dbgMsg[0]) std::cerr << " [" << dbgMsg << "]";
      std::cerr << std::endl;
      std::cerr.flush ();
      SxUtil::getGlobalObj().unlockLog ();

      // --- wall clock time 1
      t1 = SxTime::getRealTime ();
   }
}
   
   
SxTrace::~SxTrace ()
{
   if (enabled > 1)  {
      // --- wall clock time 2
      t2 = SxTime::getRealTime ();

      // --- Print time spend in procedure in milliseconds.
      //     general performance is in digit place ([xy]ms vs [z]ms)
      SxUtil::getGlobalObj().lockLog ();
      std::cerr << "___"
                << filename << ":" << line << " "
                << functionNameShort << ": END ";
      if (dbgMsg[0]) std::cerr << " [" << dbgMsg << "] ";
      std::cerr << std::fixed << std::setprecision(3)
                << ((t2 - t1) * 1.0e3) << " ms" << std::endl;

      std::cerr.flush ();
      SxUtil::getGlobalObj().unlockLog ();
   }
}

