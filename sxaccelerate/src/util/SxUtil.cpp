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

#include <SxUtil.h>
#include <SxCrashHandler.h>
#include <SxString.h>
#include <SxException.h>
#include <iostream>
#include <string.h>   /* strstr */
#ifndef WIN32
#  include <unistd.h>   /* gethostname() */
#  include <sys/stat.h>
#endif
#ifdef WIN32
#  include <Windows.h>
#endif
//#include <netdb.h>    /* gethostbyname() */ 
#ifdef MACOSX
#  include <mach-o/dyld.h>	/* _NSGetExecutablePath */
#endif


//#ifdef USE_ESSL
////#  define ESVERR
//#  include <essl.h>
//   extern "Fortran" {
//      int enotrm (int &, int &);
//      typedef int (*FN) (int &, int &);
//   }
//#endif /* USE_ESSL */   

#ifdef USE_OPENMP
#  include <omp.h>
#  include <pthread.h>
#endif /* USE_OPENMP */

char __flags_ = 0x0;
long sxChunkSize = 5000L;

SxUtil &SxUtil::getGlobalObj ()
{
   static SxUtil sxUtil;
   return sxUtil;
}

SxUtil::SxUtil ()
{
#  ifndef SX_IOS
      // sync should be by default?
      std::cout.sync_with_stdio();
      std::cerr.sync_with_stdio();
#  endif /* SX_IOS */

   __flags_ = 0x0;
   nProcs    = 1;
   sxChunkSize = 5000L;

   // --- initializing crash handler
#  ifndef SX_IOS
      SxCrashHandler::getGlobalInstance ();
      atexit (SxCrashHandler::trapCrash);
#  endif /* SX_IOS */
   
//#  ifdef USE_ESSL
//      FN iusadr;
//      int ierno,inoal,inomes,itrace,irange,irc,dummy;
//      char storarea[8];
//#  endif /* USE_ESSL */ 

// char name[10240];

//#  ifdef USE_ESSL
//      iusadr = enotrm;
//      // --- init ESSL error option table
//      dummy  = 0;
//      einfo (0, dummy, dummy);
//      // --- make error code 2015 a recoverable error and
//      //     suppress printing all error messages for it
//      ierno = 2015;
//      inoal = 0;
//      inomes = -1;
//      itrace = 0;
//      irange = 2015;
//      errset (ierno,inoal,inomes,itrace,iusadr,irange);
//
//      errsav (ierno, storarea);
//#  endif /* USE_ESSL */

#  ifdef USE_OPENMP
      char *str = getenv ("SX_THREADS");
      if (str)  nProcs = atoi (str);
      if (nProcs == 0)    nProcs = getNProcs ();
      if (nProcs < 1)     nProcs = 1;
      if (nProcs > 1000)  nProcs = 1000;
      omp_set_num_threads (nProcs);

      str = getenv ("SX_CHUNK_SIZE");
      if (str)  sxChunkSize = atol (str);
      if (sxChunkSize < 1)  sxChunkSize = 0L;
#  endif /* USE_OPENMP */


#  ifdef HPUX
      // --- HPUX needs to have a standard pthread stack size of 4MB
      //     see also: /opt/mlib/newconfig/RelNotes/B6061-96016.pdf
      size_t newSize = 4*1024*1024;  // 4 MB
      size_t oldSize = 0;
      pthread_default_stacksize_np (newSize, &oldSize);      
#  endif  /* HPUX */

#  ifdef USE_SX_LOG
#     ifdef WIN32
         InitializeCriticalSection (&logMutex);
#     else
         pthread_mutexattr_init(&logMutexAttr);
         pthread_mutexattr_settype(&logMutexAttr, PTHREAD_MUTEX_RECURSIVE);
         pthread_mutex_init (&logMutex, &logMutexAttr);
#     endif
#  endif
}


SxUtil::~SxUtil ()
{
#  ifdef USE_SX_LOG
#     ifdef WIN32
         DeleteCriticalSection (&logMutex);
#     else
         pthread_mutex_destroy (&logMutex);
         pthread_mutexattr_destroy(&logMutexAttr);
#     endif
#  endif
}

#ifdef USE_SX_LOG
void SxUtil::lockLog ()
{
#  ifdef WIN32
      EnterCriticalSection (&logMutex);
#  else
      pthread_mutex_lock (&logMutex);
#  endif
}

void SxUtil::unlockLog ()
{
#  ifdef WIN32
      LeaveCriticalSection (&logMutex);
#  else
      pthread_mutex_unlock (&logMutex);
#  endif
}
#endif


void setNProcs (int n)
{
#  ifndef SX_IOS
      SxUtil::getGlobalObj().nProcs = n;
#  endif /* SX_IOS */

#  ifdef USE_OPENMP
      omp_set_num_threads (n);
#  endif /* USE_OPENMP */
}

int getNProcs ()
{
#  ifdef USE_OPENMP
      return omp_get_num_procs ();
#  else  /* USE_OPENMP */
      return 1;
#  endif /* USE_OPENMP */
}

void sxGetExecPath (char *resPath, size_t len)
{
   if (!resPath)  abort ();

#  if defined(LINUX)
      pid_t pid = getpid ();
      SxString link = "/proc/" + SxString(pid) + "/exe";
      ssize_t idx = readlink (link.ascii(), resPath, len); // returns int
      if (idx < 0)  {
         SX_THROW("Can't access '"+link+"' file system: readlink() failed: "
                  + sxstrerror());
         resPath[0] = '\0';
      }  else if ((size_t)idx >= len)  {
         SX_THROW ("Can't access '"+link+"' file system: readlink() "
                   "buffer is too small");
      }  else  {
         resPath[idx] = '\0';
      }
#  elif defined(MACOSX)
      uint32_t size = static_cast<uint32_t>(len);
      _NSGetExecutablePath (resPath, &size);
#  elif defined(FREEBSD)
      size_t size = len;
      int mib[4];
      mib[0] = CTL_KERN;
      mib[1] = KERN_PROC;
      mib[2] = KERN_PROC_PATHNAME;
      mib[3] = -1;
      sysctl (mib, 4, resPath, &size, NULL, 0);
#  elif defined(WIN32)
      DWORD size = (DWORD)len;
      GetModuleFileName (NULL, resPath, size);
#  elif defined(SX_ANDROID)
      // FIXME: implement something like https://github.com/gpakosz/whereami 
      SX_THROW ("Can't access path to executable");
#  elif defined(SX_IOS)
      // FIXME: implement something like https://github.com/gpakosz/whereami 
      SX_THROW ("Can't access path to executable");
#  else
#     error "Unsupported platform"
#  endif
}

void sxchmod (const char *pathname_, int mode_)
{
   SxString path(pathname_);
   SxString error;
#  ifdef WIN32
      typedef int mode_t;
#  endif

   if (mode_ < 0 || mode_ > UINT16_MAX) {
      error = "Invalid mode value " + SxString(mode_) + ". "
              "Enter a file mode that fits to octal mask 0177777.";
   }  else  {
      mode_t mode = static_cast<mode_t>(mode_);
      if (chmod (path.ascii(), mode) != 0) {
         error = "chmod() failed: " + sxstrerror ();
      }
   }

   if (error != "")  {
      SX_THROW("Can't change permissions of '"
               + path + "' to mode 0"
               + SxString::sprintf("%06o", mode_)
               + ". " + error);
   }
}
