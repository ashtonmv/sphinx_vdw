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

#ifndef _SX_LOG_H_
#define _SX_LOG_H_

#include <SxError.h>
#include <SxConfig.h>
#include <SxUtil.h>
#include <SxUniqueList.h>

#ifndef SX_LOG_ID
#  error "Compiler was invoked without -DSX_LOG_ID"
#endif

#ifndef SX_LOG_HASH
#  error "Compiler was invoked without -DSX_LOG_HASH"
#endif

#ifdef WIN32
#  define _WINSOCKAPI_ // prevent winsock.h to be included from windows.h
                       // otherwise later include of winsock2.h would show
                       // redefinition errors
#  include <windows.h> // DWORD, GetLastError()
#else
#   include <pthread.h>
#endif /* WIN32 */

#ifdef SX_ANDROID
#  include <android/log.h>
#endif


/** \brief Helper class for SX_LOG

    This class provides management function to provide basic functionality
    used in the SX_LOG macro.

    \author Sixten Boeck */
class SX_EXPORT_UTIL SxLog
{
   public:

     ~SxLog ();

      static void enable  (const char *soName);
      static void disable (const char *soName);
      static bool isEnabled (const uint32_t &);

      static uint32_t djb2 (const char *str);

      static void lock ();
      static void unlock ();

   protected:

      //SxMap<uint32_t,SxString> compIDs;  // hashes of shared object names
      SxUniqueList<uint32_t> compIDs;  // hashes of shared object names

      SxLog ();
      static SxLog &getGlobalObj ();
};



// --------------------------------------------------------------------------
// in windows: in order to display chinese characters, set console font to: 'SimSun-ExtB'

#if defined(SX_ANDROID)
#  define SX_LOG_MSG(TAG,MSG)   __android_log_print (ANDROID_LOG_INFO, TAG,\
                                                 "%s\n", (MSG).ascii())
#elif defined(SX_IOS)
void sxLogMsg (const char *, const char *); // defined in SxLog-ios.mm
#  define SX_LOG_MSG(TAG,MSG)   sxLogMsg (TAG, (MSG + "\n").ascii())
#elif defined(LINUX)
#  define SX_LOG_MSG(TAG,MSG)   std::cout << "[" << TAG << "] " << MSG << std::endl
#elif defined(WIN32)
#  define SX_LOG_MSG(TAG,MSG)   std::cout << "[" << TAG << "] " << MSG << std::endl
#elif defined(MACOSX)
#  define SX_LOG_MSG(TAG,MSG)   std::cout << "[" << TAG << "] " << MSG << std::endl
#else
#  error "Undefined platform"
#endif

// --------------------------------------------------------------------------
// SX_LOG(msg)
// SX_LOG(msg << "abc" << "def" << i);
// iterator.SX_LOG(it, "val=" << *it).foo ();
// --------------------------------------------------------------------------

// Note: sxLog() is also defined in SxIterator.h
template<class Fn>
void sxLog (const char *logId, uint32_t logHash,
            const char *file, long line, const char *func, Fn f)
{ 
#  ifdef USE_SX_LOG
      static const int dummy = 0;
      f (logId, logHash, file, line, func, dummy);
#  else
      SX_UNUSED (logId, logHash, file, line, func, f);
#  endif
}

// --- compare component hashes:
//     SX_COMPONENT_ID      defined in SxUtil.h, SxIPC.h, SxFoo.h
//     SX_EXPORT_COMPONENT  -D in Makefile / MSVC project
#define _SxLog1(MSG)                                                         \
    sxLog(SX_LOG_ID,SX_LOG_HASH,SX_FILE,__LINE__,SX_FUNC,                    \
         [&](const char *__sxLogId,uint32_t __sxLogHash,                     \
             const char *__sx_file,long __sx_line,                           \
             const char *__sx_func,auto& __sx_dummy)                         \
         {                                                                   \
            if (!SxLog::isEnabled (__sxLogHash))  return;                    \
            SX_UNUSED(__sx_dummy);                                           \
            SxLog::lock ();                                                  \
            SX_LOG_MSG(__sxLogId,MSG);                                       \
            int dbg = SxDebug::isEnabled (__sx_file, __sx_func, __sxLogId,   \
                                          NULL);                             \
            if (dbg != 0)                                                    \
               SX_PRINT_MSG(MSG,__sx_file,__sx_line,__sx_func,dbg)           \
            SxLog::unlock ();                                                \
    })
#define _SxLog2(CAPTURE,MSG)                                                 \
    sxLog(SX_LOG_ID,SX_LOG_HASH,SX_FILE,__LINE__,SX_FUNC,                    \
         [&](const char *__sxLogId,uint32_t __sxLogHash,                     \
             const char *__sx_file,long __sx_line,                           \
             const char *__sx_func,auto& CAPTURE)                            \
         {                                                                   \
            if (!SxLog::isEnabled (__sxLogHash))  return;                    \
            SX_UNUSED(CAPTURE);                                              \
            SxLog::lock ();                                                  \
            SX_LOG_MSG(__sxLogId,MSG);                                       \
            int dbg = SxDebug::isEnabled (__sx_file, __sx_func, __sxLogId,   \
                                          NULL);                             \
            if (dbg != 0)                                                    \
               SX_PRINT_MSG(MSG,__sx_file,__sx_line,__sx_func,dbg)           \
            SxLog::unlock ();                                                \
    })

#ifdef USE_SX_LOG
#   define SX_LOG(...) SX_VMACRO(_SxLog, __VA_ARGS__)
#else
#   define SX_LOG(...) ((void)0)
#endif



#endif /* _SX_LOG_H_ */

