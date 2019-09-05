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
#ifndef _SX_ERROR_H_
#define _SX_ERROR_H_

#include <SxConfig.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <signal.h>
#include <SxUtil.h>
#include <SxMacroLib.h>

#ifdef WIN32

#  define _WINSOCKAPI_ // prevent winsock.h to be included from windows.h
                       // otherwise later include of winsock2.h would show
                       // redefinition errors
#  include <windows.h> // DWORD, GetLastError()
#endif /* WIN32 */


// --- call stack
#ifdef WIN32
#  pragma warning (push)
#  pragma warning (disable:4090)
#  pragma warning (disable:4091)
#  include <TlHelp32.h>
#  include <DbgHelp.h>
#  pragma warning (pop)
#  pragma comment(lib, "Dbghelp.lib")
#endif


#ifdef SX_IOS
// CFErrorRef, CFStringRef, CFStringGetLength(), CFErrorGetCode(), ...
#  include <CoreFoundation/CoreFoundation.h>
#endif /* SX_IOS */

#ifndef WIN32
#  define SX_FILE                                                        \
         (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/')+1  : __FILE__)
#else
#  define SX_FILE                                                        \
         (strrchr(__FILE__, '\\') ? strrchr(__FILE__, '\\')+1 : __FILE__)
#endif

/** \brief Error handling functions

    This class is a container of various static functions which are
    called in DEBUG mode in case the application crashes in order to
    provide the developer with detailed information about the issue
    such as:
       - function, source line of error
       - call stack
    \author Sixten Boeck */
class SX_EXPORT_UTIL SxDebug
{
   public:


/** \brief Test if debugging is enabled.
    \return 1 if debugging is enabled, 0 otherwise.

    \par Environment variables:
    \li SX_DEBUG_FILE specifies a comma separated list of filenames
    \li SX_DEBUG_FUNCTION specifies a comma separated list of function names

    \par Disable all debug messages:
\verbatim
> export SX_DEBUG_FILE=none
\endverbatim

    \par Example file main.cpp:
\code
#include <SxError.h>

class SxMyDate
{
   public:
      int year;
      SxMyDate ();
      void setYear (int year);
      int getYear () const;
};

SxMyDate::SxMyDate ()
   : year (1970)
{
   SX_DBG_MSG ("default year="<<year);
}

void SxMyDate::setYear (int year_)
{
   year = year_;
   SX_DBG_MSG ("year="<<year);
}

int SxMyDate::getYear () const
{
   SX_DBG_MSG ("year="<<year);
   return year;
}

int main (int, char **)
{
   SX_TRACE ();

   SxMyDate date;
   date.setYear (2010);
   date.getYear ();
   return 0;
}
\endcode

   \par Output:
   \li with > export SX_DEBUG_FILE=
   \li or with > export SX_DEBUG_FILE=main.cpp
   \li or with > export SX_DEBUG_FILE=main.cpp,SxError.h
\verbatim
main.cpp main: BEGIN
main.cpp SxMyDate: default year=1970
main.cpp setYear: year=2010
main.cpp getYear: year=2010
main.cpp main: END 1.359 ms
\endverbatim

   \li with > export SX_DEBUG_FUNCTION=setYear
\verbatim
main.cpp setYear: year=2010
\endverbatim

   \li with > export SX_DEBUG_FUNCTION=setYear:getYear
\verbatim
main.cpp setYear: year=2010
main.cpp getYear: year=2010
\endverbatim */
      static int isEnabled (const char *file, const char *func,
                            const char *feature);

      static int isEnabled (const char *file, const char *func,
                            const char *logId, const char *feature);

      static bool dump (); ///< dump current call stack

      static std::ostream &tag (std::ostream &s,
                                const char *file, long line, const char *func,
                                int tron);
#  ifdef WIN32
      static std::wostream &tag(std::wostream &s,
                                const char *file, long line, const char *func,
                                int tron);
#  endif
      static void log (int level, const char *tag, const char *func,
                       const char *fmt, ...);


   protected:

      static void stackDump (int skipLvls=2);

#     ifdef WIN32
         static BOOL CALLBACK filterCB (PVOID inParam,
                                        const PMINIDUMP_CALLBACK_INPUT,
                                        PMINIDUMP_CALLBACK_OUTPUT);

         static DWORD CALLBACK writeCB (LPVOID inParam);

         static void enumerateThreads (DWORD (WINAPI *cb)(HANDLE),
                                       DWORD skipThreadId);
#     endif
};

// --------------------------------------------------------------------------

extern void SX_EXPORT_UTIL DECL_NO_RETURN sxExit (int = 1) ATTR_NO_RETURN;
extern void SX_EXPORT_UTIL DECL_NO_RETURN sxQuit ();

// --- GCC >= 2.0 provides these two magic variables which hold the name
//     of the current function, as a string.
#define SX_FUNCTION __FUNCTION__
#ifdef WIN32
#  define SX_FUNC   __FUNCTION__  // __FUNCSIG__
#else
#  define SX_FUNC   __PRETTY_FUNCTION__
#endif


// --------------------------------------------------------------------------
// SX_EXIT
// --------------------------------------------------------------------------

#ifdef NDEBUG
#   define SX_EXIT                                                           \
         { std::cerr << "\n+---------------------------------------------"   \
                     << "--------------------------------\n";                \
           std::cerr << "| Abnormal program stop in file " << __FILE__       \
                     << " at line " << __LINE__ << ".\n";                    \
           std::cerr <<"+---------------------------------------------"      \
                     << "--------------------------------\n";                \
           sxExit(1); }
#else
#   define SX_EXIT  sxExit(1);
#endif


// --------------------------------------------------------------------------
// SX_CHECK(cond,var1,var2,...,varN)
// --------------------------------------------------------------------------

#define _SxCheck1(cond)                                                     \
           if ( !(cond) ) {                                                 \
              std::cerr << std::endl << "ASSERTATION FAILED in "            \
                        << __FILE__ << ":" << __LINE__ << std::endl         \
                        << "                   " << SX_FUNC                 \
                        << std::endl << #cond << std::endl; sxExit(1); }
#define _SxCheck2(cond,var1)                                                \
           if ( !(cond) ) {                                                 \
              std::cerr << std::endl << "ASSERTATION FAILED in "            \
                        << __FILE__ << ":" << __LINE__ << std::endl         \
                        << "                   " << SX_FUNC                 \
                        << std::endl << #cond << ", "                       \
                        << #var1 << "=" << (var1)                           \
                        << std::endl; sxExit(1); }
#define _SxCheck3(cond,var1,var2)                                           \
           if ( !(cond) ) {                                                 \
              std::cerr << std::endl << "ASSERTATION FAILED in "            \
                        << __FILE__ << ":" << __LINE__ << std::endl         \
                        << "                   " << SX_FUNC                 \
                        << std::endl << #cond << ", "                       \
                        << #var1 << "=" << (var1) << ", "                   \
                        << #var2 << "=" << (var2)                           \
                        << std::endl; sxExit(1); }
#define _SxCheck4(cond,var1,var2,var3)                                      \
           if ( !(cond) ) {                                                 \
              std::cerr << std::endl << "ASSERTATION FAILED in "            \
                        << __FILE__ << ":" << __LINE__ << std::endl         \
                        << "                   " << SX_FUNC                 \
                        << std::endl << #cond << ", "                       \
                        << #var1 << "=" << (var1) << ", "                   \
                        << #var2 << "=" << (var2) << ", "                   \
                        << #var3 << "=" << (var3)                           \
                        << std::endl; sxExit(1); }
#define _SxCheck5(cond,var1,var2,var3,var4)                                 \
           if ( !(cond) ) {                                                 \
              std::cerr << std::endl << "ASSERTATION FAILED in "            \
                        << __FILE__ << ":" << __LINE__ << std::endl         \
                        << "                   " << SX_FUNC                 \
                        << std::endl << #cond << ", "                       \
                        << #var1 << "=" << (var1) << ", "                   \
                        << #var2 << "=" << (var2) << ", "                   \
                        << #var3 << "=" << (var3) << ", "                   \
                        << #var4 << "=" << (var4)                           \
                        << std::endl; sxExit(1); }

#ifdef NDEBUG
#   define SX_CHECK(...)   ((void)0)
#else
#   define SX_CHECK(...)   SX_VMACRO(_SxCheck, __VA_ARGS__)
#endif

#   define SX_TAG   SxString(SxString(SX_FILE)+":"+SxString(__LINE__)\
                             +":"+SX_FUNC).ascii()
// --------------------------------------------------------------------------
// Deprecate macros:
//    SX_CHECK_VAR
//    SX_CHECK_VARS
//    SX_CHECK_3VARS
//    SX_CHECK_RANGE
// Note: please use variadic SX_CHECK(...) instead
// --------------------------------------------------------------------------


#ifdef NDEBUG
#   define SX_CHECK_VAR(expr,var)              ((void)0)
#   define SX_CHECK_VARS(expr,var1,var2)       ((void)0)
#   define SX_CHECK_3VARS(expr,var1,var2,var3) ((void)0)
#   define SX_CHECK_RANGE(idx,size)            ((void)0)
#else
#   define SX_CHECK_VAR(expr,var)  _SxCheck2(expr,var)
#   define SX_CHECK_VARS(expr,var1,var2)                                     \
              _SxCheck3(expr,var1,var2)
#   define SX_CHECK_3VARS(expr,var1,var2,var3)                               \
              _SxCheck4(expr,var1,var2,var3)
#   define SX_CHECK_RANGE(idx,size)                                          \
              _SxCheck3((idx) >= 0 && (idx) < (size), idx, size)
#endif


// --------------------------------------------------------------------------
// SX_CHECK_ERROR(errCode)
// --------------------------------------------------------------------------
#define SX_CHECK_ERR(expr)                                                    \
    {                                                                         \
       int errcode = expr;                                                    \
       if (errcode < 0)  {                                                    \
          int err = errno;                                                    \
          std::cerr << "ERROR CODE " << errcode << " in "                     \
                    << __FILE__ << ", line " << __LINE__ << "!\n";            \
          if (errcode == -1)  errcode = err;                                  \
          std::cerr << strerror (errcode) << std::endl;                       \
                                                                              \
          SX_EXIT;                                                            \
       }                                                                      \
    }



// --------------------------------------------------------------------------
// SX_BREAK
// --------------------------------------------------------------------------

// Note: sxBreak is also defined in SxIterator.h
inline void SX_EXPORT_UTIL sxBreak (const char *file, long line, const char *func)
{
#  ifdef SX_ANDROID
      SX_EXIT; // not yet implemented
#  else
#     ifdef NDEBUG
         SX_UNUSED (file, line, func);
#     else
         char buf[BUFSIZ];
         do {
            std::cout << "BREAK in " << file << ":" << line << ": "
                                     << func << std::endl;
            std::cout << "Hit <Enter> to continue, (t)race, or dump (c)ore\n";
            std::cout.flush ();
            fgets(buf, sizeof(buf), stdin);
            if      (buf[0] == 't')  SxDebug::dump ();
            else if (buf[0] == 'c')  SX_EXIT;
         } while (buf[0] == 't');

#     endif
#  endif /* SX_ANDROID */
}

#define SX_BREAK sxBreak(SX_FILE,__LINE__,SX_FUNC)





// --------------------------------------------------------------------------
// SX_QUIT
// --------------------------------------------------------------------------
#ifdef NDEBUG
#   define SX_QUIT                                                           \
       { std::cout << "\n+---------------------------------------------"     \
                   << "--------------------------------\n";                  \
         std::cout << "| Program stopped in file " << __FILE__               \
                   << " at line " << __LINE__ << "\n"                        \
                   << "| Invalid or inconsistent input as described "        \
                   << "above.\n";                                            \
         std::cout << "+---------------------------------------------"       \
                   << "--------------------------------\n";                  \
         sxExit(1); }
#else
#   define SX_QUIT                                                           \
       { std::cout << "\n+---------------------------------------------"     \
                   << "--------------------------------\n";                  \
         std::cout << "| Program stopped in file " << __FILE__               \
                   << " at line " << __LINE__ << "\n"                        \
                   << "| Invalid or inconsistent input as described "        \
                   << "above.\n";                                            \
         std::cout << "+---------------------------------------------"       \
                   << "--------------------------------\n";                  \
         sxQuit (); }
#endif




// --------------------------------------------------------------------------
// SX_DBG_MSG(msg)
// SX_DBG_MSG(msg << "abc" << "def" << i);
// iterator.SX_DBG_MSG(it, "val=" << *it).foo ();
// --------------------------------------------------------------------------

#ifdef WIN32
#  define SX_PRINT_MSG(MSG,...)  \
   SxDebug::tag(std::wcerr,__VA_ARGS__) << MSG << std::endl;
#else
#  define SX_PRINT_MSG(MSG,...)  \
   SxDebug::tag(std::cerr,__VA_ARGS__) << MSG << std::endl;
#endif

// --------------------------------------------------------------------------

// Note: sxDbgMsg() is also defined in SxIterator.h
template<class Fn>
void sxDbgMsg (const char *file, long line, const char *func, Fn f)
{ 
#  ifdef NDEBUG
      SX_UNUSED (file, line, func, f);
#  else
      static const int dummy = 0;
      f (file, line, func, dummy);
#  endif
}


#define _SxDbgMsg1(MSG)                                                      \
    sxDbgMsg(SX_FILE,__LINE__,SX_FUNC,                                       \
             [&](const char *__sx_file,long __sx_line,                       \
                 const char *__sx_func, auto& __sx_dummy)                    \
             {                                                               \
                int dbg = SxDebug::isEnabled (__sx_file, __sx_func, NULL);   \
                if (dbg == 0)  return;                                       \
                SX_UNUSED(__sx_dummy);                                       \
                SxUtil::getGlobalObj().lockLog ();                           \
                SX_PRINT_MSG(MSG,__sx_file,__sx_line,__sx_func,dbg)          \
                SxUtil::getGlobalObj().unlockLog ();                         \
    })
#define _SxDbgMsg2(CAPTURE,MSG)                                              \
    sxDbgMsg(SX_FILE,__LINE__,SX_FUNC,                                       \
             [&](const char *__sx_file,long __sx_line,                       \
                 const char *__sx_func, auto& CAPTURE)                       \
             {                                                               \
                int dbg = SxDebug::isEnabled (__sx_file, __sx_func, NULL);   \
                if (dbg == 0)  return;                                       \
                SX_UNUSED(CAPTURE);                                          \
                SxUtil::getGlobalObj().lockLog ();                           \
                SX_PRINT_MSG(MSG,__sx_file,__sx_line,__sx_func,dbg)          \
                SxUtil::getGlobalObj().unlockLog ();                         \
    })

#define SX_DBG_MSG(...) SX_VMACRO(_SxDbgMsg, __VA_ARGS__)



// --------------------------------------------------------------------------
// SX_TRACE ()
// SX_TRACE ("some dbg msg")
//
// Function-like Macro SX_TRACE ();
// http://gcc.gnu.org/onlinedocs/cpp/Function_002dlike-Macros.html
// --------------------------------------------------------------------------

class SX_EXPORT_UTIL SxTrace
{
   public:
      SxTrace (int enabled_, const char *file, long line,
               const char *func, const char *funcShort,
               const char *dbgMsg=NULL);
      ~SxTrace ();

   protected:
      int enabled;
      char filename[80];
      long line;
      char functionName[80];
      char functionNameShort[80];
      char dbgMsg[80];
      double t1;
      double t2;
};

#ifdef NDEBUG
#   define SX_TRACE(...) ((void)0)
#else
#   define SX_TRACE(...) \
           static int _tron = SxDebug::isEnabled (SX_FILE, SX_FUNCTION, \
                                                  NULL);                \
           SxTrace sxTrace (_tron,SX_FILE,__LINE__,SX_FUNC,SX_FUNCTION, \
			   {__VA_ARGS__});
#endif



// --------------------------------------------------------------------------

class SxString;
/** Get an error message from a Windows-specific system error code resp. a POSIX
    errno value. Beware that intermediate library calls, syscalls etc. may
    overwrite the values. It is better to store the value immediately after a
    failure in a local variable 'err' and call sxstrerror (err). */
extern SxString SX_EXPORT_UTIL sxstrerror ();
/** Get an error message from an errno error code (formally errno_t) */
extern SxString SX_EXPORT_UTIL sxstrerror (int errnum);

#ifdef WIN32
/** Get an error message for a Windows-specific system error code (last-error
    code;
    <http://msdn.microsoft.com/en-us/library/ms680347(v=vs.85).aspx>,
    <http://msdn.microsoft.com/en-us/library/ms681381(v=vs.85).aspx>)
\code
   BOOL status = WriteFile (...);
   if (!status)  {
      DWORD err = ::GetLastError ();
      cout << "ERROR: " << sxstrerror (err) << endl;
   }
\endcode
*/
extern SxString SX_EXPORT_UTIL sxstrerror (DWORD sysErrCode);
#endif /* WIN32 */

#ifdef SX_IOS
extern SxString SX_EXPORT_UTIL sxstrerror (CFErrorRef error);
#endif /* SX_IOS */

/// This auxiliary function states that we ran out of memory and sxExits
extern void SX_EXPORT_UTIL DECL_NO_RETURN sxOutOfMemoryHandler () ATTR_NO_RETURN;

#endif /* _SX_ERROR_H_ */
