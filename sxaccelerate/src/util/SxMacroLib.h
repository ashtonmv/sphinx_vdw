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

#ifndef _SX_MACRO_LIB_H_
#define _SX_MACRO_LIB_H_

#include <SxConfig.h>

// --------------------------------------------------------------------------
// SX_UNUSED (var1, var2, ...)
//
// Avoid compiler warning of an unused variable
// Usage:
//    void foo (int dummy) {
//       SX_UNUSED(dummy);
//    }
// --------------------------------------------------------------------------
template<class ...Args> void SX_UNUSED(Args&& ... args)
{
    (void)(sizeof...(args));
}




// --------------------------------------------------------------------------
// SX_UNIQUE_ID (abc)
//
// Define unique identifiers
// Usage:
//     int SX_UNIQUE_ID(abc)  => int abc##__LINE##
// --------------------------------------------------------------------------
#define SX_UNIQUE_ID_NAME(var,x)  var ## x
#define SX_UNIQUE_ID_WRAP(var,x)  SX_UNIQUE_ID_NAME(var,x)
#define SX_UNIQUE_ID(var)         SX_UNIQUE_ID_WRAP(var,__LINE__)



// --------------------------------------------------------------------------
// SX_BLOCK (SxMutex)
//
// Define unique blocks
// Usage:
//    SERIALIZE {
//       ...
//    }
//    define SX_SERIALIZE  SX_BLOCK(SxMutex)
// --------------------------------------------------------------------------
#define SX_BLOCK(type) \
   bool SX_UNIQUE_ID(cond)=true; for (type SX_UNIQUE_ID(var); SX_UNIQUE_ID(cond); SX_UNIQUE_ID(cond) = false)

#define SX_BLOCK_ARG(type,arg) \
   bool SX_UNIQUE_ID(cond)=true; for (type SX_UNIQUE_ID(var)(arg); SX_UNIQUE_ID(cond); SX_UNIQUE_ID(cond) = false)




// --------------------------------------------------------------------------
// SX_VMACRO (macroSet, __VA_ARGS__)
//
// Environment to provide variadic macros
// Usage:
//    #define _MyMacro0() ...
//    #define _MyMacro1(a) ...
//    #define _MyMacro2(a,b) ...
//    #define _MyMacro3(a,b,c) ...
//    
//    #define MyMacro(...) SX_VMACRO(_MyMacro, __VA_ARGS__)
//
// Note:    
//    MSVC expands __VA_ARGS__ according to C99. Hence, proper variadic macro
//    unrolling doesn't work:
//    http://connect.microsoft.com/VisualStudio/feedback/details/380090/variadic-macro-replacement
// --------------------------------------------------------------------------
#define _SxVMacroGlue(x,y)  x y
#define _SxVMacroGetNArgs(_1,_2,_3,_4,_5,_6,_7,_8,_9,nArgs,...) nArgs
#define _SxVMacroExpandArgs(args)  _SxVMacroGetNArgs args
#define _SxVMacroCountArgs(...)                                               \
           _SxVMacroExpandArgs((__VA_ARGS__,9,8,7,6,5,4,3,2,1,0))

#define _SxVMacroOvl2(name,nArgs) name##nArgs
#define _SxVMacroOvl1(name,nArgs) _SxVMacroOvl2(name,nArgs)
#define _SxVMacroOvl(name,nArgs)  _SxVMacroOvl1(name,nArgs)

#define SX_VMACRO(name,...)                                                   \
          do  {                                                               \
             _SxVMacroGlue(_SxVMacroOvl(name,_SxVMacroCountArgs(__VA_ARGS__)),\
                           (__VA_ARGS__));                                    \
          } while (0)

#define SX_VMACRO_FUNC(name,...)                                              \
          _SxVMacroGlue(_SxVMacroOvl(name,_SxVMacroCountArgs(__VA_ARGS__)),   \
                        (__VA_ARGS__))






// --------------------------------------------------------------------------
//
// Extract an individual argument from __VA_ARGS__
//
// Usage:
//    #SX_VMACRO_ARG2(__VA_ARGS__)
// --------------------------------------------------------------------------
#define SX_VMACRO_ARG0(ARG,...)                                ARG
#define SX_VMACRO_ARG1(_0,ARG,...)                             ARG
#define SX_VMACRO_ARG2(_0,_1,ARG,...)                          ARG
#define SX_VMACRO_ARG3(_0,_1,_2,ARG,...)                       ARG
#define SX_VMACRO_ARG4(_0,_1,_2,_3,ARG,...)                    ARG
#define SX_VMACRO_ARG5(_0,_1,_2,_3,_4,ARG,...)                 ARG
#define SX_VMACRO_ARG6(_0,_1,_2,_3,_4,_5,ARG,...)              ARG
#define SX_VMACRO_ARG7(_0,_1,_2,_3,_4,_5,_6,ARG,...)           ARG
#define SX_VMACRO_ARG8(_0,_1,_2,_3,_4,_5,_6,_7,ARG,...)        ARG
#define SX_VMACRO_ARG9(_0,_1,_2,_3,_4,_5,_6,_7,_8,ARG,...)     ARG




// --------------------------------------------------------------------------
// ATTR_NO_RETURN
//
// ATTR_NO_RETURN indicates that the function never returns
// this removes many false alarms about possibly uninitialized vars
// --------------------------------------------------------------------------
#ifdef ATTR_NO_RETURN
#  error ATTR_NO_RETURN macro name collision - find a better name!
#endif
#ifdef WIN32
#  define ATTR_NO_RETURN
#  define DECL_NO_RETURN __declspec(noreturn)
#else
#  define DECL_NO_RETURN
#  ifdef __GNUC__
#    define ATTR_NO_RETURN __attribute__ ((noreturn))
#  else
#    define ATTR_NO_RETURN
#  endif
#endif


#endif /* _SX_MACRO_LIB_H_ */
