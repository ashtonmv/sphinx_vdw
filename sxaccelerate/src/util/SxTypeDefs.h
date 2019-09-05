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

/** Portable C types
    \file SxTypeDefs.h
    \ingroup group_os
    \author Vaclav Bubnik, bubnik@mpie.de */

#ifndef _SX_TYPE_DEFS_H_
#define _SX_TYPE_DEFS_H_

#include <SxConfig.h>

// --- INTn_MIN, INTn_MAX C99 standard, in C++ on request.
#ifndef __STDC_LIMIT_MACROS
#  define __STDC_LIMIT_MACROS 1
#endif

/*
#ifdef WIN32
#ifdef INT8_MIN
#  undef INT8_MIN
#endif
#ifdef INT16_MIN
#  undef INT16_MIN
#endif
#ifdef INT32_MIN
#  undef INT32_MIN
#endif
#ifdef INT8_MAX
#  undef INT8_MAX
#endif
#ifdef INT16_MAX
#  undef INT16_MAX
#endif
#ifdef INT32_MAX
#  undef INT32_MAX
#endif
#ifdef UINT8_MAX
#  undef UINT8_MAX
#endif
#ifdef UINT16_MAX
#  undef UINT16_MAX
#endif
#ifdef UINT32_MAX
#  undef UINT32_MAX
#endif
#endif */  /* WIN32 */

#include <stdint.h>

#ifdef MSVC
//#  undef UINT8_MAX
//#  define UINT8_MAX  255U
// --- [util/SxWinConfig.h:45]: typedef long ssize_t;
#else
/* size_t ISO/ANSI standard C
   ssize_t not a standard */
#  include <sys/types.h>
#endif /* MSVC */

#include <climits>
#include <ctype.h>

/** Casting between pointer-to-function and pointer-to-object.

    There is no valid cast between pointer to function and pointer to object:
    http://www.trilithium.com/johan/2004/12/problem-with-dlsym/
    
    Use it to remove warning:
       ISO C++ forbids casting between pointer-to-function and pointer-to-object
*/
template<class T>
inline T sxFnCast (void *ptr_)
{
   return reinterpret_cast<T> (reinterpret_cast<size_t> (ptr_));
}

/** numeric limits */
template<class T>
inline T sxmin ()
{
   return 0;
}
/** numeric limits */
template<class T>
inline T sxmax ()
{
   return 0;
}

// --- signed char
template<> inline signed char sxmin<signed char> ()
{
   return SCHAR_MIN;
}
template<> inline signed char sxmax<signed char> ()
{
   return SCHAR_MAX;
}

// --- int
template<> inline int sxmin<int> ()
{
   return INT_MIN;
}
template<> inline int sxmax<int> ()
{
   return INT_MAX;
}

// --- long
template<> inline long sxmin<long> ()
{
   return LONG_MIN;
}
template<> inline long sxmax<long> ()
{
   return LONG_MAX;
}

// --- long long
template<> inline long long sxmin<long long> ()
{
   return LLONG_MIN;
}
template<> inline long long sxmax<long long> ()
{
   return LLONG_MAX;
}

template<> inline unsigned char sxmax<unsigned char> ()
{
   return UCHAR_MAX;
}
template<> inline unsigned int sxmax<unsigned int> ()
{
   return UINT_MAX;
}
template<> inline unsigned long sxmax<unsigned long> ()
{
   return ULONG_MAX;
}
template<> inline unsigned long long sxmax<unsigned long long> ()
{
   return ULLONG_MAX;
}

#endif /* _SX_TYPE_DEFS_H_ */

