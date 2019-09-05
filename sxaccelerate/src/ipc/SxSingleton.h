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

#ifndef _SX_SINGLETON_H_
#define _SX_SINGLETON_H_

#include <SxPtr.h>

/** \brief ...

    \b SxClass = S/PHI/nX ...

    ....

    \author Sixten Boeck, boeck@mpie.de */
template<class T>
class SxSingleton
{
   public:

      static T &getObj ();
      static SxPtr<T> &getObjPtr ();

      static void destroyKeyCB (void *);

};


#include <SxSingleton.hpp>

#endif /* _SX_SINGLETON_H_ */
