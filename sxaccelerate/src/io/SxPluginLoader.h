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

#ifndef _SX_PLUGIN_LOADER_H_
#define _SX_PLUGIN_LOADER_H_

#ifdef WIN32
#  include <windows.h>
#endif

#include <SxString.h>
#include <SxIO.h>
#include <SxConfig.h>

/** \brief Load/unload program parts dynamically at run-time.

    \b SxPluginLoader = SPHInX Plugin Loader

    State-of-the-art operating systems allow loading and unloading
    of program parts at run-time. Such loadable part is usually called
    \e plugin. By using a plugin-mechanism
    -# the executable can become smaller (only required parts are loaded)
    -# program might required less memory (only required objects are created)
    -# functionality of the program can be changed during run-time

    A plugin has to be compiled as plugin:
    - under UNIX: lib*.so
    - under Windows: lib*.dll

    \author Sixten Boeck, boeck@mpie.de */
class SX_EXPORT_IO SxPluginLoader
{
   public:

      SxPluginLoader ();
     ~SxPluginLoader ();


      bool openSharedObject (const SxString &soName);
      void closeSharedObject ();
      typedef const char *(*SymPtr)(void);
      SymPtr getSymbol (const SxString &name);

   protected:

#     ifdef WIN32
         typedef HINSTANCE SoHandleType;
#     else
         typedef void *    SoHandleType;
#     endif

      SoHandleType soHandle;
      
};

#endif /* _SX_PLUGIN_LOADER_H_ */
