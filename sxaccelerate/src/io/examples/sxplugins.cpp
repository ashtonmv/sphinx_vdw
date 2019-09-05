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
#include <SxPluginManager.h>
// --- define the common plugin base (our plugin API)
#include "plugins/SxDemoPluginBase.h"


/**  \author  Sixten Boeck  */
int main ()
{
   // --- prepare plugin manager and load list of plugins
   //     note: only the plugin list is loaded, not the plugins themselves
   SxPluginManager pluginMgr;
   pluginMgr.updatePluginList (SxList<SxString>() << "plugins/.libs");
   pluginMgr.print ();

   // --- load and execute plugin A
try {
   SxPtr<SxDemoPluginBase> plugin = pluginMgr.create<SxDemoPluginBase> ("SxPluginA");
   plugin->foo ();

   // --- unload plugin A (thanks to SxPtr), load and execute plugin A
   plugin = pluginMgr.create<SxDemoPluginBase> ("SxPluginB");
   plugin->foo ();
} catch (SxException e)  {
   e.print ();
}

   return 0;
}
