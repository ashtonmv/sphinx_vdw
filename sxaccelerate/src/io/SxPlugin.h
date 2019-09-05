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

#ifndef _SX_PLUGIN_H_
#define _SX_PLUGIN_H_

#include <SxString.h>
#include <SxIO.h>
#include <stdio.h>

class SxPlugin;
class SxPluginManager;

extern "C" SX_EXPORT_IO const char *getPluginIdentifier ();
extern "C" SX_EXPORT_IO void sxUnrefPlugin (SxPluginManager *, SxPlugin *);



/** \brief Dynamically loaded plugin

    \b SxPlugin = SFHIngX Plugin

    Classes derived from SxPlugin can be dynamically loaded from the
    SxPluginLoader. In order to support dynamic loading the corresponding
    object files have to be compiled with the -shared or -fPIC option.
    See for example the GNUmakefile in src/modules.

    Note that SxPlugins have to be registered using the REGISTER_SX_PLUGIN
    macro. It furthermore requires the definition of an static constructor.
\code
   class SxMyPlugin : public SxPlugin
   {
      public:

         SxMyPlugin (void *parent) : SxPlugin (parent) { }
         virtual ~SxMyPlugin ();
   };
   REGISTER_SX_PLUGIN("SxMyPlugin", 1,0,0, "Potentials", "MyPlugin",
                      "sphinx:1.1-1.5,someOtherModule:1.2|1.3,"
                      "anotherModule:1.5");
\endcode
    \author Sixten Boeck, boeck@mpie.de */
class SxPlugin
{
   public:
      SxPlugin (SxPluginManager *parent_, 
                void * /*userData*/ = NULL,
                const SxString &pluginId_="undefined")
         : parent(parent_),
           pluginId (pluginId_) { } 

      virtual ~SxPlugin () 
      { 
         if (parent)  sxUnrefPlugin (parent, this); 
      }

      void setPluginId (const SxString &id)  { pluginId = id; }
      SxString getPluginId () const { return pluginId; }

   protected:
      SxPluginManager *parent;
      SxString pluginId;
};


#define REGISTER_SX_PLUGIN(name,major,minor,patch,type,identifier,req)        \
   REGISTER_SX_PLUGIN_VER(name,major,minor,patch,type,identifier,req)

#define REGISTER_SX_PLUGIN_VER(name,major,minor,patch,type,identifier,req)    \
   extern "C"  {                                                              \
      SxPlugin    *create (SxPluginManager *p, void *d)   \
                                                 { return new name (p,d);   } \
      const char  *getPluginType ()              { return type;             } \
      const char  *getPluginIdentifier ()        { return identifier;       } \
      unsigned int getMajor ()                   { return major;            } \
      unsigned int getMinor ()                   { return minor;            } \
      unsigned int getPatch ()                   { return patch;            } \
      const char  *getRequirements ()            { return req;              } \
   }                                                                          \
   void sx_##name##_plugin ()


#endif /* _SX_PLUGIN_H_ */
