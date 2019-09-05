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

#ifndef _SX_PLUGIN_MANAGER_H_
#define _SX_PLUGIN_MANAGER_H_

#include <SxPluginLoader.h>
#include <SxIO.h>
#include <SxPtr.h>
#include <SxMap.h>
#include <SxList.h>
#include <SxString.h>
#include <SxPlugin.h>
#include <SxException.h>

//SX_EXPORT_IO void sxUnrefPlugin (void *, void *);


/** \brief ...

    \b SxPluginManager = SFHIngX Plugin Manager

    ....

    \author Sixten Boeck, boeck@mpie.de */
class SX_EXPORT_IO SxPluginManager
{
   public:

      enum LoadPolicy { UpdateOnly, Preload };

      SxPluginManager ();
     ~SxPluginManager ();

      void updatePluginList (const SxList<SxString> &folders,
                             const SxString &type="",
                             enum LoadPolicy loadPolicy=UpdateOnly);
      void loadPlugins (const SxList<SxString> &folders,
                        const SxString &type="");

      bool contains (const SxString &identifier);
      SxList<SxString> getPluginIds () const;

      template<class T>
         SxPtr<T> create (const SxString &identifier, // #exception
                          void *userData=NULL) const
         {
            cleanUp ();

            if ( !absPaths.containsKey(identifier) )  
               SX_THROW ("Plugin with identifier '" 
                        + identifier + "' not found in Plugin manager");

            SxPtr<SxPluginLoader> loader;
            if (!plugins.containsKey (identifier))  {
               loader = SxPtr<SxPluginLoader>::create ();
               SxString file = absPaths(identifier);
               loader->openSharedObject(file);
               plugins(identifier) = loader;
               pluginRefCounter(identifier) = 0;
            }  else  {
               loader = plugins(identifier);
            }

            T *(*instantiate)(void *, void *, const SxString &) = NULL;
            instantiate = (T *(*)(void *,void *, const SxString &))(
                             loader->getSymbol("create")
                          );

            const char *(*getId)() = NULL;
            getId = (const char *(*)())(
                      loader->getSymbol ("getPluginIdentifier")
                    );

            SxString pluginId = (*getId)();

            T *plugin = (*instantiate)(const_cast<SxPluginManager *>(this),
                                       userData, pluginId);
            pluginRefCounter(identifier)++;

            plugin->setPluginId (identifier);


            SxPtr<T> res;
            res.initFromCPointer (plugin);

            return res;
         }

      void unref (SxPlugin *);  // called from destructor of SxPlugin

      void print ();

   protected:

      SxMap<SxString,SxString>                       absPaths;
      mutable SxMap<SxString,SxPtr<SxPluginLoader> > plugins;
      mutable SxMap<SxString,int>                    pluginRefCounter;
      SxMap<SxString,LoadPolicy>                     pluginLoadPolicy;
      mutable SxPtr<SxPluginLoader>                  toBeUnloaded;

   private:
      void cleanUp () const;
};


#endif /* _SX_PLUGIN_MANAGER_H_ */
