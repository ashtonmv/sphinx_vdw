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

#include <SxPluginManager.h>
#include <SxPlugin.h>
#include <SxFSAction.h>

void sxUnrefPlugin (SxPluginManager *parent, SxPlugin *plugin)
{
   parent->unref (plugin);
}


SxPluginManager::SxPluginManager ()
{
   // empty
}


SxPluginManager::~SxPluginManager ()
{
   cleanUp ();
   int i=0;
   SxList<SxString> ids = plugins.getKeys ();
   SxList<SxString>::ConstIterator it;
   for (it = ids.begin(); it != ids.end(); ++it, ++i)  {

      // --- we cannot close a shared object if it is still in use!
      SX_CHECK (pluginRefCounter(*it) == 0, 
                pluginRefCounter(*it), ids(i));

      plugins(*it)->closeSharedObject ();
   }
}


void SxPluginManager::updatePluginList (const SxList<SxString> &folders,
                                        const SxString &type,
                                        enum LoadPolicy loadPolicy)
{
   cleanUp ();
   // --- common interface to plugins
   const char *(*getType)() = NULL;
   const char *(*getId)()   = NULL;
// const char *(*getReq)()  = NULL;

   SxPluginLoader loader;
   SxPtr<SxPluginLoader> preLoader;

   SxDir dir;
   SxList<SxString>::ConstIterator folderIt;
   SxList<SxFileInfo>::ConstIterator fileIt;
   SxList<SxFileInfo> files;
   SxString pluginId, pluginType, filename;
   SxString soPattern;
#  ifdef WIN32
      soPattern = "*.dll";
#  else
      soPattern = "*.so";
#  endif
   for (folderIt = folders.begin(); folderIt != folders.end(); ++folderIt)  {
      
      files = SxFSAction::ls ((*folderIt) + 
                              SxString (SxFileInfo::getSeparator ()) + 
                              soPattern);
      if (files.getSize() == 0)  
           files = SxFSAction::ls (*folderIt +
                                   SxString (SxFileInfo::getSeparator ()) +
                                   ".libs" +
                                   SxString (SxFileInfo::getSeparator ()) +
                                   soPattern);
      for (fileIt = files.begin(); fileIt != files.end(); ++fileIt)  {
         if (SxFSAction::test_f ((*fileIt).getAbsPath ()))  {
            if (fileIt->getName().tokenize('.').getSize() == 2)  { // skip symlinks
               //FIXME: better use fileIt->isSymLink (makes problems under macosx)
               filename = fileIt->getAbsPath();
               sxprintf ("checking %s... ", filename.ascii()); 
               fflush (stdout);

               if (loader.openSharedObject (filename))  {

                  // --- retrieve common symbols
                  getType = (const char *(*)())loader.getSymbol ("getPluginType");
                  getId   = (const char *(*)())loader.getSymbol ("getPluginIdentifier");
//                getReq  = (const char *(*)())loader.getSymbol ("getRequirements");

                  pluginId   = (*getId)();
                  pluginType = (*getType)();

                  loader.closeSharedObject ();

                  if ( type == "" || type == pluginType )  {

                     // if (checkDependencies ( (*getRequirements)() ))  {
                     absPaths(pluginId)         = filename;
                     pluginLoadPolicy(pluginId) = loadPolicy;

                     // --- preload plugin if required
                     if (loadPolicy == Preload)  {
                        if (!plugins.containsKey (pluginId))  {
                           preLoader = SxPtr<SxPluginLoader>::create ();
                           preLoader->openSharedObject(filename);
                           plugins(pluginId) = preLoader;
                           pluginRefCounter(pluginId) = 0;
                        } 
                        sxprintf ("Loaded\n");
                     }  else  {
                        sxprintf ("OK\n");
                     }

                     // }  else  {
                     //    sxprintf ("wrong dependencies\n");
                     // }
                  }
               }
            } 
         }

      }
   }
}

bool SxPluginManager::contains (const SxString &identifier)
{
   return absPaths.containsKey (identifier);
}

SxList<SxString> SxPluginManager::getPluginIds () const
{
   return absPaths.getKeys ();
}


void SxPluginManager::unref (SxPlugin *p)
{
   SX_CHECK (p);
   cleanUp ();
   SxString id = p->getPluginId();
// sxprintf ("SxPluginManager::unref %p (%s)\n", p, id.ascii());

   pluginRefCounter(id)--;

   if (pluginLoadPolicy(id) == Preload)  {
      if (pluginRefCounter(id) == -1)  pluginRefCounter(id) = 0;
      return;
   }
   
   if (pluginRefCounter(id) == 0)  {
      // --- we cannot unload the shared object right away as this
      //     function was called from the plugin's destructor.
      //     the shared object must be be unloaded as long we are 
      //     still inside the plugin's destructor!
      toBeUnloaded = plugins(id);
      plugins.removeKey (id);
      pluginRefCounter.removeKey (id);
   }
}

void SxPluginManager::print ()
{
   cleanUp ();
   SxList<SxString> ids = absPaths.getKeys ();
   SxList<SxString>::ConstIterator it;
   SxString state, policy;
   int nRefs;
   for (it = ids.begin(); it != ids.end(); ++it)  {
      if (pluginRefCounter.containsKey (*it))  {
         nRefs = pluginRefCounter(*it);
         if      (nRefs == 0)  
            state = "preloaded";
         else if (pluginRefCounter(*it) >  0)  
            state = SxString("loaded (") 
                  + SxString (pluginRefCounter(*it))
                  + ")";
         else
            SX_EXIT;
      }  else  {
         state = "available";
      }
      policy = (pluginLoadPolicy(*it) == Preload) ? "+" : " ";
      sxprintf ("%-10s (%s)| %-10s| %s\n",
                (*it).ascii(),
                policy.ascii(),
                state.ascii(),
                absPaths(*it).ascii());
   }

}


void SxPluginManager::cleanUp () const
{
   if (toBeUnloaded != SxPtr<SxPluginLoader> ())  {
      toBeUnloaded->closeSharedObject ();
      toBeUnloaded = SxPtr<SxPluginLoader> ();
   }
}
