#ifndef _SX_DEMO_PLUGIN_BASE_H_
#define _SX_DEMO_PLUGIN_BASE_H_

#include <SxPlugin.h>

class SxDemoPluginBase : public SxPlugin
{
   public:
      SxDemoPluginBase (void *parent_) 
         : SxPlugin ((SxPluginManager *)parent_) { }
      virtual ~SxDemoPluginBase () { }
      
      virtual void foo ()=0;
};

#endif /* _SX_DEMO_PLUGIN_BASE_H_ */
