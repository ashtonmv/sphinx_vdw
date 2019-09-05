#ifndef _SX_PLUGIN_B_H_
#define _SX_PLUGIN_B_H_

#include <SxDemoPluginBase.h>

class SxPluginB : public SxDemoPluginBase
{
   public:
      SxPluginB (SxPluginManager *parent, void *data);
      virtual ~SxPluginB ();
      
      virtual void foo ();
};

REGISTER_SX_PLUGIN (SxPluginB, 1,0,0, "SxPluginB", "SxPluginB", "Test");

#endif /* _SX_PLUGIN_B_H_ */
