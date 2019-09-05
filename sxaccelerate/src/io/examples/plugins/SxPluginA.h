#ifndef _SX_PLUGIN_A_H_
#define _SX_PLUGIN_A_H_

#include <SxDemoPluginBase.h>

class SxPluginA : public SxDemoPluginBase
{
   public:
      SxPluginA (SxPluginManager *parent, void *data);
      virtual ~SxPluginA ();
      
      virtual void foo ();
};

REGISTER_SX_PLUGIN (SxPluginA, 1,0,0, "SxPluginA", "SxPluginA", "Test");

#endif /* _SX_PLUGIN_A_H_ */
