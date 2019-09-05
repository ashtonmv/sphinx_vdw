#include <SxPluginB.h>
#include <stdio.h>

SxPluginB::SxPluginB (void *parent_) : SxDemoPluginBase (parent_)
{
   //empty
}


SxPluginB::~SxPluginB ()
{
   printf ("About to destroy plugin B.\n");
}

void SxPluginB::foo ()
{
   printf ("This is plugin B.\n");
}