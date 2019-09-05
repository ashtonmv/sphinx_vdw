#include <SxPluginB.h>
#include <stdio.h>

SxPluginB::SxPluginB (SxPluginManager *parent_, void *)
   : SxDemoPluginBase (parent_)
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