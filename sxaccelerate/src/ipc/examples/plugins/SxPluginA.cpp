#include <SxPluginA.h>
#include <stdio.h>

SxPluginA::SxPluginA (void *parent_) : SxDemoPluginBase (parent_)
{
   //empty
}


SxPluginA::~SxPluginA ()
{
   printf ("About to destroy plugin A.\n");
}

void SxPluginA::foo ()
{
   printf ("This is plugin A.\n");
}