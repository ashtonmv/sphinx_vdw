#include <SxPluginA.h>
#include <stdio.h>

SxPluginA::SxPluginA (SxPluginManager *parent_, void *)
   : SxDemoPluginBase (parent_)
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