#ifndef _SX_PARALLELHIERARCHY_H_
#define _SX_PARALLELHIERARCHY_H_

#include <SxTaskGroup.h>
#include <SxMpiTaskGroup.h>
#include <SxMap.h>
#include <SxString.h>
#include <SxSymbolTable.h>
#include <SxMpiComm.h>

class SX_EXPORT_MPI SxParallelHierarchy
{
   public:
      SxParallelHierarchy();
      SxParallelHierarchy(const SxString &hierarchyFile);

      /// tear down the parallel hierarchy manually
      void destroy();

      /// access a task group specified by the string label
      SxTaskGroup * getTaskGroup(const SxString &label);

      /// write the current hierarchy in sx-format
      static void write (const SxString &filename);

      /// print basic information about the hierarchy
      void info();

      /// store the instances of SxMpiTaskGroup via SxPtr to allow to modify the instances later
      static SxMap< SxString, SxPtr<SxMpiTaskGroup> > mpiTaskGroups;
      /// store plain C pointers to the SxTaskGroup objects in the map mpiTaskGroups
      static SxMap< SxString, SxTaskGroup* > taskGroups;
      static bool initialized;

   private:
      /// set up the parallel hierarchy from an input file
      void buildHierarchy(const SxString &hierarchyFile);

      /// add one node to the hierarchy from parser tree
      void build_hierarchy_recursively(SxSymbolTable * node, SxSymbolTable * parent);
      /// set up the taskGroups pointer map and create child-parent pointer links
      void build_hierarchy_pointers ();

      /// access an MPI task group specified by the string label
      SxMpiTaskGroup * getMpiTaskGroup(const SxString &label);
};

#endif
