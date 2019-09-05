//#include <SxMpi.h>
#include <SxTaskGroup.h>
#include <SxMpiTaskGroup.h>
#include <SxParallelHierarchy.h>

#include <SxMpiTest.h>

#include <SxLoopMPI.h>
#include <SxString.h>

#include <SxPtr.h>

#define SUCCESS true
#define FAILURE false


/// Check a taskGroup for completeness (not for consistency)
bool check_task_group(SxTaskGroup*);

/// Basic test of SxParallelHierarchy, construction and manual destruction.
/// The existence of the default task group "top" is checked as well.
bool test_01()
{
   cout << "| running " << __FUNCTION__ << " ..." << endl;
   SxParallelHierarchy ph;
   ph.info();
   ph.write("ph-default.sx");
   SxTaskGroup *tg = ph.getTaskGroup("top");
   SX_CHECK( check_task_group(tg) );
   ph.destroy();
   ph.info();
   return SUCCESS;
}

/// Basic test of SX_MPI_LEVEL with automatic extension of SxParallelHierarchy.
bool test_02()
{
   cout << "| running " << __FUNCTION__ << " ..." << endl;
   SxParallelHierarchy ph;
   SX_MPI_LEVEL("topX");
   SX_MPI_LEVEL("topY");
   SX_MPI_LEVEL("topZ");
   ph.info();
   ph.write("ph-sx_mpi_level.sx");
   ph.destroy();
   return SUCCESS;
}

/// Usage examples of SxPtr in combination with the task group maps (obsolete).
bool test_03()
{
   cout << "| running " << __FUNCTION__ << " ..." << endl;
   SxPtr<SxMpiTaskGroup> tg = SxPtr<SxMpiTaskGroup>::create ();
   SxMap< SxString, SxPtr<SxMpiTaskGroup> > mtgMap;
   mtgMap(SxString("TEST")) = tg;
   SxMap< SxString, SxTaskGroup* > tgMap;
   tgMap(SxString("TEST")) = mtgMap(SxString("TEST")).getPtr ();
   return SUCCESS;
}

/// Basic tests of extending a parallel hierarchy with SxAutoLevel.
bool test_04()
{
   cout << "| running " << __FUNCTION__ << " ..." << endl;
   SxParallelHierarchy ph;
   {
      // A "context" can be passed to SxAutoLevel.
      // It is prefixed to globalName (== key to the SxTaskGroup maps).
      SxAutoLevel("auto")
            .append("foo-1", 12)
            .append("foo-2", 6);
   }
   ph.info();

   // The context is then used to fetch a task group from the parallel hierarchy as follows:
   SxTaskGroup * tg = ph.getTaskGroup("auto.foo-1");
   SX_CHECK( check_task_group(tg) );

   //
   tg = tg->getChild("foo-2");
   SX_CHECK( check_task_group(tg) );

   ph.write("ph-sx_auto_level.sx");

   ph.destroy();
   return SUCCESS;
}

//bool test_05()
//{
//   cout << "| running " << __FUNCTION__ << " ..." << endl;
//   SxParallelHierarchy ph;
//   SxTaskGroup *tg = ph.getTaskGroup("top");
//   SX_CHECK( check_task_group(tg) );
//   ph.destroy();
//   return SUCCESS;
//}



bool check_task_group(SxTaskGroup * tg)
{
   // --- check basic properties (numbers)
   SX_CHECK ( tg->getName().getSize () > 0 );
   SX_CHECK ( tg->getNmembers() > 0 );
   SX_CHECK ( tg->getNsiblings() > 0 );
   SX_CHECK ( tg->getNworkers() > 0 );
   SX_CHECK ( tg->getSiblingRank() >= 0 );
   SX_CHECK ( tg->getMemberRank() >= 0 );
   // --- check parents and children
   SX_CHECK ( tg->getParent() );
   SX_CHECK ( tg->getParentId().getSize () > 0);
   if ( tg->hasChildren() )
   {
      for (SxMap<SxString, SxTaskGroup*>::Iterator it = tg->children.begin();
            it != tg->children.end(); ++it )
      {
         SxTaskGroup * tg2 = it.getValue();
         SX_CHECK( tg2 );
      }
      SX_CHECK( tg->children.getSize() == tg->childId.getSize() );
   }
   return SUCCESS;
}



int main(int argc, char * argv[])
{
   SxLoopMPI::init(argc, argv);

   SX_CHECK( test_01() );
   SX_CHECK( test_02() );
   SX_CHECK( test_03() );
   SX_CHECK( test_04() );
//   SX_CHECK( test_05() );

   return 0;
}
