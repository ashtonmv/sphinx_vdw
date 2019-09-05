#include <SxMpiComm.h>

#include <SxTaskGroup.h>
#include <SxMpiTaskGroup.h>
#include <SxParallelHierarchy.h>

#include <SxMap.h>
#include <SxString.h>
#include <SxList.h>
#include <SxParser.h>
#include <SxSymbolTable.h>
#include <SxError.h>
#include <SxLoopMPI.h>


SxMap<SxString, SxTaskGroup*> SxParallelHierarchy::taskGroups
   = SxMap<SxString, SxTaskGroup*>();
SxMap< SxString, SxPtr<SxMpiTaskGroup> > SxParallelHierarchy::mpiTaskGroups
   = SxMap< SxString, SxPtr<SxMpiTaskGroup> >();
bool SxParallelHierarchy::initialized = false;


SxParallelHierarchy::SxParallelHierarchy()
{
   if (initialized) return;

   if (SxFile ("parallelHierarchy.sx").exists ())
   {
      buildHierarchy("parallelHierarchy.sx");
   }
   else
   {
      // --- build default top-level group
      SxPtr<SxMpiTaskGroup> mpiTg = SxPtr<SxMpiTaskGroup>::create ();

      SxMpiComm mpiComm;
      SxMpiComm parentIntraComm;
      SxString parentId;

      int mpiSize = SxLoopMPI::nr ();
      SxString name = "top";
      int siblings = 1;

      mpiTg->setName(name);
      mpiTg->setNsiblings(siblings);
      mpiTg->setNmembers(mpiSize);

      parentId = "ROOT";
      // --- better set parent pointer of "top" to itself?
//      mpiTg->setParent(NULL);
      mpiTg->setParent( mpiTg.getPtr() );
      // ---
      if (mpiSize > 1) {
         mpiComm = parentIntraComm = SxMpi::mpiCommWorld();
         mpiTg->setIntraCommunicator(mpiComm);
      }
      mpiTg->setParentId(parentId);

      if (mpiSize > 1) {
         mpiComm = SxMpi::joinMpiComm(mpiTg->getMasterCommunicator(), parentIntraComm);
         mpiTg->setInterCommunicator(mpiComm);

         mpiTg->setSiblingRank( SxMpi::rank(parentIntraComm) / (SxMpi::size(parentIntraComm)/siblings) );
         mpiTg->setMemberRank( SxMpi::rank(mpiTg->getIntraCommunicator()) );
      } else {
         mpiTg->setSiblingRank (0);
         mpiTg->setMemberRank (0);
      }

      mpiTaskGroups(name) = mpiTg;

      // complete the setup
      build_hierarchy_pointers ();
   }
}

SxParallelHierarchy::SxParallelHierarchy(const SxString &hierarchyFile)
{
   if (initialized) return;
   buildHierarchy(hierarchyFile);
}

void SxParallelHierarchy::destroy()
{
   taskGroups.removeAll();
   mpiTaskGroups.removeAll();
   initialized = false;
   cout << SX_SEPARATOR;
   cout << "| parallel hierarchy cleared by ::destroy()" << endl;
   cout << SX_SEPARATOR;
}

void SxParallelHierarchy::buildHierarchy(const SxString &hierarchyFile)
{
   if (initialized) return;

   SxParser parser;
   SxParser::Table table = parser.read(hierarchyFile);
   SxSymbolTable * top = const_cast<SxSymbolTable*> (table.getPtr());

   // (0) brief sloppy consistency check of the symbol table
   try {
      if (top->containsGroup("level")) {
         SxSymbolTable * child = top->getGroup("level");
         if ((child->get("name")->toString() == SxString("top-level")) &&
               (child->get("members")->toInt() > 0) &&
               (child->get("members")->toInt() != SxLoopMPI::nr ())) {
            cout << SX_SEPARATOR;
            cout << "| Error: Parallel hierarchy is incompatible with actual"
                    " MPI resources" << endl;
            cout << "| Number of top-level members in hierarchy: "
                 << child->get("members")->toInt() << endl;
            cout << "| Number of MPI processes:                  "
                 << SxLoopMPI::nr () << endl;
            cout << "| Check hierarchy file '" << hierarchyFile
                 << "' or mpirun parameters." << endl;
            cout << SX_SEPARATOR;
            SX_QUIT;
         }
      }
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }

   // (1) recursively parse the symbol table and build the parallel hierarchy accordingly
   build_hierarchy_recursively (top, NULL);

   // (2) set-up pointers between the task groups
   build_hierarchy_pointers ();

   info ();
}


// TODO integrate build_hierarchy_pointers() into build_hierarchy_recursively()
void SxParallelHierarchy::build_hierarchy_pointers ()
{
   if (initialized)
      return;

   for (SxMap< SxString, SxPtr<SxMpiTaskGroup> >::Iterator it = mpiTaskGroups.begin(); it != mpiTaskGroups.end(); ++it)
   {
      taskGroups(it.getKey()) = it.getValue().getPtr();
   }

   // set up SxTaskGroup* connections between parents and children
   for (SxMap<SxString, SxTaskGroup*>::Iterator it = taskGroups.begin(); it != taskGroups.end(); ++it)
   {
      SxTaskGroup * parentTg;
      parentTg = it.getValue();
      SX_CHECK(parentTg);
      // autoLvl routines make the links themselves
      if (parentTg->hasAutoLvlName())
         continue;
      if (parentTg->childId.getSize () > 0)
      {
         for (SxList<SxString>::Iterator it2 = parentTg->childId.begin();
               it2 != parentTg->childId.end(); ++it2)
         {
            SxTaskGroup * childTg = taskGroups(*it2);
            SX_CHECK(childTg);
            parentTg->addChild(childTg);
            childTg->setParent(parentTg);
         }
         for (SxMap<SxString, SxTaskGroup*>::Iterator it2 = parentTg->children.begin(); it2 != parentTg->children.end(); ++it2)
         {
            SX_CHECK( it2.getValue() );
         }
      }
   }

   initialized = true;
   return;
}


void SxParallelHierarchy::info()
{
   cout << SX_SEPARATOR;
   cout << "| SxParallelHierarchy::info() ..." << endl;
   if (taskGroups.getSize() > 0)
   {
      for (SxMap<SxString, SxTaskGroup*>::Iterator it = taskGroups.begin(); it != taskGroups.end(); ++it)
      {
         SxString id = it.getKey();
         SxTaskGroup *tg  = it.getValue();
         SX_CHECK( (id == tg->getName()) || (id == tg->getAutoLvlName()) );
         sxprintf("| %s : members=%d, siblings=%d", tg->getName().ascii(), tg->getNmembers(), tg->getNsiblings());
         if (tg->hasChildren())
         {
            sxprintf("; children :");
            for (SxMap<SxString,SxTaskGroup*>::Iterator it2 = tg->children.begin(); it2 != tg->children.end(); ++it2)
            {
               SxString id2 = it2.getKey();
#ifndef NDEBUG
               SxTaskGroup *tg2 = it2.getValue();
               SX_CHECK( id2 == tg2->getName() );
#endif
               sxprintf(" %s", id2.ascii());
            }
         }
         else
         {
            sxprintf("; no children");
         }
         sxprintf("\n");
      }
   }
   else
   {
      sxprintf("| no entries\n");
   }
   cout << SX_SEPARATOR;
}


SxTaskGroup * SxParallelHierarchy::getTaskGroup(const SxString &name)
{
   if (taskGroups.containsKey(name))
      return taskGroups(name);
   else
   {
      sxprintf("\n%s", SX_SEPARATOR.ascii());
      sxprintf("| In %s :\n", SX_FUNC);
      sxprintf("|    SxTaskGroup \"%s\" is not defined.\n", name.ascii());
      sxprintf("%s", SX_SEPARATOR.ascii());
      SX_EXIT;
      return NULL;
   }
}


void SxParallelHierarchy::build_hierarchy_recursively(SxSymbolTable * node, SxSymbolTable * parent)
{
   if (node == NULL)
      return;

   SxString name = "";
   int siblings = 0;
   int members = 0;

   if (node->contains("name", true)) {
      name = node->get("name")->toString();
      siblings = node->get("siblings")->toInt();
      members = node->get("members")->toInt();

      // catch default values, see /src/share/sphinx/parallelHierarchy.sx
      if (siblings < 1)
         siblings = SxLoopMPI::nr ();
      if (members < 1)
         members = SxLoopMPI::nr ();

      SxPtr<SxMpiTaskGroup> mpiTg = SxPtr<SxMpiTaskGroup>::create ();

      SxMpiComm mpiComm;
      SxMpiComm parentIntraComm;
      SxString parentId;

      mpiTg->setName(name);
      mpiTg->setNsiblings(siblings);
      mpiTg->setNmembers(members);

      if (parent->contains("name", true)) {
         SxMpiTaskGroup * parentTg;
         parentId = parent->get("name")->toString();
         parentTg = mpiTaskGroups(parentId).getPtr();
         parentIntraComm = parentTg->getIntraCommunicator();
         mpiComm = SxMpi::divideMpiComm(parentIntraComm, siblings);
      }
      else {
         // --- better set parent pointer of "top" to itself?
         parentId = "ROOT";
         //mpiTg->setParent(NULL);
         // ---
         mpiTg->setParent( mpiTg.getPtr() );
         mpiComm = parentIntraComm = SxMpi::mpiCommWorld();
      }
      mpiTg->setIntraCommunicator(mpiComm);
      mpiTg->setParentId(parentId);

      mpiComm = SxMpi::joinMpiComm(mpiTg->getMasterCommunicator(), parentIntraComm);
      mpiTg->setInterCommunicator(mpiComm);

      mpiTg->setSiblingRank( SxMpi::rank(parentIntraComm) / (SxMpi::size(parentIntraComm)/siblings) );
      mpiTg->setMemberRank( SxMpi::rank(mpiTg->getIntraCommunicator()) );

      // work with strings for now and build the parent-child pointer maps later
      // as the memory addresses inside the map are not yet fixed
      // TODO eliminate childIds
      for (SxList<SxSymbolTable*>::ConstIterator child = node->children.begin();
            child != node->children.end(); child++) {
         if ((*child)->name == "level")
            mpiTg->addChildId( (*child)->get("name")->toString() );
      }

      mpiTaskGroups(name) = mpiTg;
   }

   for (SxList<SxSymbolTable*>::ConstIterator child = node->children.begin();
         child != node->children.end(); child++) {
      if ((*child)->name == "level")
         build_hierarchy_recursively(*child, node);
   }

   return;
}


SxMpiTaskGroup * SxParallelHierarchy::getMpiTaskGroup(const SxString &name)
{
   SX_CHECK( mpiTaskGroups.containsKey(name) );
   return mpiTaskGroups(name).getPtr ();
}


void SxParallelHierarchy::write (const SxString &filename)
{
   if (!initialized)
      return;
   // find the root node
   SxTaskGroup* root = taskGroups.begin ().getValue ();
   SX_CHECK (root);
   while (root->getParentId () != "ROOT") root = root->getParent ();
   if (!root->master ()) return; // only the master node must write this

   // --- printout
   ofstream out(filename.ascii ());
   out << "format parallelHierarchy;" << endl;
   root->write (out);
   out.close ();
}
