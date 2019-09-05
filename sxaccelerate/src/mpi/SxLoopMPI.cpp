#include <SxLoopMPI.h>
#include <SxString.h>
#include <stdio.h>
#include <SxError.h>

#include <SxParallelHierarchy.h>
#include <SxTaskGroup.h>
#include <SxStack.h>
#include <SxNArray.h>

#ifdef USE_LOOPMPI
#include <mpi.h>
#include <SxCLI.h>
#include <SxFileIO.h>
#endif // USE_LOOPMPI

void SxLoopMPICrash ()
{
#ifdef USE_LOOPMPI
   MPI_Abort (MPI_COMM_WORLD, 5);
#endif
}


bool SxLoopMPI::doParallelFft=false;
bool SxLoopMPI::initialized=false;
SxStack<SxTaskGroup*> SxLoopMPI::mpiLevel;

SxLoopMPI::TaskSelection SxLoopMPI::MPIsource;
SxLoopMPI::TaskSelection SxLoopMPI::MPItarget;

// reader function of SxParser
extern SxArray<char> (*SxParser_readFile)(const SxString &, int64_t);

#ifdef USE_LOOPMPI
SxArray<char> readFile_MPI(const SxString &filename, int64_t nBytes)
{
   SxArray<char> res;
   SX_MPI_MASTER_ONLY  {
      // read file in master task
      res = SxFileIO::readBinary (filename, nBytes);
      if (res.getSize () != int(res.getSize ()))  {
         cout << "Too large file '" << filename << "' for MPI broadcast: "
              << res.getSize () << " bytes" << endl;
         SxLoopMPICrash ();
         SX_EXIT;
      }
   }
   SX_MPI_SOURCE(TopLevel);
   SX_MPI_TARGET(TopLevel, TaskGroupAll);
   // broadcast content
   int size = SxLoopMPI::bcast ((int)res.getSize (), 0);
   if (res.getSize () == 0) res.resize (size);
   SxLoopMPI::bcast (res.elements, size, 0);
   return res;
}
#endif

#ifdef USE_LOOPMPI
void SxLoopMPI::init (int argc, char **argv)
{
#  ifdef USE_MPI
      MPI_Init (&argc, &argv);
      SxCLI::atExit.append (&SxLoopMPI::finalize); // to be run if SxCLI exits
      SxParser_readFile = readFile_MPI;
#  else
      SX_UNUSED (argc);
      SX_UNUSED (argv);
#  endif
   initialized = true;
   if (me ()>0)  // redirecting stdout if we are not master
   {
      SxCLI::Log::disable ("not MPI master");
      const char *sxlogstdout = ::getenv ("SX_LOG_STDOUT");
      if (sxlogstdout && SxString(sxlogstdout) == "ALL")  {
         SxString file = ".stdoutDump_"+SxString(me ());
         freopen (file.ascii(), "w", stdout);
      } else {
         freopen ("/dev/null", "w", stdout);
      }
   }
   nr ();
   sxprintf ("MPI parallel run: %d/%d\n", me (), nr ());
//   sxExitCB = SxLoopMPICrash;
//   sxQuitCB = SxLoopMPICrash;
}
#else
void SxLoopMPI::init (int, char **) {}
#endif // USE_LOOPMPI

void SxLoopMPI::finalize ()
{
#ifdef USE_LOOPMPI
   if (initialized)  {
      // shut down MPI
      (void)MPI_Finalize ();
      SxParser_readFile = SxFileIO::readBinary;
   }
   initialized = false;
#endif // USE_LOOPMPI
}
// automatically call SxLoopMPI::finalize
SxLoopMPI::ShutDown SxLoopMPI::shutDown;


int SxLoopMPI::me ()
{
#ifdef USE_LOOPMPI
   static int myId=-1;
   if (!initialized) {
      sxprintf ("Warning: MPI init has not been called!\n");
      return 0;
   }
   if (myId<0) // call MPI_Comm_rank only once
   {
      MPI_Comm_rank(MPI_COMM_WORLD, &myId);
   }
   return myId;
#else
   return 0;
#endif // USE_LOOPMPI
}

int SxLoopMPI::nr ()
{
#ifdef USE_LOOPMPI
   static int nrProc=-1;
   if (!initialized) {
      sxprintf ("Warning: MPI init has not been called!\n");
      return 1;
   }
   if (nrProc<0) // call MPI_Comm_size only once
   {
      MPI_Comm_size(MPI_COMM_WORLD, &nrProc);
   }
   return nrProc;
#else
   return 1;
#endif // USE_LOOPMPI
}

#ifdef USE_LOOPMPI
bool SxLoopMPI::myWork (int i)
{
   if (mpiLevel.isEmpty())
      return (i%nr ()) == me ();
   else
   {
      return (mpiLevel.top())->myWork(i);
   }
}
#else
bool SxLoopMPI::myWork (int)
{
   return true;
}
#endif // USE_LOOPMPI

#ifdef USE_LOOPMPI
int SxLoopMPI::whoseWork (int i)
{
   if (mpiLevel.isEmpty())
      return (i%nr ());
   else
      return mpiLevel.top ()->whoseWork (i);
}
#else
int SxLoopMPI::whoseWork (int)
{
   return 0;
}
#endif // USE_LOOPMPI

int SxLoopMPI::invMap ()
{
   return me ();
}

bool SxLoopMPI::lor (bool in)
{
#ifdef USE_LOOPMPI
   if (!initialized) {
      sxprintf ("Warning: MPI init has not been called!\n");
      return in;
   }
   int send, recv;
   send=in;
   MPI_Allreduce(&send,&recv,1,MPI_INT,MPI_LOR,MPI_COMM_WORLD);
   return recv;
#else
   return in;
#endif // USE_LOOPMPI
}

#ifdef USE_LOOPMPI
int SxLoopMPI::sum (int in)
{
   int res = 0;
   sum (&in, &res, 1);
   return res;
}

double SxLoopMPI::sum (double in)
{
   double res = 0.;
   sum (&in, &res, 1);
   return res;
}

double SxLoopMPI::bcast (double in, int source)
{
   double msg=in;
   SxLoopMPI::bcast(&msg,1,source);
   return msg;
}

int SxLoopMPI::bcast (int in, int source)
{
   int msg=in;
   SxLoopMPI::bcast(&msg,1,source);
   return msg;
}
#endif // USE_LOOPMPI

#ifdef USE_LOOPMPI
void SxLoopMPI::barrier ()
{
   if (!initialized) {
      sxprintf ("Warning: MPI init has not been called!\n");
      return;
   }
   MPI_Barrier(MPI_COMM_WORLD);
}
#endif // USE_LOOPMPI


#ifdef USE_LOOPMPI
void SxLoopMPI::setParallelFft (bool mode)
{
   SxLoopMPI::doParallelFft=mode;
   SxLoopMPI::barrier ();
}
#else
void SxLoopMPI::setParallelFft (bool) {}
#endif // USE_LOOPMPI

bool SxLoopMPI::parallelFft ()
{
#ifdef USE_LOOPMPI
   return SxLoopMPI::doParallelFft;
#else
   return false;
#endif // USE_LOOPMPI
}


#ifdef USE_MPI
MPI_Comm SxLoopMPI::MpiCommWorld()
{
   return MPI_COMM_WORLD;
}

MPI_Info SxLoopMPI::MpiInfoNull()
{
   return MPI_INFO_NULL;
}
#endif



void SxLoopMPI::pushMpiLevel(SxTaskGroup * level)
{
//   sxprintf("%s/n", __PRETTY_FUNCTION__);
   SX_CHECK( level );
	mpiLevel.push(level);
}

SxTaskGroup * SxLoopMPI::popMpiLevel()
{
//   sxprintf("%s/n", __PRETTY_FUNCTION__);
	return mpiLevel.pop();
}

static SxString autoLvlGlobal (SxParallelHierarchy &ph, const SxString &context, const SxString &level)
{
   SxString globalName = context + "." + level;
   if (ph.taskGroups.containsKey (globalName))
   {
      cout << SX_SEPARATOR;
      cout << "| SxTaskGroup global name collision detected for " << globalName << "." << endl;
      {
         int id;
         for (id = 1; id < 100; ++id)  {
            if (!ph.mpiTaskGroups.containsKey (globalName + "-" + id))
               break;
         }
         if (id >= 100)  {
            // We arrive here if this level has been used in 100
            // different contexts (=outer MPI levels).
            // Is this sensible?
            cout << "| ERROR: auto-level explosion for " << level << endl;
            cout << "| MPI hierarchy needs to be checked" << endl;
            SX_EXIT;
         }
         globalName = globalName + "." + id;
      }
      cout << "| Using the modified identifier " << globalName << " instead." << endl;
      cout << SX_SEPARATOR;
   }
   return ( globalName );
}


SxMPILevel::SxMPILevel(const SxString &level)
{
   SxParallelHierarchy ph;
   SxTaskGroup *current = NULL;
   SxString context = "SxMpiLevel";

   if (SxLoopMPI::mpiLevel.isEmpty ())  {
      // --- try via global MPI level table
      if (ph.taskGroups.containsKey (level))  {
         SxLoopMPI::pushMpiLevel( ph.getTaskGroup(level) );
         return;
      }
      // otherwise, get the top level of the parallel hierarchy
      if (!ph.taskGroups.containsKey ("top"))  {
         cout << SX_SEPARATOR;
         cout << "| Error: Incomplete parallel hierarchy." << endl;
         cout << "| Automatic hierarchy requires a top level named \"top\"\n"
              << "| and levels named \"all\" within each incomplete level."
              << endl;
         cout << SX_SEPARATOR;
         SX_EXIT;
      }
      current = ph.getTaskGroup ("top");
   } else {
      current = SxLoopMPI::mpiLevel.top ();
   }

   if (current->getName () == level)  {
      // we are already in the right level. Push it again
      SxLoopMPI::pushMpiLevel (current);
      return;
   }

   // --- create requested level automatically if necessary
   if (!current->children.containsKey (level))  {
      cout << SX_SEPARATOR;
      cout << "| Current MPI level: " << current->getName () << endl;
//      cout << "| Children: " << current->children.getKeys () << endl;
      cout << "| WARNING: automatically introduced MPI level '" << level << "'";
      if (ph.taskGroups.containsKey(current->getName ()))  {
           cout << " below level '" << current->getName () << "' ";
      } 
      cout << endl << "|    (no MPI parallelization)" << endl;
      cout << SX_SEPARATOR;

      SX_CHECK (dynamic_cast<SxMpiTaskGroup*>(current));
      SxMpiComm parentIntraComm = dynamic_cast<SxMpiTaskGroup*>(current) ->getIntraCommunicator ();
      int siblings = 1;

      SxPtr<SxMpiTaskGroup> autoGroup = SxPtr<SxMpiTaskGroup>::create ();

      autoGroup->setNmembers (current->getNmembers ());
      autoGroup->setNsiblings (siblings);

      autoGroup->setSiblingRank(0);
      if (current->getNmembers () > 1)  {
         autoGroup->setIntraCommunicator (SxMpi::divideMpiComm(parentIntraComm, siblings));
         // --- masterCommunicator is set by the previous line -- the call order matters!
         autoGroup->setInterCommunicator (SxMpi::joinMpiComm(autoGroup->getMasterCommunicator(), parentIntraComm));

         autoGroup->setMemberRank( SxMpi::rank(autoGroup->getIntraCommunicator()) );
      } else {
         autoGroup->setMemberRank (0);
      }

      autoGroup->setName(level);

      SxString globalName = autoLvlGlobal (ph, context, level);
      autoGroup->setAutoLvlName(globalName);

      ph.mpiTaskGroups(globalName) = autoGroup;
      ph.taskGroups(globalName) = ph.mpiTaskGroups(globalName).getPtr();

      SxTaskGroup *autoGroupPtr = ph.taskGroups(globalName);
      autoGroupPtr->setParent(current);
      autoGroupPtr->setParentId(current->getName());
      current->addChild(autoGroupPtr);
      current->addChildId(level);
   }
   SxLoopMPI::pushMpiLevel (current->children(level));
}


SxMPILevel::~SxMPILevel()
{
   if (!SxLoopMPI::mpiLevel.isEmpty ())
      SxLoopMPI::popMpiLevel();
}

void SxAutoLevel::finalize ()
{
#ifdef USE_LOOPMPI
   SxParallelHierarchy ph;
   SxTaskGroup* current = SxLoopMPI::mpiLevel.isEmpty () 
                        ? ph.taskGroups("top") 
                        : SxLoopMPI::mpiLevel.top ();
   int iStart;
   for (iStart = 0;  iStart < levels.getSize (); iStart++)  {
      if (current->children.containsKey (levels(iStart).name))
         current = current->children (levels(iStart).name);
      else
         break;
   }
   if (iStart == levels.getSize ()) return; // done

   SxArray<int> siblings (levels.getSize ());

   // current task group's master decides on the parallelization pattern within the group
   if (current->master ())  {
      int nMember = current->getNmembers ();
      cout << SX_SEPARATOR;
      cout << "| automatically creating MPI levels" << endl;
      cout << "|   available MPI tasks within current level (" << current->getName ()
           << "): " << nMember << endl;
      for (int iLevel = iStart; iLevel < levels.getSize (); ++iLevel)  {
         cout << "|   level " << levels(iLevel).name << " requests " << levels(iLevel).work
              << " items" << endl;
      }
      // --- work space
      SxArray2<int> optSiblings((int)levels.getSize (), nMember);
      SxArray2<double> optEff((int)levels.getSize (), nMember);
      optSiblings.set (-1);

      // --- now optimize the efficiency
      getEfficiency (nMember, iStart, optSiblings, optEff);

      // --- printout
      cout << "| final decision" << endl;
      for (int iLevel = iStart; iLevel < levels.getSize (); iLevel++)  {
         siblings(iLevel) = optSiblings (iLevel, nMember - 1);
         cout << "|   " << levels(iLevel).name << " with " << siblings(iLevel)
              << " siblings ("
              << 0.01 * round(optEff(iLevel, nMember - 1) * 1e4) << "%)" << endl;
         nMember /= siblings(iLevel);
      }
   }

   // parent for the auto-level groups
   SxMpiTaskGroup *parent = dynamic_cast<SxMpiTaskGroup*>(current);
   if (parent->getNmembers () > 1)  {
      SxMpi::broadcast (siblings.elements, (int)siblings.getSize (), 
                        parent->getIntraCommunicator ());
   }

   // --- now create the Mpi levels
   for (int iLevel = iStart; iLevel < levels.getSize (); ++iLevel)  {
      SX_CHECK (parent);

      SxPtr<SxMpiTaskGroup> mpiTg = SxPtr<SxMpiTaskGroup>::create ();

      int nMembers  = parent->getNmembers () / siblings(iLevel),
          nSiblings = siblings(iLevel);
      if (iLevel + 1 == levels.getSize ())  {
         // final level: no siblings, all become workers
         nMembers = parent->getNmembers ();
         nSiblings = 1;
      }
      mpiTg->setNsiblings (nSiblings);
      mpiTg->setNmembers  (nMembers);

      // --- communicators/ranks
      SxMpiComm parentIntraComm = parent->getIntraCommunicator ();
      if (parent->getNmembers () > 1)  {
         mpiTg->setIntraCommunicator (SxMpi::divideMpiComm(parentIntraComm, nSiblings));
         mpiTg->setInterCommunicator (SxMpi::joinMpiComm(mpiTg->getMasterCommunicator(), parentIntraComm));
         mpiTg->setSiblingRank (SxMpi::rank(parentIntraComm) / (SxMpi::size(parentIntraComm)/nSiblings) );
         mpiTg->setMemberRank (SxMpi::rank(mpiTg->getIntraCommunicator()) );
      } else {
         mpiTg->setSiblingRank (0);
         mpiTg->setMemberRank (0);
      }

      mpiTg->setParent (parent);
      mpiTg->setParentId (parent->getName());
      mpiTg->setName (levels(iLevel).name);

      SxString globalName = autoLvlGlobal (ph, context, levels(iLevel).name);
      mpiTg->setAutoLvlName(globalName);

      ph.mpiTaskGroups (globalName) = mpiTg;
      SxMpiTaskGroup *tg = ph.mpiTaskGroups(globalName).getPtr ();

      ph.taskGroups (globalName) = tg;

      parent->addChildId (levels(iLevel).name);
      parent->children (levels(iLevel).name) = tg;

      cout << "| created " << levels(iLevel).name
           << "/" << globalName
           << " at " << parent->getName ()
           << " (" << nSiblings << " siblings, "
           << tg->getNmembers () << " members each)" << endl;

      parent = tg;
   }

   cout << SX_SEPARATOR;
#endif
}

double SxAutoLevel::getEfficiency (int nMember, int iLevel,
                            SxArray2<int> &optSiblings,
                            SxArray2<double> &optEff)
{
   // no more levels to parallelize?  ==> 1/N efficiency
   if (iLevel >= levels.getSize ()) {
      return 1./double(nMember);
   }

   if (optSiblings(iLevel, nMember - 1) > 0) {
      return optEff(iLevel, nMember - 1);
   }

   // no parallelization possible?  ==> 100% efficiency
   if (nMember == 1)  {
      for ( ; iLevel < levels.getSize (); ++iLevel)  {
         optSiblings(iLevel,0) = 1;
         optEff(iLevel,0) = 1.;
      }
      return 1.;
   }

   Workload &wl = levels(iLevel);
   int nSiblings;
   double workPerSibling, optEfficiency = -1.;
   nSiblings = nMember;
   SX_CHECK (wl.work >= 1);
   if (wl.work < nSiblings) nSiblings = wl.work;
   SX_CHECK (nSiblings >= 1);
   for ( ; nSiblings > 0; nSiblings --)  {
      // find divisor
      while (nMember % nSiblings > 0) nSiblings--;
      if (wl.work % nSiblings == 0)
         workPerSibling = wl.work / nSiblings;
      else
         workPerSibling = ceil (double(wl.work) / nSiblings);
      cout << "|   " << wl.name << " => consider " << nSiblings << " siblings with each "
           << (nMember / nSiblings) << " members" << endl;
      double time = wl.serial * wl.work 
                  + (1. - wl.serial) * workPerSibling
                    / getEfficiency (nMember / nSiblings, iLevel + 1, optSiblings, optEff);
      double efficiency = wl.work / (time * nSiblings);
      if (optEfficiency < efficiency)  {
         optSiblings(iLevel, nMember - 1) = nSiblings;
         optEfficiency = efficiency;
      }
   }
   optEff(iLevel, nMember - 1) = optEfficiency;

   cout << "|   efficiency for " // << level " << wl.name << " with "
        << nMember << " MPI processes: "
        << round (optEfficiency*1e4)*0.01 << "% with " << optSiblings(iLevel, nMember - 1)
        << " siblings" << endl;
   
   return optEfficiency;
}

#ifdef USE_MPI
SxMpiTaskGroup* SxLoopMPI::getLevel(const TaskSelection &lvl)
{
   // --- top level
   if (lvl.levelType == TopLevel)
      return dynamic_cast<SxMpiTaskGroup*>
            (SxParallelHierarchy::taskGroups ("top"));

   // --- current level
   if (lvl.levelType == CurrentLevel)  {
      if (SxLoopMPI::mpiLevel.isEmpty ())
         return dynamic_cast<SxMpiTaskGroup*>
               (SxParallelHierarchy::taskGroups ("top"));
      else
         return getLevel ();
   }

   // --- get the current MPI level and search by name
   SX_CHECK (lvl.levelType == SpecifiedName);
   SxMpiTaskGroup *current = SxLoopMPI::mpiLevel.isEmpty ()
                           ? dynamic_cast<SxMpiTaskGroup*>
                             (SxParallelHierarchy::taskGroups ("top"))
                           : SxLoopMPI::getLevel ();
   SX_CHECK(current);

   // name matches current level ?
   if (current->getName () == lvl.level)
      return current;

   SxString lvlName = lvl.level;
   // try children of current level
   if (current->children.containsKey (lvlName))
      return dynamic_cast<SxMpiTaskGroup*>(current->children(lvlName));
   // try parents of current level
   for ( ; current->getParentId () != "ROOT" ; 
         current = dynamic_cast<SxMpiTaskGroup*>(current->getParent ()))
   {
      if (current->getName () == lvlName) return current;
   }

   // --- program should never arrive here
   SX_CHECK(false);
   return NULL;
}

SxMpiTaskGroup* SxLoopMPI::getTarget(SxTaskGroup *source)
{
   SX_CHECK (source);
   if (source->getName () == MPItarget.level) 
      return dynamic_cast<SxMpiTaskGroup*>(source);

   // search the source's parents
   SxMpiTaskGroup *current = NULL;
   if (MPItarget.levelType == CurrentLevel)  {
      if (mpiLevel.isEmpty ())
         current = dynamic_cast<SxMpiTaskGroup*>
                   (SxParallelHierarchy::taskGroups("top"));
      else
         current = getLevel ();
   }

   while (true) {
      if (source->getParentId () == MPItarget.level)  {
         return dynamic_cast<SxMpiTaskGroup*>(source->getParent ());
      }
      if (MPItarget.levelType == CurrentLevel && source == current)  {
         return current;
      }
      if (source->getParentId () == "ROOT")  {
         // top of the level tree. Is this what we wanted?
         if (MPItarget.levelType == TopLevel)
            return dynamic_cast<SxMpiTaskGroup*>(source);
         // else ...
         cout << "Invalid source/target relation for MPI sum" << endl;
         cout << "Source level is " << source->getName () << endl;
         cout << "Target level is ";
         if (MPItarget.levelType == SpecifiedName)
            cout << MPItarget.level;
         else
            cout << "current level (" << current->getName () << ')';
         cout << endl << "Source must be sublevel of target" << endl;
         SX_EXIT;
      }
      source = source->getParent ();
   }

   // --- program should never arrive here
   SX_CHECK(false);
   return NULL;
}


template<class T>
void SxLoopMPI::sum(T* in, T* out, ssize_t nElem)
{
   if (!initialized)
      sxprintf ("Warning: MPI init has not been called!\n");

   if (!initialized || nr () == 1)  {
      if (in != out) memcpy (out, in, nElem * sizeof(T));
      return;
   }

   MPI_Datatype dt = getMpiDataType (in);

   // nElem must not exceed the value limit of int
   SX_CHECK_VARS (nElem < (1<<30), nElem, (1<<30));
   int N = int(nElem);
   void *inPtr = in;

   // --- global sum
   if (MPIsource.levelType == TopLevel)  {
      SX_CHECK (MPItarget.levelType == TopLevel);
      SX_CHECK (MPIsource.who == TaskGroupAll);
      if (MPItarget.who == TaskGroupAll)  {
         if (in == out) inPtr = MPI_IN_PLACE;
         MPI_Allreduce(inPtr, out, N, dt, MPI_SUM, MPI_COMM_WORLD);
      } else {
         if (in == out && me () == 0) inPtr = MPI_IN_PLACE;
         MPI_Reduce(inPtr, out, N, dt, MPI_SUM, 0, MPI_COMM_WORLD);
      }
      return;
   }

   if (MPIsource.who == TaskGroupAll)  {
      // --- sum via the target all-rank communicator
      // the actual source level is irrelevant
      SxMpiTaskGroup *target = getLevel (MPItarget);
      if (!target)  {
         // try to get the target via the source level
         SxTaskGroup *source = getLevel (MPIsource);
         if (source)  {
            target = getTarget(source);
         } else {
            // try direct access via parallel hierarchy
            if (SxParallelHierarchy::taskGroups.containsKey(MPItarget.level))
               target = dynamic_cast<SxMpiTaskGroup*>
                        (SxParallelHierarchy::taskGroups(MPItarget.level));
            else
               SX_EXIT;
         }
      }
      // --- do the actual summation
      MPI_Comm comm = target->getIntraCommunicator ().comm;
      if (MPItarget.who == TaskGroupAll)  {
         if (in == out) inPtr = MPI_IN_PLACE;
         MPI_Allreduce (inPtr, out, N, dt, MPI_SUM, comm);
      } else {
         if (in == out && target->master ()) inPtr = MPI_IN_PLACE;
         MPI_Reduce (inPtr, out, N, dt, MPI_SUM, 0, comm);
      }
   } else {
      SX_CHECK (MPIsource.who == TaskGroupMaster);
      // --- sum recursively to the master of each task group
      //     via the inter-master communicator of its sublevel
      
      // --- get source and target (both are needed)
      SxMpiTaskGroup *source = getLevel (MPIsource);
      SxMpiTaskGroup *target = NULL;
      if (!source)  {
         target = getLevel (MPItarget);
         if (!target) { SX_EXIT; }
         if (target->children.containsKey (MPIsource.level))
            source = dynamic_cast<SxMpiTaskGroup*>
                     (target->getChild (MPIsource.level));
         else
            SX_EXIT;
      } else {
         target = getTarget (source);
      }
      SX_CHECK (source != target);

      // --- recursive reduction onto master
      for (; source != target; 
           source = dynamic_cast<SxMpiTaskGroup*>(source->getParent ()))
      {
         MPI_Comm comm;
         if (!source->hasChildren ())  {
            SX_CHECK_VAR (source->getNsiblings () == 1, source->getNsiblings ());
            if (source->getNmembers () == 1) continue; // nothing to sum
            // bottom level => each rank would be a master
            comm = source->getIntraCommunicator ().comm;
            if (inPtr == out && source->master ())
               inPtr = MPI_IN_PLACE;
         } else {
            if (source->getNsiblings () == 1) continue; // nothing to do
            if (!source->master ()) break;
            comm = source->getInterCommunicator ().comm;
            if (inPtr == out && dynamic_cast<SxMpiTaskGroup*>
                                (source->getParent ())->master () )
            {
               inPtr = MPI_IN_PLACE;
            }
         }

         MPI_Reduce(inPtr, out, N, dt, MPI_SUM, 0, comm);
         inPtr = out;
      }
      if (target->master () && inPtr != out)  {
         // no summation occured, need to copy from in to out
         memcpy (out, in, nElem * sizeof(T));
      }
      // now result is on target master. If necessary, broadcast result
      if (MPItarget.who == TaskGroupAll)
         MPI_Bcast (out, N, dt, 0, target->getIntraCommunicator ().comm);
   }
}

// --- instantiate explicitly the relevant cases
template void SxLoopMPI::sum<int> (int*, int*, ssize_t);
template void SxLoopMPI::sum<double> (double*, double*, ssize_t);



template<class T>
void SxLoopMPI::bcast(T* inout, ssize_t nElem, int sourceId)
{
   if (!initialized)
      sxprintf ("Warning: MPI init has not been called!\n");

   if (!initialized || nr () == 1)
      return;

   MPI_Datatype dt = getMpiDataType (inout);
   // nElem must not exceed the value limit of int
   SX_CHECK_VARS (nElem < (1<<30), nElem, (1<<30));
   int N = int(nElem);

   if (MPItarget.levelType == TopLevel && MPIsource.levelType == TopLevel) {
      // ignore settings in MPI{source,target}.who
      MPI_Bcast (inout, N, dt, sourceId, MPI_COMM_WORLD);
      return;
   }

   SX_CHECK(MPIsource.who == TaskGroupMaster);
   SX_CHECK(   MPItarget.levelType == CurrentLevel 
            || MPItarget.levelType == LevelSiblings);

   SxMpiTaskGroup* source = getLevel (MPIsource);
   if (!source)  {
      cout << "Failed to find source level '" << MPIsource.level
           << "'for MPI broadcast" << endl;
      cout << "The broadcasting source level should be close to the"
           << "current SX_MPI_LEVEL" << endl;
      SX_EXIT;
   }

   MPI_Comm comm;
   if (MPItarget.levelType == LevelSiblings)  {
      if (!source->hasChildren ())  {
         SX_CHECK (source->getNsiblings () == 1);
         // bottom level
         if (source->getNmembers () == 1) return; // nothing to do
         comm = source->getIntraCommunicator ().comm;
      } else {
         if (MPItarget.who == TaskGroupAll)  {
            comm = dynamic_cast<SxMpiTaskGroup*>(source->getParent ())
                   ->getIntraCommunicator ().comm;
            sourceId *= source->getNmembers (); // only masters send data
         } else {
            SX_CHECK (MPItarget.who == TaskGroupMaster);
            if (source->getNsiblings () == 1) return; // nothing to do
            comm = source->getInterCommunicator ().comm;
            if (!source->master ()) {
#ifndef NDEBUG
               if (dt == MPI_DOUBLE)  {
                  // invalidate vector on non-master
                  T nan = T(sqrt(-1.));
                  for (ssize_t i = 0; i < nElem; ++i) inout[i] = nan;
               }
#endif
               return; // nothing to broadcast or receive on non-master
            }
         }
      }
   } else {
      SX_CHECK (MPItarget.levelType == CurrentLevel);
      SX_CHECK (MPItarget.who == TaskGroupAll);
      SX_CHECK_VAR (sourceId == 0, sourceId); // only master to others allowed
      if (source->getNmembers () == 1) return; // nothing to do
      comm = source->getIntraCommunicator ().comm;
   }
   MPI_Bcast (inout, N, dt, sourceId, comm);
}

// --- instantiate explicitly the relevant cases
template void SxLoopMPI::bcast<char> (char*, ssize_t, int);
template void SxLoopMPI::bcast<int> (int*, ssize_t, int);
template void SxLoopMPI::bcast<double> (double*, ssize_t, int);
#endif
