#include <SxMpiComm.h>
#include <SxMpiTaskGroup.h>
#include <SxString.h>
#include <SxError.h>


SxMpiTaskGroup::SxMpiTaskGroup()
{
   intraComm = masterComm = interComm = SxMpi::mpiCommNull();
}

SxMpiTaskGroup::SxMpiTaskGroup(const SxMpiTaskGroup &tg)
   : SxTaskGroup(tg)
{
   if (verboseConstructors)
      sxprintf("%s called for %s\n", SX_FUNC, tg.name.ascii());
   intraComm = tg.intraComm;
   masterComm = tg.masterComm;
   interComm = tg.interComm;
}

SxMpiTaskGroup &SxMpiTaskGroup::operator=(const SxMpiTaskGroup &tg)
{
   if (verboseConstructors)
      sxprintf("%s called for %s\n", SX_FUNC, tg.name.ascii());
   SxTaskGroup::operator=(tg);
   intraComm = tg.intraComm;
   masterComm = tg.masterComm;
   interComm = tg.interComm;
   return *this;
}


void SxMpiTaskGroup::setIntraCommunicator(const SxMpiComm comm)
{
   // set the communicator which include all ranks of the taskGroup
   intraComm = comm;

   // set the communicator which includes only rank 0 of the group
   int nRanks;
   int ranks[1];
   nRanks = 1;
   ranks[0] = 0;
   masterComm = SxMpi::extractMpiComm(intraComm, ranks, nRanks);

   return;
}

void SxMpiTaskGroup::setInterCommunicator(const SxMpiComm comm)
{
   // set the communicator involving all the master ranks from all instances
   interComm = comm;
   return;
}



SxMpiComm SxMpiTaskGroup::getIntraCommunicator()
{
   return intraComm;
}

SxMpiComm SxMpiTaskGroup::getInterCommunicator()
{
   return interComm;
}

SxMpiComm SxMpiTaskGroup::getMasterCommunicator()
{
   return masterComm;
}


/*
int SxMpiTaskGroup::size()
{
   return SxMpi::size(intraComm);
}

int SxMpiTaskGroup::rank()
{
   return SxMpi::rank(intraComm);
}
*/

bool SxMpiTaskGroup::master()
{
   return (getMemberRank () == 0);
}


void SxMpiTaskGroup::printInfo()
{
   static bool firstLog = true;
   SxString logFile = SxString(".ph-")
         + SxString(SxMpi::size (), 4) + SxString(".log");
   FILE *fp;
   if (firstLog)
   {
      fp = fopen(logFile.ascii(), "w");
      firstLog = false;
   }
   else
      fp = fopen(logFile.ascii(), "a");
   //
   // --- check integer and string variables, first ---
   sxfprintf(fp, "%s\n", name.ascii());
   sxfprintf(fp, "   nMembers=%d memberRank=%d nSiblings=%d siblingRank=%d\n",
                     nMembers,   memberRank,   nSiblings,   siblingRank);
   sxfprintf(fp, "   parentId : %s\n", parentId.ascii());
   sxfprintf(fp, "   childId  :");
   for (SxList<SxString>::Iterator it = childId.begin(); it != childId.end(); ++it)
      sxfprintf(fp, " %s", (*it).ascii());
   sxfprintf(fp, "\n");
   // --- follow pointers to parent and child task groups
   sxfprintf(fp, "   parent   : %s\n", parent ? parent->getName().ascii() : "NULL");
   sxfprintf(fp, "   children :");
   for (SxMap<SxString, SxTaskGroup*>::Iterator it = children.begin(); it != children.end(); ++it)
      sxfprintf(fp, " %s|%s", it.getKey().ascii(), it.getValue()->getName().ascii());
   sxfprintf(fp, "\n");
   //
   fclose(fp);
}

void SxMpiTaskGroup::barrier()
{
#ifdef USE_MPI
//   sxprintf("===> %s\n", SX_FUNC);
   MPI_Barrier(intraComm.comm);
#endif
}




void SxMpiTaskGroup::sum(double *d, const ssize_t &n)
{
   sum(d, (int)n, INTRACOMM, SxMpi::mpiCommNull());
}

void SxMpiTaskGroup::sum(double *d, const ssize_t &n, const SxString &_childId)
{
   SxMpiTaskGroup * tgChild;
   tgChild = dynamic_cast<SxMpiTaskGroup*> (this->getChild(_childId));
   sum(d, (int)n, CHILDCOMM, tgChild->getInterCommunicator());
}


#ifdef USE_MPI
void SxMpiTaskGroup::sum(double *d, const int &n, const CommType &flag, const SxMpiComm &childComm)
{
   if (SxMpi::size () == 1) return;
   SxMpiComm tmpComm = SxMpi::mpiCommNull();
   double *buf = NULL;

   // (0) select communicator, allocate and initialize buffers
   if      (flag == INTRACOMM) {
      // objective: sum buffer d within the task group
      tmpComm = intraComm;
      //
      buf = new double[n];
      for (int i=0; i<n; ++i) buf[i] = 0.0;
   }
   else if (flag == INTERCOMM) {
      // objective: sum buffer d over masters of the task group siblings
      tmpComm = interComm;
      //
      buf = new double[n];
      for (int i=0; i<n; ++i) buf[i] = 0.0;
   }
   else if (flag == MASTERCOMM) {
      tmpComm = masterComm;
   }
   else if (flag == CHILDCOMM) {
      // objective: sum buffer d over masters of the child task groups
      tmpComm = childComm;
      //
      buf = new double[n];
      for (int i=0; i<n; ++i) buf[i] = 0.0;
   }
   else {
      sxprintf("Error %s : Unknown communicator type.\n", SX_FUNC);
      SX_EXIT;
   }

   // (1) reduce operations
   if (tmpComm != SxMpi::mpiCommNull())
   {
      if      (flag == INTRACOMM) {
         for (int i=0; i<n; ++i) buf[i] = d[i];
         SxMpi::allreduce(buf, d, n, SxMpi::SUM, tmpComm);
      }
      else if (flag == INTERCOMM) {
         SxMpi::allreduce(d,  buf, n, SxMpi::SUM, tmpComm);
      }
      else if (flag == MASTERCOMM) {
         // nothing to do
      }
      else if (flag == CHILDCOMM) {
         SxMpi::reduce(d, buf, n, SxMpi::SUM, tmpComm);
      }
   }

   // (2) broadcast operations to all local task group members
   {
      if      (flag == INTRACOMM)
         delete [] buf;
      else if (flag == INTERCOMM) {
         // one of the group ranks now has the result, the others have
         // buf zeroed ---> sum operation is a broadcast
         SxMpi::allreduce(buf, d, n, SxMpi::SUM, intraComm);
         delete [] buf;
      }
      else if (flag == MASTERCOMM) {
         SxMpi::broadcast(d, n, intraComm);
      }
      else if (flag == CHILDCOMM) {
         // one of the group ranks now has the result, the others have
         // buf zeroed ---> sum operation is a broadcast
         SxMpi::allreduce(buf, d, n, SxMpi::SUM, intraComm);
         delete [] buf;
      }
   }
}
#else
void SxMpiTaskGroup::sum(double *, const int &, const CommType &, const SxMpiComm &)
{ /* empty */}
#endif /* USE_MPI */
