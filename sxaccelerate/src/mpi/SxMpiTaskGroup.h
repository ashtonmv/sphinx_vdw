#ifndef _SX_MPITASKGROUP_H_
#define _SX_MPITASKGROUP_H_

#include <SxMpiComm.h>
#include <SxString.h>
#include <SxTaskGroup.h>

class SX_EXPORT_MPI SxMpiTaskGroup : public SxTaskGroup
{
   public:
      SxMpiTaskGroup();
      SxMpiTaskGroup(const SxMpiTaskGroup &tg);
      ~SxMpiTaskGroup() {}

      SxMpiTaskGroup &operator=(const SxMpiTaskGroup &tg);

      void setIntraCommunicator(const SxMpiComm comm);
      void setInterCommunicator(const SxMpiComm comm);

      SxMpiComm getIntraCommunicator();
      SxMpiComm getInterCommunicator();
      SxMpiComm getMasterCommunicator();

      //int size();
      //int rank();
      bool master();
      void printInfo();
      void barrier();

      // sum over all ranks of the task group
      void sum(double *d, const ssize_t &n);
      // sum over all master ranks of all childs labeled childId
      void sum(double *d, const ssize_t &n, const SxString &childId);


   private:
      /* MPI communicators */
      // communicator which includes all ranks/members belonging to the current class instance
      SxMpiComm intraComm;
      // communicator which includes only rank 0 (/ the "master") of the current class instance
      SxMpiComm masterComm;
      // communicator which includes all masters of all siblings
      SxMpiComm interComm;

      typedef enum
      {
         INTRACOMM,
         INTERCOMM,
         MASTERCOMM,
         CHILDCOMM
      }
      CommType;

      void sum(double *d, const int &n, const CommType &flag, const SxMpiComm &childComm);
};

#endif
