#ifndef _SX_MPI_NEW_H_
#define _SX_MPI_NEW_H_

#include <SxConfig.h>
#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef WIN32
#  if defined(_EXPORT_sxmpi)
#     define SX_EXPORT_MPI __declspec(dllexport)
#  else
#     define SX_EXPORT_MPI __declspec(dllimport)
#  endif
#else
#  define SX_EXPORT_MPI
#endif


class SX_EXPORT_MPI SxMpiComm
{
   public:
#ifdef USE_MPI
      MPI_Comm comm;
#endif
      SxMpiComm();

//      ~SxMpiComm();
//      SxMpiComm(const SxMpiComm &mpic);
//      SxMpiComm &operator=(const SxMpiComm &mpic);

      bool operator== (const SxMpiComm &other) const;
      bool operator!= (const SxMpiComm &other) const;
};

class SX_EXPORT_MPI SxMpi
{
   public:
      static void init(int * argc, char *** argv);
      static void finalize();
      static void stop();

      static SxMpiComm mpiCommWorld();
      static SxMpiComm mpiCommNull();
      static SxMpiComm mpiCommSelf();

      static void barrier(SxMpiComm comm = mpiCommWorld());
      static int     size(SxMpiComm comm = mpiCommWorld());
      static int     rank(SxMpiComm comm = mpiCommWorld());

      // sub-divides an MPI communicator into nParts
      static SxMpiComm divideMpiComm (SxMpiComm _commIn, int nParts);
      // creates a new MPI communicator by extracting a subset of ranks from an existing one
      static SxMpiComm extractMpiComm (SxMpiComm _commIn, int * ranks, int nRanks);
      // join separate communicators commChild, the ranks of which are part of commParent
      static SxMpiComm joinMpiComm (SxMpiComm _commChild, SxMpiComm _commParent);

      typedef enum { SUM } mpiOp;

      static void    reduce(double *in, double *out, const int &nelem, const mpiOp &_mpiOp, const SxMpiComm &_comm);
      static void    reduce(int    *in, int    *out, const int &nelem, const mpiOp &_mpiOp, const SxMpiComm &_comm);

      static void allreduce(double *in, double *out, const int &nelem, const mpiOp &_mpiOp, const SxMpiComm &_comm);
      static void allreduce(int    *in, int    *out, const int &nelem, const mpiOp &_mpiOp, const SxMpiComm &_comm);

      static void broadcast(double *inout,           const int &nelem,                      const SxMpiComm &_comm);
      static void broadcast(int    *inout,           const int &nelem,                      const SxMpiComm &_comm);
};

#ifdef USE_MPI
inline MPI_Datatype getMpiDataType (double *) { return MPI_DOUBLE; }
inline MPI_Datatype getMpiDataType (int *)    { return MPI_INT; }
inline MPI_Datatype getMpiDataType (char *)    { return MPI_CHAR; }
#endif

#endif
