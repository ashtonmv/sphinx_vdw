#include <SxMpiComm.h>
#include <SxString.h>
#include <SxError.h>
#include <SxConfig.h>
#ifdef USE_MPI
#include <mpi.h>
#include <typeinfo>
#endif


SxMpiComm::SxMpiComm()
{
#ifdef USE_MPI
      comm = MPI_COMM_NULL;
#endif
}

// --- Do not use the MPI free,dup,compare functions!
//     If we stick to the plain pointer (the native MPI_Comm datatype)
//     everything should be fine -- and much simpler.

//SxMpiComm::~SxMpiComm()
//{
//#ifdef USE_MPI
//   MPI_Comm_free(&comm);
//#endif
//}
//
//SxMpiComm::SxMpiComm(const SxMpiComm &mpic)
//{
//#ifdef USE_MPI
//   MPI_Comm_dup(mpic.comm, &comm);
//#endif
//}
//
//SxMpiComm &SxMpiComm::operator=(const SxMpiComm &mpic)
//{
//#ifdef USE_MPI
//   MPI_Comm_dup(mpic.comm, &comm);
//#endif
//   return *this;
//}

#ifdef USE_MPI
bool SxMpiComm::operator== (const SxMpiComm &other) const
{
   return (comm == other.comm);
//   int result;
//   MPI_Comm_compare(comm, other.comm, &result);
//   return (result == MPI_IDENT);
}
#else
bool SxMpiComm::operator== (const SxMpiComm &) const { return true; }
#endif

bool SxMpiComm::operator!= (const SxMpiComm &other) const
{
   return !(*this == other);
}


#ifdef USE_MPI
void SxMpi::init(int *argc, char ***argv)
{
   MPI_Init(argc, argv);
}
#else
void SxMpi::init(int *, char ***) {}
#endif

void SxMpi::finalize()
{
#ifdef USE_MPI
   MPI_Finalize();
#endif
}

void SxMpi::stop()
{
   finalize();
   exit(1);
}

#ifdef USE_MPI
void SxMpi::barrier(SxMpiComm comm /* = MPI_COMM_WORLD */)
{
   MPI_Barrier(comm.comm);
}
#else
void SxMpi::barrier(SxMpiComm) {}
#endif

#ifdef USE_MPI
int SxMpi::size(SxMpiComm comm /* = MPI_COMM_WORLD */)
{
   int n;
   MPI_Comm_size(comm.comm, &n);
   return n;
}
#else
int SxMpi::size(SxMpiComm) { return 1; }
#endif

#ifdef USE_MPI
int SxMpi::rank(SxMpiComm comm /* = MPI_COMM_WORLD */)
{
   int r;
   MPI_Comm_rank (comm.comm, &r);
   return r;
}
#else
int SxMpi::rank(SxMpiComm) { return 0; }
#endif

SxMpiComm SxMpi::mpiCommWorld()
{
   SxMpiComm c;
#ifdef USE_MPI
   c.comm = MPI_COMM_WORLD;
#endif
   return c;
}

SxMpiComm SxMpi::mpiCommNull()
{
   SxMpiComm c;
#ifdef USE_MPI
   c.comm = MPI_COMM_NULL;
#endif
   return c;
}

SxMpiComm SxMpi::mpiCommSelf()
{
   SxMpiComm c;
#ifdef USE_MPI
   c.comm = MPI_COMM_SELF;
#endif
   return c;
}


#ifdef USE_MPI
SxMpiComm SxMpi::divideMpiComm(SxMpiComm _commIn, int nParts)
{
   SxMpiComm _newComm;

   MPI_Comm commIn;
   MPI_Comm newComm;
   MPI_Group inGroup, newGroup;
   int inRank, inSize;
   int ranges[1][3];
   int n;

   commIn = _commIn.comm;

   MPI_Comm_rank(commIn, &inRank);
   MPI_Comm_size(commIn, &inSize);
   MPI_Comm_group(commIn, &inGroup);

   // set up n and ranges
   if (nParts > inSize)
   {
      sxprintf("ERROR: %s : (nParts < inSize)!\n", SX_FUNC);
      SX_EXIT;
   }
   else
   {
      int partSize = inSize / nParts;
      int nGroup = inRank / partSize;
      ranges[0][0] = partSize * nGroup;
      //
      if ((inSize % nParts) == 0)
      {
         ranges[0][1] = ranges[0][0] + partSize - 1;
      }
      else
      {
         sxprintf("ERROR: %s : ((inSize MOD nParts) != 0)!\n", SX_FUNC);
         SX_EXIT;
      }
      //
      ranges[0][2] = 1;
      n = 1;
   }

   MPI_Group_range_incl(inGroup, n, ranges, &newGroup);
   MPI_Comm_create(commIn, newGroup, &newComm);

   _newComm.comm = newComm;

   return _newComm;
}
#else
SxMpiComm SxMpi::divideMpiComm(SxMpiComm, int) { return SxMpiComm (); }
#endif /* USE_MPI */



#ifdef USE_MPI
SxMpiComm SxMpi::extractMpiComm(SxMpiComm _commIn, int * ranks, int nRanks)
{
   SxMpiComm _newComm;

   MPI_Comm commIn;
   MPI_Comm newComm;
   MPI_Group inGroup, newGroup;

   commIn = _commIn.comm;

   if (commIn != MPI_COMM_NULL)
   {
      MPI_Comm_group(commIn, &inGroup);
      MPI_Group_incl(inGroup, nRanks, ranks, &newGroup);
      MPI_Comm_create(commIn, newGroup, &newComm);
   }
   else
   {
      newComm = MPI_COMM_NULL;
   }

   _newComm.comm = newComm;

   return _newComm;
}
#else
SxMpiComm SxMpi::extractMpiComm(SxMpiComm, int *, int) { return SxMpiComm (); }
#endif /* USE_MPI */



#ifdef USE_MPI
SxMpiComm SxMpi::joinMpiComm(SxMpiComm _commChild, SxMpiComm _commParent)
{
   SxMpiComm _commOut;

   MPI_Comm commOut, commChild, commParent;
   MPI_Group groupIn, groupOut, groupParent;
   int inRank, inSize, parentSize;
   int outSize;
   int * rankMask;
   int * rankMaskTmp;
   int * ranksOut;

   commChild = _commChild.comm;
   commParent = _commParent.comm;

   commOut = MPI_COMM_NULL;
   groupIn = groupOut = groupParent = MPI_GROUP_NULL;
   inRank = inSize = parentSize = -1;

   if (commParent != MPI_COMM_NULL)
   {
      MPI_Comm_size(commParent, &parentSize);
      MPI_Comm_group(commParent, &groupParent);
   }
   else
   {
      sxprintf("ERROR: %s : (commParent == MPI_COMM_NULL)\n", SX_FUNC);
      SX_EXIT;
   }

   rankMask = new int[parentSize];
   rankMaskTmp = new int[parentSize];

   for (int i=0; i<parentSize; ++i)
   {
      rankMask[i] = 0;
      rankMaskTmp[i] = 0;
   }

   if (commChild != MPI_COMM_NULL)
   {
      MPI_Comm_size(commChild, &inSize);
      MPI_Comm_rank(commChild, &inRank);
      MPI_Comm_group(commChild, &groupIn);

      int * ranksIn;
      int * ranksParent;
      ranksIn = new int[inSize];
      ranksParent = new int[inSize];

      for (int i=0; i<inSize; i++)
         ranksIn[i] = i;

      MPI_Group_translate_ranks(groupIn, inSize, ranksIn, groupParent, ranksParent);

      for (int i=0; i<inSize; ++i)
         rankMaskTmp[ranksParent[i]] = 1;

      delete [] ranksIn;
      delete [] ranksParent;
   }

   MPI_Reduce(rankMaskTmp, rankMask, parentSize, MPI_INT, MPI_LOR, 0, commParent);
   MPI_Bcast(rankMask, parentSize, MPI_INT, 0, commParent);

   outSize = 0;
   for (int i=0; i<parentSize; ++i)
      outSize += rankMask[i];

   ranksOut = new int[outSize];
   for (int i=0, j=0; i<parentSize; ++i)
      if (rankMask[i])
         ranksOut[j] = i, ++j;

   MPI_Group_incl(groupParent, outSize, ranksOut, &groupOut);
   MPI_Comm_create(commParent, groupOut, &commOut);

   _commOut.comm = commOut;


   return _commOut;
}
#else
SxMpiComm SxMpi::joinMpiComm(SxMpiComm, SxMpiComm) { return SxMpiComm (); }
#endif /* USE_MPI */



// --- MPI wrapper template implementations -----------------------------------
//
//   - we cannot put these implementations into <SxMpi.h> due to the MPI
//     functions and datatypes from <mpi.h> which we do *not* want to be
//     visible/necessary anywhere where <SxMpi.h> is included
//   - caller functions are implemented further below, these caller functions
//     appear with their declarations in <SxMpi.h>
//
#ifdef USE_MPI

template <class T>
   void SxMpi_reduce(T *in, T *out, const int &nelem, const SxMpi::mpiOp &_mpiOp, const SxMpiComm &_comm)
{
   MPI_Datatype dt = getMpiDataType (in);
   MPI_Op op;
   switch (_mpiOp)
   {
      case SxMpi::SUM:
         op = MPI_SUM;
         break;
      default:
         SX_EXIT;
         break;
   }
   //
   int status = MPI_Reduce(in, out, nelem, dt, op, 0, _comm.comm);
   if (status != MPI_SUCCESS) { SX_EXIT; }
}
#else
template <class T>
void SxMpi_reduce(T *, T *, const int &, const SxMpi::mpiOp &, const SxMpiComm &) { }
#endif


#ifdef USE_MPI
template <class T>
   void SxMpi_allreduce(T *in, T *out, const int &nelem, const SxMpi::mpiOp &_mpiOp, const SxMpiComm &_comm)
{
   MPI_Datatype dt = getMpiDataType (in);
   MPI_Op op;
   switch (_mpiOp)
   {
      case SxMpi::SUM:
         op = MPI_SUM;
         break;
      default:
         SX_EXIT;
         break;
   }
   //
   int status = MPI_Allreduce(in, out, nelem, dt, op, _comm.comm);
   if (status != MPI_SUCCESS) { SX_EXIT; }
}
#else
template <class T>
void SxMpi_allreduce(T *, T *, const int &, const SxMpi::mpiOp &,
                     const SxMpiComm &)
{}
#endif


#ifdef USE_MPI
template <class T>
   void SxMpi_broadcast(T *inout, const int &nelem, const SxMpiComm &_comm)
{
   MPI_Datatype dt = getMpiDataType (inout);
   int status = MPI_Bcast(inout, nelem, dt, 0, _comm.comm);
   if (status != MPI_SUCCESS) { SX_EXIT; }
}
#else
template <class T> void SxMpi_broadcast(T *, const int &, const SxMpiComm &) {}
#endif

// --- caller functions for the MPI template wrappers -------------------------
void SxMpi::reduce(double *in, double *out, const int &nelem, const SxMpi::mpiOp &_mpiOp, const SxMpiComm &_comm)
   { SxMpi_reduce(in, out, nelem, _mpiOp, _comm); }
void SxMpi::reduce(int    *in, int    *out, const int &nelem, const SxMpi::mpiOp &_mpiOp, const SxMpiComm &_comm)
   { SxMpi_reduce(in, out, nelem, _mpiOp, _comm); }

void SxMpi::allreduce(double *in, double *out, const int &nelem, const SxMpi::mpiOp &_mpiOp, const SxMpiComm &_comm)
   { SxMpi_allreduce(in, out, nelem, _mpiOp, _comm); }
void SxMpi::allreduce(int    *in, int    *out, const int &nelem, const SxMpi::mpiOp &_mpiOp, const SxMpiComm &_comm)
   { SxMpi_allreduce(in, out, nelem, _mpiOp, _comm); }

void SxMpi::broadcast(double *inout, const int &nelem, const SxMpiComm &_comm)
   { SxMpi_broadcast(inout, nelem, _comm); }
void SxMpi::broadcast(int    *inout, const int &nelem, const SxMpiComm &_comm)
   { SxMpi_broadcast(inout, nelem, _comm); }
