#include <SxString.h>
#include <SxParallelHierarchy.h>
#include <SxTaskGroup.h>
#include <SxLoopMPI.h>
#include <math.h>

#include <SxMpiTestbed.h>


// a trivial data parallel dummy compute kernel
double computeKernel()
{
   double pi;
   pi = 4.0*atan(-1.0);
   return ( sin(pi)*sin(pi) + cos(pi)*cos(pi) );
}



int main(int argc, char * argv[])
{
   const int nk = 5;
   const int nj = 5;

   SxLoopMPI::init (argc, argv);

   SxParallelHierarchy ph;

   SxTaskGroup * tgTopLevel;
   tgTopLevel = ph.getTaskGroup("top-level");

   tgTopLevel->printInfo();

   // 1: use a flat parallelization `a la LoopMpi
   {
      SxTaskGroup * tgTopAll;
      tgTopAll = tgTopLevel->getChild("top-all");

      tgTopAll->printInfo();

      double sum = 0.0;

      for (int ik = 0; ik < nk; ik++)
      {
         if (! tgTopAll->myWork(ik))
            continue;

         for (int ij = 0; ij < nj; ij++)
         {
            sum += computeKernel();
         }
      }

      sum = tgTopAll->sum(sum);
      if (tgTopAll->master())
         sxprintf("(1)  sum = %f\n", sum);
   }

   return 0;

//   // 2: use a hierarchical parallelization but no hierarchical communication
//   {
//      SxTaskGroup * tgKpoints;
//      SxTaskGroup * tgStates;
//
//      tgKpoints = tgTopLevel->getChild("k-points");
//      tgStates = tgKpoints->getChild("states");
//
//      tgKpoints->printInfo();
//      tgStates->printInfo();
//
//      double sum = 0.0;
//
//      for (int ik = 0; ik < nk; ik++)
//      {
//         if (! tgKpoints->myWork(ik))
//            continue;
//
//         for (int ij = 0; ij < nj; ij++)
//         {
//            if (! tgStates->myWork(ij))
//               continue;
//
//            sum += computeKernel();
//         }
//      }
//
//      SxTaskGroup * tgTopAll;
//      tgTopAll = tgTopLevel->getChild("top-all");
//
//      sum = tgTopAll->sum(sum);
//      if (tgTopAll->master())
//         sxprintf("(2)  sum = %f\n", sum);
//   }
//
//
//   // 3: use a hierarchical parallelization with hierarchical communication
//   {
//      SxTaskGroup * tgKpoints;
//      SxTaskGroup * tgStates;
//
//      tgKpoints = tgTopLevel->getChild("k-points");
//      tgStates = tgKpoints->getChild("states");
//
//      double sum = 0.0;
//
//      for (int ik = 0; ik < nk; ik++)
//      {
//         if (! tgKpoints->myWork(ik))
//            continue;
//
//         for (int ij = 0; ij < nj; ij++)
//         {
//            if (! tgStates->myWork(ij))
//               continue;
//
//            sum += computeKernel();
//         }
//      }
//      sum = tgStates   -> sum(sum            );
//      sum = tgKpoints  -> sum(sum, "states"  );
//      sum = tgTopLevel -> sum(sum, "k-points");
//
//      if (tgTopLevel->master())
//         sxprintf("(3)  sum = %f\n", sum);
//   }
//
//
//   // 4: use a hierarchical parallelization with hierarchical communication
//   {
//      SxTaskGroup * tgKpoints;
//      SxTaskGroup * tgStates;
//      tgKpoints = tgTopLevel->getChild("k-points");
//      tgStates = tgKpoints->getChild("states");
//
//      double sum = 0.0;
//
//      for (int ik = 0; ik < nk; ik++)
//      {
//         if (! tgKpoints->myWork(ik))
//            continue;
//
//         for (int ij = 0; ij < nj; ij++)
//         {
//            if (! tgStates->myWork(ij))
//               continue;
//
//            sum += computeKernel();
//         }
//      }
//      sum = tgTopLevel->sum( tgKpoints->sum( tgStates->sum(sum), "states" ), "k-points");
//
//      if (tgTopLevel->master())
//         sxprintf("(4)  sum = %f\n", sum);
//   }
//
//
//
//   // 5: work with the taskGroups' index generators
//   {
//      SxTaskGroup * tgKpoints;
//      tgKpoints = tgTopLevel->getChild("k-points");
//      tgKpoints->setIndexBounds(0, nk);   // can be done only once
//
//      SxTaskGroup * tgStates;
//      tgStates = tgKpoints->getChild("states");
//      tgStates->setIndexBounds(0, nj);
//
//      double sum = 0.0;
//
//      for (int ik = tgKpoints->getIdxFrom(); ik < tgKpoints->getIdxTo(); ik++)
//      {
//         for (int ij = tgStates->getIdxFrom(); ij < tgStates->getIdxTo(); ij++)
//         {
//            sum += computeKernel();
//         }
//      }
//      sum = tgTopLevel->sum( tgKpoints->sum( tgStates->sum(sum), "states" ), "k-points");
//
//      if (tgTopLevel->master())
//         sxprintf("(5)  sum = %f\n", sum);
//   }
//
//
//
//   // 6: work with index iterators, use a sum-up method to go from the node to the root
//   {
//      SxTaskGroup * tgKpoints;
//      tgKpoints = tgTopLevel->getChild("k-points");
//      tgKpoints->setIndexBounds(0, nk);   // can be done only once
//
//      SxTaskGroup * tgStates;
//      tgStates = tgKpoints->getChild("states");
//      tgStates->setIndexBounds(0, nj);
//
//      double sum = 0.0;
//
//      for (SxList<int>::ConstIterator ik = tgKpoints->getIndexIterator(SxTaskGroup::FIRST);
//            ik != tgKpoints->getIndexIterator(SxTaskGroup::LAST); ++ik)
//      {
//         for (SxList<int>::ConstIterator ij = tgStates->getIndexIterator(SxTaskGroup::FIRST);
//               ij != tgStates->getIndexIterator(SxTaskGroup::LAST); ++ij)
//         {
//            // useful operations with *ij and *ik are possible, of course
//            sum += computeKernel();
//         }
//      }
//
//      // instead of using
////      sum = tgTopLevel->sum( tgKpoints->sum( tgStates->sum(sum), "states" ), "k-points");
//      // we can simply use
//      sum = tgStates->sumUp(sum);
//
//      if (tgTopLevel->master())
//         sxprintf("(6)  sum = %f\n", sum);
//   }
//
//
//
//   // 7: use hybrid parallelism -- threads & mpi tasks
//
//
//
//   // no need to explicitly call mpi finalize -- happens automatically w/ the destructor
//   return 0;
}
