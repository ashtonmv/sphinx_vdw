#ifndef _SX_TASKGROUP_H_
#define _SX_TASKGROUP_H_

#include <SxList.h>
#include <SxMap.h>
#include <SxString.h>
#include <SxMpiComm.h>

class SX_EXPORT_MPI SxTaskGroup
{
   public:
      SxTaskGroup();
      virtual ~SxTaskGroup() {};

      SxTaskGroup(const SxTaskGroup &tg);
      SxTaskGroup &operator=(const SxTaskGroup &tg);

      void setName (const SxString &name_);
      const SxString getName ();

      void setAutoLvlName (const SxString &name_);
      const SxString getAutoLvlName ();
      bool hasAutoLvlName();

      void setNmembers (const int &n);
      int getNmembers () { return nMembers; }

      void setNsiblings (const int &n);
      int getNsiblings () { return nSiblings; }

      void setSiblingRank (const int &n);
      int getSiblingRank() { return siblingRank; }

      void setMemberRank (const int &n);
      int getMemberRank() { return memberRank; }


      void addChild(SxTaskGroup * child);
      SxTaskGroup * getChild(const SxString &childId);

      void setParent(SxTaskGroup * _parent);
      SxTaskGroup * getParent();


      /// String IDs are used to temporarily store parent-child relations
      /// during the first phase of the hierarchy build-up process.
      /// The second phase then sets up pointers based on these string IDs.
      /// See SxParallelHierarchy::buildHierarchy() for details.
      void setParentId(const SxString &parentId);
      const SxString &getParentId() const;
      SxString parentId;
      void addChildId(const SxString &childId);
      SxList<SxString> childId;

      bool hasChildren();
      bool isLeaf();


      virtual bool master();
      virtual bool myWork(int ik);
      virtual int whoseWork (int ik);
      bool notMyWork(int ik);

      virtual void printInfo() = 0;

      virtual void barrier() = 0;


// --- summation routines
      double sum(const double &d);
      double sum(const double &d, const SxString &childId);
      virtual void sum(double *d, const ssize_t &n) = 0;
      virtual void sum(double *d, const ssize_t &n, const SxString &childId) = 0;
      //
      double sumUp(const double &d,                   const SxString upTo="ROOT");
      void   sumUp(      double *d, const ssize_t &n, const SxString upTo="ROOT");
      void   sumUp(        void *c, const ssize_t &n, const SxString upTo="ROOT");
//      void sumUp(SxComplex<double> *c, const ssize_t &n);
// ---

   protected:
      static const bool verboseConstructors = false;

      int nMembers;     // number of processing units within the task group
      int nSiblings;    // number of siblings (task group instances with the same parent task group)
      int memberRank;   // (MPI/OpenMP/PThread) rank of the processing unit within the task group
      int siblingRank;  // rank of the task group within the siblings

      SxString name;    // (local) name of the task group
      SxString autoLvlName; // copy of the global name used by the task group maps in SxParallelHierarchy


   public:
      /// pointer to the parent task group ()
      SxTaskGroup * parent;

      /// pointer map to all the children task groups
      SxMap<SxString, SxTaskGroup*> children;

      /// write out level in parallelHierarchy.sx format
      void write (ostream &out, const SxString &indent = "") const;
};

#endif
