#include <SxList.h>
#include <SxMap.h>
#include <SxString.h>
#include <SxTaskGroup.h>
#include <SxError.h>


SxTaskGroup::SxTaskGroup()
{
   nMembers = nSiblings = memberRank = siblingRank = -1;
}


SxTaskGroup::SxTaskGroup(const SxTaskGroup &tg)
{
   if (verboseConstructors)
      sxprintf("%s called for %s\n", SX_FUNC, tg.name.ascii());
   //
   nMembers = tg.nMembers;
   nSiblings = tg.nSiblings;
   memberRank = tg.memberRank;
   siblingRank = tg.siblingRank;
   name = tg.name;
   autoLvlName = tg.autoLvlName;
   parentId = tg.parentId;
   parent = tg.parent;
   childId = tg.childId;
   children = tg.children;
}

SxTaskGroup &SxTaskGroup::operator=(const SxTaskGroup &tg)
{
   if (verboseConstructors)
      sxprintf("%s called for %s\n", SX_FUNC, tg.name.ascii());
   //
   nMembers = tg.nMembers;
   nSiblings = tg.nSiblings;
   memberRank = tg.memberRank;
   siblingRank = tg.siblingRank;
   name = tg.name;
   autoLvlName = tg.autoLvlName;
   parentId = tg.parentId;
   childId = tg.childId;
   parent = tg.parent;
   children = tg.children;
   return *this;
}


void SxTaskGroup::setNmembers(const int &n)
{
   nMembers = n;
}

void SxTaskGroup::setNsiblings(const int &n)
{
   nSiblings = n;
}

void SxTaskGroup::setSiblingRank(const int &n)
{
   siblingRank = n;
}

void SxTaskGroup::setMemberRank(const int &n)
{
   memberRank = n;
}

void SxTaskGroup::setName(const SxString &name_)
{
   name = name_;
}

const SxString SxTaskGroup::getName()
{
   return name;
}

void SxTaskGroup::setAutoLvlName(const SxString &name_)
{
   autoLvlName = name_;
}

const SxString SxTaskGroup::getAutoLvlName()
{
   return autoLvlName;
}

bool SxTaskGroup::hasAutoLvlName()
{
   return (autoLvlName.getSize() > 0);
}


/* virtual */ bool SxTaskGroup::master()
{
   if (memberRank == 0)
      return true;
   else
      return false;
}

namespace {
   void errorTruncatedHierarchy ()
   {
      cout << "Fatal problem: lowest MPI level does not exhaust MPI tasks"
           << endl;
      cout << "This may lead to wrong results" << endl;
      cout << "Check the parallel hierarchy" << endl;
   }
}

/* virtual */ bool SxTaskGroup::myWork(int ik)
{
   if (hasChildren())
      return ((ik % nSiblings) == siblingRank);
   else
      return ((ik % nMembers) == memberRank);
}

bool SxTaskGroup::notMyWork(int ik)
{
   return (! myWork(ik));
}

/* virtual */ int SxTaskGroup::whoseWork(int ik)
{
   if (hasChildren())
      return ik % nSiblings;
   else
      return ik % nMembers;
}


void SxTaskGroup::setParentId(const SxString &_parentId)
{
   parentId = _parentId;
   return;
}

const SxString& SxTaskGroup::getParentId() const
{
   return parentId;
}

void SxTaskGroup::addChildId(const SxString &_childId)
{
   if (! childId.contains(_childId))
      childId.append(_childId);
   return;
}

bool SxTaskGroup::hasChildren()
{
   return (children.getSize () > 0);
}

bool SxTaskGroup::isLeaf()
{
   return (! hasChildren());
}

void SxTaskGroup::addChild(SxTaskGroup * child)
{
   SX_CHECK( child );
   children( child->getName() ) = child;
}

SxTaskGroup * SxTaskGroup::getChild(const SxString &_childId)
{
   SX_CHECK( children.containsKey(_childId) );
   return children(_childId);
}

void SxTaskGroup::setParent(SxTaskGroup * _parent)
{
   parent = _parent;
   return;
}

SxTaskGroup * SxTaskGroup::getParent()
{
   SX_CHECK(parent != NULL);
   return parent;
}

/* virtual */ void SxTaskGroup::printInfo()
{
   return;
}


double SxTaskGroup::sum(const double &d)
{
   double e = d;
   sum(&e, 1);
   return e;
}

double SxTaskGroup::sum(const double &d, const SxString &childIdIn)
{
   double e = d;
   sum(&e, 1, childIdIn);
   return e;
}


double SxTaskGroup::sumUp(const double &d, const SxString upTo)
{
   double e = d;
   sumUp(&e, 1, upTo);
   return e;
}

void SxTaskGroup::sumUp(void *c, const ssize_t &n, const SxString upTo)
{
   sumUp((double*)c, 2*n, upTo);
}

void SxTaskGroup::sumUp(double *d, const ssize_t &n, const SxString upTo)
{
   sum(d, n);

   SxTaskGroup * tg;
   tg = this;

   SxString parName, nodeName;
   parName = tg->getParentId();

   while (parName != "ROOT")
   {
      nodeName = tg->getName();
      tg = tg->getParent();
      tg->sum(d, n, nodeName);
      if (nodeName == upTo) break;
      parName = tg->getParentId();
   }
}

void SxTaskGroup::write (ostream &out, const SxString &indent) const
{
   out << indent << "level {" << endl;
   out << indent << "   name = \"" << name << "\";";
   out << endl;
   out << indent << "   siblings = " << nSiblings << ";" << endl;
   out << indent << "   members  = " << nMembers << ";" << endl;
   SxMap<SxString, SxTaskGroup*>::ConstIterator it;
   for (it = children.begin (); it != children.end (); ++it)
      it.getValue ()->write (out, indent + "   ");
   out << indent << "}" << endl;
}
