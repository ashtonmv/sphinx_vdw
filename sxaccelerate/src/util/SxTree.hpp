// ---------------------------------------------------------------------------
//
//           The general purpose cross platform C/C++ framework
//
//                       S x A c c e l e r a t e
//
//           Home:       https://www.sxlib.de
//           License:    Apache 2
//           Authors:    see src/AUTHORS
//
// ---------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Tree::Iterator-class
//----------------------------------------------------------------------------

// Constructors
// -- Constructor
template<class T>
SxTree<T>::Iterator::Iterator () : ptr (NULL)
{
}

// -- Copyconstructor
template<class T>
SxTree<T>::Iterator::Iterator (const Iterator &rhs) : ptr (rhs.ptr)
{
   // empty
}

template<class T>
SxTree<T>::Iterator::Iterator (SxTree<T> *rhs) : ptr (rhs)
{
   // empty
}

// Operations
// -- Assignment
template<class T>
typename SxTree<T>::Iterator &
SxTree<T>::Iterator::operator= (const Iterator &rhs)
{
   ptr = rhs.ptr;

   return *this;
}

template<class T>
typename SxTree<T>::Iterator &SxTree<T>::Iterator::operator= (SxTree<T> *rhs)
{
   ptr = rhs;
   return *this;
}

// -- Equality
template<class T>
bool SxTree<T>::Iterator::operator== (const typename SxTree<T>::Iterator &rhs) const
{
   return ptr == rhs.ptr;
}

template<class T>
bool SxTree<T>::Iterator::operator!= (const typename SxTree<T>::Iterator &rhs) const
{
   return ptr != rhs.ptr;
}

// -- Dereferencation
template<class T>
SxTree<T> &SxTree<T>::Iterator::operator* ()
{
   SX_CHECK (ptr);
   return *ptr;
}

template<class T>
SxTree<T> *SxTree<T>::Iterator::operator-> ()
{
   SX_CHECK (ptr);
   return ptr;
}

// -- In- and Decrementation
template<class T>
typename SxTree<T>::Iterator &SxTree<T>::Iterator::operator++ ()
{
   return *this = ptr->rightSibling;
}

template<class T>
typename SxTree<T>::Iterator &SxTree<T>::Iterator::operator++ (int)
{
   Iterator ret (*this); ++(*this); return ret;
}

template<class T>
void SxTree<T>::Iterator::operator+= (unsigned int steps)
{
   for (unsigned int iSteps = 0; iSteps < steps; ++iSteps)  {
      ++(*this);
   }
   return ; //*this;
}

template<class T>
typename SxTree<T>::Iterator &SxTree<T>::Iterator::operator-- ()
{
   return *this = ptr->leftSibling;
}

template<class T>
typename SxTree<T>::Iterator &SxTree<T>::Iterator::operator-- (int)
{
   Iterator ret (*this); --(*this); return ret;
}

template<class T>
void SxTree<T>::Iterator::operator-= (unsigned int steps)
{
   for (unsigned int iSteps = 0; iSteps < steps; ++iSteps)  {
      --(*this);
   }
   return ; //*this;
}

//----------------------------------------------------------------------------
// Tree::ConstIterator-class
//----------------------------------------------------------------------------

// Constructors
// -- Constructor
template<class T>
SxTree<T>::ConstIterator::ConstIterator () : ptr (NULL)
{
   // empty
}

// -- Copyconstructor
template<class T>
SxTree<T>::ConstIterator::ConstIterator (const typename SxTree<T>::ConstIterator &rhs)
: ptr (rhs.ptr)
{
   // empty
}

template<class T>
SxTree<T>::ConstIterator::ConstIterator (const typename SxTree<T>::Iterator &rhs)
: ptr (rhs.getConstPtr ())
{
   // empty
}

template<class T>
SxTree<T>::ConstIterator::ConstIterator (const SxTree<T> *rhs)
: ptr (rhs)
{
   // empty
}

// Operations
// -- Assignment
template<class T>
typename SxTree<T>::ConstIterator &
SxTree<T>::ConstIterator::operator= (const typename SxTree<T>::ConstIterator &rhs)
{
   ptr = rhs.ptr;

   return *this;
}

template<class T>
typename SxTree<T>::ConstIterator &
SxTree<T>::ConstIterator::operator= (const SxTree<T> *rhs)
{
   ptr = rhs;
   return *this;
}

// -- Equality
template<class T>
bool SxTree<T>::ConstIterator::operator== (const typename SxTree<T>::ConstIterator &rhs) const
{
   return ptr == rhs.ptr;
}

template<class T>
bool SxTree<T>::ConstIterator::operator!= (const typename SxTree<T>::ConstIterator &rhs) const
{
   return ptr != rhs.ptr;
}

// -- Dereferencation
template<class T>
const SxTree<T> & SxTree<T>::ConstIterator::operator* ()
{
   SX_CHECK (ptr);
   return *ptr;
}

template<class T>
const SxTree<T> * SxTree<T>::ConstIterator::operator-> ()
{
   SX_CHECK (ptr);
   return ptr;
}

// -- In- and Decrementation
template<class T>
typename SxTree<T>::ConstIterator &SxTree<T>::ConstIterator::operator++ ()
{
   SX_CHECK (ptr);
   return *this = ptr->rightSibling;
}

template<class T>
typename SxTree<T>::ConstIterator &SxTree<T>::ConstIterator::operator++ (int)
{
   Iterator ret (*this); ++(*this); return ret;
}

template<class T>
void SxTree<T>::ConstIterator::operator+= (unsigned int steps)
{
   for (unsigned int iSteps = 0; iSteps < steps; ++iSteps)  {
      ++(*this);
   }
   return; //*this;
}

template<class T>
typename SxTree<T>::ConstIterator &SxTree<T>::ConstIterator::operator-- ()
{
   return *this = ptr->leftSibling;
}

template<class T>
typename SxTree<T>::ConstIterator &SxTree<T>::ConstIterator::operator-- (int)
{
   Iterator ret (*this); --(*this); return ret;
}

template<class T>
void SxTree<T>::ConstIterator::operator-= (unsigned int steps)
{
   for (unsigned int iSteps = 0; iSteps < steps; ++iSteps)  {
      --(*this);
   }
   return ; //*this;
}

//----------------------------------------------------------------------------
// Tree-class
//----------------------------------------------------------------------------

// Constructors
// -- Constructor
template<class T>
SxTree<T>::SxTree () : children (),
                       data (SxPtr<T>::create ((T) 0)),
                       nElements (1), parent (NULL),
                       leftSibling (NULL), rightSibling (NULL) 
{
   // empty
}

// -- Copyconstructor
template<class T>
SxTree<T>::SxTree (const SxTree<T> &rhs, CopyType ct)
: children (),
  data (SxPtr<T>::create ((T) 0)),
  nElements (1), parent (NULL),
  leftSibling (NULL), rightSibling (NULL)
{
   if (ct == deepCopy)  {
      createDeepCopy (rhs);
   }  else  {
      createShallowCopy (rhs);
   }
}

// shallow copy
template<class T>
void SxTree<T>::createShallowCopy (const SxTree<T> &rhs)
{

   // Replicating the data shallow
   data = rhs.data;
   nElements = rhs.nElements;

   // Traversing the children
   SxPtr<SxTree<T> > tmp;
   typename SxList<SxPtr< SxTree<T> > >::ConstIterator it;
   for(it = rhs.children.begin ();
       it != rhs.children.end ();
       ++it)  {

      // Creating an object
      tmp = SxPtr<SxTree<T> >::create ();

      // Linking
      tmp->parent = Iterator (this);

      // Connecting the new Element to the left and right siblings
      tmp->rightSibling = NULL;
      // creating a reference to the left sibling if one exits
      if (children.getSize () > 0)  {
         tmp->leftSibling  = Iterator (&(*children.last ()));
         children.last ()->rightSibling = Iterator (&(*tmp));
      }  else  {
         tmp->leftSibling = NULL;
      }

      // Appending the object to the list
      children.append (tmp);

      // Recalling the function to step one level down
      tmp->createDeepCopy(**it);

   }
}

/*
SxTree SxTree<T>::getDeepCopy () const
{
    // (1) local variable
    SxTree res;
    // (2) treat current level
    res.data = COPY data;
    // (3) 1 level down
    for (it = getChildren())  {
       res.insert (it.getDeepCopy ());
    }
    // (3.1) exit criteria(!)
    // (4) return current level
    return res;
}
*/


// deep copy
template<class T>
void SxTree<T>::createDeepCopy (const SxTree<T> &rhs)
{

   // Replicating the data
   data = SxPtr<T>::create (*rhs.data);
   //data = rhs.data;
   nElements = rhs.nElements;

   // Traversing the children
   typename SxList<SxPtr< SxTree<T> > >::ConstIterator it;
   for(it = rhs.children.begin ();
       it != rhs.children.end ();
       ++it)  {

      // Creating an object
      SxPtr<SxTree<T> > tmp = SxPtr<SxTree<T> >::create ();

      // Linking
      tmp->parent = Iterator (this);

      // Connecting the new Element to the left and right siblings
      tmp->rightSibling = NULL;
      // creating a reference to the left sibling if one exits
      if (children.getSize () > 0)  {
         tmp->leftSibling  = Iterator (&(*children.last ()));
         children.last ()->rightSibling = Iterator (&(*tmp));
      }  else  {
         tmp->leftSibling = NULL;
      }

      // Appending the object to the list
      children.append (tmp);

      // Recalling the function to step one level down
      tmp->createDeepCopy(**it);

   }
}

// -- Default constructor
template<class T>
SxTree<T>::SxTree (const T &t) : children(),
                                 data (SxPtr<T>::create (t)),
                                 nElements (1), parent (NULL),
                                 leftSibling (NULL), rightSibling (NULL) 
{
   // empty
}

// -- Destructor
template<class T>
SxTree<T>::~SxTree ()
{
   // empty
}

//----------------------------------------------------------------------------
// Operations
//----------------------------------------------------------------------------
// -- Assignment
template<class T>
SxTree<T> &SxTree<T>::operator= (const SxTree<T> &rhs)
{

  // Selftest
  if (this == &rhs)  return *this;

  // Removing the former children
  removeChildren ();

  // the default CopyType is shallowCopy
  // if you want to change the default copy behavior to deepCopy uncomment 
  // the following line and comment the line beneath
  // createDeepCopy (rhs); // make a deep copy
  createShallowCopy(rhs); // make a shallow copy

  return *this;
}

template<class T>
SxTree<T> &SxTree<T>::operator= (const T &t)
{
   // Selftest
   if (data.getPtr () == &t) return *this;

   data = SxPtr<T>::create (t);

   // Update the parent pointer of the children
   for (typename SxList<SxPtr<SxTree<T> > >::Iterator
         it  = children.begin ();
         it != children.end ();
         ++it)  {
      (*it)->parent = Iterator (this);
   }

   return *this;
}

// -- Input
template<class T>
SxTree<T> &SxTree<T>::operator<< (const T &t)
{
   insert (t);

   return *this;
}

template<class T>
SxTree<T> &SxTree<T>::operator<< (const SxTree<T> &tree)
{
   insert (tree);

   return *this;
}

//----------------------------------------------------------------------------
// Methods
//----------------------------------------------------------------------------
// -- Accessfunctions
template<class T>
SxTree<T> &SxTree<T>::setValue (const T &t)
{
   return *this = t;
}

// -- Insertion
template<class T>
SxTree<T> &SxTree<T>::insert (const T &t)
{
   SxPtr<SxTree<T> > p2tree = SxPtr<SxTree<T> >::create (t);
   p2tree->parent = this;

   // Connecting the new Element to the left and right siblings
   p2tree->rightSibling = NULL;
   // creating a reference to the left sibling if one exits
   if (children.getSize () > 0)  {
      p2tree->leftSibling  = Iterator (&(*children.last ()));
      //p2tree->leftSibling  = Iterator (&(*(*--children.end ())));//STL
      children.last ()->rightSibling = Iterator (&(*p2tree));
      //(*--children.end ())->rightSibling = Iterator (&(*p2tree));//STL
   }  else  {
      p2tree->leftSibling = NULL;
   }

   children.append (p2tree);
   //children.insert (children.end (), p2tree);//STL
   ++nElements;

   return *p2tree;
}

template<class T>
SxTree<T> &SxTree<T>::insert (const SxTree<T> &tree)
{

   // Copying the tree
   SxPtr<SxTree<T> > pCopy = SxPtr<SxTree<T> >::create (tree);

   pCopy->parent = this;

   // Connecting the new Element to the left and right siblings
   pCopy->rightSibling = NULL;
   // creating a reference to the left sibling if one exits
   if (children.getSize () > 0)  {
      pCopy->leftSibling  = Iterator (&(*children.last ()));
      //pCopy->leftSibling  = Iterator (&(*(*--children.end ())));//STL
      children.last ()->rightSibling = Iterator (&(*pCopy));
      //(*--children.end ())->rightSibling = Iterator (&(*pCopy));//STL
   }  else  {
      pCopy->leftSibling = NULL;
   }

   // Adding the tree to the children list
   children.append(pCopy);
   //children.insert (children.end (), pCopy);//STL

   nElements += tree.nElements;

   return *this;
}

template<class T>
SxTree<T> &SxTree<T>::appendChild (const T &t)
{
   return insert (t);
}

template<class T>
SxTree<T> &SxTree<T>::appendTree (const SxTree<T> &tree)
{
   return insert (tree);
}

template<class T>
SxTree<T> & SxTree<T>::prependChild (const T & t)
{
   SxPtr<SxTree<T> > p2tree = SxPtr<SxTree<T> >::create ();
   p2tree->data   = SxPtr<T>::create (t);
   p2tree->nElements = 1;
   p2tree->parent = this;
   p2tree->leftSibling = NULL;
   p2tree->rightSibling = begin ();
   // Adding a leftSibling reference to the former start element
   if (begin () != Iterator (NULL))  {
      p2tree->rightSibling->leftSibling = &(*p2tree);
   }
   ++nElements;
   children.prepend (p2tree);
   //children.insert (children.begin (),p2tree);//STL

   return *p2tree;
}

template<class T>
SxTree<T> & SxTree<T>::prependTree (const SxTree<T> & tree){

   // Copying the tree
   SxPtr<SxTree<T> > pCopy = SxPtr<SxTree<T> >::create (tree);

   // Setting the parent Iterator of the copied tree
   pCopy->parent = this;

   // Setting the rightSibling Iterator of the copied tree
   pCopy->rightSibling = begin ();

   // Adding a leftSibling reference to the former start element which
   // references the copied tree
   if (begin () != Iterator (NULL))  {
      pCopy->rightSibling->leftSibling = &(*pCopy);
   }

   // Prepending the tree to the children list
   children.prepend (pCopy);

   // Updating the size
   nElements += tree.nElements;

   return *this;
}

// -- Output
template<class T>
std::ostream & SxTree<T>::print (std::ostream &os, int indent) const{

   for(int i = 0; i < indent-1; ++i)  { os << " ";  }
   os << "|";
   for(int i = 0; i < indent; ++i)  { os << " ";  }
   os << "+--";

   // testing if there is anything to print
   if (data)  os << *data << " (";
   typename SxTree<T>::ConstIterator helpIter;
   helpIter = getParent ();
   if (parent != Iterator (NULL))  os << " ^" << helpIter->getValue ();
   helpIter = getLeftSibling ();
   if (leftSibling != Iterator (NULL))  os << " <" << helpIter->getValue ();
   helpIter = getRightSibling ();
   if (rightSibling != Iterator (NULL))  os << " >" << helpIter->getValue ();
   os << ")" << endl;

   for(typename SxList<SxPtr< SxTree<T> > >::ConstIterator it
       = children.begin ();
       it != children.end ();
       ++it)  {
      (*(*it)).print (os, indent+1);
  }

   return os;
}

// -- Removing of Elements
template<class T>
void SxTree<T>::removeChildren ()
{
   // Saving the number of elements which will be removed temporarily
   unsigned int nRemovedElements = getSize ()-1;

   // Changing the total number of the children for the parents
   for (Iterator it = Iterator (this);
         it != Iterator (NULL);
         it = it->parent)  {
      it->nElements -= nRemovedElements;
   }

   // Removing all children
   //children.clear();//STL
   children.removeAll();
}

template<class T>
void SxTree<T>::removeChild (int index)
{

   SX_CHECK (index >= 0 && index < children.getSize (), 
   index, children.getSize ());

   // Creating an Iterator which references the object that will be deleted
   Iterator victimIter;
   victimIter = begin ();
   victimIter += index;

   // Disconnecting the references from the left and right sibling
   if (victimIter->leftSibling != Iterator (NULL))  {
      victimIter->leftSibling->rightSibling = victimIter->rightSibling;
   }
   if (victimIter->rightSibling != Iterator (NULL))  {
      victimIter->rightSibling->leftSibling = victimIter->leftSibling;
   }

   // Saving the number of elements which will be removed temporarily
   unsigned int nRemovedElements = victimIter->getSize ();

  // Changing the total number of the children for the parents
   for (Iterator it = Iterator (this);
         it != Iterator (NULL);
         it = it->parent)  {
      it->nElements -= nRemovedElements;
   }

   // Removing the child
   children.remove (index);
}

template<class T>
void SxTree<T>::remove (Iterator victimIter)
{

   SX_CHECK (victimIter->parent != Iterator (NULL));

   // Calculating the index
   Iterator it;
   int index = 0;
   for(it = victimIter->parent->begin (); it != victimIter; ++it){
      ++index;
   }

   // Disconnecting the references from the left and right sibling
   if (victimIter->leftSibling != Iterator (NULL))  {
      victimIter->leftSibling->rightSibling = victimIter->rightSibling;
   }
   if (victimIter->rightSibling != Iterator (NULL))  {
      victimIter->rightSibling->leftSibling = victimIter->leftSibling;
   }

   // Saving the number of elements which will be removed temporarily
   unsigned int nRemovedElements = victimIter->getSize ();

   // Changing the total number of the children for the parents
   for (it = victimIter->parent; it != Iterator (NULL); it = it->parent)  {
      it->nElements -= nRemovedElements;
   }

   // Removing the element
   victimIter->parent->children.remove (index);
}

template<class T>
typename SxTree<T>::Iterator SxTree<T>::find (const T &t)
{
   if (*data == t) return Iterator (this);

   for (Iterator it = begin (); it != end (); ++it)  {
      Iterator childIter = it->find (t);
      if (childIter != Iterator (NULL))  return childIter;
   }

   return NULL;
}

template<class T>
typename SxTree<T>::ConstIterator SxTree<T>::find (const T &t) const
{
   if (*data == t) return Iterator (this);

   for (Iterator it = begin (); it != end (); ++it)  {
      Iterator childIter = it->find (t);
      if (childIter != Iterator (NULL))  return childIter;
   }

   return NULL;
}

// -- Beginning and End of the tree
template<class T>
typename SxTree<T>::Iterator SxTree<T>::begin ()
{
   return getFirstChild ();
}

template<class T>
typename SxTree<T>::ConstIterator SxTree<T>::begin () const
{
   return getFirstChild ();
}

template<class T>
typename SxTree<T>::Iterator SxTree<T>::end ()
{
   return Iterator (NULL);
}

template<class T>
typename SxTree<T>::ConstIterator SxTree<T>::end () const
{
   return ConstIterator (NULL);
}

template<class T>
typename SxTree<T>::Iterator SxTree<T>::getFirstChild ()
{
   if (children.getSize() > 0)  {
      return Iterator (&(*(*children.begin ())));
   } else {
      return Iterator (NULL);
   }
}

template<class T>
typename SxTree<T>::ConstIterator SxTree<T>::getFirstChild () const
{
   if (children.getSize() > 0)  {
      return ConstIterator (&(*(*children.begin ())));
   } else {
      return ConstIterator (NULL);
   }
}

template<class T>
typename SxTree<T>::Iterator SxTree<T>::getLastChild ()
{
   if (children.getSize() > 0)  {
      return Iterator (&(*children.last ()));
      //return Iterator (&(*(*--children.end ())));//STL
   } else {
      return Iterator (NULL);
   }
}

template<class T>
typename SxTree<T>::ConstIterator SxTree<T>::getLastChild () const
{
   if (children.getSize () > 0)  {
      return ConstIterator (&(**children.fromLast ()));
      //return Iterator (&(*(*--children.end ())));//STL
   } else  {
      return ConstIterator (NULL);
   }
}

template<class T>
std::ostream & operator<< (std::ostream & os, const SxTree<T> &rhs){
   return rhs.print (os);
}

