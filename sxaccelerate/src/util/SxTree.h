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

#ifndef _SX_TREE_H_
#define _SX_TREE_H_

// Including headerfiles
#include <SxPtr.h>
#include <SxList.h>
#include <SxString.h>
#include <iostream>
#include <fstream>

/** 
  
  \author Thomas Uchdorf, uchdorf@mpie.de

*/
// Tree-class
/** The Tree-class is a container which provides storage for random
  datatypes in the graphical structure of a tree.\n
  A tree is a certain kind of graph. Graphs consists of nodes and edges. Nodes
  contain data and references to other nodes. Connections between nodes created
  by the references represent edges.\n
  The characteristic attribute of a tree is that it has one root node of which
  every other node could be accessed. Each child node represents a tree of its
  own. The number of child nodes for this tree implementation is just limited
  by physical matters (i.e. memory). 
 
  \par Memory model
  This implementation of a tree generally uses shallow copy with reference
  counting by the SxPtr for copy-, assignment- and insertion-
  operations. That means when you replicate a tree the new tree is just linked
  to the memory of the old one. The memory of a tree is only discarded when
  there is no tree left using it. This guarantees that no superfluous memory is
  allocated.
  If you definitly want to occupy new memory for a copy which means you wish to
  make a deep copy, you can do that by using the copy constructor with a certain
  flag shown in the following example. 
  \par Example:
  - deep copy
  \code
    // Creating an empty tree of integers
    SxTree<int> tree;
    
    // Inserting elements to tree
    ...

    // Asuming the output of tree would look like this:
    // |+--0 ()
    // | +--1 ( ^0 >5)
    //  |  +--2 ( ^1 >3)
    //  |  +--3 ( ^1 <2 >4)
    //  |  +--4 ( ^1 <3)
    // | +--5 ( ^0 >6)
    // | +--6 ( ^0 <5 >7)
    // | +--7 ( ^0 <6)
    // Size: 8
 
    // Creating a deep copy of tree
    SxTree<int> copyTree(tree, SxTree<int>::deepCopy);

    // Removing element 1 of the deep copy
    SxTree<int>::Iterator deepCopyIter;
    deepCopyIter = copyTree.begin ();
    SxTree<int>::remove (deepCopyIter);

    // Result:
    // copyTree loses the sub tree with the root element 1 and other trees 
    // stay unaffected
  \endcode 
  - shallow Copy
  \code
    // Creating an empty tree of integers
    SxTree<int> tree;

    // Inserting elements to tree
    ...

    // Asuming the output of tree would look like this:
    // |+--0 ()
    // | +--1 ( ^0 >5)
    //  |  +--2 ( ^1 >3)
    //  |  +--3 ( ^1 <2 >4)
    //  |  +--4 ( ^1 <3)
    // | +--5 ( ^0 >6)
    // | +--6 ( ^0 <5 >7)
    // | +--7 ( ^0 <6)
    // Size: 8

    // Creating a shallow copy of tree
    SxTree<int> copyTree (tree);
    
    // Removing element 1 of the shallow copy
    SxTree<int>::Iterator shallowCopyIter;
    shallowCopyIter = copyTree.begin ();
    SxTree<int>::remove (shallowCopyIter);

    // Result:
    // shallowCopy loses the sub tree with the root node 1 and the other trees
    // stay unaffected   
  \endcode
  \sa SxTree::Iterator
  \sa SxTree::ConstIterator
*/ 
template<class T>
class SxTree
{
  public:
     /** \brief Specifies whether a shallow or deep copy should be used in
     the copy constructor (see SxTree::SxTree) */
     enum CopyType{ shallowCopy = 0, deepCopy};
     
  //Iterator-class
  /** The Iterator-class supplies the user with fundamental mechanisms for
    traversing a tree recursivly by a for-loop (see
    SxTree::Iterator::operator++ and SxTree::Iterator::operator--).
    Its objects represent references to trees. 
    \sa SxTree
    \sa SxTree::ConstIterator
  */
  class Iterator
  {
     public:

        // Constructors      
        // --Constructor
        /** \brief Creates an Iterator which references NULL */
        Iterator ();

        // -- Copyconstructor
        /** \brief Creates an Iterator which resembles another Iterator */
        Iterator (const typename SxTree<T>::Iterator &rhs);
        /** \brief Creates an Iterator from a pointer */
        Iterator (SxTree<T> *rhs);

        // Operations
        // -- Assignment
        /** \brief Assigns Iterators (just addresses!) */
        typename SxTree<T>::Iterator &operator= (const typename SxTree<T>::Iterator &rhs);
        /** \brief Assigns pointer to Iterator (just addresses!) */
        Iterator &operator= (SxTree<T> *rhs);

        // -- Equality        
        /** \brief Tests equality of two Iterators (just addresses!) */
        bool operator== (const typename SxTree<T>::Iterator &rhs) const;
        /** \brief Tests inequality of two Iterators (just addresses!) */
        bool operator!= (const typename SxTree<T>::Iterator &rhs) const;

        // -- Dereferencation
        /** \brief Dereferences with * */
        SxTree<T> &operator* ();
        /** \brief Dereferences with -> */
        SxTree<T> * operator-> ();

        // -- In- and Decrementation
        /** \brief Increments (prefix) 
         \par Example:
         \code
         // --- Traversing the entire tree recursivly by a for-loop
         void func (SxTree<T> &tree)
         {
           doSomethingWithTheTree (tree);
         
           for (SxTree<T>::Iterator it = tree.begin ();
                it != tree.end ();
                ++it)  {
             func (*it);
           }
         }
         
         int main (int argc, char *argv[])
         {
           // Create a tree of integers
           SxTree<int> rootNode;
 
           // Fill the tree with elements
           ...

           // Calling the function starting at the rootNode
           func (rootNode);

           ...
         }
         \endcode
         */
        typename SxTree<T>::Iterator &operator++ ();
        /** \brief Increments (postfix) */
        typename SxTree<T>::Iterator &operator++ (int);
        /** \brief Increments with a certain step width */
        void operator+= (unsigned int steps);
        /** \brief Decrements (prefix)
         \par Example:
         \code
         // --- Traversing the entire tree recursivly by a for-loop in reverse
         // order
         void func (SxTree<T> &tree)
         {
           doSomethingWithTheTree (tree);
         
           for (SxTree<T>::Iterator it = tree.getLastChild ();
                it != SxTree<T>::Iterator (NULL);
                --it)  {
             func (*it);
           }
         }

         int main (int argc, char *argv[])
         {
           // Create a tree of integers
           SxTree<int> rootNode;
 
           // Fill the tree with elements
           ...

           // Calling the function starting at the rootNode
           func (rootNode);

           ...
         }
         \endcode
         */        Iterator &operator-- ();
        /** \brief Decrements (postfix) */
        typename SxTree<T>::Iterator &operator-- (int);      
        /** \brief Decrements with a certain step width */
        void operator-= (unsigned int steps);
        /** \brief Returns a pointer to the constant tree */
        const SxTree<T> *getConstPtr() const { return (const SxTree<T> *) ptr; }
        
     protected:

        // Members
        /** \brief Points to a tree */
        SxTree<T> *ptr;
        
     private:  
  };
  
  //ConstIterator-class
  /** The ConstIterator-class supplies the user with fundamental mechanisms for
    traversing a constant tree recursivly by a for-loop (see 
    SxTree::ConstIterator::operator++ and SxTree::ConstIterator::operator--).
    Its objects represent references to unalienable trees.
  */
  class ConstIterator
  {
     public:

        // Constructors      
        // --Constructor
        /** \brief Creates a ConstIterator which references NULL */
        ConstIterator ();

        // -- Copyconstructor
        /** \brief Creates a ConstIterator which resembles another
        ConstIterator */
        ConstIterator (const typename SxTree<T>::ConstIterator &rhs);
        /** \brief Creates a ConstIterator which resembles another Iterator */
        ConstIterator (const typename SxTree<T>::Iterator &rhs);
        /** \brief Creates a ConstIterator from a pointer */
        ConstIterator (const SxTree<T> *rhs);

        // Operations
        // -- Assignment
        /** \brief Assigns ConstIterators (just addresses!) */
        typename SxTree<T>::ConstIterator &operator= (const typename SxTree<T>::ConstIterator &rhs);
        /** \brief Assigns pointer to ConstIterator (just addresses!) */
        typename SxTree<T>::ConstIterator &operator= (const SxTree<T> *rhs);

        // -- Equality        
        /** \brief Tests equality of two ConstIterators (just addresses!) */
        bool operator== (const typename SxTree<T>::ConstIterator &rhs) const;
        /** \brief Tests inequality of two ConstIterators (just addresses!) */
        bool operator!= (const typename SxTree<T>::ConstIterator &rhs) const;

        // -- Dereferencation
        /** \brief Dereferences with * */
        const SxTree<T> & operator* ();
        /** \brief Dereferences with -> */
        const SxTree<T> * operator-> ();

        // -- In- and Decrementation
        /** \brief Increments (prefix) 
         \par Example:
         \code
         // --- Traversing the entire tree recursivly by a for-loop
         void func (const SxTree<T> &tree)
         {
           doSomethingWithTheTree (tree);
         
           for (SxTree<T>::ConstIterator it = tree.begin ();
                it != tree.end ();
                ++it)  {
             func (*it);
           }
         }
         
         int main (int argc, char *argv[])
         {
           // Create a tree of integers
           SxTree<int> rootNode;
 
           // Fill the tree with elements
           ...

           // Calling the function starting at the rootNode
           func (rootNode);

           ...
         }
         \endcode
         */
        typename SxTree<T>::ConstIterator &operator++ ();
        /** \brief Increments (postfix) */
        typename SxTree<T>::ConstIterator &operator++ (int);
        /** \brief Increments with a certain step width */
        void operator+= (unsigned int steps);
        /** \brief Decrements (prefix)
         \par Example:
         \code
         // --- Traversing the entire tree recursivly by a for-loop in reverse
         // order
         void func (const SxTree<T> &tree)
         {
           doSomethingWithTheTree (tree);
         
           for (SxTree<T>::ConstIterator it = tree.getLastChild ();
                it != SxTree<T>::ConstIterator (NULL);
                --it)  {
             func (*it);
           }
         }
         
         int main (int argc, char *argv[])
         {
           // Create a tree of integers
           SxTree<int> rootNode;
 
           // Fill the tree with elements
           ...

           // Calling the function starting at the rootNode
           func (rootNode);

           ...
         }
         \endcode
         */
        typename SxTree<T>::ConstIterator &operator-- ();
        /** \brief Decrements (postfix) */
        typename SxTree<T>::ConstIterator &operator-- (int);      
        /** \brief Decrements with a certain step width */
        void operator-= (unsigned int steps);

     protected:

        // Members
        /** \brief Points to a constant tree */
        const SxTree<T> *ptr;
        
     private:  
  };
  
  // Constructors
  // -- Constructor
  /** \brief Creates a tree object */
  /** \par Example:
    \code
    // Creating an empty tree of integers
    SxTree<int> treeOfInt;
    
    // --- Output of the empty tree and its size
    std::cout << treeOfInt;
    std::cout << "Size: " << treeOfInt.getSize () << endl;
    // prints:
    // |+--0 ()
    // Size: 1
    \endcode
  */
  SxTree ();

  // -- Copy constructor
  /** \brief Creates a tree object reflecting an existing one */
  /** \par Example:
    \code
    // Creating an empty tree of integers and a sub tree with root element 1
    SxTree<int> treeOfInt, subTreeOfInt(1);

    // Inserting elements
    subTreeOfInt << 2 << 3 << 4;
    treeOfInt << subTreeOfInt << 5 << 6 << 7;

    // Creating a deep copy
    SxTree<int> deepCopyOfTreeOfInt (treeOfInt, SxTree<int>::deepCopy);

    // Creating a shallow copy
    SxTree<int> shallowCopyOfTreeOfInt (treeOfInt, SxTree<int>::shallowCopy);

    // Simpler way of creating a shallow copy
    SxTree<int> anotherShallowCopyOfTreeOfInt (treeOfInt);
    
    // --- Output of all the trees and their sizes
    std::cout << treeOfInt;
    std::cout << "Size: " << treeOfInt.getSize () << endl;
    std::cout << deepCopyOfTreeOfInt;
    std::cout << "Size: " << deepCopyOfTreeOfInt.getSize () << endl;
    std::cout << shallowCopyOfTreeOfInt;
    std::cout << "Size: " << shallowCopyOfTreeOfInt.getSize () << endl;   

    // prints three times:
    // |+--0 ()
    // | +--1 ( ^0 >5)
    // |  +--2 ( ^1 >3)
    // |  +--3 ( ^1 <2 >4)
    // |  +--4 ( ^1 <3)
    // | +--5 ( ^0 <1 >6)
    // | +--6 ( ^0 <5 >7)
    // | +--7 ( ^0 <6)
    // Size: 8
    \endcode
  */ 
  SxTree (const SxTree<T> &rhs, CopyType ct = shallowCopy);

  // -- Default constructor
  /** \brief Creates a tree object with a certain value in the root node */ 
  /** \par Example:
    \code
    // Creating a tree of integers with 1 as root element
    SxTree<int> treeOfInt(1);
    
    // --- Output of treeOfInt and its size
    std::cout << treeOfInt;
    std::cout << "Size: " << treeOfInt.getSize () << endl;
    // prints:
    // |+--1 ()
    // Size: 1
    \endcode
  */
  SxTree (const T &t);

  // -- Destructor
  /** \brief Destroys a tree object */
  ~SxTree ();
  
  // Operations
  // -- Assignment
  /** \brief Assigns a tree to another tree */
  /** \par Example:
    \code
    // Asuming our tree of integers would look like this and would be called
    // treeOfInt 
    // |+--0 ()
    // | +--1 ( ^0 >5)
    // |  +--2 ( ^1 >3)
    // |  +--3 ( ^1 <2 >4)
    // |  +--4 ( ^1 <3)
    // | +--5 ( ^0 <1 >6)
    // | +--6 ( ^0 <5 >7)
    // | +--7 ( ^0 <6)
    // Size: 8

    // Asuming another tree of integers would look like this and would be called
    // anotherTreeOfInt
    // |+--9 ()
    // | +--10 ( ^9)
    // |  +--11 ( ^10)
    // Size: 3
    
    // Assigning treeOfInt anotherTreeOfInt
    treeOfInt = anotherTreeOfInt;

    // Putting out treeOfInt
    std::cout << treeOfInt;
    std::cout << "Size: " << treeOfInt.getSize () << endl;
    
    // Prints:
    // |+--9 ()
    // | +--10 ( ^9)
    // |  +--11 ( ^10)
    // Size: 3
    \endcode
  */
  SxTree<T> &operator= (const SxTree<T> &rhs);
  /** \brief Assigns a value to the curent node of the tree */
  /** \par Example:
    \code
    // Asuming our tree of integers would look like this and would be called
    // treeOfInt 
    // |+--0 ()
    // | +--1 ( ^0 >5)
    // |  +--2 ( ^1 >3)
    // |  +--3 ( ^1 <2 >4)
    // |  +--4 ( ^1 <3)
    // | +--5 ( ^0 <1 >6)
    // | +--6 ( ^0 <5 >7)
    // | +--7 ( ^0 <6)
    // Size: 8

    // Assigning treeOfInt the element 23
    treeOfInt = 23;

    // --- Output of treeOfInt
    std::cout << treeOfInt;
    std::cout << "Size: " << treeOfInt.getSize () << endl;
    
    // Prints:
    // |+--23 ()
    // | +--1 ( ^23 >5)
    // |  +--2 ( ^1 >3)
    // |  +--3 ( ^1 <2 >4)
    // |  +--4 ( ^1 <3)
    // | +--5 ( ^23 <1 >6)
    // | +--6 ( ^23 <5 >7)
    // | +--7 ( ^23 <6)
    // Size: 8
    
    \endcode
  */ 
  SxTree<T> &operator= (const T &t);

  // -- Input
  /** \brief Inserts a new child */
  /** Works like SxTree::insert. */
  SxTree<T> &operator<< (const T &t);
  /** \brief Inserts a new tree */  
  /** Works like SxTree::insert. */
  SxTree<T> &operator<< (const SxTree<T> &tree);

  // Methods
  // -- Accessfunctions
  /** \brief Sets the value of the current node */
  /** Works like SxTree::operator=. */
  SxTree<T> &setValue (const T &t);
  /** \brief Gets the value of the current node */
  inline T &getValue ()  { return *data; }
  /** \brief Gets the value of the current node */
  inline const T &getValue () const  { return *data; }
  /** \brief Gets an Iterator which references the parent node */
  inline Iterator getParent ()  { return parent; }
  /** \brief Gets an ConstIterator which references the parent node */
  inline ConstIterator getParent () const  { return ConstIterator (parent); }
  /** \brief Gets an Iterator which references the left sibling node */
  inline Iterator getLeftSibling ()  { return leftSibling; }
  /** \brief Gets an ConstIterator which references the left sibling node */
  inline ConstIterator getLeftSibling () const  { return leftSibling; }
  /** \brief Gets an Iterator which references the right sibling node */
  inline Iterator getRightSibling ()  { return rightSibling; }
  /** \brief Gets an ConstIterator which references the right sibling node */
  inline ConstIterator getRightSibling () const  { return rightSibling; }
  /** \brief Gets the total amount of children, grandchildren and so on */
  inline unsigned int getNElements ()  { return nElements; }
  /** \brief Gets the total amount of children, grandchildren and so on */
  inline unsigned int getNElements () const  { return nElements; }
  /** \brief Gets the total amount of children, grandchildren and so on */
  /** Works like SxTree::getNElements. */
  inline unsigned int getSize ()  { return getNElements (); }
  /** \brief Gets the total amount of children, grandchildren and so on */
  /** Works like SxTree::getNElements. */
  inline unsigned int getSize () const  { return getNElements (); }
  
  // -- Insertion
  /** \brief Inserts a new child */
  /** \par Example:
    \code
    // Creating an empty tree of integers
    SxTree<int> treeOfInt;
    
    // Inserting the element 1
    treeOfInt.insert (1);
    
    // --- Output of treeOfInt and its size
    std::cout << treeOfInt;
    std::cout << "Size: " << treeOfInt.getSize () << endl;
    
    // prints:
    // |+--0 ()
    // | +--1 ( ^0)
    // Size: 2
    \endcode
  */
  SxTree<T> & insert (const T &t);
  /** \brief Inserts a new tree */
  /** \par Example:
    \code
    // Creating an empty tree of integers and a sub tree with root element 1
    SxTree<int> treeOfInt, subTreeOfInt(1);

    // Inserting elements to subTreeOfInt
    subTreeOfInt << 2 << 3 << 4;

    // Inserting subTreeOfInt
    treeOfInt.insert (subTreeOfInt);

    // Inserting more elements to treeOfInt
    treeOfInt << 5 << 6 << 7;
    
    // --- Output of treeOfInt and its size
    std::cout << treeOfInt;
    std::cout << "Size: " << treeOfInt.getSize () << endl;

    // prints:
    // |+--0 ()
    // | +--1 ( ^0 >5)
    // |  +--2 ( ^1 >3)
    // |  +--3 ( ^1 <2 >4)
    // |  +--4 ( ^1 <3)
    // | +--5 ( ^0 <1 >6)
    // | +--6 ( ^0 <5 >7)
    // | +--7 ( ^0 <6)
    // Size: 8
    \endcode
  */
  SxTree<T> & insert (const SxTree<T> &tree);
  /** \brief Inserts a child node at the end of the children list */
  /** Functions as SxTree::insert. */
  SxTree<T> & appendChild (const T &t);
  /** \brief Inserts a child tree at the end of the children list */
  /** Functions as SxTree::insert. */
  SxTree<T> & appendTree (const SxTree<T> &tree);
  /** \brief Prepends a new child to the tree */
  /** \par Example:
    \code
    // Creating an empty tree of integers
    SxTree<int> treeOfInt;
    
    // Inserting the elements 1, 2, 3
    treeOfInt << 2 << 3 << 4;
    
    // Prepending element 1
    treeOfInt.prependChild (1);
    
    // --- Output of the empty tree and its size
    std::cout << treeOfInt;
    std::cout << "Size: " << treeOfInt.getSize () << endl;
    
    // prints:
    // |+--0 ()
    // | +--1 ( ^0 >2)
    // | +--2 ( ^0 <1 >3)
    // | +--3 ( ^0 <2 >4)
    // | +--4 ( ^0 <3)
    // Size: 5
    \endcode
  */ 
  SxTree<T> & prependChild (const T &t);
  /** \brief Prepends a new child tree to the tree */
  /** \par Example:
    \code
    // Creating an empty tree of integers and a sub tree with root element 1
    SxTree<int> treeOfInt, subTreeOfInt(1);

    // Inserting elements
    subTreeOfInt << 2 << 3 << 4;
    treeOfInt << 5 << 6 << 7;

    // Prepending subTreeOfInt
    treeOfInt.prependTree (subTreeOfInt); 
    
    // --- Output of treeOfInt and its size
    std::cout << treeOfInt;
    std::cout << "Size: " << treeOfInt.getSize () << endl;

    // prints:
    // |+--0 ()
    // | +--1 ( ^0 >5)
    // |  +--2 ( ^1 >3)
    // |  +--3 ( ^1 <2 >4)
    // |  +--4 ( ^1 <3)
    // | +--5 ( ^0 <1 >6)
    // | +--6 ( ^0 <5 >7)
    // | +--7 ( ^0 <6)
    // Size: 8
    \endcode
  */
  SxTree<T> & prependTree (const SxTree<T> & tree);

  // -- Output
  /** \brief Prints the tree textually */
  /** Functions as ::operator<<. */
  std::ostream &print (std::ostream &os = std::cout, int indent = 0) const;

  // -- Removing of Elements
  /** \brief Removes all the children of the tree
  \par Example:
  \code
    ...
    // printing our example integer tree
    std::cout << exampleTree;
     
    // prints
    // |+--0
    // | +--1
    // |  +--2
    // |  +--3
    // |  +--4
    // | +--5
    // | +--6
    // |  +--7
     
    // Removing all the children of exampleTree
    exampleTree.removeChildren ();

    // printing our example integer tree
    std::cout << exampleTree;
    
    // prints
    // |+--0
  \endcode
  */
  void removeChildren ();
  /** \brief Removes a certain child of the tree specified by an index
  Be careful with the index. It starts with 0 and using an out-of-bounds-index
  leads to abortion of the program. 
  \par Example:
  \code
    ...
    // printing our example integer tree
    std::cout << exampleTree;
     
    // prints
    // |+--0
    // | +--1
    // |  +--2
    // |  +--3
    // |  +--4
    // | +--5
    // | +--6
    // |  +--7
     
    // Removing the sub tree with the root element 1
    exampleTree.removeChild (0);

    // printing our example integer tree
    std::cout << exampleTree;
    
    // prints
    // |+--0
    // | +--5
    // | +--6
    // |  +--7    
 
    // out-of-bounds-index! ==> program aborted
    exampleTree.removeChild (2);
  \endcode
  */
  void removeChild (int index);
  /** \brief Removes a certain element of the tree specified by an Iterator.
  Do not forget that root nodes can not be deleted
  \par Example:
  \code
	  ...
     // printing our example integer tree
     cout << exampleTree;

     // prints
     // |+--0
     // | +--1
     // |  +--2
     // |  +--3
     // |  +--4
     // | +--5
     // | +--6
     // |  +--7
     
    // causing an error by trying to remove rootElement
    SxTree<int>::Iterator rootIter (&exampleTree);

    SxTree<int>::remove (rootIter); // program aborted
  
    // removing element 2
    SxTree<int>::Iterator elemTwoIter;
    elemTwoIter = (exampleTree.begin ()).begin ();

    SxTree<int>::remove (elemTwoIter); // no problems
  
    // printing our example integer tree
    cout << exampleTree;
    
    // prints
    // |+--0
    // | +--1
    // |  +--3
    // |  +--4
    // | +--5
    // | +--6
    // |  +--7
  \endcode
  */
  static void remove (Iterator victimIter);

  // -- Search
  /** \brief Returns an Iterator referencing the first tree which contains the 
  searched data */
  Iterator find (const T &t);
  /** \brief Returns a ConstIterator referencing the first tree which contains 
  the searched data */  
  ConstIterator find (const T &t) const;
  /** \brief Tests whether the tree contains a certain data object or not */
  bool contains (const T &t)  { return find (t) != Iterator (NULL); }
 
  // -- Beginning and End of the tree
  /** \brief Returns an Iterator which references the first child */
  /** Functions as SxTree::getFirstChild. */
  Iterator begin ();
  /** \brief Returns a ConstIterator which references the first child */
  /** Functions as SxTree::getFirstChild. */
  ConstIterator begin () const;
  /** \brief Returns an Iterator which references the element after the last
  element */
  Iterator end ();
  /** \brief Returns a ConstIterator which references the element after the
  last element */
  ConstIterator end () const;
  /** \brief Returns an Iterator which references the first child */
  /** \par Example:
    \code
    ...
    // Asuming a tree of integers would look like this and secondElemIter
    // would reference element 1
    // |+--0
    // | +--1
    // |  +--2
    // |  +--3
    // |  +--4
    // | +--5
    // | +--6
    // |  +--7

    // Creating an Iterator which also references element 1
    Iterator firstChildIter;
    firstChildIter = secondElemIter;

    // Getting an Iterator referencing the first child of the sub tree with the
    // root element 1 
    firstChildIter.getFirstChild ();

    // Now firstChildIter references 2

    // printing the tree which is referenced by firstChildIter
    std::cout << *firstChildIter;

    // prints:
    // |+--2
   \endcode 
  */ 
  Iterator getFirstChild ();
  /** \brief Returns a ConstIterator which references the first child */
  /** \par Example:
    \code
    ...
    // Asuming a constant tree of integers would look like this and
    // secondElemIter would reference element 1
    // |+--0
    // | +--1    
    // |  +--2
    // |  +--3
    // |  +--4
    // | +--5
    // | +--6
    // |  +--7

    // Creating a ConstIterator which also references element 1
    ConstIterator firstChildIter;
    firstChildIter = secondElemIter;

    // Getting a ConstIterator referencing the first child of the sub tree with
    // the root element 1 
    firstChildIter.getFirstChild ();

    // Now firstChildIter references 2

    // printing the tree which is referenced by firstChildIter
    std::cout << *firstChildIter;

    // prints:
    // |+--2
   \endcode 
  */
  ConstIterator getFirstChild () const;
  /** \brief Returns an Iterator which references the last child */
  /** \par Example:
    \code
    ...
    // Asuming a tree of integers would look like this and secondElemIter
    // would reference element 1
    // |+--0
    // | +--1   
    // |  +--2
    // |  +--3
    // |  +--4
    // | +--5
    // | +--6
    // |  +--7

    // Creating an Iterator which also references element 1
    Iterator lastChildIter;
    lastChildIter = secondElemIter;

    // Getting an Iterator referencing the last child of the sub tree with the
    // root element 1 
    lastChildIter.getLastChild ();

    // Now lastChildIter references 4

    // printing the tree which is referenced by lastChildIter
    std::cout << *lastChildIter;

    // prints:
    // |+--4
   \endcode 
  */
  Iterator getLastChild ();
  /** \brief Returns a ConstIterator which references the last child */
  /** \par Example:
    \code
    ...
    // Asuming a const tree of integers would look like this and secondElemIter
    // would reference element 1
    // |+--0
    // | +--1   
    // |  +--2
    // |  +--3
    // |  +--4
    // | +--5
    // | +--6
    // |  +--7

    // Creating a ConstIterator which also references element 1
    ConstIterator lastChildIter;
    lastChildIter = secondElemIter;

    // Getting a ConstIterator referencing the last child of the sub tree with
    // the root element 1
    lastChildIter.getLastChild ();

    // Now lastChildIter references 4

    // printing the tree which is referenced by lastChildIter
    std::cout << *lastChildIter;

    // prints:
    // |+--4
   \endcode 
  */ 
  ConstIterator getLastChild () const;
 
  protected:

  // Methods
  void createDeepCopy (const SxTree<T> &rhs);
  void createShallowCopy (const SxTree<T> &rhs);
  
  // Members
  /** \brief List of AutoPointers which reference trees */
  SxList<SxPtr<SxTree<T> > > children;
  /** \brief Data object stored with an AutoPointer */
  SxPtr<T> data;
  /** \brief Total number of children, grandchildren and so on */
  unsigned int nElements; 
  /** \brief References the parent node */
  Iterator parent;
  /** \brief References the left node */
  Iterator leftSibling;
  /** \brief References the right node */
  Iterator rightSibling;
};

// Input-/Outputoperators
/** \brief Prints the tree textually by using the public member function 
SxTree::print */
template<class T>
std::ostream &operator<< (std::ostream & os, const SxTree<T> &rhs);

#include <SxTree.hpp>

#endif /* _SX_TREE_H_ */
