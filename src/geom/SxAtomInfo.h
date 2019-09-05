// ---------------------------------------------------------------------------
//
//      The ab-initio based multiscale library
//
//                  S / P H I / n X
//
//      Copyright:  Max-Planck-Institute for Iron Research
//                  40237 Duesseldorf, Germany
//
//      Contact:    https://sxlib.mpie.de
//      Authors:    see sphinx/AUTHORS
//      License:    see sphinx/LICENSE
//
// ---------------------------------------------------------------------------

#ifndef _SX_ATOM_INFO_H_
#define _SX_ATOM_INFO_H_

#include <SxPtr.h>
#include <SxVector.h>
#include <SxGeom.h>
#include <SxMetadata.h>

/** \brief Info on atomic species

    \b SxClass = S/PHI/nX 

    \author C. Freysoldt */
class SX_EXPORT_GEOM SxAtomInfo 
{
   public:
      /// Number of atoms per species
      SxVector<Int> nAtoms;
      /// Number of species
      int nSpecies;

      /// Metadata container
      mutable SxMetadata meta;
      
      /// Starting index for each species
      SxVector<Int> offset;

      /// Setup routine
      void set(const SxVector<Int> &nAtomsIn);
      
      /// Setup routine for offsets 
      void setupOffset ();

      /// Destructor
      inline ~SxAtomInfo ();

   protected:
      /// Constructor
      SxAtomInfo(int nSpeciesIn = -1) 
         : nAtoms(nSpeciesIn > 0 ? nSpeciesIn : 0), nSpecies(nSpeciesIn)
      { /* empty */ }
    
   public:
      /** \brief Set number of species
        @param newNSpecies new size
        @param keep        if true, keep nAtoms for current species 
                           (set to 0 otherwise)
        */
      void resize (int newNSpecies, bool keep = false);
      
      /// Get index range of atoms
      SxIdx getRange (int iSpecies) const
      {
         SX_CHECK (iSpecies < nSpecies && iSpecies >= 0,
                   iSpecies, nSpecies);
         int i0 = offset(iSpecies), i1 = i0 + nAtoms(iSpecies) - 1;
         SX_CHECK(i1 >= i0, i1, i0);
         return SxIdx(i0,i1);
      }

      /// Comparison (test for equivalence)
      bool operator== (const SxAtomInfo &in) const;
      
      /// Comparison (test for equivalence)
      bool operator!= (const SxAtomInfo &in) const
      {
         return ! operator==(in);
      }

      /// Allow SxPtr to create object
      friend class SxPtr<SxAtomInfo>;
      friend class SxConstPtr<SxAtomInfo>;

      /// Typedefs
      typedef SxConstPtr<SxAtomInfo> ConstPtr;
      typedef SxPtr<SxAtomInfo> Ptr;

      /// Create SxPtr<SxAtomInfo>
      inline static Ptr create ();
   
   private:
      /// SxAtomInfo from which this one is derived
      ConstPtr parent;
   public:
      /// Get parent 
      ConstPtr getParent () const { return parent; }

      /// Atom index mapping
      SxVector<Int> parentMap;

      /// Create a derived SxAtomInfo
      static Ptr derive(const ConstPtr &from);

      /** \brief Create a copy with different parent
        @param newParent parent of the new info
        @return the copy

        \note The newParent must be a grandparent of this info or NULL.

        The function's purpose is to shorten derivation chains
        (e.g. A->B->C->-D->E to A->E) or to orphanize a derived
        structure using the NULL ptr.
        */
      Ptr getCopy (const SxConstPtr<SxAtomInfo> &newParent) const;

      /// Check if we are derived from in
      bool isChild(const ConstPtr &in) const;
      
      /** \brief Get index map
          @param  in an SxAtomInfo from which this is somehow derived
          @return map of local index to parent index
        */
      SxVector<Int> getIdxMap(const ConstPtr &in) const;
};

SxAtomInfo::Ptr SxAtomInfo::create ()  
{ 
   return Ptr::create (); 
}

SxAtomInfo::~SxAtomInfo ()
{
#ifndef NDEBUG
   ConstPtr ptr = parent;
   while (ptr)  {
      // check for self-reference loops
      SX_CHECK(ptr.getPtr () != this);
      ptr = ptr->parent;
   }
#endif
}

#endif /* _SX_ATOM_INFO_H_ */
