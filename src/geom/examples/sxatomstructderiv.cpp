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

#include <SxAtomicStructure.h>
#include <SxRotation.h>
#include <SxOperator.h>
#include <SxNeighbors.h>

// filter: pick one species
class SxSpecFilter : public SxOperatorBase<SxAtomicStructure> SXOP_LINKFIX
{
   public:
      int iSpec;
      SxSpecFilter (int i) : iSpec(i) {}
      virtual SxAtomicStructure operator* (const SxAtomicStructure &in) const;
      virtual void applyInPlace (SxAtomicStructure &) const { SX_EXIT;}
      SXOPERATOR_GETCOPY (SxSpecFilter, SxAtomicStructure);
};

// filter: pick atoms below plane
class SxPlaneFilter : public SxOperator<SxAtomicStructure> SXOP_LINKFIX
{
   public:
      Coord planeNormal;
      double offset;
      SxPlaneFilter (Coord n, Coord p) : planeNormal(n)
      {
         planeNormal.normalize ();
         offset = planeNormal ^ p;
      }
      virtual SxAtomicStructure operator* (const SxAtomicStructure &in) const;
      virtual void applyInPlace (SxAtomicStructure &) const { SX_EXIT;}
      SXOPERATOR_GETCOPY (SxPlaneFilter,SxAtomicStructure);
};

int main ()  
{
   // define a cubic unit cell
   SxCell cell = CellMat (10., 0., 0.,
                          0., 10., 0.,
                          0., 0., 10.);

   
   // === SETTING UP AN ATOMIC STRUCTURE by hand  =============
   
   //     (see tools/SxStruct... for reading from files)
   SxAtomicStructure str(cell);
   SxList<SxString> chemNames;

   str.startCreation ();
   
   str.newSpecies ();
   chemNames << "A";
   
   // add to atoms to the structure (absolute, cartesian!)
   str.addAtom (Coord(0., 0., 0.));
   str.addAtom (Coord(5., 5., 5.));
   
   str.newSpecies ();
   chemNames << "B";
   
   str.addAtom (cell.relToCar (Coord(0.25, 0.25, 0.25)));
   str.addAtom (cell.relToCar (Coord(0.75, 0.75, 0.75)));

   // always after manual setup: end creation, update symmetries
   str.endCreation ();
   // str.updateSymmetries ();
   
   // define a pseudo-species data (for printing)
   SxSpeciesData specData(chemNames);

   cout << "Original structure:" << endl;
   str.print (specData);

   // === DERIVED STRUCTURES =============

   // --- derive a new structure (using SxPlaneFilter)
   SxPlaneFilter pf(Coord(0,0,1),    // normal
                    Coord(0,0,5.1)); // point on plane

   SxAtomicStructure derStr = pf | str;

   cout << "Derived structure:" << endl;
   derStr.print (specData);

   /** derived structures differ
       - in number of atoms
       - in atom ordering
       - or both
       but are still computationally compatible with their
       parent
   */

   // IMPORTANT: the fact that derStr is derived from str (and how)
   // is contained in atomInfo 
   SxAtomicStructure derDispl(derStr.getNAtoms (), derStr.atomInfo);
   derDispl.set(Coord(0.1,0.1,0.1));

   // change derStr
   derStr += derDispl;
   
   cout << "Derived partial structure with modifications:" << endl;
   derStr.print (specData);

   str <<= derStr; // update positions in str
   cout << "Original structure with partial modifications copied in:" << endl;
   str.print (specData);
  
   /* the following operations respect derived structures
      <<=   update
      +=    translation
      -=    translation (opposite direction)
   */
   

   // displacements
   str += derDispl;
   str -= derDispl;

   // derivation is recursive
   SxSpecFilter sf(1); // pick second species

   // derive doubleDerived from derStr
   // i.e.
   // str => derStr => doubleDerived
   // also doubleDerived is compatible to str
   SxAtomicStructure doubleDerived(sf | derStr);
   cout << "Derived partial structure doubleDerived <-- derStr <-- str:" 
        << endl;
   doubleDerived.print (specData);
   
   doubleDerived += Coord(1.,0.,0.);
   cout << "Derived partial structure doubleDerived modified" << endl;
   doubleDerived.print (specData);

   cout << "Original structure with partial modifications copied in:" << endl;
   str <<= doubleDerived; 

   str.print (specData);

   // ====== SOURCES OF DERIVED STRUCTURES
   
   // --- filters
   cout << "Filtered out first species" << endl;
   (SxSpecFilter(0) | str).print (specData);

   // --- neighbors
   SxGrid grid  (str, 10);
   SxNeighbors nn;
   Coord r0 (-7.,12.,9.);
   int nnMode = SxNeighbors::StoreAbs | SxNeighbors::IncludeZeroDistance;
   
   nn.compute (grid, str, r0, 12., nnMode);
   
   cout << "Atoms within 12 bohr of " << r0 << ":" << endl;
   nn.absPositions.print (specData);

   SxAtomicStructure delta (nn.absPositions, SxAtomicStructure::Copy);
   delta -= str; // subtract reference structure

   // now delta contains the supercell shift vectors
   for (int ia = 0; ia < nn.getSize (); ++ia)
      cout << nn.absPositions(ia)
           << " = str(" << (nn.idx(ia)+1) 
           << ") in cell " << cell.carToRel(delta(ia)) 
           << endl;
   
   // --- structure matching
   SxAtomicStructure newStr;
   newStr.startCreation ();
   newStr.newSpecies ();
   newStr.addAtom (str(0,0));
   newStr.newSpecies ();
   newStr.addAtom (str(1,0));
   newStr.endCreation ();

   // newStr is not derived from str, but some atom positions agree
   if (! newStr.atomInfo->isChild (str.atomInfo))
      cout << "newStr is not derived from str" << endl;

   // structure matching creates an atomInfo 
   SxAtomInfo::ConstPtr newAtomInfo;
   newAtomInfo = str.match (grid, newStr);
   newStr.replaceInfo (newAtomInfo);

   if (newStr.atomInfo->isChild (str.atomInfo))
      cout << "newStr is derived from str" << endl;

   
   
   
}


SxAtomicStructure SxSpecFilter::operator* (const SxAtomicStructure &in) const
{
   SxAtomInfo::Ptr aInfo = SxAtomInfo::derive (in.atomInfo);

   // new atomic structure: only atoms of 1 species
   int nAtom = in.getNAtoms(iSpec);
   
   // Keep the other species to keep species id
   aInfo->nSpecies = in.atomInfo->nSpecies;
   aInfo->nAtoms.resize (aInfo->nSpecies);
   aInfo->nAtoms.set (0);
   aInfo->nAtoms(iSpec) = nAtom;
   aInfo->setupOffset ();
   
   // create map to parent
   aInfo->parentMap.resize (nAtom);
   for (int ia = 0; ia < nAtom; ++ia)
      aInfo->parentMap(ia) = in.getIAtom(iSpec,ia);

   return SxAtomicStructure(aInfo, in);

}

SxAtomicStructure SxPlaneFilter::operator* (const SxAtomicStructure &in) const
{
   SxAtomInfo::Ptr aInfo = SxAtomInfo::derive (in.atomInfo);

   // new atomic structure: only 1 species
   aInfo->nSpecies = in.getNSpecies ();
   aInfo->nAtoms.resize (in.getNSpecies ());
   aInfo->nAtoms.set (0);
   
   // filter 
   SxList<int> atomsBelowPlane;
   for (int ia = 0; ia < in.getNAtoms (); ++ia)  {
      if ((planeNormal ^ in.constRef(ia)) < offset)  {
         atomsBelowPlane.append (ia);
         aInfo->nAtoms(in.getISpecies(ia))++;
      }
   }
   aInfo->parentMap = SxVector<Int>(atomsBelowPlane);
   aInfo->setupOffset ();

   return SxAtomicStructure(aInfo, in);

}

