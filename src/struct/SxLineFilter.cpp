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

#include <SxLineFilter.h>
#include <SxStack.h>

SxLineFilter::SxLineFilter ()
{
   // empty
}


SxLineFilter::SxLineFilter (const SxSymbolTable *structGroup,
                            const SxAtomicStructure &str)
{
   SxStack<int> lineFixIdx; // maps lineConstraint -> total atom idx
   int ia = 0;
   try  {
      SxSymbolTable *speciesGroup, *atomGroup;
      for (speciesGroup  = structGroup->getGroup("species");
           speciesGroup != NULL;
           speciesGroup  = speciesGroup->nextSibling ("species"))
      {
         for (atomGroup  = speciesGroup->getGroup("atom");
              atomGroup != NULL;
              atomGroup  = atomGroup->nextSibling("atom"), ++ia)
         {
            if (atomGroup->contains ("movableLine"))  {
               Coord line(atomGroup->get("movableLine")->toList ());
               if (line.normSqr () < 1e-12)  {
                  cout << "Invalid line constraint for atom " << (ia+1)
                       << ": line " << line << " is zero vector." << endl;
                  SX_QUIT;
               }
               line.normalize ();
               lines.addAtom (line);
               lineFixIdx << ia;
            }
         }
      }
   }  catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   lines.endCreation ();

   if (lines.getSize () == 0) return;

   // --- setup info for mapping
   SxPtr<SxAtomInfo> info = SxAtomInfo::derive (str.atomInfo);
   info->resize (1);
   info->nAtoms = lines.getSize ();
   info->setupOffset ();
   info->parentMap = lineFixIdx;
   lines.replaceInfo (info);

   // --- validation
   cout << "TODO: validation for line constraints" << endl;
}


SxLineFilter::~SxLineFilter ()
{
   // empty
}

SxAtomicStructure SxLineFilter::operator* (const SxAtomicStructure &str) const
{
   SX_CHECK (str.atomInfo == lines.atomInfo->getParent ());
   SxAtomicStructure res(str, SxAtomicStructure::Copy);
   for (int ic = 0; ic < lines.getSize (); ++ic)  {
      int ia = lines.atomInfo->parentMap(ic);
      // project vector onto constraint line
      cout << "Filtering " << ia << ": " << res.ref(ia);
      res.ref(ia) = (res.ref(ia) ^ lines(ic)) * lines(ic);
      cout << "->" << res.ref(ia) << endl;
   }
   return res;
}

