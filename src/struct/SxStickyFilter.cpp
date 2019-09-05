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

#include <SxStickyFilter.h>
#include <SxList.h>

SxStickyFilter::SxStickyFilter ()
{
   // empty
}


SxStickyFilter::SxStickyFilter (const SxSymbolTable *structGroup)
{
   SxList<SxVector3<Int> > stickyList;
   bool movableX, movableY, movableZ;
   SxSymbolTable *speciesGroup, *atomGroup;
   bool anyMovable = false;
   try  {
      for (speciesGroup  = structGroup->getGroup("species");
           speciesGroup != NULL;
           speciesGroup  = speciesGroup->nextSibling ("species"))
      {
         for (atomGroup  = speciesGroup->getGroup("atom");
              atomGroup != NULL;
              atomGroup  = atomGroup->nextSibling("atom"))
         {
            movableX = movableY = movableZ = false;

            if (atomGroup->contains ("movable"))  {
               movableX =  atomGroup->get("movable")->toAttribute();
               movableY = movableZ = movableX;
               anyMovable = true;
            }
            if ( atomGroup->contains ("movableX"))  {
               movableX =  atomGroup->get("movableX")->toAttribute();
               anyMovable = true;
            }

            if ( atomGroup->contains ("movableY")) {
               movableY =  atomGroup->get("movableY")->toAttribute();
               anyMovable = true;
            }

            if ( atomGroup->contains ("movableZ")) {
               movableZ =  atomGroup->get("movableZ")->toAttribute();
               anyMovable = true;
            }

            stickyList << SxVector3<Int> (movableX, movableY, movableZ);
         }
      }
   }  catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   sticky = stickyList;
   if (!anyMovable)  {
      // if no movable specified, make all atoms movable
      sticky.set (SxVector3<Int> (true, true, true));
   }
}


SxStickyFilter::~SxStickyFilter ()
{
   // empty
}

SxVector<TPrecTauR>
SxStickyFilter::operator* (const SxVector<TPrecTauR> &in) const
{
   SX_CHECK (in.getSize() == 3*sticky.getSize(),
             in.getSize(),     sticky.getSize());

   ssize_t nAtoms = sticky.getSize();
   SxVector<TPrecTauR> res(3 * nAtoms);

   ssize_t i, d, iDoF;
   for (i=0,iDoF=0; i < nAtoms; i++)  {
      for (d=0; d < 3; d++, iDoF++)  {
         // sticky == 0 => coordinate is fixed
         res(iDoF) = sticky(i)(d) ? in(iDoF) : 0.;
      }
   }

   return res;
}

SxAtomicStructure SxStickyFilter::operator* (const SxAtomicStructure &str) const
{
   SxAtomicStructure res;
   res.cell = str.cell;
   res.atomInfo = str.atomInfo;
   res.set ((*this * str.coords).reshape (3, str.getNAtoms ()));
   return res;
}

void SxStickyFilter::validate (const SymMat &S,
                               const SxVector<Int> &equivalentIdx)
{
   ssize_t n = sticky.getSize ();
   SX_CHECK (n == equivalentIdx.getSize(), n, equivalentIdx.getSize() );

   for (ssize_t i = 0; i < 3; ++i)  {
      for (ssize_t j = 0; j < 3; ++j)  {
         // don't check coordinates that are not intermixed
         if (fabs(S(j,i)) < 1e-7) continue;
         for (ssize_t ia = 0; ia < n; ia++)  {
            if ((bool)sticky(ia)(i) != (bool)sticky(equivalentIdx(ia))(j))  {
               cout << "Stickyness validation error" << endl;
               SxString xyz("xyz");
               int ja = equivalentIdx(ia);
               if (   sticky(ia)(0) != sticky(ia)(1)
                   || sticky(ia)(0) != sticky(ia)(2))  {
                  cout << xyz(i) << "-coordinate of ";
               }
               cout << "atom " << (ia + 1) << " is "
                    << (sticky(ia)(i) ? "fixed" : "movable") << ", but" << endl;
               if (   sticky(ja)(0) != sticky(ja)(1)
                   || sticky(ja)(0) != sticky(ja)(2))  {
                  cout << xyz(j) << "-coordinate of ";
               }
               cout << "atom " << (ja + 1) << " is "
                    << (sticky(ja)(j) ? "fixed" : "movable") << ".\n";

               cout << "Symmetry equivalent "
                    << ((ia != ja) ?  "atoms" : "directions")
                    << " should have the same stickyness." << endl;
               cout << "Please adjust the \"movable\" flags in the structure "
                       "group.\n";
               SX_QUIT;
            }
         }
      }
   }
}

