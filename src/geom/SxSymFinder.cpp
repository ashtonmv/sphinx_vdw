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

#include <SxSymFinder.h>


void SxSymFinder::compute (const SxAtomicStructure &structure)
{
   equations.removeAll ();

   SxList<SymMat> list = structure.cell.getLatticeSymmetries ();
   SxList<SymMat>::ConstIterator listIt = list.begin ();
   SxList<SymMat>::ConstIterator listEnd = list.end ();

   // ----- get affine symmetries, i.e. x -> Ux + t 
   SxRotation S;
   Coord t;
   SxAtomicStructure rotated;
   int ia;
   int count = 0;
   SymMat one (1., 0., 0.,
               0., 1., 0.,
               0., 0., 1.);

   for (; listIt != listEnd; ++listIt) {
      S = (*listIt);
      rotated = S ^ structure;
      for (ia = 0; ia < structure.getNAtoms(0); ia++)  {
         t = structure.getAtom(0,ia) - rotated.getAtom(0); 
         structure.cell.map (&t, SxCell::Origin);
         if (structure == (rotated + t))  {
            cout << "Found affine symmetry " << (++count) << endl;
            cout << S.getName () << " + " << t << endl;
            SxLinEquation eq(S - one, t);
            if (eq.isSoluble ())  {
               // --- affine symmetry can be made symmorphic
               Coord tNew;
               SxVector3<Int> nVec, nMax;
               if (eq.getDimension () < 3)  {
                  // a lattice translation in t may not result in a
                  // lattice translation for all intersections of this
                  // space with the others. Thus, we blow up the whole
                  // thing until we are sure that we got all these 
                  // intermediate positions
                  
                  // --- get space periodicity in cell basis directions
                  Coord point0, dist;
                  point0 = eq.getPoint ();
                  int loopCount = 0; // while could loop forever
                  for (int dim = 0; dim < 3; dim++)  {
                     nVec.set (0);
                     do  {
                        nVec(dim)++;
                        tNew = t + (structure.cell ^ Coord (nVec));
                        eq = SxLinEquation (S - one, tNew);
                        if (!eq.isSoluble ()) break;
                        dist = eq.getPoint () - point0;
                        structure.cell.map (&dist, SxCell::Origin);
                        if (dist.absSqr().sum() < 1e-5) break;
                     } while (loopCount < 1000);
                     nMax(dim) = nVec(dim);
                  }
                  if (loopCount >= 1000)  {
                     cout << "Too many spaces in elementary cell." << endl;
                     SX_EXIT;
                  }
               } else {
                  nMax.set (1);
               }
               // --- append all shifted spaces

               /*
               for (nVec(0) = -nMax(0); nVec(0) <= nMax(0); nVec(0)++)  {
                  for (nVec(1) = -nMax(1); nVec(1) <= nMax(1); nVec(1)++)  {
                     for (nVec(2) = -nMax(2); nVec(2) <= nMax(2); nVec(2)++)  {
                */
               for (nVec(0) = 0; nVec(0) < nMax(0); nVec(0)++)  {
                  for (nVec(1) = 0; nVec(1) < nMax(1); nVec(1)++)  {
                     for (nVec(2) = 0; nVec(2) < nMax(2); nVec(2)++)  {
                        
                        tNew = t + (structure.cell ^ Coord (nVec));
                        eq = SxLinEquation (S - one, tNew);
                        
                        if (!eq.isSoluble ()) continue;
                        if (eq.getDimension () == 0)  {
                           // substitute translation by mapped one
                           Coord point = eq.getPoint ();
                           structure.cell.map (&point, SxCell::Origin);
                           eq = SxLinEquation (one, point);
                        }
                        
                        // --- append current space
                        equations(eq) << S;

                        /*
                        // DEBUGGING
                        cout << "Can be made symmorphic by " << endl;
                        cout << eq.getName ("translation") << endl;
                         */
                     }
                  }
               }
            } // if eq soluble
         } // is affine symmetry
      } // translations
   } // lattice symmetries
            
   cout << SX_SEPARATOR;
   cout << "| Affine symmetry search complete." << endl;
   cout << "| " << equations.getSize () << " shift spaces found." << endl;
   cout << "| Checking combinations..." << endl;

   // --- get full set of equation combinations: try all (e1 && e2)
   
   SxMap<SxLinEquation, RotationList, SxNull>::Iterator it1;
   SxMap<SxLinEquation, RotationList, SxNull>::Iterator it2;
   SxMap<SxLinEquation, RotationList, SxNull>::Iterator itBegin;
   SxMap<SxLinEquation, RotationList, SxNull>::Iterator itEnd;

   it1 = equations.begin ();
   itBegin = equations.begin ();
   itEnd = equations.end ();
   
   SxLinEquation eq;
   
   for (; it1 != itEnd; ++it1)  {
      (it2 = it1)++;
      for (; it2 != itEnd; ++it2)  {
         eq = (it1.getKey() && it2.getKey());
         
         if (!eq.isSoluble ()) continue; 
         
         // --- check if eq is identical to one of the two
         if (eq == it1.getKey())  { 
            it1.getValue() << it2.getValue();
            continue;
         }
         if (eq == it2.getKey())  {
            it2.getValue() << it1.getValue();
            continue;
         }
         
         if (eq.getDimension () == 0)  {
            // substitute translation by mapped one
            Coord point = eq.getPoint ();
            structure.cell.map (&point, SxCell::Origin);
            eq = SxLinEquation (one, point);
         }
         (equations(eq) << it1.getValue()) << it2.getValue();
         
         // back & forth to update iterator
         it2--; it2++;
         if (it1 != itBegin)  {
            it1--; it1++;
         }
      } // 2nd equation
   } // 1st equation
      
   cout << SX_SEPARATOR;
   cout << "| " << equations.getSize ();
   cout << " shift spaces found after combination" << endl;
   cout << SX_SEPARATOR;
}

Coord SxSymFinder::getHighSymShift () const
{
   Coord translation(0.,0.,0.);
   ssize_t max = -1;
   
   SxMap<SxLinEquation, RotationList, SxNull>::ConstIterator it;
   SxMap<SxLinEquation, RotationList, SxNull>::ConstIterator itEnd;

   it = equations.begin ();
   itEnd = equations.end ();
   
   for (; it != itEnd; ++it)  {
      if (it.getValue().getSize () > max &&
          it.getKey().getDimension () < 3)
      {
         max = it.getValue().getSize ();
         translation = it.getKey().getPoint ();
      }
   }
   
   return translation;
}


