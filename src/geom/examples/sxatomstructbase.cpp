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
   
   // important: always end creation mode before doing anything else!
   str.endCreation ();

   // compute and store symmetries
   str.updateSymmetries ();

   // define a pseudo-species data (for printing)
   SxSpeciesData specData(chemNames);

   str.print (specData);

   // --- add a second species (normally without intermediate end/startCreation)
   str.startCreation ();
   str.newSpecies ();
   
   chemNames << "B";
   
   // let's use relative coordinates this time
   str.addAtom (cell.relToCar (Coord(0.25, 0.25, 0.25)));
   str.addAtom (cell.relToCar (Coord(0.75, 0.75, 0.75)));

   // always after manual setup: end creation, update symmetries
   str.endCreation ();
   str.updateSymmetries ();
         
   specData = SxSpeciesData (chemNames);

   str.print (specData);



   
   // === ACCESSING SINGLE ATOMS =================

   cout << "First atom 2nd species @ " << str(1,0) << endl; 
   cout << "Fourth atom whatsoever species @ " <<  str(3) << endl;

   int is = 0, ia = 1;
   Coord orig = str(is, ia);
   
   // Now let's modify the second atom
   str.ref(is, ia) = Coord(0., 0., 1.); // only ref() may be modified
   str.ref(is, ia) *= 2.;               // or +=, -=, /=, ...
   cout << "str(is,ia)=" 
        << str.getAtom (is, ia) // same as str(is,ia);
        << endl; 
   // back to orig
   str.setAtom (is, ia, orig);  // same as ref(is,ia) = orig;
   

   
   
   // === BASIC TRANSFORMATIONS (GLOBAL) =============
   
   // --- translations 
   SxAtomicStructure tStr;
   
   Coord t = cell.relToCar(Coord(0.5,0.5,0.5));
   tStr = str + t; // shift structure
   tStr = str - t; // shift structure in opposite direction
   
   // ---  rotations
   SymMat S (0., -1.,  0.,
             1.,  0.,  0.,
             0.,  0.,  1.);
   cout << "structure rotated by " << SxRotation(S).getName () << endl;
   cout << SEPARATOR << (S ^ str) << SEPARATOR << endl;

   
   // --- unit cell mapping
   // some atoms of tStr are outside the cell now. Let's map them back
   cout << "shifted" << endl;
   cout << SEPARATOR << tStr << SEPARATOR << endl;
   cout << "shifted & mapped" << endl;
   tStr.map (SxCell::Positive);
   cout << SEPARATOR << tStr << SEPARATOR << endl;
   
   
   
   
   // === STRUCTURE COMPARISON ============
   
   // is shifted structure equivalent to original one?
   cout << "Structure shifted by " << t  << " is ";
   if (tStr != str) cout << "not ";
   cout << "equivalent to original one." << endl;

   t = Coord (5., 0., 0.);
   
   cout << "Structure shifted by " << t  << " is ";
   if ((str + t) != str) cout << "not ";
   cout << "equivalent to original one." << endl;





    
   // === MORE TRANSFORMATIONS =====================
   
   // SxAtomicStructure may be used for non-coordinate-type stuff,
   // such as displacements, forces, normal modes, ...
   
   SxAtomicStructure displace;
   // alternative setup: resize and set each atom by hand. 
   displace.resize (str.getNAtoms ());
   for (int i = 0; i < str.getNAtoms (); ++i)
      displace.ref(i) = Coord (i, 0, 0);

   // --- translations part II: atom specific displacements
   tStr = str + displace; // shift like this ...
   tStr -= displace; // or like this

   // for these uses, scaling may be useful
   displace /= 10;
   tStr = str + (2.5 * displace);
   
   // last but not least: species-specific scaling (divide or multiply)

   // copy species info from str to displace 
   displace.atomInfo = str.atomInfo; 

   SxVector<Double> specScale(str.getNSpecies ());
   specScale(0) = 0.2; specScale(1) = 0.3;
   
   tStr += displace / specScale; // divide   (also as /=)
   displace *= specScale;        // multiply (also as *)

}


