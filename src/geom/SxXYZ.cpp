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

#include <SxXYZ.h>
#include <SxConstants.h>
#include <SxString.h>
#include <SxSort.h>

//#define SX_DEBUG_XYZ_FORMAT_ECHO

SxXYZ::SxXYZ ()
{
   // empty
}

SxXYZ::SxXYZ (const SxString &file_)
{
   setFilename (file_);
}

SxXYZ::~SxXYZ ()
{
   // empty
}

// Specifies the XYZ filename.
void SxXYZ::setFilename (const SxString &file_)
{
   filename = file_;
}

// Returns the number of loaded frames.
int SxXYZ::getNFrames ()
{
   return int(structures.getSize ());
}

// Returns the atomic structure.
SxAtomicStructure &SxXYZ::getStructure (int idx_)
{
   return structures(idx_);
}

// Returns unique species.
SxList<SxString> &SxXYZ::getUniqueSpecies (int idx_)
{
   return uniqueSpecies(idx_);
}

// Read all frames from a XYZ file
void SxXYZ::read ()
{
   read (0, -1, 1);
}

// Read specified frames from a XYZ file.
void SxXYZ::read (int startFrame_, int endFrame_, int step_)
{
   SX_CHECK (filename.getSize() > 0);

   SxArray<SxVector3<Double> > coords; // atomic coordinates (in Bohr)
   SxArray<SxString>           chemSymbols; // names of the elements
   SxMatrix3<Double>           cell; // Unit cell
   SxArray<ssize_t>            indices; // sorting order
   SxSort<SxString>            sort; // sorting class

   static const int BUFLEN = 10240;
   char buffer[BUFLEN]; // one line from the input file
   char name[BUFLEN]; // chemicalName or the atomic number
   int nRead;    // sscanf(): The number of items succesfully read.
   int nAtoms=0; // the number of atoms in the current frame
   int iAtom=0;  // counter [0 .. nAtoms)
   int iFrame=0; // frame counter
   int iLine=0;  // line counter
   double x,y,z; // the atomic coordinates

   SxList<SxString> uniqueSpeciesList; // {"Cl", "H", "O"}
   SxAtomicStructure structure;
   int nElements;
   int i;
   int iSpecies=-1;
   SxString lastName;
   SxList<SxAtomicStructure> tmpStructures;
   SxList<SxList<SxString> > tmpSpecies;

   // --- open the file
   FILE *fp = fopen (filename.ascii (), "r");
   if (!fp)  {
      sxprintf ("Can't open xyz file '%s'\n", filename.ascii ());
      SX_EXIT;
   }

   while (!feof (fp) && (iFrame <= endFrame_ || endFrame_ < 0))  {
      // --- read the number of atoms in the current frame
      iLine++;
      if (!fgets (buffer, BUFLEN, fp))  {
         break;
      }
      nAtoms = 0;
      nRead = sscanf (buffer, "%d", &nAtoms);
      if (nRead < 1 || nAtoms < 1)  {
         cout << "\"" << filename.ascii () << "\" Line " << iLine
              << ": The number of atoms expected." << endl;
         SX_EXIT;
      }

#  ifdef SX_DEBUG_XYZ_FORMAT_ECHO
      cout << "frame#" << iFrame
           << " nAtoms=" << nAtoms << endl;
#  endif /* SX_DEBUG_XYZ_FORMAT_ECHO */

      // --- read the title of the molecule (can be blank)
      iLine++;
      if (!fgets (buffer, BUFLEN, fp))  {
         cout << "\"" << filename.ascii () << "\" Line " << iLine
              << ": The title of the molecule or a blank line expected."
              << endl;
         SX_EXIT;
      }

      // --- atoms
      iAtom = 0;
      if (iFrame < startFrame_)  {
         // --- skip this frame: just read the lines
         // --- this approach reads a single frame and validates
         // --- the whole input file (testing)
         while (iAtom < nAtoms)  {
            iLine++;
            if (!fgets (buffer, BUFLEN, fp))  {
               cout << "\"" << filename.ascii () << "\" Line " << iLine
                  << ": 'chemName x y z' format expected. " << endl;
               SX_EXIT;
            }
            iAtom++;
         }
      }  else  {
         // --- read the atoms from this frame
         chemSymbols.resize (nAtoms);
         coords.resize (nAtoms);

         while (iAtom < nAtoms)  {
            iLine++;
            if (!fgets (buffer, BUFLEN, fp))  {
               cout << "\"" << filename.ascii () << "\" Line " << iLine
                  << ": 'chemName x y z' format expected. " << endl;
               SX_EXIT;
            }
            nRead = sscanf (buffer, "%s %lf %lf %lf", name, &x, &y, &z);
            if (nRead < 4)  {
               cout << "\"" << filename.ascii () << "\" Line " << iLine
                  << ": 'chemName x y z' format expected." << endl;
               SX_EXIT;
            }

            // --- add a new atom
            chemSymbols(iAtom) = SxString (name);
            coords(iAtom)      = SxVector3<Double> (x, y, z)*1.8897261; // Bohr

#        ifdef SX_DEBUG_XYZ_FORMAT_ECHO
            cout << "Atom " << iAtom+1 << "/" << nAtoms
               << " \"" << name << "\" ["
               << x << ", " << y << ", " << z << "]"
               << endl;
#        endif /* SX_DEBUG_XYZ_FORMAT_ECHO */
            iAtom++;
         } // end read the atoms from this frame

         // --- create frame
         nElements = int(chemSymbols.getSize ());
         if (nElements > 0)  {
            iSpecies=-1;
            lastName = SxString ();
            structure = SxAtomicStructure ();
            uniqueSpeciesList = SxList<SxString> ();

            // --- sort chemSymbols
            //   input : H H O H O O Cl
            //           0 1 2 3 4 5 6
            // indices : 6 0 1 3 2 4 5
            indices = sort.quickSortToIdx (chemSymbols);
            // --- reorder chemSymbols & coords
            //     Cl H H H O O O
            //     to have continuous species
            chemSymbols.sortByIdx (indices);
            coords.sortByIdx (indices);

            // --- create atomic structure
            structure.startCreation ();
            for (i=0; i < nElements; i++)  {
               if (lastName != chemSymbols(i))  {
                  // --- new uniqueSpeciesList
                  iSpecies++;
                  uniqueSpeciesList.append (chemSymbols(i));
                  lastName = chemSymbols(i);
                  structure.newSpecies ();
               }
               structure.addAtom (iSpecies, coords(i));
            }
            structure.endCreation ();

            tmpStructures.append (structure);
            tmpSpecies.append (uniqueSpeciesList);
         }
      }
      iFrame+=step_;
   } // end: read the frames

   // --- copy
   structures = tmpStructures;
   uniqueSpecies = tmpSpecies;

   fclose (fp);
}


// Write the data to a PDB file.
void SxXYZ::write (const SxAtomicStructure &structure_,
                   const SxArray<SxString> &chemName,
                   bool writeAMat)
{
   SX_CHECK      (filename.getSize() > 0);

   int iAtom, nAtoms;
   int iSpecies, nSpecies = structure_.getNSpecies();

   if (nSpecies == 0)  {
      // --- abortion, the structure is empty
      return;
   }

   if (nSpecies != chemName.getSize ())  {
      sxprintf ("pdb file '%s': Some species info is missing.\n",
              filename.ascii());
      SX_EXIT;
      // --- abortion, some species info is missing
      //return;
   }


   // --- open file for writing
   FILE *fp = fopen (filename.ascii(), "w");
   if (!fp)  {
      sxprintf ("Can't open xyz file '%s' for write access.\n",filename.ascii());
      SX_EXIT;
   }

   // --- write cell
   if (writeAMat)  {
      SxMatrix3<Double> cell = structure_.cell;
      SxVector3<Double> a1 = cell.col(0);
      SxVector3<Double> a2 = cell.col(1);
      SxVector3<Double> a3 = cell.col(2);
      fprintf (fp, "%15.12f %15.12f %15.12f\n", a1(0), a1(1), a1(2));
      fprintf (fp, "%15.12f %15.12f %15.12f\n", a2(0), a2(1), a2(2));
      fprintf (fp, "%15.12f %15.12f %15.12f\n", a3(0), a3(1), a3(2));
   }

   // --- write all species
   SxString name;
   Coord  coord;
	
	fprintf (fp,"%i\n\n",structure_.nTlAtoms);
	
   for (iSpecies=0; iSpecies < nSpecies; ++iSpecies)  {
      nAtoms = structure_.getNAtoms (iSpecies);
      name   = chemName(iSpecies).trim().toUpper();
      name.resize (2, true);
      name   = name.subString(0, minimum(name.getSize(), (ssize_t)2));
      //if (name != "H ")  name = name.adjustRight ();
      name = name.adjustRight ();

      // --- write all atoms of the curent species
      for (iAtom=0; iAtom < nAtoms; ++iAtom)  { 
         coord  = structure_(iSpecies, iAtom);
         coord /= 1.8897261;  // Bohr -> Angstroem
         fprintf (fp,"%-2s", name.ascii());
         fprintf (fp,"   ");
         fprintf (fp,"%15.12f %15.12f %15.12f\n", coord(0), coord(1), coord(2));
      }
   }

   fclose (fp);
}
