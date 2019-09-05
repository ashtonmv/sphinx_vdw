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

#include <SxAims.h>
#include <SxConstants.h>
#include <SxString.h>
#include <SxSort.h>

SxAims::SxAims ()
{
   // empty
}

SxAims::SxAims (const SxString &file_)
{
   setFilename (file_);
}

SxAims::~SxAims ()
{
   // empty
}

// Specifies the XYZ filename.
void SxAims::setFilename (const SxString &file_)
{
   filename = file_;
}

// Returns the atomic structure.
SxAtomicStructure &SxAims::getStructure ()
{
   return structure;
}

// Returns unique species.
SxList<SxString> &SxAims::getUniqueSpecies ()
{
   return uniqueSpecies;
}

// Read specified frames from a XYZ file.
void SxAims::read ()
{
   SX_CHECK (filename.getSize() > 0);

   SxArray<SxVector3<Double> > coords; // atomic coordinates (in Bohr)
   SxArray<SxString>           chemSymbols; // names of the elements
   SxMatrix3<Double>           cell; // Unit cell
   SxArray<ssize_t>            indices; // sorting order
   SxSort<SxString>            sort; // sorting class

   static const int BUFLEN = 10240;
   char buffer[BUFLEN]; // one line from the input file
   int nAtoms = 0;     // the number of atoms in the current frame
   int iLine=0;  // line counter
   double x,y,z; // the atomic coordinates

   int iSpecies=-1;
   SxString lastName;

   // --- open the file
   FILE *fp = fopen (filename.ascii (), "r");
   if (!fp)  {
      sxprintf ("Can't open Aims file '%s'\n", filename.ascii ());
      SX_EXIT;
   }

   int iLattice = 0; // Lattice Vector Counter
   while (!feof (fp))  {
      // --- Start reading and look for keywords
      if (!fgets (buffer, BUFLEN, fp)) 
         break; 
      iLine++;
      SxList<SxString> tokens = 
         SxString(buffer).left('\n').simplifyWhiteSpace ().tokenize(' ');
      if (tokens(0) == "lattice_vector") {
         if (iLattice > 2) {
            sxprintf ("More than 3 lattice vectors!\n");
            sxprintf ("Check Aims file '%s'\n", filename.ascii ());
            SX_QUIT;
         }
         x = tokens(1).toFloat ();
         y = tokens(2).toFloat ();
         z = tokens(3).toFloat ();
         cell(iLattice,0) = x * A2B;
         cell(iLattice,1) = y * A2B;
         cell(iLattice,2) = z * A2B;
         iLattice ++;
      }
      if (tokens(0) == "atom") {
         nAtoms++;
         x = tokens(1).toFloat ();
         y = tokens(2).toFloat ();
         z = tokens(3).toFloat ();
         coords.append(SxVector3<Double> (x, y, z) * A2B);
         chemSymbols.append(tokens(4));
      }
      if (tokens(0) == "atom_frac") {
         if (iLattice != 3) {
            sxprintf ("Cell is not specified!\n");
            sxprintf ("Check Aims file '%s'\n", filename.ascii ());
            SX_QUIT;
         }
         nAtoms++;
         x = tokens(1).toFloat ();
         y = tokens(2).toFloat ();
         z = tokens(3).toFloat ();
         coords.append(SxVector3<Double> (x, y, z) ^ cell);
         chemSymbols.append(tokens(4));
      }
   }

   // --- create atomic structure and chemList
   if (nAtoms > 0)  {
      iSpecies=-1;
      lastName = SxString ();
      structure = SxAtomicStructure ();
      uniqueSpecies = SxList<SxString> ();

      if (nAtoms > 1)  {
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
      }

      // --- create atomic structure
      structure.cell = cell.transpose ();
      structure.startCreation ();
      for (int i = 0; i < nAtoms; i++)  {
         if (lastName != chemSymbols(i))  {
            // --- new uniqueSpeciesList
            iSpecies++;
            uniqueSpecies.append (chemSymbols(i));
            lastName = chemSymbols(i);
            structure.newSpecies ();
         }
         structure.addAtom (iSpecies, coords(i));
      }
      structure.endCreation ();
   }

   fclose (fp);
}


// Write the data to an Aims file.
void SxAims::write (const SxAtomicStructure &structure_,
                    const SxArray<SxString> &chemName)
{
   FILE *fp = fopen (filename.ascii(), "w");
   if (!fp)  {
      sxprintf ("Can't open Aims file '%s' for write access.\n",
            filename.ascii());
      SX_EXIT;
   }
   
   // --- write lattice vectors
   Coord latA = structure_.cell.basis(0) / A2B;
   Coord latB = structure_.cell.basis(1) / A2B;
   Coord latC = structure_.cell.basis(2) / A2B;
   fprintf (fp,"lattice_vector   ");
   fprintf (fp,"%15.6f %15.6f %15.6f\n", latA(0), latA(1), latA(2));
   fprintf (fp,"lattice_vector   ");
   fprintf (fp,"%15.6f %15.6f %15.6f\n", latB(0), latB(1), latB(2));
   fprintf (fp,"lattice_vector   ");
   fprintf (fp,"%15.6f %15.6f %15.6f\n", latC(0), latC(1), latC(2));

   // --- write all atoms
   int nAtoms = structure_.getNAtoms ();
   for (int iAtom = 0; iAtom < nAtoms; iAtom++)  {
      int iSpecies = structure_.getISpecies(iAtom); 
      Coord coord  = structure_(iAtom);
      coord /= A2B;
      fprintf (fp,"atom   ");
      fprintf (fp,"%15.6f %15.6f %15.6f", coord(0), coord(1), coord(2));
      fprintf (fp,"   ");
      fprintf (fp,"%-2s\n", chemName(iSpecies).ascii());
   }
}
