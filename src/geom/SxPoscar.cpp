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

#include <SxPoscar.h>
#include <SxConstants.h>
#include <SxString.h>
#include <SxSort.h>

SxPoscar::SxPoscar ()
{
   // empty
}

SxPoscar::SxPoscar (const SxString &file_)
{
   setFilename (file_);
}

SxPoscar::~SxPoscar ()
{
   // empty
}

// Specifies the Poscar filename.
void SxPoscar::setFilename (const SxString &file_)
{
   filename = file_;
}

// Returns the atomic structure.
SxAtomicStructure SxPoscar::getStructure () const
{
   return structure;
}

// Return chemNames
SxArray<SxString> SxPoscar::getUniqueSpecies() const
{
   return chemNames;
}

// Read from a Poscar file.
void SxPoscar::read ()
{
   static const int BUFLEN = 10240;
   char buffer[BUFLEN];
   int iLine = 0, nRead = 0;

   FILE *fp = fopen (filename.ascii(), "r");
   if (!fp)  {
      sxprintf ("ERROR: Can't open file '%s'\n", filename.ascii());
      SX_EXIT;
   }

   // --- skip comment (line 1 of POSCAR)
   fgets (buffer, BUFLEN, fp); iLine++;

   // --- read scaling factor (line 2 of POSCAR)
   double scaling = 0;
   fgets (buffer, BUFLEN, fp); iLine++;
   nRead = sscanf (buffer, "%lf", &scaling);
   if (nRead != 1)  {
      sxprintf ("ERROR: Syntax error in line %d of file %s\n", 
              iLine, filename.ascii());
      SX_EXIT;
   }

   // --- read lattice vectors (lines 3-5 of POSCAR)
   Coord a1, a2, a3;
   fgets (buffer, BUFLEN, fp); iLine++;
   nRead = sscanf (buffer, "%lf %lf %lf", &a1(0), &a1(1), &a1(2));
   if (nRead != 3)  {
      sxprintf ("ERROR: Syntax error in line %d of file %s\n", 
              iLine, filename.ascii());
      SX_EXIT;
   }
   fgets (buffer, BUFLEN, fp); iLine++;
   nRead = sscanf (buffer, "%lf %lf %lf", &a2(0), &a2(1), &a2(2));
   if (nRead != 3)  {
      sxprintf ("ERROR: Syntax error in line %d of file %s\n", 
              iLine, filename.ascii());
      SX_EXIT;
   }
   fgets (buffer, BUFLEN, fp); iLine++;
   nRead = sscanf (buffer, "%lf %lf %lf", &a3(0), &a3(1), &a3(2));
   if (nRead != 3)  {
      sxprintf ("ERROR: Syntax error in line %d of file %s\n", 
              iLine, filename.ascii());
      SX_EXIT;
   }
  
    // if scaling (line 2 of POSCAR) is defined negative, scaling is the volume of the supercell 
   if (scaling < 0.0) {
    double volume = -scaling;
    SxMatrix3<Double> cell(a1, a2, a3);
    double vCell = cell.determinant ();
    scaling =  pow(volume/vCell,(1.0/3.0));
    }

    // SxMatrix3 has to be transposed ( = convention to read in matrix ) otherwise trouble with e.g. hcp cells
   structure.cell = scaling * (SxMatrix3<Double> (a1, a2, a3)).transpose () * A2B; // Angstroem -> Bohr

   // --- read species affiliation
   fgets (buffer, BUFLEN, fp); iLine++;
   SxList<SxString> tokens = SxString(buffer).left('\n').simplifyWhiteSpace ().tokenize(' ');
   SxArray<int> nAtoms;

   // chemNames or atomnumbers
   bool names = false;
   tokens(0).toInt (&names); // try to convert; set names=true if not possible 
   if (names) {
      SxList<SxString>::ConstIterator it;
      for (it = tokens.begin(); it != tokens.end(); ++it)
         chemNames.append(*it);
      fgets (buffer, BUFLEN, fp); iLine++;
      tokens = SxString(buffer).left('\n').tokenize(' ');
   } else  {
      int number = 0;
      SxList<SxString>::ConstIterator it;
      for (it = tokens.begin(); it != tokens.end(); ++it, number++)
         chemNames.append("POTCAR-ELEMENT-" + SxString (number));
   }

   try {
      SxList<int> nAtomList;
      SxList<SxString>::ConstIterator it;
      for (it = tokens.begin(); it != tokens.end(); ++it)
         nAtomList << it->toInt ();
      nAtoms = nAtomList;
   } catch (SxException e) {
      sxprintf ("ERROR: Syntax error in line %d of file %s\n", 
              iLine, filename.ascii());
      e.print ();
      SX_EXIT;
   }

   // --- get coordinate type (VASP reads only first character!)
   //FIXME: Selective dynamics comes on the line BEFORE the actual coord type!
   enum CoordMode { SelectiveDynamics, Direct, Cartesian } coordMode;
   fgets (buffer, BUFLEN, fp); iLine++;
   tokens = SxString(buffer).tokenize(' ');
   if (tokens.getSize() == 0)  {
      sxprintf ("ERROR: Syntax error in line %d of file %s\n", 
              iLine, filename.ascii());
      SX_EXIT;
   }
   char letter = tokens(0).toUpper()(0);
   switch (letter)  {
      case 'S' : coordMode = SelectiveDynamics;
                 break;
      case 'D' : coordMode = Direct;
                 break;
      case 'C' : coordMode = Cartesian;
                 break;
      case 'K' : coordMode = Cartesian;
                 break;
      default  : // unknown
                 sxprintf ("ERROR: Syntax error in line %d of file %s\n", 
                         iLine, filename.ascii());
                 SX_EXIT;
   }

   if (coordMode == SelectiveDynamics)  {
      fgets (buffer, BUFLEN, fp); iLine++;
      tokens = SxString(buffer).tokenize(' ');
      if (tokens.getSize() == 0)  {
         sxprintf ("ERROR: Syntax error in line %d of file %s\n", 
               iLine, filename.ascii());
         SX_EXIT;
      }
      letter = tokens(0).toUpper()(0);
      switch (letter)  {
         case 'D' : coordMode = Direct;
                    break;
         case 'C' : coordMode = Cartesian;
                    break;
         case 'K' : coordMode = Cartesian;
                    break;
         default  : // unknown
                    sxprintf ("ERROR: Syntax error in line %d of file %s\n", 
                          iLine, filename.ascii());
                    SX_EXIT;
      }
   }

   // --- read coordinates
   structure.startCreation ();
   Coord coord;
   int iAtom, iSpecies;
   for (iSpecies = 0; iSpecies < nAtoms.getSize(); ++iSpecies)  {
      structure.newSpecies();
      for (iAtom = 0; iAtom < nAtoms(iSpecies); ++iAtom)  {
         fgets (buffer, BUFLEN, fp); iLine++;
         nRead = sscanf (buffer, "%lf %lf %lf", 
                         &coord(0), &coord(1), &coord(2));
         if (nRead != 3)  {
            sxprintf ("ERROR: Syntax error in line %d of file %s\n", 
                    iLine, filename.ascii());
            SX_EXIT;
         }

         if (coordMode == Direct) coord = structure.cell ^ coord; // Angstroem -> Bohr is already in the cell
         else coord = A2B * coord * scaling; // Angstroem -> Bohr
         structure.addAtom (iSpecies, coord);
      }
   }
   structure.endCreation ();

   fclose (fp);

}


// Write the data to a POSCAR file.
void SxPoscar::write (const SxAtomicStructure &structure_,
                      const SxArray<SxString> &chemName)
{
   SX_CHECK (filename.getSize() > 0);

   int iAtom, nAtoms;
   int iSpecies, nSpecies = structure_.getNSpecies();

   if (nSpecies == 0)  {
      // --- abortion, the structure is empty
      return;
   }

   // --- open file for writing
   FILE *fp = fopen (filename.ascii(), "w");
   if (!fp)  {
      sxprintf ("Can't open POSCAR file '%s' for write access.\n",filename.ascii());
      SX_EXIT;
   }

   // --- write header
   fprintf (fp, " Generated by S/PHI/nX\n");

   // --- write cell
   SxMatrix3<Double> cell = structure_.cell;
   SxVector3<Double> a1 = cell.col(0) / A2B;  // Bohr -> Angstroem
   SxVector3<Double> a2 = cell.col(1) / A2B;
   SxVector3<Double> a3 = cell.col(2) / A2B;
   fprintf (fp, "1.0\n");                            // scaling factor
   fprintf (fp, "%15.10f %15.10f %15.10f\n", a1(0), a1(1), a1(2));  // a1
   fprintf (fp, "%15.10f %15.10f %15.10f\n", a2(0), a2(1), a2(2));  // a2
   fprintf (fp, "%15.10f %15.10f %15.10f\n", a3(0), a3(1), a3(2));  // a3

   // --- write species affiliation
   for (iSpecies=0; iSpecies < nSpecies; ++iSpecies)  {
      fprintf (fp, "%d ", structure_.getNAtoms(iSpecies));
   }
   fprintf (fp, "\n");

   // --- write all atomic coordinates
   fprintf (fp, "Cartesian\n");
   Coord  coord;
   for (iSpecies=0; iSpecies < nSpecies; ++iSpecies)  {
      nAtoms = structure_.getNAtoms(iSpecies);

      // --- write all atoms of the curent species
      for (iAtom=0; iAtom < nAtoms; ++iAtom)  { 
         coord  = structure_(iSpecies, iAtom);
         coord /= A2B;  // Bohr -> Angstroem
         fprintf (fp,"%15.10f %15.10f %15.10f\n", coord(0), coord(1), coord(2));
      }
   }
   
   // --- write species order
   SxString name;
   sxprintf ("\nSpecies order:\n");
   for (iSpecies=0; iSpecies < nSpecies; ++iSpecies)  {
      name = chemName(iSpecies).trim().toUpper();
      name.resize (2, true);
      name = name.subString(0, minimum(name.getSize(), (ssize_t)2));
      //if (name != "H ")  name = name.adjustRight ();
      name = name.adjustRight ();
      sxprintf ("   %-2s - %s\n", 
              name.ascii(), chemName(iSpecies).ascii());
   }

   fclose (fp);
}
