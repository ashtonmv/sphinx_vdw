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

#include <SxLocpot.h>
#include <SxConstants.h>
#include <SxString.h>
#include <SxSort.h>

SxLocpot::SxLocpot ()
{
   // empty
}

SxLocpot::SxLocpot (const SxString &file_)
{
   setFilename (file_);
}

SxLocpot::~SxLocpot ()
{
   // empty
}

// Specifies the Locpot filename.
void SxLocpot::setFilename (const SxString &file_)
{
   filename = file_;
}

SxCell SxLocpot::readVASPCell (FILE *fp, const SxString &file)
{
   SX_CHECK (fp);
   SxCell result;
   char buffer[10240];
   if (!fgets (buffer, 10240, fp))  {
      cout << "Unexpected end of file " << file << endl;
      SX_EXIT;
   }
   double scale;
   fscanf (fp, "%lf", &scale);
   for (int iBasis = 0; iBasis < 3; ++iBasis)
      for (int xyz = 0; xyz < 3; ++xyz)
         fscanf (fp, "%lf", &result(xyz,iBasis));
   result*= scale * A2B;
   result.setup ();
   if (feof(fp))  {
      cout << "Reading cell: Unexpected end of VASP file "
           << file << endl;
      SX_EXIT;
   }
   char c = (char)fgetc (fp); // \n
   if (c != '\n')  {
      cout << "Problem while reading cell from VASP file "
           << file << endl;
      cout << "Expected to arrive at end of line after reading cell vectors,"
           << endl << "but found '" << c << "'." << endl;
      cout << "The cell that was read now is " << endl << result << endl;
      cout << "Please check file format or report this error." << endl;
      SX_EXIT;
   }
   
   return result;
}

void SxLocpot::read ()
{

   SxArray<SxString> chemNames;

   try  {
      SxBinIO io(filename, SxBinIO::ASCII_READ_ONLY);
      cell = readVASPCell(io.fp, filename);
      char buffer[10240];
      // read first line
      fgets(buffer, 10240, io.fp);
      // --- optional vasp comment line for species information chem names
      SxList<SxString> nAt = 
         SxString(buffer).left('\n').simplifyWhiteSpace ().tokenize(' ');
      if (nAt.getSize () == 0)  {
         cout << "Problem while reading VASP file " << filename << endl;
         cout << "Syntax error: cannot read structure information" << endl;
         cout << "Please check file format or report this error." << endl;
         SX_EXIT;
      }
      bool names = false;
      nAt(0).toInt (&names); // try to convert; otherwise set names=true 
      // --- if chem names read next line, these are the atom numbers
      if (names) {
         SxList<SxString>::ConstIterator it;
         for (it = nAt.begin(); it != nAt.end(); ++it)
            chemNames.append(*it);
         fgets(buffer, 10240, io.fp);
         nAt = SxString(buffer).left('\n').simplifyWhiteSpace ().tokenize(' ');
      } else  {
         int number = 0;
         SxList<SxString>::ConstIterator it;
         for (it = nAt.begin(); it != nAt.end(); ++it, number++)
            chemNames.append("POTCAR-ELEMENT-" + SxString (number));
      }
      if (nAt.getSize () == 0)  {
         cout << "Problem while reading VASP file " << filename << endl;
         cout << "Syntax error: number of atoms line is empty" << endl;
         cout << "Please check file format or report this error." << endl;
         SX_EXIT;
      }

      // --- number of atoms
      SxList<int> nAtomList;
      SxVector<Int> nAtoms(chemNames.getSize ());
      for (SxList<SxString>::Iterator it = nAt.begin ();
            it != nAt.end (); ++it)
      {
         if ((*it).getSize () > 0)
            nAtomList << it->toInt ();
      }
      nAtoms = nAtomList;
      // one line (type of coordinates)
      fgets(buffer, 10240, io.fp);
      {
         char c;
         if (sscanf (buffer, " %c", &c) == 0)  {
            cout << "Problem while reading VASP file " << filename << endl;
            cout << "Syntax error: cannot read coordinate type" << endl;
            cout << "Please check file format or report this error." << endl;
            SX_EXIT;
         }
      }
      char letter = SxString(buffer).tokenize(' ')(0).toUpper()(0);
      bool direct = true;
      switch (letter)  {
      case 'D' : direct = true;
                 break;
      case 'C' : direct = false;
                 break;
      case 'K' : direct = false;
                 break;
      default  : // unknown
                 sxprintf ("ERROR: Syntax error in file %s\n", 
                         filename.ascii());
                 SX_EXIT;
      }

      structure.cell = cell;
      structure.startCreation ();
      Coord coord;
      int iAtom, iSpecies;
      for (iSpecies = 0; iSpecies < nAtoms.getSize(); ++iSpecies)  {
         structure.newSpecies ();
         for (iAtom = 0; iAtom < nAtoms(iSpecies); ++iAtom)  {
            // --- atom coordinates
            fgets(buffer, 10240, io.fp);
            int nRead = sscanf (buffer, "%lf %lf %lf", 
                  &coord(0), &coord(1), &coord(2));
            if (nRead != 3)  {
               sxprintf ("ERROR: Syntax error in file %s\n", 
                     filename.ascii());
               SX_EXIT;
            }

            if (direct) coord = structure.cell ^ coord; // Angstroem -> Bohr is already in the cell
            else coord = A2B * coord; // Angstroem -> Bohr
            structure.addAtom (iSpecies, coord);
            if (feof(io.fp))  {
               cout << "Reading atoms: Unexpected end of VASP file "
                  << filename << endl;
               SX_EXIT;
            }
         }
      }
      structure.endCreation ();

      // --- mesh size
      if (fscanf (io.fp, "%d%d%d", &mesh(0), &mesh(1), &mesh(2)) != 3)  {
         cout << "Error reading mesh from VASP file " << filename << endl;
         SX_EXIT;
      }
      cout << "VASP mesh: " << mesh << endl; 
      // --- read potential
      potential.resize (mesh.product ());
      SxVector3<Int> x;
      for (x(2) = 0; x(2) < mesh(2); x(2)++)
         for (x(1) = 0; x(1) < mesh(1); x(1)++)
            for (x(0) = 0; x(0) < mesh(0); x(0)++)
               fscanf (io.fp, "%lf",
                     &potential(mesh.getMeshIdx(x, SxMesh3D::Positive)));
      if (feof(io.fp))  {
         cout << "Reading potential: Unexpected end of VASP file "
            << filename << endl;
         SX_EXIT;
      }
      io.close ();
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
}
