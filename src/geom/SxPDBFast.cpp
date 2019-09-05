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

#include <SxPDBFast.h>
#include <SxConstants.h>
#include <SxString.h>
#include <SxSort.h>

//#define SX_DEBUG_PDB_FORMAT_ECHO

SxPDBFast::SxPDBFast ()
{
   // empty
}

SxPDBFast::SxPDBFast (const SxString &file_)
{
   setFilename (file_);
}

SxPDBFast::~SxPDBFast ()
{
   // empty
}


void SxPDBFast::setFilename (const SxString &file_)
{
   filename = file_;
}


SxAtomicStructure &SxPDBFast::getStructure ()
{
   return structure;
}


SxList<SxString> &SxPDBFast::getUniqueSpecies ()
{
   return uniqueSpecies;
}


void SxPDBFast::read ()
{
   SX_CHECK (filename.getSize() > 0);

   SxArray<ssize_t>           indices; // sorting order
   SxSort<SxString>           sort; // sorting class

   SxList<SxVector3<Double> > coords; // atomic coordinates (in Bohr)
   SxList<SxString>           chemSymbols; // names of the elements
   SxMatrix3<Double>          cell; // Unit cell

   static const int BUFLEN = 10240;
   char buffer[BUFLEN];
   int i;

   SxString line;
   SxString recordName;  //  1 - 6
// int serial;           //  7 - 11
   SxString name;        // 13 - 14  chem. symbol
                         // 15       remoteness indicator (alphabetic)
                         // 16       Branch designator (numeric)
// SxString altLoc;      // 17
// SxString resName;     // 18 - 20
// SxString chainId;     // 22
// int resSeq;           // 23 - 26
// SxString iCode;       // 27
   double xOrtho;        // 31 - 38
   double yOrtho;        // 39 - 46
   double zOrtho;        // 47 - 54
// double occupancy;     // 55 - 60
// double tempFactor;    // 61 - 66
// SxString seqID;       // 73 - 76
// SxString element;     // 77 - 78
// SxString charge;      // 79 - 80

   FILE *fp = fopen (filename.ascii(), "r");
   if (!fp)  {
      sxprintf ("Can't open pdb file '%s'\n", filename.ascii());
      SX_EXIT;
   }
   while ( !feof (fp) )  {
      fgets(buffer, BUFLEN, fp); // skip rest of line
#     ifndef WIN32      
         if ( feof (fp) )  break;
#     endif /* WIN32 */      
      line = SxString (buffer);
      line.resize (80, true);
      recordName = line.subString( 0,  5).trim();
      if (recordName == "ATOM")  {
//       serial     = line.subString( 6, 10).toInt();
         name       = line.subString(12, 13).trim().toUpper();
         xOrtho     = line.subString(30, 37).toDouble();
         yOrtho     = line.subString(38, 45).toDouble();
         zOrtho     = line.subString(46, 53).toDouble();
         // --- remove numbers from element name
         for (i=0; i < name.getSize(); i++) if (isdigit(name(i)))  
            name(i) = ' ';
         name = name.trim();
   
         coords.append (SxVector3<Double> (xOrtho, yOrtho, zOrtho));
         chemSymbols.append (name);
      } else if (recordName == "CRYST1")  {
         double a     = line.subString(6,14).toDouble();
         double b     = line.subString(15,23).toDouble();
         double c     = line.subString(24,32).toDouble();
         double alpha = line.subString(33,39).toDouble();
         double beta  = line.subString(40,46).toDouble();
         double gamma = line.subString(47,53).toDouble();
         if (   fabs(alpha-90.) < 1e-10 
             && fabs(beta-90.)  < 1e-10 
             && fabs(gamma-90.) < 1e-10)
         {
            cell = SxMatrix3<Double> (a,0.,0., 0.,b,0., 0.,0.,c);
         }  else {
            alpha /= RAD2DEG;
            beta  /= RAD2DEG;
            gamma /= RAD2DEG;
         // --- see PDB guide:
         //     http://autodep.ebi.ac.uk/prepare_coords.html
         double omega = a*b*c  
                      * sqrt((1. - cos(alpha)*cos(alpha)
                                 - cos(beta) *cos(beta)
                                 - cos(gamma)*cos(gamma)
                                 + 2.*(cos(alpha)*cos(beta)*cos(gamma))));
         cell = SxMatrix3<Double> (
            a,  b*cos(gamma), c*cos(beta),
            0., b*sin(gamma), c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma),
            0., 0.,           omega/(a*b*sin(gamma))
         );
         }
         cell *= 1.8897261; // Angstroem -> bohr
      }
//       } else {
//          sxprintf ("PDB Error: Record name %s not supported. Ignored.\n",
//                  recordName.ascii());
//       }
   }
   fclose (fp);


   SxArray<SxVector3<Double> > coordsArray (coords);
   SxArray<SxString>           chemArray (chemSymbols);
   // --- sort chemSymbols
   //   input : H H O H O O Cl
   //           0 1 2 3 4 5 6
   // indices : 6 0 1 3 2 4 5
   indices = sort.quickSortToIdx (chemArray);
   // --- reorder chemSymbols & coords
   //     Cl H H H O O O
   //     to have continuous species
   chemArray.sortByIdx (indices);
   coordsArray.sortByIdx (indices);

   // --- fill SxAtomicStructure
   SxList<SxString> uniqueSpeciesList; // {"Cl", "H", "O"}
   ssize_t nElements = chemArray.getSize ();
   int iSpecies=-1;
   SxString lastName;

   structure     = SxAtomicStructure ();
   uniqueSpecies = SxList<SxString> ();

   if (nElements > 0)  {
      // --- generate uniqueSpeciesList
      structure.startCreation ();
      for (i=0; i < nElements; i++)  {
         if (lastName != chemArray(i))  {
            // --- new uniqueSpeciesList
            iSpecies++;
            uniqueSpeciesList.append (chemArray(i));
            lastName      = chemArray(i);
            structure.newSpecies ();
         }
         structure.addAtom (iSpecies, coordsArray(i));
      }
      structure.endCreation ();
      structure *= 1.8897261; // Angstroem -> bohr
      uniqueSpecies = uniqueSpeciesList;
      structure.atomInfo->meta.attach (SxAtomicStructure::Elements,
                                       uniqueSpecies);
   }
}

// Write the data to a PDB file.
void SxPDBFast::write (const SxAtomicStructure &structure_,
                       const SxArray<SxString> &chemName)
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

   Coord             coord; // atomic coordinate (in Bohr)
   SxMatrix3<Double> cell; // Unit cell

   /// \todo why not structure_.getCell ();
   cell = structure_.cell;

   SxString line;
   SxString recordName;  //  1 - 6
   int serial;           //  7 - 11
   SxString name;        // 13 - 14    chem. symbol
                         // 15         remoteness indicator (alphabetic)
                         // 16         branch designator (numeric)
   SxString altLoc;      // 17
   SxString resName;     // 18 - 20
   SxString chainId;     // 22
// int resSeq;           // 23 - 16
   SxString iCode;       // 27
   SxString xOrtho;      // 31 - 38
   SxString yOrtho;      // 39 - 46
   SxString zOrtho;      // 47 - 54
// double occupancy;     // 55 - 60
// double tempFactor;    // 61 - 66
   SxString seqID;       // 73 - 76
   SxString element;     // 77 - 78
   SxString charge;      // 79 - 80

   // --- open file for writing
   FILE *fp = fopen (filename.ascii(), "w");
   if (!fp)  {
      sxprintf ("Can't open pdb file '%s' for write access.\n",filename.ascii());
      SX_EXIT;
   }

   // --- setup cell
   SxVector3<Double> a1 = cell.col(0);
   SxVector3<Double> a2 = cell.col(1);
   SxVector3<Double> a3 = cell.col(2);
   double a     = sqrt(a1^a1) / 1.8897261;  // Bohr -> Angstroem
   double b     = sqrt(a2^a2) / 1.8897261;
   double c     = sqrt(a3^a3) / 1.8897261;
   double alpha = RAD2DEG * getAngle (a3, a2);
   double beta  = RAD2DEG * getAngle (a1, a3);
   double gamma = RAD2DEG * getAngle (a1, a2);
   SxString spaceGroup (" ");
   SxString zValue (" ");

   // --- write cell
   fprintf (fp, "%-6s", "CRYST1");                                //  1 -  6
   fprintf (fp, "%9s", SxString(a).ascii());                      //  7 - 15
   fprintf (fp, "%9s", SxString(b).ascii());                      // 16 - 24
   fprintf (fp, "%9s", SxString(c).ascii());                      // 25 - 33
   fprintf (fp, "%7s", SxString(alpha).ascii());                  // 34 - 40
   fprintf (fp, "%7s", SxString(beta).ascii());                   // 41 - 47
   fprintf (fp, "%7s", SxString(gamma).ascii());                  // 48 - 54
   fprintf (fp, " ");                                             // 55     
   fprintf (fp, "%11s", spaceGroup.ascii());                      // 56 - 66
   fprintf (fp, "%4s\n",  zValue.ascii());                        // 67 - 70

   // --- write all species
   serial = 1;
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
         xOrtho = SxString (coord(0), "% 7.3f");
         yOrtho = SxString (coord(1), "% 7.3f");
         zOrtho = SxString (coord(2), "% 7.3f");
         fprintf (fp,"%-6s", "ATOM");                             //  1 -  6
         fprintf (fp,"%5d", serial);                              //  7 - 11
         fprintf (fp," ");                                        // 12
         fprintf (fp,"%-2s", name.ascii());                       // 13 - 14
         fprintf (fp," ");                                        // 15
         fprintf (fp," ");                                        // 16
         fprintf (fp," ");                                        // 17
         fprintf (fp,"   ");                                      // 18 - 20
         fprintf (fp," ");                                        // 21
         fprintf (fp," ");                                        // 22
         fprintf (fp,"    ");                                     // 23 - 26
         fprintf (fp," ");                                        // 27
         fprintf (fp,"   ");                                      // 28 - 30
         fprintf (fp,"%8s", xOrtho.ascii());                      // 31 - 38
         fprintf (fp,"%8s", yOrtho.ascii());                      // 39 - 46
         fprintf (fp,"%8s", zOrtho.ascii());                      // 47 - 54
         fprintf (fp,"\n");
         serial++;
      }
   }

   fclose (fp);
}

