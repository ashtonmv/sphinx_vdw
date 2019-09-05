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

#include <SxCLI.h>
#include <SxConfig.h>
#include <SxIO.h>
#include <SxParser.h>
#include <SxAtomicStructure.h>
#include <SxGrid.h>

int main (int argc, char **argv)
{
   // --- init SFHIngX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.preUsageMessage = "This add-on generates a set of unique displacements.";
   cli.authors = "C. Freysoldt";

   SxString inFile 
      = cli.option ("-i|--input", "input file", "supercell input file")
        .toString ("input.sx");
   
   SxString outFile
      = cli.option ("-o","filename", "output file name base")
        .toString ("input-disp-");

   double displ
      = cli.option ("-d","distance", "displacement magnitude").toDouble (0.1);

   bool printDispl
      = cli.option ("--show", "show the displacements").toBool ();

   SxMatrix3<Double> dispDirs(1., 0., 0., 0., 1., 0., 0., 0., 1.);
   if (cli.option ("--dirs","matrix", "9 numbers representing 3 "
                   " displacement directions").exists ())
      dispDirs = SxMatrix3<Double> (cli.last ().toDoubleList ()).transpose ();

   bool header = cli.option ("+h|--header", "print structure file header")
                 .toBool ();

   bool debug = cli.option ("--verbose","print some debugging info").toBool ();

   cli.finalize ();

   if (fabs(dispDirs.determinant ()) < 1e-7)  {
      cout << "Displacement directions are linear dependent!" << endl;
      SX_QUIT;
   }
   // --- normalize column vectors
   for (int i = 0; i < 3; ++i)  {
      dispDirs.setCol (i, dispDirs.col(i).normalize ());
      cout << "Displacement " << (i+1) << " along " << dispDirs.col(i) << endl;
   }
   
   // --- read input
   SxParser parser;
   SxConstPtr<SxSymbolTable> tree;
   tree = parser.read (inFile, "std/structure.std");
   SxAtomicStructure structure = SxAtomicStructure (&*tree);
   structure.epsEqual = 2e-3;

   // symmetry variables
   structure.updateSymGroup (SxAtomicStructure::SupercellCompatible);
   SxSymGroup &syms = *structure.cell.symGroupPtr;
   if (debug) syms.print ();
   int nSym = syms.getSize ();

   // primitive cell output
   cout << "PrimitiveCell: " << structure.getPrimitiveCell () << endl;

   // --- get mapping after applying symmetries
   SxGrid grid(structure, 3);
   SxArray<SxConstPtr<SxAtomInfo> > aI(nSym);
   for (int iSym = 0; iSym < nSym; ++iSym)  {
      aI(iSym) = structure.match (grid, syms(iSym) ^ structure);
   }
   
   // --- generate unique displacements
   int nAtoms = structure.getNAtoms ();
   SxAtomicStructure displacement;
   SxList<SxAtomicStructure> uniqDisps;
   SxVector<TPrecTauR> ortho;
   SxList<SxVector<TPrecTauR> >::Iterator dIt;
   SxAtomicStructure rotD;
   bool append;
   // basis vectors of space spanned by displacements + its symmetry equivalents
   SxList<SxVector<TPrecTauR> > dispSpace; 
   /* it is enough to know the inequivalent atoms!
      inequivalent dir's are determined as follows:
   SxSymType symType = SxRotation(rotMatrix).getType ();

   SxVector3<Double> axis = symType.opCoord;

   The classification (Mirror, Rotation, ..., see
   SxSymType::Classification) is also available, e.g.

   if (symType == SxSymType::Mirror)  {
   ...
   }  */

   /* --- THIS PART TRIES TO DETERMINE an optimal dispDirs
   // it is not used, since user-supplied ones are better for symmetry
   int found = 0;
   SxCell rot1, rot2; 
   if (nSym == 1) cout << "3 directions necessary" << endl;
   else  {
      for (int iSym = 0; iSym < nSym; iSym++)  {
         if (fabs(syms(iSym).rot.trace ()) < 3)//not 1 or -1.
            if (found == 0)  {
               found++;
               rot1 = syms(iSym).rot;
            }  else  {
               if (syms(iSym).rot != -1 * rot1)  {//not inverse of previous
                  found++;
                  rot2 = syms(iSym).rot;
                  break;
               }
            }
      }
   }
   int iBas;
   SxMatrix3<Double> disp;
   Coord basVec;
   if (found == 2)  {//only one dir necessary
      for (iBas = 0; iBas < 3; iBas++)  {//first try the basis vectors
         basVec = structure.cell.basis (iBas);
         disp = SxMatrix3<Double> (basVec, rot1 ^ basVec, rot2 ^ basVec);
         if (fabs(disp.determinant ()) > 1e-4) break;
      }
      if (iBas == 3)  {//basVec not usefull, try their sum
         basVec = structure.cell.basis (0) 
                + structure.cell.basis (1) + structure.cell.basis (2);
         disp = SxMatrix3<Double> (basVec, rot1 ^ basVec, rot2 ^ basVec);
         if (fabs(disp.determinant ()) < 1e-4)  {//then use a nasty vector
            basVec = Coord (SxList<double> () << 1 << ::exp (1) << PI);
            disp = SxMatrix3<Double> (basVec, rot1 ^ basVec, rot2 ^ basVec);
            if (fabs(disp.determinant ()) < 1e-4)  {
               cout << "Error in SxUniqDisp " << endl
                    << "could not find single uniq dir!" << endl;
               SX_QUIT;
            }
         }
      }
      cout << "Usefull unique dir for disp is " << basVec << endl;
   }  else if (found == 1)  {//two dirs necessary
      for (iBas = 0; iBas < 3; iBas++)  {//first try the basis vectors
         basVec = structure.cell.basis (iBas);
         basVec.normalize ();
         if (fabs (basVec ^ (rot1 ^ basVec)) < 1)  {break; }
      }
      if (iBas == 3)  {//basVec not usefull, try their sum
         basVec = structure.cell.basis (0) 
                + structure.cell.basis (1) + structure.cell.basis (2);
         basVec.normalize ();
         if (fabs (basVec ^ (rot1 ^ basVec)) ==  1)  {//then use a nasty vector
            basVec = Coord (SxList<double> () << 1 << ::exp (1) << PI);
            basVec.normalize ();
            if (fabs (basVec ^ (rot1 ^ basVec)) ==  1)  {
               cout << "Error in SxUniqDisp " << endl
                    << "could not find double uniq dir!" << endl;
               SX_QUIT;
            }
         }
      }
      disp = SxMatrix3<Double> (basVec, rot1 ^ basVec, 
                   Coord (SxList<double> () << 1 << PI << ::exp (1)));
      cout << "Usefull unique dirs for disp are " << basVec << endl << " and "
           << SxCell (disp).getReciprocalCell ().basis (2).normalize () << endl;
      //the third vector is normalized and orthogonal to the first 2.
   }  else  {//found == 0
      cout << "3 (arbitrary) directions necessary" << endl; 
   }//TODO make this much shorter!
   */

   for (int iDof = 0; iDof < 3 * nAtoms; ++iDof)  {
      // --- set up new displacement
      displacement = structure.getNewStr ();
      displacement.set(0., 0., 0.);
      displacement.ref(iDof / 3) = dispDirs.col(iDof % 3);

      // --- check for symmetry reduction
      append = true;
      ortho.copy (displacement.coords);
      for (dIt = dispSpace.begin (); dIt != dispSpace.end (); ++dIt)  {
         // Gram-Schmidt orthogonalization to existing displacement space
         // note: dot ignores shape, i.e. refers to dof representation
         ortho.plus_assign_ax(-dot(ortho,*dIt), *dIt);
         if (ortho.normSqr () < 1e-7)  { append = false; break; }
      }
      // and append the unique displacements
      if (append)  {
         if (printDispl)  {
            int ia = iDof / 3;
            cout << "Displace atom " << (ia + 1) 
                 << " along " << dispDirs.col(iDof % 3) 
                 << " => " << (structure(ia) + displ * displacement(ia))
                 << endl;
         }
         if (debug)
            cout << "Accepted";
         else 
            cout << '+';
         
         // --- add to unique displacement list
         uniqDisps << displacement;

         // --- update subspace covered by uniqDisps + its symmetry-equivalents
         for (int iSym = 0; iSym < nSym && append; ++iSym)  {
            // get rotated and reordered displacement
            rotD = SxAtomicStructure (aI(iSym), 
                                      syms.getRot(iSym) ^ displacement);
            // orthogonalize to current subspace
            for (dIt = dispSpace.begin (); dIt != dispSpace.end (); ++dIt)
               rotD.coords.plus_assign_ax(-dot(rotD.coords,*dIt), *dIt);
            double normSqr;
            if ((normSqr = rotD.coords.normSqr ()) > 1e-7)  {
               // new basis vector
               dispSpace << (rotD.coords / sqrt(normSqr));
            }
         }
      } else {
         if (debug) 
            cout << "Rejected";
         else 
            cout << '.';
      }
      if (debug) 
         cout << " displacement no. " << (iDof+1)
              << "; remaining norm: " << ortho.normSqr () << endl;
      cout.flush ();
      // --- stop if dispSpace spans all degrees of freedom
      if (dispSpace.getSize () == 3 * nAtoms)  {
         if (debug)
            cout << "Degree-of-freedom space complete - no further"
                    " displacements considered.";
         break;
      }
   }
   cout << endl;
   
   // --- output
   FILE *file;
   SxList<SxAtomicStructure>::Iterator dispIt = uniqDisps.begin ();
   SxString name;
   for (int i = 0; i < uniqDisps.getSize (); ++i, ++dispIt)  {
      name = outFile + (i+1) + ".sx";
      if ((file = fopen(name.ascii(),"w")) == NULL)  {
         cout << "Can't open '" << name << "'." << endl;
         SX_EXIT;
      }
      if (header)
         fprintf(file, "format structure;\n");
      // print displaced structure
      (structure + displ * (*dispIt)).fprint (file);
      fclose (file);
   }

   if (uniqDisps.getSize () > 1)  {  
      cout << SX_SEPARATOR <<
   "| Warning: the displacements are unique in the sense that no displacement\n"
   "| can be generated from the others. However, it may not be the minimal\n"
   "| set with this property. " << SX_SEPARATOR;
   }
        
}
