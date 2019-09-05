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
#include <SxBinIO.h>
#include <SxParser.h>
#include <SxAtomicStructure.h>
#include <SxNeighbors.h>
#include <SxEquivalence.h>

void fprint(FILE *file, const SxAtomicStructure &str, 
            const SxArray<SxString>& symbol, int ia)
{
   SX_CHECK (ia < str.getNAtoms (), ia, str.getNAtoms ());
   const char *label = str.hasLabels () ? str.getLabels ()(ia).ascii () : NULL;
   int is = 0;
   while (ia >= str.getNAtoms(is))  ia -= str.getNAtoms(is++);

   sxfprintf(file, "%s%i", symbol(is).ascii (), (ia+1));
   if (label) sxfprintf (file, "(%s)", label);
   
}

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.preUsageMessage =
      "This add-on prints structural parameters like distances and angles. "
      "Lists of atoms are comma-separated lists of indices or index-ranges "
      "a-b, e.g.\n\n"
      "1,4-9,12-13   corresponds to 1 4 5 6 7 8 9 12 13\n\n"
      "The index runs through a) all species and b) all atoms, i.e. for "
      "NH3 it would be\n"
      "1 N1\n2 H1\n3 H2\n4 H3\n";
      
   cli.authors = "C. Freysoldt";

   bool sxbFile = cli.option ("-b|--sxb", "input file is binary waves file "
                              "rather than S/PHI/nX input file")
                  .toBool ();
   SxString inFile 
      = cli.option ("-i|--input", "input file", 
                    "take original input file")
        .toString (sxbFile ? "waves.sxb" : "input.sx");
   cli.last ().defaultValue = "default: input.sx, or waves.sxb for --sxb flag";
   
   SxString outFile
      = cli.option ("-o","filename", "output file name (screen otherwise)")
        .toString ("");
   
   double maxDist  
      = cli.option ("-c|--cutoff","distance", "maximum neighbor distance")
        .toDouble (10,0);
   double minDist  
      = cli.option ("--cMin","distance", "minimum neighbor distance")
        .toDouble (0,0);

   bool angles = cli.option ("-a|--angles","printout angles")
                 .toBool ();

   SxList<int> atoms = cli.option("-n|--atoms","list",
                       "list of atom indices").toIdxList ();
   cli.last ().defaultValue = "default: all";

   SxList<int> dihedral = cli.option("-d|--dihedral","list",
      "list of atom indices for which to print dihedral angles").toIdxList ();

   bool printReduced
      = cli.option ("--printReduced", "file",
                    "print symmetry-reduced list of neighbors").toBool ();
   FILE *reducedFile = stdout;
   if (cli.last ().hasValue ())
      reducedFile = fopen (cli.last ().getValue ().ascii (), "w");
   
   SxString bondList = cli.option ("--bondlist", "file",
                                   "print list of bonded atoms").toString ("");

   int centerGroup = cli.newGroup ("reference center");
   Coord center(cli.option ("--center", "center",
                            "find neighbors of arbitrary point")
                .toList3 ());

   cli.finalize ();
   bool refCenter = cli.groupAvailable (centerGroup);

   if (!reducedFile)  {
      cout << "Cannot open file for reduced neighbor printout" << endl;
      SX_QUIT;
   }

   // --- read input
   SxAtomicStructure structure;
   SxArray<SxString> chemNames;
   if (sxbFile)  {
      try {
         SxBinIO io (inFile, SxBinIO::BINARY_READ_ONLY);
         structure.read (io);
         SxString chemNameList;
         io.read("chemNames", &chemNameList);
         chemNames = chemNameList.tokenize (',');
         io.close ();
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
   } else  {
      SxParser parser;
      SxConstPtr<SxSymbolTable> tree;
      tree = parser.read (inFile, "std/structure.std");
      structure = SxAtomicStructure (&*tree);
      chemNames = SxSpeciesData::getElements (&*tree);
   }

   
   // --- output
   FILE *file;
   if (outFile.getSize () > 0)  {
      file = fopen(outFile.ascii(),"w");
      if (file == NULL)  {
         cout << "Can't open '" << outFile << "'." << endl;
         SX_EXIT;
      }
   } else {
      file = stdout;
   }
   
   cout << "Cell volume: " << structure.cell.volume << " bohr^3" << endl;

   if (refCenter) {
      atoms << -1;
   }

   int nAtoms = int(atoms.getSize ());
   if (nAtoms == 0)  {
      nAtoms = structure.getNAtoms ();
      for (int ia = 0; ia < nAtoms; ++ia) atoms << ia;
   }

   SxVector3<Int> gMesh 
      = SxGrid::suggestMesh (structure.cell, structure.getNAtoms () / 10 + 1);
   SxGrid grid (structure, gMesh);

   SxNeighbors neighbors, nextNeighbors;
   Coord tau, d1, d2;
   int nNeighbors, nNextNeighbors;

   // lattice symmetries (for reduced neighbors)
   SxArray<SymMat> latSym;
   // identify atoms for which the reduced neighbor list has been
   // printed
   SxArray<bool> eqAtoms(structure.getNAtoms ());
   eqAtoms.set (false);
   SxSymGroup fullSym;

   if (printReduced)  {
     fullSym  = structure.getSymGroup (SxAtomicStructure::FullSymmetries);
     fullSym.print ();
   }

   SxArray<int> equivalenceMap;
   if (printReduced)  {
      latSym = structure.cell.getLatticeSymmetries ();
      SxSymGroup (latSym).fprintsx (reducedFile);
      SxAtomicStructure primStr
         = structure.getPrimStr (structure.getPrimitiveCell ());
      primStr.updateSymGroup ();

      // Generate equivalence classes of primitive structure
      SxEquivalence equiv(primStr);
      equiv.fprintsx (reducedFile);

      // --- now set up mapping from primitive to full structure
      equivalenceMap = equiv.mapToComplex (structure).equivId;
   }

   FILE* bondFile = NULL;
   if (bondList.getSize () > 0)  {
      bondFile = sxfopen (bondList, "w");
   }

   for (int i = 0; i < nAtoms; ++i)  {
      int ia = atoms(i);
      if (refCenter && ia == -1)  {
         sxfprintf(file,"center @ (%f %f %f):\n",center(0),center(1),center(2));
         neighbors.compute (grid, structure, center, 
                            maxDist, SxNeighbors::StoreAll + SxNeighbors::IncludeZeroDistance);
      } else {
         if (ia >= structure.getNAtoms () || ia < 0)  {
            cout << "Ignoring illegal atom number " << (ia+1) << endl;
            continue;
         }
         fprint(file,structure,chemNames,ia);
         sxfprintf(file, " @ (%f %f %f): \n",
                   structure(ia)(0),structure(ia)(1),structure(ia)(2));
         neighbors.compute (grid, structure, structure(ia), 
                            maxDist, SxNeighbors::StoreAll);
      }
      nNeighbors = neighbors.absPositions.getNAtoms ();
      for (int ja = 0; ja < nNeighbors; ++ja)  {
         double d = sqrt(neighbors.distSqr(ja));
         if (d >= minDist)  {
            sxfprintf(file, "%i) ", ja+1);
            tau = neighbors.absPositions(ja);
            fprint(file,structure,chemNames,neighbors.idx(ja));
            sxfprintf(file, " @ (%f %f %f)", tau(0),tau(1),tau(2));
            sxfprintf(file, ": d= %f\n", d);

            if (bondFile && ia <= neighbors.idx(ja))  {
               // print "internal" bond if not across boundary
               int jTl = neighbors.idx(ja);
               if ((tau - structure(jTl)).normSqr () < 1e-4)
                  fprintf (bondFile, "%d %d\n", ia, jTl);
            }
         }
      }
      sxfprintf(file,"\n");
      if (angles)  {
         // angle ja - ia - ka
         for (int ja = 0; ja < nNeighbors; ++ja)  {
            for (int ka = ja+1; ka < nNeighbors; ++ka)  {
               sxfprintf(file, "%i(", ja+1);
               fprint(file,structure,chemNames,neighbors.idx(ja));
               sxfprintf(file,")-.-%i(", ka+1);
               fprint(file,structure,chemNames,neighbors.idx(ka));
               sxfprintf(file, "): angle= %f\n", 
                       getAngle(neighbors.relPositions(ja),
                                neighbors.relPositions(ka)) * 180./PI );
            }
         }
      }
      if (dihedral.contains(ia))  {
         // dihedral angles  ka - ia - ja - la
         int jaIdx, kaIdx, laIdx;
         Coord b;
         for (int ja = 0; ja < nNeighbors; ++ja)  {
            // central bond
            b = neighbors.relPositions (ja);
            nextNeighbors.compute (grid, structure, neighbors.absPositions(ja),
                                   maxDist, SxNeighbors::StoreRel);
            nNextNeighbors = nextNeighbors.relPositions.getNAtoms ();

            b.normalize ();
            jaIdx = neighbors.idx(ja);
            for (int ka = 0; ka < nNeighbors; ++ka)  {
               if (ka == ja) continue;
               // bond starting from current atom
               d1 = neighbors.relPositions (ka);
               d1 -= (d1 ^ b) * b; // project out bond component
               d1.normalize ();
               kaIdx = neighbors.idx (ka);
               for (int la = 0; la < nNextNeighbors; ++la)  {
                  laIdx = nextNeighbors.idx (la);
                  
                  // bond starting from neighbor atom
                  d2 = nextNeighbors.relPositions (la);
                  // check that neighbor of neighbor isn't starting point
                  if ((d2 + neighbors.relPositions (ja)).normSqr () < 1e-6)
                     continue;
                  d2 -= (d2 ^ b) * b; // project out bond component
                  d2.normalize ();
                  
                  // printout
                  sxfprintf(file,"%i:%i(", ja+1, la+1);
                  fprint(file,structure,chemNames, laIdx);
                  sxfprintf(file, ")-%i(", ja+1);
                  fprint(file,structure,chemNames,jaIdx);
                  sxfprintf(file,")-.-");
                  sxfprintf(file, "%i(", ka+1);
                  fprint(file,structure,chemNames,kaIdx);
                  sxfprintf(file, "): dihedral angle= %f\n", 
                          acos(d1 ^ d2) * 180. / PI);
               }
            }
         }
      }
      if (printReduced && (refCenter || !eqAtoms(ia)))  {
         if (!refCenter)  {
            if ((ia == 0 || 
                     structure.getISpecies (ia) != structure.getISpecies (ia-1))
                  && (reducedFile != stdout))
            {
               if (ia > 0) fprintf (reducedFile, "   }\n");
               fprintf (reducedFile, "   species {\n      element=\"%s\";\n\n", 
                        chemNames(structure.getISpecies (ia)).ascii ());
            }
         }
         if (nNeighbors == 0) continue;

         int nSym = int(latSym.getSize ());
         SxAtomicStructure shiftedStr = structure - (refCenter ? center
                                                               : structure(ia));

         // --- find local symmetries
         SxList<int> syms;
         for (int iSym = 0; iSym < nSym; ++iSym)  {
            if ( (latSym(iSym) ^ shiftedStr) == shiftedStr)
               syms << iSym;
         }

         SxArray<bool> found (nNeighbors);
         found.set (false);
         for (int ja = 0; ja < nNeighbors; ++ja)  {
            if (found(ja)) continue;
            // known interaction as ja->ia
            if (eqAtoms(neighbors.idx(ja))) continue;
            SxList<int> symId;
            SxList<int> matSymId;
            for (int iSym = 0; iSym < syms.getSize (); ++iSym)  {
               // rotate neighbor by local symmetry
               Coord rotN = latSym(syms(iSym)) ^ neighbors.relPositions(ja);

               // --- find rotated neighbor
               int ka;
               for (ka = 0; ka < nNeighbors; ++ka)
                  if ((rotN - neighbors.relPositions(ka)).normSqr () < 1e-6)
                     break;
               if (ka == nNeighbors)  {
                  double d = sqrt(neighbors.distSqr (ja));
                  cout << "Numerical noise problem: local symmetry "
                       << latSym(syms(iSym)) << endl
                       << "maps neighbor " << (ja + 1)
                       << " @ " << neighbors.relPositions (ja)
                       << " from center (d=" << d << ") to " << rotN
                       << " (absolute: " << (structure(ia) + rotN) << ")."
                       << endl << "There is no neighbor!" << endl;
                  if (fabs( (d - maxDist)) < 3. * structure.epsEqual)  {
                     cout << "Atom sits right on cutoff sphere. "
                          << "The missing neighbor might be slightly outside. "
                          << endl << "Increasing the cutoff slightly may help."
                          << endl;
                     SX_QUIT;
                  }
                  SX_EXIT;
               }

               if (!found(ka))  {
                  symId << (syms(iSym)+1);
                  found(ka) = true;
               }
               if (ka == ja) matSymId << (syms(iSym) + 1);
            }
            // --- print out
            if (neighbors.relPositions(ja).norm () >= minDist)  {
               sxfprintf (reducedFile, "      neighbor {\n");
               if (!refCenter)  {
                  sxfprintf (reducedFile, "         // relative to\n");
                  sxfprintf (reducedFile, "         equivalenceId = %d;\n",
                             equivalenceMap(ia) + 1);
               }
               sxfprintf (reducedFile, "         // distance = % .6f\n", neighbors.relPositions(ja).norm ());
               sxfprintf (reducedFile, "         coords = [% .6f, % .6f, % .6f];\n",
                     neighbors.relPositions(ja)(0),
                     neighbors.relPositions(ja)(1),
                     neighbors.relPositions(ja)(2));
               sxfprintf (reducedFile, "         atomType=\"%s\";\n",
                     chemNames(structure.getISpecies (neighbors.idx (ja)))
                     .ascii ());

               // --- print symmetry list
               sxfprintf (reducedFile, "         // number of neighbors = %i\n",
                         int(symId.getSize ()));
               if (symId.getSize () > 1)  {
                  sxfprintf (reducedFile, "         symmetries = [%d", symId(0));
                  for (int iSym = 1; iSym < symId.getSize (); ++iSym)  {
                     if (iSym % 10) 
                        sxfprintf (reducedFile, ", %d", symId(iSym));
                     else
                        sxfprintf (reducedFile, ",\n                       %d", symId(iSym));
                  }
                  sxfprintf (reducedFile, "];\n");
               }
               if (matSymId.getSize () > 1)  {
                  sxfprintf (reducedFile, "         matSyms    = [%d", matSymId(0));
                  for (int iSym = 1; iSym < matSymId.getSize (); ++iSym)  {
                     if (iSym % 10) 
                        sxfprintf (reducedFile, ", %d", matSymId(iSym));
                     else
                        sxfprintf (reducedFile, ",\n                       %d", matSymId(iSym));
                  }
                  sxfprintf (reducedFile, "];\n");
               }
               sxfprintf (reducedFile, "      }\n");
            }
         }
         if (!refCenter)  {
            // --- mark all equivalent atoms
            for (int iSym = 0; iSym < fullSym.getSize (); ++iSym)  {
               int ja = structure.find (fullSym(iSym) ^ structure(ia), grid);
               eqAtoms(ja) = true; 
            }
         }
      }
   }

   
   if (file != stdout) fclose (file);
   if (reducedFile != stdout) {
      fprintf (reducedFile, "   }\n");
      fclose (reducedFile);
   }
   if (bondFile) fclose (bondFile);

}

