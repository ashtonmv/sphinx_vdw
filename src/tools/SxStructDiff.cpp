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

void fprint(FILE *file, const SxAtomicStructure &str,
            const SxArray<SxString>& symbol, int ia)
{
   SX_CHECK (ia < str.getNAtoms () && ia >= 0, ia, str.getNAtoms ());
   int is = 0, nAtoms;
   while (ia >= (nAtoms = str.getNAtoms(is++)) )  ia -= nAtoms; 

   sxfprintf(file, "%s%i", symbol(is).ascii (), (ia+1));

}

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.preUsageMessage =
      "This add-on prints the structure.";
   cli.authors = "C. Freysoldt";

   bool sxbFile = cli.option ("-b|--sxb", "input file is binary waves file "
                              "rather than S/PHI/nX input file")
                  .toBool ();
   SxString inFile
      = cli.option ("-i|--input", "input file",
                    "take original input file")
        .toString (sxbFile ? "waves.sxb" : "input.sx");
   cli.last ().defaultValue = "default: input.sx, or waves.sxb for --sxb flag";

   SxString refFile
      = cli.option ("-r|--reference","filename","reference file")
        .toString ("input.sx");

   double maxDist
      = cli.option ("-d|--maxdist", "distance", "max. distance (bohr) "
                    "considered a displacement").toDouble (1., 0.);

   bool printDist
      = cli.option ("--printdist", "print shift distances").toBool ();

   cli.newGroup ("boundaries");
   bool boundaryShift
      = cli.option ("--boundaryShift", "repeat shift for atoms at boundary "
                    "for all boundaries").toBool ();
   cli.option ("--boundaryCenter", "set center for boundary");
   bool bCenterShift;
   SxVector3<Double> bCenter(0,0,0);
   if ( (bCenterShift = cli.last ().exists ()) )
      bCenter = SxVector3<Double> (cli.last ().toList3 ());
   cli.last ().optional = true;

   cli.setGroup (cli.generalGroup);

   bool driftCheck
      = cli.option ("--drift", 
                    "calculate average drift for all displaced atoms").toBool();

   cli.finalize ();

   // --- read input
   SxAtomicStructure structure;
   SxArray<SxString> chemNames;
   SxArray<SxVector3<Int> > movableInfo;
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
      chemNames = SxSpeciesData::getElements(&*tree);
   }
   
   SxAtomicStructure refStructure;
   SxArray<SxString> refChemNames;
   {
      // --- read in reference structure
      SxParser parser;
      SxConstPtr<SxSymbolTable> tree;
      tree = parser.read (refFile, "std/structure.std");
      refStructure = SxAtomicStructure (&*tree);
      refChemNames = SxSpeciesData::getElements(&*tree);
   }

   // --- synchronize species ordering in structure and reference
   {
      int nSpecies = refStructure.getNSpecies ();
      SxList<int> specMap;
      SxList<SxString> specName;
      // copy reference structure species
      for (int is = 0; is < nSpecies; ++is) specName << refChemNames(is);

      // --- map species (structure) -> species (refStructure)
      ssize_t isInRef;
      for (int is = 0; is < structure.getNSpecies (); ++is)  {
         if (chemNames(is).getSize () > 0)  {
            // identify species in reference structure
            isInRef = specName.findPos (chemNames(is));
            if (isInRef < 0)  {
               // not found: new species
               isInRef = nSpecies++;
               specName << chemNames(is);
            }
         } else {
            // new species
            isInRef = nSpecies++;
            specName << SxString(nSpecies,"[species %d]");
         }
         specMap << int(isInRef);
      }
      chemNames = specName;

      if (nSpecies > refStructure.getNSpecies ())  {
         // --- reference structure: extend species list
         SxPtr<SxAtomInfo> info = SxPtr<SxAtomInfo>::create ();
         info->resize (nSpecies);
         info->nAtoms(SxIdx(0, refStructure.getNSpecies ()-1))
               <<= refStructure.atomInfo->nAtoms;
         info->setupOffset ();
         refStructure.atomInfo = info;
      }
         
      bool reorder = false;
      SxList<int>::Iterator it = specMap.begin ();
      for (int i = 0; i < structure.getNSpecies (); ++i, ++it)
         if ( (reorder = (*it != i)) ) break;

      if (reorder)  {
         // --- input structure: reorder atoms
         SxPtr<SxAtomInfo> info = SxAtomInfo::derive(structure.atomInfo);
         info->resize (nSpecies);
         it = specMap.begin ();
         for (int is = 0; is < structure.getNSpecies (); ++is, ++it)
            info->nAtoms(*it) = structure.atomInfo->nAtoms(is);
         info->setupOffset ();
         
         // create reordering map
         info->parentMap.resize (structure.getNAtoms ());
         it = specMap.begin ();
         SxVector<Int>::Iterator mapIt;
         for (int is = 0; is < structure.getNSpecies (); ++is, ++it)  {
            mapIt = info->parentMap.begin ();  
            mapIt += info->offset(*it);
            int offset = structure.atomInfo->offset(is);
            for (int i = 0; i < structure.getNAtoms (is); ++i, ++mapIt)
               *mapIt = i + offset;
         }
         // use info to remap structure
         structure = SxAtomicStructure(info, structure);
      }
   }
   if (bCenterShift) bCenter -= refStructure.cell.relToCar (Coord(0.5,0.5,0.5));
   if (boundaryShift)  {
      if (bCenterShift) refStructure -= bCenter;
      refStructure.map (SxCell::Positive);
      if (bCenterShift) refStructure += bCenter;
   }

   int nDrift = 0;
   Coord drift(0.,0.,0.); 

   SxGrid grid(refStructure, 10);
   SxNeighbors nn;

   // --- now loop over atoms in structure and look what has happened to them
   // IMPORTANT NOTE: Lines starting with "> " are parsed by SxStructPatch
   //                 DO NOT CHANGE the format without adapting SxStructPatch
   int ja = 0, is;
   SxArray<bool> refFound(refStructure.getNAtoms ());
   refFound.set (false);
   for (int ia = 0; ia < structure.getNAtoms (); ++ia)  {
      is = structure.getISpecies (ia);
      if ((ja = refStructure.find (structure(ia), grid)) >= 0) {
         // atom position coincides
         int js = refStructure.getISpecies(ja);
         if (refFound(ja))
            cout << "Warning: multiple match for atom " << (ja+1)
                 << " at " << refStructure(ja) << endl;
         refFound(ja) = true;
         // check species
         if (is != js)
            cout << "> atom " << (ja+1) << " (" << chemNames(js)
                 << ") @ " << refStructure(ja) 
                 << ": new species " << chemNames(is) << endl;

      } else {
         Coord shift;
         // look in close area
         nn.compute (grid, refStructure, structure(ia), maxDist,
                     SxNeighbors::StoreRel | SxNeighbors::StoreDistSqr
                     | SxNeighbors::IncludeZeroDistance);
         
         // fetch closest neighbor of same species
         double dMin2 = maxDist*maxDist*1.1; // little overdoing (numerics)
         for (int in = 0; in < nn.getSize (); ++in)  {
            if (dMin2 > nn.distSqr (in) 
                && refStructure.getISpecies(nn.idx(in)) == is)
            {
               ja = nn.idx(in);
               dMin2 = nn.distSqr (in);
               shift = -nn.relPositions(in);
            }
         }
         
         if (ja >= 0)  {
            // atom displaced
            if (refFound(ja))
               cout << "Warning: multiple match for atom " << (ja+1)
                    << " @ " << refStructure(ja) << endl;
            refFound(ja) = true;
            cout << "> atom " << (ja+1) << " (" << chemNames(is) << ") @ "
                 << refStructure(ja) << ": shift " << shift << endl;
            if (driftCheck)  {
               nDrift++;
               drift += shift;
            }
            if (boundaryShift)  {
               Coord relPos = refStructure.cell.carToRel (refStructure(ja) - bCenter);
               SxVector3<Int> from, to,R;
               from.set (0);
               to.set (0);
               // --- find boundaries
               for (int d = 0; d < 3; ++d)  {
                  if (fabs(relPos(d)) < 1e-4) to(d) = 1;
                  if (fabs(relPos(d)-1.) < 1e-4) from(d) = -1;
               }
               // --- replicate boundary atoms
               for (R(0) = from(0); R(0) <= to(0); ++R(0))  {
                  for (R(1) = from(1); R(1) <= to(1); ++R(1))  {
                     for (R(2) = from(2); R(2) <= to(2); ++R(2))  {

                        if (R.absSqr ().sum () == 0) continue;

                        Coord Rr = refStructure.cell.relToCar (R);
                        cout << "> atom " << (ja+1) << " ("
                             << chemNames(is) << ") @ "
                             << (refStructure(ja) + Rr)
                             << ": shift " << shift << endl;
                     }
                  }
               }
            }
            if (printDist)
               cout << "d=" << sqrt(dMin2) << endl;
         } else {
            // atom added
            cout << "> new " << chemNames(structure.getISpecies(ia))
                 << " @ " << structure(ia) << endl;
         }
      }
   }

   // --- list remaining atoms in refStructure
   for (int ia = 0; ia < refStructure.getNAtoms (); ++ia)  {
      is = refStructure.getISpecies (ia);
      if (!refFound(ia))  {
         cout << "> atom " << (ia+1) << " (" << chemNames(is) << ") @ " 
              << refStructure(ia) << ": deleted." << endl;
      }
   }

   if (driftCheck && nDrift > 0)  {
      cout << "drift (" << nDrift << " atoms) = "
           << (drift / double(nDrift)) << endl;
   }

}

