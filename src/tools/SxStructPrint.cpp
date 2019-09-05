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
#include <SxStickyFilter.h>
#include <SxPDBFast.h>
#include <SxPoscar.h>
#include <SxAims.h>

SxAtomicStructure readRaw (FILE *fp)
{
   SX_CHECK(fp);
   char buffer[10240];

   SxAtomicStructure res;
   res.startCreation ();

   enum Mode {None, StructStream, VectorStream, Xyz} mode = None;
   Coord pos;
   int lineNo = 0;
   int iSpecies;
   while (!feof(fp))  {
      if (!fgets (buffer, 10240, fp)) break;
//      cout << buffer;
      ++lineNo;
      if (buffer[0]=='{')  {
         if (mode == None) res.newSpecies ();
         cout << "Vector stream mode detected!" << endl;
         mode = VectorStream;
         char *where = buffer+1;
         int nChar;
         do {
            nChar = 0;
            if (sscanf(where,"{%lf},{%lf},{%lf}%n",
                           &pos(0), &pos(1), &pos(2), &nChar) >= 3)  {
               res.addAtom (pos);
               where += nChar;
            } else {
               cout << "Unexpected end in vector stream mode" << endl;
               break;
            }
         } while (*where++ == ',');
      } else if ((sscanf(buffer, "species %d", &iSpecies) == 1))  {
         res.newSpecies ();
         if (mode != StructStream)
            cout << "Structure stream mode detected." << endl;
         mode = StructStream;
      } else if (sscanf(buffer, "atom %*d: {%lf,%lf,%lf}",
                                &pos(0), &pos(1), &pos(2)) == 3) {
         if (mode != StructStream)
            cout << "Structure stream mode detected." << endl;
         if (mode == None) res.newSpecies ();
         mode = StructStream;
         res.addAtom (pos);
      } else if (sscanf(buffer,"%lf %lf %lf",&pos(0),&pos(1),&pos(2)) == 3)  {
         if (mode != Xyz)
            cout << "xyz style detected." << endl;
         if (mode == None) res.newSpecies ();
         mode = Xyz;
         res.addAtom(pos);
      } else if (sscanf(buffer, " %*s") > 0) {
         cout << "Failed to parse line no. " << lineNo << ':' << endl;
         cout << buffer << endl;
      } else {
         // empty line
      }
   }
   res.endCreation ();
   if (mode == VectorStream || mode == Xyz)
      res.atomInfo = SxPtr<SxAtomInfo> ();
   return res;
}

SxAtomicStructure readFHI98 (const SxString &fileName)
{
   SxBinIO io(fileName, SxBinIO::ASCII_READ_ONLY);
   FILE *fp = io.fp;

   const int BUFLEN = 1024;
   int i, ic;
   char buffer[BUFLEN];
   char *bs = NULL;
   SxString tk;
   SxMatrix3<Double> cell;

   for (i = 0; i < 3; i++)  {
         bs = fgets(buffer, BUFLEN, fp);
         tk = SxString(bs);
         SxList<SxString> list;
         list = tk.tokenize (' ');
         for (ic = 0; ic < 3; ic++) {
            try {
               cell(ic, i)=list(ic).toDouble ();
            } catch (SxException e)  {
               e.print ();
               SX_QUIT;
            }
         }
   }
   bs = fgets (buffer, BUFLEN, fp);
   tk = SxString(bs);
   int nSpecies = tk.toInt ();

   SxArray<SxString> chemName(nSpecies);
   SxAtomicStructure res;
   res.cell = cell;

   for (int is=0; is < nSpecies; is++) {
      bs = fgets (buffer, BUFLEN, fp);
      tk = SxString(bs);
      int nAtoms = tk.toInt ();

      bs = fgets (buffer, BUFLEN, fp);
      tk = SxString(bs);
      tk.resize (tk.getSize () - 1, true); // chop off '\n'
      chemName(is) = tk.trim ();

      res.newSpecies ();

      for (int ia=0; ia < nAtoms; ia++)  {
         bs = fgets(buffer, BUFLEN, fp);
         tk = SxString(bs);
         SxList<SxString> list = tk.tokenize (' ');
         Coord pos;
         for (ic = 0; ic < 3; ic++) {
            pos(ic)=list(ic).toDouble ();
         }
         res.addAtom (pos);
      }
   }
   res.endCreation ();
   res.atomInfo->meta.attach (SxAtomicStructure::Elements, chemName);
   return res;
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

   bool wrap = cli.option ("--wrap", "keep atoms within cell").toBool ();
   int WS = cli.newGroup ("Wigner-Seitz mapping");
   Coord center ( cli.option ("--center","point","center of Wigner-Seitz cell")
                  .toList3());
   cli.setGroup (cli.generalGroup);
   bool regular = cli.option ("--regular", "regularize cell")
                  .toBool ();

   int binGroup = cli.newGroup ("binary files");
   bool sxbFile = cli.option ("-b|--sxb", "input file is binary waves file "
                              "rather than S/PHI/nX input file")
                  .toBool ();
   int sxGroup = cli.newGroup ("sx files");
   cli.excludeGroup (binGroup);
   bool keepMovable
      = cli.option("-m|--keep-movable","retain mobility of atoms as "
                   "it is defined in the input structure")
        .toBool ();

   int pdbGroup = cli.newGroup ("pdb files");
   cli.excludeGroup (binGroup);
   cli.excludeGroup (sxGroup);
   bool pdbFile = cli.option ("--pdb","input file is PDB file").toBool ();

   int vaspGroup = cli.newGroup ("vasp files");
   cli.excludeGroup (sxGroup);
   cli.excludeGroup (pdbGroup);
   cli.excludeGroup (binGroup);
   bool poscarFile = cli.option ("--vasp","input file is VASP POSCAR file")
      .toBool ();

   int fhiGroup = cli.newGroup ("FHI98 data");
   cli.excludeGroup (sxGroup);
   cli.excludeGroup (pdbGroup);
   cli.excludeGroup (binGroup);
   cli.excludeGroup (vaspGroup);
   bool fhi98File = cli.option ("--fhi98", "input file is fhi98md format")
                    .toBool ();
   int fhiAimsGroup = cli.newGroup ("FHIAims data");
   cli.excludeGroup (sxGroup);
   cli.excludeGroup (pdbGroup);
   cli.excludeGroup (binGroup);
   cli.excludeGroup (vaspGroup);
   cli.excludeGroup (fhiGroup);
   bool fhiAimsFile = cli.option ("--aims", "input file is FHIAims format")
                    .toBool ();

   // this is a temporary feature. Must be removed when trajectories are
   // available
   cli.newGroup ("raw data");
   cli.excludeGroup (sxGroup);
   cli.excludeGroup (pdbGroup);
   cli.excludeGroup (binGroup);
   cli.excludeGroup (vaspGroup);
   cli.excludeGroup (fhiGroup);
   cli.excludeGroup (fhiAimsGroup);
#ifndef NDEBUG
   bool rawFile = cli.option ("--raw", "input file is raw data. Experimental.")
                  .toBool ();
#endif

   cli.setGroup(SxCLI::generalGroup);

   SxString inFile
      = cli.option ("-i|--input", "input file",
                    "take original input file")
        .toString (sxbFile ? "waves.sxb" : "input.sx");
   cli.last ().defaultValue = "default: input.sx, or waves.sxb for --sxb flag";

   SxString speciesFile
      = cli.option ("-s|--species", "input file",
                    "take species data from this S/PHI/nX input file")
        .toString ("");

   SxString outFile
      = cli.option ("-o","filename", "output file name (screen otherwise)")
        .toString ("");

   bool outSxb = cli.option ("--outsxb", "output file is binary "
                              "rather than S/PHI/nX input file").toBool ();
   if (outSxb && outFile.getSize () == 0) outFile = "structure.sxb";

   int printOptions = SxAtomicStructure::DefaultPrint;
   if (cli.option("--printsym", "print symmetries").toBool ())
      printOptions |= SxAtomicStructure::PrintSymmetries;

   cli.newGroup ("cut structure");
   bool cut = cli.option ("--cut", "cut Cell").toBool ();
   Coord lower (cli.option ("-l","lower","lower Bound").toList3 ("x,"));
   Coord upper (cli.option ("-u","upper","upper Bound").toList3 ("x,"));

   cli.finalize ();

   // --- read input
   SxAtomicStructure structure;
   if (sxbFile)  {
      try {
         SxBinIO io (inFile, SxBinIO::BINARY_READ_ONLY);
         structure.read (io);
         if (speciesFile.getSize () == 0)  {
         }
         io.close ();
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
   } else if (pdbFile) {
      SxPDBFast pdb(inFile);
      pdb.read ();
      structure = pdb.getStructure ();
   } else if (poscarFile) {
      SxPoscar poscar(inFile);
      poscar.read ();
      structure = poscar.getStructure ();
      structure.atomInfo->meta.update (SxAtomicStructure::Elements,
                                       poscar.getUniqueSpecies ());
   } else if (fhi98File)  {
      structure = readFHI98 (inFile);
   } else if (fhiAimsFile)  {
      SxAims geometry(inFile);
      geometry.read();
      structure = geometry.getStructure ();
#  ifndef NDEBUG
   } else if (rawFile)  {
      FILE *fp = fopen(inFile.ascii (), "r");
      if (fp)  {
         structure = readRaw (fp);
      } else {
         cout << "Cannot open " << inFile << "for reading. " << endl;
         SX_QUIT;
      }
#  endif
   } else {
      SxParser parser;
      SxConstPtr<SxSymbolTable> tree;
      tree = parser.read (inFile, "std/structure.std");
      structure = SxAtomicStructure (&*tree);
      if (speciesFile.getSize () == 0) structure.readElements (&*tree);
      if (keepMovable)
         structure.atomInfo->meta.attach (SxAtomicStructure::StickyFilter,
               SxStickyFilter(tree->getGroup("structure")).getStickyArray ());
   }

   if (speciesFile.getSize () > 0)  {
      SxParser parser;
      SxConstPtr<SxSymbolTable> tree;
      tree = parser.read (speciesFile, "std/structure.std");
      SxArray<SxString> chemNames = SxSpeciesData::getElements (&*tree);
      if (structure.cell.volume < 1e-7)  {
         // structure probably from raw data. Feed in our cell...
         cout << "WARNING: Taking cell from '" << speciesFile << "'.\n";
         structure.cell = SxCell(&*tree);
      }
      if (!structure.atomInfo)  {
         SxAtomicStructure specStr(&*tree);
         if (specStr.getNAtoms () == structure.getNAtoms ())  {
            structure.atomInfo = specStr.atomInfo;
         } else {
            cout << "Cannot transfer species info: number of atoms differs.\n";
            cout << inFile << ": nAtoms=" << structure.getNAtoms () << endl;
            cout << speciesFile << ": nAtoms=" << specStr.getNAtoms () << endl;
            SX_QUIT;
         }
      }
      if (chemNames.getSize () != structure.getNSpecies ())  {
         cout << SX_SEPARATOR;
         cout << "| WARNING: number of species differs" << endl;
         cout << "| nSpecies=" << structure.getNSpecies ()
              << " in '" << inFile << "'." << endl;
         cout << "| nSpecies=" << chemNames.getSize ()
              << " in '" << speciesFile << "'." << endl;
         cout << SX_SEPARATOR;
         chemNames.resize (structure.getNSpecies (), true);
      }
      structure.atomInfo->meta.update (SxAtomicStructure::Elements, chemNames);
   } else {
      if (!structure.atomInfo)  {
         SxPtr<SxAtomInfo> info = SxPtr<SxAtomInfo>::create (1);
         info->nAtoms(0) = structure.getNAtoms ();
         info->setupOffset ();
         structure.replaceInfo (info);
      }
   }

   if (regular) structure.cell = structure.cell.getRegularCell ();

   if (cli.groupAvailable (WS))  {
      structure -= center;
      structure.map (SxCell::WignerSeitz);
      structure += center;
   } else if (wrap) {
      structure %= structure.cell;
   }

   if (cut)  {
      structure %= structure.cell;
      structure = structure.cut(lower, upper);
   }

   // --- output
   if (outSxb)  {
      structure.write (outFile);
   } else {
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

      structure.fprint (file, printOptions);

      if (file != stdout) fclose (file);
   }

}
