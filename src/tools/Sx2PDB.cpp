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

#include <SxUtil.h>
#include <SxConstants.h>
#include <SxList.h>
#include <SxPrecision.h>
#include <SxVector3.h>
#include <SxMatrix3.h>
#include <SxString.h>
#include <SxPDBFast.h>
#include <SxCLI.h>
#include <SxBinIO.h>
#include <SxAtomicStructure.h>
#include <SxSpeciesData.h>

int main (int argc, char **argv)
{
   // --- parse command line options
   SxCLI cli (argc, argv);
   SxCLI::CliArg* opt;
   
   // --- input and output files
   SxString inFile 
      = cli.option("-i","file","SPHInX input file to read from")
        .toString("input.sx");
   SxString wavesFile
      = cli.option ("-w","file","if given, read from SPHInX waves file "
                    "(waves.sxb) rather than from input file")
        .toString ("");
   SxString outFile 
      = cli.option("-o","file","PDB output file to be written")
        .toString("input.pdb");
   
   // --- elements
   SxList<SxString> elements; 
   opt = &(cli.option("-elem","list","overwrites the default list of chemical "
                      "elements with user-given ':'-separated list.\n"
                      "-elem H:Ga:Mg means the 1st species is assumed to be "
                      "Hydrogen, the 2nd one Gallium and the 3rd Magnesium"));
   if (opt->exists (true))
      elements = (opt->toString ("")).tokenize (':');
   opt = &(cli.option("--colourful","overwrites the default list of chemical "
                  "elements with O:B:Na:S:P:Mg which equals red:green:blue:"
                  "yellow:oraSxe:forest green in the CPK model"));
   if (opt->toBool () && elements.getSize () == 0)
      elements << "O" << "B" << "Na" << "S" << "P" << "Mg";
   
   // --- repetition
   int xRepeat, yRepeat, zRepeat;
   SxString xyzRepeat
      = cli.option("-r|--repeat","AxBxC",
                   "with numbers A, B, C, and x beiSx an 'x'. "
                   "Repeat the cell A times in first (x) direction, B times in "
                   "second (y) direction and C times in 3rd (z) direction")
        .toString ("1x1x1");
   if (xyzRepeat.tokenize('x').getSize () == 3)  {
      SxList<SxString> xyz = xyzRepeat.tokenize ('x');
      xRepeat = xyz(0).toInt ();
      yRepeat = xyz(1).toInt ();
      zRepeat = xyz(2).toInt ();
   } else {
      cout << "Illegal format: -r " << xyzRepeat << endl;
      SX_QUIT;
   }
   cli.finalize ();
   
   SxMatrix3<TReal8> aMat;
   SxVector3<Double> coords;
   SxAtomicStructure structure;
   SxArray<SxString> chemName;
   
   if (wavesFile.getSize () == 0)  {
      // --- read from input file
      initSPHInXMath ();
      SxParser parser;
      SxParser::Table table = parser.read (inFile, "std/structure.std");
      
      try  {
         structure = SxAtomicStructure(&*table);
         chemName = SxSpeciesData::getElements (&*table);
      } catch (SxException e)  {
         e.print ();
         SX_QUIT;
      }

      
      } else {

         cout << "Read in from waves.sxb not implemented yet" << endl;
      
         //--- TODO: read in structure from waves-file
         // --- read from waves file
      /*
         try  {

         
         SxBinIO io(wavesFile, SxBinIO::BINARY_READ_ONLY);
         
         // read cell
         io.read ("cell", &aMat);
         
         // --- read tau
         nSpecies = io.getDimension ("nSpecies");
         tau.resize (nSpecies);
         int nTotalAtoms = io.getDimension("nAllAtoms");
         io.read ("nAtoms", &nAtoms, nSpecies);
         SxMatrix<Double> tauMatrix (nTotalAtoms, 3);
         io.read ("tau", &tauMatrix, nTotalAtoms, 3);
         int iTotAtoms = 0;
         for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
            for (iAtom = 0; iAtom < nAtoms(iSpecies); iAtom++, iTotAtoms++)  {
               tau(iSpecies).append (
                     SxVector3<Double> (tauMatrix(iTotAtoms,0),
                                        tauMatrix(iTotAtoms,1),
                                        tauMatrix(iTotAtoms,2)) );
            }

         // --- read elements
         if (readElements)  {
            SxString chemNames;
            io.read ("chemNames", &chemNames);
            elements = chemNames.tokenize (',');
         }
         io.close ();
      } catch (SxException e)  {
         e.print ();
         SX_QUIT;
      }
      */
   }

   if (xRepeat != 1 || yRepeat != 1 || zRepeat != 1)
      structure = structure.repeat (xRepeat, yRepeat, zRepeat);
   SxVector3<Double> coord, trans;
   SxPDBFast(outFile).write (structure, chemName);
   return 0;
}

