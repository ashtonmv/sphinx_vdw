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

#include <SxWavesJoin.h>
#include <SxCLI.h>

#ifndef SX_STANDALONE

SxWavesJoin::SxWavesJoin (const SxKPoints &kPoints,
                          const SxList<SxString> &fileNames)
{
   SxPW::gkBasisPtr = SxThis<SxWavesJoin>::getThis ();

   // --- set up a few members and auxiliaries
   double nElectrons = -1.;

   // --- set up k-points
   SxKPoints::operator= (kPoints);
   nk = nkPoints = SxKPoints::getNk ();
   SX_CHECK (nk > 0);
   int ik, jk, iFile, iSpin;
   nGkMax = -1;

   gBasisList.resize (nk);
   for (ik = 0; ik < nk; ++ik)
      gBasisList(ik) = SxPtr<SxGBasis>::create ();

   kMap.resize (nk);
   nPerK.resize (nk);
   nGIPerK.resize (nk);

   // --- open the input files and find the k-points
   files.resize (fileNames.getSize ());
   for (iFile = 0; iFile < files.getSize (); ++iFile)  {
      try {
         // open file
         files(iFile).open (fileNames(iFile), SxBinIO::BINARY_READ_ONLY);

         // --- Read k-point data
         SxKPoints kInFile;
         kInFile.read (files(iFile));
         SxVector<Int> nGk(kInFile.getNk ());
         files(iFile).read ("nGk", &nGk, kInFile.getNk ());
         SxVector<Int> nPerKFile(kInFile.getNk ());
         files(iFile).read ("nPerK", &nPerKFile, kInFile.getNk ());

         // --- Read the fermi data
         SxFermi fermiFile;
         fermiFile.read (files(iFile));
         if (iFile == 0)  {
            nSpin = files(0).getDimension ("nSpin");
            waves.resize (1);
            waves(0).resize (nSpin);
            fermi = SxFermi (nElectrons, nPerKFile(0), nSpin, *this);
         }

         // Set up the map
         for (ik = 0; ik < nk; ++ik)  {
            for (jk = 0; jk < kInFile.getNk (); ++jk)  {
               if ((getK(ik) - kInFile.getK(jk)).normSqr () < 1e-7)  {
                  /// Found a k-point
                  if (kMap(ik).iFile >= 0)  {
                     cout << "Warning: found k-point " << (ik+1)
                          << " (k=" << getK(ik) << ") multiple times."
                          << endl;
                     cout << "Taking data from '" << fileNames(kMap(ik).iFile)
                          << "'.\n";
                  } else {
                     cout << "Taking ik=" << (ik + 1) << " ("
                          << kVec(ik) << ") from '" << fileNames(iFile)
                          << "'.\n";

                     // --- set up the map
                     kMap(ik).iFile = iFile;
                     kMap(ik).ik = jk;
                     nPerK(ik) = nPerKFile(jk);
                     nGIPerK(ik) = nPerK(ik) * nGk(jk);

                     // --- set up the G+k basis
                     int start = (jk == 0) ? 0 : nGk.sum (0, jk-1);
                     gBasisList(ik)->read (files(iFile), nGk(jk), start);

                     // setup nGkMax (needed by SxPW::writeWavesFile)
                     if (nGkMax < nGk(jk)) nGkMax = nGk(jk);

                     // --- copy the fermi data
                     for (iSpin = 0; iSpin < nSpin; ++iSpin)  {
                        fermi.eps(iSpin, ik)  = fermiFile.eps (iSpin, jk);
                        fermi.focc(iSpin, ik) = fermiFile.focc (iSpin, jk);
                     }

                  }
               }
            }
         }
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
   }

   // --- check that all k-points are available
   for (ik = 0; ik < nk; ++ik)  {
      if (kMap(ik).iFile == -1)  {
         cout << "Couldn't find k-point " << (ik+1) << " (k="
              << getK(ik) << ")!\n";
         SX_QUIT;
      }
   }

   if (nPerK.minval () != nPerK.maxval ())  {
      cout << "Varying number of states is not supported!" << endl;
      SX_QUIT;
   }

}

void SxWavesJoin::write (const SxString &filename)
{
   SX_CHECK (files.getSize () > 0);
   SX_CHECK (SxKPoints::getNk () > 0);
   SxAtomicStructure structure;
   gBasisList(0)->fft3d.resize (1);
   try {
      files(0).read ("meshDim", &gBasisList(0)->fft3d(0).mesh);
      files(0).read ("eCut", &gkCut);
      structure.read (files(0));
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   /// Use standard write routine, which uses operator(i,iSpin,ik) as call-back
   writeWavesFile (filename, fermi, structure, false);
}

const SxGBasis::TPsi SxWavesJoin::operator () (int i, int iSpin, int ik) const
{
   SX_CHECK (i >= 0 && i < getNStates(ik), i, getNStates(ik));
   SX_CHECK (iSpin >= 0 && iSpin < nSpin, iSpin, nSpin);
   SX_CHECK (ik >= 0 && ik < nk, ik, nk);

   /*
   cout << "Read (" << i << "," << iSpin << "," << ik << ") from "
        << files(kMap(ik).iFile).filename << ": ik=" << kMap(ik).ik << endl;
   */
   return SxPW::readPsi (files(kMap(ik).iFile), i, iSpin, kMap(ik).ik);
}

int SxWavesJoin::getNStates (int ik) const
{
   return nPerK(ik);
}

#else /* SX_STANDALONE */


int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   cli.preUsageMessage = SxString (
      "This add-on composes a single waves file from k-points scattered "
      "over different waves files."
      "It's main use is to enable manual k-point parallelisation for "
      "bandstructure-type runs.\n"
      "This is an expert tool, so expect little support and "
      "no robustness against errors. It is up to you to make sure that "
      "the provided files contain all necessary data and are consistent "
      "(same system, same nStates, same spin, ...). There's a good "
      "chance that errors will cause a hard crash (seg fault) "
      " - I am to lazy to implement all the necessary checks."
      ).wrap ();
   cli.authors="C. Freysoldt";

   SxString input = cli.option ("-i|--input", "input file", "all-k input file")
                    .toString ("input.sx");

   SxString output = cli.option ("-o|--output", "output file",
                                 "output file for all k-points from input file"
                                 " with data from waves files")
                     .toString ("waves-joined.sxb");
   SxList<SxString> files = cli.argument ("waves files", "wave function files")
                            .toList ();
   cli.finalize ();

   SxParser::Table table;
   try {
      SxParser parser;
      table = parser.read (input);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   SxAtomicStructure str(&*table);
   SxKPoints kPoints (str.cell, &*table);

   SxWavesJoin joiner (kPoints, files);

   joiner.write (output);

   return 0;

}


#endif /* SX_STANDALONE */
