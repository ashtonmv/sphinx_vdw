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

#include <SxWavesCmp.h>
#include <SxFermi.h>
#include <SxAtomicStructure.h>
#include <SxCLI.h>
#include <SxBinIO.h>
#include <SxList.h>

#ifndef SX_STANDALONE

SxWavesCmp::SxWavesCmp ()
{
   // empty
}


SxWavesCmp::~SxWavesCmp ()
{
   // empty
}


void SxWavesCmp::readWaves1 (const SxString &file)
{
   try  {
      SxBinIO io (file, SxBinIO::BINARY_READ_ONLY);
      Gk.read (io);
      waves1.read (io);
      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}


void SxWavesCmp::readWaves2 (const SxString &file)
{
   try  {
      SxBinIO io (file, SxBinIO::BINARY_READ_ONLY);

      sxprintf ("TODO: Should check that |G+k> are the same\n");
      sxprintf ("      SxGkBasis::operator== perhaps???\n");

      waves2.read (io);
      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}




SxPW SxWavesCmp::getDifferences ()
{
   //TODO Check for identical Gk Basis of waves1 and waves2
   int i, nStates   = waves1.getNStates ();
   int iSpin, nSpin = waves1.getNSpin ();
   int ik, nk       = waves1.getNk ();
   SxPW res (nStates, nSpin, waves1.getGkBasisPtr());
   for (ik=0; ik < nk; ++ik)  {
      for (iSpin=0; iSpin < nSpin; ++iSpin)  {
         for (i=0; i < nStates; ++i)  {
            res(iSpin,ik).colRef(i) <<= waves1(iSpin,ik).colRef(i)
                                      - waves2(iSpin,ik).colRef(i);
            VALIDATE_VECTOR (waves1(i,iSpin,ik));
            VALIDATE_VECTOR (waves2(i,iSpin,ik));
            VALIDATE_VECTOR (res(i,iSpin,ik));
         }
      }
   }

   return res;
}


SxDiracVec<Double> SxWavesCmp::absSqr (int i, int iSpin, int ik) const
{
   return waves1(i,iSpin,ik).absSqr();
}


SxMatrix<Double> SxWavesCmp::getOverlap (const SxArray<int> &states,
                                         int iSpin, int ik) const
{
   int i, j, nStates = int(states.getSize());
   SxMatrix<Double> S (nStates, nStates);

   PsiGI psi1 = waves1(iSpin,ik);
   PsiGI psi2 = waves2(iSpin,ik);
   for (i=0; i < nStates; ++i)  {
      for (j=0; j < nStates; ++j)  {
         S(i,j) = (psi1.colRef(states(i)) ^ psi2.colRef(states(j))).chop()
                  .absSqr ();
      }
   }

   return S;
}

#else /* SX_STANDALONE */


int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   SxString waves1     = cli.option ("-w1|--waves1", "file",
                                     "1st SPHInX input file")
                                    .toString ("waves.sxb");
   SxString waves2     = cli.option ("-w2|--waves2", "file",
                                     "2nd SPHInX input file").toString ();
   SxString outFile    = cli.option ("-o|--output", "file",
                                     "Out wavefunction file.")
                                     .toString("waves.diff.sxb");
   SxList<int> states  = cli.option ("-n|--nStates", "indices", 
                                     "index list of Bloch states")
                                     .required (false).toIdxList ();
   SxList<int> spins   = cli.option ("-s|--spin", "indices", 
                                     "index list of spin channels "
                                     "(1-up, 2-down)").required(false)
                                     .toIdxList();
   SxList<int> kPoints = cli.option ("-k|--kPoints", "indices", 
                                     "index list of Bloch states")
                                     .required (false).toIdxList ();
   bool calcS = cli.option ("-S|--overlap",
                            "Compute overlap matrix/matrices |S(i,j)(s,k)|^2")
                            .toBool();
   bool calcD = cli.option ("-d|--difference",
                            "Compute waves1(i,s,k) - waves2(i,s,k)")
                            .toBool();
   bool calc2 = cli.option ("-2|--absSqr",
                            "Compute |waves(i,iSpin,ik)|^2)")
                            .toBool();

   cli.finalize ();

   if (!calcS && !calcD && !calc2)  {
      sxprintf ("Error: No command specified.\n");
      sxprintf ("       Commands are:  -S | --overlap\n");
      sxprintf ("                      -d | --difference\n");
      sxprintf ("                      -2 | --absSqr\n");
      return 1;
   }
   

   initSPHInXMath ();

   
   int i, nStates, iSpin, nSpin, ik, nk;

   try  {
      SxBinIO io (waves1, SxBinIO::BINARY_READ_ONLY);
      nSpin   = io.getDimension ("nSpin");
      nk      = io.getDimension ("nk");
      nStates = io.getDimension ("nAllStates") / nk;

      if (states.getSize() == 0) 
         for (i=0; i < nStates; ++i)  states << i;

      if (spins.getSize() == 0) 
         for (i=0; i < nSpin; ++i)  spins << i;

      if (kPoints.getSize() == 0)
         for (i=0; i < nk; ++i)     kPoints << i;
      
      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   cout << "states\n"  << states   << endl;
   cout << "spin\n"    << spins    << endl;
   cout << "kpoints\n" << kPoints  << endl;

   SxWavesCmp wavesCmp;
   wavesCmp.readWaves1 (waves1);

   if (calc2)  {
      SxBinIO io;
      for (ik=0; ik < kPoints.getSize(); ++ik)  {
         for (iSpin=0; iSpin < spins.getSize(); ++iSpin)  {
            for (i=0; i < states.getSize(); ++i)  {
               io.open (outFile + "-" + (states(i)+1)
                                + "-" + (spins(iSpin)+1)
                                + "-" + (kPoints(ik)+1)
                                + ".dat",
                        SxBinIO::ASCII_WRITE_ONLY);
               io.writeXYPlot (
                     toVector(
                        wavesCmp.absSqr(states(i), spins(iSpin), kPoints(ik))
                     )
               );
               io.close ();
            }
         }
      }
      return 0;
   }


   wavesCmp.readWaves2 (waves2);

   if (calcS)  {
      for (ik=0; ik < kPoints.getSize(); ++ik)  {
         for (iSpin=0; iSpin < spins.getSize(); ++iSpin)  {
            cout << SX_SEPARATOR;
            cout << "| ik="    << kPoints(ik)+1 
                 << ", iSpin=" << spins(iSpin)+1 << endl;
            cout << SX_SEPARATOR;
            cout << wavesCmp.getOverlap (states, iSpin, ik) << endl;
         }
      }
   } 

   if (calcD)  {
      SxPW wavesDiff (wavesCmp.getDifferences ());
      SxGkBasis         Gk;
      SxFermi           fermi;
      SxAtomicStructure structure;
      SxString cmdLine;
      for (i=0; i < argc; ++i) cmdLine += SxString(argv[i]) + " ";
      try  {
         // --- clone data of |G+k>, fermi, and structure from waves1-file
         SxBinIO io (waves1, SxBinIO::BINARY_READ_ONLY);
         Gk.read (io);
         fermi.read (io);
         structure.read (io);
         io.close();
      
         // --- save new difference file
         io.open (outFile, SxBinIO::BINARY_WRITE_ONLY);
         io.write ("input", cmdLine);
         fermi.write (io);
         Gk.write (io);
         wavesDiff.write (io);  // largest array must be called in the end!
         structure.write (io);
         io.setMode (SxBinIO::WRITE_DATA);
         io.write ("input", cmdLine);
         fermi.write (io);
         Gk.write (io);
         wavesDiff.write (io);
         structure.write (io);
         io.close ();
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
   }

   return 0;

}


#endif /* SX_STANDALONE */
