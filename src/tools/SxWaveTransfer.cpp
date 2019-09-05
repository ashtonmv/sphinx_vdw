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
#include <SxPW.h>
#include <SxFermi.h>
#include <SxPWOverlap.h>
#include <SxPAWOverlap.h>

// --- class for writing waves with changed mesh directly from
//     a netcdf binary
class SxChangeCutoffPW : public SxPW  {
   protected:
      SxString inFile;
      const SxGkBasis &oldGk;
      const SxGkBasis &newGk;
      bool renormalize;
   public:
      double formatVersion;
      SxChangeCutoffPW (const  SxString &filename,
                        SxPtr<SxGkBasis> oldBasis,
                        SxPtr<SxGkBasis> newBasis,
                        bool   renorm,
                        double format)
         : inFile(filename), 
           oldGk(*oldBasis),
           newGk(*newBasis),
           renormalize(renorm),
           formatVersion(format)
          { 
             nStatesPerK.resize (oldGk.getNk ());
             nkPoints = oldGk.getNk ();
             gkBasisPtr = newBasis;
             try {
                SxBinIO io(filename, SxBinIO::BINARY_READ_ONLY);
                io.read ("nPerK",&nStatesPerK,nkPoints);
                nSpin = io.getDimension("nSpin");
             } catch (SxException e)  {
                e.print ();
                SX_EXIT;
             }
             waves.resize(1);
             waves(0).resize(nSpin);
          }
      virtual void write (SxBinIO &) const;
};

// --- this is basically copied from SxPW::write
//     the current format is 1.000 - 1.001
void SxChangeCutoffPW::write (SxBinIO &io) const
{
   bool writeKwise = formatVersion > 1.000;
   bool readKwise;
   SxBinIO in;
   try  {
      in.open (inFile, SxBinIO::BINARY_READ_ONLY);
      readKwise = ! in.containsDim("nCoeff");
      if (readKwise && io.ncMode == SxBinIO::WRITE_DATA)
         cout << "Reading waves k-wise" << endl;
      int iSpin;
      int ik,    nk    = oldGk.getNk ();
      int nCoeff = 0, i, ng = -1;

      SxVector<Int> nPerK(nk), nGk(nk);
      SxVector<Int> remap, mesh(newGk(0).fft3d(0).meshSize);

      in.read ("nPerK", &nPerK, nk);
      nSpin = in.getDimension ("nSpin");

      for (ik=0; ik < nk; ik++)  {
         nGk(ik)   = ng = newGk(ik).ng;
         nCoeff += nPerK(ik) * ng * nSpin;
      }
      int offset = 0, inOffset = 0;

      // write version
      io.write ("formatVersion", formatVersion);

      // --- create dimensions
      io.addDimension ("nk", nk);
      if (writeKwise) {
         for (ik = 0; ik < nk; ++ik)
            io.addDimension ("nCoeff-"+SxString(ik+1), 
                             nPerK(ik) * newGk(ik).ng * nSpin);
      } else {
         io.addDimension ("nCoeff", nCoeff);
      }
      io.addDimension ("nSpin", nSpin);

      // --- write data
      if (!io.contains ("nGk") || (io.ncMode == SxBinIO::WRITE_DATA))
         io.write ("nGk",   nGk,   "nk");
      if (!io.contains ("nPerK") || (io.ncMode == SxBinIO::WRITE_DATA))
         io.write ("nPerK", nPerK, "nk");
      PsiG psi, psiOld;
      SxString psiVarNameWrite   = "psi",
               coeffDimNameWrite = "nCoeff",
               psiVarNameRead    = "psi";
      
      for (ik=0; ik < nk; ik++)  {
         // --- prepare map
         if (io.ncMode == SxBinIO::WRITE_DATA)  {
            mesh = -1;
            ng = oldGk(ik).ng;
            SxVector<TPrecFFTIdx>::Iterator n123It = oldGk(ik).n123(0).begin ();
            SxVector<Int>::Iterator mapIt;
            SxMesh3D oldFft = oldGk(ik).fft3d(0).mesh,
                     newFft = newGk(ik).fft3d(0).mesh;
            // write old indices to new mesh
            SxVector3<Int> G;
            for (int ig = 0; ig < ng; ++ig, ++n123It)  {
               G = oldFft.getMeshVec(*n123It,SxMesh3D::Origin);
               mesh(newFft.getMeshIdx(G, SxMesh3D::Origin)) = ig;
            }
            // now mesh(iMesh) contains the index in the old basis

            // --- collect indices
            ng = newGk(ik).ng;
            remap.resize (ng);
            mapIt = remap.begin ();
            n123It = newGk(ik).n123(0).begin ();
            for (int ig = 0; ig < ng; ++ig, ++n123It, ++mapIt)
               *mapIt = mesh(*n123It);
            // --- prepare psi
            psi.resize (ng);
            psiOld.resize(oldGk(ik).ng);
            psi.set (0.);
         }

         if (writeKwise)  {
            offset = 0;
            SxString sK(ik+1);
            psiVarNameWrite = "psi-" + sK;
            coeffDimNameWrite = "nCoeff-" + sK;
         }
         if (readKwise)  {
            inOffset = 0;
            psiVarNameRead = "psi-" + SxString(ik+1);
         }
         for (iSpin=0; iSpin < nSpin; iSpin++)  {
            for (i=0; i < nPerK (ik); i++)  {
               if (io.ncMode == SxBinIO::WRITE_DATA)  {
                  // read psi in old basis
                  in.read (psiVarNameRead, &psiOld, oldGk(ik).ng, inOffset);
                  // --- remap it
                  SxVector<Int>::Iterator mapIt = remap.begin ();
                  PsiG::Iterator psiIt = psi.begin ();
                  SX_CHECK (ng == nGk(ik), ng, nGk(ik));
                  for (int ig = 0; ig < ng; ++ig, ++psiIt, ++mapIt)
                    if (*mapIt != -1) *psiIt = psiOld(*mapIt);
                  // normalize if requested
                  if (renormalize) psi.normalize ();
               }
               // write psi
               io.write (psiVarNameWrite, psi, coeffDimNameWrite, offset);
               SX_CHECK (   io.ncMode != SxBinIO::WRITE_DATA
                         || psi.getSize () == nGk(ik),
                         psi.getSize (), nGk(ik));
               offset += nGk(ik);
               inOffset += oldGk(ik).ng;
            }
         }
         // at the end, all coefficients must be written
         SX_CHECK ( (!writeKwise) || offset == nPerK(ik) * nGk(ik) * nSpin,
                   offset, nPerK(ik) * nGk(ik) * nSpin);
      }
      // at the end, all coefficients must be written
      SX_CHECK (writeKwise || offset == nCoeff, offset, nCoeff);

      // close input file
      in.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}

int main (int argc, char** argv)
{
   // --- parse command line options
   SxCLI cli (argc, argv);
   cli.preUsageMessage = SxString(
         "This add-on manipulates wave function files. It can compare the "
         "stored basis to the one generated by the system (check), and it can "
         "reorder the stored waves according to the system-generated basis "
         "(reorder)- which is the only way to transfer waves files between "
         "different systems for older S/PHI/nX-versions.\n"
         "It can also rewrite the waves to a changed cut-off energy (change) "
         " or reduce the number of states (resize)."
         ).wrap ();

   SxString inputFile 
      = cli.option ("--S/PHI/nXinput", "file", "the S/PHI/nX input file")
        .toString ("input.sx");
   SxString wavesIn 
      = cli.option ("-i", "waves file",
                    "waves function to be read in")
        .toString ("waves.sxb");
   SxString wavesOut
      = cli.option ("-o", "waves file",
                    "waves function to be written out")
        .toString ("waves-out.sxb");

   double formatVersion 
      = cli.option("--format", "version", "request specific format version (change,reorder)")
        .toDouble(1.001,1.000-1e-9,1.001+1e-9);

   if (wavesIn == wavesOut)  {
      cout << "Input and output file must differ!" << endl;
      cli.setError ();
   }
   int energyGroup = cli.newGroup ("change cutoff energy");
   double newCutoff
      = cli.option ("-e|--eCut", "energy", "the new cut-off energy "
                    "(only for command 'change')")
        .toDouble ();
   int statesGroup = cli.newGroup ("change cutoff energy");
   int nStatesNew = cli.option ("-n|--nStates", "states", 
                    "the new number of states (only for command resize)")
                    .toInt (false,0);
   cli.setGroup(cli.generalGroup);
   SxString commandList ("check reorder change resize select propagate");
   SxString command
      = cli.argument ("command","one of: " + commandList).toString ()
        .toLower ();
   if (!cli.error && !commandList.tokenize(' ').contains(command))  {
      cout << "\nERROR: Unknown command '" << command << "'.\n";
      cli.setError ();
   }

   if (command == "resize")  {
      if (!cli.groupAvailable(statesGroup))  {
         cout << "Command 'resize': no new number of states given (-n option)!"
              << endl;
         cli.setError ();
      }
   }

   SxArray<int> states;
   cli.option ("--states", "states", "select: selected states");
   if (command == "select")  {
      if (cli.last ().exists ())  {
         states = cli.last ().toIdxList ();
      } else {
         cout << "Missing list of select states (--states)" << endl;
         cli.setError ();
      }
   }

   if (command == "change")  {
      if (!cli.groupAvailable(energyGroup))  {
         cout << "Command 'change': no new cutoff energy given (-e option)!"
              << endl;
         cli.setError ();
      }
   }

   cli.finalize ();

   initSPHInXMath ();
   SxLoopMPI::init (argc, argv);
   
   SxAtomicStructure structure;
   double eCut;
   SxMesh3D mesh;
   
   SxPtr<SxGkBasis> gkOrigPtr = SxPtr<SxGkBasis>::create ();
   SxGkBasis &gk = *gkOrigPtr;
   try  {
      SxBinIO io(wavesIn, SxBinIO::BINARY_READ_ONLY);
      structure.read (io);
      gk.read (io);
      io.read ("eCut", &eCut);
      io.read ("meshDim", &mesh);
      io.close ();
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
   int nk = gk.getNk ();

   SxPtr<SxGkBasis> gkPtr;
   SxPtr<SxGBasis> gPtr;
   bool renormalize = false;

   // --- Read and check G+k basis from input file
   if (command == "check" || command == "reorder")  {
      SxParser parser;
      SxParser::Table table = parser.read (inputFile);
      double eCutIn = SxGBasis::getECut(&*table);
      if (fabs(eCutIn - eCut) > 1e-6)  {
         cout << "Mismatch in cutoff energies" << endl;
         cout << "Use sxwavetransfer 'change' to change the cutoff energy." 
              << endl;
         SX_QUIT;
      }
      gPtr = SxPtr<SxGBasis>::create (SxGBasis::getMesh(&*table),
                                              structure, 
                                              4. * eCutIn);
      gkPtr = SxPtr<SxGkBasis>::create (*gPtr, &*table);
      if (gkPtr->getNk () != nk)  {
         cout << "Mismatch in number of k-points." << endl;
         cout << inputFile << ": " << gkPtr->getNk () << endl;
         cout << wavesIn << ": " << nk << endl;
         cout << 
            "Possible causes\n"
            "* input file doesn't match waves file\n"
            "* symmetries have changed\n";
         SX_QUIT;
      }
      
      SxGBasis *inG, *waveG;
      for (int ik = 0; ik < nk; ++ik)  {
         if ((gkPtr->getK(ik) - gk.getK(ik)).normSqr () > 1e-6)  {
            cout << "Mismatch in k-point no. " << (ik+1) << endl;
            cout << inputFile << ": " << gkPtr->getK(ik) << endl;
            cout << wavesIn << ": " << gk.getK(ik) << endl;
            cout << 
               "Possible causes\n"
               "* input file doesn't match waves file\n"
               "* symmetries have changed\n";
            SX_QUIT;
         }
         inG = &(*gkPtr)(ik);
         waveG = &gk(ik);
         if (inG->ng != waveG->ng)  {
            cout << "Mismatch in number of G-vectors for k-point no. " 
                 << (ik+1) << endl;
            cout << inputFile << ": " << inG->ng << endl;
            cout << wavesIn   << ": " << waveG->ng << endl;
            cout << 
               "Possible causes\n"
               "* input file doesn't match waves file\n"
                 << endl;
            cout << "command 'reorder' will renormalize waves" << endl;
            if (command == "check")  {
               cout << "|G+k> bases differ in number of G vectors."  << endl
                    << "You must run 'reorder' to get this fixed." << endl;
               return 1;
            }
            renormalize = true;
         } else if (command == "check") {
            SxVector<TPrecFFTIdx>::Iterator a = inG->n123(0).begin (),
                                            b = waveG->n123(0).begin ();
            for (int ig = 0; ig < inG->ng; ++ig, ++a, ++b)
               if (*a != *b) {
                  cout << "|G+k> bases differ."  << endl
                       << "You may run 'reorder' to get this fixed." << endl;
                  return 1;
               }
         }
         
      }
      if (command == "check")  {
         cout << "The G-vector ordering in the waves file is correct." << endl;
         return 0;
      }
   }

   // --- create new G+k basis with desired cutoff energy
   if (command == "change")  {
      mesh = SxGBasis::getMeshSize(newCutoff, structure.cell);
      gPtr = SxPtr<SxGBasis>::create (mesh, structure, 4. * newCutoff);
      gkPtr = SxPtr<SxGkBasis>::create (gk, *gPtr, newCutoff, true);
      renormalize = (newCutoff <= eCut);
   }


   // --- read Fermi and chemical symbols
   SxFermi fermi;
   try  {
      SxBinIO io(wavesIn, SxBinIO::BINARY_READ_ONLY);
      fermi.read (io);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   if (command == "reorder" || command == "change")  {
      SxChangeCutoffPW reorder(wavesIn, gkOrigPtr, gkPtr, renormalize,formatVersion);
      reorder.writeWavesFile (wavesOut, fermi, structure, false);
   }

   // --- reduce number of states
   if (command == "resize")  {
      // open for reading wave function coefficients k-point-wise
      SxPW waves(wavesIn, SxPW::ReadOnDemand);

      // --- check number of states
      for (int ik = 0; ik < nk; ++ik)
         if (nStatesNew > waves.getNStates (ik))  {
            cout << "Increasing the number of states is not possible" << endl;
            cout << "new nStates = " << nStatesNew << endl;
            cout << "nStates for k-point no. " << (ik+1) << ": " 
                 << waves.getNStates(ik) << endl;
            SX_QUIT;
         }
      // --- resize
      waves.changeNStates (nStatesNew);
      fermi.resize (nStatesNew);
      // output
      waves.writeWavesFile (wavesOut, fermi, structure, false);

   }

   if (command == "select")  {
      gPtr = SxPtr<SxGBasis>::create (mesh, structure, 4. * eCut);
      gk.gBasis = gPtr.getPtr ();
      // open for reading wave function coefficients k-point-wise
      SxPW waves(wavesIn, SxPW::ReadOneByOne);

      nStatesNew = int(states.getSize ());
      SxPW wavesNew(nStatesNew, waves.getNSpin (), gkPtr);

      for (int ik = 0; ik < nk; ++ik)  {
         for (int iSpin = 0; iSpin < waves.getNSpin (); ++iSpin)  {
            for (int iNew = 0; iNew < states.getSize (); ++iNew)  {
               wavesNew(iNew,iSpin,ik) <<= waves(states(iNew), iSpin, ik);
               fermi.eps(iNew,iSpin,ik) = fermi.eps(states(iNew), iSpin, ik);
               fermi.focc(iNew,iSpin,ik) = fermi.focc(states(iNew), iSpin, ik);
            }
         }
      }
      fermi.resize (nStatesNew);
      wavesNew.writeWavesFile (wavesOut, fermi, structure, false);
   }

   if (command == "propagate")  {
      SxPW waves(wavesIn, SxPW::ReadOnDemand);
      waves.setGkBasisPtr(gkOrigPtr);
      SxParser parser;
      SxParser::Table table = parser.read (inputFile);
      double eCutIn = SxGBasis::getECut(&*table);
      if (fabs(eCutIn - eCut) > 1e-6)  {
         cout << "Mismatch in cutoff energies" << endl;
         cout << "Use sxwavetransfer 'change' to change the cutoff energy." 
              << endl;
         SX_QUIT;
      }
      SxAtomicStructure newStruct(&*table);
      
      gPtr = SxPtr<SxGBasis>::create (SxGBasis::getMesh(&*table),
                                              newStruct, 
                                              4. * eCutIn);
      gkPtr = SxPtr<SxGkBasis>::create (*gPtr, &*table);
      if (gkPtr->getNk () != nk)  {
         cout << "Mismatch in number of k-points." << endl;
         cout << inputFile << ": " << gkPtr->getNk () << endl;
         cout << wavesIn << ": " << nk << endl;
         cout << 
            "Possible causes\n"
            "* input file doesn't match waves file\n"
            "* symmetries have changed\n";
         SX_QUIT;
      }
      SxPW wavesNew(waves.getNStates(), waves.getNSpin (), gkPtr);

      SxPtr<SxOverlapBase> SPtr;
      if (table->containsGroup("pseudoPot")) 
         SPtr = SxPtr<SxPWOverlap>::create ();
      else if (table->containsGroup("pawPot"))   {
         SxPtr<SxPAWPot> pawPotPtr = SxPtr<SxPAWPot>::create (&*table);
         SxPtr<SxPartialWaveBasis> pBasis
            = SxPtr<SxPartialWaveBasis>::create (pawPotPtr, structure);
         pBasis->createProjBasis (*gkPtr);
         SPtr = SxPtr<SxPAWOverlap>::create (pBasis, pawPotPtr);
      } else   {
         cout << "No known Potential Group found!" << endl;
         SX_QUIT;
      }

      for (int ik = 0; ik < waves.getNk (); ik++)  {
         SxCell oldRecCell = 
            structure.cell.getReciprocalCell();
         SxCell newRecCell = 
            newStruct.cell.getReciprocalCell();
         SxVector3<Double> oldKRel = 
            oldRecCell.carToRel(gkOrigPtr->kVec(ik));
         SxVector3<Double> newKRel = 
            newRecCell.carToRel(gkPtr->kVec(ik));
         if ((newKRel - oldKRel).norm () > 1e-6)  {
            cout << "Consistency check failed for k-point " << (ik+1);
               cout << ".\n Expected " << oldKRel << " but found ";
               cout << newKRel << "." << endl;
               SX_QUIT;
         }
         for (int iSpin = 0; iSpin < waves.getNSpin (); iSpin++)  {
            for (int iState = 0; iState < waves.getNStates (ik); iState++)  {
               PsiG psiOld = waves(iState,iSpin,ik);
               PsiG psiNew = wavesNew(iState,iSpin,ik);
               psiNew <<= (*gkPtr)(ik).transferWaves (psiOld);
            }
            SPtr->orthonormalize(&wavesNew(iSpin,ik));
         }
      }
      wavesNew.writeWavesFile (wavesOut, fermi, newStruct, false);
   }
   return 0;

}
