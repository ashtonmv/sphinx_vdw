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

#include <SxConfig.h>

#include <SxCLI.h>
#include <SxParser.h>
#include <SxHamInit.h>

void printSize (double size)
{
   cout << size;
   SxList<SxString> units;
   units << "B" << "KB" << "MB" << "GB" << "TB";
   int iUnitMax = int(units.getSize ()) - 1;

   int iUnit = 0;
   double count = size;
   while (count > 1200. && iUnit < iUnitMax)  {
      count /= 1024.;
      iUnit++;
   }
   if (iUnit == 0) return;

   int digits;
   if (count >= 100.) digits = 0;
   else if (count >= 10.) digits = 1;
   else digits = 2;
   
   sxprintf(" (%.*f %s)", digits, count, units(iUnit).ascii ());
}

void printSize (ssize_t size)
{
   printSize (double(size));
}

void printSize (size_t size)
{
   printSize (double(size));
}


int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.authors = "C. Freysoldt";
   cli.preUsageMessage =
      "This add-on estimates the data sizes from a S/PHI/nX calculation.";
   SxString inFile =
      cli.option ("-i|--input","file","S/PHI/nX input file (for nl potentials)")
      .toString ("input.sx");

   cli.finalize ();

   // --- read input
   SxParser parser;
   SxParser::Table table = parser.read (inFile);

   int nPulay = -1, blockSize = -1;
	try {

      SxSymbolTable *main = table->getGroup ("main");

      SxSymbolTable *scfDiag = main->getGroup ("scfDiag",true);

      nPulay = scfDiag->get("nPulaySteps")->toInt ();

      if (scfDiag->containsGroup ("blockCCG")) {
         
         SxSymbolTable *blockCCG = scfDiag->getGroup ("blockCCG");
         blockSize = blockCCG->get("blockSize")->toInt ();

       } else {
          blockSize = 1;
       }

   } catch (SxException e) {

       nPulay = 1;
       blockSize = 1;

   }

   SxAtomicStructure structure(&*table);
   (cout << "structure done\n").flush ();
   SxKPoints kPoints(structure.cell, &*table);
   (cout << "k done\n").flush ();
   SxSpeciesData pPot = *SxHamInit::setupPotential (&*table);
   (cout << "species data done\n").flush ();

   double eCut = SxGBasis::getECut (&*table);
   SxVector3<Int> mesh = SxGBasis::getMesh (&*table);
   
   // k-points, FFT mesh
   cout << endl << SX_SEPARATOR;
   cout << "number of k-points : " << kPoints.nk << endl;
   cout << "plane-wave cut-off : " << eCut << endl;
   cout << "FFT mesh           : ";
   sxprintf ("%dx%dx%d\n", mesh(0), mesh(1), mesh(2));
   cout << "FFT meshsize       : " << mesh.product () << endl;
   cout << SX_SEPARATOR;

   // G-vectors
   double lCut = sqrt(eCut);
   double cutVol = FOUR_PI / 3. * lCut * lCut * lCut;
   int ng = int (lround(cutVol / structure.cell.getReciprocalCell ().volume));
   
   cout << "Number of G-vectors per k-point: ~" << ng << endl;

   // states
   int nSpin         = SxHamiltonian::getNSpin (&*table);
   double nElectrons = SxHamiltonian::getNExcessElectrons (&*table);
   for (int is = 0; is < pPot.getNSpecies (); ++is)
      nElectrons += pPot.valenceCharge(is) * structure.getNAtoms (is);
   int nStates = SxHamiltonian::getNStates (&*table, nElectrons);

   cout << "Number of spin channels : " << nSpin << endl;
   cout << "Number of electrons     : " << nElectrons << endl;
   cout << "Number of states        : " << nStates << endl;

   int nAllStates = nSpin * nStates * kPoints.nk;

   double waveCoeffSize = double(ng) * nAllStates * sizeof(PrecCoeffG);
   double n123Size = double(kPoints.nk) * double(ng) * sizeof(int);
   double gkVecSize = double(3 * kPoints.nk) * double(ng) * sizeof(PrecG);

   double tauSize = double(structure.getNAtoms () * 3) * sizeof(PrecTauR);
// double inputSize = SxParser_buffer.getSize ();
   ssize_t epsFoccSize = nAllStates * (sizeof (PrecEps) + sizeof(PrecFocc));

   ssize_t other = 2 * kPoints.nk * sizeof(int)           // nPerK, nGk
               + 3 * sizeof(int)                        // meshDim
               + nAllStates * sizeof(int)               // epsSortIdx
               + 4 * kPoints.nk * sizeof(double)        // kWeights, kVec
               + structure.getNSpecies () * sizeof(int) // nAtoms
               + 10 * sizeof(double)                    // cell, eCut
               + 1200                                   // ~ netCDF header
                                                        // chemNames
              ;

   double sum = waveCoeffSize + n123Size + gkVecSize + tauSize
                + /*inputSize +*/ double(epsFoccSize + other);

   cout << SX_SEPARATOR;
   cout << "Estimated 'waves.sxb' size: " << endl;
   cout <<   "Wave function coefficients: "; printSize (waveCoeffSize);
   cout << "\neps/focc                  : "; printSize (epsFoccSize);
   cout << "\nFFT mesh mapping          : "; printSize (n123Size);
   cout << "\nG+k vectors               : "; printSize (gkVecSize);
   cout << "\natomic structure          : "; printSize (tauSize);
// cout << "\ninput file                : "; printSize (inputSize);
   cout << "\nother                     : "; printSize (other);
   cout << "\n-------------------------------------------------------";
   cout << "\nsum                       : "; printSize (sum);

   cout << endl << SX_SEPARATOR;

   cout << "Estimated 'rho.sxb' size: ";
   printSize (double(mesh.product ())*sizeof(double)*double(nSpin) // mesh
              + 9 * (int)sizeof(double)                            // cell
              + 3 * (int)sizeof(int)                               // meshDim
              + 200                                                // header
             );
   cout << endl;
   

   cout << SX_SEPARATOR;
   double fftMeshReal    = double(mesh.product ()) * sizeof(double);
   ssize_t fftMeshComplex = mesh.product () * sizeof(SxComplex16);

   ssize_t phase1D = structure.getNAtoms () * mesh.sum ()
                   * sizeof(SxComplex16); // 1D phase factors
   ssize_t structFact = structure.getNSpecies () * ng * sizeof(SxComplex16);
   double phaseFact  = double (kPoints.nk * structure.getNAtoms () * ng)
                       * sizeof(SxComplex16)
                     + double(kPoints.nk * 2 * fftMeshComplex); // FFT meshes
   double phaseFacSmall = double(kPoints.nk 
                          * (phase1D              // 1D phase factors
                             + ng * sizeof(ssize_t) * 3) // packedGRel
                        + fftMeshComplex);        // on-the-fly FFT mesh

   double gkBase = double(kPoints.nk * (structFact + ng*sizeof(double)/*g2*/))
                  + n123Size + gkVecSize;
   ssize_t rgSize = 2*fftMeshComplex         // R FFT mesh
                  + 8*structFact             // big G
                  + 8*ng*sizeof(double)*4    // gVec + g2 in big G
                  + phase1D                  // 1D phase factors
                  + 8*ng*sizeof(ssize_t) * 3 // packedGRel
                  + 8*ng*sizeof(int);        // n123 in big G

   cout << "| Major main memory items:     " << endl;
   cout << "| G+k basis:                   ";
   printSize (gkBase +  phaseFact); cout << endl;
   cout << "| G+k basis (saveMemory)       ";
   printSize (gkBase + phaseFacSmall); cout << endl;
   cout << "| R/G basis:                   ";
   printSize (rgSize); cout << endl;
   cout << "| density/potentials:          ";
   printSize (fftMeshReal * (5 * nSpin // vEff vXC vEx vCor rho 
                             + 2));    // vHartreeR vLoc
   cout << endl
        << "| waves (single state)         ";
   printSize (ng * sizeof(PrecCoeffG));
   cout << endl
        << "| waves (single k-point)       ";
   printSize (nStates * ng * sizeof(PrecCoeffG));
   cout << endl;

   cout << "| mixer (per Pulay step)       ";
   printSize (fftMeshReal * 2 * nSpin);
   cout << endl;

   cout << "| blockCCG (per block state):  ";
   
   printSize (ng*sizeof(PrecCoeffG)*5); // hPsi X Xold + 2*workspace
   cout << endl;

   cout << "+---------------------------------------------------------" << endl
        << "| memory size:" << endl
        << "| max                          ";
   printSize (  gkBase + phaseFact + double(rgSize) + fftMeshReal *(5*nSpin+2)
              + kPoints.nk * nSpin * nStates * double(ng * sizeof(PrecCoeffG))
              + fftMeshReal * 2 * nSpin * nPulay 
              + nSpin * double(ng*sizeof(PrecCoeffG))*5 * blockSize);

   cout << endl
        << "| min (saveMem,keepWaves)      ";
   printSize (  gkBase + phaseFacSmall + double(rgSize)
              + fftMeshReal *(5*nSpin+2)
              + nStates * double(ng * sizeof(PrecCoeffG))
              + fftMeshReal * 2 * nSpin * nPulay
              + double(ng*sizeof(PrecCoeffG))*5 * blockSize);
   cout << endl;

   /*
   cout << "| lcao                     ";
   printSize (double(nAO)*nAO*sizeof(PrecCoeffG) // matrix 
              * (4.5 + 2 * kPoints.nk)); // H S Sinv eig.vecs Hsym + Hkstat + Lk
   cout << endl;
   */

   cout << SX_SEPARATOR;

   
}

