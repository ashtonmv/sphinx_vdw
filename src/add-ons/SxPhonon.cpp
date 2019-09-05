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

#include <SxPhonon.h>
#include <SxParser.h>
#include <iostream> 
#include <fstream>
#include <SxCLI.h>
#include <SxDynMat.h>
#include <SxRotation.h>

#ifndef SX_STANDALONE

SxArray<SxMatrix<Complex16>::Eigensystem> SxPhonon::readInput (
      const SxString &setInput,
      const SxString &dynMatInput,
      const double &epsE)
{  
   cout << "readInput" << endl;

   //check that inputs provided
   SX_CHECK (dynMatInput != "" && setInput != "");
   
   SxArray<SxMatrix<Complex16>::Eigensystem> inEigs;
   try  {
      // setInput to parser
      SxParser parser;
      SxParser::Table parTab = parser.read (setInput);
      
      // ---read optional qMesh
      SxSymbolTable *qMeshPtr = NULL; 
      if (parTab->containsGroup ("qMesh"))  {
         qMeshPtr = parTab->getGroup ("qMesh");  
         cout << "@qMesh read" << endl;
      }  else  {
         cout << "@qMesh not read" << endl;
      }
      
      // ---read optional qPath
      SxSymbolTable *qPathPtr = NULL; 
      if (parTab->containsGroup ("qPath"))  {
         qPathPtr = parTab->getGroup ("qPath");
         cout << "@qPath read" << endl;
      }  else  {
         cout << "@qPath not read" << endl;
      }
      
      // ---read optional indices for special q-point / mode
      if (parTab->containsGroup ("eigenVecs"))  {
         SxSymbolTable *evPtr = parTab->getGroup ("eigenVecs");
         sQ = evPtr->get ("sQ")->toInt ();
         sMod = evPtr->get ("sMode")->toInt ();
      }
      cout << "sQ, sMod read" << sQ << " " << sMod << endl;

      // ---read sxb file
      SxBinIO io (dynMatInput, SxBinIO::BINARY_READ_ONLY);

      // -read primitive cell
      SxString inVar ("primCell");
      SxMatrix3<Double> cellMat;
      io.read (inVar, &cellMat);//typecast to SxMatrix first
      primStr = SxAtomicStructure (SxCell (cellMat));
      primStr.epsEqual = epsE;
      cout << "primCell read" << primStr.cell << endl;
      
      // -read primitive positions
      int nAtom = io.getDimension ("nAtom");
      SxMatrix<Double> primPos (3, nAtom);
      io.read ("primPos", &primPos, 3, nAtom);
      cout << "primPos read" << endl;
      primPos.print (true);

      // -transform them to a structure
      primStr.resize (nAtom);
      for (int iAtom = 0; iAtom < nAtom; iAtom++)  
      { primStr.setAtom (iAtom, primPos.col (iAtom).toVector3 ()); }
      cout << "primStr set" << endl;

      // -read exact qPoints
      int nQ = io.getDimension ("nQ");
      SxMatrix<Double> xQp (3, nQ);
      io.read ("exactQ", &xQp, 3, nQ);
      cout << "exactQ read" << endl;
      xQp.print (true);

      // -transform them to an SxArray<Coord>
      exactQ.resize (nQ);
      for (int iQ = 0; iQ < nQ; iQ++)  
      { exactQ(iQ) = Coord (xQp.col (iQ).toArray ()); }

      // set dimensions of eigensystem
      int nMode = 3 * nAtom;
      
      // -read eigenvals
      SxMatrix<Double> eigVals (nMode, nQ);
      io.read ("eigVals", &eigVals, nMode, nQ);
      cout << "eigVals read " << endl;
      eigVals.print (true);

      // -read eigenvecs
      SxVector<Complex16> eigVecs (sqr (nMode) * nQ);
      io.read ("eigVecs", &eigVecs, sqr (nMode) * nQ);
      cout << "eigVecs read, size " << eigVecs.getSize () << endl;

      // --set eigensystem
      inEigs.resize (nQ);
      eigVecs.reshape (sqr (nMode), nQ);
      SxMatrix<Complex16> vecs; 
      for (int iQ = 0; iQ < nQ; iQ++)  { 
         inEigs(iQ).vals = eigVals.col (iQ);       
         vecs = eigVecs.col (iQ);
         vecs.reshape (nMode, nMode); //columnwise
         inEigs(iQ).vecs = vecs; //cannot direct?!
      }  
         
      // generate qMesh (or qPath when mesh not provided) from primCell
      if (qMeshPtr != NULL)  {
         qPoints = SxKPoints (primStr.cell, qMeshPtr);
         cout << "qMesh initialised\n";
      }  else  {
         
         // at least one mesh is required
         SX_CHECK (qPathPtr != NULL);
         qPoints = SxKPoints (primStr.cell, qPathPtr);
         cout << "qPath initialised\n";
      }
      // no. qPoints to standard out
      cout << "@nq " << qPoints.nk << endl;
   }  catch (SxException e)  {
      e.print();
      SX_QUIT;
   } 
   return inEigs;
}

SxArray<Coord> SxPhonon::getLatShifts (ofstream &output) const {  

   cout << "getLatShift" << endl;

   // ---make list of Coord
   int nCell = int(exactQ.getSize ());
   
   // use reciprocal of the primitive cell to get the right shifts
   SxCell recPrimCell = primStr.cell.getReciprocalCell ();
   
   // add basis vectors in case of no lattice shifts
   SxList<Coord> qL;
   qL << recPrimCell.basis (0) << recPrimCell.basis (1) 
      << recPrimCell.basis (2);
  
   // map into cell for set from vectors  
   for (int iQ = 0; iQ < nCell; iQ++)  
   { qL << recPrimCell.getMapped (exactQ(iQ), SxCell::Origin); }
   
   // set reciprocal Cell
   SxCell recCell;
   recCell.setFromVectors (qL);
   recCell = recCell.getRegularCell (); //=compact cell
   
   //---remove noise
   SxMatrix3<Int> relS;
   for (int iDir = 0; iDir < 3; iDir++)  
   { relS.setCol (iDir, 
            round (recPrimCell.carToRel (recCell.basis (iDir)) * nCell)); }
   cout << "relS" << relS << " / " << nCell << endl;
   recCell = (recPrimCell ^ relS) / nCell;
   SX_CHECK (nCell == round (recPrimCell.volume / recCell.volume),
             nCell, recPrimCell.volume / recCell.volume);
   
   // ---use SxDynMat.getExactQ with reciprocal cells
   SxDynMat D;

   // set "supercell"
   D.str = SxAtomicStructure (recPrimCell);
   
   // provide "primCell"
   output << "#get latShifts as 'exactQ' of recCell:" << endl << "#";
   SxArray<Coord> res = D.getExactQ (recCell, output);
   cout << "latShifts" << res << endl;
   
   // set supercell from reciprocal of recCell
   SxCell supCell = recCell.getReciprocalCell ();
   
   //-- remove more noise (?)
   for (int iDir = 0; iDir < 3; iDir++)  {
      relS.setCol (iDir, round (primStr.cell.carToRel (supCell.basis (iDir)))); 
   }
   supCell = primStr.cell ^ relS;
   cout << "supCell" << supCell << endl;
   
   // -add the supercell basis to the lattice shifts
   res.resize (nCell + 3, true);
   res(nCell) = supCell.basis (0);
   res(nCell+1) = supCell.basis (1);
   res(nCell+2) = supCell.basis (2);
   
   //return the results
   return res;
}

SxMatrix<Complex16> SxPhonon::fourierReduce (
      const int iQ, 
      const Coord &latShift,
      const SxCell &supCell,
      const bool debug) const 
{  
   if (debug)  { cout << "fourierRed" << endl; }
   
   // -set correct sizes
   bool write = false;
   int nAtom = primStr.getNAtoms (),
       iAtom, iLat, iDir, jDir; //index, 
   double delta; //summed distances of atom to borders 
   Coord dC, diffC, //coord difference between mode and forc atom
         factor; //decrease factor of lattice vector for atom 
   SxVector3<Int> imag; //counts supCell images
   SxCell doubCell (2. * supCell);
   if (iQ == 0 && debug) cout << "doubCell" << doubCell;
   int facDen; //denominator of factor
   SxComplex16 facNum; //numerator of factor
   SxMatrix<Complex16> res (3 * nAtom, nAtom * 3); //output
   for (int iMode = 0; iMode < nAtom; iMode++)  {
      if (debug) cout << endl << "iMode" << iMode << endl;

      // loop over primitive atoms
      for (iAtom = 0; iAtom < nAtom; iAtom++)  {
         facNum.set (0.);
         facDen = 0;

         // set coord difference in WignerSeitz cell
         diffC =  primStr (iAtom) - primStr (iMode) + latShift;
         supCell.map (&diffC, SxCell::WignerSeitz); 
         if (debug) cout << "iAtom" << iAtom << "diffC" << diffC;

         //---loop over supCell vectors
         for (iLat = 0; iLat < 8; iLat++)  {

            // set image (shift)
            imag = SxVector3<Int> (SxList<int> () 
                             <<  (iLat % 2) << ((iLat / 2) % 2) <<  (iLat / 4));
      
            // add image vector in double WignerSeitz cell
            dC = diffC + (supCell ^ imag);
            doubCell.map (&dC, SxCell::WignerSeitz);

            // sum individual dC's, remove off-set 
            // TODO what if cell not orthorhombic?
            delta = (absSqrt (dC)).sum () - (absSqrt (diffC)).sum ();
            if (delta < primStr.epsEqual)  {
               
               //atom at border: make phase of dC
              facNum += expI (qPoints.kVec (iQ) ^ dC);
              facDen++;
            }
         } 
         if (debug)  { cout << "facDen" << facDen << endl; }

         // combine Num/Den into res, this is now Hermitian
         facNum = facNum / double (facDen);
         for (iDir = 0; iDir < 3; iDir++)  {
            for (jDir = 0; jDir < 3; jDir++)  {
               res (3 * iAtom + iDir, 3 * iMode + jDir) = facNum; 
            }
         }
      }
   }
   if (write)  { cout << endl; }

   // check if res is a square matrix
   SX_CHECK (res.nCols () == res.nRows (),
             res.nCols (), res.nRows ());
   return res;
}

SxArray<SxMatrix<Complex16>::Eigensystem> SxPhonon::transformDynMat (
      const SxArray<SxMatrix<Complex16>::Eigensystem> &inEigs, 
      const SxArray<Coord> &latShifts, //includes supCell
      const bool debug) 
{  
   cout << "transformDynMat" << endl;
   int nCell = int(exactQ.getSize ()),
       nMode = int(inEigs(0).vals.getSize ()),
       iDir, iXQ, //counter over exact qPoints 
       nAtom = primStr.getNAtoms ();
   SX_CHECK (nMode = 3 * nAtom);
   SX_CHECK ((nCell == latShifts.getSize () - 3)
             && (nCell == inEigs.getSize ()));

   //---make dynMat with reverted phase factors
   SxVector<Complex16> phaseQ (nMode), //for atoms in primCell
                       phaseQs (nAtom);
   SxArray<SxMatrix<Complex16> > dmQ (nCell);//(nAtDir, nMode);
   double asym;
   cout << "Check hermiticity input dynMats" << endl;
   for (iXQ = 0; iXQ < nCell; iXQ++)  {

      //- make dynMat from eigensystem: D=Vecs*Vals*Vecs^-1 (columnwise basis)
      dmQ(iXQ) = inEigs(iXQ).vecs * inEigs(iXQ).vals; //matrix * vector 
      dmQ(iXQ) = dmQ(iXQ) ^ inEigs(iXQ).vecs.inverse ();
      
      //--check hermiticity dynMats
      asym = (dmQ(iXQ) - dmQ(iXQ).transpose ().conj ()).norm () 
           / dmQ(iXQ).norm ();
      if (asym > .02)
      {  cout << iXQ <<"iXQ Warning! dynMatIn " 
              << int (50. * asym) << "% anti-hermitian." << endl; }

      //--inverse phase factors (in primCell)
      phaseQs = expI (- exactQ(iXQ) ^ primStr); //:nAtom 
      phaseQ.reshape (nAtom, 3); //row is fast index!

      // three copies because phase is independent of direction
      for (iDir = 0; iDir < 3; iDir++)  { phaseQ.colRef (iDir) <<= phaseQs; }
      phaseQ = phaseQ.transpose (); //:(3, nAtom)
      dmQ(iXQ) = phaseQ * dmQ(iXQ) * phaseQ.conj ();
   }         //SxVector * SxMatrix * SxVector!
   cout << "dynMats revert phase constructed" << endl;

   //---dynMat(iQ)-> dynMat(iCell), uses exactQ and is supercell periodic
   SxComplex16 phaseQl; //for each qPoint and lattice vector
   SxArray<SxMatrix<Complex16> > dmC (nCell);
   double norm = 1. / nCell, //normalization factor
          fract; 
   int iCell, jCell;
   SxMatrix<Double> dmGam;
   for (iCell = 0; iCell < nCell; iCell++)  {

      // first initialize
      phaseQl = expI (- exactQ(0) ^ latShifts(iCell));
      dmC(iCell) = phaseQl * dmQ(0) * norm;

      // then add the rest
      for (iXQ = 1; iXQ < nCell; iXQ++)  {
         phaseQl = expI (- exactQ(iXQ) ^ latShifts(iCell));
         dmC(iCell) += phaseQl * dmQ(iXQ) * norm; 
      }

      //-- check reality dynMat
      fract = dmC(iCell).imag ().norm () / dmC(iCell).norm ();
      if (fract > .01)  { 
        cout << iCell << "iCell Warning! real dynMat " 
             << int (100. * fract) << "% imaginary." << endl; 
      }

      //--check drift via eigenvalues at Gamma 
      if (iCell == 0)  { dmGam = dmC(iCell).real (); } //copy
      else  { dmGam += dmC(iCell).real (); }
   }
   cout << "drift? (3 e.v.@Gam) " << dmGam.eigenvalues ()(0) 
        << dmGam.eigenvalues ()(1)
        << dmGam.eigenvalues ()(2) << endl;
   cout << "Real dynMats constructed & reality checked." << endl;
   
   //---check for symmetry real dynMats 
   SxCell supCell (latShifts(nCell), latShifts(nCell+1), latShifts(nCell+2));
   SxArray<int> shift (nCell); //inverse lattice shift numbers
   shift.set (-1);
   cout << "Check symmetry real dynMats." << endl; 
   for (iCell = 0; iCell < nCell; iCell++)  {
      for (jCell = iCell; jCell < nCell; jCell++)  {
         if (supCell.getMapped (latShifts(jCell) + latShifts(iCell),
                  SxCell::Origin).norm () < 1.e-7)  
         {
            shift(iCell) = jCell;
            fract = (dmC(iCell).real () - dmC(jCell).transpose ().real ())
                    .norm () / (dmC(iCell).real ().norm () + dmC(jCell).real ()
                    .norm ());
            if (fract > .01)
            { cout << iCell << "iCell Warning! real dynMat " 
                   << int (100. * fract) << "% anti-symmetric." << endl;}
         }
      }
   }
   cout << "iInvShift" << shift << endl;

   //--- get newDynMat ---------------------------------------------
   SxMatrix<Complex16> redPhas, invPhas, nDM; //:(3*nAtom, nMode);
   SxArray<SxMatrix<Complex16>::Eigensystem> eig (qPoints.nk);
   cout << "Start Fourier transform" << endl; 
   SxVector<Complex16>::Iterator phasIt;

   for (int iQ = 0; iQ < qPoints.nk; iQ++)  {
      nDM = SxVector<Complex16> (3 * nAtom * nMode, 0.); //initialize@0.
      nDM.reshape (3 * nAtom, nMode); //cannot in 1 step?!
      if (debug)  {cout << "\nqPoint" << qPoints.kVec (iQ);}
      for (iCell = 0; iCell < nCell; iCell++)  {  
         if (shift(iCell) > -1) {
            if (debug)  {cout << "iCell" << iCell;}

            // get complete phasefactors, reduced at boundaries
            redPhas = fourierReduce (iQ,(latShifts(iCell)),supCell, debug);
            if (shift(iCell) == iCell)  { //no double calculation
               redPhas *= .5;
               invPhas = redPhas;
            }  else  {
               invPhas = fourierReduce 
                            (iQ, (latShifts(shift(iCell))), supCell, debug);
            }
         
            // ---check hermiticity phases
            fract = (redPhas.conj ().transpose () - invPhas).norm ();
            if (fract > .01 * nMode && redPhas.norm () < nMode / 2.)  { 
               cout << iCell << "iCell Warning! phases " 
                    << int (100. * fract / (redPhas.norm () + invPhas.norm ()))
                    << "% anti-hermitian." << endl;
                 
               // print current lattice shift, qPoint and phases for debugging
               if (debug)  {
                  cout << latShifts(iCell) << " latShifts " 
                       << latShifts(shift(iCell)) << qPoints.kVec (iQ) 
                       << "qP, jCell " << shift(iCell) << endl << "phases:"; 
                  redPhas.print (true); 
                  invPhas.conj ().transpose ().print (true);
               }
            } 

            // combine with real dynMat: not hermitian when added separately
            nDM += redPhas * dmC(iCell).real () //is real
                +  invPhas * dmC(shift(iCell)).real (); // SxMatrix * SxMatrix
         }
      }  

      //--check hermiticity new dynMats
      asym = (nDM - nDM.transpose ().conj ()).norm ();
      if (asym * 20. > nDM.norm ()) {  
         cout << "iQ" << iQ << "Warning! dynMatOut " 
              << int (50. * asym / nDM.norm ()) << "% anti-hermitian: noise?!"
              << " hermitized. norm" << nDM.norm () << endl; 

         //-hermitize!
         nDM += nDM.transpose ().conj ();
         nDM *= .5;
      }
      eig(iQ) = nDM.eigensystem ();

      // check drift again (assume first qPoint is Gamma)
      if (iQ == 0)  {
         cout << " drift? (3 e.v.@q=0) " << eig(iQ).vals (0) 
              << eig(iQ).vals (1) << eig(iQ).vals (2) << endl;
      }
   }
   cout << "New dynMats constructed & Hermiticity checked" << endl; 
   return eig;
}

void SxPhonon::printFull (
      const SxArray<SxMatrix<Complex16>::Eigensystem> &eigs) const
{
  cout << endl << "print output for sorting code (qPoints, freqs, vecs)" << endl;

  // check if matrix provided not empty
  SX_CHECK (eigs.getSize () > 0);

  FILE *file;
  int iQ, nQ=qPoints.nk, i, j;
  int nMode = int(eigs(0).vals.getSize ());

  file = fopen("qPoints","w");
  for (iQ=0; iQ<nQ; iQ++) fprintf(file,"%f  %f  %f\n",qPoints.kVec(iQ)(0),qPoints.kVec(iQ)(1),qPoints.kVec(iQ)(2));
  fclose(file);

  // --copy squared frequencies (meV)^2 (sorted)
  SxMatrix<Double> freq (nMode, nQ); // :(iMode,iQ)
  for (iQ = 0; iQ < nQ; iQ++)  { freq.colRef (iQ) <<= eigs(iQ).vals; }

  // get sqrt, conserve sign: a / sqrt (|a|)
  freq = sgnSqrt (freq); //meV

  file = fopen("freqs","w");
  for (iQ=0; iQ<nQ; iQ++) {
    for (i=0; i<nMode; i++) fprintf(file,"%f  ",freq(i,iQ));
    fprintf(file,"\n");
  }
  fclose(file);

  file = fopen("vecs","w");
  for (iQ=0; iQ<nQ; iQ++) {
    for (i=0; i<nMode; i++) {
      for (j=0; j<nMode; j++) {
        fprintf(file,"%15.10f %15.10f    ",eigs(iQ).vecs(i,j).re,eigs(iQ).vecs(i,j).im);
      }
      fprintf(file,"\n");
    }
  }
  fclose(file);
}

void SxPhonon::printFreq (
      const SxArray<SxMatrix<Complex16>::Eigensystem> &eigs,
      const SxString &sxbOut,
            ofstream &output) const
{  
   cout << "printFreq" << endl;
   
   // check if matrix provided not empty
   SX_CHECK (eigs.getSize () > 0);

   // --check if square matrices
   int nQ = int(eigs.getSize());
   for (int iQ = 0; iQ < nQ; iQ++)  {
      if (eigs(iQ).vecs.nCols () != eigs(iQ).vecs.nRows ())  {
         cout << "Error in SxPhonon::getFreq! "
              << "Non-square dynMat for q-points no. " << iQ << endl
              << "no. rows is " << eigs(iQ).vecs.nRows () 
              << "no. cols is" << eigs(iQ).vecs.nCols () << endl;
         SX_QUIT;
      }
   }

   // --copy squared frequencies (meV)^2 (sorted)
   int nMode = int(eigs(0).vals.getSize ());
   SxMatrix<Double> freq (nMode, nQ); // :(iMode,iQ)
   for (int iQ = 0; iQ < nQ; iQ++)  { freq.colRef (iQ) <<= eigs(iQ).vals; }

   // get sqrt, conserve sign: a / sqrt (|a|)
   freq = sgnSqrt (freq); //meV

   // --count number of negative frequencies^2
   int nImagFreq = 0;
   SxVector<Double>::Iterator freqIt;
   for (freqIt = freq.begin (); freqIt != freq.end (); freqIt++)  
   { if (*freqIt < -1.e-2)  { nImagFreq++; } }
   
   // --find frequency and eigenvector of sQ, sMod OR the smallest ones
   int rQ = sQ, 
       rMod = sMod;
   double sFreq;
   SxVector<Complex16> sVec;
   if ((sMod >= 0) && (sQ >= 0))  {

      // check that indices are not too large
      if (sMod >= freq.nRows ())  {
         cout << "Oops sMod too large, 0 taken instead." << endl;
         rMod = 0;
      }
      if (sQ >= freq.nCols ())  {
         cout << "Oops sQ too large, 0 taken instead." << endl;
         rQ = 0;
      }
      sFreq = freq(rMod, rQ);
   }  else  {// find minimum 
      int iMinFreq;
      sFreq = freq.minval (&iMinFreq);
      rQ = iMinFreq % nQ;
      rMod = iMinFreq / nQ; //should be 0
   }
   sVec = eigs(rQ).vecs.colRef (rMod); //colwise
   
   // write data of special q-point to output
   double vol = primStr.cell.volume * cube (scal);
   output << "#There are " << nImagFreq << "imaginary frequencies for vol "
          << vol << endl << "#The frequency and eigenvector of q-point " << rQ 
          << " and mode " << rMod << endl << "#are resp. " 
          << sFreq << "(meV) and " << sVec << endl;

   // ---print frequencies to sxb file
   try  {
      SxBinIO io (sxbOut, SxBinIO::BINARY_WRITE_ONLY);

      // ---first write header frequencies: SxMatrix (nMode, nQ)
      vol *= 3. / nMode; //vol. per atom
      io.write ("volume", vol);
      io.addDimension ("nMode", nMode);
      io.addDimension ("nQ", nQ);
      io.write ("frequencies", SxMatrix<Double> (nMode, nQ), "nMode", "nQ");

      // ---then write data
      io.setMode (1); //SxBinIO::WRITE_DATA
      io.write ("volume", vol);
      io.write ("frequencies", freq, "nMode", "nQ");
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   // --write phonon dispersions to output
   if (nQ > 500)  {
      output << "#no Phonon dispersion: too many qPoints " << endl; 
   }  else  {
      int iMode;
      output << "#Phonon dispersion (meV) for provided qPoints " << endl; 
      for (int iQ = 0; iQ < nQ; iQ++)  {
         output << iQ << " ";
         for (iMode = 0; iMode < nMode; iMode++) {
            output << freq(iMode, iQ) << " ";
         }
         output << endl;
      }
   }
}

SxArray<SxAtomicStructure> SxPhonon::readRefine (
      const SxString &input,
      const SxArray<Coord> &latShifts) //includes supCell
{ 
   cout << "readRefine" << endl;
  
   // check that input not empty
   if (input == "")  {
      cout << "Error in SxPhonon::readRefine! empty input provided" << endl;
      SX_EXIT;
   }
   
   // output of read process
   SxAtomicStructure forc, dTau;
   try  {
      // input to parser
      SxParser parser;
      SxParser::Table parTab = parser.read (input);
      
      // ---read supStr, forces and dTau
      SxAtomicStructure supStr;
      supStr = SxAtomicStructure (parTab->getGroup ("structure")); 
      forc = SxAtomicStructure (parTab->getGroup ("structure"), "force");
      dTau = SxAtomicStructure (parTab->getGroup ("structure"), "dTau");
      cout << "strs read" << endl 
           << supStr << forc << dTau << endl;

      //check if structure has right size (eigenvecs)
      int nCell = int(exactQ.getSize ()),
          nPrimAt = primStr.getNAtoms (); 
      if (supStr.getNAtoms () != nCell * nPrimAt)  {
         cout << "Error in SxPhonon::readRefine, wrong no. atoms!\n"
              << supStr.getNAtoms () << "is not" << nCell << "*" << nPrimAt;
         SX_EXIT;
      } 

      //get scaling factor
      scal = cbrt (supStr.cell.volume 
                 / primStr.cell.volume / PrecTauR(exactQ.getSize ()));
      cout << "volume scaling" << scal << endl;

      // we have no info on the spins and so can't check the primCell
      // we check the superstructure with the scaled primStr & the latShifts
      SxAtomicStructure newPrimStr = supStr.getPrimStr (primStr.cell * scal);

    //SX_CHECK (newPrimStr == (primStr * scal),
    //          newPrimStr, (primStr * scal)); FIXME typecast to double
      for (int iCell = 0; iCell < nCell; iCell++)  
      { SX_CHECK (supStr + latShifts(iCell) * scal == supStr); }
      //Note: primStr, qP and exQ are NOT reset, only the scaling factor is saved

      // ---read recMasses (no class member)
      SxSymbolTable *specPtr = NULL;
      SxList<double> recMassL;
      SxVector<Double> recMass;
      for (specPtr = parTab->getGroup ("pseudoPot")->getGroup ("species");
           specPtr != NULL; specPtr = specPtr->nextSibling ("species"))  {
         recMassL << (specPtr->get ("reciprocalMass")->toReal ()); 
      }
      recMass = recMassL;

      // check if recMass has right size
      if ( recMass.getSize () != supStr.getNSpecies () )  {
         cout << "Mass mismatch in SxPhonon::read !" << endl
              << "No. masses provided is " << recMass.getSize () << endl
              << "No. species is " << supStr.getNSpecies () << endl;
         SX_EXIT;
      }
 
      // scale with masses BEFORE reorder (symmetries unchanged)
      forc *= sqrt(recMass); 
      dTau /= sqrt(recMass); 
    
      // ---get iPrimAt together
      supStr.epsEqual = primStr.epsEqual * scal; 
      SX_EXIT; SxAtomicStructure splitStr; // next line needs splitSpecies
      //SxAtomicStructure splitStr = supStr.splitSpecies (primStr.cell * scal);
      //TODO: CHECK not just reorder, but also new species info  

      forc.replaceInfo (supStr.atomInfo);
      forc = SxAtomicStructure (splitStr.atomInfo, forc);
      dTau.replaceInfo (supStr.atomInfo);
      dTau = SxAtomicStructure (splitStr.atomInfo, dTau);

      // correct drift TODO apply force filter?
      forc -= (forc.sum () / supStr.getNAtoms ()); 
   }  catch (SxException e)  {
      e.print();
      SX_QUIT;
   } 
   return SxList<SxAtomicStructure> () << forc << dTau; 
}

SxArray<SxVector<Double> > SxPhonon::decompose ( 
      const SxArray<SxAtomicStructure> &strs,
      const SxArray<SxMatrix<Complex16> > &eigVecs,
      const SxArray<Coord> &latShifts) const
{  
   cout << "decompose, experimental!" << endl;
   
   //---reorder dTau, forc to (3 * nAtom, nCell)
   int nCell = int(latShifts.getSize ()) - 3, 
       nMode = int(eigVecs(0).nCols ()), //colWise
       nAtom = nMode / 3; 
   SxVector<Double> forc (strs(0).coords), 
                    dTau (strs(1).coords), //3*nCell*nAtom
                    ref;
   dTau.reshape (3, nCell * nAtom);
   dTau = dTau.transpose (); //:(nCell*nAtom,3)
   forc.reshape (3, nCell * nAtom);
   forc = forc.transpose (); //:(nCell*nAtom,3)
   for (int iDir = 0; iDir < 3; iDir++)  {
      ref = dTau.colRef (iDir); //:(nCell*nAtom)
      ref.reshape (nCell, nAtom);
      ref << ref.transpose (); //:(nAtom, nCell) CHECK non-square
      ref = forc.colRef (iDir);
      ref.reshape (nAtom, nCell); //:(nAtom, nCell)
      ref << ref.transpose ();
   }
   dTau = dTau.transpose ();//:(3,nAtom*nCell)
   dTau.reshape (3 * nAtom, nCell);
   forc = forc.transpose ();
   forc.reshape (3 * nAtom, nCell);

   int iCell;
   SxVector<Complex16> phaseQl (nCell), 
                       resF (nMode),
                       resT (nMode);
   SxArray<SxVector<Double> > vals (nCell);
   for (int iQ = 0; iQ < nCell; iQ++)  {

      // phase factor for lattice shift
      for (iCell = 0; iCell < nCell; iCell++)  
      { phaseQl(iCell) = expI (exactQ(iQ) ^ latShifts(iCell)); } //COMPLEX 
      resF = (forc ^ phaseQl) ^ eigVecs(iQ); //colWise
      resT = (dTau ^ phaseQl) ^ eigVecs(iQ); //(mat ^ vec) ^ mat !

      // -check that eigenvals are real
      cout << "real? vals.imag.norm:" 
           <<  sqrt ((resF / resT).imag ().absSqr ().sum ()) * sqr (AU2MEV);
      vals(iQ) = (- resF / resT) * sqr (AU2MEV); //meV^2//colWise
   }
   return vals;
}

#else /* SX_STANDALONE */

int main (int argc, char **argv)  {

   // command line interface: collect identical options
   SxCLI cli (argc, argv);
   SxString dynMatInput = cli.option("-d|--dynmat","file", 
                      "binary dynamical matrix file").toString ("dynmat.sxb"),
            setInput = cli.option ("-s|--settings", "file", 
                                   "phonon settings").toString ("phononSet.sx");
   SxArray<SxString> refineInput = cli.option ("-r|--refine", "file",
                                "(multiple) refinement data file(s)").toList ();
   double epsE = cli.option ("-e|--epsEqual", "accuracy structure")
                                .toDouble (5.e-3);
   bool printFull = cli.option ("-p|--printFull", "print qpoints, freqs, vecs files needed for sorting code")
                                .toBool ();
   bool debug = cli.option ("-v|--verbose", "print some extra information")
                                .toBool ();
   cli.finalize (); 

   // open output file
   ofstream output;
   output.open ("phonon.out");
   output << "#Input files are " << dynMatInput << " " << setInput 
          << " " << refineInput << endl << "#epsE" << epsE << endl;

   // read the input
   SxPhonon P;
   SxArray<SxMatrix<Complex16>::Eigensystem> inEigs = //:nExactQ,(nMode, nMode)
      P.readInput (setInput, dynMatInput, epsE); 

   // get the lattice shifts + the supercell (the last three)
   SxArray<Coord> latShifts = P.getLatShifts (output);//:nCell; lattice vectors
   
   // fourier transform from exactQ to provided qPoints
   SxArray<SxMatrix<Complex16>::Eigensystem> eigs 
      = P.transformDynMat (inEigs, latShifts, debug);//AUTO resize

   // print the frequencies
   SxString sxbOut ("phonon.sxb");
   P.printFreq (eigs, sxbOut, output);

   // print output for sorting code
   if (printFull) P.printFull(eigs);
   
   if (refineInput.getSize () == 0)  {//no refinement

      // close the output
      output.close ();

      // end the program
      return 0;
   }

   // --refine the frequencies
   int iCell, nCell = int(latShifts.getSize ()) - 3;
   SxArray<SxAtomicStructure> forcdTauStr;
   SxArray<SxMatrix<Complex16> > eigVecs (nCell);
   SxArray<SxVector<Double> > eigVals; //:nCell
   for (iCell = 0; iCell < nCell; iCell++)  //nasty!
   { eigVecs(iCell) = inEigs(iCell).vecs; }

   for (int iVol = 0; iVol < refineInput.getSize (); iVol++)  { 
   
      // read forces + displacements for refinement (volumes) 
      forcdTauStr = P.readRefine (refineInput(iVol), latShifts);
    
      // decompose to the frequencies of the phonon modes (meV^2)
      eigVals = P.decompose (forcdTauStr, eigVecs, latShifts);//:(nExactQ,nMode) 

      for (iCell = 0; iCell < nCell; iCell++)  //nasty!
      { inEigs(iCell).vals = eigVals(iCell); }
    
      // transform to the new q-points
      eigs = P.transformDynMat (inEigs, latShifts, debug);
    
      // print the frequencies
      P.printFreq (eigs, sxbOut + SxString (iVol), output);
   }
   // close the output
   output.close ();

   // end the program
   return 0;
}
#endif /* SX_STANDALONE */

