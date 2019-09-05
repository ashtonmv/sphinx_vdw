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

#include <SxSecondaryStructure.h>
#include <SxParser.h>
#include <SxElemDB.h>
#include <SxHessianOps.h>

//----------------------------------------------------------------------------
//    Constructors and Destructors
//----------------------------------------------------------------------------

SxSecondaryStructure::SxSecondaryStructure ()
{
   //empty
} 

SxSecondaryStructure::SxSecondaryStructure 
(  const SxString &structureFilename, const SxString &periodicityFilename
 ,  const SxString &chainTypeIn,int turns, bool toSymmetrize)
{
   setPeriodicity (periodicityFilename);
   setChainType (chainTypeIn);
   setStructure (structureFilename);

   isInInputOrder = true;

   nPeptides = (int)periodicity(0).getSize ();
   nKPoints = nPeptides;
   nAtomsPP = (int)periodicity.getSize (); 
   nAtomsPU = nAtomsPP; 
   nDoF = nPeptides*nAtomsPP*3;
   nDoFPP = nAtomsPP*3;
   nDoFPU = nDoFPP;
   nTurns = turns;
   nUnits = nPeptides;
   twistRad = 2.*PI*(double)nTurns/(double)nPeptides;
   twistDeg = 360.*(double)nTurns/(double)nPeptides;
   zPitch = tauCart.cell(2)(2)/(double)nPeptides;

   //--- irreducible unit for translational symmetry operations 
   //    contains 4 peptide units for beta-sheets
   if (chainType.contains ("Sheet")) { 
      zPitch *= 2.;
      nDoFPU *= 4;
      nUnits = nUnits/4;
      nAtomsPU *= 4;
   }
   
   symCenter = getSymCenter ();
   
   // --- getting cylindrical coordinates
   tauCyl.copy (tauCart);
   tauCartShifted.copy (tauCart);
   updateCylCoords ();
   updateBMatrix ();
   //if (systemClass == SxString("bulk")) vibRots = 3;
   if (toSymmetrize) symmetrizeCoords ();
   if ( (chainType == "helix" ) ||  (chainType == "FES" ) ) vibRots = 2;
}

//--ugly: the constructor should be merged with that above 
SxSecondaryStructure::SxSecondaryStructure 
(  const SxAtomicStructure &tauIn, const SxSpeciesData &speciesData, 
   const SxString &periodicityFilename
 ,  const SxString &chainTypeIn,int turns, bool toSymmetrize)
{
   int counter = 0;
   setPeriodicity (periodicityFilename);
   setChainType (chainTypeIn);
   SxAtomicStructure inputStructure; inputStructure.copy (tauIn);
  
   tauCart = SxAtomicStructure(inputStructure.cell);
   tauCart.startCreation ();
   masses.startCreation ();
   SxVector3<Double> vec;
   
   for (int is = 0; is < inputStructure.getNSpecies (); is++) {
      for (int ia = 0; ia < inputStructure.getNAtoms (is); ia++) {
         vec = inputStructure.ref (counter);
         tauCart.addAtom (vec);
         for (int i = 0; i < 3; i++)
            vec(i) = speciesData.ionicMass(is);
         masses.addAtom (vec);
         counter++;
      }
   }

   
   masses.endCreation ();
   tauCart.endCreation ();
   //tauCart.copy (tauIn);
   //setStructure (structureFilename);

   isInInputOrder = true;
   nPeptides = (int)periodicity(0).getSize ();
   nKPoints = nPeptides;
   nAtomsPP = (int)periodicity.getSize (); 
   nAtomsPU = nAtomsPP; 
   nDoF = nPeptides*nAtomsPP*3;
   nDoFPP = nAtomsPP*3;
   nDoFPU = nDoFPP;
   nTurns = turns;
   nUnits = nPeptides;
   twistRad = 2.*PI*(double)nTurns/(double)nPeptides;
   twistDeg = 360.*(double)nTurns/(double)nPeptides;
   zPitch = tauCart.cell(2)(2)/(double)nPeptides;
   
   //--- irreducible unit for translational symmetry operations 
   //    contains 4 peptide units for beta-sheets
   
   if (chainType.contains ("Sheet")) {
      zPitch *= 2.;
      nDoFPU *= 4;
      nUnits = nUnits/4;
      nAtomsPU *= 4;
   }

   
   symCenter = getSymCenter ();
   
   // --- getting cylindrical coordinates
   tauCyl.copy (tauCart);
   tauCartShifted.copy (tauCart);
   updateCylCoords ();
   updateBMatrix ();
   //if (systemClass == SxString("bulk")) vibRots = 3;
   if (toSymmetrize) symmetrizeCoords ();
   if ( (chainType == "helix" ) ||  (chainType == "FES" ) ) vibRots = 2;
   
}

 

SxSecondaryStructure::~SxSecondaryStructure ()
{
   // empty
}

//----------------------------------------------------------------------------
//    Member Functions
//----------------------------------------------------------------------------

double SxSecondaryStructure::getZPitch () 
{
   return zPitch;
}

int SxSecondaryStructure::getNPeptides () 
{
   return nPeptides;
}

int SxSecondaryStructure::getNAtomsPP () 
{
   return nAtomsPP;
}

int SxSecondaryStructure::getNDoF () 
{
   return nDoF;
}

int SxSecondaryStructure::getNDoFPIU () 
{
   if (chainType.contains ("Sheet")) 
    return (3*2*nAtomsPP);
   else return (3*nAtomsPP);
}

SxMatrix<Double> SxSecondaryStructure::getFullCartHessian () 
{
   return hessianCart;
}

void SxSecondaryStructure::switchToPeptideOrder ()
{
   if (isInInputOrder) {
      peptideOrderTau ();
      isInInputOrder = false;
   }
}

void SxSecondaryStructure::peptideOrderTau ()
{
   int pep, atoms, counter;

   SxAtomicStructure newTauCart; newTauCart.copy (tauCart);
   SxAtomicStructure newTauCartShifted; newTauCartShifted.copy (tauCartShifted);
   SxAtomicStructure newTauCyl; newTauCyl.copy (tauCyl);
   SxAtomicStructure newMasses; newMasses.copy (masses);

   counter = 0;
   for (pep = 0; pep < nPeptides; pep++) {
         for (atoms = 0; atoms < nAtomsPP; atoms++) {
            newTauCart.setAtom(0, counter,
               tauCart.ref((int)periodicity(atoms)(pep) - 1));
            newTauCartShifted.setAtom(0, counter,
               tauCartShifted.ref((int)periodicity(atoms)(pep) - 1));
            newTauCyl.setAtom(0, counter,  
               tauCyl.ref((int)periodicity(atoms)(pep) - 1));
            newMasses.setAtom(0, counter,  
               masses.ref((int)periodicity(atoms)(pep) - 1));
            counter++;
         }
   }

   masses.copy (newMasses); tauCart.copy (newTauCart); tauCyl.copy (newTauCyl);
   tauCartShifted.copy (newTauCartShifted);
}

void SxSecondaryStructure::switchToInputOrder ()
{
   if (!isInInputOrder) {
      inputOrderTau ();
      isInInputOrder = true;
   }
}

void SxSecondaryStructure::inputOrderTau ()
{
   int i;
   
   SxAtomicStructure newTauCart; newTauCart.copy (tauCart);
   SxAtomicStructure newTauCartShifted; newTauCartShifted.copy (tauCart);
   SxAtomicStructure newTauCyl; newTauCyl.copy (tauCyl);
   SxAtomicStructure newMasses; newMasses.copy (masses);
   SxArray<int> per(2);
   
   for (i = 0; i < nPeptides*nAtomsPP; i++) {
      per = getPeriodicIndices(i);
      newTauCart.setAtom (0, i, tauCart.ref (per(0) + per(1)*nAtomsPP));
      newTauCartShifted.setAtom (0, i, 
            tauCartShifted.ref (per(0) + per(1)*nAtomsPP));
      newTauCyl.setAtom (0, i, tauCyl.ref (per(0) + per(1)*nAtomsPP));
      newMasses.setAtom (0, i, masses.ref (per(0) + per(1)*nAtomsPP));
   }

   masses.copy (newMasses); tauCart.copy (newTauCart); tauCyl.copy (newTauCyl);
   tauCartShifted.copy (newTauCartShifted);
}   

void SxSecondaryStructure::setPeriodicity(const SxString &fileName)
{
   FILE *fp = NULL;
   const int BUFLEN = 10240;
   int  j, nP;
   char buffer[BUFLEN];
   char *bs = NULL;
   bool isLine;
   SxVector <Double> line;

   if ( !(fp = fopen (fileName.ascii (), "r")) ) {
      sxprintf ("Can't open file %s", fileName.ascii ());
      SX_EXIT;
   }
   
   periodicity.resize (0);
   isLine = true;
   nP = 0;

   while (isLine) {
      bs = fgets (buffer, 40*3*10, fp);
      SxString tk(bs);
      SxList<SxString> list;
      list = tk.tokenize (' ');
      if (list.getSize() != 0) {
         if (nP == 0) {
            nP = (int)list.getSize ();
            line.resize (nP);
         }
         else {
            if (list.getSize () != nP) {
               cout << "Error while reading periodicity-file:" 
                    << "Wrong number of peptides in one line"
                    << endl;
               fflush (stdout);
               SX_QUIT;
            }
         }

         for (j = 0; j < nP; j++)  {
            line(j) = list(j).toDouble ();
         }

         SxVector<Double> newLine;
         newLine.copy (line);
         periodicity.append (newLine);
      } else {
         isLine = false;
      }
   }

   fclose (fp); 
}


void SxSecondaryStructure::setChainType (const SxString &chainTypeIn) 
{
   chainType = SxString (chainTypeIn);
}

void SxSecondaryStructure::setStructure (const SxString &fileName) 
{
   
   // --- loading in cartesian structure, masses and cell 
   //     from FHI98 structure mass
   SxBinIO io;
   SxList<SxString> speciesNameList;
   SxList<SxList<SxVector3<Double> > > tauList;
   int i, is, ia;
   SxVector3<Double> massVec;
   SxMatrix3<Double> cell;
   (cout << "getting element table ..." << endl).flush ();
   SxElemDB elemDB;
   (cout << "ready" << endl).flush ();
   
   io.open (fileName, SxBinIO::ASCII_READ_ONLY);
   tauList = io.loadStructureFHI98 (); io.close ();
   io.open (fileName, SxBinIO::ASCII_READ_ONLY);
   speciesNameList = io.loadSpeciesFHI98 ();io.close ();
   io.open (fileName, SxBinIO::ASCII_READ_ONLY);
   cell = io.loadCellFHI98 ();io.close ();

   tauCart = SxAtomicStructure (cell);
   masses = SxAtomicStructure (cell);
   tauCart.startCreation ();
   masses.startCreation ();
   
   for (is = 0; is < tauList.getSize (); is++) {
      for (ia = 0; ia < tauList(is).getSize (); ia++) {
         tauCart.addAtom (tauList(is)(ia));
         for (i = 0; i < 3; i++)
            massVec(i) = elemDB.getAtomicWeight(speciesNameList(is));
         masses.addAtom (massVec);
      }
   }
   tauCart.endCreation ();
   masses.endCreation ();
}

void SxSecondaryStructure::setHessian 
    (const SxMatrix<Double> &hessian, int /*supercellCutoff*/, bool applyAveraging) 
{
   int i;
   double sym;

   //--- for projection to internal coordinates
   switchToInputOrder ();
   SxArtifactFilter f;
   f.set (tauCart, SxString("z"), false);
   
   hessianCart = SxMatrix<Double> (nDoF, nDoF);
   hessianCyl = SxMatrix<Double> (nDoF, nDoF);
   
   blockArrayCart =  getInitialBlockArray ();
   blockArrayCyl = getInitialBlockArray ();
   
   switchToPeptideOrder ();
   
   SxMatrix<Double> hessianIn (nDoF, nDoF);
   SxMatrix<Double> difference (nDoF, nDoF);
   hessianIn.copy (hessian);
   
   for (i = 0; i < 2; i++) {

      blockArrayCart = getPeptideOrderedBlocks (hessianIn);
      hessianCart = getHessianFromBlocks(blockArrayCart);
      //--- in principle one could treat an FES as an helix here 
      //    should be unified if there's time left
      if ((chainType == "helix")) 
         hessianCyl = getCylHessian (hessianCart, tauCyl);
      if ((chainType == "FES") || (chainType.contains("Sheet")))
         hessianCyl = getFESHessian (hessianCart);

      if (chainType.contains ("Sheet")) 
         hessianCyl = applyBSheetSymmetry (hessianCyl);

      blockArrayCyl = getBlocks (hessianCyl);
     // blockArrayCyl = applySupercellCutoff (blockArrayCyl, supercellCutoff);
      blockArrayCyl = applySymmetry (blockArrayCyl, applyAveraging);
      hessianCyl = getHessianFromBlocks (blockArrayCyl);

      // --- (cart -> fes) = (fes -> cart) ^ -1
      if ((chainType == "FES")  || (chainType.contains("Sheet")))
         hessianCart = getFESHessian (hessianCyl);
      if (chainType == "helix") 
         hessianCart = getCartHessian (hessianCyl, tauCartShifted);
   
         
      blockArrayCart = getBlocks (hessianCart);
      hessianCart = getHessianFromBlocksInputOrder (blockArrayCart);
      
      //--- Does the Hessian still change ?
      //    (difference should be close to zero in the second iterative step)
      difference = hessianCart - hessianIn;
      sym = difference.absSqr ().sum ();
      if (i == 1) 
         cout << "Difference (should be close to zero): " << sym << endl;

      hessianCart = getSymmetrizedMatrix (hessianCart);
      //--- applying artifact filter commented out (should be user defined)
      //switchToInputOrder ();
      //hessianCart = (f | hessianCart);
      //switchToPeptideOrder ();
      hessianIn.copy (hessianCart);
   }

   switchToInputOrder ();
   if (chainType == "helix") 
      hessianCyl = getCylHessian (hessianCart, tauCyl);
   if (chainType == "FES") 
      hessianCyl = getFESHessian (hessianCart);


   
   if (! (chainType.contains("Sheet")) ) {
      computeKMesh ();
   }
      computeNormalDisplacements ();
}

void SxSecondaryStructure::setHessianFromRefinement 
    (const SxMatrix<Double> &hessianIn, const SxMatrix<Double> &basisHessianIn, 
     int supercellCutoff, bool applyAveraging) 
{
   int nDoFPIU = getNDoFPIU ();
   SxMatrix<Double>  deviation(nDoF, nDoF);
   SxMatrix<Double>  basisHessian(nDoF, nDoF);
   SxMatrix<Double>  hessian(nDoF, nDoF);
   SxMatrix<Double> normalProjections (nDoF, nDoF);
   SxMatrix<Double> reducedHessian (nDoFPIU, nDoFPIU);
   SxMatrix<Double> gammaHessian (nAtomsPU * 3, nAtomsPU *3);
   SxMatrix<Double> block (nDoFPIU, nDoFPIU);
   SxMatrix<Double>::Eigensystem eigBasis;
   SxVector<Double> column (nDoF);
   SxVector<Double> orth (nDoF);
   SxVector<Double> longColumn (nDoF);
   SxVector<Double> workAroundPsi (nDoFPIU);
   SxArray<SxArray<SxMatrix<Double> > > blockArray;
   int i, j, k;
   bool toChange;

   basisHessian.copy (basisHessianIn);
   hessian.copy (hessianIn);
  
   SxArray<int> fromRefineHessian (1);
   fromRefineHessian(0) = -1;

   //--- this flag should be taken out of the code one day
   //    however it remains for safety reasons
   bool oldCalc = false;

   blockArray = getInitialBlockArray();
   basisHessian = getSymmetrizedMatrix(basisHessian);
   blockArray = getPeptideOrderedBlocks (basisHessian);

   gammaHessian.set (blockArray(0)(0));
   for (i = 0; i < nDoFPIU; i++) {
      for (j = 0; j < nDoFPIU; j++) {
         reducedHessian(i, j) = gammaHessian (i, j);
      }
   }
   
   eigBasis = reducedHessian.eigensystem ();
   
   //--- this procedure assures that the results do not suffer 
   ///   from the fact, that different platforms use different 
   //    math-libs
   if (!oldCalc) {
      workAroundPsi.resize (nDoFPIU);
      for (i = 0; i < nDoFPIU; i++) {
         workAroundPsi.copy(eigBasis.vecs.colRef(i));
         if (workAroundPsi(0) <= 0.)
            workAroundPsi = -workAroundPsi;
         eigBasis.vecs.colRef(i).set(workAroundPsi);
      }
   }
   
   // --- a matrix with normed displacement vectors according to the 
   //     input basis is generated 
   for (i = 0; i < nDoF; i++) {
      if (i < nDoFPIU) {
         for (j = 0; j < nDoFPIU; j++) 
            column(j) = eigBasis.vecs.colRef(i)(j);
         
         for (j = nDoFPIU; j < nDoF; j++) 
            column(j) = 0.;
         
         column = column/sqrt(column.absSqr().sum());
         for (j = 0; j < i; j++) {
            orth.copy (deviation.colRef(j));
            column = column - (column^orth).chop()*orth;
            column = column/sqrt(column.absSqr().sum());
         }
         deviation.colRef(i).set(column);
      } 
      else 	deviation.colRef(i).set (0.);
   }
   
   for (i = 0; i < nDoF; i++) {
      for (j = 0; j < nDoF; j++) 
         longColumn(j) = hessian (i, j);
      
      longColumn = getPeptideOrderedVec(longColumn);
      
      for (j = 0; j < nDoF; j++) 
         hessian(i, j) = longColumn(j);
   }
         
   if (fromRefineHessian(0) != -1) {
      
      for (i = 0; i < nDoFPIU; i++) {
         toChange = true;
         for (j = 0; j < fromRefineHessian.getSize (); j++) {
            if (fromRefineHessian(j) == (i + 1)) toChange = false; 
         }
         if (toChange) {
            eigBasis.vecs.colRef(i).set (0.);
            eigBasis.vecs.colRef(i)(i) = 1.;
         }
      }
   }
   
   blockArray = getInitialBlockArray ();
   blockArray = getBlocks(deviation);
   //applyScrewSymmetry(&blockArray, 0, periodicity);
   for (i = 1; i < nUnits; i++) 
      blockArray(i)(i).copy(blockArray(0)(0));
   
   deviation = getHessianFromBlocks(blockArray);
   
   for (j = 0; j < nDoF; j++) {
      //column.copy(blockArray(0)(i).colRef(j));
      
      for (k = 0; k < nDoF; k++) {
         column(k) = hessian(j, k);
      }
      
      if (fromRefineHessian(0) != -1) {
         SxVector<Double> toReplace (nDoF);
         toChange = true;
         for (int ii = 0; ii < fromRefineHessian.getSize (); ii++) {
            if (fromRefineHessian(ii) == (j + 1)) toChange = false; 
         }
         if (toChange) {
            toReplace.copy(deviation.colRef(j));
            toReplace = getInputOrderedVec(toReplace);
            toReplace = basisHessian^toReplace;
            toReplace = getPeptideOrderedVec(toReplace);
            column.copy (toReplace);
         }
      }
      
      for (k = 0; k < nDoF; k++) 
         orth(k) = (deviation.colRef(k)^column).chop ();
      for (k = 0; k < nDoF; k++) 
         normalProjections(j, k) = orth(k);
   }
   hessian = deviation ^ (normalProjections^(deviation.transpose ()));
   
   blockArray = getInitialBlockArray ();
   blockArray = getBlocks(hessian);

   hessian = getHessianFromBlocksInputOrder(blockArray);

   setHessian (hessian, supercellCutoff, applyAveraging); 
}


void SxSecondaryStructure::symmetrizeCoords () 
{
   double r, z1, z2, delta, zLat, scorr, 
          phi, phi1, phi2, globalPhi;
   int i, j;

   switchToPeptideOrder ();


   if (chainType == "helix") {
      zLat = zPitch*(double) nPeptides;
   
      //--- Symmetrizing r coordinate
      for (i = 0; i < nAtomsPP; i++) {
         r = 0.;
         for (j = 0; j < nPeptides; j++) 
            r += tauCyl.ref(j*nAtomsPP + i)(1)/(double)nPeptides;
         for (j = 0; j < nPeptides; j++) 
            tauCyl.ref(j*nAtomsPP + i)(1) = r;
      }
   
      //--- Symmetrizing z coordinate;

      //--- if zPitch is not initialized by lattice parameter it 
      //    is obtained by lattice paramter
      if ((zPitch - 0.) < 1e-15)  {
         for (i = 0; i < nAtomsPP; i++) {
            z1 = 0.;
            for (j = 0; j < nPeptides ; j++) 
               z1 += tauCyl.ref(j*nAtomsPP + i)(2)/nPeptides;
            for (j = 0; j < nPeptides ; j++) 
               tauCyl.ref(j*nAtomsPP + i)(2) = z1;
         }
      }   else {

         for (int it = 0; it < 5; it++) {
            for (i = 0; i < nAtomsPP; i++) {
               delta = 0.;
            
               for (j = 0; j < nPeptides ; j++) {
                  z1 = tauCyl.ref (j*nAtomsPP + i)(2);
                  z2 = tauCyl.ref(((j+1) % nPeptides)*nAtomsPP + i)(2);
            
                  if (fabs (z1 - z2) < 1.5*zPitch) scorr = 0.;
                  else {
                     if (z2 > z1) scorr = -1.;
                     else  scorr = 1. ;
                  }
                  if (z1 < z2 + scorr*zLat) {
                     delta = zPitch - ((z2 + scorr*zLat) - z1); 
                     z2 += delta/2.;
                     z1 -= delta/2.;
                  } else {
                     delta = zPitch - (z1 - (z2 + scorr*zLat)); 
                     z2 -= delta/2.;
                     z1 += delta/2.;
                  }
                  tauCyl.ref(j*nAtomsPP + i)(2) = z1;
                  tauCyl.ref(((j+1) % nPeptides)*nAtomsPP + i)(2) = z2;
               }
            }
         }
      }
      
      if ((zPitch - 0.) < 1e-15)  {
         for (i = 0; i < nAtomsPP; i++) {
            phi = 0.;
            for (j = 0; j < nPeptides ; j++) 
               phi += tauCyl.ref(j*nAtomsPP + i)(0)/nPeptides;
            for (j = 0; j < nPeptides ; j++) 
               tauCyl.ref(j*nAtomsPP + i)(0) = phi;
         }
      }
   
      else {
         phi1 = tauCyl.ref(0)(0);
         globalPhi = 2.*PI/(double)nPeptides*double(nTurns);
         phi2 = tauCyl.ref(nAtomsPP)(0);
         
         phi2 += globalPhi;
         if (phi2 > PI) phi2 -= 2*PI;
         if (fabs(phi1 - phi2) <  0.5*globalPhi) globalPhi *= -1.;

         for (i = 0; i < nAtomsPP; i++) {
            for (j = 0; j < nPeptides - 1; j++) {
               if (j == 0) 
                  tauCyl.ref(i)(0) 
                     = tauCyl.ref((nPeptides - 1)*nAtomsPP + i)(0)   
                     + globalPhi;
               else 
                  tauCyl.ref(j*nAtomsPP + i)(0)
                     =  tauCyl.ref((j-1)*nAtomsPP + i)(0)
                     +  globalPhi;
            
               if (tauCyl.ref(j*nAtomsPP + i)(0) > PI) 
                  tauCyl.ref(j*nAtomsPP + i)(0) -= 2.*PI;   

               if (tauCyl.ref(j*nAtomsPP + i)(0) < PI) 
                  tauCyl.ref(j*nAtomsPP + i)(0) += 2.*PI;   
            }
         }
      }
   }
   
   if (chainType == "pSheet") {
      zLat = zPitch*(double) (nPeptides/2);
   
      //--- Symmetrizing r coordinate
      for (i = 0; i < nAtomsPP; i++) {
         r = 0.;
         for (j = 0; j <= (nPeptides - 1)/2; j++) 
            r += tauCyl.ref(j*2*nAtomsPP + i)(1)/(double)nPeptides*2.;
         for (j = 0; j <= (nPeptides - 1)/2; j++) 
            tauCyl.ref(j*2*nAtomsPP + i)(1) = r;
      }
      //--- Symmetrizing r coordinate
      for (i = 0; i < nAtomsPP; i++) {
         r = 0.;
         for (j = 0; j <= (nPeptides - 1)/2; j++) 
            r += tauCyl.ref((j*2 + 1)*nAtomsPP + i)(1)/(double)nPeptides*2.;
         for (j = 0; j <= ( nPeptides -1)/2; j++) 
            tauCyl.ref((j*2 + 1)*nAtomsPP + i)(1) = r;
      }
      
      //--- Symmetrizing z coordinate
      //    accounts for asymmetric pitch 
      double zPitch1 = 0.;
      for (i = 0; i < nAtomsPP; i++) {
         zPitch1 +=  
            (tauCyl.ref(i)(2) -tauCyl.ref(nAtomsPP*2 + i)(2))/(double)nAtomsPP;
      }
      double zPitch2 = 0.;
      for (i = 0; i < nAtomsPP; i++) {
         zPitch2 +=  
            (tauCyl.ref(i + nAtomsPP)(2) 
             - tauCyl.ref(nAtomsPP*3 + i)(2))/(double)nAtomsPP;
      }
      double signum = -1.;
      if (zPitch1 > 0.) signum = 1.;
      double corr = (2.*zPitch - (fabs(zPitch1 + zPitch2)))/2.;

      zPitch1 = signum* ( corr/2. + fabs(zPitch1));
      zPitch2 = signum* ( corr/2. + fabs(zPitch2));
      
      for (int it = 0; it < 5; it++) {
         for (i = 0; i < nAtomsPP; i++) {
            delta = 0.;
            
            for (j = 0; j <= (nPeptides - 1)/2 ; j++) {
               z1 = tauCyl.ref ((j*2)*nAtomsPP + i)(2);
               z2 = tauCyl.ref((((j+1)*2) % (nPeptides))*nAtomsPP + i)(2);
            
               if (fabs (z1 - z2) < 1.5*zPitch1) scorr = 0.;
               else {
                  if (z2 > z1) scorr = -1.;
                  else  scorr = 1. ;
               }
               if (z1 < z2 + scorr*zPitch1*(double)(nPeptides/2)) {
                  delta = zPitch 
                     - ((z2 + scorr*zPitch1*(double)(nPeptides/2)) - z1); 
                  z2 += delta/2.;
                  z1 -= delta/2.;
               } else {
                  delta = zPitch 
                     - (z1 - (z2 + scorr*zPitch1*(double)(nPeptides/2))); 
                  z2 -= delta/2.;
                  z1 += delta/2.;
               }
               tauCyl.ref(j*2*nAtomsPP + i)(2) = z1;
               tauCyl.ref((((j+1)*2) % (nPeptides))*nAtomsPP + i)(2) = z2;
            }
            
            for (j = 0; j <= (nPeptides/2 - 1); j++) {
               z1 = tauCyl.ref (((j*2) + 1)*nAtomsPP + i)(2);
               z2 = tauCyl.ref(((((j+1)*2) % nPeptides) + 1)*nAtomsPP + i)(2);
            
               if (fabs (z1 - z2) < 1.5*zPitch2) scorr = 0.;
               else {
                  if (z2 > z1) scorr = -1.;
                  else  scorr = 1. ;
               }
               if (z1 < z2 + scorr*zPitch2*(double)(nPeptides/2)) {
                  delta = zPitch 
                     - ((z2 + scorr*zPitch2*(double)(nPeptides/2)) - z1); 
                  z2 += delta/2.;
                  z1 -= delta/2.;
               } else {
                  delta = zPitch 
                     - (z1 - (z2 + scorr*(double)(nPeptides/2))); 
                  z2 -= delta/2.;
                  z1 += delta/2.;
               }
               tauCyl.ref((j*2 + 1)*nAtomsPP + i)(2) = z1;
               tauCyl.ref(((((j+1)*2) % nPeptides) + 1)*nAtomsPP + i)(2) = z2;
            }
         }
      }
      
      //--- Symmetrizing phi ordinate
   
      phi1 = tauCyl.ref(0)(0);
      globalPhi = PI;
      phi2 = tauCyl.ref(2*nAtomsPP)(0);
         
      phi2 += globalPhi;
      if (phi2 > PI) phi2 -= 2*PI;
      if (fabs(phi1 - phi2) <  0.5*globalPhi) globalPhi *= -1.;
      
      for (i = 0; i < nAtomsPP; i++) {
         for (j = 0; j <= (nPeptides - 1)/2; j++) {
            tauCyl.ref(j*2*nAtomsPP + i)(0)
               =  (double)j*globalPhi +  tauCyl.ref(i)(0);
            
            if (tauCyl.ref(j*2*nAtomsPP + i)(0) > PI) 
               tauCyl.ref(j*2*nAtomsPP + i)(0) -= 2.*PI;   

            if (tauCyl.ref(j*2*nAtomsPP + i)(0) < PI) 
               tauCyl.ref(j*2*nAtomsPP + i)(0) += 2.*PI;   
         }
         for (j = 0; j <= (nPeptides - 1)/2; j++) {
            tauCyl.ref((j*2 + 1)*nAtomsPP + i)(0)
               =  (double)j*globalPhi +  tauCyl.ref(nAtomsPP + i)(0);
            
            if (tauCyl.ref((j*2 + 1)*nAtomsPP + i)(0) > PI) 
               tauCyl.ref((j*2 + 1)*nAtomsPP + i)(0) -= 2.*PI;   

            if (tauCyl.ref((j*2 + 1)*nAtomsPP + i)(0) < PI) 
               tauCyl.ref((j*2 + 1)*nAtomsPP + i)(0) += 2.*PI;   
         }
      }
   
   }
   if (chainType == "FES") {
      //todo
   }
      updateCartCoords ();
}

void SxSecondaryStructure::updateCylCoords ()
{
   int i;
   double x, y, z, r, phi;
   symCenter = getSymCenter ();
   
         for (i = 0; i < tauCart.nTlAtoms; i++) {
            x = tauCart.ref(i)(0) - symCenter(0);
            y = tauCart.ref(i)(1) - symCenter(1);
            z = tauCart.ref(i)(2);
            phi = atan2(y, x);
            //phi += PI; 
            r = sqrt(x*x + y*y);
            tauCyl.ref(i)(0) = phi;
            tauCyl.ref(i)(1) = r;
            tauCyl.ref(i)(2) = z;   
         }
}
   
void SxSecondaryStructure::updateCartCoords ()
{
   int i;
   double x, y, z, r, phi;
   
   for (i = 0; i < tauCart.nTlAtoms; i++) {
      
      phi = tauCyl.ref(i)(0);
      r   = tauCyl.ref(i)(1);
      z   = tauCyl.ref(i)(2);
      x = cos(phi)*r;
      y = sin(phi)*r;
            
      tauCart.ref(i)(0) = x + symCenter(0);
      tauCart.ref(i)(1) = y + symCenter(1);
      tauCart.ref(i)(2) = z;

      tauCartShifted.ref(i)(0) = x;
      tauCartShifted.ref(i)(1) = y;
      tauCartShifted.ref(i)(2) = z;

   }
}

void SxSecondaryStructure::updateBMatrix ()
{

   double phi, r;

   B = SxMatrix<Double> (nDoFPP, nDoFPP);
   B.set (0.);
  
   switchToPeptideOrder ();
      
   for (int x = 0; x < nAtomsPP; x++) {
      phi = tauCyl.ref(x)(0);
      r =   tauCyl.ref(x)(1);
      B(3*x + 0, 3*x + 0) = -sin(phi)/r;
      B(3*x + 1, 3*x + 1) = sin(phi);
      B(3*x + 1, 3*x + 0) = cos(phi);
      B(3*x + 0, 3*x + 1) = cos(phi)/r;
      B(3*x + 2, 3*x + 2) = 1.;
   }
}

SxVector<Double> SxSecondaryStructure::
displCartToCyl (const SxVector<Double> &in)
{
   SxVector<Double> returnValue(nDoF);
   double dx, dy, dz;
   double dr, dphi, phi;
   double x, y, r;
   switchToPeptideOrder ();
   for (int i = 0; i < nDoF/3; i++) {
      dx = in(i*3 + 0);
      dy = in(i*3 + 1);
      dz = in(i*3 + 2);
      x = tauCartShifted.ref(i)(0);
      y = tauCartShifted.ref(i)(1);
      
      r = sqrt(x*x + y*y);
      dr = sqrt((x + dx)*(x + dx) + (y + dy)*(y + dy)) - r;
     // - x*y*dx*dy/r*r*r; 

      //L = 1./(1 + y*y/x/x);
      phi = atan2(y, x);
      dphi = atan2(y + dy, x+ dx) - phi ;
     // - (L/x/x- 2.*L*L*y*y/x/x/x/x)*dx*dy; 

      returnValue(i*3 + 0) = dphi;
      returnValue(i*3 + 1) = dr;
      returnValue(i*3 + 2) = dz;
   }
   return returnValue;
}

SxVector<Double> SxSecondaryStructure::
displCylToCart  (const SxVector<Double> &in)
{
   SxVector<Double> returnValue(nDoF);
   double dx, dy, dz;
   double dr, dphi;
   double r, phi;
   switchToPeptideOrder ();
   for (int i = 0; i < nDoF/3; i++) {
      dphi = in(i*3 + 0);
      dr = in(i*3 + 1);
      dz = in(i*3 + 2);
      r = tauCyl.ref(i)(1);
      phi = tauCyl.ref(i)(0);
      
      dx = cos(phi + dphi)*(r + dr) - cos(phi)*r;
     // - sin(phi)*dphi*dr;
      dy = sin(phi + dphi)*(r + dr) - sin(phi)*r;
     // + cos(phi)*dphi*dr;

      returnValue(i*3 + 0) = dx;
      returnValue(i*3 + 1) = dy;
      returnValue(i*3 + 2) = dz;
   }
   return returnValue;
}

SxArray<double> SxSecondaryStructure::getSymCenter () 
{
   int i;
   int size = (tauCart.nTlAtoms);
   SxArray<double> returnValue(2);
   
   returnValue(0) = returnValue (1) =0.;
   
   for (i = 0; i < size; i++) {
      returnValue(0) += tauCart.ref(i)(0);
      returnValue(1) += tauCart.ref(i)(1);
   }
   
   returnValue (0) = returnValue(0)/(double)(size);
   returnValue (1) = returnValue(1)/(double)(size);
   return returnValue;
}
  

SxArray<SxArray<SxMatrix<Double> > > 
SxSecondaryStructure::getInitialBlockArray ()

{
   int i, j;
   SxArray<SxArray<SxMatrix<Double> > >  returnValue (nUnits);
   
   for (i = 0; i < nUnits; i++) {
      returnValue(i) = SxArray<SxMatrix<Double> > (nUnits);
      for (j = 0; j < nUnits; j++) {
         returnValue(i)(j) = SxMatrix<Double> (nDoFPU, nDoFPU);
         returnValue(i)(j).set(0.);
      }
   }
   return returnValue;
}


SxArray<SxArray<SxMatrix<Double> > >  
   SxSecondaryStructure:: getBlocks (const SxMatrix<Double> &Hessian)
{
   int i, j,  cartI, cartJ, atomIndexI, atomIndexJ, iPep, jPep;
   SxArray<SxArray<SxMatrix<Double> > > returnValue (nUnits);
   
   for (iPep = 0; iPep < nUnits; iPep++) {
      returnValue(iPep) = SxArray<SxMatrix<Double> > (nUnits);
      for (jPep = 0; jPep < nUnits; jPep++) {
         returnValue(iPep)(jPep) = SxMatrix<Double> (nAtomsPU*3, nAtomsPU*3);
         returnValue(iPep)(jPep).set(0.);
         for (i = 0; i < nAtomsPU; i++) {
            for (j = 0; j < nAtomsPU; j++) {
               atomIndexI = iPep*nAtomsPU + i;
               atomIndexJ = jPep*nAtomsPU + j;
               for (cartI = 0; cartI < 3; cartI++) {
                  for (cartJ = 0; cartJ < 3; cartJ++) {
                     returnValue(iPep)(jPep)(i*3 + cartI, j*3 + cartJ) 
                        = Hessian (atomIndexI*3 + cartI, atomIndexJ*3 + cartJ);
                  }
               }
            }
         }
      }
   }
   return returnValue;
}


SxArray<SxArray<SxMatrix<Double> > > SxSecondaryStructure::
getPeptideOrderedBlocks (const SxMatrix<Double> &Hessian)
{
   int i, j,  cartI, cartJ, atomIndexI, atomIndexJ, iPep, jPep;
   SxArray<SxArray<SxMatrix<Double> > > returnValue (nUnits);

   for (iPep = 0; iPep < nUnits; iPep++) {
      returnValue(iPep) = SxArray<SxMatrix<Double> > (nUnits);
      for (jPep = 0; jPep < nUnits; jPep++) {
         returnValue(iPep)(jPep) = SxMatrix<Double> (nAtomsPU*3, nAtomsPU*3);
         returnValue(iPep)(jPep).set(0.);
         for (i = 0; i < nAtomsPU; i++) {
            for (j = 0; j < nAtomsPU; j++) {
               atomIndexI = (int)periodicity(i%nAtomsPP)
                                (i/nAtomsPP + (nPeptides/nUnits)*iPep) - 1;
               atomIndexJ = (int)periodicity(j%nAtomsPP)
                                (j/nAtomsPP + (nPeptides/nUnits)*jPep) - 1;
                  
               for (cartI = 0; cartI < 3; cartI++) {
                  for (cartJ = 0; cartJ < 3; cartJ++) {
                     returnValue(iPep)(jPep)(i*3 + cartI, j*3 + cartJ) 
                        = Hessian (atomIndexI*3 + cartI, atomIndexJ*3 + cartJ);
                  }
               }
            }
         }
      }
   }
   return returnValue;
}

SxVector<Double> SxSecondaryStructure::getPeptideOrderedVec 
                 (const SxVector<Double> &in) 
{
   int iPOrder, atoms, pep, counter;
   SxVector<Double> returnValue(nDoF);
   counter = 0;
   for (pep = 0; pep < nPeptides; pep++) {
      for (atoms = 0; atoms < nAtomsPP; atoms++) {
         iPOrder = (int)periodicity(atoms)(pep) - 1;
         returnValue (counter*3 + 0) = in(iPOrder*3 + 0);
         returnValue (counter*3 + 1) = in(iPOrder*3 + 1);
         returnValue (counter*3 + 2) = in(iPOrder*3 + 2);
         counter++;
      }
   }
   return returnValue;   
}
     
SxVector<Double> SxSecondaryStructure::getInputOrderedVec 
                 (const SxVector<Double> &in) 
{
   int i, iPOrder;
   SxVector<Double> returnValue(nDoF);
   SxArray<int> per(2);

   for (i = 0; i < nPeptides*nAtomsPP; i++) {
      per = getPeriodicIndices(i);
      iPOrder = per(0) + per(1)*nAtomsPP;
      returnValue (i*3 + 0) = in(iPOrder*3 + 0);
      returnValue (i*3 + 1) = in(iPOrder*3 + 1);
      returnValue (i*3 + 2) = in(iPOrder*3 + 2);
   }
   return returnValue;   
}
                    

SxArray<SxArray<SxMatrix<Double> > > SxSecondaryStructure::applySupercellCutoff
(const SxArray<SxArray<SxMatrix<Double> > > &blockArray, int supercellCutoff)
{
   int i;
   SxMatrix<Double> block(nAtomsPP*3, nAtomsPP*3);
   SxMatrix<Double> diag (nAtomsPP*3, nAtomsPP*3);
   SxMatrix<Double> offDiag(nAtomsPP*3, nAtomsPP*3);
   SxArray<SxArray<SxMatrix<Double> > > returnValue (blockArray);
   SxMatrix<Double> Binv;

   Binv = B.inverse ();

   for (i = 0; i < nUnits; i++) {
      if ((i > supercellCutoff) && ((nUnits - i) > supercellCutoff))
         returnValue(0)(i).set (0.);
   }
   return returnValue;
}
   

SxMatrix<Double> SxSecondaryStructure::getHessianFromBlocks
(const SxArray<SxArray<SxMatrix<Double> > > &blockArray) 
{
   int i, j, iPep, jPep, iHessian, jHessian;
   int nBlocks = (int)blockArray.getSize ();
   SxMatrix<Double> returnValue(nBlocks*nAtomsPU*3, nBlocks*nAtomsPU*3);
   for (iPep = 0; iPep < nBlocks; iPep++) {
      cout << iPep << endl;
      cout << iPep << endl;
      for (jPep = 0; jPep < nBlocks; jPep++) {
         for (i = 0; i < nAtomsPU*3; i++) {
            for (j = 0; j < nAtomsPU*3; j++) {
               iHessian = (iPep*nAtomsPU)*3 + i;
               jHessian = (jPep*nAtomsPU)*3 + j;
               
               returnValue(iHessian, jHessian) =
                  blockArray(iPep)(jPep)(i, j);
            }
         }
      }
   }
   return returnValue;
}

SxMatrix<Double> SxSecondaryStructure::applyBSheetSymmetry 
     (const SxMatrix<Double> &in)
{
   int j, shiftedIndex;
   SxVector<Double> col (nDoF);
   SxVector<Double> shiftedCol (nDoF);
   SxMatrix<Double> returnValue (nDoF, nDoF);
   
   returnValue.copy (in);
   for (int i = 0; i < 6*nAtomsPP; i++) {
      col.set (in.colRef (i));
      for (j = 0; j < nDoF; j++) {
         shiftedIndex = ( (j + 6*nAtomsPP) % nDoF );
         shiftedCol (shiftedIndex) = col(j);
      }
      returnValue.colRef (i + 6*nAtomsPP).set (shiftedCol);
   }
   return returnValue;
}

SxArray<SxArray<SxMatrix<Double> > > SxSecondaryStructure::applySymmetry
     (const SxArray<SxArray<SxMatrix<Double> > > &blockArray, 
      bool applyAveraging)
{
   int i, j;

   SxArray<SxArray<SxMatrix<Double> > > symmetrizedArray;
   SxArray<SxArray<SxMatrix<Double> > > returnValue;
   symmetrizedArray = getInitialBlockArray ();
   returnValue = getInitialBlockArray ();
   
   SxMatrix<Double> block(nDoFPU, nDoFPU);
   //--- applying screw symmetry by copying blocks from representative 
   //    peptide unit to the neighbors (filling hessian matrix)
  
   if (!applyAveraging) {
      block.copy (blockArray(0)(0));
      
      block = getSymmetrizedMatrix(block);
      returnValue(0)(0).copy (block);

      for (i = 0; i < nUnits; i++) {
          for (j = 0; j < nUnits; j++) {
             if (!((i == 0) && (j == 0))) {
                returnValue(i)(j).copy(blockArray(i)(j));
             }
          }
      } 

      for (i = 0; i < nUnits; i++) {
         for (j = 0; j < nUnits; j++) {
            if (j >= i)
               symmetrizedArray(i)(j).copy(returnValue(0)(j-i));
             else symmetrizedArray(i)(j).copy(returnValue(0)((j-i) 
                      + nUnits));
         }
      }

      //--- applying screw symmetry, by averaging over all peptide units 
   } else {
      
      for (j = 0; j < nUnits; j++) 
         returnValue(0)(j) 
            = blockArray(0)(j)/(double)(nUnits);          
      
      for (i = 1; i < nUnits; i++) {
         for (j = 0; j < nUnits; j++) {
            if (j >= i) 
               returnValue(0)(j-i) = returnValue(0)(j-i) 
                                     + blockArray(i)(j)
                                     / (double)nUnits;
             else 
               returnValue(0)(j-i + nUnits) 
                                     = returnValue(0)(j-i+nUnits)
                                     + blockArray(i)(j)
                                     / (double)nUnits;
         }
      } 
   
      block.copy (returnValue(0)(0));
       block = getSymmetrizedMatrix(block);
      returnValue(0)(0).copy (block);
      
      for (i = 0; i < nUnits; i++) {
         for (j = 0; j < nUnits; j++) {
            if (j >= i) symmetrizedArray(i)(j).
               copy(returnValue(0)(j-i));
             else 
                symmetrizedArray(i)(j).
                   copy(returnValue(0)((j-i) + nUnits));
         }
      } 
   }
   
   for (i = 0; i < nUnits; i++) {
         for (j = 0; j < nUnits; j++) {
            returnValue(i)(j).
               copy((symmetrizedArray)(i)(j));
         }
   } 
   return returnValue;
}
         
SxMatrix<Double> SxSecondaryStructure::getHessianFromBlocksInputOrder
(const SxArray<SxArray<SxMatrix<Double> > > &blockArray)
{
   int i, j, cartI, cartJ;
   SxArray<int> perI, perJ;
   SxMatrix<Double> returnValue(nDoF, nDoF);
   
   for (i = 0; i < nAtomsPU*nUnits; i++) {
      for (j = 0; j < nAtomsPU*nUnits; j++) {

         perI = getPeriodicIndices (i);
         perJ = getPeriodicIndices (j);
         for (cartI = 0; cartI < 3; cartI++) {
            for (cartJ = 0; cartJ < 3; cartJ++) {
               returnValue(3*i + cartI, 3*j + cartJ)
                  = blockArray(perI(1))(perJ(1))
                  (perI(0)*3 + cartI, perJ(0)*3 + cartJ);

            }
         }
      }
   }
   return returnValue;
}


SxMatrix<Double> SxSecondaryStructure::getCylHessian 
(const SxMatrix<Double> &HC, const SxAtomicStructure &t)
{
      
   int i, j, size;
   SxString xphi = SxString("xphi");
   SxString yphi = SxString("yphi");
   SxString xr = SxString("xr");
   SxString yr = SxString("yr");
   size = (int)sqrt((double)HC.getSize ());
   
   SxMatrix<Double> rV(size, size);
   double Hxx, Hyy, Hzz, Hxy, Hyx, Hxz, Hzx, Hyz, Hzy;
   double Hpp, Hrr, Hpr, Hpz, Hrp, Hzp, Hrz, Hzr;


   if (chainType == "helix") {
      for (i = 0; i < size/3; i++) {
         for (j = 0; j < size/3; j++) {

            Hxx = HC(i*3+0, j*3+0); Hxy = HC(i*3+0, j*3+1); 
            Hxz = HC(i*3+0, j*3+2);
            
            Hyx = HC(i*3+1, j*3+0); Hyy = HC(i*3+1, j*3+1); 
            Hyz = HC(i*3+1, j*3+2);
            
            Hzx = HC(i*3+2, j*3+0); Hzy = HC(i*3+2, j*3+1); 
            Hzz = HC(i*3+2, j*3+2);


            Hpp =   getT(xphi, j, t)*getT(xphi, i, t)*Hxx 
               + getT(xphi, j, t)*getT(yphi, i, t)*Hyx    
               + getT(yphi, j, t)*getT(xphi, i, t)*Hxy 
               + getT(yphi, j, t)*getT(yphi, i, t)*Hyy;

            Hrr =   getT(xr, j, t)*getT(xr, i, t)*Hxx 
               + getT(xr, j, t)*getT(yr, i, t)*Hyx    
               + getT(yr, j, t)*getT(xr, i, t)*Hxy 
               + getT(yr, j, t)*getT(yr, i, t)*Hyy;

            Hpr =   getT(xr, j, t)*getT(xphi, i, t)*Hxx 
                  + getT(xr, j, t)*getT(yphi, i, t)*Hyx    
                  + getT(yr, j, t)*getT(xphi, i, t)*Hxy 
                  + getT(yr, j, t)*getT(yphi, i, t)*Hyy;

            Hrp =   getT(xphi, j, t)*getT(xr, i, t)*Hxx 
                  + getT(xphi, j, t)*getT(yr, i, t)*Hyx    
                  + getT(yphi, j, t)*getT(xr, i, t)*Hxy 
                  + getT(yphi, j, t)*getT(yr, i, t)*Hyy;

         
            Hpz = getT(xphi, i, t)*Hxz + getT(yphi, i, t)*Hyz;
            Hrz = getT(xr, i, t)*Hxz + getT(yr, i, t)*Hyz;
            Hzp = getT(xphi, j, t)*Hzx + getT(yphi, j, t)*Hzy;
            Hzr = getT(xr, j, t)*Hzx + getT(yr, j, t)*Hzy;

            rV(i*3+0, j*3+0)=Hpp; rV(i*3+0, j*3+1)=Hpr; rV(i*3+0, j*3+2)=Hpz;
            rV(i*3+1, j*3+0)=Hrp; rV(i*3+1, j*3+1)=Hrr; rV(i*3+1, j*3+2)=Hrz; 
            rV(i*3+2, j*3+0)=Hzp; rV(i*3+2, j*3+1)=Hzr; rV(i*3+2, j*3+2)=Hzz;
         }
      }
   }
   
   if (chainType == "pSheet") {
      for (i = 0; i < size/3; i++) {
         for (j = 0; j < size/3; j++) {
            Hxx = HC(i*3+0, j*3+0); Hxy = HC(i*3+0, j*3+1); 
            Hxz = HC(i*3+0, j*3+2);
            
            Hyx = HC(i*3+1, j*3+0); Hyy = HC(i*3+1, j*3+1); 
            Hyz = HC(i*3+1, j*3+2);
            
            Hzx = HC(i*3+2, j*3+0); Hzy = HC(i*3+2, j*3+1); 
            Hzz = HC(i*3+2, j*3+2);
         }
      }
            
   }
      
return rV;
}

SxMatrix<Double> SxSecondaryStructure::getFESHessian 
(const SxMatrix<Double> &HC)
{
   int i, j, size;
   bool mirroredI, mirroredJ;
   mirroredI = mirroredJ = false;
   size = nDoF;
   SxMatrix<Double> rV(size, size);
   double Hxx, Hyy, Hzz, Hxy, Hyx, Hxz, Hzx, Hyz, Hzy;
   
   //--- atoms per symmetry unit are equal to the number of 
   //    atoms per irreducible unit for FES
   //    beta-sheets contain a additional point group symmetry 

   if (chainType.contains("Sheet")) {
     cout << "TODO: POINT GROUP SYMMETRY OF ANTI PARALLEL BETASHEETS" << endl;
     //return HC;
   }
     
   for (i = 0; i < size/3; i++) {
            for (j = 0; j < size/3; j++) {

            Hxx = HC(i*3+0, j*3+0); Hxy = HC(i*3+0, j*3+1); 
            Hxz = HC(i*3+0, j*3+2);
            
            Hyx = HC(i*3+1, j*3+0); Hyy = HC(i*3+1, j*3+1); 
            Hyz = HC(i*3+1, j*3+2); 
            
            Hzx = HC(i*3+2, j*3+0); Hzy = HC(i*3+2, j*3+1); 
            Hzz = HC(i*3+2, j*3+2);

            if (chainType == "FES") {
               if ( ((i / nAtomsPP) % 2) == 0) mirroredI = false;
               else mirroredI = true;
               if ( ((j / nAtomsPP) % 2) == 0) mirroredJ = false;
               else mirroredJ = true;
            }
            
            if (chainType.contains ("Sheet"))  {
               if ( ((i / nAtomsPP) % 4) <= 1) mirroredI = false;
               else mirroredI = true;
               if ( ((j / nAtomsPP) % 4) <= 1) mirroredJ = false;
               else mirroredJ = true;
            }
               

            if ( (!mirroredI) && (mirroredJ))
            {
               Hxx = - Hxx; Hxy = -Hxy; Hxz = Hxz;
               Hyx = - Hyx; Hyy = -Hyy; Hyz = Hyz; 
               Hzx = - Hzx; Hzy = -Hzy; Hzz = Hzz;
            }

            if ( (mirroredI) && (mirroredJ) )
            {
               Hxx =   Hxx; Hxy =  Hxy; Hxz = -Hxz;
               Hyx =   Hyx; Hyy =  Hyy; Hyz = -Hyz; 
               Hzx = - Hzx; Hzy = -Hzy; Hzz =  Hzz;
            }
            
            if ( (mirroredI) && (!mirroredJ) )
            {
               Hxx =  -Hxx; Hxy = -Hxy; Hxz = -Hxz;
               Hyx =  -Hyx; Hyy = -Hyy; Hyz = -Hyz; 
               Hzx =   Hzx; Hzy =  Hzy; Hzz =  Hzz;
            }

            rV(i*3+0, j*3+0) = Hxx; rV(i*3+0, j*3+1) = Hxy; 
            rV(i*3+0, j*3+2) = Hxz;
            
            rV(i*3+1, j*3+0) = Hyx; rV(i*3+1, j*3+1) = Hyy; 
            rV(i*3+1, j*3+2) = Hyz;
            
            rV(i*3+2, j*3+0) = Hzx; rV(i*3+2, j*3+1) = Hzy; 
            rV(i*3+2, j*3+2) = Hzz;
         }
   }
return rV;
}
   
SxMatrix<Double> SxSecondaryStructure::getCartHessian 
    (const SxMatrix<Double> &HC, const SxAtomicStructure &t)
{
   int i, j;
   int size = (int)sqrt((double)HC.getSize ());
   SxString phix = SxString("phix");
   SxString phiy = SxString("phiy");
   SxString rx = SxString("rx");
   SxString ry = SxString("ry");
   SxMatrix<Double> rV(size, size);
   double Hxx, Hyy, Hzz, Hxy, Hyx, Hxz, Hzx, Hyz, Hzy;
   double Hpp, Hrr, Hpr, Hpz, Hrp, Hzp, Hrz, Hzr;

   for (i = 0; i < size/3; i++) {
         for (j = 0; j < size/3; j++) {

            Hpp = HC(i*3+0, j*3+0); Hpr = HC(i*3+0, j*3+1); 
            Hpz = HC(i*3+0, j*3+2);
            
            Hrp = HC(i*3+1, j*3+0); Hrr = HC(i*3+1, j*3+1); 
            Hrz = HC(i*3+1, j*3+2);
            
            Hzp = HC(i*3+2, j*3+0); Hzr = HC(i*3+2, j*3+1); 
            Hzz = HC(i*3+2, j*3+2);

            Hxx =   getT(phix, j, t)*getT(phix, i, t)*Hpp 
                  + getT(phix, j, t)*getT(rx, i, t)*Hrp    
                  + getT(rx, j, t)*getT(phix, i, t)*Hpr 
                  + getT(rx, j, t)*getT(rx, i, t)*Hrr;

            Hyy =   getT(phiy, j, t)*getT(phiy, i, t)*Hpp 
                  + getT(phiy, j, t)*getT(ry, i, t)*Hrp    
                  + getT(ry, j, t)*getT(phiy, i, t)*Hpr 
                  + getT(ry, j, t)*getT(ry, i, t)*Hrr;

            Hxy =   getT(phiy, j, t)*getT(phix, i, t)*Hpp 
                  + getT(phiy, j, t)*getT(rx, i, t)*Hrp    
                  + getT(ry, j, t)*getT(phix, i, t)*Hpr 
                  + getT(ry, j, t)*getT(rx, i, t)*Hrr;

            Hyx =   getT(phix, j, t)*getT(phiy, i, t)*Hpp 
                  + getT(phix, j, t)*getT(ry, i, t)*Hrp    
                  + getT(rx, j, t)*getT(phiy, i, t)*Hpr 
                  + getT(rx, j, t)*getT(ry, i, t)*Hrr;

         
            Hxz = getT(phix, i, t)*Hpz + getT(rx, i, t)*Hrz;
            Hyz = getT(phiy, i, t)*Hpz + getT(ry, i, t)*Hrz;
            Hzx = getT(phix, j, t)*Hzp + getT(rx, j, t)*Hzr;
            Hzy = getT(phiy, j, t)*Hzp + getT(ry, j, t)*Hzr;

            rV(i*3+0, j*3+0)=Hxx; rV(i*3+0, j*3+1)=Hxy; rV(i*3+0, j*3+2)=Hxz;
            rV(i*3+1, j*3+0)=Hyx; rV(i*3+1, j*3+1)=Hyy; rV(i*3+1, j*3+2)=Hyz; 
            rV(i*3+2, j*3+0)=Hzx; rV(i*3+2, j*3+1)=Hzy; rV(i*3+2, j*3+2)=Hzz;
         }
   }
return rV;
}

SxMatrix<Double> SxSecondaryStructure::getSymmetrizedMatrix 
(const SxMatrix<Double> &in) 
{
   int i, j, size;
   double avgValue;
   size = (int) sqrt((double) (in.getSize ()));
   SxMatrix<Double> returnValue (size, size);
   
   for (i = 0; i < size; i++) {
      for (j =0; j < size; j++) {
         avgValue = (in(i, j) + in(j, i))/2.;
         returnValue(i, j) = avgValue;
         returnValue(j, i) = avgValue;
      }
   }
   return returnValue;
}

double SxSecondaryStructure::getT
        (const SxString &whichT, int atomIndex, const SxAtomicStructure &tau) 
{
   double T;
   double phi, r, x, y/*, z*/;
   SxVector3<Double> tripel = tau(atomIndex);
   phi = tripel(0);r = tripel(1);//z = tripel(2);
   x = phi;y = r;

   T = 0;

   // should be switch/case
   if (whichT == SxString("xphi")) { 
      T = -sin(phi)*r; return T;
   }
   if (whichT == SxString("yphi")) {
      T = cos(phi)*r; return T;
   }
   if (whichT == SxString("xr")) {
      T = cos(phi); return T;
   }
   if (whichT == SxString("yr")) {
      T = sin(phi); return T;
   }
   if (whichT == SxString("rx")) {
      T = x/sqrt(x*x + y*y); return T;
   }
   if (whichT == SxString("ry")) {
      T = y/sqrt(x*x + y*y); return T;
   }
   if (whichT == SxString("phix")) { 
      T = -y/(x*x*(1. + y*y/(x*x))); return T;
   }
   if (whichT == SxString("phiy")) {
      T = 1./(x*(1. + y*y/(x*x))); return T;
   }
   
   cout << "Wrong coded String in getT" << endl;
   return 0;
}

SxArray<int> SxSecondaryStructure::getPeriodicIndices (int atomIndex ) 
{
   int i, j;
   bool found =false;
   SxArray<int> returnValue(2);
   
   for (i = 0; i < nAtomsPP; i++) {
      for (j = 0; j < nPeptides; j++) {
         if (((int)periodicity(i)(j) - 1) == atomIndex) {
            returnValue(0) = i;
            returnValue(1) = j;
            i = nAtomsPP;
            j = nPeptides;
            found = true;   
         }
      }
   }
   // cout << "ret " << returnValue(0) << " " <<  returnValue(1) << endl;
   if (found == false) {
      cout << "Periodic Index not found in Periodicity File. Atom Index: "
         << atomIndex << endl;
      SX_EXIT;
   }

   //--- in case of beta-sheets the irreducible unit contains 4 peptide units
   if (chainType.contains ("Sheet")) {
      returnValue(0) += (returnValue(1) % 4) * nAtomsPP;
      returnValue(1) = returnValue(1)/4;
   }

   return returnValue;
}

double SxSecondaryStructure::getPDRes () 
{
   int ggT = 0;
   if (((int)nPeptides % (int)nTurns) == 0) {
      ggT = (int) nTurns;
   }

   int nPoints = nPeptides;
   int fac = 2;
   while (nPoints < 500) {
      if (fac != ggT) nPoints*=fac;
      fac++;
   }

   //cout << "nKPoints hard coded !!!!!!!!!: 2*990"<< endl;
   //nPoints = 2*990;
   return (360./((double)(nPoints)));
}

SxAtomicStructure SxSecondaryStructure::getCoords () 
{
   switchToInputOrder (); 
   return tauCart;
}
    
void SxSecondaryStructure::resizeKMesh (int nK)
{
   if (((nK % nPeptides) != 0) || (nK < nPeptides)) {
      cout << "Ill-defined nKPoints: must be an integer multiple of " 
           << "the number of peptides " << endl; 
           SX_QUIT;
   }
   nKPoints = nK;
   computeKMesh ();
}

void SxSecondaryStructure::computeKMesh ()
{
   double KPoint;
   int i;
   SxArray<SxMatrix<Complex16> > kFoldDynamicals ((int) nKPoints);
   SxHessianOps hOps;
   
   switchToPeptideOrder ();

   SxVector<Double> massVecPP (nDoFPU);
   for (i = 0; i < nDoFPU; i++) massVecPP(i) = masses.coordRef()(i);
   
   for (i = 0; i < nKPoints; i++) {
      kFoldDynamicals(i) = SxMatrix<Complex16> (nDoFPU, nDoFPU);
      kFoldDynamicals(i).set (0.);
   }

      cout << "Dispersion Relation in Cyl. Coordiantes !" << endl;
      cout << "Calculating the reduced Hessians ... start" << endl;   
      
      for (int K = 0; K < nKPoints; K++) {
         if ((K % 5) == 0) cout << ".";
         KPoint =  360./(double)nKPoints*(double)K;
         kFoldDynamicals(K).set (getK (KPoint));
      }

      kMeshFreqs = SxArray<SxVector<Complex16> > (nKPoints);
      for (i = 0; i < nKPoints; i++) {
         kMeshFreqs(i) = SxVector<Double>(nDoFPU);
         hOps = SxHessianOps(kFoldDynamicals(i), massVecPP);
         kMeshFreqs(i).set (hOps.getEigenFrequencies ());
      }
}

SxMatrix<Complex16> SxSecondaryStructure::getK (double phaseInDeg)
{
   double weight = 1.;
   double phase;
   int i;
   SxMatrix<Complex16> dynCyl (nDoFPU, nDoFPU);
   SxMatrix<Complex16> returnValue (nDoFPU, nDoFPU);
   SxMatrix<Complex16> block (nDoFPU, nDoFPU);

   phase = phaseInDeg/360.*2.*PI;
   
   dynCyl.set (0.);
   returnValue.set (0.);
   switchToPeptideOrder ();
      
   for (i = -(nUnits)/2 ; i <= nUnits/2; i++) {
      block.set (0.);
      if (i >= 0) block.set (blockArrayCyl(0)(i));
      else    block.set (blockArrayCyl(0)(nPeptides + i)) ;   
         
         if (((nUnits % 2) == 0) && (abs(i) == nUnits/2)) 
            weight = 0.5;
         else weight = 1.0;
                
         dynCyl = dynCyl + block
            *    ((cos(phase * (double)i) + I*sin(phase * (double)i))) 
            * weight;
   }


   if (chainType == ("helix"))   
      returnValue = B.transpose () ^ ( (dynCyl^B));
  
   return returnValue; 
}

void SxSecondaryStructure::computeNormalDisplacements ()
{
   int nDoFPIU = getNDoFPIU ();
   SxMatrix<Double> reducedHessian (nDoFPIU, nDoFPIU);
   SxMatrix<Double> gammaHessian (nAtomsPU * 3, nAtomsPU *3);
   SxVector<Double> column(nDoF);
   SxVector<Double> orth(nDoF);
   SxVector<Double> workAroundPsi(nDoF);
   SxMatrix<Double>::Eigensystem eig;
 
   normalDisplacements = SxMatrix<Double> (nDoF, nDoF);
   normalCurvatures = SxVector<Double> (nDoF);

   normalDisplacements.set (0.); normalCurvatures.set (0.);
   
   gammaHessian = getKcart(0);
   for (int i = 0; i < nDoFPIU; i++) {
      for (int j = 0; j < nDoFPIU; j++) {
         reducedHessian(i, j) = gammaHessian (i, j);
      }
   }
   
   eig = SxMatrix<Double>::Eigensystem ();
   eig = reducedHessian.eigensystem ();


   workAroundPsi.resize(nDoFPIU);
   for (int i = 0; i < nDoFPIU; i++) {
      workAroundPsi.set(eig.vecs.colRef(i));
      if (workAroundPsi(0) <= 0.) 
         workAroundPsi = -workAroundPsi;
      eig.vecs.colRef(i).set(workAroundPsi);
   }

   normalCurvatures.set (0.);
   // --- refinement calculation is based on eigenvalues of hessians,
   //     these must be limited
   for (int i = 0; i < nDoFPIU; i++) {
      if (fabs(eig.vals(i).re) < 0.001) normalCurvatures(i) = 0.001;
      else normalCurvatures(i) = fabs(eig.vals(i).re);
   }

   for (int i = 0; i < nDoF; i++) {
            
      if (i < nDoFPIU) {
         column.set(0.);
         for (int j = 0; j < nDoFPIU; j++) 
            column(j) = eig.vecs.colRef(i)(j);
              
         //--- orthonormalisation
         column = column/sqrt(column.absSqr().sum());
         for (int j = 0; j < i; j++) {
            orth.set (normalDisplacements.colRef(j));
            column = column - (column^orth).chop()*orth;
            column = column/sqrt(column.absSqr().sum());
         }
               normalDisplacements.colRef(i).set(column);
      } 
            else 
               normalDisplacements.colRef(i).set(0.);
         }
      
         // --- the deviation vectors are normed according 
         //     to the reciprocal curvatures of the PES
         //     this results in more accurate results 
         //     for the low frequency branches
      
         for (int i = 0; i < nDoFPIU; i++) {
            column.set (normalDisplacements.colRef(i));
            column = column / sqrt(column.absSqr().sum())
                            / sqrt(fabs(eig.vals(i).re));
            normalDisplacements.colRef(i).set (getInputOrderedVec(column));
         }
}

SxVector<Double> SxSecondaryStructure::getDisplacementVector (int DoF) 
{
   SxVector<Double> returnValue (nDoF);
   returnValue.set (normalDisplacements.colRef (DoF));
   return returnValue;
}

SxVector<Double> SxSecondaryStructure::getNormedDisplacementVector (int DoF) 
{
   SxVector<Double> returnValue(nDoF);
   double norm;
   returnValue.set(normalDisplacements.colRef(DoF));
   norm = sqrt(returnValue.absSqr().sum());
   if (norm > 1e-15)
      returnValue = returnValue/norm;
   return returnValue;
}

SxMatrix<Double> SxSecondaryStructure::getKcart (int i)
{
   SxArray<SxArray<SxMatrix<Double> > > blockArray (nPeptides);
   SxMatrix<Double> returnValue = SxMatrix<Double> (nDoFPP, nDoFPP);

   //--- a bit inefficient but not a bottleneck yet
   blockArray = getPeptideOrderedBlocks (hessianCart);
   returnValue.copy (blockArray(0)(i));
   return returnValue;
}
   
   
SxComplex16 SxSecondaryStructure::getFrequency 
                      (int branch, double phase, const SxString &interpolation) 
{
   int i;

   if (interpolation == ("spline")) {
      // ugly: splines are recalculated for each phi
      // but: is not a bottleneck
      
      SxArray<double> x (nKPoints + 1);
      SxArray<double> y (nKPoints + 1);
      SxArray<double> y2 (nKPoints + 1);
      double dy, dx;
         
      for (i = 0; i <= nKPoints; i++) {
         if (i == nKPoints) 
            y (i) = kMeshFreqs(0)(branch);
         else y(i) = kMeshFreqs(i)(branch);

         x(i) = 360./(double)(nKPoints)*(double)i;

         //if (y(i) > 10.) y(i) += 10;
         //if (j == 0)  cout << x(i) << " " << y(i) << endl;
      }
      if (branch <= 1) {
         dy = kMeshFreqs(1)(branch);
         dx = 360./nKPoints;
         if (branch == 0) {
            if ((phase < twistDeg) 
                  || (phase > (360.-twistDeg)) 
                  || (fabs(twistDeg - 180.) < 1e-2)) { 
               SxArray<double> x0 ((int)nTurns + 1);
               SxArray<double> y0 ((int)nTurns + 1);
               SxArray<double> y20 ((int)nTurns + 1);
               double shiftedPhase = phase;
               if (shiftedPhase > twistDeg) 
                  shiftedPhase = -(shiftedPhase - 360.);
               for (i = 0; i <= nTurns; i++) {
                  x0(i) = x(i);
                  y0(i) = y(i);
                  //cout << "x " << x(i) << endl;
                  //cout << "y " << y(i) << endl;
               } 
               spline (x0, y0, (int)y0.getSize (), dy/dx, 0, &y20);
               return  splint (x0, y0, y20, (int)x0.getSize (), shiftedPhase);
            } else {
               int nK = nKPoints - 2*(int)nTurns + 1;
               int offset = (int)nTurns;
               SxArray<double> x0 (nK);
               SxArray<double> y0 (nK);
               SxArray<double> y20 (nK);
               for (i = 0; i < nK; i++) {
                  x0(i) = x(i + offset);
                  y0(i) = y(i + offset);
               } 
               spline (x0, y0, (int)y0.getSize (), 0, 0, &y20);
               return  splint (x0, y0, y20, (int)x0.getSize (), phase);
            }
         }
            
         if (branch == 1) 
            spline (x, y, (int)y.getSize (), dy/dx, -dy/dx, &y2);

      } else spline (x, y, (int)y.getSize (), 0., 0., &y2);
      return  splint (x, y, y2, (int)x.getSize (), phase);
   }
  
   if (interpolation == SxString("linear")) {
      double kSpacing = 360./(double) nKPoints;
      int kIndex = (int)(phase/kSpacing);
      int k1 = kIndex;
      if (k1 == nKPoints) k1 = 0;
      int k2 = kIndex + 1;
      if (k2 >= nKPoints) k2 = 0;

      double dY = kMeshFreqs(k2)(branch) - kMeshFreqs(k1)(branch);
     
      return (dY * (phase - (double)kIndex*kSpacing)/
               kSpacing + kMeshFreqs(k1)(branch));
   }

   if (interpolation == SxString("fourier")) {
      int nKPointsInterpol = (int)(360./getPDRes());
      if (nKPoints != nKPointsInterpol) {
         nKPoints = nKPointsInterpol;
         computeKMesh ();
      }
      return  (kMeshFreqs((int)(phase/getPDRes() + 0.5))(branch));
   } 
      
   if (interpolation == ("splineForLowest")) {
      /*
*/
         //if (y(i) > 10.) y(i) += 10;
         //if (j == 0)  cout << x(i) << " " << y(i) << endl;
  //    }
      if (branch >= 1) {

         //--- in case of systems where the second branch (with index 1)
         //    gets zero at 180 degrees (i.e. e.g for FES)
         //    a spline interpolation close to 180 degrees
         //    for all other cases (except the lowest branch) the fourier
         //    interpolation is better 
         if( (fabs((twistDeg - 180.0)) < 1e-2) 
            && (branch == 1) 
            && (fabs(180.0 - phase) <  1./(double)nPeptides * 360.0))  {
            
            if (nKPoints !=  nPeptides) { 
               nKPoints = nPeptides;
               computeKMesh ();
            }
            return  getFrequency (1, phase, SxString ("spline"));
                  
         }
         else  {
            int nKPointsInterpol = (int)(360./getPDRes());
      
            if (nKPoints != nKPointsInterpol) {
               nKPoints = nKPointsInterpol;
               computeKMesh ();
            }
            return  (kMeshFreqs((int)(phase/getPDRes() + 0.5))(branch));
         }
      }
    
      else {
         //--- stable estimate for 1st derivative at gamma point
         //    dy/dx = y(10.0 deg)/10.0 deg
         //   for smaller delta x's there's the risk of numerical
         //   noise
         int nKPointsInterpol = (int)(360./getPDRes());
         if (nKPoints != nKPointsInterpol) {
            nKPoints = nKPointsInterpol;
            computeKMesh ();
         }
      
         SxArray<double> x (nKPoints + 1);
         SxArray<double> y (nKPoints + 1);
         SxArray<double> y2 (nKPoints + 1);
         
         for (i = 0; i <= nKPoints; i++) {
            if (i == nKPoints) 
               y (i) = kMeshFreqs(0)(branch);
            else y(i) = kMeshFreqs(i)(branch);
            x(i) = 360./(double)(nKPoints)*(double)i;
         } 

         int kLinear = (int) (10.0/getPDRes ()); 
         double dy = kMeshFreqs (kLinear)(0);
         double dx = (double) kLinear*getPDRes ();
         
         //--- estimate for frequency at high symmetry point (180 deg)
         //    improves the spline interpolation for supercell with odd 
         //    number of peptides         
         double antiSymmFreq = kMeshFreqs(nKPoints/2)(0);
      
         //--- spline interpolation used for the lowest branch
         //    boundary conditions are taken partially from fourier
         //    interpolation
         int nKbyNPep = nKPoints/nPeptides;
            
         if ((phase < twistDeg) 
               || (phase > (360.-twistDeg)) 
               || (fabs(twistDeg - 180.) < 1e-2)) { 
            SxArray<double> x0 ((int)nTurns + 1);
            SxArray<double> y0 ((int)nTurns + 1);
            SxArray<double> y20 ((int)nTurns + 1);
            double shiftedPhase = phase;
            if (shiftedPhase > twistDeg) shiftedPhase = -(shiftedPhase - 360.);
            for (i = 0; i <= nTurns; i++) {
               x0(i) = x(i*nKbyNPep);
               y0(i) = y(i*nKbyNPep);
               //cout << "x " << x(i) << endl;
               //cout << "y " << y(i) << endl;
            } 
            spline (x0, y0, (int)y0.getSize (), dy/dx, 0, &y20);
            return  splint (x0, y0, y20, (int)x0.getSize (), shiftedPhase);
         } else {
            int nK = nPeptides - 2*(int)nTurns + 1;
            int offset = (int)nTurns*nKbyNPep;
            SxArray<double> x0 (nK);
            SxArray<double> y0 (nK);
            SxArray<double> y20 (nK);
            
            for (i = 0; i < nK; i++) {
               x0(i) = x(i*nKbyNPep + offset);
               y0(i) = y(i*nKbyNPep + offset);
            } 

            if ((nPeptides % 2) == 1) {
               SxArray<double> xExtended (nK + 1);
               SxArray<double> yExtended (nK + 1);
               SxArray<double> y2Extended (nK + 1);
               for (i = 0; i < nK/2; i++) {
                  xExtended(i) = x0(i);
                  yExtended(i) = y0(i);
               } 
               xExtended(nK/2) = 180.0;
               yExtended(nK/2) = antiSymmFreq;
               for (i = (nK/2 + 1); i <= nK; i++) {
                  xExtended(i) = x0(i - 1);
                  yExtended(i) = y0(i - 1);
               } 
               x0  = SxArray<double> (nK + 1);
               y0  = SxArray<double> (nK + 1);
               y20 = SxArray<double> (nK + 1);
               
               for (i = 0; i <= nK; i++) {
                  x0(i) = xExtended(i);
                  y0(i) = yExtended(i);
               }
            }
            
            spline (x0, y0, (int)y0.getSize (), 0, 0, &y20);
            return  splint (x0, y0, y20, (int)x0.getSize (), phase);
         }
      }
   }   
   return 0.;
}

SxArray<double> SxSecondaryStructure::
splineItLowestBranch (const SxArray<double> &in)
{
   
   //--- stable estimate for 1st derivative at gamma point
   //    dy/dx = y(10.0 deg)/10.0 deg
   //   for smaller delta x's there's the risk of numerical
   //   noise
   int nK = (int)in.getSize ();
   double scale = (double)(nK)/180.;
   SxArray<double> x (nK);
   SxArray<double> y (nK);
   SxArray<double> returnValue (nK);
         
   for (int i = 0; i < nK; i++) {
         y (i) = in(i);
         x(i) = (double)i/scale;
   } 

   int kLinear = (int) (10.0*scale); 
   double dy = in(kLinear);
   double dx = (double) 10.0;
         
   //--- estimate for frequency at high symmetry point (180 deg)
   //    improves the spline interpolation for supercell with odd 
   //    number of peptides         
   double antiSymmFreq = in(nK-1);
      
   //--- spline interpolation used for the lowest branch
   //    boundary conditions are taken partially from fourier
   //    interpolation
   int nKbyNPep = nK/nPeptides*2;
   //int twistIndex = (int)(twistDeg*scale);
   int splinePoints = (int)nTurns + 1;
   
   SxArray<double> x0ST (splinePoints);
   SxArray<double> y0ST (splinePoints);
   SxArray<double> splineST (splinePoints);

   //--- generating spline coefficients for phase angles < twist angle
   
   cout << "FIRST HALF "<< endl;
   for (int i = 0; i < splinePoints; i++) {
      x0ST(i) = x(i*nKbyNPep);
      cout << "X" << x0ST(i) << endl;
      y0ST(i) = y(i*nKbyNPep);
      cout << "Y" << y0ST(i) << endl;
   } 

   spline (x0ST, y0ST, (int)y0ST.getSize (), dy/dx, 0, &splineST);
           
   //--- generating spline coefficients for phase angles > twist angle
   
   splinePoints = (nPeptides - 2*(int)nTurns + 2)/2;
   int offset = (int)nTurns*nKbyNPep;
   SxArray<double> x0LT (splinePoints);
   SxArray<double> y0LT (splinePoints);
   SxArray<double> splineLT (splinePoints);
   
   cout << "SECOND HALF "<< endl;
   for (int i = 0; i < splinePoints; i++) {
      x0LT(i) = x(i*nKbyNPep + offset);
      cout << "X" << x0LT(i) << endl;
      y0LT(i) = y(i*nKbyNPep + offset);
      cout << "Y" << y0LT(i) << endl;
   } 
   
   if ((nPeptides % 2) == 1) {
      SxArray<double> xExtended (splinePoints + 1);
      SxArray<double> yExtended (splinePoints + 1);
      SxArray<double> y2Extended (splinePoints + 1);
      for (int i = 0; i < splinePoints; i++) {
         xExtended(i) = x0LT(i);
         yExtended(i) = y0LT(i);
      } 
      xExtended(splinePoints) = 180.0;
      yExtended(splinePoints) = antiSymmFreq;
      
      x0LT  = SxArray<double> (splinePoints + 1);
      y0LT  = SxArray<double> (splinePoints + 1);
      splineLT = SxArray<double> (splinePoints + 1);
               
      for (int i = 0; i <= splinePoints; i++) {
         x0LT(i) = xExtended(i);
         y0LT(i) = yExtended(i);
      }
   }
            
   spline (x0LT, y0LT, (int)y0LT.getSize (), 0, 0, &splineLT);

   for (int i = 0; i < nK; i++) {
      if (x(i) <= twistDeg) {
         returnValue (i) = 
           splint (x0ST, y0ST, splineST, (int)x0ST.getSize (), x(i));
      } else
         returnValue (i) = 
           splint (x0LT, y0LT, splineLT, (int)x0LT.getSize (), x(i));
   }

   //-- for pi-helix it is more convenient to take the region 120 
   //   degs to 240 degs from fourier interpolation
   
   if (fabs ((twistDeg - 80.0) ) < 5.0) {
      for (int i = nK - 1; fabs (180. - (double)i/scale) < 60.0 ; i--) {
         returnValue(i) = in (i);
      }
   }
   
   return returnValue;
}

SxVector<Complex16> 
SxSecondaryStructure::getFrequencies 
(const SxString &interpolation, double dphi) 
{
   int size = 0;
   //double phase;
   SxVector<Complex16> returnValue;
   int counter = 0;
   int i =0;
   
   // --- ugly
   if (interpolation == "none") {
      size = nKPoints*nDoFPP;
      returnValue.resize (size);
      for (i = 0; i < nDoFPP; i++) {
         for (int phaseId = 0; phaseId < nKPoints; phaseId++) {
            //phase = (double)phaseId * 360. /(double)nKPoints;
            returnValue(i*nKPoints + phaseId) = kMeshFreqs(phaseId)(i);
         }
      }
   } else {
      size = (int) ((360./dphi) + .5) * nDoFPP;
      returnValue.resize (size);
      for (i = 0; i < nDoFPP; i++) {
         for (counter = 0; counter < ((int) ((360./dphi) + .5)); counter++) { 
            returnValue (i*((int) ((360./dphi) + .5)) + counter) 
               = getFrequency(i, dphi*(double)counter, interpolation);
           
            // --- important change (2.9.2005)
            //     for the region 120 degs < phi < 180 degs 
            //     the fourier interpolation is pursuingly 
            //     applied for each branch
            //     to also account for oscillations in this region 
            //     of the first vibrational branch in the ces of pi-helix
            
            if ((fabs(180. - counter*dphi) < 60.) && (twistDeg < 90)) {
             returnValue (i*((int) ((360./dphi) + .5)) + counter) 
               = getFrequency(i, dphi*(double)counter, "fourier");
            }
  
         }
      }
   }
   return returnValue;
}

int SxSecondaryStructure::getDoFInputOrder (int DoFPeptideOrder) 
{
  int nA = DoFPeptideOrder/3;
  int nXYZ = DoFPeptideOrder - 3*nA;
  int nP = nA/nAtomsPP;
  int returnValue = (int)(periodicity(nA - nP*nAtomsPP)(nP) - 1)*3 + nXYZ;
  return returnValue;
}
  
void SxSecondaryStructure::printPhononDispersionRelation 
         (const SxString &fileName, const SxString &interpolation, double /*dphi*/)
{
   int i;
   double phase;
   SxComplex16 freq;
   
   FILE *fp = NULL;
   if ( !(fp = fopen (fileName.ascii (), "w")) ) {
      sxprintf ("Can't open file %s",fileName.ascii ()); 
      SX_EXIT;
   }

   // --- ugly
   if (interpolation == "none") {
      for (i = 0; i < nDoFPP; i++) {
         for (int phaseId = 0; phaseId <= nKPoints; phaseId++) {
            phase = (double)phaseId * 360. /(double)nKPoints;
            if (phaseId ==  nKPoints) 
               fprintf (fp, "%f %f\n", phase, kMeshFreqs(0)(i).re);
            else 
               fprintf (fp, "%f %f\n", phase, kMeshFreqs(phaseId)(i).re);
            
         }
         fprintf (fp, "\n");
      }
   } else {
      for (i = 0; i < nDoFPP; i++) {

         for (int K = 0; K < (int)(360./getPDRes ()); K++) {  
            phase = getPDRes () * (double) K;
            freq   
               = getFrequency(i, phase, interpolation); 

            
            // --- important change (2.9.2005)
            //     for the region 120 degs < phi < 180 degs 
            //     the fourier interpolation is pursuingly 
            //     applied for each branch
            //     to also account for oscillations in this region 
            //     of the first vibrational branch in the case of pi-helix
            
            if ((fabs(180. - phase) < 60.) && (twistDeg < 90.)) 
            freq  
              = getFrequency (i, phase, "fourier");

            
            
            if (freq.im > freq.re) freq.re = -freq.im;
            fprintf 
               (fp, "%f %f\n", phase, freq.re);
         }
            fprintf (fp, "\n");
      }
   }
   fclose (fp);
}

SxComplex16 SxSecondaryStructure::getEValue (double phaseInDeg, const SxVector<Complex16> &mode, const SxArray<SxArray<SxMatrix<Double> > > &blocks)
{
   double phase;
   int iPep;
   SxMatrix<Complex16> block (nDoFPU, nDoFPU);
   SxMatrix<Complex16> dynamical (nDoFPU, nDoFPU);
   switchToPeptideOrder ();
   SxVector<Double> massVecPP (nDoFPU);
   for (int i = 0; i < nDoFPU; i++) massVecPP(i) = masses.coordRef()(i);

   SxComplex16 w = 0.;
   for (int i = 0; i < nPeptides; i++) {
      if (i > nPeptides/2) iPep = i - nPeptides;
      else iPep = i;
      phase = phaseInDeg*DEG2RAD*(double)iPep;
      block.copy (blocks(0)(i));
                
      block =  block
         *    ((cos(phase ) + I*sin(phase ))); 

      block = B.transpose () ^ ( (block^B));
   
      for (int k = 0; k < nDoFPU; k++) {
         for (int l = 0; l < nDoFPU; l++) {
            dynamical(k, l) = block(k, l) 
               /sqrt (massVecPP(k))/ sqrt(massVecPP(l)); 
         }
      }
     // SxHessianOps hessianOps (block, massVecPP);
      SxComplex16 l 
         = AU2CM*AU2CM*(mode.adjoint ()^dynamical^mode).chop ();
      w += l;
     // w += (mode.adjoint ()^dynamical^mode).chop ();
   }
   return (sqrt(w));
}
 
SxComplex16 SxSecondaryStructure::getIaCoeff (double phaseInDeg, const SxVector<Complex16> &mode, int iPepIn)
{
   double phase;
   int iPep;
   SxMatrix<Complex16> block (nDoFPU, nDoFPU);
   switchToPeptideOrder ();
   SxVector<Double> massVecPP (nDoFPU);
   for (int i = 0; i < nDoFPU; i++) massVecPP(i) = masses.coordRef()(i);

   if (iPepIn > nPeptides/2) iPep = iPepIn - nPeptides;
   else iPep = iPepIn;

   phase = phaseInDeg*DEG2RAD*(double)iPep;
   block.copy (blockArrayCyl(0)(iPepIn));
                
   
   block =  block
      *    ((cos(phase ) + I*sin(phase ))); 


   if (chainType == ("helix")) 
      block = B.transpose () ^ ( (block^B));
   
   SxHessianOps hessianOps (block, massVecPP);
/*
   if ((iPep == 2) && (int (phaseInDeg + 0.5) == 60)) {
      cout << mode << endl;
      cout << AU2CM*AU2CM*(mode.adjoint ()^hessianOps.dynamical^mode).chop () 
         << endl;
      SX_EXIT;
   }
   */
   return (AU2CM*AU2CM*(mode.adjoint ()^hessianOps.dynamical^mode).chop ());

}

SxComplex16 SxSecondaryStructure::getIaCoeffdWdPhi (double phaseInDeg, const SxVector<Complex16> &mode, int iPepIn, double threshold)
{
   double phase;
   int iPep;
   SxMatrix<Complex16> block (nDoFPU, nDoFPU);
   switchToPeptideOrder ();
   int counter = 0.;
   SxVector<Double> massVecPP (nDoFPU);
   for (int i = 0; i < nDoFPU; i++) massVecPP(i) = masses.coordRef()(i);

   if (iPepIn > nPeptides/2) iPep = iPepIn - nPeptides;
   else iPep = iPepIn;

   phase = phaseInDeg*DEG2RAD*(double)iPep;
   block.copy (blockArrayCyl(0)(iPepIn));
                
 /*
   block =  block
      *    ((I*cos(phase ) - sin(phase ))) 
      * (double)iPep;
*/
  
  // cout << "iPep = " << iPepIn; 
   for (int k = 0; k < nDoFPP; k++) {
      for (int l = 0; l < nDoFPP; l++) {
            if ((block (k, l).abs()  <= threshold)) {
               block (k , l) = 0.;
               counter ++;
            }
      }
   }
   
  // cout << " " << counter << endl;
         
   if ((iPep == 1) || (iPep == -1))
      block = 1.0*block;
      
   if (chainType == ("helix")) 
      block = B.transpose () ^ ( (block^B));
   
   SxHessianOps hessianOps (block, massVecPP);
   return 
      (0.5*AU2CM*AU2CM*DEG2RAD
       *(double)iPep*((I*cos(phase ) - sin(phase ))) 
          *(mode.adjoint ()^hessianOps.dynamical^mode).chop ());
}

SxComplex16 SxSecondaryStructure::getIaCoeffdWdPhiThinned (double phaseInDeg, const SxVector<Complex16> &mode, const SxMatrix<Double> &im, int iPepIn, double threshold)
{
   double phase;
   int iPep;
   SxMatrix<Complex16> block (nDoFPU, nDoFPU);
   switchToPeptideOrder ();
   int counter = 0.;
   SxVector<Double> massVecPP (nDoFPU);
   for (int i = 0; i < nDoFPU; i++) massVecPP(i) = masses.coordRef()(i);

   if (iPepIn > nPeptides/2) iPep = iPepIn - nPeptides;
   else iPep = iPepIn;

   phase = phaseInDeg*DEG2RAD*(double)iPep;
   block.copy (blockArrayCyl(0)(iPepIn));
                
 /*
   block =  block
      *    ((I*cos(phase ) - sin(phase ))) 
      * (double)iPep;
*/
  
   cout << "iPep = " << iPepIn; 
   for (int k = 0; k < nDoFPP; k++) {
      for (int l = 0; l < nDoFPP; l++) {
            if (im (k, l) <= threshold) {
               block (k , l) = 0.;
               counter ++;
            }
      }
   }
   
   cout << " " << counter << endl;
         
   if ((iPep == 1) || (iPep == -1))
      block = 1.0*block;
      
   if (chainType == ("helix")) 
      block = B.transpose () ^ ( (block^B));
   
   SxHessianOps hessianOps (block, massVecPP);
   return 
      (0.5*AU2CM*AU2CM*DEG2RAD
       *(double)iPep*((I*cos(phase ) - sin(phase ))) 
          *(mode.adjoint ()^hessianOps.dynamical^mode).chop ());
}

SxMatrix<Double> SxSecondaryStructure::getImportanceMatrix (double phaseInDeg, const SxVector<Complex16> &mode, int iPepIn)
{
   double phase;
   int iPep;
   double correctValue, eraseValue;
   SxMatrix<Complex16> block (nDoFPU, nDoFPU);
   SxMatrix<Complex16> tBlock (nDoFPU, nDoFPU);
   SxMatrix<Double> returnValue (nDoFPU, nDoFPU);
   switchToPeptideOrder ();
   SxVector<Double> massVecPP (nDoFPU);
   for (int i = 0; i < nDoFPU; i++) massVecPP(i) = masses.coordRef()(i);

   if (iPepIn > nPeptides/2) iPep = iPepIn - nPeptides;
   else iPep = iPepIn;

   phase = phaseInDeg*DEG2RAD*(double)iPep;
                
 /*
   block =  block
      *    ((I*cos(phase ) - sin(phase ))) 
      * (double)iPep;
*/
  
   block.copy (blockArrayCyl(0)(iPepIn));
   if (chainType == ("helix")) 
      block = B.transpose () ^ ( (block^B));
   SxHessianOps hessianOps (block, massVecPP);
   cout << "iPep = " << iPepIn; 
   correctValue =  
      (0.5*AU2CM*AU2CM*DEG2RAD
           *(double)iPep*((I*cos(phase ) - sin(phase ))) 
             *(mode.adjoint ()^hessianOps.dynamical^mode).chop ());
      
   
   for (int k = 0; k < nDoFPP/3; k++) {
      cout << "Erasing ..." << k << endl;
      for (int l = 0; l < nDoFPP/3; l++) {
            
            block.copy (blockArrayCyl(0)(iPepIn));
           
            for (int i = 0; i < 3; i++) {
               for (int j = 0; j < 3; j++) {
                  block (3*k + i , 3*l + j) = 0.;
               }
            }
            
            if (chainType == ("helix")) 
               block = B.transpose () ^ ( (block^B));
        
            SxHessianOps hOps (block, massVecPP);
            eraseValue = 
               (0.5*AU2CM*AU2CM*DEG2RAD
                *(double)iPep*((I*cos(phase ) - sin(phase ))) 
                *(mode.adjoint ()^hOps.dynamical^mode).chop ());

            for (int i = 0; i < 3; i++) {
               for (int j = 0; j < 3; j++) {
                  returnValue (3*k + i, 3*l + j) 
                     = fabs((eraseValue - correctValue));
               }
            }
      }
   }
/*
   for (int k = 0; k < nDoFPP; k++) {
      cout << "Erasing ..." << k << endl;
      for (int l = 0; l < nDoFPP; l++) {
            
            block.copy (blockArrayCyl(0)(iPepIn));
           
            block (k , l) = 0.;
            
            if (chainType == ("helix")) 
               block = B.transpose () ^ ( (block^B));
        
            SxHessianOps hOps (block, massVecPP);
            eraseValue = 
               (0.5*AU2CM*AU2CM*DEG2RAD
                *(double)iPep*((I*cos(phase ) - sin(phase ))) 
                *(mode.adjoint ()^hOps.dynamical^mode).chop ());

               returnValue (k, l) = fabs((eraseValue - correctValue));
               //cout << correctValue << " " << eraseValue << " " << returnValue(k, l) << endl;
      }
   }
  */    
   return returnValue;
      
}

SxMatrix<Double> SxSecondaryStructure::getImportanceMatrixFreq
(int /* branchIndex */, int /* iPep */)
{
   //--- store old k-mesh setting and set k-mesh size to nPeptides
   //int nKPointsOld = nKPoints;
   resizeKMesh (nPeptides);
   
   //--- store old cartesian hessian 
   SxMatrix<Double> hessianCartPOrder (nDoF, nDoF);
   hessianCartPOrder = getHessianFromBlocks (blockArrayCart);
   SxMatrix<Double> hessianDel (nDoF, nDoF);
   
   hessianDel.copy (hessianCartPOrder);
   hessianCyl = getCylHessian (hessianDel, tauCyl);
   blockArrayCyl = getBlocks (hessianCyl);

   computeKMesh ();

   for (int i  = 0; i < nDoFPP; i++) {
     for (int j  = 0; j < nKPoints; j++) {
        cout << kMeshFreqs(j)(i) << endl;
     }
   }
        
  SX_EXIT;
  return SxMatrix<Double> ();
}

int SxSecondaryStructure::
    getInChainDistance (int iPep, int jPep, int iAtoms, int jAtoms)
{
          int n1 = iPep*nAtomsPP + iAtoms;
          int n2 = jPep*nAtomsPP + jAtoms;
          int dist1 = abs (n1 - n2);
          n2 = (jPep - nPeptides)*nAtomsPP + jAtoms;
          int dist2 = abs (n1 - n2);
          if (dist2 < dist1 ) 
          dist1 = dist2;
          n2 = (jPep + nPeptides)*nAtomsPP + jAtoms;
          dist2 = abs (n1 - n2);
          if (dist2 < dist1 ) 
          dist1 = dist2;
         return dist1;
}

bool SxSecondaryStructure::inShell 
(int iPep, int jPep, int iAtoms, int jAtoms, int rad) 
{
    if ((getInChainDistance (iPep, jPep, iAtoms, jAtoms) 
                            <= rad ) 
                          ) 
    return true;
    else return false;
}

void SxSecondaryStructure::rescaleForceConstants
(int modus, double scaling, int rad)
{
   
   //--- store old cartesian hessian 
   SxMatrix<Double> hessianCartPOrder (nDoF, nDoF);
   hessianCartPOrder = getHessianFromBlocks (blockArrayCart);
   
   //-- hessian to be manipulated
   SxMatrix<Double> hessianDel (nDoF, nDoF);
   hessianDel.copy (hessianCartPOrder) ;
   int ii, jj;
  
   //-- some control output
   /*
   double sum = 0.;
   for (int x = 0; x < 3; x++) {
      for (int y = 0; y < 3; y++) {
      
         for (int iPep = 0; iPep < nPeptides; iPep++) {
            for (int iAtoms = 0; iAtoms < nAtomsPP; iAtoms++) {
               sum = 0.;
               for (int jPep = 0; jPep < nPeptides; jPep++) {
                  for (int jAtoms = 0; jAtoms < nAtomsPP; jAtoms++) {
               
                     ii = (iPep*nAtomsPP + iAtoms)*3 + x;
                     jj = (jPep*nAtomsPP + jAtoms)*3 + y;

                     sum += hessianDel(ii, jj);
                  }
               }
               cout << sum << endl;
            }
         }
      }
   }

   cout << "3x3" << endl;
   cout << hessianDel(120, 120) << " " << hessianDel(120, 121) << " "
        << hessianDel(120, 122) << endl;
   cout << hessianDel(121, 120) << " " << hessianDel(121, 121) << " "
        << hessianDel(121, 122) << endl;
   cout << hessianDel(122, 120) << " " << hessianDel(122, 121) << " "
        << hessianDel(122, 122) << endl;

  */
   
  // --- pertubing hessian matrix -> correction of diagonal entries
  int neighs = 0; 
  SxMatrix<Double> pivot (nAtomsPP*nPeptides, nAtomsPP*nPeptides);
  pivot.set (0.);
   for (int x = 0; x < 3; x++) {
      for (int y = 0; y < 3; y++) {
         for (int iPep = 0; iPep < nPeptides; iPep++) {
            for (int iAtoms = 0; iAtoms < nAtomsPP; iAtoms++) {
               neighs = 0;
               for (int jPep = 0; jPep < nPeptides; jPep++) {
                  for (int jAtoms = 0; jAtoms < nAtomsPP; jAtoms++) {
                     if ( ((inShell (iPep, jPep, iAtoms, jAtoms, rad) && 
                           (modus == 0))
                           ||
                           
                           ((!(inShell (iPep, jPep, iAtoms, jAtoms, rad))) &&
                           (modus == 1)))
                           
                           &&
                           (!((iAtoms == jAtoms) && (iPep == jPep)))
                           )
                     {
                        neighs++;
                        ii = (iPep*nAtomsPP + iAtoms)*3;
                        jj = (jPep*nAtomsPP + jAtoms)*3;
                        hessianDel(ii +x ,ii + y) = 
                           hessianDel(ii + x, ii + y) 
                           - (scaling - 1.) *(0.5*(
                                    hessianDel (ii + x,jj + y) +
                                   hessianDel (ii + x,jj + y)))   
                                 ;
                     }
                  }
               }
            }
         }
      }
   }
  
  // --- pertubing hessian matrix -> rescaling of force constants 
   for (int x = 0; x < 3; x++) {
      for (int y = 0; y < 3; y++) {
         for (int iPep = 0; iPep < nPeptides; iPep++) {
            for (int iAtoms = 0; iAtoms < nAtomsPP; iAtoms++) {
               neighs = 0;
               for (int jPep = 0; jPep < nPeptides; jPep++) {
                  for (int jAtoms = 0; jAtoms < nAtomsPP; jAtoms++) {
                     if ( ((inShell (iPep, jPep, iAtoms, jAtoms, rad) && 
                           (modus == 0))
                           
                           || 
                           ((!(inShell (iPep, jPep, iAtoms, jAtoms, rad))) &&
                           (modus == 1)))
                           
                           &&
                           (!((iAtoms == jAtoms) && (iPep == jPep)))
                           )
                     {


                        neighs++;
                        ii = (iPep*nAtomsPP + iAtoms)*3;
                        jj = (jPep*nAtomsPP + jAtoms)*3;
                        
                        hessianDel(ii + x,jj + y) 
                           = scaling* 
                                 (hessianCartPOrder (ii + x,jj + y)); 
                                 
                        
                     }
                  }
               }
            }
         }
      }
   }
  
  
   blockArrayCart = getBlocks(hessianDel);
   hessianDel = getHessianFromBlocksInputOrder (blockArrayCart);
   switchToInputOrder ();
   
   //-- projecting out external forces
   SxArtifactFilter f;
   f.set (tauCart, SxString("z"), false);
   hessianDel = (f | hessianDel);

   //-- updating the hessian with manipulated one
   setHessian (hessianDel, 100, false);
}

 
   
  

void SxSecondaryStructure::printFC (const SxSymbolTable *symbolTable)  
{
   SxString fileName = 
            symbolTable->get("file")->toString ();

   SxComplex16 sum, coeff, df, f, ia1, ia2, ia, iaFull;

   FILE *fp = NULL;
   if ( !(fp = fopen (fileName.ascii (), "w")) ) {
      sxprintf ("Can't open file %s",fileName.ascii ()); 
      SX_EXIT;
   }
   
   //--- store old cartesian hessian 
   SxMatrix<Double> hessianCartStore (nDoF, nDoF);
   hessianCartStore = getFullCartHessian();
   
   rescaleForceConstants (0, 1.2, nAtomsPP + nAtomsPP/2);
   printPhononDispersionRelation ("disp_mod_nn.out", "splineForLowest", getPDRes ());
  
   setHessian (hessianCartStore, 100, false);

   rescaleForceConstants (1, 1.2, nAtomsPP + nAtomsPP/2);
   printPhononDispersionRelation ("disp_mod_hb.out", "splineForLowest", getPDRes ());
   
   SX_EXIT;  
}

void SxSecondaryStructure::spline 
(const SxArray<double> &x, const SxArray<double> &y, int n, 
 double yp1, double ypn, SxArray<double> *y2) 
{
   int i, k;
   double p, qn, sig, un;

   SxArray<double> u (x.getSize ());

   if (yp1 > 0.99e30)
      (*y2)(0) = u(0) = 0.;
   else {
      (*y2)(0) = -0.5;
      u(0) = (3.0/(x(1) - x(0)))
           * ((y(1) - y(0))/(x(1) - x(0)) - yp1);
   }

   for (i = 1; i <= (n-2); i++) {
      sig = (x(i) - x(i-1))/(x(i+1) - (x(i-1)));
      p   = sig* (*y2)(i-1) + 2.;
      (*y2)(i) = (sig - 1.)/p;
      u(i)  = (y(i+1) - y(i))/(x(i+1) - x(i)) - (y(i) - y(i-1))/(x(i) - x(i-1));
      u(i)  = (6.0*u(i)/(x(i+1) - x(i-1)) -sig*u(i-1))/p;
   }

   if (ypn > 0.99e30)
      qn = un = 0.0;
   else {
      qn = 0.5;
      un = (3.0/(x(n-1) - x(n-2)))*(ypn - (y(n-1) - y(n-2))/(x(n-1) - x(n-2)));
   }
   (*y2)(n-1) = (un - qn*u(n-2))/(qn* (*y2)(n-2) + 1.);

   for (k = (n-2); k>= 0; k--) 
      (*y2)(k) = (*y2)(k)*(*y2)(k + 1) + u(k);

}
      
double SxSecondaryStructure::splint 
(const SxArray<double> &xa, const SxArray<double> &ya, 
 const SxArray<double> &y2a, int n, double x)
{
   int klo, khi, k;
   double h,b,a;

   klo = 0;
   khi = n-1;

   while ( (khi-klo) > 1) {
      k = (khi + klo) >> 1;
      if (xa(k) > x) khi = k;
      else klo = k;
   }

   h = (xa(khi) - xa(klo));

   if (h == 0.) {
      cout << "Bad xa input to routine splint" << endl;
      SX_EXIT;
   }

   a = (xa(khi) - x)/h;
   b = (x - xa(klo))/h;

   double y = a*ya(klo) + b*ya(khi)+((a*a*a - a)*y2a(klo) 
            + (b*b*b - b)*y2a(khi))*(h*h)/6.0;
   return y;
}


void SxSecondaryStructure::printBlock 
(const SxMatrix<Double> &toPrint,const SxString &fileName, bool atomWise)
{
   SxBinIO io;
   int i, j;
   SxMatrix<Double> block(toPrint);
   if (atomWise){
      block = SxMatrix<Double> (nDoFPU/3, nDoFPU/3);
      for (i = 0; i < (nDoFPU/3); i++) {
         for (j = 0; j < (nDoFPU/3); j++) 
            block(i, j) = (toPrint(3*i, 3*j));
      }
   }
   
   
/*
   for (int x = 0; x < nDoFPP; x++) {
      for (int y = 0; y < nDoFPP; y++) {
         if (((x/3) > (y/3)) && (i >= 1))  block(x, y) = 0.;
      }
   }
*/
   io.open (fileName, SxBinIO::ASCII_WRITE_ONLY);
   io.write ("dummy", block, "dummy", "dummy");
   io.close();
}

void SxSecondaryStructure::printBlocks (const SxString &fileName)
{
   SxBinIO io;
   SxString fileNameId;
   SxMatrix<Double> block (nDoFPP, nDoFPP);
   double threshold = 0.;
   double absValue;
   SxHessianOps hOps;
   SxVector<Double> massVecPP (nDoFPU);

   for (int i = 0; i < nDoFPU; i++) massVecPP(i) = masses.coordRef()(i);
   
   for (int i = 0; i < blockArrayCart.getSize (); i++) {
      fileNameId = SxString (fileName + SxString (i) + SxString (".out"));
      io.open (fileNameId, SxBinIO::ASCII_WRITE_ONLY);
      block.copy (blockArrayCart(0)(i));
      block = B.transpose () ^ ( (block^B));
      switchToPeptideOrder ();
      hOps.set (block, massVecPP);
      block.set (hOps.getDynamical ());
     // block.set (hOps.getHessian ());
      for (int x = 0; x < nDoFPP; x++) {
        for (int y = 0; y < nDoFPP; y++) {
           absValue = fabs (block (x, y));
           if (absValue > threshold) block(x, y) = absValue;
           else block(x, y) = 0.;
           if (((x/3) > (y/3)) && (i >= 1))  block(x, y) = 0.;
        }
      }
      
      io.write ("dummy", block, "dummy", "dummy");
      io.close();
   }
   
   blockArrayCart = getPeptideOrderedBlocks (hessianCart);
   hessianCart = getHessianFromBlocks(blockArrayCart);
   hessianCyl = getCylHessian (hessianCart, tauCyl);
   blockArrayCyl = getBlocks (hessianCyl);
   for (int i = 0; i < blockArrayCart.getSize (); i++) {
      for (int j = 0; j < blockArrayCart.getSize (); j++) {
         block.copy (blockArrayCyl(i)(j));
         block = B.transpose () ^ ( (block^B));
         blockArrayCyl (i)(j).copy (block);
      }
   }
  
   
   hessianCyl = getHessianFromBlocks (blockArrayCyl);
  
   for (int x = 0; x < nDoF; x++) {
      for (int y = 0; y < nDoF; y++) {
         absValue = fabs (block (x, y));
         if (absValue > threshold) hessianCyl(x, y) = absValue;
         else block(x, y) = 0.;
         
      }
   }
  
   fileNameId = SxString (fileName +  SxString ("_full.out"));
   io.open (fileNameId, SxBinIO::ASCII_WRITE_ONLY);
   io.write ("dummy", hessianCyl, "dummy", "dummy");
   io.close();
   
   
}
