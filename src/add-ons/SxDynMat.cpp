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

#include <SxDynMat.h>
#include <SxCLI.h>
#include <SxParser.h>
#include <iostream> 
#include <fstream>

#ifndef SX_STANDALONE

SxCell SxDynMat::editInput (const SxAtomicStructure &initForc, 
                          const SxAtomicStructure &initSpin,
                          const SxVector<Double> &recM)
{  cout << "editInput" << endl;
  
   /* set primitive cell with the symmetries allowed by the initial spin 
      structure (including the identity)*/
   Coord initSpinAvg = initSpin.sum () / initSpin.getNAtoms ();
   if ((initSpin - initSpinAvg).absSqr ().sum () > 1e-6 * initSpin.getNAtoms ())  {
      cout << "Spin structures are currently not supported!" << endl;
      SX_EXIT;
   }
   //SxCell primCell = str.getSpinSyms (initSpin, false);
   SxPtr<SxSymGroup> syms;
   syms = syms.create (str.getSymGroup(SxAtomicStructure::SupercellCompatible));
   SxCell primCell = syms->getPrimCell ();
   syms->setToPrimitive (); // remove translational symmetries
   primCell.symGroupPtr = syms;

    
   /* split the species to sets of translation equivalent atoms. The cell,
    * symmetries and the map to the parent are transferred*/
   SxAtomicStructure splitStr = str.splitSpecies (primCell);

   /* set the no. atoms in the primitive cell and the no. species and print to
      standard out.*/
   int nPrim = splitStr.getNSpecies (),
       nSpec = str.getNSpecies ();
   cout << nPrim << "nPrim, nSpec" << nSpec << endl;

   //--- enlarge recMass to nPrimAt (used in hermit)
   recMass.resize (nPrim);
   int iSpec = 0; //counters over direction and species
   for (int iPrim = 0; iPrim < nPrim; ++iPrim)  {

      // go to next species if it has changed
      if (str.getISpecies (splitStr.getIAtom (iPrim, 0)) != iSpec)  
      { iSpec++; }
      
      // transfer the correct recM
      recMass (iPrim) = recM(iSpec);
   } 

   int nDof = (int)dTauStr.getSize (), //no. degrees of freedom
       nAtom = str.getNAtoms (); //no. atoms
   Coord drift; //the drift
   for (int iDof = 0; iDof < nDof; iDof++)  {

      // correct for initial forces
      forcStr(iDof).replaceInfo (str.atomInfo);
      forcStr(iDof) -= initForc;

      // reorder the atoms, but keep the atomInfo to have the correct nSpec
      forcStr(iDof) = SxAtomicStructure (splitStr.atomInfo, forcStr(iDof));
      forcStr(iDof).replaceInfo (str.atomInfo);
      
      // correct drift BEFORE the mass rescaling TODO apply a force filter instead
      drift = forcStr(iDof).sum ();
      forcStr(iDof) -= (drift / nAtom); //the drift is equally ditributed
    
      // scale the forces with the masses (this doesn't influence the symmetries)
      forcStr(iDof) *= sqrt (recM);

      // determine the displacement structures
      dTauStr(iDof).replaceInfo (str.atomInfo);
      dTauStr(iDof) = (dTauStr(iDof) - str).map (str.cell, SxCell::Origin);

      /* reorder the atoms AFTER the mass scaling, but keep the atomInfo to have 
         the correct nSpec*/
      dTauStr(iDof) = SxAtomicStructure (splitStr.atomInfo, dTauStr(iDof));
      dTauStr(iDof).replaceInfo (str.atomInfo);
            
      // check the sizes of the displacement, which could be small
      if (sqrt (dTauStr(iDof).absSqr ().sum ()) < str.epsEqual / 3.)  {
         cout << "Displacement error in SxDynMat::editInput !" << endl
              << "Total abs. dTau for (displ.) structure no. " << iDof << endl
              << "is only " << sqrt (dTauStr(iDof).absSqr ().sum ()) << endl;
         SX_EXIT;
      } 

      // scale dTauStr with recMass AFTER the size check (this doesn't influence
      // the symmetries)
      dTauStr(iDof) /= sqrt (recM); //should epsEqual also be rescaled?
   }

   // set splitStr with the correct symmetries, but keep the atomInfo!
   splitStr.replaceInfo (str.atomInfo);
   str = splitStr;
   str.cell.symGroupPtr = primCell.symGroupPtr; //reference!
   return primCell;
}

SxCell SxDynMat::readInput (const SxString &input, 
                                  ofstream &output,
                            const double &epsE)  
{  
   // check that input not empty
   if (input == "")  {
      cout << "Error in SxDynMat::readInput! empty input provided" << endl;
      SX_EXIT;
   } 
   cout << "read input " << input << endl;

   // output of read process
   SxVector<Double> recM;
   SxAtomicStructure initForc, initSpin;
   try  {
      
      // input to parser
      SxParser parser;
      SxParser::Table parTab = parser.read (input);
      SxString speciesId = "structure";
      if (parTab->containsGroup("pseudoPot"))  speciesId = "pseudoPot"; 
      
      // ---read coordinates, forces and possibly spins. The first set should be
      // the initial structure!
      int nAtom = -1, //no. atoms
          nSpec = -1, //no. species
          nCurSpec, //current no. species
          nCurAtom; //current no. atoms
      SxCell inCell; //initial cell
      SxSymbolTable *strucPtr = NULL;
      SxList<SxAtomicStructure> strL, inForcL, inSpinL; //:nDof + 1
      for (strucPtr = parTab->getGroup ("structure"); strucPtr != NULL;
           strucPtr = strucPtr->nextSibling ("structure"))  
      {
         strL << SxAtomicStructure (strucPtr);
         inForcL << SxAtomicStructure (strucPtr, "force"); //set to 0 if not given
         inSpinL << SxAtomicStructure (strucPtr, "spin"); //set to 0 if not given

         // --some consistency checks
         if (strL.getSize () == 1)  {
            nSpec = strL (0).getNSpecies ();
            nAtom = strL (0).getNAtoms ();
         }
         nCurSpec = strL.last ().getNSpecies ();
         nCurAtom = strL.last ().getNAtoms ();
         if (nCurSpec != nSpec || nCurAtom != nAtom)  {
            cout << "Structure mismatch in SxDynMat::read!" << endl
                 << "Reference structure has " << nSpec << "species and "
                 << nAtom << "atoms." << endl
                 << "Structure no. " << strL.getSize () << " has " 
                 << nCurSpec << "species and " << nCurAtom << "atoms." << endl;
             SX_EXIT;
         }
         
         // ---read the cell of the first structure and check the others
         if (strL.getSize () == 1) inCell = SxCell (strucPtr);
         else if (strucPtr->contains ("cell"))  {
            if ((inCell - strL.last ().cell).absSqr ().sum () > 1.e-4) {
               cout << "Cell mismatch in SxDynMat::read !" << endl
                    << "Reference structure has " << inCell << endl
                    << "Structure no. " << strL.getSize ()
                    << " has " << strL.last ().cell << endl;
               SX_EXIT;
            }
         }
      }

      // ---read recMass
      SxSymbolTable *specPtr = NULL;
      SxList<double> recML; //list of recMass
      for (specPtr = parTab->getGroup (speciesId)->getGroup ("species");
           specPtr != NULL; specPtr = specPtr->nextSibling ("species"))  
      {
         recML << (specPtr->get ("reciprocalMass")->toReal ()); 
      }
      recM = recML;

      if ( recM.getSize () != nSpec )  {
         cout << "Mass mismatch in SxDynMat::read !" << endl
              << "No. masses provided is " << recM.getSize () << endl
              << "No. species is " << nSpec << endl;
         SX_EXIT;
      }

      // ---set the initial coordinates, forces and spins
      str = strL.first (); //initialise str
      str.epsEqual = epsE;
      initForc = inForcL.first ();
      initForc.replaceInfo (str.atomInfo); //share
      initSpin = inSpinL(0); //ignore other spins
      initSpin.epsEqual = epsE; //share
      initSpin.replaceInfo (str.atomInfo); //share

      backgroundForces.copy(initForc);
       
      // ---write some data to stdout as check
      output << "Init. forces.absSqr ().sum () is " 
             << initForc.absSqr ().sum () << endl;

      // --make array of displaced coordinates and forces
      strL.remove (0);
      dTauStr = strL; //:nDof does epsEqual need to be set?
      inForcL.remove (0);
      forcStr = inForcL; //:nDof
   }  catch (SxException e)  {
      e.print();
      SX_QUIT;
   } 

   // edit the input
   SxCell primCell = editInput (initForc, initSpin, recM);
   SX_CHECK (primCell.symGroupPtr);
      
   // write output to file
   int nDof = (int)dTauStr.getSize ();
   output << "undisplaced structure is " << str << endl
          << "supercell is " << str.cell << endl
          << "primitive cell is " << primCell << endl
          << "nSym: " << primCell.symGroupPtr->getSize () << endl
          << "(mass scaled) displacements are " ;
   for (int iDof = 0; iDof < nDof; iDof++)  { 
      output << dTauStr(iDof).coords << endl; 
   }  // this is much data! use print instead to write it to standard out?

   output << endl << "reciprocal masses are " << recM << endl
          << "nDof (provided)" << nDof << endl;
   return primCell;
}
      
SxArray<Coord> SxDynMat::getExactQ (const SxCell &primCell, 
                                    ofstream &output) const
{  cout << "getExactQ" << endl;

   // get the reciprocal cells
   SxCell recPrimCell = primCell.getReciprocalCell (), //including symmetries?
          recCell = str.cell.getReciprocalCell ();
   SxVector3<Int> mult; //lattice multiplicities in 3 directions 

   // set the no. q-points
   int iQ = 0, jQ, 
       nQ = int(str.cell.volume / primCell.volume + .5);

   // print the relative basis to standard out for checks
   cout << "relative basis " //factors 10 against rounding problems
        << SxVector3<Int> (recPrimCell.carToRel (10*nQ * recCell.basis (0)))/10. 
        << SxVector3<Int> (recPrimCell.carToRel (10*nQ * recCell.basis (1)))/10.
        << SxVector3<Int> (recPrimCell.carToRel (10*nQ * recCell.basis (2)))/10.
        << "/" << nQ << endl;

   mult = SxVector3<Int> (nQ, 1, 1);//1st guess of multiplicities
   
   // determine the exact q-points as the multiples of recCell vectors in recPrimCell
   SxVector3<Int> ind; //indices
   SxArray<Coord> res (nQ); //result
   while (iQ < nQ)  {
      ind = SxVector3<Int> (iQ % mult(0), (iQ / mult(0)) % mult(1),
                            iQ / mult(0) / mult(1));
      for (jQ = 0; jQ < iQ; jQ++)  { 
         if (recPrimCell.getMapped ((recCell ^ ind) - res(jQ), SxCell::Origin )
               .norm () < 1e-6)  { 
            break; 
         }
      }
      if (jQ == iQ)  {//not found
         res (iQ) = recCell ^ ind;
         iQ++;
      }  else  {//decrease mult in one direction, increase the next
         if (mult(0) == nQ)  {
            mult(0) = iQ;
            mult(1) = nQ / iQ;
         }  else if (mult(1) == nQ / mult(0))  {
            mult(1) = iQ / mult(0);
            mult(2) = nQ / iQ;
         }  else  {
            cout << "Error in SxDynMat::getExactQ!" << endl
                 << "could not find " << nQ << "exact q-points" << endl
                 << "found exactQ: " << res << endl;
            SX_QUIT;
         }
      }
   }

   // --check that we have a complete set
   ind = SxVector3<Int> ( nQ % mult(0), (nQ / mult(0)) % mult(1),
                          nQ / mult(0) / mult(1) );
   Coord diff;
   for (jQ = 0; jQ < nQ; jQ++)  { 
      diff = recPrimCell.getMapped ((recCell ^ ind) - res(jQ), SxCell::Origin);
      if (diff.norm () < 1e-6)  { break; } 
   }
   if ( jQ == nQ )  {//not found
      cout << "Error in SxDynMat::getExactQ!" << endl
           << "no complete set found: recCell ^ " << ind << " misses" << endl;
      SX_EXIT;
   }

   // write result to output file
   output << "exact q-points are " << res  << "in recCel: " << recPrimCell
          << endl; 

   // return q-points
   return res;
}
 
bool SxDynMat::symMod (const int iSym, 
                       const int iDof)  
{ 
   // current symmetry operation
   SxSymOp symOp (str.cell.symGroupPtr->operator() (iSym));

   // get symmetry reordering info
   SxAtomicStructure mapStr = symOp ^ str; 
   SxGrid grid (mapStr, 10); 
   SxConstPtr<SxAtomInfo> reorder = mapStr.match (grid, str);
   SX_CHECK (reorder);
   // this goes wrong when epsEqual is not copied (implicitely), then symOp is
   // no more a symmetry

   // get mapped dTauStr (with the atomInfo of str)
   mapStr = symOp.rot ^ dTauStr(iDof);
   SxAtomicStructure dTauMap (reorder, mapStr);

   // share the atom info for the subtraction!
   dTauMap.replaceInfo (dTauStr(iDof).atomInfo);

   // ---get sign of mapping to remove the inversion duality of the displacements
   double dTauDiff = (dTauStr(iDof) - dTauMap).absSqr ().sum (),
          dTauInvDiff = (dTauStr(iDof) + dTauMap).absSqr ().sum ();
   int sign = 1;
   if (dTauInvDiff < dTauDiff)  {
      dTauDiff = dTauInvDiff;
      sign = -1;
   }

   // ---symmetrize with symOp when mapped to equiv. dTauStr
   SxAtomicStructure forcMap;
   if (dTauDiff < 1.e-5)  {//this is a sqr!
      
      SxAtomicStructure temp (forcStr(iDof), SxAtomicStructure::Copy);

      // --add transformed (and reordered) forces then half them (they have the
      // atomInfo of str)
      mapStr = sign * symOp.rot ^ forcStr(iDof);
      mapStr = SxAtomicStructure (reorder, mapStr);
      mapStr.replaceInfo (forcStr(iDof).atomInfo);//otherwise they will be reordered again!
      forcStr(iDof) += mapStr;
      forcStr(iDof) /= 2.; 

      // symmetrize dTauStr as well 
      temp = SxAtomicStructure (dTauStr(iDof), SxAtomicStructure::Copy);
      dTauStr(iDof) += (sign * dTauMap);
      dTauStr(iDof) /= 2.; 
      return true; //symmetrized
   }  else  { return false; } //not symmetrized
}

int SxDynMat::genNewMods (const SxArray<SxArray<SxSymOp> > &genSyms, //:nDof,nSym
                          const bool debug)
{  
   cout << "genNewMod";
   int iSym, nSym, //counts the symmetries
       nDof = (int)genSyms.getSize (), //no. degrees of freedom
       nMode = nDof; //no. modes, this will increase
   SxConstPtr<SxAtomInfo> reorder;
   SxAtomicStructure mapStr;
   for (int iDof = 0; iDof < nDof; iDof++)  {
      nSym = (int)genSyms(iDof).getSize ();
      dTauStr.resize (nMode + nSym, true);
      forcStr.resize (nMode + nSym, true);
      for (iSym = 0; iSym < nSym; iSym++)  {

         // get symmetry reordering info
         mapStr = genSyms(iDof)(iSym) ^ str;
         SxGrid grid (mapStr, 10); 
         reorder = mapStr.match (grid, str);
         SX_CHECK (reorder);

         // create empty new dTauStr
         dTauStr(nMode + iSym) = dTauStr(iDof).getNewStr ();
         
         // add new dTauStr
         mapStr = genSyms(iDof)(iSym).rot ^ dTauStr(iDof);
         dTauStr(nMode + iSym) = SxAtomicStructure (reorder, mapStr);

         //output to standard out in debug mode
         if (debug)  {
            cout << "newDTau"; 
            dTauStr(nMode + iSym).coords.transpose ().print (true);
         }
  
         // create empty new forcStr
         forcStr(nMode + iSym) = forcStr(iDof).getNewStr ();
         
         // add new forcStr
         mapStr = genSyms(iDof)(iSym).rot ^ forcStr(iDof);
         forcStr(nMode + iSym) = SxAtomicStructure (reorder, mapStr);
      }
      nMode += nSym; //the identity must be excluded
   }
   cout << nMode << endl;
   return nMode;
}

int SxDynMat::applySymOps (ofstream &output, const bool debug)  
{    
   int nDof = (int)dTauStr.getSize (),
       iSym, nSym = str.cell.symGroupPtr->getSize ();
   //these symmetries are only the ones within the primitive cell and include
   //the identity

   SxList<SxSymOp> kernL; // symmetries of dTau(iDof) (kernel)
   SxArray<SxArray<SxSymOp> > genSyms (nDof); //for all dTau

   for (int iDof = 0; iDof < nDof; iDof++)  {
      kernL.removeAll ();
      cout << endl;
      for (iSym = 0; iSym < nSym; iSym++)  {
      
         // symmetrize force mode iDof if dTau invariant 
         if (symMod (iSym, iDof))  {
    
            // add symOp to the kernel of dTau(iDof)
            kernL << str.cell.symGroupPtr->operator()(iSym);
         }
      }
      // print some info about the symmetries 
      output << "iDof" << iDof << "nSym" << kernL.getSize () << endl;
      cout << "kernel" << kernL;
      genSyms(iDof) = *str.cell.symGroupPtr / SxArray<SxSymOp> (kernL);
   }//RIGHT division (symmetry operations do not commute)!
   //genSyms excludes the idenity

   return genNewMods (genSyms, debug);//:nMode
}

void SxDynMat::fourier (const SxArray<Coord> &exactQ)  
{
   cout << "fourier " << endl;

   // --sum contributions from different primitive cells with the apropriate phase factors
   int iQ, nQ = (int)exactQ.getSize (),
       nAtom = str.getNAtoms (),
       iPrimAt, nPrimAt = nAtom / nQ, //no. atoms in the primitive cell (differs from nSpec)
       nMode = (int)dTauStr.getSize (),
       iDir, begin, end;
   SxVector<Complex16> phaseQ (nAtom); //phase factor for certain Q
   SxMatrix<Complex16> dTauPhasM (3, nAtom), //dTau of certain mode with phase factors
                       forcPhasM (3, nAtom), //idem for forces
                       oldPhas (3, nAtom);  //old phase factor
   for (iQ = 0; iQ < nQ; iQ++)  {

      // standard phase factor (complex)
      phaseQ = expI (exactQ(iQ) ^ str);

      for (int iMode = 0; iMode < nMode; iMode++)  {

         // multiply with phase factors & put dTau, forc in matrices
         dTauPhasM = dTauStr(iMode).coords * phaseQ; //:3*nAtom
         forcPhasM = forcStr(iMode).coords * phaseQ; //Matrix * Vector

         // sum over cells
         for (iPrimAt = 0; iPrimAt < nPrimAt; iPrimAt++)  {
            begin = iPrimAt * nQ;
            end = begin + nQ - 1;
            for (iDir = 0; iDir < 3; iDir++)  {
               dTau(iQ)(iDir + 3 * iPrimAt, iMode) 
                  = dTauPhasM.row (iDir).sum (begin, end);
               forces(iQ)(iDir + 3 * iPrimAt, iMode) 
                  = forcPhasM.row (iDir).sum (begin, end);
            }
         }
      }
   }
}   

SxArray<SxMatrix<Complex16>::Eigensystem> SxDynMat::getDynMat (
       const SxArray<Coord> &exactQ, const bool debug) 
{ 
   // check that exactQ not empty
   int nQ = (int)exactQ.getSize ();
   if (nQ == 0)  {
      cout << "ERROR in SxDynMat::reduceToPrimCell! "
           << "no qPoints provided, qPoints = " << exactQ << endl;
      SX_QUIT;
   }
   cout << "getDynMat" << endl;
 
   // --- resize matrices
   dTau.resize (nQ);
   forces.resize (nQ);
   int nPrimAt = str.getNAtoms () / nQ,
       nMode = (int)dTauStr.getSize ();
   for (int iQ = 0; iQ < nQ; iQ++)  {
      dTau(iQ) = SxMatrix<Complex16> (3 * nPrimAt, nMode);
      forces(iQ) = SxMatrix<Complex16> (3 * nPrimAt, nMode);//row is fast index! 
   }

   //make fourier transform of forces and displacements
   fourier (exactQ);
   cout << "transformed" << endl;
   
   SxMatrix<Complex16> dynMatQ;
   SxArray<SxMatrix<Complex16>::Eigensystem> res (nQ);
   SxMatrix<Complex16> inverse, identity;
   bool gamma;
   cout << "calculate inverse of displacements."
        << " This goes wrong if they don't span the phase space!\n"
        << "If so then check the symmetries and the data provided\n";
   for (int iQ = 0; iQ < nQ; iQ++)  {
      inverse = dTau(iQ).inverse ();//this is a right inverse

      // --check inverse
      identity = dTau (iQ) ^ inverse;
      if (identity.trace ().re < double(dTau(iQ).nRows ()) * .9999)  {
         cout << "Warning from SxDynMat::getDynMat: inverse might be wrong!" 
              << endl << identity << " should be the identity! trace:" 
              << identity.trace () << endl << "did you provide sufficient data?";
         if (debug)  {
            cout << "dTau,inverse: ";
            dTau (iQ).print (true);
            inverse.print (true);
         }
         SX_EXIT;
      }
      dynMatQ = - (forces(iQ) ^ inverse);

    // --check and correct Hermiticity of dynMat
      gamma =  (exactQ(iQ).norm () < 1.e-6);
      if (! hermit (dynMatQ, gamma, debug))  { cout << "for iQ" << iQ << endl; }

      dynMatQ *= sqr (AU2MEV); //rescale to meV^2
      res(iQ) = dynMatQ.eigensystem (); //reference
   }
   // --print extra information in debug mode
   if (debug)  {
      cout << "eigVec orthogonal? (dot products)";
      int jMode;
      for (int iMode = 0; iMode < 3*nPrimAt - 1; iMode++)  {
         for (jMode = iMode + 1; jMode < 3*nPrimAt; jMode++)  {
            cout << " " << (res(0).vecs.col (iMode) ^ res(0).vecs.col (jMode));
         }
         cout << endl;
      }
      cout << endl;
   }
   cout << "dynMats constructed & Hermiticity checked" << endl;
   return res;
}

bool SxDynMat::hermit (SxMatrix<Complex16> &dynMatQ, 
                       bool driftCor, 
                       const bool debug) 
{ 
   // --check hermiticy of dynMatQ
   SX_CHECK (dynMatQ.nCols () == dynMatQ.nRows ());
   SX_CHECK (dynMatQ.nCols () % 3 == 0);
   bool herm = true;
   SxMatrix<Complex16> asymMat = dynMatQ - dynMatQ.transpose ().conj ();
   double asym = asymMat.norm ();
   if (asym * 10. > dynMatQ.norm ())  
   {  cout << "Warning! dynMat " << (asym / dynMatQ.norm () / 2.) 
           << " fraction anti-hermitian. norm:" << dynMatQ.norm ();
      herm = false; 
      if (debug)  {
         cout << "asymMat:"; 
         asymMat.print (true);
      }
   }

   // make hermitian
   dynMatQ += dynMatQ.transpose ().conj ();
   dynMatQ *= .5; 
   
   //--- set & correct drift (self-consistently); for the dynMat at Gamma only!
   if (driftCor)  { 

      // check reality dynMat
      SX_CHECK (100. * dynMatQ.imag ().norm () < dynMatQ.norm ());
      
      bool repeat = true;
      int index, iMode, nMode = (int)dynMatQ.nCols ();
      SxVector<Double> drift (3 * nMode, 0.);
      SxVector<Complex16>::Iterator dmIt;
      while (repeat)  { 

         //--- get drift (columnwise 3-fold sum)
         dmIt = dynMatQ.begin ();
         for (iMode = 0; iMode < nMode; iMode++)  {
            for (index = 0; index < nMode; index++, dmIt++)  {
               drift (3 * iMode + (index % 3)) += (*dmIt).re 
                                              / sqrt (recMass(index / 3));
            } 
         }

         // check if further repeat needed 
         if (drift.norm () * sqr (nMode) < 1.e-3 * dynMatQ.norm ())  
         { repeat = false; }

         // rescale drift per atom
         drift *= 3. / nMode;
         if (debug)  { 
            cout << dynMatQ.norm () << ":norm dynMat, drift:" 
                 << drift << "corrected" << endl; 
         }

         //--- correct drift 
         dmIt = dynMatQ.begin ();
         for (iMode = 0; iMode < nMode; iMode++)  {
            for (index = 0; index < nMode; index++, dmIt++)  {
               *dmIt -= drift (3 * iMode + (index % 3)) 
                      * sqrt (recMass(index / 3));
            } 
         } 
         // make again Hermitian 
         dynMatQ += dynMatQ.transpose ().conj ();
         dynMatQ *= .5; 
      }
   }
   return herm;
}
  
void SxDynMat::print (const SxArray<SxMatrix<Complex16>::Eigensystem> &eig,
                      const SxArray<Coord> &qPoints,
                      const SxCell &primCell) const
{  cout << "print" << endl;

   // ---check the no. q-points
   int nQ = (int)qPoints.getSize ();
   if (nQ == 0 || nQ != eig.getSize ())  {
      cout << "Size error in SxDynMat::print! "
           << "provided no. q-points is " << nQ
           << "no. dynamical matrices is " << eig.getSize () << endl;
      SX_EXIT;
   }
   
   // set the primitive structure
   SxAtomicStructure primStr = str.getPrimStr (primCell);
   primStr.map (SxCell::Origin);
  
   // ---print primStr, qPoints & eig to sxb file
   try  {
      SxBinIO io ("dynmat.sxb", SxBinIO::BINARY_WRITE_ONLY);

      // ---first write all headers
      io.addDimension ("nDir", 3);
      io.write ("primCell", primCell, "nDir"); //SxMatrix3
       
      // primPos: SxVector (3*nAtom)
      io.addDimension ("nAtom", primStr.getNAtoms ());
      io.write ("primPos", primStr.coords, "nDir", "nAtom");
      
      // qPoints: SxMatrix (3, nQ)
      io.addDimension ("nQ", nQ);
      io.write ("exactQ", SxMatrix<Double> (3, nQ), "nDir", "nQ");
      
      // eigensystem: values SxMatrix (nMode, nQ)
      int nMode = (int)eig(0).vals.getSize ();
      io.addDimension ("nMode", nMode);
      io.write ("eigVals", SxMatrix<Double> (nMode, nQ), "nMode", "nQ");

      // eigensystem: vectors SxVector (sqr(nMode)*nQ)
      io.addDimension ("vecSize", sqr (nMode) * nQ);
      io.write ("eigVecs", SxVector<Complex16> (sqr (nMode) * nQ), "vecSize");
      
      // ---then write all data
      io.setMode (1);
      io.write ("primCell", primCell, "nDir");
      io.write ("primPos", primStr.coords, "nDir", "nAtom");
      
      // -make a matrix of the q-points
      SxMatrix<Double> qP (3, nQ);
      for (int iQ = 0; iQ < nQ; iQ++)  { 
         qP.colRef (iQ) <<= qPoints(iQ);
      }
      io.write ("exactQ", qP, "nDir", "nQ");
      
      // -make a matrix of the eigenvalues
      SxMatrix<Double> vals (nMode, nQ);
      for (int iQ = 0; iQ < nQ; iQ++)  {
         vals.colRef (iQ) <<= eig(iQ).vals;//vals are vectors
      }
      io.write ("eigVals", vals, "nMode", "nQ");
      
      // make one vector of the eigenvectors
      SxMatrix<Complex16> vecs (sqr(nMode), nQ);
      for (int iQ = 0; iQ < nQ; iQ++)  {

         // print the imaginary parts to standard out
         cout << "imag%vec(iQ):" << int (100. * eig(iQ).vecs.imag ().norm () / nMode); 
         
         vecs.colRef (iQ) <<= eig(iQ).vecs; //vecs are columnwise vectors
      }
      cout << endl;
      vecs.reshape (sqr(nMode) * nQ);
      io.write ("eigVecs", vecs, "vecSize");
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}

void SxDynMat::printHesse (const SxArray<SxMatrix<Complex16>::Eigensystem> &eig,
                      const SxArray<Coord> &exactQ, const SxCell &primCell, const bool hesse,
                      const bool sxhesse, const bool bgForces, const double disp) const
{
  int nQ = int(exactQ.getSize ());
  int iQ, iPrim, jPrim, nDof;
  FILE *file;
  int i, j;
  Coord deltaR;

  SxArray<SxMatrix<Complex16> > dynMatQ (nQ);//(nAtDir, nMode);
  for (iQ = 0; iQ < nQ; iQ++)  {
     //- make dynMat from eigensystem: D=Vecs*Vals*Vecs^-1 (columnwise basis)
     dynMatQ(iQ) = eig(iQ).vecs * eig(iQ).vals; //matrix * vector 
     dynMatQ(iQ) = dynMatQ(iQ) ^ eig(iQ).vecs.inverse ();
  }
  
  nDof = 3*str.getNAtoms ();
  SxMatrix<Double> realH(nDof,nDof);
  SxComplex16 phase, Dij;
  SxAtomicStructure primStr = str.getPrimStr (primCell);
  primStr.map (SxCell::Origin);
  SxGrid grid (primStr, 3); // auxiliary for atom finding
  SX_LOOP(ia)  {
     iPrim = primStr.find (str.getAtom(ia),grid);
     SX_LOOP(ja)  {
        jPrim = primStr.find (str.getAtom(ja),grid);
        deltaR = str.getAtom(ia) - str.getAtom(ja);
        SX_LOOP2(iDir(3),jDir(3))  {
           Dij = 0.;
           for (iQ = 0; iQ < nQ; iQ++) {
              phase = expI (-exactQ(iQ) ^ deltaR);
              Dij += phase * dynMatQ(iQ)(3*iPrim+iDir,3*jPrim +jDir) * (1./nQ);
           }
           // dynMatQ is at this stage in units of meV (i.e. not in the sphinx default units) due to the scaling at the end of getDynMat
           // here we rescale back to a.u. for consistency with Thermodynamics scripts and also with masses to get from dynMat to Hesse
           realH(3*ia+iDir, 3*ja+jDir) = Dij*(1./sqrt(recMass(iPrim)*recMass(jPrim)))/sqr (AU2MEV);
      }
    }
  }
  
  if (hesse) {
    // export hesse
    cout << endl 
         << "Exporting HesseMatrix_sphinx in units of hartree/bohrradius^2" 
         << endl;
    file = fopen("HesseMatrix_sphinx", "w");
    for (i=0; i<nDof; i++) {
      for (j=0; j<nDof; j++) {
        fprintf(file,"%25.15f ",realH(i,j));
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }

  if (sxhesse) {
    // export sxhesse.sx
    cout << endl 
         << "Exporting sxhesse.sx in units of hartree and hartree/bohr and "
         << "using disp of " << disp << " bohr" << endl;
    SxAtomicStructure force (str, SxAtomicStructure::Copy);
    SxVector3<TPrecTauR> dispVec;
    SxVector<TPrecTauR> U(nDof), F(nDof);
    file = fopen("sxhesse.sx","w");
    fprintf(file,"format phononDat;\n\n");
    fprintf(file,"pseudoPot  {\n");
    for (i=0; i<str.getNSpecies(); i++) {
      iPrim = primStr.find (str.getAtom(i,0),grid);
      fprintf(file,"  species { reciprocalMass = %.10f ; }\n", recMass(iPrim));
    }
    fprintf(file,"}\n\n");
    // equilibrium structure first
    force *= 0.;
    // add background forces if requested
    if (bgForces) force += backgroundForces;
    fprintf(file,"# equilibrium structure\n");
    str.fprint(file,force);
    for (i=0; i<str.getNAtoms(); i++) {
     for (int iDir=0; iDir<3; iDir++) {
      U = 0.;
      U(i*3+iDir) = disp;
      F = -realH^U;
      for (j=0; j<str.getNAtoms(); j++) {
        force.ref(j) = SxVector3<TPrecTauR> (F(3*j+0), F(3*j+1), F(3*j+2));
      }
      // add background forces if requested
      if (bgForces) force += backgroundForces;
      dispVec = SxVector3<TPrecTauR> (0,0,0);
      dispVec(iDir) += disp;
      const_cast<SxAtomicStructure&>(str).ref(i) += dispVec;
      // displaced structures
      fprintf(file,"# structure with displacement of atom %d in %c direction\n",i+1,'x'+iDir);
      str.fprint(file,force);
      const_cast<SxAtomicStructure&>(str).ref(i) -= dispVec;
     }
    }
    fclose(file);
  }
}

#else /* SX_STANDALONE */

int main (int argc, char **argv)
{
   // ---command line interface
   SxCLI cli (argc, argv);
   double epsE = cli.option ("-e|--epsEqual", "accuracy structure")
                             .toDouble (5.e-3);
   SxString input = cli.option ("-i|--input", "file", 
                                "input").toString ("phononDat.sx");
   bool debug = cli.option ("-v|--verbose", "print some extra info").toBool ();
   bool sxhesse = cli.option ("-HS|--HesseSX", "print sxhesse.sx file").toBool ();
   bool hesse = cli.option ("-H|--Hesse", "print Hesse matrix to file").toBool ();
   bool bgForces = cli.option ("-k|--keepBG", "keep background forces in sxhesse.sx").toBool ();
   double disp = cli.option ("-d|--disp", "use disp for displacement in sxhesse.sx")
                             .toDouble (0.01);
   cli.finalize ();
   
   // -open output file
   ofstream output;
   output.open ("dynmat.out");
   cout << "dynmat.out opened" << endl;

   // -read input
   output << "Input file is " << input << " epsEqual is " << epsE << endl;
   SxDynMat D;
   SxCell primCell = D.readInput (input, output, epsE);

   // - get the exact q-points of the supercell
   SxArray<Coord> exactQ = D.getExactQ (primCell, output);

   // symmetrize & generate mapped displacements, forces 
   int nMode = D.applySymOps (output, debug), 
       nAtom = D.str.getNAtoms (), 
       nCell = (int)exactQ.getSize ();

   // test minimum no. modes (equivalent displacements could be included)
   if (nMode < 3 * nAtom / nCell)  {
      cout << "Error in SxDynMat" << endl 
           << "Too few modes generated!: " << nMode << endl
           << "There should be at least: " << 3 * nAtom / nCell << endl;
      SX_EXIT;
   }

   // get the dynamical matrices (eigensystems) for the exact q-points
   SxArray<SxMatrix<Complex16>::Eigensystem> eig = D.getDynMat (exactQ, debug);

   // --write the frequencies to ASCII output
   int iMod, nImagFreq = 0;
   SxComplex16 eigVal;
   nMode = (int)eig(0).vals.getSize ();
   output << "exactFreq"; 
   for (int iQ = 0; iQ < eig.getSize (); iQ++)  {
      output << endl << iQ << " ";
      for (iMod = 0; iMod < nMode; iMod++)  {
         eigVal = eig(iQ).vals (iMod);
         if (eigVal.re < -1.e-4)  { nImagFreq++; }
         output << (eigVal.re / sqrt (eigVal.abs ())) << " ";
      }
   }
   output << "\n #There are " << nImagFreq << " imaginary frequencies\n";

   // close output
   output.close ();

   // write the dynamical matrices to binary output
   D.print (eig, exactQ, primCell);

   // write hesse matrix to file if requested
   if (hesse||sxhesse) D.printHesse (eig, exactQ, primCell, hesse, sxhesse, bgForces, disp);

   // end program
   return 0;
}
#endif /* SX_STANDALONE */
