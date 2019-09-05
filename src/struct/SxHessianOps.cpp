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

#include <SxHessianOps.h>

//----------------------------------------------------------------------------
//    Constructors and Destructors
//----------------------------------------------------------------------------

SxHessianOps::SxHessianOps ()
{
   //empty
} 

SxHessianOps::SxHessianOps 
(  const SxHessianOps &in)
{
   nDoF = in.getNDoF ();
   hessian = SxMatrix<Complex16> (nDoF, nDoF); hessian.set (in.hessian);
   dynamical = SxMatrix<Complex16> (nDoF, nDoF); dynamical.set (in.dynamical);
   massVec = SxVector<Double> (nDoF); massVec.set (in.massVec);
   
   //--- using the constructors as shown in the below two lines
   //   doesn't work properly ... therefore commented out 
   //   however, updating by re-calculating the eigensystems 
   //   is a bit inefficient and should be replaced soon 
   
  // eigDyn   = SxMatrix<Double>::Eigensystem (in.eigDyn);
  // eigHessian 
  //    = SxMatrix<Double>::Eigensystem (in.eigHessian);
   
   needsUpdate = true;
   update ();
}

SxHessianOps::SxHessianOps 
(  const SxMatrix<Complex16> &hessianIn, const SxVector<Double> &massVecIn)
{
   int i, j;
   nDoF = (int)massVecIn.getSize ();
   setMasses (massVecIn);
   setHessian (hessianIn);
   dynamical = SxMatrix<Complex16>(nDoF, nDoF);
   for (i = 0; i < nDoF; i++) {
      for (j = 0; j < nDoF; j++)  
         dynamical(i, j) = hessian(i, j)/sqrt (massVec(i)*massVec(j));
   }
   update ();
}

SxHessianOps::~SxHessianOps ()
{
   // empty
}

void SxHessianOps::set (const SxMatrix<Complex16> &inH, const SxVector<Double> &inM)
{
   int i, j;
   nDoF = (int)inM.getSize ();
   setMasses (inM);
   setHessian (inH);
   dynamical = SxMatrix<Complex16>(nDoF, nDoF);
   for (i = 0; i < nDoF; i++) {
      for (j = 0; j < nDoF; j++)  
         dynamical(i, j) = hessian(i, j)/sqrt (massVec(i)*massVec(j));
   }
   update ();
}

void SxHessianOps::setFromRefinement 
(const SxHessianOps &basisHOps, const SxMatrix<Double> &responsesIn)
{
   int i, j;
   nDoF = basisHOps.getNDoF ();
   SxMatrix<Double> T (nDoF, nDoF);
   SxMatrix<Double> responses (nDoF, nDoF);
   SxVector<Double> col (nDoF);
   SxVector<Double> newCol (nDoF);
  
   setMasses (basisHOps.massVec);
   T.set (basisHOps.eigHessian.vecs);
   
   //--- projection on displacement vectors
   for (i = 0; i < nDoF; i++) {
      col.set (responsesIn.colRef(i));
      for (j = 0; j < nDoF; j++) {
							newCol(j) = 
                        (T.colRef(j)^col).chop ();
      }
      responses.colRef(i).set(newCol);
   }
 
   //--- transformation to cartesian basis
   setHessian (T ^ (responses ^ T.transpose ()));
   dynamical = SxMatrix<Complex16> (nDoF, nDoF);
   for (i = 0; i < nDoF; i++) {
      for (j = 0; j < nDoF; j++)  
         dynamical(i, j) = hessian(i, j)/sqrt (massVec(i)*massVec(j));
   }
   update ();
}


SxHessianOps &SxHessianOps::operator= (const SxHessianOps &in)
{
   if ( &in == this )  return *this;

   nDoF = in.getNDoF ();
   hessian = SxMatrix<Complex16> (nDoF, nDoF); hessian.copy (in.hessian);
   dynamical = SxMatrix<Complex16> (nDoF, nDoF); dynamical.copy (in.dynamical);
   massVec = SxVector<Double> (nDoF); massVec.copy (in.massVec);
   
   //--- using the constructors as shown in the below two lines
   //   doesn't work properly ... therefore commented out 
   //   however, updating by re-calculating the eigensystems 
   //   is a bit inefficient and should be replaced soon 
   
  // eigDyn   = SxMatrix<Double>::Eigensystem (in.eigDyn);
  // eigHessian 
  //    = SxMatrix<Double>::Eigensystem (in.eigHessian);
   
   
   needsUpdate = true;
   update ();
   return *this;
}

//----------------------------------------------------------------------------
//    Member Functions
//----------------------------------------------------------------------------

SxMatrix<Complex16> SxHessianOps::getHessian () {
   return hessian;
}

SxMatrix<Complex16> SxHessianOps::getDynamical () {
   return dynamical;
}

int SxHessianOps::getNDoF () const
{
   return nDoF;
}

SxMatrix<Double> SxHessianOps::getEigenVelocities ()
{
   int i, j;
   SxMatrix<Double> returnValue (nDoF, nDoF);
   update ();
   
   for (i = 0; i < nDoF; i++) {
      returnValue.colRef(i).set (eigDyn.vecs.colRef(i));
        for (j = 0; j < nDoF; j++) 
           returnValue.colRef(i)(j) = 
            returnValue.colRef(i)(j)/sqrt (massVec(j));
   }

   return returnValue;
}

SxVector<Complex16> SxHessianOps::getEigenFrequencies ()
{
   int i;
   double eV;
   SxVector<Complex16> returnValue (nDoF);
   
   update ();
   returnValue.set (0.);
   
   for (i = 0; i < nDoF; i++) {
      eV = eigDyn.vals(i);
      if (eV >= 0.)
         returnValue(i).re = sqrt(eV)*5123.75;
      else 
         returnValue(i).im = sqrt(-eV)*5123.75;
   } 
   return returnValue;
}

SxVector<Double> SxHessianOps::getReducedMasses ()
{
   int i, j;
   double norm;
   SxVector<Double> returnValue (nDoF);
   SxVector<Double> eigenVector (nDoF);

   update ();
   for (i = 0; i < nDoF; i++) {
      norm = 0.; 
      eigenVector.set(eigDyn.vecs.colRef(i));
      for (j = 0; j < nDoF; j++) 
         norm += (eigenVector(j)*eigenVector(j)/massVec(j));
      returnValue(i) = 1./norm;
   }
   return returnValue;
}

void SxHessianOps::setMasses (const SxVector<Double> &massVecIn)
{
   needsUpdate = true;
   massVec = SxVector<Double> (nDoF);
   massVec.copy (massVecIn);
}

void SxHessianOps::setHessian (const SxMatrix<Complex16> &hessianIn)
{
   needsUpdate = true;
   hessian = SxMatrix<Complex16> (nDoF, nDoF);
   hessian.copy (hessianIn);
}


void SxHessianOps::update () 
{
   SxMatrix<Double> a (nDoF, nDoF);
   SxMatrix<Double> b (nDoF, nDoF);
   if (needsUpdate) {
      eigDyn = dynamical.eigensystem ();
      eigDyn.vecs.set (getOrthonormalizedMatrix(eigDyn.vecs));
      eigHessian = hessian.eigensystem ();
      eigHessian.vecs.set (getOrthonormalizedMatrix(eigHessian.vecs));
      
      needsUpdate = false;
   }
}
      
void SxHessianOps::setCurvaturesDynBasis 
                     (SxList<SxList<double> > &toReplace) 
{
   int i, j, index;
   SxMatrix<Double> diagonal (nDoF, nDoF);
   SxVector<Double> reducedMasses = getReducedMasses ();
   cout << reducedMasses << endl;
   SxVector<Double> column(nDoF); SxVector<Double> orth(nDoF);
   diagonal.set (0.);
   for (i = 0; i < nDoF; i++)  diagonal(i, i) = eigDyn.vals (i);
   needsUpdate = true;

   //--- to avoid inconsitencies a symmetrization of the 
   //    new curvatures according to the degenracy of 
   //    the old curvatures is necessary 
   bool changed = true;
   double k1, k2;
   while (changed) { 
      changed = false;
      for (i = 0; i < (toReplace.getSize () - 1); i++) {
         index = (int)(toReplace(i)(0) - 1);
         k1 = toReplace(i)(1);
         k2 = toReplace(i + 1)(1);
         if (   (fabs (diagonal (index, index) 
                    - diagonal(index + 1, index + 1)) < 1e-7) 
             && (fabs (k1 - k2) > 1e-7) ) {
            changed = true;  
            k1 = (k1 + k2)/2.;
            toReplace(i)(1) = k1;
            toReplace(i + 1)(1) = k1;
         }
      }
   }
   

   for (i = 0; i < toReplace.getSize (); i++) {
      index = (int)(toReplace(i)(0) - 1);
      if ((index < 0) || (index >= nDoF)) {
         cout << "Ill-defined index for replacing curvature: "
              << index << endl;
         cout << "Should be in the DoF-range : 1 ..." << nDoF << endl;
         SX_QUIT;
      }
      diagonal (index, index) = toReplace(i)(1);
   } 
     
   SxMatrix<Double> BW (nDoF, nDoF);
   SxMatrix<Double> BWT (nDoF, nDoF);

   BW.set(eigDyn.vecs);

   BWT.set(BW); BWT = BWT.inverse ();
   dynamical = BW ^ (diagonal ^ BWT);
   
   for (i = 0; i < nDoF; i++) {
      for (j = 0; j < nDoF; j++) {
         hessian(i, j) = sqrt(massVec(i)*massVec(j))*dynamical(i, j);
      }
   }

}

void SxHessianOps::setFrequencies 
                     (SxList<SxList<SxComplex16> > &toReplace) 
{
   int i, index;
   SxVector<Double> reducedMasses = getReducedMasses ();
   SxList<SxList<double> > curvatures;
   SxList<double> listEntry;
   double k;
   SxComplex16 freq;

   for (i = 0; i < toReplace.getSize (); i++) {
      //rescale to curvatures consumerable by setCurvaturesDynBasis
      freq = toReplace(i)(1);
      index = (int)(toReplace(i)(0).re - 1);
      k = (freq.abs () * freq.abs ())
          /5123.75/5123.75;
      //k = k/reducedMasses(index)/reducedMasses(index);
      
      listEntry.resize (0);
      listEntry.append (index + 1);
      listEntry.append (k);
      curvatures.append (listEntry);
   } 
   setCurvaturesDynBasis (curvatures);
}

void SxHessianOps::write (const SxString &fileName)
{
   SxBinIO io;
   io.open (fileName, SxBinIO::SX_WRITE_ONLY);
   io.write ("matrix", (SxMatrix<Double>)hessian, "dummy", "dummy");
   io.close();

   SxVector<Complex16> vals = getEigenFrequencies ();
   io.open (fileName, SxBinIO::SX_APPEND_ONLY);
   io.write ("frequencies (unit 1/cm)", vals, "dummy", 0);
   io.close();

   SxString fileNameRaw = SxString (fileName + SxString (".out"));
   io.open (fileNameRaw, SxBinIO::ASCII_WRITE_ONLY);
   io.write ("dummy", (SxMatrix<Double>)hessian, "dummy", "dummy");
   io.close();

}

SxVector<Double> SxHessianOps::getDisplacementVector (int DoF) 
{
   SxVector<Double> returnValue (nDoF);
   returnValue.set (eigHessian.vecs.colRef (DoF));

   
   //--- displacement length is normed according to curvature of PES 
   double curv = fabs(eigHessian.vals(DoF).re);
   
   //--- threshold value (if you change this, please also change the according 
   //    value in getHessianFromRefinement)
   
   if (curv < 1e-9) curv = 1e-9;

   returnValue = returnValue/sqrt(curv);
   return returnValue;
}

SxMatrix<Double> SxHessianOps::getSymmetrizedMatrix 
(const SxMatrix<Double> &h) 
{
   int i, j;
   int size = (int)h.nRows ();
   double avg;
   SxMatrix<Double> returnValue (size, size);
   for (i = 0; i < size; i++) {
      for (j = 0; j <= i; j++) {
         avg = (h(i, j) + h(j, i))/2.;
         returnValue(i, j) = returnValue(j, i) = avg;
      }
   }
   return returnValue;
}

SxMatrix<Complex16> SxHessianOps::getOrthonormalizedMatrix 
(const SxMatrix<Complex16> &h) 
{
   int i, j, k;
   bool found = false;
   int size = (int)h.nRows ();
   SxMatrix<Complex16> returnValue (size, size);
   SxVector<Complex16> col (size);
   SxVector<Complex16> orth (size);

   for (i = 0; i < size; i++) {
      col.set (h.colRef(i));
      //--- for portability purposes (not all platforms return the samle signum)
      k = 0;
      found = false;
      while (!found) {
         if (fabs(col(k).re) > 1e-15) {
            if (col(k).re <= 0.)
               col = -col;
            found = true;
         } else k++;
      }
         
		
      for (j = 0; j < i; j++) {
         orth.set (returnValue.colRef(j));
         col = col - (col^orth).chop()*orth;
         col = col/sqrt(col.absSqr().sum());
      }
      returnValue.colRef(i).set (col);
   }
     return returnValue;
}
         
void SxHessianOps::printMolden (const SxString &fileName, 
                               const SxAtomicStructure &tau,
                               int REPLICZ)
                               
{
   FILE *fp = NULL;
   int i,  zReplic, j;
   SxAtomicStructure tauLoc; tauLoc.copy (tau);
   double zExt = tau.cell(2)(2);
    
   SxString line, numberOne, numberTwo, numberThree, numberZero, toSave,
            innerLine;
   SxVector<Double> help(nDoF);
   SxList<SxString> chemName;
   SxVector<Complex16> eigenFrequencies = getEigenFrequencies ();
   SxMatrix<Double> eigenVelocities = getEigenVelocities ();
   
   //--- parsing in chemical elements
   SxElemDB elemDB;
   
   for (i = 0; i < nDoF/3; i++) 
      chemName.append (elemDB.getChemSymbol (massVec(i*3)));
   
   line = "[Molden Format]\n";
   
   line += "[FREQ]\n";
   for (i = 0; i < nDoF; i++) {
      if (eigenFrequencies(i).re > eigenFrequencies(i).im)
         line += SxString ((float) eigenFrequencies(i).re);
      else {
         line += SxString ((float) eigenFrequencies(i).im);
         line += "*I";
      }
      line += "\n";
   }
   
   line += "[FR-COORD]\n";

   // --- the input-variable zReplic contains how often the structure should
   //     be replicated in the z-Direction, this is a very special option
   //     for being able to vizualize infinite helical molecules 

   for (zReplic = 1; zReplic <= REPLICZ; zReplic++) {  
      for (i = 0; i < nDoF/3; i++) {
            // --- has to be replaced by something like element(is) 
            //     (chemical formula)
            numberZero = chemName(i);
            numberOne = SxString(tauLoc.ref(i)(0));
            numberTwo = SxString(tauLoc.ref(i)(1));
            numberThree = SxString(tauLoc.ref(i)(2) 
                  + (double)(zReplic - 1.)*zExt) ;
            // NgString numberThree (tau(is)(ia)(2)); 
            
            toSave = 
               numberZero + SxString(" ") + numberOne + SxString(" ")
               + numberTwo + SxString(" ")
               + numberThree + SxString("\n");
            line += toSave;
         }
   } 
   
   
   line += "[FR-NORM-COORD]\n";
   
   for (i = 0; i < nDoF; i++) {
      line += SxString("dummy ") + SxString(i) + SxString("\n");
      help.set (eigenVelocities.colRef(i));
     
      innerLine = SxString();
      for (zReplic = 1; zReplic <= REPLICZ; zReplic++) {  
         for (j = 0; j < nDoF; j++) {
            innerLine += SxString(help(j)); 
            innerLine += SxString(" ");
            j++;
            innerLine += SxString(help(j)) + SxString(" ");
            j++;
            innerLine += SxString(help(j)) + SxString("\n");
         }
      }
      line += innerLine;
   }
   
   if ( !(fp = fopen (fileName.ascii (), "w")) ) {
      sxprintf ("Can't open file %s", fileName.ascii ()); 
      SX_EXIT; 
   }     
   fprintf (fp, "%s\n", line.ascii ());
   
   fclose (fp);
}
 
