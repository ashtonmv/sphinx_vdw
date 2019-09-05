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
//   AUTHOR: Lars Ismer, ismer@fhi-berlin.mpg.de

#include <SxArtifactFilter.h>
SxArtifactFilter::SxArtifactFilter ()
{
     // empty
}


void SxArtifactFilter::set (const SxAtomicStructure &tauIn, 
                                    const SxString &freezeModeIn,
                                    bool massWeightingIn) 
{
   needsUpdateDM = true;
   setStructure (tauIn);
   setFreezeMode (freezeModeIn);
   setMassWeighting (massWeightingIn);
   DMatrix = getDMatrix ();
   // empty
}

SxArtifactFilter::~SxArtifactFilter ()
{
   // empty
}

SxVector<TPrecTauR> 
SxArtifactFilter::operator*(const SxVector<TPrecTauR> &in) const 
{
  if (freezeMode != SxString ("notr")) {
   if ((in.nRows () > 1) && (in.nCols () > 1)) {
      SxMatrix<Double> returnValue(in.nRows (), in.nCols ());
      returnValue = 
         DMatrix^(DMatrix.transpose ()^(in^DMatrix))^DMatrix.transpose ();
      return returnValue;
   } else {
      SxVector<Double> returnValue(in.getSize ());
      returnValue = 
         DMatrix^(DMatrix.transpose ()^in);
      return returnValue;
   }
  } else {
   //  cout << "Artifact Filter has no effect !" << endl;
     return in;
  }
}   

void SxArtifactFilter::setFreezeMode (const SxString &freezeModeIn)
{
   freezeMode = freezeModeIn;
   //needsUpdateDM = true;
}

void SxArtifactFilter::setStructure (const SxAtomicStructure &structureIn)
{
   tau.copy (structureIn);
   //needsUpdateDM = true;
}

void SxArtifactFilter::setMassWeighting (bool isActive)
{
   massWeighting = isActive;
   //needsUpdateDM = true;
}

void SxArtifactFilter::optimizeDMatrix (const SxMatrix<Double> &/*m*/) 
{
   SX_EXIT;
   /*
   cout << "Optimizing D-Matrix ..." << endl;

   int freeRots = 0;
   if (freezeMode == SxString("z")) freeRots = 1;
   if (freezeMode == SxString("xyz")) freeRots = 3;
   if (freezeMode == SxString("trans"))  freeRots = 0;
   */

}
   
   
   
SxMatrix<Double> SxArtifactFilter::getDMatrix () 
{
      int i, j, k;
      int freeRots = 0;
      int size = tau.getNAtoms ()*3;
      SxVector<Double> tauVec(size);;
      tauVec.copy (tau.coordRef ());
      SxVector3<Double> COM;
      SxMatrix<Double> inertiaTensor (3, 3);
      SxMatrix<Double> D (size, size);
      SxMatrix<Double>::Eigensystem IEig;	
      
      SxVector<Double> masses (size);
      SxVector<Double> column (size);
      SxVector<Double> orth(size);
    	
      masses = getMasses ();

      if (freezeMode == SxString("z")) freeRots = 1;
      if (freezeMode == SxString("xyz")) freeRots = 3;
      if (freezeMode == SxString("trans"))  freeRots = 0;
      
      D.set (0.);
      
      if (freezeMode != SxString ("notr")) {	
         
         // --- columns of the D-Matrix according to translational degrees 
         //     of freedom
         for (j = 0; j < 3; j++) {
            column.copy (D.colRef(j));
            for (i = 0; i < size/3; i++) column(i*3 + j) = sqrt (masses(i*3));
            D.colRef(j).set(column);
         }
         
         if (freezeMode == SxString ("xyz") || freezeMode == SxString ("z")) {
			
            // --- structure needs to be shifted to center of mass 
            //     to calculate inertia tensor
            COM = translateToCOM(&tauVec);
            inertiaTensor = getInertiaTensor(tauVec);
            IEig = inertiaTensor.eigensystem();
	
            SxVector3<Double> R;
            SxVector3<Double> IV;
			
            int vp1, vp2;

            // --- if freezing only around the z-axis is used to generate 
            //     rotational movements (see below) for xyz the eigenvectors
            //     of the inertia tensor are used
            
            if (freezeMode == SxString ("z")) {
               IV (0) = 0.; IV(1) = 0.; IV(2) = 0.;
               IEig.vecs.colRef(0).set (IV);
               IV (0) = 0.; IV(1) = 0.; IV(2) = 0.;
               IEig.vecs.colRef(1).set (IV);
               IV (0) = 1.0; IV(1) = 0.; IV(2) = 0.;
               IEig.vecs.colRef(2).set (IV);
            }
	
            // --- the columns of D according to rotational movements are 
            //     generated
            for (k = 0; k < freeRots; k++) {
               column.copy(D.colRef(k + 3));
               for (i = 0; i < size/3; i++) {
                  R(0) = tauVec(i*3);
                  R(1) = tauVec(i*3 + 1);
                  R(2) = tauVec(i*3 + 2);
	   
                  IV(0) = IEig.vecs.colRef(0)(k);
                  IV(1) = IEig.vecs.colRef(1)(k);
                  IV(2) = IEig.vecs.colRef(2)(k);
                  for (j = 0; j < 3; j++) {
                     vp1 = j + 1;
                     vp2 = vp1 + 1;
                     vp1 = vp1 - (vp1/3)*3;
                     vp2 = vp2 - (vp2/3)*3;
					
                     column(i*3 + j) 
                        = (R(vp1)*IV(vp2) - R(vp2)*IV(vp1))*sqrt(masses(3*i));
                  }
               }
               D.colRef(3 + k).set (column);
            }
         }
        
         // --- Normalisation of the translational and rotational columns
         //     of D
         for (i = 0; i < (3 + freeRots); i++) {
            column.copy (D.colRef (i));
            column = column/sqrt(column.absSqr().sum());
            D.colRef(i).set (column);
         }
         
         //--- gram-schmidt orthonormalisation yields the other 
         //    columns of D, which are then purely internal movements
         SxMatrix<Double> id(size, size);
         int shift = 1;
         id = id.identity ();
	
         for (i = 1; i < (3 + freeRots); i++) {
            column.set (D.colRef(i));
            for (j = 0; j < i; j++) {
               orth.set (D.colRef(j));
               column = column - (column^orth).chop()*orth;
               column = column/sqrt(column.absSqr().sum());
            }
            D.colRef(i).set (column);
         }
	
         for (i = 0; i < (size - (3 + freeRots)); i++) {
            column.copy (id.colRef(i + shift));
            for (j = 0; j < (3 + freeRots + i); j++) {
               column = column - (column^D.colRef(j)).chop()*D.colRef(j);
               if (sqrt(column.absSqr().sum()) < 1.e-1) {
                  i --;
                  j = 3 + freeRots + i;
                  shift += 1;
               } else
                  column = column/sqrt(column.absSqr().sum());
            }
            D.colRef(i + 3 + freeRots).set(column);
         }
			
         for (i = 0; i < size; i++) {
            column.copy (D.colRef (i));
            column = column/sqrt(column.absSqr().sum());
            D.colRef(i).set (column);
         }

         for (i = 1; i < size; i++) {
            column.copy (D.colRef(i));
            for (j = 0; j < i; j++) {
               orth.copy (D.colRef(j));
               column = column - (column^orth).chop()*orth;
               column = column/sqrt(column.absSqr().sum());
            }
            D.colRef(i).set(column);
         }
         
         //--- the D matrix further used is build up by the 
         //    the columns which belong to purely internal 
         //    movement  DMatrix.colRef(i - (3+freeRots)).copy (D.colRef(i));

         
         SxMatrix<Double> D_ = SxMatrix<Double> (size, size - (3 + freeRots));
         for (i = 3 + freeRots; i < size ; i++)  
            D_.colRef(i - (3+freeRots)).set (D.colRef(i));
         return D_;
   } else {
      SxMatrix<Double> D_ = SxMatrix<Double> (size, size - (3 + freeRots));
      D_.set(0.);
      return D_;
   }
}

SxVector<Double> SxArtifactFilter::getMasses ()
{
   int i, size;
   size = tau.getNAtoms () * 3;
   SxVector<Double> masses(size);
   if (massWeighting) masses.copy (massVec); 
   else {for (i = 0; i < size; i++)  masses(i) = 1.;}
   return masses;
}
   
void SxArtifactFilter::setMasses (const SxVector<Double> &massesIn) 
{
   massVec.resize (massesIn.getSize ());
   massVec.copy (massesIn);
}

void SxArtifactFilter::setMasses (const SxSpeciesData &speciesData) 
{
   int is, ia, counter, i;
   massVec = SxVector<Double>(tau.getNAtoms () * 3);
   counter = 0;
	for (is = 0; is < tau.getNSpecies (); is++) {
		for (ia = 0; ia < tau.getNAtoms (is); ia++) {
			for (i = 0; i < 3; i++) {
				massVec(counter) = speciesData.ionicMass(is);
				counter++;
			}
		}
	}

}

SxVector3<Double> SxArtifactFilter::translateToCOM (SxVector<Double> *tauVec) 
{
 	int i, size;
	size = (int)(*tauVec).getSize ();
	SxVector3<Double> COM;
   SxVector<Double> masses =  getMasses ();
	double M = 0.;
	
	COM(0) = COM(1) = COM (2) = M +0.;
   for (i = 0; i < size;) {
      M += masses(i);
      COM(0) += (*tauVec)(i)*masses(i); i++;
      COM(1) += (*tauVec)(i)*masses(i); i++;
      COM(2) += (*tauVec)(i)*masses(i); i++;
   }
   COM (0) /= M; COM (1) /= M; COM (2) /= M;
   for (i = 0; i < size;) {
      (*tauVec)(i) -= COM(0);	i++;
      (*tauVec)(i) -= COM(1);	i++;
      (*tauVec)(i) -= COM(2);	i++;
   }
   return COM;
}

SxMatrix<Double> 
    SxArtifactFilter::getInertiaTensor (const SxVector<Double> &tauVec)
{
	int i;
	int size = tau.getNAtoms () * 3;
	SxMatrix<Double> Ten(3, 3);
   SxVector<Double> masses = getMasses ();
   
	Ten.set (0.);
	for (i = 0; i < (size/3.); i++) {
		Ten(0, 0) += 
			masses(3*i)*
			( tauVec(3*i + 1)*tauVec(3*i + 1) 
			+ tauVec(3*i + 2)*tauVec(3*i + 2));
		Ten(1, 1) += 
			masses(3*i)*
			( tauVec(3*i + 0)*tauVec(3*i + 0) 
			+ tauVec(3*i + 2)*tauVec(3*i + 2));
		
		Ten(2, 2) += 
			masses(3*i)*
			( tauVec(3*i + 0)*tauVec(3*i + 0) 
			+ tauVec(3*i + 1)*tauVec(3*i + 1));

		Ten(0, 1) -= 
			masses(3*i)*
			( tauVec(3*i + 0)*tauVec(3*i + 1));
		Ten(1, 0) = Ten(0, 1);
		
		Ten(0, 2) -= 
			masses(3*i)*
			(tauVec(3*i + 0)*tauVec(3*i + 2));
		Ten(2, 0) = Ten(0, 2);
		
		Ten(1, 2) -= 
			masses(3*i)*
			(tauVec(3*i + 1)*tauVec(3*i + 2));
		Ten(2, 1) = Ten(1, 2);
	}

	return Ten;
}

void SxArtifactFilter::validate (const SxVector<Int> &/*equivalentIdx*/)
{
   // empty
}

