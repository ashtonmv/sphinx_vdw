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

#include <SxTaylorExpPotential.h>
#include <fstream>

SxTaylorExpPotential::SxTaylorExpPotential ()
{  //SVN_HEAD; // this line MUST NOT be removed!!!
   sxprintf ("This is TaylorExpPotetial\n");
}




SxTaylorExpPotential::SxTaylorExpPotential (const SxAtomicStructure &str,
      const SxSymbolTable *table) {
    speciesData = SxSpeciesData (&*table);
    SxBinIO io;
    int size;
    contributions = 0;
    nAtoms = str.getNAtoms();

    SxSymbolTable *cmd;
    cmd = table -> getGroup("taylorExpPotential");
    
   //--- only squared matrices allowed for reading in via parser!
   hessianGroup = cmd -> getGroup ("equTau");
   if (hessianGroup -> contains ("file")) {
      SxString fn = hessianGroup -> get("file") -> toString ();
      equTau = loadStructureFHI98(fn);
   } 

   else equTau = SxAtomicStructure(cmd -> getGroup("equTau"));
   size = equTau.nTlAtoms * 3;
   
   hessianGroup = cmd -> getGroup ("hessian") ;
   contributions += 1;
   if (hessianGroup -> contains ("file") ) {
      hessian.resize (size*size);
      hessian.reshape(size, size);
      SxString fName = hessianGroup -> get("file") -> toString ();
      io.open (fName, SxBinIO::ASCII_READ_ONLY);
      io.read ("dummy", &hessian, 
            size, size, 0, 0);
      io.close();
   }

   else {
      hessian = 
         SxMatrix<Double> (hessianGroup -> get("matrix") -> toList ());
      size =  int(lround(sqrt((double)hessian.getSize ())));
      hessian.reshape (size, size);
   }
   symmetrizeHessian ();
   if ( cmd -> containsGroup ("ttensor")) {
      contributions += 2;
      ttensor.resize (size);
      ifstream filestr;
      filestr.open((cmd -> getGroup("ttensor") -> get("file") -> toString()).ascii());
      int i,j,k;
      for (i=0;i<size;i++) {
         ttensor(i).resize(size);
         for (j=0; j<size; j++) {
            ttensor(i)(j).resize(size);
            for (k=0; k<size; k++) {
               filestr >> ttensor(i)(j)(k);
            }
         }
      }
      filestr.close();      
   }

 	getForces(str,table);
   double disp=0.01;
 exportHesse(getNumericalHesseMatrix(str,disp),"HesseMatrix");
   exportTTensor(getNumericalTTensor(str, disp),"TTensor");
}

  


SxTaylorExpPotential::SxTaylorExpPotential (const SxSymbolTable *cmd)
{
   contributions = 0;
   //nAtoms = str.getNAtoms();
   SxBinIO io;
   int size;
   //--- only squared matrices allowed for reading in via parser!
   hessianGroup = cmd -> getGroup ("equTau");
   if (hessianGroup -> contains ("file")) {
      SxString fn = hessianGroup -> get("file") -> toString ();
      equTau = loadStructureFHI98(fn);
   }

   else equTau = SxAtomicStructure(cmd -> getGroup("equTau"));
   size = equTau.nTlAtoms * 3;
   
   hessianGroup = cmd -> getGroup ("hessian") ;
   contributions += 1;
   if (hessianGroup -> contains ("file") ) {
      hessian.resize (size*size);
      hessian.reshape(size, size);
      SxString fName = hessianGroup -> get("file") -> toString ();
      io.open (fName, SxBinIO::ASCII_READ_ONLY);
      io.read ("dummy", &hessian, 
            size, size, 0, 0);
      io.close();
   }

   else {
      hessian = 
         SxMatrix<Double> (hessianGroup -> get("matrix") -> toList ());
      size =  int(round(sqrt((double)hessian.getSize ())));
      hessian.reshape (size, size);
   }
   symmetrizeHessian ();

   if ( cmd -> containsGroup ("ttensor")) {
      contributions += 2;
      ttensor.resize (size);
      ifstream filestr;
      filestr.open((cmd -> getGroup("ttensor") -> get("file") -> toString()).ascii());
      int i,j,k;
      for (i=0;i<size;i++) {
         ttensor(i).resize(size);
         for (j=0; j<size; j++) {
            ttensor(i)(j).resize(size);
            for (k=0; k<size; k++) {
               filestr >> ttensor(i)(j)(k);
            }
         }
      }
      filestr.close();      
   }
  }

SxAtomicStructure SxTaylorExpPotential::loadStructureFHI98 (const SxString &fileName)
{
   SxBinIO io;
   SxList<SxList<SxVector3<Double> > > tauList;
   SxList<SxString> speciesNameList;
   SxMatrix3<Double> cell;
   int i, j;
   SxAtomicStructure tau;
   
   io.open (fileName, SxBinIO::ASCII_READ_ONLY);
   tauList = io.loadStructureFHI98 ();
   io.close ();
   io.open (fileName, SxBinIO::ASCII_READ_ONLY);
   speciesNameList = io.loadSpeciesFHI98 ();
   io.close ();
   io.open (fileName, SxBinIO::ASCII_READ_ONLY);
   cell = io.loadCellFHI98 ();
   io.close ();

   tau = SxAtomicStructure (cell);
   tau.startCreation ();
   
   //Should be iterators (but not a bottleneck here)
   for (j = 0; j < tauList.getSize (); j++) {
      for (i = 0; i < tauList(j).getSize (); i++) {
         tau.addAtom (tauList(j)(i));
      }
   }
   tau.endCreation ();
   return tau;
}  

SxTaylorExpPotential::~SxTaylorExpPotential ()
{
   // empty
}


bool SxTaylorExpPotential::isRegistered (const SxSymbolTable *) const
{
   return false;
}

void SxTaylorExpPotential::execute (const SxSymbolTable *, bool /*calc*/)
{
   cout << "This is TaylorExpPotential." << endl;
   return;
}

double SxTaylorExpPotential::getEnergy () const
{
   return totalEnergy;
}

SxSpeciesData SxTaylorExpPotential::getSpeciesData () const
{
   return speciesData;
}

SxAtomicStructure SxTaylorExpPotential::getForces 
(const SxAtomicStructure &tauIn, const SxSymbolTable *)
{  //SVN_HEAD; // this line MUST NOT be removed
   double epsilon=1e-10;
   SxAtomicStructure tau;
   tau.copy (tauIn);
   SxVector<Double> delta = -(tau.coordRef () - equTau.coordRef ());
   SxAtomicStructure force (tau);
   force.set (hessian^delta);
   totalEnergy = 0.5*(delta.transpose() ^ force.coordRef()).chop ();

   if ( containsContribution ("ttensor")) {
      int dof = 3*nAtoms;
      int i, j, k;
      double energyFromTTensor = 0.;
		SxVector<Double> forceFromTTensor(dof), x(dof);
      forceFromTTensor.set (0.);
      x = -delta;

      SxVector<Int> ind1(2), ind2(3);
      for (i = 0; i < dof; i++) {
         if (fabs(x(i)) > epsilon) {
            ind1(0) = i;
            ind2(0) = i;
            for (j = i; j < dof; j++) {
               if (fabs (x(j)) > epsilon) {
                  ind1(1) = j;
                  ind2(1) = j;
                  forceFromTTensor = forceFromTTensor +
                                     e(ind1) * x(i)*x(j) * ttensor(i)(j);
                  for (k = j; k < dof; k++) {
                     if (fabs(x(k)) > epsilon) {
                        
                        ind2(2)=k;
                        energyFromTTensor = energyFromTTensor + 
                                   e(ind2) * x(i)*x(j)*x(k) * ttensor(i)(j)(k);
               
 //               cout << i << " " << j << " " << k << endl;
 //               cout << x(i) << " " << x(j) << " " << x(k) << endl;
 //               cout << e(ind2) << endl;
 //               cout << ttensor(i)(j)(k) << endl;
 //               cout << e(ind2) * x(i)*x(j)*x(k) * ttensor(i)(j)(k) << endl << endl;


                     }
                  }
               }
            }
         }
      }
 //     cout << "energyFromTTensor " << energyFromTTensor << endl;
      totalEnergy += energyFromTTensor;
      for (i=0; i<equTau.nTlAtoms; i++) {
			for (j=0; j<3; j++) {
				force.ref(i)(j) += forceFromTTensor(3*i+j);
         }
      }
   }

// cout << "TOTAL ENERGY: " << totalEnergy << endl;
// cout << "FORCES: " << endl << force << endl;

   return force;
}

SxArray<SxAtomicStructure> SxTaylorExpPotential::getNumericalHesseMatrix
   (const SxAtomicStructure &tauIn, const double &dx)
{
   SxArray<SxAtomicStructure> Hesse;
   Hesse.resize(3*nAtoms);
   SxAtomicStructure A, B;
	SxAtomicStructure devCoord; devCoord.copy    (tauIn); 
   SxSymbolTable *table = NULL;

	for (int i = 0; i < nAtoms; i++) {
		for (int j = 0; j < 3; j++) {
			devCoord.ref(i)(j) += dx;
			A = getForces (devCoord,table);
			devCoord.ref(i)(j) -= 2.*dx;
			B = getForces (devCoord,table);
         Hesse(3*i+j) = (1./dx/2.)*(A-B);
         devCoord.ref(i)(j) += dx;
		}
	}
	return Hesse;
}

SxArray<SxArray<SxAtomicStructure> > SxTaylorExpPotential::getNumericalTTensor 
   (const SxAtomicStructure &tauIn, const double &dx) 
{
   double threshold = 1e-9;

   SxArray<SxArray<SxAtomicStructure> > TTensor;
   TTensor.resize(3*nAtoms);
   SxArray<SxAtomicStructure> A, B, Delta;
   SxAtomicStructure devCoord; devCoord.copy    (tauIn);
   int i,j,k,l,m;

   for (i = 0; i < nAtoms; i++) {
      cout << "atom " << i+1 << endl;
   	for (j = 0; j < 3; j++) {
   		devCoord.ref(i)(j) += dx;
   		A = getNumericalHesseMatrix (devCoord,dx);
   		devCoord.ref(i)(j) -= 2.*dx;
   		B = getNumericalHesseMatrix (devCoord,dx);
         Delta = A-B;
         for (k=0;k<Delta.getSize();k++) {
   			Delta(k) = (1./dx/2.) * Delta(k);
            for (l=0; l<nAtoms; l++) {
               for (m=0; m<3; m++) {
                  if (Delta(k).ref(l)(m)<threshold) Delta(k).ref(l)(m)=0.;
               }
            }
         }
         
         TTensor(3*i+j) = Delta;
         devCoord.ref(i)(j) += dx;
   	}
   }
  	return TTensor;
}

void SxTaylorExpPotential::symmetrizeHessian () 
{
   int i, j, size;
	size = (int)hessian.nRows ();
   double sym;

   for (i = 0; i < size; i++) {
      for (j =0; j < size; j++) {
         sym = (hessian(i, j) + hessian(j, i))/2.;
         hessian(i, j) = sym;
         hessian(j, i) = sym;
      }
   }
}

bool SxTaylorExpPotential::containsContribution (const SxString &contr)
{
   int c;
	if (contr=="hessian") c=1;
	if (contr=="ttensor") c=2;
	if (contr=="quadtensor") c=4;
	if (contr=="quintensor") c=8;

   double fracpart, intpart;
   modf (contributions/c,&intpart);
   fracpart = modf (intpart/2,&intpart);
   if ( fracpart > .25 ) return true; 
                    else return false;
}

void SxTaylorExpPotential::exportHesse
   (const SxArray<SxAtomicStructure> &hesse, const SxString &fileName)
{
   int i, j;
   ofstream filestr;
   filestr.open(fileName.ascii(), ifstream::trunc);
   for (i=0; i<3*nAtoms; i++) {
      for (j=0; j<nAtoms; j++) {
         filestr << -hesse(i)(j)(0) << " " << -hesse(i)(j)(1) << " " << -hesse(i)(j)(2) << " ";
      }
      filestr << endl;
   }
   filestr.close();
}

void SxTaylorExpPotential::exportTTensor 
   (const SxArray<SxArray<SxAtomicStructure> > &TTensor, 
    const SxString &fileName)
{
   int i, j, k;
   ofstream filestr;
   filestr.open(fileName.ascii(), ifstream::trunc);
	for (k=0;k<3*nAtoms;k++) {
      for (i=0; i<3*nAtoms; i++) {
         for (j=0; j<nAtoms; j++) {
            filestr << TTensor(k)(i)(j)(0) << " " << TTensor(k)(i)(j)(1) 
               << " " << TTensor(k)(i)(j)(2) << " ";
         }
      }
      filestr << endl;
   }
   filestr.close();
}

double SxTaylorExpPotential::e (const SxVector<Int> &v)
{
   SxList<int> count, unionList;
   int i,j;
   bool isIn;
   
   for (i=0; i<v.getSize(); i++) {
      isIn=false;
		for (j=0; j<unionList.getSize(); j++) {
			if (v(i)==unionList(j)) {
				count(j)++;
            isIn=true;
            break;
         }
      }
      if (isIn==false) {
			unionList.append(v(i));
         count.append(1);
      }
   }

   double fac=1.;
   for (i=0; i<count.getSize(); i++) {
		fac *= 1./(factorial(count(i)));
   }

   return fac;
}

int SxTaylorExpPotential::factorial(int i) {
    if (i <= 1)
        return i;
    return (i * factorial(i - 1));
}
