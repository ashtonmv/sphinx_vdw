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
// Author: Lars Ismer, ismer@fhi-berlin.mpg.de

#include <SxVDW.h>

SxVDW::SxVDW () 
{ //SVN_HEAD;
   
}

SxVDW::SxVDW 
(const SxAtomicStructure &t, const SxSymbolTable *cmd, 
 const SxSpeciesData &speciesData) 
{ 
   int i;
   SxAtomicStructure tau;
   tau.copy(t);

   nAtoms = tau.nTlAtoms;
   
	superCell             = tau.cell;
   //---workaround: will be solved if connected to sxatomicstructure
   if (superCell.absSqr().sum() < 1e-10) { 
      superCell (0, 0) = superCell(1, 1) = superCell(2, 2) = 100;
   }
      
	coord                 = SxArray<SxVector3<Double> > (nAtoms);
	species               = SxArray<SxString> (nAtoms);

   //SxSpeciesData speciesData = SxSpeciesData(cmd -> topLevel ());
   
   for (i = 0; i < nAtoms; i++) {
      coord(i) = tau.ref (i);
      species(i) = speciesData.chemName(tau.getISpecies (i));
   }
      
	energyContrib   = SxArray<double>             (nAtoms);
   
	
	dist                  = SxArray<SxList<double> > (nAtoms);
	expTerm               = SxArray<SxList<double> > (nAtoms);
	neighbours            = SxArray<SxList<int> >    (nAtoms);
	supercellIds          = SxArray<SxList<int> >    (nAtoms);
	sNeighbours            = SxArray<SxList<int> >    (nAtoms);
	sSupercellIds          = SxArray<SxList<int> >    (nAtoms);

	borderCrossers            = SxArray<SxList<int> >    (27);
	output = false;

	Forces = SxArray<SxVector3<Double> > (nAtoms);

   potentialType = cmd -> getGroup ("vdwCorrection") 
                       -> get("potentialType") -> toString ();
      
}

SxVDW::~SxVDW ()
{
   // empty
}

void SxVDW::resize 
(SxList<SxList<SxVector3<Double> > >  tau, SxList<SxString> speciesNameList, SxMatrix3<Double> aMat, SxString pt)
{
	int i, j;
	int counter;
  	
	potentialType = pt;
	nAtoms = 0;
	for (i = 0; i < tau.getSize (); i++) {
			nAtoms += (int)tau(i).getSize ();
	}

	superCell             = aMat;
	
	coord                 = SxArray<SxVector3<Double> > (nAtoms);
	species               = SxArray<SxString> (nAtoms);
	
	counter = 0;
	for (i = 0; i < tau.getSize (); i++) {
		for (j = 0; j < tau(i).getSize (); j++) {
			species(counter) = speciesNameList(i);
			counter++;
		}
	}

	energyContrib   = SxArray<double>             (nAtoms);
	
	dist                  = SxArray<SxList<double> > (nAtoms);
	expTerm               = SxArray<SxList<double> > (nAtoms);
	neighbours            = SxArray<SxList<int> >    (nAtoms);
	supercellIds          = SxArray<SxList<int> >    (nAtoms);
	sNeighbours            = SxArray<SxList<int> >    (nAtoms);
	sSupercellIds          = SxArray<SxList<int> >    (nAtoms);

	borderCrossers            = SxArray<SxList<int> >    (27);
	output = false;

	Forces = SxArray<SxVector3<Double> > (nAtoms);

}

void SxVDW::updateHybridisation ()
{
	int i;
	for (i = 0; i < nAtoms; i++) {
		if ( (species(i) == SxString ("C")) 
			  || (species(i).contains(SxString("C_")))) {
			if (neighbours(i).getSize() == 3) {
				species(i) = SxString("C_sp2");
			}
			if (neighbours(i).getSize() == 4) {
				species(i) = SxString("C_sp3");
			}
         if ( (neighbours(i).getSize () < 3) || 
              (neighbours(i).getSize () > 4) ) {
         cout << "NgVDW: WARNING !!!!!!!: The Hybridisation of the "
              << "Carbon Atom with id " << i << " is neither sp2 nor sp3 "
              << "is treated as sp3" << endl;
         cout << "Check atomic basis and supercell geometry !!!" << endl;
         species(i) = SxString("C_sp3");
         }

		}
	}
}

void SxVDW::update (SxArray<SxVector3<Double > > newcoord) 
{
	if (output) cout << " ---- Updating Coords" << endl;
	updateCoord (newcoord);
	//printParameters();
	//printCoord();
	//printSpecies ();
	if (output) cout << " ---- Updating Border Crossers" << endl;
	updateBorderCrossers();
	if (output) cout << " ---- Updating Neighbourhood and Distances" << endl;
	updateNeighbours (SxString ("covalentCutoff"));
	//printNeighbours ();
	updateHybridisation ();
	//printSpecies ();
	updateNeighbours (SxString ("vanDerWaalsCutoff"));
	//printNeighbours();
	//updateSecondNeighbours();
}

SxVector3<Double> SxVDW::getLatticeVec(int supercellId) 
{
	int a1Coeff, a2Coeff, a3Coeff;
	int help = supercellId/9;
   a1Coeff = a2Coeff = a3Coeff = 0; 
	
	if (help == 0) a1Coeff = 0; 
 	if (help == 1) a1Coeff = 1;
	if (help == 2) a1Coeff = -1;

	supercellId = supercellId - help*9;

  	help = supercellId/3;
	
	if (help == 0) a2Coeff = 0; 
 	if (help == 1) a2Coeff = 1;
	if (help == 2) a2Coeff = -1;

	supercellId = supercellId - help*3;

  	help = supercellId;
	
	if (help == 0) a3Coeff = 0; 
 	if (help == 1) a3Coeff = 1;
	if (help == 2) a3Coeff = -1;

	//cout << "Supercell: " << a1Coeff 
	//	                   << " " << a2Coeff << " " << a3Coeff << endl; 

	SxVector3 <Double> returnValue = a1Coeff*superCell(0) + a2Coeff*superCell(1) 
		     + a3Coeff*superCell(2);
   //cout << returnValue << endl;
	return returnValue;
		     
}

double SxVDW::getDist (int i, int j, int supercellId) {
	SxVector3<Double> distVec;
	
	distVec = coord(i) - coord(j) - getLatticeVec(supercellId);
	
	return 
		sqrt (distVec(0)*distVec(0) + distVec(1)*distVec(1) + distVec(2)*distVec(2));
}
		
void SxVDW::updateCoord (SxArray<SxVector3<Double> > newcoord) 
{
	
	for ( int i = 0; 
		   i < nAtoms; 
			i++) {
		coord(i) = newcoord(i);
	}
}

	    	  
void SxVDW::updateBorderCrossers () 
{
	double maxDist = getParam("a", 0, 0);
	//double length, projection;
	int j, i, superCellId;
	//int a0Flip, a1Flip, a2Flip;
	//SxArray<int> acoeff(3);

	// comment 20.02.2003: not so sure about the following comment
	// only works for orthorombic supercells, but easy to generalize


	bool smallSuperCell = false;
	smallSuperCell = true;
	
	//-- should be placed elsewhere (Initialisiation)
	for (j = 0; j < 3; j++) {
		if (sqrt(superCell(j).absSqr().sum ()) < 2.*maxDist) smallSuperCell = true;
	}
	
	for (superCellId = 0; superCellId < 27; superCellId++) {
		borderCrossers(superCellId).resize(0);
	}
	
	for (i = 0; i < nAtoms; i++) {
	if (smallSuperCell) {
		for (superCellId = 0; superCellId < 27; superCellId++) {
			borderCrossers(superCellId).append(i);
	}

	} else {
   //--- did not work, so commented out (has importance for efficience)
		/*
		acoeff(0) = acoeff(1) = acoeff(2) = 0;
		for (j = 0; j < 3; j++) {
			lenth = sqrt(superCell(j).absSqr ());
			projection = (  superCell(j)(0)*coord(i)(0)
					       +  superCell(j)(1)*coord(i)(1)
							 +  superCell(j)(2)*coord(i)(2) )
				         /length;
			if (projection < maxDist) {
				acoeff(j) = 1;
			} else {
			if ((length - projection) < maxDist) {
				acoeff(j) = -1;
			} else {
				acoeff(j) = 0;
			}
			}
		}
		
		
		for (a0Flip = 0; abs(a0Flip) < 2 ; a0Flip += acoeff(0)) {
			for (a1Flip = 0; abs(a1Flip) < 2 ; a1Flip += acoeff(1)) {
				for (a2Flip = 0; abs(a2Flip) < 2; a2Flip += acoeff(2)) {

					//cout << a0Flip << " " << a1Flip << " " << a2Flip << endl;
					superCellId = 9*(a0Flip*(3*a0Flip - 1)/2)
						+ 3*(a1Flip*(3*a1Flip - 1)/2)
						+ (a2Flip*(3*a2Flip - 1)/2);
					borderCrossers(superCellId).append (i);
					
					if (acoeff(2) == 0) a2Flip = 2; 
				}
				if (acoeff(1) == 0) a1Flip = 2; 
			}
			if (acoeff(0) == 0) a0Flip = 2; 
		}
		
    */
	}
	
	}
}
  	

void SxVDW::updateNeighbours (SxString modus) 
{
   double distance;
	SxList<int>::Iterator itBC;
	neighbours   = SxArray<SxList<int> >    (nAtoms);
	dist         = SxArray<SxList<double> > (nAtoms);
	supercellIds = SxArray<SxList<int> >    (nAtoms);
	int i, j, supercells;
	/*
	for (supercells = 0; supercells < 27; supercells++) {
		cout << endl << "Border Crossers in cell" << getLatticeVec (supercells) <<
			endl;
		for ( itBC  = borderCrossers (supercells).begin();
					itBC != borderCrossers (supercells).end();
					++itBC) {
				j = *itBC;
					cout << j << " ";
		}
		cout << endl;
	}
*/
	
	for (i = 0; i < nAtoms; i++) {
		dist(i).         resize (0);
		neighbours(i).   resize (0);
		supercellIds(i). resize (0);
		
	for (supercells = 0   ; supercells < 27; supercells++) {
		if (supercells == 0) {
			for (j = 0; j < nAtoms; j++) {
				if ( ((i != j)) && ((distance = getDist(i, j, supercells)) < getParam (modus, i , j))) {
					neighbours(i).append (j);
					dist(i).append (distance);
					supercellIds(i).append (supercells);
				} 
			}
		} else {
			for ( itBC  = borderCrossers (supercells).begin();
					itBC != borderCrossers (supercells).end();
					++itBC) {
				j = *itBC;
				if ((distance = getDist(i, j, supercells)) < 
						getParam (modus, i , j)) { 
					neighbours(i).append (j);
					dist(i).append (distance);
					supercellIds(i).append (supercells);
				}

			}
		}
	}
	}
}
 
double SxVDW::getDampingFunction (double R, double Rm) 
{
	double fd = 0.;
	double cdamp = getParam("cdamp", 0, 0);
	double beta = getParam("beta", 0, 0);
   double dstar = getParam("dstar", 0, 0);

	if (potentialType == SxString("WTYang_I")) 
		fd = ::pow ((1 - exp(-cdamp*::pow((R/Rm), 3.))), 2.);
	if (potentialType == SxString("WTYang_II")) 
		fd = (1/(1 + exp(-beta*(R/Rm - 1))));
   if (potentialType.contains("Elstner"))
      fd = ::pow((1 - exp(-dstar*::pow((R/Rm), 7.))), 4);

	return fd;
}

double SxVDW::getDampingDerivative (double R, double Rm) 
{
	double fd = 0.;
	double e = 0.;

	double cdamp = getParam("cdamp", 0, 0);
	double beta = getParam("beta", 0, 0);
   double dstar = getParam("dstar", 0, 0);
	
	if (potentialType == SxString("WTYang_I")) {
		e = exp(-cdamp*(::pow( (R/Rm), 3.)));
		fd = 6.*cdamp*R*R/Rm/Rm/Rm*(1. - e)*e;
	}
	
	if (potentialType == SxString("WTYang_II")) {
		e = exp(-beta*(R/Rm - 1.));
		fd = beta/Rm*e/(1. + e)/(1. + e);
	}
	
   if (potentialType.contains("Elstner")) {
		e = exp(-dstar*(::pow( (R/Rm), 7.)));
		fd = ::pow((1-e), 3.)*28.*dstar/::pow(Rm, 7.)*::pow(R, 6.)*e;
	}
	return fd;
}

double SxVDW::getDampingSecondDerivative (double R, double Rm) {
	double fd, e, ePrime, ePrimePrime, a;
	double cdamp = getParam("cdamp", 0, 0);
	double beta = getParam("beta", 0, 0);
	double dstar = getParam("dstar", 0, 0);
   fd = e = ePrime = ePrimePrime = a = 0.;
	
	if (potentialType == SxString("WTYang_I")) {
		e = exp(-cdamp*(::pow( (R/Rm), 3.)));
		ePrime = -e*cdamp*3.*R*R/Rm/Rm/Rm;
		fd = 12*cdamp*R/Rm/Rm/Rm*(1.-e)*e
			+ 6.*cdamp*R*R/Rm/Rm/Rm*(ePrime - 2.*e*ePrime);
	}
	
	if (potentialType == SxString("WTYang_II")) {
		e = exp(-beta*(R/Rm - 1.));
		fd = - ::pow((beta/Rm), 2.)*e/(1. + e)/(1. + e)*(1. - 2*e/(1 + e));
	}

   if (potentialType.contains("Elstner")) {
		a = dstar/::pow (Rm, 7.);
      e = exp(-dstar*(::pow( (R/Rm), 7.)));
      ePrime = -7.*a*::pow(R, 6.)*e;
      ePrimePrime = 49.*a*a*::pow(R, 12.)*e - 42*a*::pow(R, 5.)*e; 
      fd = -(4.*ePrimePrime * ::pow((1 - e), 3.) 
            - 12*ePrime*ePrime*(::pow((1-e), 2)));
	}
   
	return fd;
}
	
double SxVDW::getTotalEnergy () {
	double R, Rm, fd, eVDW, C6;
	int i, j, neighj;
	eVDW = 0.;
	for (i = 0; i < nAtoms; i++) {
		for (j = 0; j < neighbours(i).getSize (); j++) {
			R = dist(i)(j);
			neighj = neighbours(i)(j);
			Rm = getRm (i, neighj);
			C6 = getC6 (i, neighj);
			fd = getDampingFunction (R, Rm);
			/*
			cout << "R: " << R << endl;
			cout << "Rm: " << Rm << endl;
			cout << "C6: " << C6 << endl;
			cout << "fd: " << fd << endl;
			cout << "eVDW: " << eVDW << endl;
			*/
			eVDW += -fd*C6/(::pow(R, 6.))/2.;
		}
	}
	
	return eVDW;
					
}		

	

SxVector3<Double>  SxVDW::getForceOnAtom (int i) {
	
	SxVector3<Double>  returnValue;
	
   double R, Rm, C6, fd, fdPrime, derivative;
			 
	int neighj;
	
	returnValue(0) = 0.;
	returnValue(1) = 0.;
	returnValue(2) = 0.;

	
	//cout << "Atom "<< i << endl;

	for (int j = 0; j < neighbours(i).getSize (); j++) {
		//if (neighbours (i)(j)) {
	      neighj = neighbours(i)(j);	
			//cout << "Bonding Partner: " << neighj << endl;;
			R = dist(i)(j);
			Rm = getRm(i, neighj);
			C6 = getC6(i, neighj);
			fdPrime = getDampingDerivative (R, Rm);
			fd = getDampingFunction (R, Rm);

			derivative  
				= fdPrime*C6/::pow(R, 6.) - 6.*fd*C6/::pow(R, 7.); 

			returnValue += derivative * (coord(i) - coord(neighj) 
					       - getLatticeVec (supercellIds(i)(j)))
				            / dist(i)(j);
		   	
		//}	
	}
	
	return returnValue;
}

void SxVDW::updateForces () {
	for (int i = 0; i < nAtoms; i++) {
		for (int j = 0; j < 3; j++) {
		Forces(i)(j) = 0;
		}
	}
		

	for (int i = 0; i < nAtoms; i++) {
		Forces(i) = getForceOnAtom (i);
	}
}


SxArray<SxVector3<Double> > SxVDW::getForces () {
	
	totalEnergy = 0.;
	if (output) cout << "...2 Body Contributions..." << endl;
	updateForces ();
	return Forces;
}
			
bool SxVDW::areNeighbors (int atom1, int atom2) {
	int i;
   bool isNeighbor = false;
	
	for (i = 0; i < neighbours(atom1).getSize (); i++) {
	  if (atom2 == neighbours(atom1)(i)) {
		  isNeighbor = true;
		  i = (int)neighbours(atom1).getSize ();
	  } 
	}

	return isNeighbor;
}
	
int SxVDW::getNeighborIndex (int atom1, int atom2) {
	int i;
   int neighbor = 0;
	
	
	for (i = 0; i < neighbours(atom1).getSize (); i++) {
	  if (atom2 == neighbours(atom1)(i)) {
		  
		  neighbor = i;
		  i = (int)neighbours(atom1).getSize ();
	  } 
	}

	return neighbor;
}
SxMatrix3<Double> SxVDW::getInteraction (int atom1, int atom2, int neighborIndex) 
{
	SxMatrix3<Double> returnValue;
	SxMatrix3<Double> dtwor;
	SxMatrix3<Double> twodr;
	int i, j;
	
	SxVector3<Double> coord1, coord2, deltaCoord;
	coord1 = coord(atom1);
	coord2 = coord(atom2);
	double dg, g, fd, fdP, fdPP;
	
	double R = dist(atom1)(neighborIndex);
	double Rm = getRm (atom1, atom2);
	double C6 = getParam(SxString("C6"), atom1, atom2);

	
	deltaCoord = coord1 - coord2 
					       - getLatticeVec (supercellIds(atom1)(neighborIndex));

	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			if (i == j) {

				dtwor(i, j) = -1/R + ::pow((deltaCoord(i)), 2.)/::pow(R, 3.);
			} else {
				dtwor(i, j) = (deltaCoord(i))* ( deltaCoord(j))
					          /::pow(R, 3.);
			}
		}
	}

	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
				twodr(i, j) =  - (deltaCoord(i))
					             *(deltaCoord(j))/::pow(R, 2.);
		}
	}


	fd = getDampingFunction (R, Rm);
	fdP = getDampingDerivative (R, Rm);
	fdPP = getDampingSecondDerivative (R, Rm);

	g = fdP*C6/::pow(R, 6.) - 6.*fd*C6/::pow(R, 7.);
	dg = fdPP*C6/::pow(R, 6.) - 12.*fdP*C6/::pow(R, 7.) + 42.*fd*C6/::pow(R, 8.);

	returnValue =  twodr*dg + dtwor*g;
	//cout << g << endl;
	
	//cout << returnValue << endl << endl;

	return returnValue;
}

SxMatrix<Double> SxVDW::getHessian () {
	int i, j, k, l;
	int size = nAtoms*3;
	SxMatrix<Double> returnValue (size, size);
	SxMatrix3<Double> interaction;

	returnValue.set (0.);

	for (k = 0; k < 3; k++) {
		for (l = 0; l < 3; l++) {
			interaction(k)(l) = 0.;
		}
	}

	for (i = 0; i < nAtoms; i++) {
		for (j = 0; j < nAtoms; j++) {
			
			interaction = interaction*0.;
			if (i == j) {
				for (k = 0; k < neighbours(i).getSize (); k++) {
					if (i != neighbours(i)(k)) {
						interaction = 
							interaction + getInteraction(i, neighbours(i)(k), k);
					}
				}
			} else {
				if (areNeighbors (i, j)) {
					for (k = 0; k < neighbours(i).getSize (); k++) {
						if ( neighbours(i)(k) == j ){
							interaction = interaction - getInteraction(i, j, k);
							//interaction = (-1.)*interaction;
						}
					}
				} else {
					interaction = 0.*interaction;
				}
			}
		
			for (k = 0; k < 3; k++) {
				for (l = 0; l < 3; l++) {
					returnValue(3*i + k, 3*j +l) = interaction(k, l);
				}
			}
		}
	}

	return returnValue;
}
					

SxMatrix<Double>  SxVDW::getNumericalHessian (double dx) {
	int i, j, k, l;
	SxArray<SxVector3<Double> > undevCoord (nAtoms);
  	SxArray<SxVector3<Double> > devCoord   (nAtoms);
	SxArray<SxVector3<Double> > devForcesLeft(nAtoms);
	SxArray<SxVector3<Double> > devForcesRight(nAtoms);

	SxMatrix<Double> hessian (nAtoms*3, nAtoms*3);
   //update(coord);
	
	//undevForces = getForces();

	for (i = 0; i < nAtoms; i++) {
      undevCoord(i) = coord(i);
		devCoord(i) = coord(i);
	}

	for (i = 0; i < nAtoms; i++) {
		cout << "Atom " << i << " of " << nAtoms << endl;
		for (j = 0; j < 3; j++) {
			devCoord(i)(j) += dx;
			update(devCoord);
			devForcesRight = getForces ();
			devCoord(i)(j) -= 2.*dx;
			update(devCoord);
			devForcesLeft = getForces ();
			devCoord(i)(j) += dx;
			for (k = 0; k < nAtoms; k++) {
				for (l = 0; l < 3; l++) {
					hessian(i*3 + j, k*3 +l) 
						= -(devForcesRight(k)(l) - devForcesLeft(k)(l))/2./dx;
				}
			}
		}
	}
	update(undevCoord);
	return hessian;
}

SxArray<SxVector3<Double> > SxVDW::getNumericalForces (double dx) {
	SxArray<SxVector3<Double> > forces     (nAtoms);
	SxArray<SxVector3<Double> > dummy      (nAtoms); 
	SxArray<SxVector3<Double> > undevCoord (nAtoms);
  	SxArray<SxVector3<Double> > devCoord   (nAtoms);
	double undevEnergy;

   update(coord);
	
	totalEnergy = getTotalEnergy ();
	undevEnergy = totalEnergy;

	for (int i = 0; i < nAtoms; i++) {
      undevCoord(i) = coord(i);
		devCoord(i) = coord(i);
	}

	for (int i = 0; i < nAtoms; i++) {
		cout << "Atom " << i << " of " << nAtoms << endl;
		for (int j = 0; j < 3; j++) {
			forces(i)(j) = 0.;
			devCoord(i)(j) += dx;
			update(devCoord);
			totalEnergy = getTotalEnergy ();
			forces(i)(j) -= (totalEnergy - undevEnergy)/dx/2.;
			devCoord(i)(j) -= 2.*dx;
			update(devCoord);
			totalEnergy = getTotalEnergy();
			forces(i)(j) += (totalEnergy - undevEnergy)/dx/2.;
			devCoord(i)(j) += dx;
		}
	}
  	totalEnergy = undevEnergy;
	return forces;
	
}

double SxVDW::getRm (int atom1, int atom2)
{
	double Rm = 0.;
   if (potentialType.contains ("WTYang")) {
      if (species(atom1).contains(SxString("N"))) Rm += 2.93;
      if (species(atom1).contains(SxString("H"))) Rm += 2.27;
      if (species(atom1).contains(SxString("O"))) Rm += 2.87;
      if (species(atom1).contains(SxString("C"))) Rm += 3.21;
	
      if (species(atom2).contains(SxString("N"))) Rm += 2.93;
      if (species(atom2).contains(SxString("H"))) Rm += 2.27;
      if (species(atom2).contains(SxString("O"))) Rm += 2.87;
      if (species(atom2).contains(SxString("C"))) Rm += 3.21;
   }
   
   if (potentialType.contains ("Elstner")) {
       Rm = 3.8 * 1.8897;
   }

	return Rm;
}

double SxVDW::getC6 (int atom1, int atom2) 
{
	double C6one = 0., None = 0., C6two = 0., Ntwo = 0., C6 = 0.;
   double conv = 1./0.5976;
   double conv2 = 6.7481;
   double Pone = 0., Ptwo = 0.;
      
   if (potentialType.contains("WTYang")) {
   //--- C6 coefficients from Q. Wu and W. Yang, JCP 116, 515 (2002), Table II
      if (species(atom1).contains("C_sp2")) {
         C6one = 27.32;
         None = 2.02;
      }

      if (species(atom1).contains("C_sp3")) {
         C6one = 22.05;
         None = 2.02;
      }
	
      if (species(atom1).contains("H")) {
         C6one = 2.845;
         None = 0.53;
      }
	
      if (species(atom1).contains("O")) {
         C6one = 13.07;
         None = 2.65;
      }

      if (species(atom1).contains("N")) {
         C6one = 19.48;
         None = 2.52;
      }
	
      if (species(atom2).contains("C_sp2")) {
         C6two = 27.32;
         Ntwo = 2.02;
      }

      if (species(atom2).contains("C_sp3")) {
         C6two = 22.05;
         Ntwo = 2.02;
      }
	
      if (species(atom2).contains("H")) {
         C6two = 2.845;
         Ntwo = 0.53;
      }
	
      if (species(atom2).contains("O")) {
         C6two = 13.07;
         Ntwo = 2.65;
      }

      if (species(atom2).contains("N")) {
         C6two = 19.48;
         Ntwo = 2.52;
      }

	C6 = 2 * 
		::pow (C6one*C6one*C6two*C6two*None*Ntwo, 1./3.)/
		(::pow (C6one*Ntwo*Ntwo, 1./3.) + (::pow (C6two*None*None, 1./3.)));
   }

   if (potentialType.contains("Elstner_I")) {
   //--- C6 coefficients from M. Elstner et al., JCP 114, 5149 (2001), 
   //    Table II (first row), effective number of electrons (None)
   //    are calculated according to Eq. 6   
      
      
      if (species(atom1).contains("C_sp2")) {
         C6one = 18.56*conv;
         Pone = 1.352*conv2;
      }

      if (species(atom1).contains("C_sp3")) {
         C6one = 12.93*conv;
         Pone = 1.061*conv2;
      }
	
      if (species(atom1).contains("H")) {
         C6one = 1.61*conv;
         Pone = 0.387*conv2;
      }
	
      if (species(atom1).contains("O")) {
         C6one = 5.71*conv;
         Pone = 0.569*conv2;
      }

      if (species(atom1).contains("N")) {
         C6one = 13.16*conv;
         Pone = 1.030*conv2;
      }
      
      if (species(atom2).contains("C_sp2")) {
         C6two = 18.56*conv;
         Ptwo = 1.352*conv2;
      }

      if (species(atom2).contains("C_sp3")) {
         C6two = 12.93*conv;
         Ptwo = 1.061*conv2;
      }
	
      if (species(atom2).contains("H")) {
         C6two = 1.61*conv;
         Ptwo = 0.387*conv2;
      }
	
      if (species(atom2).contains("O")) {
         C6two = 5.71*conv;
         Ptwo = 0.569*conv2;
      }

      if (species(atom2).contains("N")) {
         C6two = 13.16*conv;
         Ptwo = 1.030*conv2;
      }

      //--- Slater-Kirkwood interpolation formula 
      //    Note: Elstner paper (Eq. 7) is wrong here
      //    Formula is correct in T. Halgren, JACS 114, 7827 (1992) 
      C6 = (2.*C6one*C6two*Pone*Ptwo)
         / (Ptwo*Ptwo*C6one + Pone*Pone*C6two);
   }
   if (potentialType.contains("Elstner_II")) {
   //--- C6 coefficients from M. Elstner et al., JCP 114, 5149 (2001), 
   //    Table II (second row), effective number of electrons (None)
   //    are calculated according to Eq. 6   
      
      
      if (species(atom1).contains("C_sp2")) {
         C6one = 16.07*conv;
         Pone = 1.352*conv2;
      }

      if (species(atom1).contains("C_sp3")) {
         C6one = 12.37*conv;
         Pone = 1.061*conv2;
      }
	
      if (species(atom1).contains("H")) {
         C6one = 1.53*conv;
         Pone = 0.387*conv2;
      }
	
      if (species(atom1).contains("O")) {
         C6one = 4.15*conv;
         Pone = 0.569*conv2;
      }

      if (species(atom1).contains("N")) {
         C6one = 11.55*conv;
         Pone = 1.030*conv2;
      }
      
      if (species(atom2).contains("C_sp2")) {
         C6two = 16.07*conv;
         Ptwo = 1.352*conv2;
      }

      if (species(atom2).contains("C_sp3")) {
         C6two = 12.37*conv;
         Ptwo = 1.061*conv2;
      }
	
      if (species(atom2).contains("H")) {
         C6two = 1.53*conv;
         Ptwo = 0.387*conv2;
      }
	
      if (species(atom2).contains("O")) {
         C6two = 4.15*conv;
         Ptwo = 0.569*conv2;
      }

      if (species(atom2).contains("N")) {
         C6two = 11.55*conv;
         Ptwo = 1.030*conv2;
      }
   
      //--- Slater-Kirkwood interpolation formula 
      //    Note: Elstner paper (Eq. 7) is wrong here
      //    Formula is correct in T. Halgren, JACS 114, 7827 (1992) 
      C6 = (2.*C6one*C6two*Pone*Ptwo)
         / (Ptwo*Ptwo*C6one + Pone*Pone*C6two);
   }
 
   return C6;
}	

	


double SxVDW::getParam (SxString name, int atom1, int atom2) 
{
	
	if (name == SxString ("covalentCutoff")) {
		return 3.0;
	}

	if (name == SxString ("vanDerWaalsCutoff")) {
		return 3. * getRm(atom1, atom2);
	}
	
	if (name == SxString ("cdamp")) return 3.54;
	
	if (name == SxString ("dstar")) return 3.00;
	
   if (name == SxString ("C6")) {
		return getC6(atom1, atom2);
   }
	
	if (name == SxString ("beta")) 
		return 23.0;
	
	if (name == SxString ("epsilon")) return 1.;
	
	if (name == SxString ("lamda")) return 21.0;
	if (name == SxString ("sigma")) return 1.;
	if (name == SxString ("gamma")) return 1.2;
	if (name == SxString ("constant")) return 21.0;

	
	return 1.;
}

void SxVDW::printParameters () {
	cout << "\nEmpirical Potential" << endl;
	cout << "\nNrOfAtoms " << nAtoms << endl;
	cout << "Omega" << superCell.determinant() << endl;
}

void SxVDW::printCoord () {
	
	cout << "\nCoordinates: \n";
	for (int i = 0; i < nAtoms; i++) {
		cout << coord(i) << endl;
	}
}

void SxVDW::printSpecies () {
	
	cout << "\nCoordinates: \n";
	for (int i = 0; i < nAtoms; i++) {
		cout << species(i) << endl;
	}
}

void SxVDW::printNeighbours () {
	
	for (int i = 0; i < nAtoms; i++) {
		cout << "Atom " << i << " has "<<  neighbours(i).getSize()
		     << " neighbours: ";
		for (int j = 0; j < neighbours(i).getSize(); j++) {
			/*
				if (neighbours(i)(j) == 1)
				cout << j << " ";
		}
		*/
			cout << neighbours(i)(j) << " ";
	}
	cout << endl;
}
}




