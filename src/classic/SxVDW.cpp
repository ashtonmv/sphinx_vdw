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
// Authors: Lars Ismer (ismer@fhi-berlin.mpg.de)
//			Michael Ashton (ashton@mpie.de)

#include <SxVDW.h>
#include <SxElemDB.h>

SxVDW::SxVDW () 
{ //SVN_HEAD;
   
}

SxVDW::SxVDW 
(const SxAtomicStructure &t, const SxSymbolTable *table, const SxRho &dens) 
{ 
    int i;
    SxAtomicStructure tau;
    SxParser parser;

    SxParser::Table elemTable = parser.read ("species/elements.sx");
    SxElemDB elemDB(elemTable);

    tau.copy(t);

	rho = dens;
    nAtoms = t.nTlAtoms;
	superCell             = tau.cell;
   //---workaround: will be solved if connected to sxatomicstructure
    if (superCell.absSqr().sum() < 1e-10) { 
        superCell (0, 0) = superCell(1, 1) = superCell(2, 2) = 100;
    }
      
	coord                 = SxArray<SxVector3<Double> > (nAtoms);
	species               = SxArray<SxString> (nAtoms);
	polarizability        = SxArray<double> (nAtoms);
	C6                    = SxArray<double> (nAtoms);
	vdwRadius             = SxArray<double> (nAtoms);
	effectiveVolume       = SxArray<double> (nAtoms);
	hirshfeldVolume       = SxArray<double> (nAtoms);
	freeVolume            = SxArray<double> (nAtoms);


    SxSpeciesData speciesData = SxSpeciesData(table->topLevel());
   
    for (i = 0; i < nAtoms; i++) {
        coord(i) = tau.ref (i);
        species(i) = speciesData.chemName(tau.getISpecies (i));
	    polarizability(i) = elemDB.getPolarizability(species(i));
 	    C6(i) = elemDB.getC6(species(i));
 	    vdwRadius(i) = elemDB.getVdwRadius(species(i));
		freeVolume(i) = 1.; // Need to figure out how to read this in.
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

    SxSymbolTable *vdwGroup = table -> getGroup ("vdwCorrection");

    correctionType = vdwGroup -> get("correctionType") -> toString ();
    if (vdwGroup -> contains("combinationRule")) {
        combinationRule = vdwGroup -> get("combinationRule") -> toString ();
    } else combinationRule = SxString("Default");
    if (vdwGroup -> contains("damping")) {
        damping = vdwGroup -> get("damping") -> toString ();
    } else damping = SxString("Fermi");

}

SxVDW::~SxVDW ()
{
   // empty
}

int SxVDW::getnAtoms() {
    return nAtoms;
}

void SxVDW::resize 
(SxList<SxList<SxVector3<Double> > >  tauarr, SxList<SxString> speciesNameList, SxMatrix3<Double> aMat)
{
	int i, j;
	int counter;
	nAtoms = 0;
	for (i = 0; i < tauarr.getSize (); i++) {
			nAtoms += (int)tauarr(i).getSize ();
	}

	superCell             = aMat;
	
	coord                 = SxArray<SxVector3<Double> > (nAtoms);
	species               = SxArray<SxString> (nAtoms);
	
	counter = 0;
	for (i = 0; i < tauarr.getSize (); i++) {
		for (j = 0; j < tauarr(i).getSize (); j++) {
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
	updateNeighbours (SxString ("vdwCutoff"));
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
	
	// This caused a crash.
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
 
double SxVDW::getDampingFunction (double R, double Rij) 
{
	double fd = 0.;
	double s6 = getParam("s6", 0, 0);
	double d = getParam("d", 0, 0);
    double sR = getParam("sR", 0, 0);
	double vdwCutoff = getParam("vdwCutoff", 0, 0);


	// Fermi damping with cutoff (set by default to~ 100 Bohr)

	if (R < vdwCutoff) {
		fd = s6 / (1 + exp(-d * (R / (sR * Rij) - 1)));
	}

	return fd;
}

double SxVDW::getDampingDerivative (double R, double Rij) 
{
	double fdprime = 0.;

	double s6 = getParam("cdamp", 0, 0);
	double d = getParam("beta", 0, 0);
    double sR = getParam("dstar", 0, 0);
	double vdwCutoff = getParam("vdwCutoff", 0, 0);
	double scale = 1.0 / (sR * Rij);
	double x = R * scale;
	double chi = exp(-d * (x - 1.0));

	// Again fermi damping with a cutoff, vdwCutoff

	if (R < vdwCutoff) {
		fdprime = d * scale * chi / ::pow(1.0 + chi, 2);
	}

	return fdprime;
}

/*
double SxVDW::getDampingSecondDerivative (double R, double Rm) {
	double fd, e, ePrime, ePrimePrime, a;
	double cdamp = getParam("cdamp", 0, 0);
	double beta = getParam("beta", 0, 0);
	double dstar = getParam("dstar", 0, 0);
   fd = e = ePrime = ePrimePrime = a = 0.;
	
	if (correctionType == SxString("WTYang_I")) {
		e = exp(-cdamp*(::pow( (R/Rm), 3.)));
		ePrime = -e*cdamp*3.*R*R/Rm/Rm/Rm;
		fd = 12*cdamp*R/Rm/Rm/Rm*(1.-e)*e
			+ 6.*cdamp*R*R/Rm/Rm/Rm*(ePrime - 2.*e*ePrime);
	}
	
	if (correctionType == SxString("WTYang_II")) {
		e = exp(-beta*(R/Rm - 1.));
		fd = - ::pow((beta/Rm), 2.)*e/(1. + e)/(1. + e)*(1. - 2*e/(1 + e));
	}

   if (correctionType.contains("Elstner")) {
		a = dstar/::pow (Rm, 7.);
      e = exp(-dstar*(::pow( (R/Rm), 7.)));
      ePrime = -7.*a*::pow(R, 6.)*e;
      ePrimePrime = 49.*a*a*::pow(R, 12.)*e - 42*a*::pow(R, 5.)*e; 
      fd = -(4.*ePrimePrime * ::pow((1 - e), 3.) 
            - 12*ePrime*ePrime*(::pow((1-e), 2)));
	}
   
	return fd;
}
*/

void computeHirshfeldVolume () {
	for (int i = 0; i<nAtoms; i++) {
		hirshfeldVolume(i) = 1.;
	}
}

void computeEffectiveVolume () {
	if (correctionType == SxString("TS")) {
		computeHirshfeldVolume ();
		for (int i = 0; i<nAtoms; i++) {
			effectiveVolume(i) = hirshfeldVolume(i) / freeVolume(i);
		}
	}
	else {
		for (int i = 0; i<nAtoms; i++) {
			effectiveVolume(i) = 1.;
		}
	}
}

void SxVDW::compute () {
	// update attributes "totalEnergy" (double) and
	// "Forces" (SxArray<SxVector3<Double>>)
	int i, j, neighj;
	double R, Rij, C6ij, fd, fdPrime, derivative;

	computeEffectiveVolume();

	// Reset vdW energy to 0
	totalEnergy = 0.;
	for (i = 0; i < nAtoms; i++) {
		// Reset force array for atom i to 0
		for (int j = 0; j < 3; j++) {
			Forces(i)(j) = 0;
		}
		for (j = 0; j < neighbours(i).getSize (); j++) {
			R = dist(i)(j);
			neighj = neighbours(i)(j);
			Rij = getRij (i, neighj);
			C6ij = getC6ij (i, neighj);
			fd = getDampingFunction (R, Rij);
			fdPrime = getDampingDerivative (R, Rij);

			// Update Energy
            totalEnergy += -fd*C6ij/(::pow(R, 6.))/2.;

			// Update Forces
			derivative = fdPrime * C6ij/::pow(R, 6.) -
						 6. * fd * C6ij/::pow(R, 7.);

			Forces(i) += derivative * (coord(i) - coord(neighj) 
					     - getLatticeVec (supercellIds(i)(j))) / R;

		}
	}
}

double SxVDW::getTotalEnergy () {
	double R, Rij, fd, eVDW, C6ij;
	int i, j, neighj;
	eVDW = 0.;
	for (i = 0; i < nAtoms; i++) {

		for (j = 0; j < neighbours(i).getSize (); j++) {
			R = dist(i)(j);
			neighj = neighbours(i)(j);
			Rij = getRij (i, neighj);
			C6ij = getC6ij (i, neighj);
			fd = getDampingFunction (R, Rij);
            eVDW += -fd * C6ij/(::pow(R, 6.))/2.;
		}
	}
	return eVDW;
					
}		

SxVector3<Double>  SxVDW::getForceOnAtom (int i) {
	
	SxVector3<Double>  returnValue;
	
    double R, Rij, C6ij, fd, fdPrime, derivative;
			 
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
			Rij = getRij(i, neighj);
			C6ij = getC6ij(i, neighj);
			fdPrime = getDampingDerivative (R, Rij);
			fd = getDampingFunction (R, Rij);

			derivative  
				= fdPrime*C6ij/::pow(R, 6.) - 6.*fd*C6ij/::pow(R, 7.); 

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
	//updateForces ();
	compute (); // inherently updates Forces
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

/*
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
	double Rij = getRij (atom1, atom2);
	double C6ij = getC6ij (atom1, atom2);

	
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


	fd = getDampingFunction (R, Rij);
	fdP = getDampingDerivative (R, Rij);
	fdPP = getDampingSecondDerivative (R, Rij);

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
*/
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

double SxVDW::getRij (int atom1, int atom2)
{
    /*
    Return the combined vdw radius of atom1 and atom2.
    If the TS correction is chosen, this value is modified
    according to each atom's Hirshfeld volume ratio.
    */
    double Ri = vdwRadius(atom1);
    double Rj = vdwRadius(atom2);

    if (correctionType == SxString("TS")) {
        Ri = Ri * getEffectiveVolume(atom1);
        Rj = Rj * getEffectiveVolume(atom2);
    }

    double Rij = Ri + Rj;
	return Rij;
}

double SxVDW::getC6ij (int atom1, int atom2) 
{
    /*
    Return the combined C6 dispersion coefficient between atom1 and atom2.
    If the TS correction is chosen, this value is modified
    according to each atom's Hirshfeld volume ratio.
    */
    double C6i, C6j, C6ij, alphai, alphaj;
    
    C6i = C6(atom1);
    C6j = C6(atom2);
    alphai = polarizability(atom1);
    alphaj = polarizability(atom2);

    if (correctionType == SxString("TS")) {
        double vEffi = effectiveVolume(atom1);
        double vEffj = effectiveVolume(atom2);

        C6i = C6i * ::pow(vEffi, 2);
        C6j = C6j * ::pow(vEffj, 2);
        alphai = alphai * vEffi;
        alphaj = alphaj * vEffj;
    }

    if (correctionType == SxString("TS") or combinationRule == SxString("Tang")) {
        // from Tang,K. T. Phys. Rev. 1969, 177, 108
        C6ij = (2 * C6i*C6j /
            ( C6j * alphai / alphaj + C6i * alphaj / alphai ) );
    }
    else if (combinationRule == SxString("GB")) {
        // from Gould & Bucko (https://arxiv.org/pdf/1604.02751.pdf) page 18
        double alphaij = ::pow ( alphai * alphaj, 0.5 );
        C6ij = 1.43 * ::pow (alphaij, 1.45);
    }
    else  C6ij = ::pow (C6(atom1) * C6(atom2), 0.5); // Default

    return C6ij;
}	


double SxVDW::getParam (SxString name, int atom1, int atom2) 
{
	
	if (name == SxString ("covalentCutoff")) {
		return 3.0;
	}

	if (name == SxString ("vdwCutoff")) {
		return 94.5;  // 50 Angstroms, but in Bohr
	}
	
	if (name == SxString ("cdamp")) return 3.54;
	
	if (name == SxString ("dstar")) return 3.00;
	
   if (name == SxString ("C6ij")) {
		return getC6ij(atom1, atom2);
   }
	
	if (name == SxString ("beta")) 
		return 23.0;
	
	if (name == SxString ("epsilon")) return 1.;
	
	if (name == SxString ("lamda")) return 21.0;
	if (name == SxString ("sigma")) return 1.;
	if (name == SxString ("gamma")) return 1.2;
	if (name == SxString ("constant")) return 21.0;

	if (name == SxString ("s6")) return 0.75;
	if (name == SxString ("d")) return 20.0;
	if (name == SxString ("sR")) return 1.0;

	
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




