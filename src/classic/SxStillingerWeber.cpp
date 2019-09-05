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

#include <SxStillingerWeber.h>
SxStillingerWeber::SxStillingerWeber (const SxAtomicStructure &/*str*/, 
                          const SxSymbolTable *table)
{
   sxprintf ("This is SxStillingerWeber::SxStillingerWeber\n");
   speciesData = SxSpeciesData (&*table);
}

SxStillingerWeber::SxStillingerWeber ()
{  //SVN_HEAD; // this line MUST NOT be removed!!!
   sxprintf ("This is SxStillingerWeber::SxStillingerWeber\n");
}

SxStillingerWeber::~SxStillingerWeber ()
{
   // empty
}


bool SxStillingerWeber::isRegistered (const SxSymbolTable *) const
{
   return false;
}

void SxStillingerWeber::execute (const SxSymbolTable *, bool /*calc*/)
{
   cout << "This is StilWeb. " << endl;
   return;
}

void SxStillingerWeber::resize (int nA, double scale, SxMatrix3<Double> aMat)
{
	
	nAtoms                = nA;
	scalingFactor         = scale;
	superCell             = scalingFactor*aMat;
	
	coord                 = SxArray<SxVector3<Double> > (nAtoms);
	directEnergyContrib   = SxArray<double>             (nAtoms);
	indirectEnergyContrib = SxArray<double>             (nAtoms);
	
	dist                  = SxArray<SxList<double> > (nAtoms);
	expTerm               = SxArray<SxList<double> > (nAtoms);
	neighbours            = SxArray<SxList<int> >    (nAtoms);
	supercellIds          = SxArray<SxList<int> >    (nAtoms);
	sNeighbours            = SxArray<SxList<int> >    (nAtoms);
	sSupercellIds          = SxArray<SxList<int> >    (nAtoms);

	borderCrossers            = SxArray<SxList<int> >    (27);
	output = false;

	Forces = SxArray<SxVector3<Double> > (nAtoms);
	directForces = SxArray<SxVector3<Double> > (nAtoms); 
	indirectForces = SxArray<SxVector3<Double> > (nAtoms); 

}

void SxStillingerWeber::update (SxArray<SxVector3<Double > > newcoord) 
{
	if (output) cout << " ---- Updating Coords" << endl;
	updateCoord (newcoord);
	//printParameters();
	//printCoord();
	//updateDistances();
	if (output) cout << " ---- Updating Border Crossers" << endl;
	updateBorderCrossers();
	if (output) cout << " ---- Updating Neighbourhood and Distances" << endl;
	updateNeighbours ();
	//printNeighbours();
	//updateSecondNeighbours();
	//printSecondNeighbours();
}

SxVector3<Double> SxStillingerWeber::getLatticeVec(int supercellId) 
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

double SxStillingerWeber::getDist (int i, int j, int supercellId) {
	SxVector3<Double> distVec;
	
	distVec = coord(i) - coord(j) - getLatticeVec(supercellId);
	
	return 
		sqrt (distVec(0)*distVec(0) + distVec(1)*distVec(1) + distVec(2)*distVec(2));
}
		
void SxStillingerWeber::updateCoord (SxArray<SxVector3<Double> > newcoord) 
{
	
	for ( int i = 0; 
		   i < nAtoms; 
			i++) {
		coord(i) = newcoord(i);
	}
}

	    	  
void SxStillingerWeber::updateBorderCrossers () 
{
	double maxDist = getParam("a", 0, 0);
	double length, projection;
	int j, i, superCellId, a0Flip, a1Flip, a2Flip;
	SxArray<int> acoeff(3);

	// comment 20.02.2003: not so sure about the following comment
	// only works for orthorombic supercells, but easy to generalize


	bool smallSuperCell = false;
	//smallSuperCell = true;
	
	//-- should be placed elsewhere (Initialisiation)
	double abssqr;
	for (j = 0; j < 3; j++) {
		abssqr = superCell(j)(0)*superCell(j)(0) 
			    + superCell(j)(1)*superCell(j)(1)
				 + superCell(j)(2)*superCell(j)(2);
		
		if (sqrt(abssqr) < 2.*maxDist) smallSuperCell = true;
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

		acoeff(0) = acoeff(1) = acoeff(2) = 0;
		for (j = 0; j < 3; j++) {
			length = sqrt(superCell(j).absSqr ().sum ());
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
		

	}
	
	}
}
		
	
	
			
					

					
							
					  	

void SxStillingerWeber::updateNeighbours () 
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
				if ( ((i != j)) && ((distance = getDist(i, j, supercells)) < getParam ("a",i , j))) {
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
				if ((distance = getDist(i, j, supercells)) < getParam ("a",i , j)) { 
					neighbours(i).append (j);
					dist(i).append (distance);
					supercellIds(i).append (supercells);
				}

			}
		}
	}
	}
}
 
void  SxStillingerWeber::updateSecondNeighbours() {
	sNeighbours   = SxArray<SxList<int> >    (nAtoms);
	sSupercellIds = SxArray<SxList<int> >    (nAtoms);
	int i, j, k, neighj, neighk, supercellId, position;
	
	for (i = 0; i < nAtoms; i++) {
		sNeighbours(i).resize (0);
		sSupercellIds(i).resize (0);

		for (j = 0; j < neighbours(i).getSize (); j++) {
			sNeighbours(i).append (neighbours(i)(j));
			sSupercellIds(i).append (supercellIds(i)(j));
		}
		for (j = 0; j < neighbours(i).getSize (); j++) {
			neighj = neighbours(i)(j);
		
			for (k = 0; k < neighbours(neighj).getSize (); k++) {
				
				neighk = neighbours(neighj)(k);
				supercellId = supercellIds(neighj)(k);
				
				if (!(sNeighbours(i).contains(neighk))) {
					if (!(neighk == i) && (supercellId == 0)) {
					sNeighbours(i).append (neighk);
					sSupercellIds(i).append (supercellId);
					}
				} else {
					position = (int)sNeighbours(i).findPos(neighk);
					if (sSupercellIds(i)(position) != supercellId) {
					  cout << "Warning: Atom with position " 
					       << coord(neighk) << " is neighbour "
							 << "of Atom with position " << coord(i)
							 << " in at least to supercells: \n" 
							 << getLatticeVec(supercellId) << endl
							 << getLatticeVec(sSupercellIds(i)(position)) << endl;
						sNeighbours(i).append (neighk);
						sSupercellIds(i).append (supercellId);
					}
				}
			}
		}
	}
}
					
							 
					  		 
					
	 
SxVector3<Double>  SxStillingerWeber::getIndirectCenterForceOnAtom (int i) {
	
	SxVector3<Double>  returnValue, deltaIJ, deltaIK, radialContribForceJ,
	                   radialContribForceK, angularContribForceJ, angularContribForceK;
 
	
	double /* distBySigma, */indirectPotentialContrib,
	       aij, gammaij, contribJ, contribK, gammaik, 
	       aik, sigmaik, sigmaij, /* epsilonik, */derivativeJ, 
			 theta, constantik,
			 derivativeK, distij, distik, expij, expik;
	
	int neighj, neighk, superCellIdij, superCellIdik;		 
	
	SxList<int>::Iterator itNeighJ;
	SxList<int>::Iterator itNeighK;
	SxList<double>::Iterator itDistJ;
	SxList<double>::Iterator itDistK;
	SxList<double>::Iterator itExpJ;
	SxList<double>::Iterator itExpK;
	SxList<int>::Iterator itSCIDK;
	SxList<int>::Iterator itSCIDJ;


	
	returnValue(0) = 0.;
	returnValue(1) = 0.;
	returnValue(2) = 0.;
   
	itNeighJ = neighbours(i).begin();
	itExpJ = expTerm(i).begin();
	itDistJ  = dist(i).begin();
	itSCIDJ = supercellIds(i).begin();


	for (int j = 0; j < neighbours(i).getSize(); j++) {
      
		neighj = *itNeighJ; 
		aij = getParam ("a", i, neighj);
		gammaij = getParam ("gamma", i, neighj);
		sigmaij = getParam ("sigma", i, neighj);
	
		superCellIdij = *itSCIDJ;

			distij = *itDistJ;
			expij = *itExpJ;
				 
			
			contribJ = (-1.)*(gammaij/sigmaij)
				* pow (( distij/sigmaij - aij ), -2.); 
 
			itNeighK = neighbours(i).begin();
			itDistK  = dist(i).begin();
			itExpK   = expTerm(i).begin();
			itSCIDK  = supercellIds(i).begin();

			//for (int k = 0; k < neighbours(i).getSize (); k++) {
				for (int k = 0; k < j; k++) {

			//for (int k = 0; k < neighbours(i).getSize(); k++) {
				//if ((neighbours (i)(k)) && ( k!=j ) ){	
   
				neighk = *itNeighK;
				superCellIdik = *itSCIDK;
				if ((neighj != neighk) || (superCellIdij != superCellIdik)) {
					aik   = getParam ("a", i, neighk);
					gammaik   = getParam ("gamma", i, neighk);
					sigmaik   = getParam ("sigma", i, neighk);
					//epsilonik = getParam ("epsilon", i, neighk);
					constantik = getParam ("constant", i, neighk);
					distik = *itDistK;
					expik = *itExpK;

					//distBySigma = distik/sigmaik;

					deltaIJ = coord(i) - coord(neighj) 
						                - getLatticeVec (superCellIdij);
					deltaIK = coord(i) - coord(neighk) 
						                - getLatticeVec (superCellIdik);

					theta 
						= (deltaIJ(0)*deltaIK(0) + deltaIJ(1)*deltaIK(1) 
								+ deltaIJ(2)*deltaIK(2))
						      /distij/distik;

					indirectPotentialContrib 
						= constantik 
						* (pow (expij, gammaij))
						* (pow (expik, gammaik))
						* pow((theta + 1./3.), 2.);

					//cout << "Indir. Pot. Contrib. "
					//	  << "(" << i << ", " << neighj <<", " << neighk << "):" 
					//	  << indirectPotentialContrib << endl;

					/*
					if ((i < neighj) && (neighj < neighk)) {
					indirectEnergyContrib(i) += 
						indirectPotentialContrib;
					}
					*/
					//if ((i < neighj) && (neighj < neighk)) 
					
					totalEnergy += indirectPotentialContrib;

					derivativeJ = indirectPotentialContrib * contribJ;
				       			
					radialContribForceJ 
						= deltaIJ/distij*derivativeJ;

					contribK = (-1.)*(gammaik/sigmaik)
						* pow ((distik/sigmaik - aik), -2.); 
 

					derivativeK = indirectPotentialContrib * contribK;

					radialContribForceK 
						= deltaIK/distik*derivativeK;
					/*
					cout << "Rad. Forc."
						  << "(" << i << ", " << neighj <<", " << neighk << "):" 
						  << "("<< radialContribForceK << ", " 
						  << radialContribForceJ << ")" << endl;
				
					cout << "Dist/Delta"
						  << "(" << i << ", " << neighj <<", " << neighk << "):" 
						  << "("<< dist(i)(k) << ", " 
						  << deltaIK << ")" << endl;
				
					cout << "derivativeK"
						  << "(" << i << ", " << neighj <<", " << neighk << "):" 
						  << "("  
						  << derivativeK << ")" << endl;

*/
				if (theta != -(1./3.)) {	
					angularContribForceJ =
						(deltaIK/distik 
						 - theta*deltaIJ/distij)/distij
						* 2. * indirectPotentialContrib
						/ (theta + 1./3.);
	
					angularContribForceK =
						(deltaIJ/distij 
						 - theta*deltaIK/distik)/distik
						* 2. * indirectPotentialContrib
						/ (theta + 1./3.);
				}
				else {
					angularContribForceJ(0)=0.;angularContribForceJ(1)=0.;
					angularContribForceJ(2)=0.;
					angularContribForceK =	angularContribForceJ;
				}	

					/*cout << "BAng. Forc. "
						  << "(" << i << ", " << neighj <<", " << neighk << "):" 
						  << "("<< angularContribForceK << ", " 
						  << angularContribForceJ << ")" << endl;
					
					cout << "Dist/Delta"
						  << "(" << i << ", " << neighj <<", " << neighk << "):" 
						  << "("<< dist(i)(k) << ", " << dist(i)(j) <<  
						   deltaIK << "," << deltaIJ<< ")" << endl;
					
					cout << "Indir. Pot. Contrib, theta"
						  << indirectPotentialContrib << "," << theta << ")" << endl;

				
					cout << "derivativeK"
						  << "(" << i << ", " << neighj <<", " << neighk << "):" 
						  << "("  
						  << derivativeK << ")" << endl;


					cout << "TAng. Forc. "
						  << "(" << i << ", " << neighj <<", " << neighk << "):" 
						  << "("<< angularContribForceK << ", " 
						  << angularContribForceJ << ")" << endl;

					cout << "TRad. Forc. "
						  << "(" << i << ", " << neighj <<", " << neighk << "):" 
						  << "("<< radialContribForceK << ", " 
						  << radialContribForceJ << ")" << endl;
*/
				indirectForces(neighj) += angularContribForceJ + 
						  radialContribForceJ;

				indirectForces(neighk) += angularContribForceK + 
						  radialContribForceK;

				returnValue -= 
						(angularContribForceJ + angularContribForceK 
						 + radialContribForceJ + radialContribForceK);
					
			
				//}
			}
				++itExpK;
				++itDistK;
				++itSCIDK;	
				++itNeighK;	
			}
			++itExpJ;
			++itDistJ;
	++itSCIDJ;		
	++itNeighJ;

	}
	return returnValue;
}

void SxStillingerWeber::updateExponentialTerms () {
	double a, sigma;
   expTerm.resize (nAtoms);
	for (int i = 0; i < nAtoms; i++) {
		expTerm(i).resize (0);
		for (int j = 0; j < neighbours(i).getSize (); j++) {
			//if (neighbours (i)(j)) {
				a     = getParam ("a", i, j);
				sigma = getParam ("sigma", i, j);

				//cout << "Exp. Term (" << i << "," << j <<") :";
				
				expTerm(i). append  
				(exp (1. / (dist(i)(j)/sigma - a)));				
				
				//cout << expTerm (i)(j) << endl;
				
			}
			/*
				else {
				
			expTerm (i)(j) = 0.; 
		
			}
			*/
		//}	
	}	
}
	
SxVector3<Double>  SxStillingerWeber::getDirectForceOnAtom (int i) {
	
	SxVector3<Double>  returnValue;
	
	double distBySigma, directPotentialContrib,
	       a, A, B, sigma, epsilon, derivative;
			 
	int neighj;
	
	returnValue(0) = 0.;
	returnValue(1) = 0.;
	returnValue(2) = 0.;

	
	//cout << "Atom "<< i << endl;

	for (int j = 0; j < neighbours(i).getSize (); j++) {
		//if (neighbours (i)(j)) {
	      neighj = neighbours(i)(j);	
			//cout << "Bonding Partner: " << neighj << endl;;
			a = getParam ("a", i, neighj);
			A = getParam ("A", i, neighj);
			B = getParam ("B", i, neighj);
			sigma = getParam ("sigma", i, neighj);
			epsilon = 	getParam ("epsilon", i, neighj);

		
			distBySigma = dist(i)(j)/sigma;
			//cout << "(eps, a, A, B, sigma): " << epsilon << " "
			//	 << a <<" "<< A <<" "<< B <<" " << sigma << endl;


			directPotentialContrib  
				= epsilon*A
				* (B*(pow (distBySigma, -4.)) - 1.)
				* expTerm(i)(j);
	
			//cout << "Exponential Term : " << expTerm(i)(j) << endl;
			
			//cout << "Factor : " << (B*(pow (distBySigma, -4.)) - 1.)
 //<< endl;
		
			
	//		cout << "Direct Potential Contrib : " << directPotentialContrib << endl;
//			cout << "Total Energy : " << totalEnergy << endl;

			directEnergyContrib (i) +=  directPotentialContrib;
	   	totalEnergy += directPotentialContrib/2.;
			
			derivative  
				= epsilon*A*expTerm (i)(j)
				* ( (-1.) * pow(distBySigma - a, -2.)/sigma 
					  	* (B*pow (distBySigma, -4.) - 1.)
				  		- 4.*B*pow (distBySigma, -4.)/dist(i)(j));

			//cout << i << " " << neighj << " " << returnValue << endl;
			
			returnValue -= derivative * (coord(i) - coord(neighj) 
					       - getLatticeVec (supercellIds(i)(j)))
				            / dist(i)(j);
		   	
		//}	
	}
	
	return returnValue;
}

void SxStillingerWeber::updateDirectForces () {
	for (int i = 0; i < nAtoms; i++) {
		for (int j = 0; j < 3; j++) {
		directForces(i)(j) = 0;
		}
	}
		

	for (int i = 0; i < nAtoms; i++) {
		directForces(i) = getDirectForceOnAtom (i);
	}
}

void SxStillingerWeber::updateIndirectForces () {
	for (int i = 0; i < nAtoms; i++) {
		for (int j = 0; j < 3; j++) {
		indirectForces(i)(j) = 0;
		}
	}

	for (int i = 0; i < nAtoms; i++) {
		indirectForces(i) += getIndirectCenterForceOnAtom (i);
	}

}

SxArray<SxVector3<Double> > SxStillingerWeber::getForces () {
	
   if (output) cout << "...updating Exponential terms..." << endl;
	updateExponentialTerms ();
	
	totalEnergy = 0.;
	if (output) cout << "...2 Body Contributions..." << endl;
	updateDirectForces ();

	if (output) cout << "...3 Body Contributions..." << endl;
	updateIndirectForces ();
	
	for (int i = 0; i < nAtoms; i++) {
		//directEnergyContrib(i) = 0.;
		//indirectEnergyContrib(i) = 0.;
		
		//cout << "Direct Force on Atom " << i << ": " << directForces(i) << endl;
		//cout << "Indirect Force on Atom " << i << ": " << indirectForces(i) << endl;

		Forces(i) = (directForces(i) + indirectForces(i))/scalingFactor;
	}

	return Forces;
}
			
SxArray<SxVector3<Double> > SxStillingerWeber::getNumericalForces (double dx) {
	SxArray<SxVector3<Double> > forces     (nAtoms);
	SxArray<SxVector3<Double> > dummy      (nAtoms); 
	SxArray<SxVector3<Double> > undevCoord (nAtoms);
  	SxArray<SxVector3<Double> > devCoord   (nAtoms);
	double undevEnergy;

   update(coord);
	updateExponentialTerms ();
	
	dummy = getForces ();
	undevEnergy = totalEnergy;

	for (int i = 0; i < nAtoms; i++) {
      undevCoord(i) = coord(i);
		devCoord(i) = coord(i);
	}

	for (int i = 0; i < nAtoms; i++) {
		for (int j = 0; j < 3; j++) {
			forces(i)(j) = 0.;
			devCoord(i)(j) += dx;
			update(devCoord);
			updateExponentialTerms();
			dummy = getForces();
			forces(i)(j) -= (totalEnergy - undevEnergy)/dx/2./scalingFactor;
			devCoord(i)(j) -= 2.*dx;
			update(devCoord);
			updateExponentialTerms();
			dummy = getForces();
			forces(i)(j) += (totalEnergy - undevEnergy)/dx/2./scalingFactor;
			devCoord(i)(j) += dx;
		}
	}
  	totalEnergy = undevEnergy;
	return forces;
}

double SxStillingerWeber::getParam (SxString name, int /*i*/, int /*j*/) {
	
	if (name == SxString ("a")) {
		return 1.8;
	}
	if (name == SxString ("A")) return 7.049556277;
	if (name == SxString ("B")) return 0.6022245584;
	if (name == SxString ("epsilon")) return 1.;
	if (name == SxString ("lamda")) return 21.0;
	if (name == SxString ("sigma")) return 1.;
	if (name == SxString ("gamma")) return 1.2;
	if (name == SxString ("constant")) return 21.0;

	
	return 1.;
}

void SxStillingerWeber::printParameters () {
	cout << "\nEmpirical Potential" << endl;
	cout << "\nScaling-Factor " << scalingFactor << endl;
	cout << "\nNrOfAtoms " << nAtoms << endl;
	cout << "Omega" << superCell.determinant() << endl;
}

void SxStillingerWeber::printCoord () {
	
	cout << "\nCoordinates: \n";
	for (int i = 0; i < nAtoms; i++) {
		cout << coord(i) << endl;
	}
}

void SxStillingerWeber::printNeighbours () {
	
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
void SxStillingerWeber::printSecondNeighbours () {
	
	for (int i = 0; i < nAtoms; i++) {
		cout << "Atom " << i << " has "<<  sNeighbours(i).getSize()
		     << " second neighbours: ";
		for (int j = 0; j < sNeighbours(i).getSize(); j++) {
			/*
				if (neighbours(i)(j) == 1)
				cout << j << " ";
		}
		*/
			cout << sNeighbours(i)(j) << " ";
	}
	cout << endl;
}
}

double SxStillingerWeber::getEnergy () const
{
	return totalEnergy;
}


SxSpeciesData SxStillingerWeber::getSpeciesData () const
{
   return speciesData;
}

SxAtomicStructure SxStillingerWeber::getForces (const SxAtomicStructure &tau,
                                                const SxSymbolTable *)
{  //SVN_HEAD; // this line MUST NOT be removed

   int nAtoms_, i;
	double scale = 1.0;

	
	nAtoms_ = tau.getNAtoms ();
	SxArray<SxVector3<Double> > coords(nAtoms_);
	SxArray<SxVector3<Double> > forces(nAtoms_);
	SxArray<SxVector3<Double> > numForces(nAtoms_);

	
	for (i = 0; i < nAtoms_; i++) coords(i) = tau(i);

	
	resize (nAtoms_, scale, tau.cell );

	update (coords);
	forces = getForces ();

   cout << forces << endl;

//	numForces = getNumericalForces (0.001);

//	for (i = 0; i < nAtoms_; i++) {
//		cout << forces(i) << endl;
//		cout << numForces(i) << endl;
//		cout << endl;
//	}

   // --- place holder
   //SxAtomicStructure str (tau, true);
	
   SxAtomicStructure str;
	str.copy (tau);

	for (i = 0; i < nAtoms_; i++)  
		str.ref(i) = forces(i);

   return str;

}


