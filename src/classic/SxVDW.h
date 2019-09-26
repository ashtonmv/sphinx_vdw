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

#ifndef _SX_VDW_H_
#define _SX_VDW_H_

#include <SxClassic.h>
#include <SxAtomicStructure.h>
#include <SxRho.h>
#include <SxSpeciesData.h>
#include <SxSymbolTable.h>

/** This class is for the empirical van der waals correction  */
class SX_EXPORT_CLASSIC SxVDW
{
   public:
	 
	  SxString correctionType;
      SxString combinationRule;
      SxString damping;
	  SxArray<SxVector3<Double> > coord;
	  SxArray<SxString> species;
	  
	  RhoR rho;

	  SxMatrix3<Double> superCell;
	  SxArray<double>  energyContrib;

	  SxArray<SxList<int> > neighbours;
      SxArray<SxList<int> > sNeighbours;

	  SxArray<SxList<int> >  borderType;
	  SxArray<SxList<int> > supercellIds;
	  SxArray<SxList<int> > sSupercellIds;
	  
      SxArray<SxList<int> > borderCrossers;
	  SxArray<SxVector3<Double> > Forces;

	  SxArray<SxList<double> > dist;
	  SxArray<SxList<double> > expTerm;

	  SxArray<double> polarizability;
	  SxArray<double> C6;
	  SxArray<double> vdwRadius;
	  SxArray<double> effectiveVolume;
	  SxArray<double> hirshfeldVolume;
	  SxArray<double> freeVolume;

	  bool output;
	  int nAtoms;
	  double totalEnergy;

	   SxVDW ();
      SxVDW (const SxAtomicStructure &, const SxSymbolTable *, const RhoR &);
     ~SxVDW ();

	  
	  double getDist(int, int, int);
	  void resize 
		  (SxList<SxList<SxVector3<Double> > > , SxList<SxString>, SxMatrix3<Double>);
	  void update(SxArray<SxVector3<Double> >);
	  void updateNeighbours (SxString);

	  void updateCoord (SxArray<SxVector3<Double> >);
	  void updateBorderCrossers ();
	  
	  void updateDampingTerms (); 
	  void updateHybridisation (); 
	  
      void computeHirshfeldVolume ();
      void computeEffectiveVolume ();
	  void compute ();

	  SxVector3<Double>  getLatticeVec(int);
	  SxVector3<Double>  getForceOnAtom (int); 
	  SxMatrix3<Double> getInteraction(int, int, int);
	  void updateForces (); 
	  SxArray<SxVector3<Double> > getForces ();
	  SxArray<SxVector3<Double> > getNumericalForces (double);
	  SxMatrix<Double> getNumericalHessian (double);
	  SxMatrix<Double> getHessian ();

	  bool areNeighbors(int,int);
	  int getNeighborIndex(int,int);
	  double getDampingFunction(double, double);
	  double getDampingDerivative(double, double);
	  double getDampingSecondDerivative(double, double);
	  
	  double getC6ij (int, int); 
	  double getRij (int, int); 
	  double getParam (SxString, int, int);
	  double getTotalEnergy ();

      int getnAtoms ();

	  void printParameters ();
	  void printNeighbours ();
	  
	  void printCoord ();
	  void printSpecies ();
};


#endif /* _SX_VDW_H_ */
