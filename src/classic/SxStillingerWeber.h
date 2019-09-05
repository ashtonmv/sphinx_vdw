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

#ifndef _SX_STILLINGER_WEBER_H_
#define _SX_STILLINGER_WEBER_H_

#include <SxClassic.h>
#include <SxPotential.h>

/** \brief Stillinger-Weber empirical potential

    \b SxStillingerWeber = S/PHI/nX Stillinger-Weber Empirical Potential

    ...

    \author Lars Ismer, ismer@fhi-berlin.mpg.de */
class SX_EXPORT_CLASSIC SxStillingerWeber : public SxPotential
{
   public:

      SxStillingerWeber ();
      SxStillingerWeber (const SxAtomicStructure &, const  SxSymbolTable *);
      virtual ~SxStillingerWeber ();

      virtual bool isRegistered (const SxSymbolTable *) const;
      virtual void execute (const SxSymbolTable *, bool calc=true);
      
      virtual SxAtomicStructure getForces (const SxAtomicStructure &,
                                           const SxSymbolTable * =NULL);

      virtual SxSpeciesData getSpeciesData () const;
      virtual PrecEnergy getEnergy () const;
	
      SxArray<SxVector3<Double> > coord;

	  SxMatrix3<Double> superCell;
	  SxArray<double>  directEnergyContrib;
	  SxArray<double>  indirectEnergyContrib;

	  SxArray<SxList<int> > neighbours;
	   SxArray<SxList<int> > sNeighbours;

	  SxArray<SxList<int> >  borderType;
	  SxArray<SxList<int> > supercellIds;
	  SxArray<SxList<int> > sSupercellIds;


	  SxArray<SxList<int> > borderCrossers;
	  SxArray<SxVector3<Double> > directForces;
	  SxArray<SxVector3<Double> > indirectForces;
	  SxArray<SxVector3<Double> > Forces;

	  SxArray<SxList<double> > dist;
	  SxArray<SxList<double> > expTerm;
	  
	SxSpeciesData speciesData;


	  bool output;
	  int nAtoms;
	  double scalingFactor;
	  double totalEnergy;

	  
	  double getDist(int, int, int);
	  
	  void resize(int, double, SxMatrix3<Double>);
	  void update(SxArray<SxVector3<Double> >);
	  void updateDistances ();
	  void updateNeighbours ();
	  void updateSecondNeighbours ();

	  void updateBondType ();
	  void updateCoord (SxArray<SxVector3<Double> >);
	  void updateBorderCrossers ();
	  
	  void updateExponentialTerms (); 
	  
	  SxVector3<Double>  getLatticeVec(int);
	  SxVector3<Double>  getIndirectCenterForceOnAtom (int); 
	  SxVector3<Double>  getDirectForceOnAtom (int); 
	  void updateDirectForces (); 
	  void updateIndirectForces (); 
	  SxArray<SxVector3<Double> > getForces ();
	  SxArray<SxVector3<Double> > getNumericalForces (double);

	  double getParam (SxString, int, int); 

	  void printParameters ();
	  void printNeighbours ();
	  
	  void printSecondNeighbours ();
	  void printCoord ();


};

#endif /* _SX_STILLINGER_WEBER_H_ */
