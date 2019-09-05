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

#ifndef _SX_EAM_H_
#define _SX_EAM_H_

#include <SxClassic.h>
#include <SxPotential.h>

/** \brief embedded atom empirical potential

    \b SxEAM = S/PHI/nX Embedded Atom Empirical Potential

    ...

    \author Lars Ismer, ismer@mpie.de */
class SX_EXPORT_CLASSIC SxEAM : public SxPotential
{
   public:

      SxEAM ();
      SxEAM (const SxAtomicStructure &, const  SxSymbolTable *);
      virtual ~SxEAM ();

      virtual bool isRegistered (const SxSymbolTable *) const;
      virtual void execute (const SxSymbolTable *, bool calc=true);
      
      virtual SxAtomicStructure getForces (const SxAtomicStructure &,
                                           const SxSymbolTable * =NULL);

      virtual SxSpeciesData getSpeciesData () const;
      virtual PrecEnergy getEnergy () const;
      virtual PrecEnergy getPotentialEnergy () const;
      
      SxAtomicStructure getNumericalForces (const SxAtomicStructure &,
            const double &);
      SxArray<SxAtomicStructure> getNumericalHesseMatrix (const SxAtomicStructure &,
            const double &);
      SxArray<SxArray<SxAtomicStructure> > getNumericalTTensor (const SxAtomicStructure &,
            const double &);
      void exportForces (const SxAtomicStructure &, const SxString &);
      void exportHesse (const SxArray<SxAtomicStructure> &, const SxString &);
      void exportTTensor (const SxArray<SxArray<SxAtomicStructure> > &, const SxString &);
      
      double getTotalEnergy ();
      double getV (const double &, int);
      double getRho (const double &i, int);
      double getF (const double &i, int);
      double getdV (const double &, int);
      double getdRho (const double &i, int);
      double getdF (const double &i, int);


      void setupSplines (const SxSymbolTable *);
      void updateNeighs ();

      // functions for link cell method
      void setLinkCellMethodParams (const SxAtomicStructure &);
      void updateNeighsLinkCellMethod ();
      //
      void setFixedNeighborsMethodParams (const SxAtomicStructure &);
      void updateNeighsFixedNeighborsMethod ();
      
      void update (const SxAtomicStructure &);

      void spline 
         (const SxArray<double> &, const SxArray<double> &, int , 
          double, double,  SxArray<double> *);

      double splint 
         (const SxArray<double> &, const SxArray<double> &, 
          const SxArray<double> &, int n, double x);
 
      
      double getDist (int, int, const SxVector3<Double> &);
      
      
      int getInteractionType (int, int);
      
      SxList<double> readData (const SxString &);

      
      int nPoints, nAtoms, n, nSpecies;
      double cutoff;
      double totalEnergy;
      SxAtomicStructure tau;
      SxList<SxArray<SxArray<double> > > V;
      SxList<SxArray<SxArray<double> > > F;
      SxList<SxArray<SxArray<double> > > rho;
      
      SxList<SxArray<SxArray<double> > > dV;
      SxList<SxArray<SxArray<double> > > dF;
      SxList<SxArray<SxArray<double> > > drho;


      SxSpeciesData speciesData;
     
      SxArray<SxArray<int> > neigh;
      SxArray<SxArray<int> > pw;
      SxArray<SxArray<int> > weight;
      SxArray<SxArray<double> > dist;
      SxArray<double > rhoSum;
      SxArray<SxArray<SxVector3<Double> > > cell;

      // constants for link cell method
      double length, puffer;
      bool noLinkCell;
      SxVector3<Double> meshOrigin;
      SxArray<SxVector3<Int> > origMesh;
      SxArray<SxArray<SxVector3<Int> > > neighbors;

      // constants for fixedNeighborsMethod
      double cutoffPlusRadius, radius;
      SxArray<SxArray<int> > potentialNeighborsNr;
      SxArray<SxArray<SxVector3<Double> > > potentialNeighborsCell;
      SxAtomicStructure startTau;
      SxArray<SxString> species;


};

#endif /* _SX_EAM_H_ */
