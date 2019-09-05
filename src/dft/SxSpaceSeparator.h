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

#ifndef _SX_SPACESX_SEPARATOR_H_
#define _SX_SPACESX_SEPARATOR_H_

#include <SxBinIO.h>
#include <SxMesh3D.h>
#include <SxNeighbors.h>
#include <SxRBasis.h>
#include <SxProjector.h>
#include <SxRho.h>

class SxSpaceSeparator
{
   public:

      /// Constructor
      SxSpaceSeparator () { /*empty*/ };

      SxSpaceSeparator (const SxAtomicStructure &structureIn);

      ~SxSpaceSeparator () { /*empty*/ };

      SxArray<SxDiracVec<Double> > voronoi (const SxRBasis &R);

      SxArray<SxDiracVec<Double> > bader (const SxDiracVec<Double> &rhoIn);

      SxArray<SxDiracVec<Double> > baderTrinkle (
            const SxDiracVec<Double> &rhoIn);

      ssize_t getNextPerGrad (ssize_t index, SxVector3<Double> &corr);

      SxArray<int> getNextAtom (ssize_t index);

      double getShareFactor(ssize_t index, 
            const SxDiracVec<Int> &atomicSpace, int iAtom);

      int findMaximum (ssize_t index, SxDiracVec<Int> &currentMap, 
                       SxVector3<Double> &corr, bool setTrajectory);
      int getCriticalType (ssize_t index);

      bool isBorderPoint (ssize_t index, SxDiracVec<Int> &currentMap);

      SxDiracVec<Double> getFilter (
            const SxDiracVec<Int> atomicSpace, int iAtom);

      SxList<SxDiracVec<Double> > computeGradRho ();
      SxList<SxDiracVec<Double> > computeGradGradRho ();
      SxDiracVec<Double> getGradRho(int dir);
      SxDiracVec<Double> getGradGradRho(int dir);

      ssize_t getNMaxima() {return maximaIndices.getSize();};
      SxArray<ssize_t> getMaximaIndices() {
         if (getNMaxima () > 0) return SxArray<ssize_t>(maximaIndices);
         else return SxArray<ssize_t> (0);
      };
      ssize_t getNMinima() {return minimaIndices.getSize ();};
      SxArray<ssize_t> getMinimaIndices() {
         if (getNMinima () > 0) return SxArray<ssize_t>(minimaIndices);
         else return SxArray<ssize_t> (0);

      };
      ssize_t getNRing() {return ringIndices.getSize ();};
      SxArray<ssize_t> getRingIndices() {
         if (getNRing () > 0) return SxArray<ssize_t>(ringIndices);
         else return SxArray<ssize_t> (0);
      };
      ssize_t getNBond() {return bondIndices.getSize ();};
      SxArray<ssize_t> getBondIndices() {
         if (getNBond () > 0) return SxArray<ssize_t>(bondIndices);
         else return SxArray<ssize_t> (0);
      };

      void plotDensity(const SxString &filename);
      void plotGradDensity(const SxString &filename);

      SxVector<Double> PAWDist;

   protected:

      enum pointType {Maximum, Interior, Boundary};

      SxDiracVec<Double> rho;

      SxList<SxDiracVec<Double> > gradRhoComponents;

      SxList<SxDiracVec<Double> > gradGradRhoComponents;

      SxAtomicStructure structure;

      SxMesh3D mesh;

      SxStack<ssize_t> maximaIndices;
      SxStack<ssize_t> ringIndices;
      SxStack<ssize_t> bondIndices;
      SxStack<ssize_t> minimaIndices;
      
      SxArray<ssize_t> findNeighbors(ssize_t index, int shell);

      double ramp (double u);

      SxDiracVec<Double> getProbabilityFlux (ssize_t index,
            const SxDiracVec<Double> &rhoIn);

      void printWeights (const SxArray<SxDiracVec<Double> > &weights, ssize_t idx);

};

namespace Timer {
   enum SpaceSeparatorTimer {
      Voronoi,
      Bader,
      BaderTrinkle
   };
}

SX_REGISTER_TIMERS(Timer::SpaceSeparatorTimer)
{
   using namespace Timer;
   regTimer (Voronoi,     "Voronoi");
   regTimer (Bader,       "Bader bond critical");
   regTimer (BaderTrinkle,"Bader-Trinkle");
}

#endif /* _SX_SPACESX_SEPARATOR_H_ */
