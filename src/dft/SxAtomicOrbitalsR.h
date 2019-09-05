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

#ifndef _SX_ATOMIC_ORBITALSR_H_
#define _SX_ATOMIC_ORBITALSR_H_

#include <SxDirac.h>
#include <SxArray.h>
#include <SxBinIO.h>
#include <SxAtomicStructure.h>
#include <SxRadRBasis.h>

class SX_EXPORT_DFT SxAtomicOrbitalsR
{
   public:

      SxAtomicOrbitalsR ();
      SxAtomicOrbitalsR (const SxArray<SxArray<SxDiracVec<Double> > >  &in,
                         SxConstPtr<SxRadRBasis> radRBasisPtrIn,
                         bool splineRepIn);
      
      /// \brief Copy Constructor
      SxAtomicOrbitalsR (const SxAtomicOrbitalsR &);

      /// \brief Destructor
      virtual ~SxAtomicOrbitalsR ();

      void set (double val);

      /// \brief Assignment operator
      void operator= (const SxAtomicOrbitalsR &);

      void operator+= (const SxAtomicOrbitalsR &);
      void operator-= (const SxAtomicOrbitalsR &);

      SxAtomicOrbitalsR operator+ (const SxAtomicOrbitalsR &) const;
      SxAtomicOrbitalsR operator- (const SxAtomicOrbitalsR &) const;

      SxAtomicOrbitalsR operator* (double skalar) const;
      SxAtomicOrbitalsR operator* (const SxArray<SxArray<double> > &skalar) const;
      SxAtomicOrbitalsR operator* (const SxAtomicOrbitalsR &in) const;

      /// \brief Get number of Species
      int getNSpecies () const;

      /** \brief Gets the maximum angular component for a certain species.

          This function returns the index of the maximum angular component
          of a specified species.

          \b Note: The returned value is not the number of components. Thus,
          if used in for() loops the '<=' operator has to be used, e.g.
          \code
             // --- Loop over all angular components
             for (int l=0; l <= orbitals.getLMax(); l++)  {
                ...
             }
          \endcode
       */
      int getNOrbTypes (int iSpecies) const;

      SxArray<SxQuantumNumbers> getReducedOrbitalMap () const;
      
      SxArray<SxQuantumNumbers> getOrbitalMap (SxAtomicStructure structure) const;

      //// \brief get muSet(is)(iot)
      SxDiracVec<Double> &operator() (int is, int iot);
      const SxDiracVec<Double> &operator() (int is, int iot) const;

      SxDiracVec<Double> operator() (int is, int ia, int iot, int l, int m) const;
      
      void print (const SxString fileIn) const;

      double getNormSqr (int is, int iot) const;
      double getNormSqrSum () const;

      double dot (const SxAtomicOrbitalsR &orbitalsIn, int is, int iot, int js, int jot) const;
      double dot (const SxAtomicOrbitalsR &orbitalsIn) const;
      
      double sum (const SxAtomicOrbitalsR &orbitalsIn, int is, int iot, int js, int jot) const;
      double sum (const SxAtomicOrbitalsR &orbitalsIn) const;
      
      void normalize ();

      void setBasis (SxConstPtr<SxRadRBasis> radRBasisPtrIn);

      SxConstPtr<SxRadRBasis> getBasisPtr () const { return radRBasisPtr; };

      void toSpline () const;

      void toVec () const;

      SxDiracVec<Double> toVec (int is, int iot) const;

      bool isSpline () const {return splineRep;};

      void createFuncLMap ();

      SxRadRBasis::TPsi &getFuncL (int is, int l, int ifl);

      const SxRadRBasis::TPsi &getFuncL (int is, int l, int ifl) const;

      int getFuncPerL (int is, int l) const;

      SxArray<SxArray<SxMatrix<Double> > > orthogonalize ();

      SxArray<SxArray<SxMatrix<Double> > > getOverlap () const;

      SxMatrix<Double> getOverlap (int is, int l) const;

      SxArray<SxArray<SxArray<int> > > funcLMap;

      int getLMax () const;

      int getLMax (const int iSpecies) const;

      void orthogonalizeOn(SxAtomicOrbitalsR &basis);

      void rotate(SxArray<SxArray<SxMatrix<Double> > > rotMat);

      int getNOrbitals (SxAtomicStructure structure) const;

      int getIOT (int is, int n, int l) const;

      int getOrbitalIdx (int is, int ia, int iot, int l, int m, SxArray<SxQuantumNumbers> map) const;

   protected:

      mutable bool splineRep;

      SxConstPtr<SxRadRBasis>  radRBasisPtr;

      /** \brief The sampling points of each orbitals.

        <b>Strorage order</b>: 
        -# iSpecies
        -# l
        -# r */
      mutable SxArray<SxArray<SxDiracVec<Double> > >  muSet;  // :is,:iot,:r

};

inline SX_EXPORT_DFT SxAtomicOrbitalsR operator* (double skalar, const SxAtomicOrbitalsR &in)
{
   return in * skalar;
}

inline SX_EXPORT_DFT SxAtomicOrbitalsR operator* (const SxArray<SxArray<double> > &skalar, const SxAtomicOrbitalsR &in)
{
   return in * skalar;
}

#endif /* _SX_ATOMIC_ORBITALSR_H_ */
