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

#ifndef _SX_ATOMIC_ORBITALSG_H_
#define _SX_ATOMIC_ORBITALSG_H_

#include <SxDirac.h>
#include <SxArray.h>
#include <SxBinIO.h>
#include <SxAtomicStructure.h>
#include <SxRadGBasis.h>

class SX_EXPORT_DFT SxAtomicOrbitalsG
{
   public:

      SxAtomicOrbitalsG ();
      SxAtomicOrbitalsG (const SxArray<SxArray<SxDiracVec<Double> > >  &in,
                         SxConstPtr<SxRadGBasis> radGBasisPtrIn,
                         bool splineRepIn);
      
      /// \brief Copy Constructor
      SxAtomicOrbitalsG (const SxAtomicOrbitalsG &);

      /// \brief Destructor
      virtual ~SxAtomicOrbitalsG ();

      void set (double val);

      /// \brief Assignment operator
      void operator= (const SxAtomicOrbitalsG &);

      void operator+= (const SxAtomicOrbitalsG &);
      void operator-= (const SxAtomicOrbitalsG &);

      SxAtomicOrbitalsG operator+ (const SxAtomicOrbitalsG &) const;
      void sumMPI() const;
      SxAtomicOrbitalsG operator- (const SxAtomicOrbitalsG &) const;

      SxAtomicOrbitalsG operator* (double skalar) const;
      SxAtomicOrbitalsG operator/ (double skalar) const;
      SxAtomicOrbitalsG operator* (const SxArray<SxArray<double> > &skalar) const;
      SxAtomicOrbitalsG operator* (const SxAtomicOrbitalsG &in) const;

      /// \brief Get number of Species
      int getNSpecies () const;

      int getNOrbTypes () const;
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

      double dot (const SxAtomicOrbitalsG &orbitalsIn, int is, int iot, int js, int jot) const;
      double dot (const SxAtomicOrbitalsG &orbitalsIn) const;

      double sum (const SxAtomicOrbitalsG &orbitalsIn, int is, int iot, int js, int jot) const;
      double sum (const SxAtomicOrbitalsG &orbitalsIn) const;

      void normalize ();

      void setBasis (SxConstPtr<SxRadGBasis> radGBasisPtrIn);

      SxConstPtr<SxRadGBasis> getRadGBasisPtr () const { return radGBasisPtr; };

      void toSpline () const;

      void toVec () const;

      SxDiracVec<Double> toVec (int is, int iot) const;

      bool isSpline () const {return splineRep;};

      void createFuncLMap ();

      SxRadGBasis::TPsi &getFuncL (int is, int l, int ifl);

      const SxRadGBasis::TPsi &getFuncL (int is, int l, int ifl) const;

      int getFuncPerL (int is, int l) const;

      SxArray<SxArray<SxMatrix<Double> > > orthogonalize ();

      SxArray<SxArray<SxMatrix<Double> > > getOverlap () const;

      SxMatrix<Double> getOverlap (int is, int l) const;

      SxArray<SxArray<SxArray<int> > > funcLMap;

      int getLMax () const;

      int getLMax (const int iSpecies) const;

      void orthogonalizeOn(SxAtomicOrbitalsG &basis);

      void rotate(SxArray<SxArray<SxMatrix<Double> > > rotMat);

      int getNOrbitals (SxAtomicStructure structure) const;

      int getIOT (int is, int n, int l) const;

      int getOrbitalIdx (int is, int ia, int iot, int l, int m, SxArray<SxQuantumNumbers> map) const;

   protected:

      mutable bool splineRep;

      SxConstPtr<SxRadGBasis>  radGBasisPtr;

      /** \brief The sampling points of each orbitals.

        <b>Strorage order</b>: 
        -# iSpecies
        -# l
        -# r */
      mutable SxArray<SxArray<SxDiracVec<Double> > >  muSet;  // :is,:iot,:r

};

inline SX_EXPORT_DFT SxAtomicOrbitalsG operator* (double skalar, const SxAtomicOrbitalsG &in)
{
   return in * skalar;
}

inline SX_EXPORT_DFT SxAtomicOrbitalsG operator* (const SxArray<SxArray<double> > &skalar, const SxAtomicOrbitalsG &in)
{
   return in * skalar;
}

#endif /* _SX_ATOMIC_ORBITALSG_H_ */
