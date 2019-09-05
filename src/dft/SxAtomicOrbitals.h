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

#ifndef _SX_ATOMIC_ORBITALS_H_
#define _SX_ATOMIC_ORBITALS_H_

#include <SxPrecision.h>
#include <SxRadBasis.h>
#include <SxRadRBasis.h>
#include <SxRadGBasis.h>
#include <SxGkBasis.h>
#include <SxPsiSet.h>
#include <SxDirac.h>
#include <SxOrbital.h>
#include <SxPW.h>
#include <SxArray.h>
#include <SxDFT.h>
#include <SxNaturalCubicSpline.h>
#include <SxPAWPot.h>
#include <SxPseudoPot.h>
#include <SxHDF5.h>

/** \brief Container of atomic orbitals \f$
              \langle r | \mu_{i_s,i_a,n,l,m} \rangle
           \f$

    \b SxAtomicOrbitals = S/PHI/nX Atomic Orbitals

    This class is a container for atomic orbitals sampled on a radial
    grid. The atomic orbitals can be specified by the index of the
    species \f$i_s\f$ and the angular momentum \em l.

    \sa      \ref page_dirac
    \ingroup group_dft
    \ingroup group_dirac
    \author  Sixten Boeck
  */
class SX_EXPORT_DFT SxAtomicOrbitals : public SxPsiSet
{
   public:

      SxAtomicOrbitals ();
      /** \brief Constructor

          \param radFunc   The function values sampled on the radial
                           grid in :iSpecies,:l storage order
          \param radBasis  The reference to the radial basis
       */
      SxAtomicOrbitals (const SxArray<SxArray<SxRadBasis::TPsi> >  &,
                        const SxConstPtr<SxRadBasis> radBasisPtrIn);
      
      SxAtomicOrbitals (const SxArray<SxDiracMat<Double> > &in,
                        const SxConstPtr<SxRadBasis> radBasisPtrIn);


      /// \brief Copy Constructor
      SxAtomicOrbitals (const SxAtomicOrbitals &);

      /// \brief Read In Constructor
      SxAtomicOrbitals (SxBinIO &io);

      /// \brief Destructor
      virtual ~SxAtomicOrbitals ();

      /// \brief Assignment operator
      void operator= (const SxAtomicOrbitals &);

      void operator+= (const SxAtomicOrbitals &);

      SxAtomicOrbitals operator+ (const SxAtomicOrbitals &) const;
      SxAtomicOrbitals operator- (const SxAtomicOrbitals &) const;

      SxAtomicOrbitals operator* (double skalar) const;
      SxAtomicOrbitals operator* (const SxArray<SxArray<double> > &skalar) const;
      SxAtomicOrbitals operator* (const SxAtomicOrbitals &in) const;

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
      int getNOrbTypes () const;

      //// \brief get muSet(is)(iot)
      SxRadBasis::TPsi &operator() (int is, int iot);
      const SxRadBasis::TPsi &operator() (int is, int iot) const;

      int getIOT (int is, int n, int l) const;
      /** \brief Extract a single orbital

          In principle only \em iSpecies and \em l are needed to extract
          a single orbital. However, for some projections (e.g. onto 
          plane-waves) a Dirac vector also need to know it's atomic index and
          the quantum numbers. This function extracts an orbital and 
          initializes the necessary auxillary Dirac data. */
      SxRadBasis::TPsi operator() (int is, int ia, int n, int l, int m) const;

      void set (double val);

      SxQuantumNumbers getQuantumNumbers (int iSpecies, int idx) const;

      SxArray<SxQuantumNumbers> getReducedOrbitalMap () const;
      
      SxArray<SxQuantumNumbers> getOrbitalMap (SxAtomicStructure structure) const;

      int getNOrbitals (SxAtomicStructure structure) const;

      void createFuncLMap ();

      SxRadBasis::TPsi getFuncL (int is, int l, int ifl);

      const SxRadBasis::TPsi getFuncL (int is, int l, int ifl) const;

      int getLMax () const;

      int getLMax (const int iSpecies) const;

      SxArray<SxArray<int> > getFuncPerL () const;

      int getFuncPerL (int is, int l) const;

      void write (SxBinIO &io) const;

      void write (SxString filename) const;

      void read (SxBinIO &io);

      void read (const SxString &file);

      void readSiesta (const SxString &file, 
                       const SxDiracVec<Double> &radFunc);
      
      void writeSiesta (const SxArray<SxString> &file, 
                        const SxArray<SxDiracVec<Double> > &basis);

      /// Set the radial basis pointer
      void setBasis (const SxConstPtr<SxRadBasis> );

      const SxRadBasis& getRadBasis () const { return *radBasisPtr; };
      SxConstPtr<SxRadBasis> getBasis () const { return radBasisPtr; };
      SxConstPtr<SxRadBasis> getRadBasisPtr () const { return radBasisPtr; };

      void print (SxString file) const;
      void print () const;

      void setup (const SxSymbolTable *table);

      void addOrbital (const SxDiracVec<Double> &orbitalIn);

      /** \brief Read an orbital from a file
        @param fp an open file, containing lines of "r psi(r)"
        @param is  species
        @param n   main index
        @param l   l channel
        @param iot  if set, read to orbital with index iot.
                    Append orbital otherwise.
        @param ignoreComments if # line comments should be ignored 
        */
      void readOrbital (FILE *fp, int is, int n, int l, int iot=-1, 
                        bool ignoreComments = true);
      /** \brief Read orbitals from a file
        @param fp   an open file, containing lines of "r psi(r)" preceded by
                    comment lines that give l= (and optionally n=)
        @param is   species
       */
      void readOrbitals (FILE *fp, int is);

      double getNormSqr(int is, int iot) const;
      double getNormSqrSum() const;
      double dot (const SxAtomicOrbitals &in) const;
      double dot (const SxAtomicOrbitals &in, int is, int iot, int js, int jot) const;

      void normalize ();

      SxArray<SxArray<SxMatrix<Double> > > orthogonalize ();

      SxArray<SxArray<SxMatrix<Double> > > getOverlap () const;

      SxMatrix<Double> getOverlap (int is, int l) const;

      SxArray<SxArray<SxArray<int> > > funcLMap;

      int getOrbitalIdx (int is, int ia, int iot, int l, int m, SxArray<SxQuantumNumbers> map) const;

      void refine (SxArray<int> &factor);

      SxDiracVec<Double> compressWave(SxString &file, 
            int iState, 
            int iSpecies, 
            int l);
      
      #ifdef USE_HDF5
      void writeHDF5(const SxString &name);
      #endif

   protected:

      /** \brief The pointer to the \f$|r\rangle\f$ basis. */
      SxConstPtr<SxRadBasis>               radBasisPtr;

      /** \brief The sampling points of each orbitals.

        <b>Strorage order</b>: 
        -# iSpecies
        -# l
        -# r */
      SxArray<SxArray<SxRadBasis::TPsi> >  muSet;  // :is,:iot,:r

      void registerMemoryObservers ();
   public:
      /// Get the muSet
      const SxArray<SxArray<SxRadBasis::TPsi> > & getMuSet () const
      {
         return muSet;
      }

      /// Inquire l-number of some orbital type
      int getL (int is, int iot) const
      {
         SX_CHECK (is >= 0 && is < muSet.getSize (), is, muSet.getSize ());
         SX_CHECK (iot >= 0 && iot < muSet(is).getSize (),
                   iot, muSet(is).getSize ());
         return muSet(is)(iot).handle->auxData.l;
      }

};

inline SX_EXPORT_DFT SxAtomicOrbitals operator* (double skalar, const SxAtomicOrbitals &in)
{
   return in * skalar;
}

inline SX_EXPORT_DFT SxAtomicOrbitals operator* (const SxArray<SxArray<double> > &skalar, const SxAtomicOrbitals &in)
{
   return in * skalar;
}

#endif /* _SX_ATOMIC_ORBITALS_H_ */
