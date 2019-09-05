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

#ifndef _SX_PSEUDO_POT_H_
#define _SX_PSEUDO_POT_H_

#include <SxSymbolTable.h>
#include <SxArray.h>
#include <SxDirac.h>
#include <SxDFT.h>
#include <SxSpeciesData.h>

/** \brief Container for pseudo potentials

    \b SxPseudoPot = S/PHI/nX Pseudo Potential Container

    ....

    \author Sixten Boeck */
class SX_EXPORT_DFT SxPseudoPot : public SxSpeciesData
{
   public:
            /** Number of species in the unit cell */
      int                                      nSpecies;
      /** Radius of artificial Gaussian charge */
      SxArray<double>                          rGauss;          // :iSpecies
      /** Maximum angular momentum used in pseudopotential. Please mind,
          this is the maximum component! Hence a for loop for over all
          components would include ( <= ) lMax!
          \verbatim
             for (int l=0; l <= lMax(iSpecies); l++)  {
                ...
             }   
          \endverbatim
        */
      SxVector<Int>                            lMax;            // :iSpecies
      /// Get max. l of all species
      int getLMax () const;
      /** Local compontent in pseudopotential */
      SxArray<int>                             lLoc;            // :iSpecies
      /// \brief Whether to use real-space projectors for certain species
      SxArray<bool>                            realSpace;
      /** Name of pseudopotential file */
      SxArray<SxString>                        pseudoPotFiles;  // :iSpecies
      /** The pseudopotential data */
      SxArray<SxDiracVec<TReal8> >             rad; // :iSpecies,:r
      SxArray<SxArray<SxDiracVec<TReal8> > >   pseudoPsi; // :iSpecies,:l,:r
      SxArray<SxArray<SxDiracVec<TReal8> > > getPseudoPsi ();
      SxDiracVec<TReal8> getPseudoPsi (int is, int l);
      SxArray<SxArray<SxDiracVec<TReal8> > >   pseudoPot; // :iSpecies,:l,:r
      /** Obsolete. Occupency numbers used in LCAO initialization */
      SxArray<SxArray<bool> >                  pseudoFocc; // :iSpecies,:l
      /** Occupency numbers used in LCAO initialization */
      SxArray<SxArray<SxVector<TReal8> > >     foccAtomicRho;  // :iSpecies,:l,:iSpin
      bool                                     nlcc;
      SxArray<SxArray<SxList<Real8> > >        radCore;    // :iSpecies,:2,:r
      // TODO: radCore needs to have an own logDr 
      /** The logDr parameter */
      SxArray<Real8>                           logDr;     // :iSpecies

      
      SxPseudoPot ();
      SxPseudoPot (const SxSymbolTable *);

   public:

      ~SxPseudoPot ();

      /** \brief re-confines pseudo-wavefunctions (without affecting
                 pseudo-potentials) by multiplying them by a Gaussian
                 distribution, and then normalizes them again
          \param sigma standard deviation of the Gaussian distribution. Thus,
                       increasing sigma reduces the effect of this routin */
      void reConfinePsi (double sigma);

      /// Check if real-space projectors are used
      bool useRealSpace () const {
         for (int is = 0; is < nSpecies; ++is) 
            if (realSpace(is)) return true;
         return false;
      }

      /** \brief Get initial density for single atom in G-basis
          @param G the G basis
          @param iSpecies which species
          @param iSpin which spin channel (-1 = sum)
        */
      SxDiracVec<TPrecCoeffG>
      getAtomRhoG(const SxGBasis &G, int iSpecies, int iSpin = -1) const;

   protected:
       /** \todo to be changed completely */
      void initPseudoPotentials ();

      void initPseudoPotSiesta (SxArray<SxString> psiFiles); 

//      void write (SxBinIO &);
//      void read  (SxBinIO &);

      void print () const;

};

#endif /* _SX_PSEUDO_POT_H_ */
