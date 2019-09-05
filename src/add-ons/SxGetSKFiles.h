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
#ifndef _SX_GET_SK_FILES_H_
#define _SX_GET_SK_FILES_H_

#ifdef WITH_TB
#include <SxConfig.h>
#include <SxAtomicStructure.h>
#include <SxExt.h>

/** \brief calculates Tight Binding Slater-Koster data files 

    \b SxGetTBSKFiles = S/PHI/nX Tight Binding Slater-Koster Files generator 

    This is not a final work. It is supposed to give us a feeling. 
    
    \author   Hazem Abu-Farsakh, hazem@phys.upb.de */

class SX_EXPORT_EXT SxGetSKFiles 
{
   public:

      SxGetSKFiles (const SxSymbolTable *tableIn,
                    double maxDIn, double stepIn, double sigmaIn);

      ~SxGetSKFiles ();

      /** \brief check the number of atoms and the cell, and set initial 
                 atomic positions */
      void setInitialStr ();

      /** \brief update structure and basis and calculate 
           Hamiltonian and Overlap*/
      void calculate (const SxSymbolTable *tableIn);

      /** \brief update structure and set Hamiltonian and Overlap to a number*/
      void setTo (const double num);

      /** \brief used to calculate atomic eigenvalues */
      void calcEigVals (const SxSymbolTable *tableIn);

      /** \brief print file headers */
      void printHeaders (FILE *fp);

      /** \brief update file */
      void updateFile (FILE *fp);

      /** \brief the atomic structure information */
      SxAtomicStructure structure;

      /** \brief number of species */
      int nSpecies;

      /** \brief the grid spacing between two structural steps */
      double strStep;

      /** \brief the number of structural steps */
      int nStrSteps;

      /** \brief maximum distance to reach between the two atoms (Bohr) */
      double maxDist;

      /** \brief the standard deviation (in Bohr) of a Gaussian distribution 
           that is used to re-confine the pseudo wavefunctions more.
           Thus, increasing sigma reduces the effect of this Gaussian. */
      double sigma;

      /** \brief object of the SxTBSKData class */
      SxTBSKData SKData;

      /** \brief Hamiltonian elements. (coupling type)(iOrb, jOrb) */
      SxArray<SxMatrix<Double> >  H;

      /** \brief Overlap elements. (coupling type)(iOrb, jOrb) */
      SxArray<SxMatrix<Double> >  S;

      /** \brief orbitals (l quantum numbers) per species. 
                 :(iSpecies)(iOrbital) */
      SxArray<SxList<int> > lOrbs;

      /** \brief elements name. :(iSpecies)*/
      SxArray<SxString>     elementName;

      /** \brief chemical symbols. :(iSpecies)*/
      SxArray<SxString>     chemName;

      /** \brief valence charges. :(iSpecies)*/
      SxArray<double>       valenceCharge;

      /** \brief set of atomic eigenvalues. :(iSpecies)(il QN)*/
      SxArray<SxVector<Double> > eigVals;
};
#endif /* WITH_TB */
#endif /* _SX_GET_SK_FILES_H_ */
