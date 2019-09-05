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
#ifndef _SX_TB_ATOM_ATOM_H_
#define _SX_TB_ATOM_ATOM_H_

#include <SxAtomicStructure.h>
#include <SxTBSpecies.h>

/** \brief Tight Binding Hamiltonian and overlap between orbitals of two atoms 

    \b SxTBAtomAtom = S/PHI/nX Tight Binding Atom - Atom interaction

    \author Hazem Abu-Farsakh, h.farsakh@mpie.de
*/
class SxTBAtomAtom 
{
   public:

   /** constructor */
   SxTBAtomAtom ();

   /** destructor */
   ~SxTBAtomAtom ();
   
   /** functions */
   /** \brief reads slater-kostrer file 

       reads epsAtom (only in case of the same species), and   
       the slater-koster parameters (Vss-sigma, Vsp-sigma, ... etc,
       Sss-sigma, ... etc), and coefficients used for calclating repulsive energy
       and other related details 
       \param iS            first Species 
       \param jS            second Species
       \param inSpeciesData tb species data 
       \param noRepul       true will not read repulsive coefficients from file, 
                            and will cause to set repulsive energy and its 
                            gradient to 0 */
   void  init (int iS, int jS, const SxTBSpecies &inSpeciesData, 
                               const bool        noERepul = false); 

   /** \brief calculates the atom - atom Hamiltonian and Overlap matrices 
       \param dVecIn  distance vector between the two atoms */
   void   setDistVec (const SxVector3<Double> &dVecIn);

   /** \brief calculates the atom - atom Hamiltonian and Overlap matrices 
              in case if the two atoms are the same 
              
       (set the sub matrices in the diagonal part of the complete matrix)*/
   void  setDiagElements ();

   /** \brief calculates the atom - atom Hamiltonian and Overlap matrices 
              in case if the two atoms are different. 
              
         setDistVec should be used before calling this function */
   void  setElements ();

   /** \brief 
       \return (iOrb, jOrb) element of the atom-atom Hamiltonian matrix*/
   double getHam (int iOrb, int jOrb) const;    

   /** \brief 
       \return (iOrb, jOrb) element of the atom-atom Overlap matrix*/
   double getOvl (int iOrb, int jOrb) const;    

   /** \brief 
       \return atom-atom Hamiltonian matrix*/
   SxMatrix<Double> getAtomAtomHam () const;

   /** \brief 
       \return atom-atom Overlap matrix */
   SxMatrix<Double> getAtomAtomOvl () const;

   /** \brief 
       \return the total number of orbitals of the first atom */
   int getNOrb1 () const;

   /** \brief 
       \return the total number of orbitals of the second atom */
   int getNOrb2 () const;

   /** \brief get number of valence electrons 
       \return the number of valence electrons. Can be used ONLY if the two
        atoms are of the same species */
   int getNElect () const;

   /** \brief 
       \return the atom - atom repulsive energy */
   double getERepulsive (double distIn) const;

   /** \brief 
       \return the atom - atom repulsive energy gradient */
   double getRepulsiveEGrad (double distIn) const;

   /** \brief 
       \return the cutoff radius for neighbors, based on the maximum distance
        between atoms in the slater-koster files */
   double getNeighborsCutoff () const;

   /** \brief 
       \return the Hubbard parameter (or called chemichal hardness 
       for spin-unpolarized atom), for a certain orbital.
       works only if the two atoms are of the same species type 
       because they exist in the homo-atomic slater-koaster files. 
       \param lOrb  orbital quantum number*/
   double getHubbardU (int iOrb) const;

   int    lMax1, lMax2, lOrb1, lOrb2;
   bool   isSameSpecies;
   double dist, dist2;

   bool   withoutERepulsive;

   protected:

   int    nCol, nPoints;
   double fraction, delta, cutoff;
   int    index;
   
   /** \brief atom - atom Hamiltonain */
   SxMatrix<Double>  ham;

   /** \brief atom - atom Overlap */
   SxMatrix<Double>  ovl;

   /** \brief atom - atom slater-koster parameters for the Hamiltonian */
   SxMatrix<Double>  slkoHam;

   /** \brief atom - atom slater-koster parameters for the Overlap */
   SxMatrix<Double>  slkoOvl;

   /** \brief number of intervals to read in the repulsive energy coefficients*/
   int     nLines;

   /** \brief number of valence elctrons */
   int     nValElect;

   /** \brief distance cutoff for the repulsive energy */
   double  dCut;

   /** \brief repulsive enerrgy data type provided in sk-files*/
   SxString repType; 

   /** \brief stored as data on a grid of distances. (index)(distance,eRep) */
   SxArray<SxVector<Double> > repData;

   /** \brief parameters for a function to calculate the repulsive energy, when
              the distance between atoms is less than a certain distancce */
   SxVector<Double> func;
   
   /** \brief intervals of distances (for the repulsive enerrgy) */
   SxMatrix<Double> interval;

   /** \brief coefficients of a function used to calculate the repulsive
              eneryg for each interval*/
   SxMatrix<Double> repulsiveCoef1;

   /** \brief coefficients of the last distance interval */
   SxVector<Double> repulsiveCoef2;

   /** \brief dCos are the directional cosines, dSin are the 
              corresponding sines*/
   SxVector3<Double> dCos, dSin;

   /** \brief atomic eigenvalues 
     
     used only if the two atoms are of the same species type */
   SxVector<Double>  epsAtom;

   /** \brief hubbard parameter
     
     used only if the two atoms are of the same species type. :(orbital)*/
   SxVector<Double>  hubbardU;
   
   /** \brief Hamiltonian slater-koster parameters 

        vCoeff(0) = Vss-sigma, vCoeff(1) = Vsp-sigma, vCoeff(2) = Vps-sigma, 
        vCoeff(3) = Vpp-pi,    vCoeff(4) = Vpp-sigma */
   SxVector<Double>  vCoeff;

   /** \brief Overlap slater-koster parameters 

        sCoeff(0) = Sss-sigma, sCoeff(1) = Ssp-sigma, sCoeff(2) = Sps-sigma, 
        sCoeff(3) = Spp-pi,    sCoeff(4) = Spp-sigma */
   SxVector<Double> sCoeff;

   /** \brief name of sk file*/ 
   SxString fileName;
   SxString fileName1;
   SxString fileName2;
   
   /** functions */
   /** \brief calculates and returns the Hamiltonian or the Overlap
              matrix elements for atom - atom interaction  

       \return the Hamiltonian or the Overlap matrix elements 
               for atom - atom interaction 
       \param coeff  slater-koster parameters vCoeff (in case of Hamiltonian)
                     or sCoeff (in case of overlap)*/ 
   SxMatrix<Double> getSubMatrix (const SxVector<Double> &coeff) const;

   /** \brief 
       \return the total number of orbitals 
       \param  lMax maximum angular quantum number */
   int getlOrb (int lMax) const;
   
   /** \brief spline function for a uniform grid 
       \return the values of the required slater-koster parameters 
               for a certain distance betweeen the interacting atoms 
       \param slko slater-koster parameters
       \param ind  index of the distance location in the file 
                    (ind <= nPoints-2)
       \param frac fraction of the location between ind and ind+1
                    (1. >= frac >= 0.)*/
   SxVector<Double> spline (const SxMatrix<Double> &slko, 
                            int ind, double frac) const;

   /** if the file that contains the first compination of species can't 
       be found, this bool is used to transpose the slater-koster matrix 
       elements when calculating the H and S matrices from slater-koster
       file of the opposite compination */
   bool isTransposed;
};


#endif /* _SX_TB_ATOM_ATOM */
