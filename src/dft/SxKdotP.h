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
#ifndef _SX_KDOTP_H_
#define _SX_KDOTP_H_

#include <SxRBasis.h>
#include <SxBinIO.h>
#include <SxVector3.h>
#include <SxRho.h>
#include <SxHamiltonian.h>
#include <SxFermi.h>
#include <SxTypes.h>

/** \brief SxKdotP n band Hamiltonian

    \b SxKdotP = SFHIngX implementation of k dot p theory
    for general, user-defined n-band Hamiltonians

    \author Oliver Marquardt, marquardt@mpie.de
  */
class SX_EXPORT_DFT SxKdotP : public SxHamiltonian,
                   public SxRho
{
   public:
   
   SxKdotP ();
   /// Constructor
   SxKdotP (const SxPWSet &,
               const RhoR &);
   
   PsiG psiIn;
   enum Preconditioner { Payne, Arias };

   /// Destructor
   virtual ~SxKdotP () {/*empty*/}

   static const SxSymbolTable * getHamiltonianGroup (const SxSymbolTable *);

   short nEmptyStates;
   const SxFermi *fermiPtr;
   const SxPWSet *wavesPtr;
   
   virtual void read (const SxSymbolTable *); // --- read input parameter
   virtual void update (const SxFermi &);
   virtual void compute (const SxFermi &, bool, bool);
   virtual SxRho &getRho ();
   virtual void normalizeRho();
   virtual void validateNStates (int);

   virtual void computeRho (const SxPsiSet &, const SxFermi &); 
   virtual PrecEnergy getEnergy (const SxPsiSet &, const SxFermi &);
   virtual void writeRho (const SxString &) const;
   virtual void readRho (const SxBinIO &);
   virtual PrecEnergy getETrial ();
   SxMeshR returnMesh (SxString inString);
   virtual void printEnergies () const;
   virtual double getEnergy () const { return 0.; }

   virtual int getNEmptyStates ();
 
   // --- calculate d^2/dx_i dx_j etc.
   SxArray<SxArray<PsiR> > secondDerivative (const PsiG &);
   // --- calculate d/dx_i
   SxArray<PsiR> firstDerivative (const PsiG &);
   SxArray<PsiR> rFirstDerivative (const PsiR &);
   void formDerivatives ();
                     
   PsiG operator* (const PsiG &);
   /// \brief Apply NxN Hamiltonian
   // hook in for SxHamiltonian
   virtual PsiG apply (const PsiG &psi) const
   {
      // note: can operator* be turned into const?
      return const_cast<SxKdotP*>(this)->SxKdotP::operator* (psi);
   }
 
   PsiR zero, one, psiR;
   SxMeshR vChg;
   PsiG gZero, gOne;
   const SxGBasis *gBasisPtr;
   SxArray<PsiG> gVecs, gG2;
   SxArray<SxMeshR> materials;
   SxArray<SxMeshR> eIJ; // strain fields
   SxMeshR vP, vExt, outMesh, help, help2, chargePotential;


   // --- for bandstructure calculations
   SxVector3<Double> bsStart, bsEnd, k, wgt;
   SxList<SxVector3<Double> > bsPts; // list of all points of bandstructure
   bool bs, old, onlyBS, moreIO, accurateInterfaces, firstStep, derivatives;
   int stepsBs, counter, nCharges;
   SxString bsMat, inFile;
   /** \brief
      used to calculate <psiI|H|psiJ> for optical momentum matrix
      */    
//   double psiIHPsiJ (const PsiG &, const PsiG &);
   void writeMeshASCII (const SxString &name, const PsiR &data) const;

   void showBandstructure(SxString, SxVector3<Int>);
   /** \brief
     Evaluate parse tree.
     */
   PsiR evaluateTree(int, int, int, bool);
   /** \brief
     Evaluate parse tree for a band structure plot.
     */
   SxComplex16 evaluateTreeBS(int, int, int, int, SxVector3<Double>);

   /** \brief
     Transform Hamiltonian element to parse tree
     */
   void buildTree(int, int, SxString, SxString);
   bool containsOperator(int, int, int);
   bool isRealNumber(int, int, int);
   void printElement(int, int, int);
   void countOperatorMultiplications(SxString);
   void determineDerivatives(int);

   // resolve expression and return position of right branch in tree list
   int resolveExpr(SxString, int, int);
   int resolveOpKey(SxString);
   /** \brief
     A specified parameter can be separately extracted in an .sxb file,
     e.g. for viewing using the Phinax viewer.
     */
   int outputPar, nOpMult, iMult;
   SxString outPar;
   SxArray<int> opKey;
   SxArray<SxString> opMultStr;
   SxArray<PsiR> kIPar, kJPar, kIJPar, par, kKPar;

   /** \brief
     Routines that mark elements as operator or potential.
     To make sure that no element in the parse tree is
     counted twice (i.e. double multiplication of the
     element with the wave function, unnecessary marks
     are removed after parsing the tree.
     */
   void correctOperators(int, int);
   void removeOpMarks(int, int, int);
   /** \brief
     Determine wether element is an operator or a potential term
     */
   int whatIsElement(int, int, int);
   int firstOutsideBracket(SxString, char);
   /** \brief
     The parse tree: expression, left branch, right branch
     */
   SxArray<SxArray<SxList<SxString> > >expression;
   SxArray<SxArray<SxList<int> > > leftPtr, rightPtr;
   /** \brief
     routines search for unknown expressions in the hamiltonian
     and replace these if given previously in the input file
     */
   SxString containsUnknown(SxList<SxString>);
   SxString replaceUnknown(SxString, SxString);
   /** \brief
     create Hamiltonian matrix for band structure plot
     */
   SxMatrix<Complex16> hMatrixBS(SxVector3<Double> k, SxVector3<Int> r);

   SxArray<SxMeshR> pMem; // for speed-optimized calculation: parameters in memory
   bool speedOpt;
   SxMeshR parameters (int iParam);
   int iEpsR;
   SxMeshR epsilonR, V, totalCharge;
   SxMatrix<Double> matParam; // all parameters of all materials
   SxMatrix<Double> bowParam; // all bowings of involved materials
   int nMat, nParam, nComp, precMaterial, gSize, rSize, bsMatIdx;
   SxArray<SxArray<PsiR> > sdC;
   SxArray<PsiR> fdC;
   SxList<SxString> matNames; // names of involved materials
   SxList<SxString> paramNames; // names of involved parameters
   SxList<SxString> chargeFiles; // names of external charge files

   SxList<double> paramVals, chargePrefactors;
   SxList<double> bowingVals; // bowing for material pairs
   SxMap<SxString,int> paramMap; // each parameter gets an integer index

   PsiG D;
   
   PrecEnergy eKin, eExt, eTrial;
   PsiR H_Psi;
   inline const SxPWSet &getWavesRef () const {
      SX_CHECK (wavesPtr);
      return *wavesPtr;
   }
 
   virtual void set (const SxPsiSet &, const SxFermi &);
   /** \brief The general n-band k.p-Hamiltonian
    */
   PsiR hamNxN (const PsiG &);
   PsiR returnAccurate(int);
   
   virtual SxDiracVec<TPrecCoeffG::TReal> preconditioner (const PsiG &psi) const {
      return preconditioner (psi, Payne);
   }

   SxDiracVec<TPrecCoeffG::TReal> 
      preconditioner (const PsiG &, 
            Preconditioner type) const;
};

namespace Timer {
   enum KPNxNHamTimer {
      sndDer,
      fstDer,
      NxNham,
      init,
      evalTree,
      derivatives,
      kk2s,
      extField,
      complexNr,
      doubleNr,
      known,
      constant
   };
}

SX_REGISTER_TIMERS (Timer::KPNxNHamTimer)
{
   using namespace Timer;
   regTimer (sndDer,      "second derivative");
   regTimer (fstDer,      "first derivative");
   regTimer (NxNham,      "Hamiltonian routine");
   regTimer (init,        "Hamiltonian init");
   regTimer (evalTree,    "Evaluate tree");
   regTimer (derivatives, "Set up derivatives");
   regTimer (kk2s,        "kk2s");
   regTimer (extField,    "ext");
   regTimer (complexNr,   "complex");
   regTimer (doubleNr,    "double");
   regTimer (known,       "known");
   regTimer (constant,    "constant");
}

#endif /* _SX_KDOTP_H_ */
