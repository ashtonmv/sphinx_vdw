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

#ifndef _SX_AOMAT_TB_H_
#define _SX_AOMAT_TB_H_

#include <SxExt.h>
#include <SxAtomicStructure.h>
#include <SxPotential.h>
#include <SxSpeciesData.h>
#include <SxGBasis.h>
#include <SxNeighbors.h>
#include <SxCLI.h>
#include <SxYlm.h>
#include <SxKPoints.h>

enum AoMatTBTimer { HPsiTime, SPsiTime, InitTimer, FTimer, GradTimer, orthoTimer, computeTimer, CholeskyTimer};
SX_REGISTER_TIMERS(AoMatTBTimer)
{
   regTimer (HPsiTime,  "H|Psi>");
   regTimer (SPsiTime,  "S|Psi>");
   regTimer (InitTimer, "Initialization");
   regTimer (FTimer,    "calcFunctional");
   regTimer (GradTimer, "calcGradient");
   regTimer (orthoTimer, "Orthonormalize");
   regTimer (computeTimer, "Compute");
   regTimer (CholeskyTimer, "Cholesky");
}

class SX_EXPORT_EXT SxAOMatTB 
{
   protected:
      /// orbitalTypeList per species
      SxArray<SxArray<int> > orbitalTypeList;

      /// List of neighbor idx for each atom
      SxArray<SxArray<int> > neighborIdx;

      /// List of neighbor coords
      SxArray<SxArray<Coord> > neighborCoord;

      /// List of Hamiltonian matrices for all possible neighbors indexed
      SxArray<SxMatrix<Double> > hamiltonians;

      /// List of Overlap matrices for all possible neighbours indexed
      SxArray<SxMatrix<Double> > overlaps;

      /// Index List for all neighbors of a specific atom
      SxArray<SxArray<int> > idxList; 

      /// Orientation of the atoms
      SxArray<int> orientations;

      /// Species data
      SxSpeciesData speciesData;

      /// actual kPoint
      mutable Coord kPoint;

      /// Generating neighbor list
      SxArray<Coord>  genNeighbors;
      /// Symmetry list
      SxArray<SymMat> syms;

      /// Map of neighbor origin: which genNeighbor (iTlAtom,iNeighbor)
      SxArray<SxArray<int> >  genNeighborId;
      /// Map of neighbor origin: which sym (iTlAtom,iNeighbor)
      SxArray<SxArray<int> >  genSymId;

   public:

      /// Constructor
      SxAOMatTB () { /*empty*/ };
      /// Constructor using input file and atomic strusture
      SxAOMatTB (const SxSymbolTable *tablePtr);
      /// Destructor
      ~SxAOMatTB () {/*empty*/ };

      void info() const;

      int getOffset (int iSpecies) const;
      /// gives number of all orbitals for this structure
      int getNOrbs () const;
      /// gives number of Orbitals for one atom of one Species
      int getNOrbs (int iSpecies) const;

      /// orbitalTypeList per species
      const SxArray<SxArray<int> > &getOrbitalTypeList() const {return orbitalTypeList;};
  
      /// Implementation of H |Psi>
      SxVector<Complex16> applyH (const SxVector<Complex16> &psi) const;
      
      /// Implementation of S |Psi>
      SxVector<Complex16> applyS (const SxVector<Complex16> &psi) const;

      /// set kPoint for Hamiltonian
      void setKPoint (const Coord &kPointIn) const {kPoint = kPointIn;};

      /// Atomic structure
      SxAtomicStructure structure;

      void directMin (const SxKPoints &kPoints);

      int findNeighbor (SxBinIO &io, int iSpecies, int jSpecies, Coord dist);
      SxMatrix3<Double> readR (SxBinIO &io, int iNeighbor);
      SxMatrix<Double>  readS (SxBinIO &io, int iNeighbor);
      SxMatrix<Double>  readH (SxBinIO &io, int iNeighbor);
};

class SX_EXPORT_EXT SxWaves  {
   
   protected:
      
   SxArray<SxVector<Complex16> > waves;
   SxConstPtr<SxAOMatTB> tbHamPtr;

   public:

   SxWaves () {/*empty*/};
   SxWaves (const SxWaves &in) 
   {
      tbHamPtr = in.tbHamPtr;
      int nStates = in.getNStates();
      waves.resize(nStates);
      for (int iState = 0; iState < nStates; iState++)  {
         waves(iState) = 1.0 * in(iState);
      }
   };

   SxWaves (int nStates, int dimIn, SxConstPtr<SxAOMatTB> tbHamPtrIn)
      :tbHamPtr(tbHamPtrIn)
   {
      SX_CHECK (nStates > 0);
      SX_CHECK (dimIn > 0);
      waves.resize(nStates);
      for (int iState = 0; iState < nStates; iState++)
         waves(iState).resize(dimIn);
   };
   SxWaves (SxArray<SxVector<Complex16> >&wavesIn, SxConstPtr<SxAOMatTB> tbHamPtrIn)
      :tbHamPtr(tbHamPtrIn)
   {
      waves = wavesIn;
   }
   ~SxWaves () {/*empty*/};

   void randomize ()
   {
      for (int iState = 0; iState < waves.getSize(); iState++)  {
         waves(iState).randomize();
         waves(iState).normalize();
      }
   };

   void occupySpecies (int iSpecies)
   {
      SxVector<Complex16> filter (tbHamPtr->getNOrbs ());
      int offset = tbHamPtr->getOffset(iSpecies);
      int nSpeciesOrbs = tbHamPtr->getNOrbs(iSpecies) * tbHamPtr->structure.getNAtoms(iSpecies);
      if (getNStates() > nSpeciesOrbs) {
         cout << "To many wave states, orthogonalization would produce zero vectors." << endl;
         SX_EXIT;
      }
      filter.set(0.0);
      for (int i = 0; i < nSpeciesOrbs; i++) filter(offset + i) = 1.0;
      for (int iState = 0; iState < getNStates(); iState++) waves(iState) *= filter;
   }

   void orthonormalize ()
   {
      SX_CHECK(getNStates () > 0);
      waves(0).print();
      tbHamPtr->applyS(waves(0)).print();
      waves(0) /= sqrt( dot(waves(0),tbHamPtr->applyS(waves(0)) ) );
      for (int iState = 1; iState < waves.getSize(); iState++)  {
         for (int jState = 0; jState < iState; jState++)   {
            waves(iState) -= dot(waves(jState),tbHamPtr->applyS(waves(iState))) 
                           * waves(jState);
         }
         waves(iState) /= sqrt( dot(waves(iState),tbHamPtr->applyS(waves(iState)) ) );
      }
   };

   bool checkNorm () const
   {
      bool result = true;
      for (int iState = 0; iState < getNStates(); iState++)  {
         SxComplex16 norm = dot(waves(iState),tbHamPtr->applyS(waves(iState)));
         cout << norm << endl;
         if (fabs(norm.abs() - 1) > 1e-12) result = false;
      }
      return result;
   };

   SxVector<Complex16>& operator()(int iState) { return waves(iState);};
   const SxVector<Complex16>& operator()(int iState) const { return waves(iState);};

   void operator= (const SxWaves &in)   
   {
      tbHamPtr = in.tbHamPtr;

      int nStates = in.getNStates ();
      waves.resize(nStates);
      for (int iState = 0; iState < nStates; iState++)  {
         waves(iState) = 1.0 * in.waves(iState);
      }
   };

   SxWaves operator+ (const SxWaves &in) const
   {
      SX_CHECK (getNStates () == in.getNStates ());
      SX_CHECK (getDim () == in.getDim ());
      
      SxArray<SxVector<Complex16> > result (getNStates ());
      for (int iState = 0; iState < getNStates (); iState++)  {
         result(iState) = waves(iState) + in(iState);
      }

      return SxWaves(result,tbHamPtr);
   };

   SxWaves operator- (const SxWaves &in) const
   {
      SX_CHECK (getNStates () == in.getNStates ());
      SX_CHECK (getDim () == in.getDim ());
      
      SxArray<SxVector<Complex16> > result (getNStates ());
      for (int iState = 0; iState < getNStates (); iState++)  {
         result(iState) = waves(iState) - in(iState);
      }

      return SxWaves(result,tbHamPtr);
   };

   SxWaves operator* (const SxVector<Double> &in) const
   {
      SX_CHECK(in.getSize () == getNStates ());
      SxArray<SxVector<Complex16> > result (getNStates ());
      for (int iState = 0; iState < getNStates (); iState++)  {
         result(iState) = waves(iState) * in(iState);
      }

      return SxWaves(result,tbHamPtr);
   };

   SxWaves operator/ (const SxVector<Double> &in) const
   {
      SX_CHECK(in.getSize () == getNStates ());
      SxArray<SxVector<Complex16> > result (getNStates ());
      for (int iState = 0; iState < getNStates (); iState++)  {
         result(iState) = waves(iState) / in(iState);
      }

      return SxWaves(result,tbHamPtr);
   };

   int getNStates () const { return (int)waves.getSize();};
   int getDim () const { return (int)waves(0).getSize();};

};

inline SxWaves operator*(const SxVector<Double> &in, const SxWaves &waves) { return waves * in;}

// Auxillary class for Minimization
class SX_EXPORT_EXT SxHeS {
   
   protected:
   SxConstPtr<SxAOMatTB> hamTBPtr;
   SxWaves psi;
   SxArray<SxVector<Complex16> > HPsi;
   SxArray<SxVector<Complex16> > SPsi;
   SxArray<SxVector<Complex16> > HeS2Psi;
   SxVector<Double> epsilon;
   SxMatrix<Complex16> kappa;
   
   // Constructor only for operator usage
   SxHeS (SxConstPtr<SxAOMatTB> hamTBPtrIn, int nStates)
   {
      SX_CHECK (hamTBPtrIn.getPtr());
      
      hamTBPtr = hamTBPtrIn;
      HPsi.resize(nStates);
      SPsi.resize(nStates);
      HeS2Psi.resize(nStates);
      epsilon.resize(nStates);
      kappa.reformat(nStates,nStates);
   };

   public:
   //Constructor
   SxHeS () {/*empty*/};
   SxHeS (const SxHeS &in) {setup (in.hamTBPtr, in.psi, in.epsilon);};
   SxHeS (SxConstPtr<SxAOMatTB> hamTBPtrIn, const SxWaves &psiIn, const SxVector<Double> &epsilonIn) 
   {
      setup (hamTBPtrIn, psiIn, epsilonIn);
   };

   ~SxHeS () {/*empty*/};

   void setup (SxConstPtr<SxAOMatTB> hamTBPtrIn, const SxWaves &psiIn, 
               const SxVector<Double> &epsilonIn)
   {
      SX_CHECK (hamTBPtrIn.getPtr());
      SX_CHECK (psiIn.getNStates () > 0);
      SX_CHECK (hamTBPtrIn->getNOrbs () == psiIn.getDim (),
                hamTBPtrIn->getNOrbs (), psiIn.getDim ());
      hamTBPtr = hamTBPtrIn;
      psi = psiIn;
      epsilon = 1.0 * epsilonIn;
      int nStates = psi.getNStates ();
      HPsi.resize(nStates);
      SPsi.resize(nStates);
      HeS2Psi.resize(nStates);
      kappa.reformat(nStates,nStates);

      compute ();
      updateKappa ();
   }

   int getNStates () const {return psi.getNStates ();};

   const SxWaves& getPsi() const {return psi;};
   const SxVector<Complex16>& getSPsi (int iState) const {return SPsi(iState);};
   const SxVector<Complex16>& getPsi (int iState) const {return psi(iState);};
   const SxVector<Complex16>& getHeS2Psi (int iState) const {return HeS2Psi(iState);};

   SxMatrix<Complex16> getCholesky ()
   {
      int nStates = getNStates ();
      SxMatrix<Complex16> S (nStates, nStates);

      for (int iState = 0; iState < nStates; iState++)  {
         S(iState,iState) = dot(psi(iState), SPsi(iState));
         for (int jState = iState + 1; jState < nStates; jState++)  {
            S(iState,jState) = dot(psi(iState), SPsi(jState));
            S(jState,iState) = S(iState,jState).conj ();
         }
      }

      SxMatrix<Complex16> result = S.choleskyDecomposition ();

      return result;
   }

   void normalize ()
   {
      for (int iState = 0; iState < getNStates (); iState++)  
         normalize (iState);
   }

   void normalize (int iState)
   {
      double factor = sqrt(dot(psi(iState),SPsi(iState)).re);

      psi(iState)     /= factor;
      HPsi(iState)    /= factor;
      SPsi(iState)    /= factor;
      HeS2Psi(iState) /= factor;

   }

   void compute ()
   {
      SX_START_TIMER(computeTimer);
      int nStates = psi.getNStates ();
      for (int iState = 0; iState < nStates; iState++)  {
         HPsi(iState) = hamTBPtr->applyH(psi(iState));
         SPsi(iState) = hamTBPtr->applyS(psi(iState));
         SxVector<Complex16> HeSPsi = HPsi(iState) - epsilon(iState) * SPsi(iState);
         HeS2Psi(iState) = hamTBPtr->applyH(HeSPsi) - epsilon(iState) * hamTBPtr->applyS(HeSPsi);
         SX_CHECK (dot(psi(iState), HPsi(iState)).im < 1e-14,
                   dot(psi(iState), HPsi(iState)).im);
         SX_CHECK (dot(psi(iState), SPsi(iState)).im < 1e-14,
                   dot(psi(iState), SPsi(iState)).im);
         SX_CHECK (dot(psi(iState), HeS2Psi(iState)).im < 1e-14,
                   dot(psi(iState), HeS2Psi(iState)).im);
      }

      SX_STOP_TIMER(computeTimer);

   }

   double getFunctional ()
   {
      SX_START_TIMER(FTimer);
      SxComplex<double> result = 0.0;
      int nStates = psi.getNStates ();
      for (int iState = 0; iState < nStates; iState++)  {
         SX_CHECK((HPsi(iState) - hamTBPtr->applyH(psi(iState))).norm () < 1e-14,
                  (HPsi(iState) - hamTBPtr->applyH(psi(iState))).norm ());
         SX_CHECK((SPsi(iState) - hamTBPtr->applyS(psi(iState))).norm () < 1e-14,
                  (SPsi(iState) - hamTBPtr->applyS(psi(iState))).norm ());
         result += dot(psi(iState),HeS2Psi(iState));
         SX_CHECK (fabs(dot(psi(iState),HeS2Psi(iState)).im) < 1e-14,
                   dot(psi(iState),HeS2Psi(iState)).im);
      }
    
      SX_CHECK (fabs(result.im) < 1e-14, fabs(result.im));
      SX_STOP_TIMER(FTimer);

      return result.re;
   }

   SxWaves getGradient () const
   {
      SX_START_TIMER(GradTimer);
      int nStates = psi.getNStates ();
      SxArray<SxVector<Complex16> > grad (nStates);
      
      SxArray<SxVector<Complex16> > SPsiSum (nStates);
      for (int iState = 0; iState < nStates; iState++)  {
         SPsiSum(iState).resize(SPsi(iState).getSize());
         SPsiSum(iState).set(0.0);
         for (int jState = 0; jState < nStates; jState++)  {
            SPsiSum(iState) += kappa(iState,jState) * SPsi(jState);
         }
      }

      for (int iState = 0; iState < nStates; iState++)  {
         grad(iState) = HeS2Psi(iState) - SPsiSum(iState);
      }

      SxWaves result (grad, hamTBPtr);
      SX_STOP_TIMER(GradTimer);

      return result;
   };

   SxMatrix<Complex16> setEnergyByCholesky ()
   {
      SX_START_TIMER(CholeskyTimer);
      cout << "START ENERGY BY CHOLESKY" << endl;
      int nStates = getNStates ();
      SxMatrix<Complex16> H (nStates, nStates);

      for (int iState = 0; iState < nStates; iState++)  {
         H(iState,iState) = dot(psi(iState), HPsi(iState));
         for (int jState = iState + 1; jState < nStates; jState++)  {
            H(iState,jState) = dot(psi(iState), HPsi(jState));
            H(jState,iState) = H(iState,jState).conj ();
         }
      }
      SxMatrix<Complex16> LInv = getCholesky ().inverse ();
      H = LInv ^ H ^ LInv.adjoint ();

      SxMatrix<Complex16>::Eigensystem eig = H.eigensystem ();
      SxMatrix<Complex16> rot = LInv.adjoint () ^ eig.vecs;

      SxHeS orig = *this;

      epsilon = eig.vals.real ();
      
      for (int i = 0; i < nStates; i++)  {
         psi(i).set(0.0);
         HPsi(i).set(0.0);
         SPsi(i).set(0.0);
         for (int j = 0; j < nStates; j++)  {
            // Konvention in eig.vecs : Summation has to be over first index!
            psi(i)     += rot(j,i) * orig.psi(j);
            HPsi(i)    += rot(j,i) * orig.HPsi(j);
            SPsi(i)    += rot(j,i) * orig.SPsi(j);
         }
         SxVector<Complex16> HeSPsi = HPsi(i) - epsilon(i) * SPsi(i);
         HeS2Psi(i) = hamTBPtr->applyH(HeSPsi) - epsilon(i) * hamTBPtr->applyS(HeSPsi);
         SX_CHECK (fabs(dot(psi(i), HPsi(i)).im) < 1e-14,
                   dot(psi(i), HPsi(i)).im);
         SX_CHECK (fabs(dot(psi(i), SPsi(i)).im) < 1e-14,
                   dot(psi(i), SPsi(i)).im);
         SX_CHECK (fabs(dot(psi(i), HeS2Psi(i)).im) < 1e-14,
                   dot(psi(i), HeS2Psi(i)).im);
      }
      normalize();
      updateKappa ();

      SX_STOP_TIMER(CholeskyTimer);

      return rot;
   }

   SxVector<Double> getEnergy () const
   {
      int nStates = psi.getNStates ();
      SxVector<Double> result (nStates);
      for (int iState = 0; iState < nStates; iState++)  {
         result(iState) = getEnergy(iState);
      }
      return result;
   }

   double getEnergy (int iState) const
   {
      double result = dot(psi(iState),HPsi(iState)).re / dot(psi(iState),SPsi(iState)).re;

      return result;
   }

   void updateKappa ()
   {
      int nStates = psi.getNStates ();
      for (int iState = 0; iState < nStates; iState++)  {
         for (int jState = 0; jState < nStates; jState++)  {
            kappa(iState,jState) = dot(psi(jState),HeS2Psi(iState));
         }
      }
   }
   
   void updateEpsilon ()
   {
      setEnergyByCholesky ();
      updateKappa ();
   }

   void setEpsilon (const SxVector<Double> &in)
   {
      epsilon = 1.0*in;
      for (int iState = 0; iState < getNStates (); iState++)  {
         SxVector<Complex16> HeSPsi = HPsi(iState) - epsilon(iState) * SPsi(iState);
         HeS2Psi(iState) = hamTBPtr->applyH(HeSPsi) - epsilon(iState) * hamTBPtr->applyS(HeSPsi);
      }
      updateKappa ();
   }

   const SxVector<Double>& getEpsilon () const { return epsilon; }
   
   SxVector<Double> getResiduum () const
   {
      SX_CHECK (psi.getNStates () > 0);
      int nStates = psi.getNStates ();
      SxVector<Double> result(nStates);
      for (int iState = 0; iState < nStates; iState++)  {
         result(iState) = getResiduum(iState);
      }
      return result;
   }
   
   double getResiduum (int iState) const
   {
      return (HPsi(iState) - getEnergy(iState) * SPsi(iState)).norm ();
   }

   void orthonormalize (bool calcHeS2Psi)
   {
      SX_START_TIMER(orthoTimer);
      int nStates = getNStates ();
      normalize (0);
      for (int iState = 1; iState < nStates; iState++)  {
         for (int jState = 0; jState < iState; jState++)   {
            SxComplex<double> factor = dot(psi(jState),SPsi(iState));
            psi(iState) -= factor * psi(jState);
            HPsi(iState) -= factor * HPsi(jState);
            SPsi(iState) -= factor * SPsi(jState);
            if(!calcHeS2Psi) HeS2Psi(iState) -= factor * HeS2Psi(jState);
         }
         normalize(iState);
      }
      if (calcHeS2Psi)  {
         for (int iState = 0; iState < nStates; iState++)  {
            SxVector<Complex16> HeS = HPsi(iState) - epsilon(iState) * SPsi(iState);
            HeS2Psi(iState) = hamTBPtr->applyH(HeS) - epsilon(iState) * hamTBPtr->applyS(HeS);
         }
      }
      SX_STOP_TIMER(orthoTimer);
   }

   void operator=(const SxHeS &in)
   {
      hamTBPtr = in.hamTBPtr;
      psi = in.psi;
      int nStates = psi.getNStates ();
      epsilon = 1.0 * in.epsilon;
      kappa = 1.0 * in.kappa;
      HPsi.resize(nStates);
      SPsi.resize(nStates);
      HeS2Psi.resize(nStates);
      for (int iState = 0; iState < nStates; iState++)  {
         HPsi(iState) = 1.0 * in.HPsi(iState);
         SPsi(iState) = 1.0 * in.SPsi(iState);
         HeS2Psi(iState) = 1.0 * in.HeS2Psi(iState);
      }
   }

   void operator-=(const SxHeS &in)
   {
      SX_CHECK (hamTBPtr.getPtr () == in.hamTBPtr.getPtr());
      SX_CHECK (psi.getDim () == in.psi.getDim());
      SX_CHECK (psi.getNStates () == in.psi.getNStates());
      SX_CHECK (fabs(in.epsilon.norm () - epsilon.norm ()) < 1e-12);
      
      int nStates = psi.getNStates ();
      for (int iState = 0; iState < nStates; iState++)  {
         psi(iState)     -= in.psi(iState);
         HPsi(iState)    -= in.HPsi(iState);
         SPsi(iState)    -= in.SPsi(iState);
         HeS2Psi(iState) -= in.HeS2Psi(iState);
      }
   }

   void operator+=(const SxHeS &in)
   {
      SX_CHECK (hamTBPtr.getPtr () == in.hamTBPtr.getPtr());
      SX_CHECK (psi.getDim () == in.psi.getDim());
      SX_CHECK (psi.getNStates () == in.psi.getNStates());
      SX_CHECK (fabs(in.epsilon.norm () - epsilon.norm ()) < 1e-12);
      
      int nStates = psi.getNStates ();
      for (int iState = 0; iState < nStates; iState++)  {
         psi(iState)     += in.psi(iState);
         HPsi(iState)    += in.HPsi(iState);
         SPsi(iState)    += in.SPsi(iState);
         HeS2Psi(iState) += in.HeS2Psi(iState);
      }
   }

   void operator*=(const SxVector<Double> &skalar)
   {
      int nStates = psi.getNStates ();
      for (int iState = 0; iState < nStates; iState++)  {
         psi(iState) *= skalar(iState);
         HPsi(iState) *= skalar(iState);
         SPsi(iState) *= skalar(iState);
         HeS2Psi(iState) *= skalar(iState);
      }
   }

   void operator/=(const SxVector<Double> &skalar) 
   {
      int nStates = psi.getNStates ();
      for (int iState = 0; iState < nStates; iState++)  {
         psi(iState) /= skalar(iState);
         HPsi(iState) /= skalar(iState);
         SPsi(iState) /= skalar(iState);
         HeS2Psi(iState) /= skalar(iState);
      }
   }

   SxHeS operator*(const SxVector<Double> &skalar) const
   {
      int nStates = getNStates ();
      SxHeS result(this->hamTBPtr, nStates);

      result.epsilon = 1.0 * epsilon;
      result.kappa = 1.0 * kappa;
      result.psi = skalar * psi;
      for (int iState = 0; iState < nStates; iState++)  {
         result.HPsi(iState) = skalar(iState) * HPsi(iState);
         result.SPsi(iState) = skalar(iState) * SPsi(iState);
         result.HeS2Psi(iState) = skalar(iState) * HeS2Psi(iState);
      }
      
      return result;
   }

   SxHeS operator/(const SxVector<Double> &skalar) const
   {
      int nStates = getNStates ();
      SxHeS result(this->hamTBPtr, nStates);

      result.epsilon = 1.0 * epsilon;
      result.kappa = 1.0 * kappa;
      result.psi = psi / skalar;
      for (int iState = 0; iState < nStates; iState++)  {
         result.HPsi(iState) = HPsi(iState) / skalar(iState);
         result.SPsi(iState) = SPsi(iState) / skalar(iState);
         result.HeS2Psi(iState) = HeS2Psi(iState) / skalar(iState);
      }
      
      return result;
   }

   SxHeS operator+(const SxHeS &in) const
   {
      SX_CHECK (in.hamTBPtr.getPtr () == hamTBPtr.getPtr ());
      SX_CHECK (psi.getDim () == in.psi.getDim());
      SX_CHECK (psi.getNStates () == in.psi.getNStates());
      SX_CHECK (fabs(in.epsilon.norm () - epsilon.norm ()) < 1e-12);
      
      int nStates = getNStates ();
      SxHeS result(this->hamTBPtr, nStates);

      result.epsilon = 1.0 * epsilon;
      result.kappa = 1.0 * kappa;
      result.psi = psi + in.psi;
      for (int iState = 0; iState < nStates; iState++)  {
         result.HPsi(iState) = HPsi(iState) + in.HPsi(iState);
         result.SPsi(iState) = SPsi(iState) + in.SPsi(iState);
         result.HeS2Psi(iState) = HeS2Psi(iState) + in.HeS2Psi(iState);
      }
      
      return result;
   }

   SxHeS operator-(const SxHeS &in) const
   {
      SX_CHECK (in.hamTBPtr.getPtr () == hamTBPtr.getPtr ());
      SX_CHECK (psi.getDim () == in.psi.getDim());
      SX_CHECK (psi.getNStates () == in.psi.getNStates());
      SX_CHECK (fabs(in.epsilon.norm () - epsilon.norm ()) < 1e-12);
      
      int nStates = getNStates ();
      SxHeS result(this->hamTBPtr, nStates);

      result.epsilon = 1.0 * epsilon;
      result.kappa = 1.0 * kappa;
      result.psi = psi - in.psi;
      for (int iState = 0; iState < nStates; iState++)  {
         result.HPsi(iState) = HPsi(iState) - in.HPsi(iState);
         result.SPsi(iState) = SPsi(iState) - in.SPsi(iState);
         result.HeS2Psi(iState) = HeS2Psi(iState) - in.HeS2Psi(iState);

         SX_CHECK (dot(result.psi(iState), result.HPsi(iState)).im < 1e-14,
                   dot(result.psi(iState), result.HPsi(iState)).im);
         SX_CHECK (dot(result.psi(iState), result.SPsi(iState)).im < 1e-14,
                   dot(result.psi(iState), result.SPsi(iState)).im);
         SX_CHECK (dot(result.psi(iState), result.HeS2Psi(iState)).im < 1e-14,
                   dot(result.psi(iState), result.HeS2Psi(iState)).im);
      }
      
      return result;
   }
};

inline SxHeS operator* (const SxVector<Double> &skalar, const SxHeS &in) {return in * skalar;}

#endif /* _SX_AOMAT_TB_H_ */
