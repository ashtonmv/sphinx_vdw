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

#ifndef _SX_KDOTP8X8_H_
#define _SX_KDOTP8X8_H_

#include <SxDFT.h>
#include <SxRBasis.h>
#include <SxBinIO.h>
#include <SxVector3.h>
#include <SxRho.h>
#include <SxHamiltonian.h>
#include <SxFermi.h>
#include <SxTypes.h>

/** \brief SxKdotP 8 band Hamiltonian

    \b SxKdotP8x8 = S/PHI/nX implementation of k dot p theory
    using a six band Hamiltonian
    from \ref Pryor97

    \author Oliver Marquardt, marquardt@mpie.de
  */
class SX_EXPORT_DFT SxKdotP8x8 : public SxHamiltonian,
                                 public SxRho
{
   public:
   
   SxKdotP8x8 ();// { }
   /// Constructor
   SxKdotP8x8 (const SxPWSet &,
               const RhoR &);
   
   PsiG psiIn;
   enum Preconditioner { Payne, Arias };

   /// Destructor
   virtual ~SxKdotP8x8 () {/*empty*/}

   static const SxSymbolTable * getHamiltonianGroup (const SxSymbolTable *);

   Real8 dOmega, eps; 
   double ekt, precFactor, f;
   int nEmptyStates, step;
   short Hamiltonian, nStates;
   const SxFermi *fermiPtr;
   const SxPWSet *wavesPtr;
   /** \brief used parameters:
     
       The Luttinger parameter \f$ \gamma_{1,2,3} \f$, valence band edges
       \f$ E_v\f$ and spin orbit coupling \f$ \Delta_{so} \f$ are different in
       two materials.
     */    
   
   SxMeshR deltaSpinOrbital;
   SxArray<SxMeshR> gamma;
   PsiR Tkin;
   SxArray<SxMeshR> delta;  // Delta1-3 in WZ 
   SxArray<SxArray<PsiR> > H_Strain; // ???
   SxMeshR recMassPar, recMassOrth, P1, P2;  // parameter set in WZ
   SxArray<SxMeshR> A, eVBOMod;
   SxMeshR A1,A2, L1, L2, M1, M2, M3, N1, N2, N3;
   SxMeshR reciprocalEffMass, externalPotential, eVBO, eCBO, eP, p0;
   SxMeshR pO, pP, externalField;
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
   virtual void printEnergies () const;
   virtual double getEnergy () const { return 0.; }

   virtual int getNEmptyStates (); // TODO: generalize
 
   SxArray<SxArray<PsiR> > secondDerivative (const PsiG &);
   SxArray<PsiR> firstDerivative (const PsiG &);
                     
   double interpolate(double, double, double, double);   
   PsiG operator* (const PsiG &);
   /// \brief Apply 8x8 Hamiltonian
   // hook in for SxHamiltonian
   virtual PsiG apply (const PsiG &psi) const
   {
      // note: can operator* be turned into const?
      return const_cast<SxKdotP8x8*>(this)->SxKdotP8x8::operator* (psi);
   }
   PsiR zero, psiR;
   bool strained, electrons, bs;
   const SxGBasis *gBasisPtr;
   SxArray<SxArray<PsiR> > sdC;
   SxArray<double> aP, gammaP;
   double aP1, aP2, rMassP;
   double kx,ky,kz;
   bool kxSet, kySet, kzSet;
   SxArray<PsiR> fdC;
   const SxGkBasis *Gk;
   SxArray<PsiG> gVecs;
   PsiG g2;
   PsiG D;
   double s32, s2, s3, s6;
   double p, bow;
   double mixing(double, double, double, double, double);
   SxVector3<Double> recMesh;
   
   PrecEnergy eKin, eExt, eTrial;
   PsiR H_Psi;
   inline const SxPWSet &getWavesRef () const {
      SX_CHECK (wavesPtr);
      return *wavesPtr;
   }
 
   virtual void set (const SxPsiSet &, const SxFermi &);
   /** \brief Zincblend Hamiltonian is calculated as follows:

       \f[ \hat{H} = \left[\begin{array}{cccccccc}
       A & 0 & V^\star & 0 & \sqrt{3} V & -\sqrt{2} U & -U & \sqrt{2}V^\star\\
       0 & A & -\sqrt{2} U & -\sqrt{3} V^\star & 0 & -V & \sqrt{2} V & U\\
       V & - \sqrt{2} U & -P + Q & -S^\star & R & 0
       & \sqrt{\frac{3}{2}} S & -\sqrt{2} Q\\
       0 & -\sqrt{3} V & -S & -P-Q & 0 & R & -\sqrt{2} R & \frac{1}{\sqrt{2}} S\\
       \sqrt{3} V^\star & 0 & R^\star & 0 & -P-Q & S^\star
       & \frac{1}{\sqrt{2}} S^\star & \sqrt{2} R^\star\\
       -\sqrt{2} U & -V^\star & 0 & R^\star & S & -P +Q
       & \sqrt{2}Q & \sqrt{\frac{3}{2}}S^\star\\
       -U & \sqrt{2}V^\star & \sqrt{\frac{3}{2}}S^\star & -\sqrt{2}R^\star
       & \frac{1}{\sqrt{2}}S & \sqrt{2} Q & -P-\Delta & 0\\
       \sqrt{2} V & U & -\sqrt{2}Q & \frac{1}{\sqrt{2}}S^\star & \sqrt{2}R
       &\sqrt{\frac{3}{2}} S & 0 & -P-\Delta
       \end{array}\right] \f]
       with:
       \f[ A = E_c - \frac{\hbar^2}{2m_0}
       (\partial_x^2 + \partial_y^2 + \partial_z^2)\f]
       \f[ P = -E_v-\gamma_1\frac{\hbar^2}{2m_0}(\partial_x^2 + \partial_y^2
       + \partial_z^2)\f]
       \f[ Q = \gamma_2\frac{\hbar^2}{2m_0}(\partial_x^2 + \partial_y^2
       - 2 \partial_z^2)\f]
       \f[ R = \sqrt{3}\frac{\hbar^2}{2m_0}[\gamma_2(\partial_x^2-\partial_y^2)
       -2 i \gamma_3\partial_x\partial_y]\f]
       \f[ S = -\sqrt{3}\gamma_3\frac{\hbar^2}{2m_0}\partial_z( \partial_x
       - i \partial_y)\f]
       \f[ U = \frac{-i}{\sqrt{3}}P_0\partial_z\f]
       \f[ V = \frac{-i}{\sqrt{6}}P_0(\partial_x -i\partial y)\f]
       where \f$ E_v\f$ is the unstrained local valence-band and \f$ E_c\f$
       the conduction band edge and \f$\gamma_i\f$ are the three modified
       Luttinger parameters:
       \f[ \gamma_1 = \gamma_1^L - \frac{E_p}{3E_g + \Delta}\f]
       \f[ \gamma_2 = \gamma_2^L - \frac{1}{2}\frac{E_p}{3E_g + \Delta}\f]
       \f[ \gamma_3 = \gamma_3^L - \frac{1}{2}\frac{E_p}{3E_g + \Delta}\f]
       where \f$Eg=E_c-E_v\f$ and
       \f[ E_p = \frac{2m_0 P^2}{\hbar^2}\f]
       with \f$ P_0\f$ being the coupling between conduction and valence bands.
     */
   PsiR ham8x8_ZB  (const PsiG &);

   /** \brief Wurtzite Hamiltonian (from PRB 54, 2491 (1996)) is calculated as follows:

     \f[ \hat{H} = \left[\begin{array}{cccccccc}
     T & 0 & -P_{xy} & P_z & P_{xy}^\dagger & 0 & 0 & 0\\
     0 & T & 0 & 0 & 0 & -P_{xy} & P_z & P_{xy}^\dagger\\
     -P_{xy}^\dagger & 0 & F & -H^\dagger & -K^\dagger & 0 & 0 & 0\\
     P_z & 0 & -H & \lambda & -H^\dagger & \Delta & 0 & 0\\
     P_{xy} & 0 & -K & H & G & 0 & \Delta & 0\\
     0 & -P_{xy}^\dagger & 0 & \Delta & 0 & G & -H^\dagger & -K^\dagger\\
     0 & P_z & 0 & 0 & \Delta & -H & \lambda & H^\dagger\\
     0 & P_{xy} & 0 & 0 & 0 & K & H & F
     \end{array}\right] \f]
     
     with:
     \f[ T = E_c + \frac{\hbar^2 k_z^2}{2m_\parallel^e} + \frac{\hbar^2(k_x^2 + k_y^2)}{2m_\perp^e} + \alpha_\parallel\varepsilon_{zz} + \alpha_\perp(\varepsilon_{xx} + \varepsilon_{yy})\f]
     \f[ F = \Delta_1 + \Delta_2 + \lambda + \theta\f]
     \f[ G = \Delta_1 - \Delta_2 + \lambda + \theta\f]
     \f[ \Delta = \sqrt{2}\Delta_2\f]
     \f[ \lambda = \tilde{A}_1 k_z^2 +\tilde{A}_2 (k_x^2 + k_y^2) + D_1\varepsilon_{zz} + D_2(\varepsilon_{xx} + \varepsilon_{yy})\f]
     \f[ \theta = \tilde{A}_3 k_z^2 + \tilde{A}_4 (k_x^2 + k_y^2) + D_3\varepsilon_{zz} + D_4(\varepsilon_{xx} + \varepsilon_{yy})\f]
     \f[ H = i(A_6k_z k_+ + A_7k_+ + D_6\varepsilon_{z+})\f]
     \f[ K = A_5 k_+^2 + D_5\varepsilon_+\f]
     \f[ P_{xy} = -\frac{i}{\sqrt{2}} \left(k_x + i\cdot k_y\right) \cdot P_p\f]
     \f[ P_{z} = -\frac{i}{\sqrt{2}} \left(k_z) \cdot P_o\right)\f]
     with:
     \f[ \varepsilon_+ = \varepsilon_{xx} - \varepsilon_{yy} +2i\varepsilon_{xy} \f]
     \f[ \varepsilon_{z+} = \varepsilon_{xz} + i \varepsilon_{yz} \f]
     \f[ k_+ = k_x + i k_y \f]
     \f[ k_\perp^2 = k_x^2 + k_y^2 \f]
   i*/
   PsiR ham8x8_WZ  (const PsiG &);

   
    /** \brief calculation of \f$ \hat{H} |\psi \rangle \f$i

       Calculates the 8 elements vector \f$ \hat{H} | \psi \rangle \f$. The spin is used here
       for creating a 8 band wave function.
   */    
  

/*   virtual PsiG preconditioner (const PsiG &psi) const {return one;};
   virtual SxDiracVec<TPrecCoeffG::TReal>
       ariasPreconditioner (const PsiG &psi) const;*/

   virtual SxDiracVec<TPrecCoeffG::TReal> preconditioner (const PsiG &psi) const {
      return preconditioner (psi, Payne);
   }

   SxDiracVec<TPrecCoeffG::TReal> 
      preconditioner (const PsiG &, 
            Preconditioner type) const;
   SxVector<Double> eKinS;
   SxVector3<Int> posBS;
   void calculateStrain(
                SxArray<SxString> &,
                SxVector<Double> &,
                SxVector<Double> &,
                SxVector<Double> &,
                SxVector<Double> &,
                SxVector<Double> &,
                SxVector<Double> &,
                double,
                bool);
};

namespace Timer {
   enum KP8x8HamTimer {
      secondDer,
      firstDer,
      hamiltonian,
      hamElements,
      hamSetup,
      op
   };
}

SX_REGISTER_TIMERS (Timer::KP8x8HamTimer)
{
   using namespace Timer;
   regTimer (secondDer,  "second derivative");
   regTimer (firstDer,   "first derivative");
   regTimer (hamiltonian,"Hamiltonian routine");
   regTimer (hamElements,"Hamiltonian elements");
   regTimer (hamSetup,   "Hamiltonian setup");
   regTimer (op,         "operator*");
}

#endif /* _SX_KDOTP8X8_H_ */
