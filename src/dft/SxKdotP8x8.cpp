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


/* Implementation of 8 band k dot p perturbation theory.
The Hamiltonian is taken from:

ref1:
Craig Pryor: Eight-band calculations of strained InAs/GaAs quantum dots
compared with one-, four- and six-band approximations,
Phys. Rev. \textbf{B}, Vol. 57, 12 (1998)
*/

#include <SxKdotP8x8.h>
#include <SxProjector.h>
#include <SxGBasis.h>

SxKdotP8x8::SxKdotP8x8 (const SxPWSet &waves,
                        const RhoR &rhoRIn)
   : SxHamiltonian (),
     SxRho (rhoRIn),
     wavesPtr (&waves)
// SxKdotP8x8
{
   SX_CHECK (rBasisPtr);
   s32 = sqrt(1.5);
   s2  = sqrt(2.);
   s3  = sqrt(3.);
   s6  = sqrt(6.);
   sdC.resize(1);
   fdC.resize(1);
   
   nStates = (short)waves.getNStates();
   step = 0;

   gBasisPtr = &(rBasisPtr->getGBasis ());
   
   SX_CHECK(rBasisPtr->getNComp() == 8);
   SX_CHECK(gBasisPtr->getNComp() == 8);
   
   gVecs.resize(1);
   D.resize(1);
   
}

SxArray<PsiR> SxKdotP8x8::firstDerivative (const PsiG &psiG)
{
//   SX_CLOCK(Timer::derivative)
   const SxRBasis &R = *rBasisPtr;
   SxArray<PsiR> gradR(3);
   int i;

   for (i=0; i < 3; i++)  {
      gradR(i) = (R | ( (-gVecs(i)) * I * psiG ) );
   }
   return gradR;
}

SxArray<SxArray<PsiR> > SxKdotP8x8::secondDerivative (const PsiG &psiG)
{
//   SX_CLOCK(Timer::derivative);
   const SxRBasis &R = *rBasisPtr;
   SxArray<SxArray<PsiR> > grad2R(3);
   int i;

	   
   for (i = 0; i < 3; i++)  {
      grad2R(i).resize(3);
      grad2R(i)(i) =
         - ( R | ((  gVecs(i) * gVecs(i)) * psiG) );
   }
      
   grad2R(0)(1) =
      - ( R | ((  gVecs(0) * (gVecs(1)) * psiG) ) );
   grad2R(0)(2) =
      - ( R | ((  gVecs(0) * (gVecs(2)) * psiG) ) );
   grad2R(1)(2) =
      - ( R | ((  gVecs(1) * (gVecs(2)) * psiG) ) );
   
   grad2R(1)(0) = grad2R(0)(1);
   grad2R(2)(0) = grad2R(0)(2);
   grad2R(2)(1) = grad2R(1)(2);
   return grad2R;
}

double SxKdotP8x8::interpolate (double y1, double y2, double comp, double bowing)
{
   return ( y2 - comp * (y2 - y1) - y2 * bowing * comp * (1- comp));
}

double SxKdotP8x8::mixing (double x1, double y1, double x2, double y2, double x)
{
   double m;
   double y2I = interpolate(y2, y1, p, bow);
   m = (y1 + (x+x2) * (y2I - y1) / (x2 - x1));
   return m;
}

PrecEnergy SxKdotP8x8::getETrial()
{
   return eTrial;
}

void SxKdotP8x8::update (const SxFermi &fermi)
{
//   SX_CLOCK (Timer::hUpdate);
   fermiPtr = &fermi;
   // --- external potential
   eExt = 0.;
}

void SxKdotP8x8::compute (const SxFermi &fermi,
                       bool tauChanged,
                       bool rhoChanged)
{
   if (tauChanged && rhoChanged) {/* just to get rid of compiler-warning*/};
   update(fermi);
}

PsiG SxKdotP8x8::operator* (const PsiG &psiG)
{
   SX_CHECK (psiG.getBasisPtr());

   const SxGBasis &G = *dynamic_cast<const SxGBasis *>(psiG.getBasisPtr());

   // if not yet done: initialize gVecs for first and second derivatives
   // NEEDS A BETTER SOLUTION (?)
   if (gVecs.getSize() == 1)  {
      int i, iComp;
      g2.resize(psiG.getSize());
      g2.setBasis(psiG.getBasisPtr());
      for (iComp = 0; iComp < 8; iComp++)
         g2.setComponent(iComp, G.g2);
      gVecs.resize(3);
      for (i = 0; i < 3; i++)  {
         gVecs(i).resize(psiG.getSize());
         gVecs(i).setBasis(psiG.getBasisPtr());
         for (iComp = 0; iComp < 8; iComp++)
            gVecs(i).setComponent(iComp, G.gVec.colRef(i));
      }
   }
   
   // return gradient: ham8x8_ZB and _WZ is the realspace HPsi
   step++;
   if (Hamiltonian == 0) // use zincblende hamiltonian
      return (G | ham8x8_ZB (psiG));
   else //if (Hamiltonian == 1) // use wurtzite hamiltonian
      return (G | ham8x8_WZ (psiG));
}

// the hamiltonian itself is implemented here

PsiR SxKdotP8x8::ham8x8_WZ (const PsiG &psiG)
{
   SX_CLOCK(Timer::hamiltonian);
   SX_CHECK (psiG.getBasisPtr());
   const SxRBasis &R = *rBasisPtr;
   if (D.getSize() == 1)  { // setting up kinetic operator for preconditioner
      cout << "preparing preconditioner...";

      SxArray<SxArray<PsiR> > sd;
      sd.resize(3);
      for (int i = 0; i < 3; i++)  {
         sd(i).resize(3);
         for (int j = 0; j < 3; j++)
            sd(i)(j) = gVecs(i) * gVecs(j);
      }
      SxArray<PsiG> T(8), lambda(8), F(8);
      for (int iComp = 0; iComp < 8; iComp++)  {
         T(iComp) = aP2 * (sd(0)(0) + sd(1)(1)).getComponent(iComp)
            + aP1 * sd(2)(2).getComponent(iComp);
         lambda(iComp) = -aP(1) * (sd(0)(0) + sd(1)(1)).getComponent(iComp)
            - aP(0) * sd(2)(2).getComponent(iComp);
         F(iComp) = -(aP(1) + aP(3)) * (sd(0)(0) + sd(1)(1)).getComponent(iComp)
            - (aP(0) + aP(2)) * sd(2)(2).getComponent(iComp);
      }
 
      D.resize(psiG.getSize());
      D = 0. * psiG;
      D.setComponent(0, T(0));
      D.setComponent(1, T(1));
      D.setComponent(2, F(2));
      D.setComponent(3, lambda(3));
      D.setComponent(4, F(4));
      D.setComponent(5, F(5));
      D.setComponent(6, lambda(6));
      D.setComponent(7, F(7));
      cout << " done." << endl;
   }


   PsiR H = (R | psiG);

   SxArray<PsiR> fd; // first derivatives dPsiG/d(x,y,z)

   SxArray<SxArray<PsiR> > sd; // second derivatives d/d(x,y,z)(dPsiG/d(x,y,z)
   {
   SX_CLOCK(Timer::firstDer);
   fd = firstDerivative(psiG);
   }
   {
   SX_CLOCK(Timer::secondDer);
   sd = secondDerivative(psiG);
   }

   // --- for given k-point
   if (kxSet)  {
      sd(0)(0) = -kx * kx * (R | psiG);
      sd(0)(1) = -kx * fd(1);
      sd(1)(0) = -kx * fd(1);
      sd(2)(0) = -kx * fd(2);
      sd(0)(2) = -kx * fd(2);
      fd(0) = -I*kx * (R | psiG);
   }
   if (kySet)  {
      sd(1)(1) = -ky * ky * (R | psiG);
      sd(0)(1) = -ky * fd(0);
      sd(1)(0) = -ky * fd(0);
      sd(2)(1) = -ky * fd(2);
      sd(1)(2) = -ky * fd(2);
      fd(1) = -I*ky * (R | psiG);
   }
   if (kzSet)  {
      sd(2)(2) = -kz * kz * (R | psiG);
      sd(2)(1) = -kz * fd(1);
      sd(1)(2) = -kz * fd(1);
      sd(2)(0) = -kz * fd(0);
      sd(0)(2) = -kz * fd(0);
      fd(2) = -I*kz * (R | psiG);
   }
   if (kxSet && kySet)  {
      sd(0)(1) = -kx * ky * (R | psiG);
      sd(1)(0) = -kx * ky * (R | psiG);
   }
   if (kxSet && kzSet)  {
      sd(0)(2) = -kx * kz * (R | psiG);
      sd(2)(0) = -kx * kz * (R | psiG);
   }
   if (kySet && kzSet)  {
      sd(1)(2) = -ky * kz * (R | psiG);
      sd(2)(1) = -ky * kz * (R | psiG);
   }
   
    




   int iComp;

   SxArray<PsiR> T(2), F(2), M(2), lambda(6), theta(6), // diagonal
      K(6), K_conj(6), L(6), L_conj(6), // off-diagonal 6 band
      Pp(8), Pm(8), Pz(8), // cb-vb coupling
      Delta(6); // deltaSO/CR

   {
      SX_CLOCK(Timer::hamElements);
 
   for (iComp = 0; iComp < 8; iComp++)  {
      psiR <<= H.getComponent(iComp);
      for (int i = 0; i < 3; i++)  {
         fdC(i) <<= fd(i).getComponent(iComp);
         for (int j = 0; j < 3; j++)
            sdC(i)(j) <<= sd(i)(j).getComponent(iComp);
      }

      if (iComp < 2) {  // conduction band
            T(iComp) =
                -A2 * (sdC(0)(0)
                   + sdC(1)(1))
               - A1 * sdC(2)(2)
               + (eCBO  + externalField)* psiR;
        
      }
      else {  // valence band
         lambda(iComp - 2) = .5 * ( A(0) * sdC(2)(2)
               + A(1) * (sdC(0)(0)
                  + sdC(1)(1)) )
               + (eVBOMod(iComp - 2) + externalField) * psiR;
 
            theta(iComp - 2) = .5 * ( A(2) * sdC(2)(2)
                  + A(3) * (sdC(0)(0)
                  + sdC(1)(1)) );

            if (iComp == 2)
               F(0) = lambda(iComp - 2) + theta(iComp - 2);
            if (iComp == 4)
               M(0) = lambda(iComp - 2) + theta(iComp - 2);
            if (iComp == 7)
               F(1) = lambda(iComp - 2) + theta(iComp - 2);
            if (iComp == 5)
               M(1) = lambda(iComp - 2) + theta(iComp - 2);


            // --- off-diagonal vb Hamiltonian elements
            K(iComp - 2) = .5 * A(4) * (sdC(0)(0)
                  + 2. * I * sdC(0)(1)
                  - sdC(1)(1));
            K_conj(iComp - 2) = .5 * A(4) * (sdC(0)(0)
                  - 2. * I * sdC(0)(1)
                  - sdC(1)(1));
            
            L(iComp - 2) = .5 * A(5) * (sdC(0)(2)
                  + I * sdC(1)(2));
            L_conj(iComp - 2) = .5 * A(5) * (sdC(0)(2)
                  - I * sdC(1)(2));
            
            Delta(iComp - 2) = s2 * delta(2) * psiR;


      }

      // --- cb <-> vb interaction elements
      Pp(iComp) = -I * (1./s2) * (fdC(0)
            + I * fdC(1)) * pP;
      Pm(iComp) = -I * (1./s2) * (fdC(0)
            - I * fdC(1)) * pP;
      Pz(iComp) = -I * fdC(2) * pO;
 
   }
 
   }
   /* basis notation:
      sU, sD, u1, u3, u2, u5, u6, u4
    */
   {
      SX_CLOCK(Timer::hamSetup);
   H.setComponent(0,T(0) - Pp(2) + Pz(3) + Pm(4));
   H.setComponent(1,T(1) - Pp(5) + Pz(6) + Pm(7));
   H.setComponent(2,-Pm(0) + F(0) - L_conj(1) - K_conj(2));
   H.setComponent(3,Pz(0) - L(0) + lambda(1) + L_conj(2) + Delta(3));
   H.setComponent(4,Pp(0) - K(0) + L(1) + M(0) + Delta(4));
   H.setComponent(5,-Pm(1) + Delta(1) + M(1) - L_conj(4) - K_conj(5));
   H.setComponent(6,Pz(1) + Delta(2) - L(3) + lambda(4) + L_conj(5));
   H.setComponent(7,Pp(1) - K(3) + L(4) + F(1));

   }
   
   return H;
 

   /*
   // --- debugging part: show bandstructure
   if (bs)  {
      const SxRBasis &R = *rBasis;
      SxFFT3d fftMesh;
      fftMesh.mesh = R.getMesh();
      int fftIdx;// = fftMesh.getFFTIdx(posBS); TODO!!!
      SxComplex16 t, f, g, lambda, theta, h, k, pp, pm, pz;
      double ev, ec, a1, a2;
      SxArray<Real8> a(7), d(3);
      int i, col, row;
      for (i = 0; i < 7; i++)
         a(i) = A(i)(fftIdx);

      a1 = A1(fftIdx);
      a2 = A2(fftIdx);
      ev = eVBO(fftIdx);
      ec = eCBO(fftIdx);
      for (i = 0; i < 3; i++)
      d(i)  = delta(i)(fftIdx);

      cout << "a1: " << a1 << endl;
      cout << "a2: " << a2 << endl;
      cout << "ev: " << ev << endl;
      cout << "ec: " << ec << endl;
      for (i = 0; i < 3; i++)
         cout << "delta(" << i << "): " << d(i) << endl;

      for (i = 0; i < 7; i++)
         cout << "A(" << i << "): " << a(i) << endl;

      double kx,kx2,ky,ky2,kz,kz2;
      int idx;

      cout << "#bandstructure:" << endl;
      for (idx = -10; idx <  10; idx++)  {*/
/*         if (idx < 0)  {// L to gamma
            kx = 0.05 * idx;
            ky = 0.05 * idx;
            kz = 0.05 * idx;
            kx2 = kx * kx;
            ky2 = ky * ky;
            kz2 = kz * kz;
         }
         else  {
            kx = 0.1 * idx;
            ky = 0.;
            kz = 0.;
            kx2 = kx * kx;
            ky2 = 0.;
            kz2 = 0.;
         }*/
   
   /*
         kx = 0.; kx2 = 0.;
         ky = 0.; ky2 = 0.;
         kz = 0.; kz2 = 0.;
         
         // diagonal terms
         
         // cb
         t = a2 * (kx2 + ky2) + a1 * kz2 + ec;
         // vb
         lambda = -.5 * a(0) * (kx2 + ky2) - .5 * a(1) * kz2 + ev;
         theta = -.5 * a(2) * (kx2 + ky2) - .5 * a(3) * kz2;
         f = lambda + theta;
         g = lambda + theta;
         
         // off-diagonal vb terms
         k = .5 * a(4) * (kx2 + 2. * I * kx * ky - ky2);
         h = .5 * a(5) * (kx * kz + I * ky * kz);
         pp = (1./s2) * pP(0) * (kx + I * ky);    
         pm = (1./s2) * pP(0) * (kx - I * ky);    
         pz = pO(0) * kz;    
         SxMatrix<Complex16> Mat(8,8);
         Mat(0,0) = t;
         Mat(0,1) = 0.;
         Mat(0,2) = -pp;
         Mat(0,3) = pz;
         Mat(0,4) = pm;
         Mat(0,5) = 0.;
         Mat(0,6) = 0.;
         Mat(0,7) = 0.;
         
         Mat(1,0) = 0.;
         Mat(1,1) = t;
         Mat(1,2) = 0.;
         Mat(1,3) = 0.;
         Mat(1,4) = 0.;
         Mat(1,5) = -pp;
         Mat(1,6) = pz;
         Mat(1,7) = pm;
         
         Mat(2,0) = -pm;
         Mat(2,1) = 0.;
         Mat(2,2) = f;
         Mat(2,3) = - h.conj();
         Mat(2,4) = - k.conj();
         Mat(2,5) = 0.;  // = 0 !
         Mat(2,6) = 0.;
         Mat(2,7) = 0.;
         
         Mat(3,0) = pz;
         Mat(3,1) = 0.;
         Mat(3,2) = -h;
         Mat(3,3) = lambda;
         Mat(3,4) = h.conj(); 
         Mat(3,5) = 0.;
         Mat(3,6) = 0.;
         Mat(3,7) = 0.;
         
         Mat(4,0) = pp;
         Mat(4,1) = 0.;
         Mat(4,2) = -k;
         Mat(4,3) = h;  
         Mat(4,4) = g;
         Mat(4,5) = 0.;
         Mat(4,6) = 0.;
         Mat(4,7) = 0.;
         
         Mat(5,0) = 0.;
         Mat(5,1) = -pm;
         Mat(5,2) = 0.;
         Mat(5,3) = 0.;
         Mat(5,4) = 0.;
         Mat(5,5) = g;
         Mat(5,6) = -h.conj();
         Mat(5,7) = -k.conj();
         
         Mat(6,0) = 0.;
         Mat(6,1) = pz;
         Mat(6,2) = 0.;
         Mat(6,3) = 0.;
         Mat(6,4) = 0.;
         Mat(6,5) = -h;
         Mat(6,6) = lambda;
         Mat(6,7) = h.conj();
         
         Mat(7,0) = 0.;
         Mat(7,1) = pp;
         Mat(7,2) = 0.;
         Mat(7,3) = 0.;
         Mat(7,4) = 0.;
         Mat(7,5) = -k;
         Mat(7,6) = h;
         Mat(7,7) = f;

	 if (strained)  {
	    for (col = 0; col < nSpin; col++)
	       for (row = 0; row < nSpin; row++)  {
	          Mat(col,row) += (H_Strain(col)(row))(fftIdx)
                                + externalPotential(fftIdx);
                  }
         }
         cout << kx << " " << (Mat.eigenvalues())(0).re * EV
            << " \t" <<  (Mat.eigenvalues())(1).re * EV;
         cout << " \t" << (Mat.eigenvalues())(2).re * EV
            << " \t" <<  (Mat.eigenvalues())(3).re * EV;
         cout << " \t" << (Mat.eigenvalues())(4).re * EV
            << " \t" <<  (Mat.eigenvalues())(5).re * EV;
         cout << " \t" << (Mat.eigenvalues())(6).re * EV
            << " \t" <<  (Mat.eigenvalues())(7).re * EV << endl;
      }
   SX_QUIT; 
   }
   
   int col, row;
   
   for (col = 0; col < nSpin; col++)  {
      for (row = 0; row < nSpin; row++) {
         if (strained)  {
            PsiR psiR;// = allPsiG(row).toR(); TODO!!!
              H_WZ(col) += H_Strain(col)(row) * psiR;
         }

      }


   }
*/
}



PsiR SxKdotP8x8::ham8x8_ZB (const PsiG &psiG)
{
   SX_CHECK (psiG.getBasisPtr());
   const SxRBasis &R = *rBasisPtr;

   if (D.getSize() == 1)  { // setting up kinetic operator for preconditioner
      cout << "preparing preconditioner...";
      SxArray<SxArray<PsiR> > sd;
      sd.resize(3);
      for (int i = 0; i < 3; i++)  {
         sd(i).resize(3);
         for (int j = 0; j < 3; j++)
            sd(i)(j) = gVecs(i) * gVecs(j);
      }
      SxArray<PsiG> T(8), P(8), Q(8);
      for (int iComp = 0; iComp < 8; iComp++)  {
         T(iComp) = 0.5 * rMassP * 
            (sd(0)(0) + sd(1)(1) + sd(2)(2)).getComponent(iComp);
         P(iComp) = -0.5 * gammaP(0) *
            (sd(0)(0) + sd(1)(1) + sd(2)(2)).getComponent(iComp);
         Q(iComp) = -0.5 * gammaP(1) * 
            (sd(0)(0) + sd(1)(1) - 2. * sd(2)(2)).getComponent(iComp);
      }
 
      D.resize(psiG.getSize());
      D = 0. * psiG;
      D.setComponent(0, T(0));
      D.setComponent(1, T(1));
      D.setComponent(2, -P(2) + Q(2));
      D.setComponent(3, -P(3) - Q(3));
      D.setComponent(4, -P(4) - Q(4));
      D.setComponent(5, -P(5) + Q(5));
      D.setComponent(6, -P(6));
      D.setComponent(7, -P(7));
      cout << " done." << endl;
   }

   PsiR H = (R | psiG);

   SxArray<PsiR> fd; // first derivatives dPsiG/d(x,y,z)

   SxArray<SxArray<PsiR> > sd; // second derivatives d/d(x,y,z)(dPsiG/d(x,y,z)
   fd = firstDerivative(psiG);
   sd = secondDerivative(psiG);

   // N = R in reference for Hamiltonian!
   SxArray<PsiR> T(2), U(8), V(8), V_conj(8), P(6), Q(6), N(6), S(6),
                 N_conj(6), S_conj(6);
 
   // --- calculating matrix elements for all six states
  
   int iComp; 
   for (iComp = 0; iComp < 8; iComp++)  {
      psiR <<= H.getComponent(iComp);
      for (int i = 0; i < 3; i++)  {
         fdC(i) <<= fd(i).getComponent(iComp);
         for (int j = 0; j < 3; j++)
            sdC(i)(j) <<= sd(i)(j).getComponent(iComp);
      }


      if (iComp < 2) {
         T(iComp) = -0.5 * reciprocalEffMass 
            * (sdC(0)(0) + sdC(1)(1) + sdC(2)(2))
                      + psiR * (eCBO + externalField);
      }
      else  {
         P(iComp - 2) = - (eVBO + externalField) * psiR - gamma(0) * 0.5 
            * (sdC(0)(0) + sdC(1)(1) + sdC(2)(2));
         Q(iComp - 2) = -gamma(1) * 0.5
            * (sdC(0)(0) + sdC(1)(1) - 2. * sdC(2)(2));
         
         if (iComp > 5)
            P(iComp - 2) += deltaSpinOrbital * psiR;
      
         // determine R
         if ((iComp == 4 ) || (iComp == 5) || (iComp == 6))
            N(iComp - 2) = s3 * 0.5
               * (gamma(1) * (sdC(0)(0) - sdC(1)(1))
                     - 2. * I * gamma(2) * sdC(1)(0));
      
         if ((iComp == 2 ) || (iComp == 3) || (iComp == 7))
            N_conj(iComp - 2) = s3 * 0.5
               * (gamma(1) * (sdC(0)(0) - sdC(1)(1))
                     + 2. * I * gamma(2) * sdC(1)(0));
      
         // determine S
         if (iComp != 3 )
            S(iComp - 2) = -s3 * gamma(2)
               * (sdC(0)(2) - I * sdC(1)(2));
      
         if (iComp != 4 )
            S_conj(iComp - 2) = -s3 * gamma(2)
               * (sdC(0)(2) + I * sdC(1)(2));
      }
      // determine U
      U(iComp) = ( -I / s3 ) * fdC(2) * p0;
     
      // determine V
      V(iComp)   = -((I / s6)) * (fdC(0) - I * fdC(1)) * p0;

      V_conj(iComp) = -((I / s6)) * (fdC(0) + I * fdC(1)) * p0;
     
   }
   H.setComponent(0, T(0)
         + V_conj(2) + s3 * V(4) - s2 * U(5) - U(6) + s2 * V_conj(7));
   H.setComponent(1, T(1)
         - s2 * U(2) - s3 * V_conj(3) - V(5) + s2 * V(6) + U(7));
   H.setComponent(2, V(0) - s2 * U(1)
         - P(0) + Q(0) - S_conj(1) + N(2) + s32 * S(4) - s2 * Q(5));
   H.setComponent(3, -s2 * V(1)
         - S(0) - P(1) - Q(1) + N(3) - s2 * N(4) + 1./s2 * S(5));
   H.setComponent(4, s3 * V_conj(0)
         + N_conj(0) - P(2) - Q(2) + S_conj(3) + 1./s2 * S_conj(4)
         + s2 * N_conj(5) );
   H.setComponent(5, -s2 * U(0) - V_conj(1)
         + N_conj(1) + S(2) - P(3) + Q(3) + s2 * Q(4) + s32 * S_conj(5));
   H.setComponent(6, -U(0) + s2 * V_conj(1)
         + s32 * S_conj(0) - s2 * N_conj(1) + 1./s2 * S(2)
         + s2 * Q(3) -P(4) );
   H.setComponent(7, s2 * V(0) + U(1)
         - s2 * Q(0) + 1./s2 * S_conj(1) + s2 * N(2) + s32 * S(3) -P(5));
/*
   // --- debugging part: show bandstructure
   if (bs)  {
      const SxRBasis &R = *rBasis;
      SxFFT3d fftMesh;
      fftMesh.mesh = R.getMesh();
      int fftIdx;// = fftMesh.getFFTIdx(posBS); TODO!!!
      SxComplex16 a,p,q,r,s,u,v,z;
      double ev, ec, l1, l2, l3, m, d, po;

      m = (reciprocalEffMass(fftIdx));
      ev = eVBO(fftIdx);
      ec = eCBO(fftIdx);
      l1 = gamma(0)(fftIdx);
      l2 = gamma(1)(fftIdx);
      l3 = gamma(2)(fftIdx);
      d  = deltaSpinOrbital(fftIdx);
      po = p0(fftIdx);

      cout << "m: " << m << endl;
      cout << "ev: " << ev << endl;
      cout << "ec: " << ec << endl;
      cout << "delta: " << d << endl;
      cout << "l1: " << l1 << endl;
      cout << "l2: " << l2 << endl;
      cout << "l3: " << l3 << endl;
      cout << "p0: " << po << endl;
      cout << "v: " << (-1./sqrt(6.))*po << endl;

      double kx,kx2,ky,ky2,kz,kz2;
      int idx;

      cout << "#bandstructure:" << endl;
      for (idx = -100; idx <  100; idx++)  {
         if (idx < 0)  { // L = .5,.5,.5 -> gamma = 0.,0.,0.
            kx = fabs(0.005 * idx);
            ky = fabs(0.005 * idx);
            kz = fabs(0.005 * idx);
            kx2 = kx * kx;
            ky2 = ky * ky;
            kz2 = kz * kz;
         } else  { // gamma -> L = 0.,.5,0.
            kx = 0.;
            ky = 0.01 * idx;
            kz = 0.;
            kx2 = 0.;
            ky2 = ky * ky;
            kz2 = 0.;
         }

         a = ec + (.5*m) * (kx2 + ky2 + kz2);
         p = -ev + .5 * l1 * (kx2 + ky2 + kz2);
         q = .5 * l2 * (kx2 + ky2 - 2. * kz2);
         z = ev - d -.5 * l1 * (kx2 + ky2 + kz2);
         r = -(sqrt(3.)/2.) * (l2 * (kx2 - ky2) - I * 2. * l3 * kx * ky);
         s = sqrt(3.) * l3 * kz* ( kx - I * ky);
         u = (1./sqrt(3.)) * po * kz;
         v = (1./sqrt(6.)) * po * (kx - I * ky);
         
         SxMatrix<Complex16> Mat(8,8);
         Mat(0,0) = a;
         Mat(0,1) = 0.;
         Mat(0,2) = v.conj();
         Mat(0,3) = 0.;
         Mat(0,4) = s3 * v;
         Mat(0,5) = -s2 * u;
         Mat(0,6) = -u;
         Mat(0,7) = s2 * v.conj();
         
         Mat(1,0) = 0.;
         Mat(1,1) = a;
         Mat(1,2) = -s2 * u;
         Mat(1,3) = -s3 * v.conj();
         Mat(1,4) = 0.;
         Mat(1,5) = -v;
         Mat(1,6) = s2 * v;
         Mat(1,7) = u;
         
         Mat(2,0) = v;
         Mat(2,1) = -s2 * u;
         Mat(2,2) = -p + q;
         Mat(2,3) = -s.conj();
         Mat(2,4) = r;
         Mat(2,5) = 0.;  // = 0 !
         Mat(2,6) = s32 * s;
         Mat(2,7) = - s2 * q;
         
         Mat(3,0) = 0.;
         Mat(3,1) = -s3 * v;
         Mat(3,2) = - s;
         Mat(3,3) = - p - q;
         Mat(3,4) = 0.; 
         Mat(3,5) = r;
         Mat(3,6) = - s2 * r;
         Mat(3,7) = (1./s2) * s;
         
         Mat(4,0) = s3 * v.conj();
         Mat(4,1) = 0.;
         Mat(4,2) = r.conj();
         Mat(4,3) = 0.;  
         Mat(4,4) = - p - q;
         Mat(4,5) = s.conj();
         Mat(4,6) = (1./s2) * s.conj();
         Mat(4,7) = s2 * r.conj();
         
         Mat(5,0) = -s2 * u;
         Mat(5,1) = -v.conj();
         Mat(5,2) = 0.;
         Mat(5,3) = r.conj();
         Mat(5,4) = s;
         Mat(5,5) = - p + q;
         Mat(5,6) = s2 * q;
         Mat(5,7) = s32 * s.conj();
         
         Mat(6,0) = -u;
         Mat(6,1) = s2 * v.conj();
         Mat(6,2) = s32 * s.conj();
         Mat(6,3) = - s2 * r.conj();
         Mat(6,4) = (1./s2) * s;
         Mat(6,5) = s2 * q;
         Mat(6,6) = z;
         Mat(6,7) = 0.;
         
         Mat(7,0) = s2 * v;
         Mat(7,1) = u;
         Mat(7,2) = -s2 * q;
         Mat(7,3) = (1./s2) * s.conj();
         Mat(7,4) = s2 * r;
         Mat(7,5) = s32 * s;
         Mat(7,6) = 0.;
         Mat(7,7) = z;
      
         cout << idx << " " << (Mat.eigenvalues())(0).re * EV
            << " \t" <<  (Mat.eigenvalues())(1).re * EV;
         cout << " \t" << (Mat.eigenvalues())(2).re * EV
            << " \t" <<  (Mat.eigenvalues())(3).re * EV;
         cout << " \t" << (Mat.eigenvalues())(4).re * EV
            << " \t" <<  (Mat.eigenvalues())(5).re * EV;
         cout << " \t" << (Mat.eigenvalues())(6).re * EV
            << " \t" <<  (Mat.eigenvalues())(7).re * EV << endl;
      }
   SX_QUIT; 
   }
   
   int col, row;
   
   for (col = 0; col < nSpin; col++)  {
      for (row = 0; row < nSpin; row++) {
         if (strained)  {
            PsiR psiR;// = allPsiG(row).toR(); TODO!!!
              H(col) += H_Strain(col)(row) * psiR;
         }

      }

   }
*/
   
   return H;
}

   // --- determine H Psi in k-Space

const SxSymbolTable *
SxKdotP8x8::getHamiltonianGroup (const SxSymbolTable *table)
{
   SX_CHECK (table);
   SxSymbolTable *hGroup = NULL;
   try  { hGroup = table->getGroup("KP8x8Hamiltonian"); }
   catch (const SxException&) { /*EXIT;*/ }
   return hGroup;
}

int SxKdotP8x8::getNEmptyStates () {
   return nEmptyStates;
}

void SxKdotP8x8::validateNStates (int nStatesIn)
{
   int nStatesMax = wavesPtr->getGkBasis().nGkMin;

   if (nStatesIn > nStatesMax)  {
      sxprintf ("Error: input parameter(s) \"nEmptyStates\" or "
              "\"nEmptyStatesChi\" in (PW)Hamil-\n"
              "       tonian group chosen too large. Please reduce it/them or "
              "try with a\n"
              "       larger cutoff \"eCut\". The number of states you want to "
              "compute must\n"
              "       not exceed %d.\n", nStatesMax);
      SX_QUIT;
   }
}

void SxKdotP8x8::read (const SxSymbolTable *table)
{
   gamma.resize(3);
   delta.resize(3);
   Hamiltonian = 0;
   double evHi, evLo, ecLo, ecHi, deltaSOHi, deltaSOLo, eMassHi, eMassLo;
   double mPHi, mPLo, mOHi, mOLo, egHi, egLo, eP1Hi, eP1Lo, eP2Hi, eP2Lo;
   double pOHi, pOLo, pPHi, pPLo; // P orthogonal and parallel in WZ 
   double factor; // prefactor in strain field
   deltaSOHi = deltaSOLo = eMassHi = eMassLo = mPHi = mPLo = mOHi = mOLo = 0.;
   eP1Hi = eP1Lo = eP2Hi = eP2Lo = pOHi = pOLo = pPHi = pPLo = factor = 0.; 
   SxVector<Double> AHi, ALo, DHi, DLo, deltaHi, deltaLo, aHi, aLo;
   SxVector<Double> ZBParHi, ZBParLo; // zincblend strain parameter set
   SxArray<SxString> fileNames; // files providing strain field
   AHi.resize(7);
   A.resize(7);
   eVBOMod.resize(6);
   ALo.resize(7);
   DHi.resize(6);
   DLo.resize(6);
   aHi.resize(2);
   aLo.resize(2);
   double ePHi, ePLo;
   ePHi = ePLo = 0.;
   eTrial = 1.e50;
   int i, j, idx;
   electrons = false;
   strained = false;
   bs = false;
   p = 1.;
   bow = 0.;
   double maximum, minimum;
   SxVector3<Double> luttingerHi, luttingerLo, modLutHi, modLutLo;
   SxMeshR mix;
   SxArray<SxMeshR> vExtR;
   cout << "reading K*P-input values\n";
   SxSymbolTable *hamiltonian =  table->getGroup("KP8x8Hamiltonian");
   SxSymbolTable *vExt = hamiltonian->getGroup("vExt");
   SxString vExtFile   = vExt->get("file")->toString();
   SxBinIO ioStructure (vExtFile, SxBinIO::BINARY_READ_ONLY);
   SxMatrix3<Double> cell;
   SxVector3<Int>    dim;
   vExtR.resize (1);
   vExtR = ioStructure.readMesh (&cell, &dim);
   externalPotential = vExtR(0);
   ioStructure.close();
   int size = (int)externalPotential.getSize();
   cout << "real-space size: " << size << endl;
   { // initialization of component second derivative
      cout << "cRSize: " << size << endl;
      sdC.resize(3);
      fdC.resize(3);
      psiR.resize(size);
      for (i = 0; i < 3; i++)  {
         sdC(i).resize(3);
         fdC(i).resize(size);
         for (j = 0; j < 3; j++)
            sdC(i)(j).resize(size);
      }
      
   }

   minimum = - externalPotential.minval();
   maximum = - externalPotential.maxval();
   cout << "external potential from: " << minimum << " to " << maximum << "\n";
   mix.resize(size);
   zero.resize(size);
   for (idx = 0; idx < size; idx++) {
      mix(idx) = mixing(minimum, 1., maximum, 0., externalPotential(idx));
      zero(idx) = 0.;
   }
   cout << "mixing: " << mix.minval() << " to " << mix.maxval() << endl;

   ekt          = hamiltonian->get("ekt")->toReal() / HA2EV;
   ecHi         = hamiltonian->get("cbOffsetHi")->toReal() / HA2EV;
   ecLo         = hamiltonian->get("cbOffsetLo")->toReal() / HA2EV;
   evHi         = hamiltonian->get("vbOffsetHi")->toReal() / HA2EV;
   evLo         = hamiltonian->get("vbOffsetLo")->toReal() / HA2EV;
   if (hamiltonian->contains("eTrial"))  
      eTrial       = hamiltonian->get("eTrial")->toReal() / HA2EV;
   if (hamiltonian->contains("composition"))
      p = hamiltonian->get("composition")->toReal ();
   if (hamiltonian->contains("bowing"))
      bow = hamiltonian->get("bowing")->toReal ();
   if (hamiltonian->contains("luttingerHi"))  {
      luttingerHi  =
         SxVector3<Double>(hamiltonian->get("luttingerHi")->toList());
      luttingerLo  =
         SxVector3<Double>(hamiltonian->get("luttingerLo")->toList());
      deltaSOHi    = hamiltonian->get("soHi")->toReal () / HA2EV;
      deltaSOLo    = hamiltonian->get("soLo")->toReal () / HA2EV;
      eMassHi      = hamiltonian->get("cbMassHi")->toReal ();
      eMassLo      = hamiltonian->get("cbMassLo")->toReal ();
      ePHi         = hamiltonian->get("ePHi")->toReal() / HA2EV;
      ePLo         = hamiltonian->get("ePLo")->toReal() / HA2EV;
      // calculate modified Luttinger parameter
      modLutHi(0) = luttingerHi(0) - (ePHi) / (3. * (ecHi - evHi) + deltaSOHi);
      modLutHi(1) = luttingerHi(1) - 0.5 * (ePHi) / (3. * (ecHi - evHi) + deltaSOHi);
      modLutHi(2) = luttingerHi(2) - 0.5 * (ePHi) / (3. * (ecHi - evHi) + deltaSOHi);
      modLutLo(0) = luttingerLo(0) - (ePLo) / (3. * (ecLo - evLo) + deltaSOLo);
      modLutLo(1) = luttingerLo(1) - 0.5 * (ePLo) / (3. * (ecLo - evLo) + deltaSOLo);
      modLutLo(2) = luttingerLo(2) - 0.5 * (ePLo) / (3. * (ecLo - evLo) + deltaSOLo);
      gammaP.resize(3);
      for (i = 0; i < 3; i++)
         gammaP(i) = luttingerLo(i);
      cout << "modified Luttinger parameter: \nHi: [ "
         << modLutHi(0) << ", " << modLutHi(1) << ", " << modLutHi(2) << "]\n"
         << "Lo: [ "
         << modLutLo(0) << ", " << modLutLo(1) << ", " << modLutLo(2) << "]"
         << endl;
         
      
      Hamiltonian = 0;
      if (hamiltonian->containsGroup("strain"))  {
         strained = true;
         ZBParHi.resize(4);
         ZBParLo.resize(4);
         fileNames.resize(6);
         SxSymbolTable *strain = hamiltonian->getGroup("strain");
         ZBParHi(0)   = strain->get("acHi")->toReal () / HA2EV;
         ZBParHi(1)   = strain->get("avHi")->toReal () / HA2EV;
         ZBParHi(2)   = strain->get("bHi")->toReal () / HA2EV;
         ZBParHi(3)   = strain->get("dHi")->toReal () / HA2EV;
         ZBParLo(0)   = strain->get("acLo")->toReal () / HA2EV;
         ZBParLo(1)   = strain->get("avLo")->toReal () / HA2EV;
         ZBParLo(2)   = strain->get("bLo")->toReal () / HA2EV;
         ZBParLo(3)   = strain->get("dLo")->toReal () / HA2EV;
         fileNames(0) = strain->get("sigmaXX")->toString();
         fileNames(1) = strain->get("sigmaXY")->toString();
         fileNames(2) = strain->get("sigmaXZ")->toString();
         fileNames(3) = strain->get("sigmaYY")->toString();
         fileNames(4) = strain->get("sigmaYZ")->toString();
         fileNames(5) = strain->get("sigmaZZ")->toString();
      }

      if (hamiltonian->containsGroup("externalField"))  {
         SxSymbolTable *extField = hamiltonian->getGroup("externalField");
         SxString fileName = (extField->get("file")->toString() );
         SxBinIO io (fileName, SxBinIO::BINARY_READ_ONLY);
         vExtR.resize (1);
         vExtR = io.readMesh (&cell, &dim);
         externalField = vExtR(0);
         io.close();
      }
      else 
         externalField = zero;
   }  else if (hamiltonian->contains("AHi"))  {
      AHi          = SxVector<Double>(hamiltonian->get("AHi")->toList ());
      ALo          = SxVector<Double>(hamiltonian->get("ALo")->toList ());
      mOHi         = hamiltonian->get("mOHi")-> toReal (); // effective electron
      mOLo         = hamiltonian->get("mOLo")-> toReal (); // mass orthogonal

      mPHi         = hamiltonian->get("mPHi")-> toReal (); // effective electron
      mPLo         = hamiltonian->get("mPLo")-> toReal (); // mass parallel
      if (hamiltonian->contains("eP1Hi"))  {
         eP1Hi         = hamiltonian->get("eP1Hi")-> toReal () / HA2EV;
         eP2Hi         = hamiltonian->get("eP2Hi")-> toReal () / HA2EV;
         eP1Lo         = hamiltonian->get("eP1Lo")-> toReal () / HA2EV;
         eP2Lo         = hamiltonian->get("eP2Lo")-> toReal () / HA2EV;
      }
      deltaHi      =
         SxVector3<Double>(hamiltonian->get("deltaHi")->toList ());
      deltaLo      =
         SxVector3<Double>(hamiltonian->get("deltaLo")->toList ());
      for (i = 0; i < 3; i++)  {
         deltaHi(i) *= -1./HA2EV;
         deltaLo(i) *= -1./HA2EV;
         cout << "Delta[" << i << "]: " << deltaHi(i) << " to " << deltaLo(i)
            << endl;
      }
      for (i = 1; i < 3; i++)  {
         deltaHi(i) *= 1./3.;
         deltaLo(i) *= 1./3.;
      }
      cout << "Delta(0): " << deltaHi(0) << ", " << deltaLo(0) << endl;
      Hamiltonian = 1;
      if (hamiltonian->containsGroup("strain"))  {
         strained = true;
         fileNames.resize(6);
         SxSymbolTable *strain = hamiltonian->getGroup("strain");
         factor       = strain->get("factor")->toReal();
         DHi          = SxVector<Double>(strain->get("DHi")->toList ());
         DLo          = SxVector<Double>(strain->get("DLo")->toList ());
         aHi          = SxVector<Double>(strain->get("aHi")->toList ());
         aLo          = SxVector<Double>(strain->get("aLo")->toList ());
         fileNames(0) = (strain->get("sigmaXX")->toString());
         fileNames(1) = (strain->get("sigmaXY")->toString());
         fileNames(2) = (strain->get("sigmaXZ")->toString());
         fileNames(3) = (strain->get("sigmaYY")->toString());
         fileNames(4) = (strain->get("sigmaYZ")->toString());
         fileNames(5) = (strain->get("sigmaZZ")->toString());
         for (i = 0; i < 6; i++)  {
            DHi(i) /= HA2EV;  
            DLo(i) /= HA2EV;  
         }
         for (i = 0; i < 2; i++)  {
            aHi(i) /= HA2EV;  
            aLo(i) /= HA2EV;  
         }
      }
      if (hamiltonian->containsGroup("externalField"))  {
         SxSymbolTable *extField = hamiltonian->getGroup("externalField");
         SxString fileName = (extField->get("file")->toString() );
         SxBinIO io (fileName, SxBinIO::BINARY_READ_ONLY);
         vExtR.resize (1);
         vExtR = io.readMesh (&cell, &dim);
         externalField = vExtR(0)*(1./HA2EV);
         io.close();
      }
      else 
         externalField = zero;
   }
   // --- piezoelectric field 
   if (hamiltonian->containsGroup("peField"))  {
      SxSymbolTable *peField = hamiltonian->getGroup("peField");
      SxString fileName = (peField->get("file")->toString() );
      SxBinIO io (fileName, SxBinIO::BINARY_READ_ONLY);
      vExtR.resize (1);
      vExtR = io.readMesh (&cell, &dim);
      externalField += (vExtR(0)*(1./HA2EV));
      io.close();
      cout << "piezoelectric field, max val: " <<
         externalField.maxval() << ", min val: " <<
         externalField.minval() << endl;
   } 
   if (hamiltonian->contains("nEmptyStates"))
      nEmptyStates = hamiltonian->get("nEmptyStates")->toInt();
   if (hamiltonian->contains("electrons"))
         electrons = hamiltonian->get("electrons")->toAttribute();
   if (hamiltonian->contains("bandstructure"))
         bs = hamiltonian->get("bandstructure")->toAttribute();
   if (hamiltonian->contains("position"))
         posBS = SxVector3<Double>(hamiltonian->get("position")->toList());
   else {posBS(0) = 0.; posBS(1) = 0.; posBS(2) = 0.;}
   cout << "reading of K*P-input was successful...\n";
   egHi = ecHi - evHi; // the band gap
   egLo = ecLo - evLo;
   cout << "Band gaps: " << egLo * HA2EV << " and " << egHi * HA2EV << endl;
   if (minimum - maximum < 1e-20)
      cout << "WARNING: external potential seems to be zero!\n";
   cout << "reading external potential was successfull...\n";
   
   
   eVBO.resize(externalPotential.getSize());
   eCBO.resize(externalPotential.getSize());
   if (Hamiltonian == 0)  {
      cout << "Calculating in Zincblend material\n";
      p0.resize(externalPotential.getSize());
      eP.resize(externalPotential.getSize());
      pO.resize(1);
      pP.resize(1);
   }

   if (Hamiltonian == 1)  {
      cout << "Calculating in Wurtzite material\n";
      p0.resize(1);
      eP.resize(1);
      pO.resize(externalPotential.getSize()); // refers to P1 in ref. [Winkelnkemper]
      pP.resize(externalPotential.getSize()); // refers to P2 in ref. [Winkelnkemper]
   }

   // --- calculation done as in Phys.Rev. B 74, 155322 (2006) [Winkelnkemper]
   // --- delta(0) = deltaCR, delta(1) = delta(2) = deltaSO;
   if (Hamiltonian == 1)  {
      cout << "moHi: " << mOHi << ", mpHi : " << mPHi << endl;
      cout << "moLo: " << mOLo << ", mpLo : " << mPLo << endl;
      pOHi = 
         sqrt(fabs(.5 * (1./mOHi - 1.) * ( 3. * egHi *
                   (deltaHi(0) + egHi)
                   + deltaHi(1) * (2. * deltaHi(0) + 3. * egHi))
                   / (2. * deltaHi(1) + 3. * egHi)
                   )
               );
      pOLo = 
          sqrt(fabs(.5 * (1./mOLo - 1.) * ( 3. * egLo *
                   (deltaLo(0) + egLo)
                   + deltaLo(1) * (2. * deltaLo(0) + 3. * egLo))
                   / (2. * deltaLo(1) + 3. * egLo)
                   )
               );
     
      pPHi = sqrt(fabs(.5 * (1./mPHi - 1.)
               * (egHi *
                  ( 3. * egHi * ( deltaHi(0) + egHi)
                    + deltaHi(1) * (2. * deltaHi(0) + 3. * egHi)))
               / ( deltaHi(0) * deltaHi(1) + 3. * deltaHi(1) * egHi
                  + 2. * deltaHi(1) * egHi + 3. * sqr(egHi))
               ) );
      pPLo = sqrt(fabs(.5 * (1./mPLo - 1.)
               * ( egLo *
                  ( 3. * egLo * ( deltaLo(0) + egLo)
                    + deltaLo(1) * (2. * deltaLo(0) + 3. * egLo)))
               / ( deltaLo(0) * deltaLo(1) + 3. * deltaLo(1) * egLo
                  + 2. * deltaLo(1) * egLo + 3. * sqr(egLo))
               ) );
      if (hamiltonian->contains("eP1Hi"))  {
         pPHi = sqrt(.5 * eP2Hi);
         pOHi = sqrt(.5 * eP1Hi);
         pPLo = sqrt(.5 * eP2Lo);
         pOLo = sqrt(.5 * eP1Lo);
      }
      
      cout << "ePOHi: " << 2. * sqr(pOHi) * HA2EV << " eV" << endl;
      cout << "ePOLo: " << 2. * sqr(pOLo) * HA2EV << " eV" << endl;
      cout << "ePPHi: " << 2. * sqr(pPHi) * HA2EV << " eV" << endl;
      cout << "ePPLo: " << 2. * sqr(pPLo) * HA2EV << " eV" << endl;
   }
   if (Hamiltonian == 0)  {
      reciprocalEffMass.resize(externalPotential.getSize());
      deltaSpinOrbital.resize(externalPotential.getSize());
   }  else  {
      reciprocalEffMass.resize(1);
      deltaSpinOrbital.resize(1);
   }
   for (i = 0; i < 3; i++)  {
      if (Hamiltonian == 0)
         gamma(i).resize(externalPotential.getSize() );
      else 
         gamma(i).resize(1);
      if (Hamiltonian == 1)
         delta(i).resize(externalPotential.getSize() );
      else 
         delta(i).resize(1);
   }
   
      // --- resize WZ parameter meshes
   if (Hamiltonian == 1)  {
      P1.resize(mix.getSize());
      P2.resize(mix.getSize());
      for (i = 0; i < 7; i++)  {
         if (i < 6)
            eVBOMod(i).resize(mix.getSize());
         A(i).resize(mix.getSize());
         A1.resize(mix.getSize());
         A2.resize(mix.getSize());
         L1.resize(mix.getSize());
         L2.resize(mix.getSize());
         M1.resize(mix.getSize());
         M2.resize(mix.getSize());
         M3.resize(mix.getSize());
         N1.resize(mix.getSize());
         N2.resize(mix.getSize());
         N3.resize(mix.getSize());

         }
   } else  {
      P1.resize(1);
      P2.resize(1);
   }
   sxprintf("ecHi %f, evHi %f, ecLo %f, evLo %f\n\n", ecHi, evHi, ecLo, evLo);
   cout << "mix: " << mix.minval() << " to " << mix.maxval() << endl;
   double alphaHi, alphaLo;

   alphaHi = 1./eMassHi - (ePHi/3.) * (2./(ecHi -evHi)
         + 1. / (ecHi - evHi + deltaSOHi));
   alphaLo = 1./eMassLo - (ePLo/3.) * (2./(ecLo -evLo)
         + 1. / (ecLo - evLo + deltaSOLo));
   cout << "s: " << alphaHi << " to " << alphaLo << endl;
   cout << "eP: " << ePHi << " to " << ePLo << endl;
   rMassP = 1./eMassLo;
   for (idx = 0; idx < externalPotential.getSize(); idx++)  {
      eCBO(idx) = (mix(idx) * (ecHi - ecLo) + ecLo);
      eVBO(idx) = (mix(idx) * (evHi - evLo) + evLo);
      if (Hamiltonian == 0)  {
         reciprocalEffMass(idx) = 
            (mix(idx) * (alphaHi - alphaLo) + alphaLo);
         eP(idx)   =  (mix(idx) * (ePHi - ePLo) + ePLo);
         p0(idx) =  sqrt(fabs(eP(idx)/2.));
         deltaSpinOrbital(idx) = (mix(idx)
               * (deltaSOHi - deltaSOLo) + deltaSOLo);
         for (i = 0; i< 3; i++)  {
            (gamma(i))(idx) =  
               (mix(idx) * (modLutHi(i) - modLutLo(i)) + modLutLo(i));
         }
      } else if (Hamiltonian == 1)  {
         // --- parameter meshes for WZ are calculated here
         pO(idx)   = (mix(idx) * (pOHi - pOLo) + pOLo);
         pP(idx)   = (mix(idx) * (pPHi - pPLo) + pPLo);
         for (i = 0; i < 3; i++)
            (delta(i))(idx) = mix(idx) * (deltaHi(i) - deltaLo(i)) + deltaLo(i);

         double hi, lo; // will be re-used for every parameter!!!


         // --- A'1
         hi = .5 / mPHi - sqr(pPHi) / (ecHi - evHi);
         lo = .5 / mPLo - sqr(pPLo) / (ecLo - evLo);
         aP1 = lo;
         A1(idx) = mix(idx) * (hi - lo) + lo;
         
         // --- A'2
         hi = .5 / mOHi - sqr(pOHi) / (ecHi - evHi);
         lo = .5 / mOLo - sqr(pOLo) / (ecLo - evLo);
         aP2 = lo;
         A2(idx) = mix(idx) * (hi - lo) + lo;
         
         if (idx == 0)  {
            aP.resize(7);
            AHi(0) = AHi(0) + 2. * sqr(pOHi) / (ecHi - evHi);
            ALo(0) = ALo(0) + 2. * sqr(pOLo) / (ecLo - evLo);
            aP(0) = ALo(0);

            // A(1) remains as is
            aP(1) = ALo(1);
         
            AHi(2) = AHi(2) - 2. * sqr(pOHi) / (ecHi - evHi);
            ALo(2) = ALo(2) - 2. * sqr(pOLo) / (ecLo - evLo);
            aP(2) = ALo(2);

            AHi(3) = AHi(3) + sqr(pPHi) / (ecHi - evHi);
            ALo(3) = ALo(3) + sqr(pPLo) / (ecLo - evLo);
            aP(3) = ALo(3);

            AHi(4) = AHi(4) + sqr(pPHi) / (ecHi - evHi);
            ALo(4) = ALo(4) + sqr(pPLo) / (ecLo - evLo);
            aP(4) = ALo(4);

            AHi(5) = AHi(5) + sqrt(2.) * (pOHi * pPHi) / (ecHi - evHi);
            ALo(5) = ALo(5) + sqrt(2.) * (pOLo * pPLo) / (ecLo - evLo);
            aP(5) = ALo(5);
         }
          for ( i = 0; i < 7; i++)  {
            A(i)(idx) = -(mix(idx) * (AHi(i) - ALo(i)) + ALo(i));
            if (idx == 0)
               cout << "A(" << i+1 << ")= "
                  << AHi(i) << " ... " << ALo(i) << endl;
          }
      }
   }
   if (electrons) {
      cout << "calculating electron states\n";
      precFactor = .001;
      if (eTrial == 1.e50)
         eTrial = eCBO.minval() + externalField.minval();
   } else {
      precFactor = 1.e-10;
      cout << "calculating hole states\n";
      if (eTrial == 1.e50)
         eTrial = eVBO.maxval() + externalField.maxval();
   }  
   cout << "eTrial: " << eTrial * HA2EV << " eV" << endl; 

   cout << "eCBO: " << eCBO.minval()*HA2EV << " to " << eCBO.maxval()*HA2EV << endl;
   if (Hamiltonian == 0)
      cout << "p0: " << p0.maxval() << " to " << p0.minval() << endl;
   
   cout << "conduction band offset: high pot.: " << ecHi * HA2EV << " low pot.: "
        << ecLo * HA2EV << " eV\n"
        << "\nvalence band offset: high pot.: " << evHi * HA2EV 
        << " low pot.: " << evLo * HA2EV << " eV" << endl;
   if (Hamiltonian == 0)
      cout << "modified luttinger parameters:" << endl 
         << "gammaHi: " << "[ " << gamma(0).maxval() << ", "
         << gamma(1).maxval() << ", " << gamma(2).maxval() << " ]" << endl
         << "gammaLo: " << "[ " << gamma(0).minval() << ", "
         << gamma(1).minval() << ", " << gamma(2).minval() << " ]" << endl;
   if (Hamiltonian == 1)
      cout << "effective Masses: " << A1.maxval() << " to " << A1.minval()
         << endl << "and " << A2.maxval() << " to " << A2.minval() << endl
         << "Delta from: [" << delta(0).maxval() << ", " << delta(1).maxval()
         << ", " << delta(2).maxval() << "]" << endl
         << "to: [" << delta(0).minval() << ", " << delta(1).minval()
         << ", " << delta(2).minval() << "]" << endl;
   
   if (Hamiltonian == 1)  {
      eVBOMod(0) = eVBO + delta(0) + delta(1);
      eVBOMod(1) = eVBO;
      eVBOMod(2) = eVBO + delta(0) - delta(1);
      eVBOMod(3) = eVBO + delta(0) - delta(1);
      eVBOMod(4) = eVBO;
      eVBOMod(5) = eVBO + delta(0) + delta(1);
   }
   if (strained)
      calculateStrain (
            fileNames,
            aHi, aLo,
            DHi, DLo,
            ZBParHi, ZBParLo,
            factor,
            Hamiltonian);
   // --- calculation for given k-Point
   double kr, phi;
   kxSet = kySet = kzSet = false;
   kx = ky = kz = kr = phi = 0.;
   if (hamiltonian->contains("kx"))  {
      kx  = hamiltonian->get("kx")->toReal();
      kxSet = true; 
   }
   if (hamiltonian->contains("ky"))  {
      ky  = hamiltonian->get("ky")->toReal();
      kySet = true; 
   }
   if (hamiltonian->contains("kz"))  {
      kz  = hamiltonian->get("kz")->toReal();
      kzSet = true; 
   }
   if (hamiltonian->contains("kr"))  {
      kr  = hamiltonian->get("kr")->toReal();
      phi = hamiltonian->get("phi")->toReal();
      kxSet = true;
      kySet = true;
      kx = cos(phi) * kr;
      ky = sin(phi) * kr;
   }
   if (kxSet || kySet || kzSet)
      cout << "calculation for fixed k: (";
   if (kxSet)
      cout << kx << ", ";
   else cout << "-, ";
   if (kySet)
      cout << ky << ", ";
   else cout << "-, ";
   if (kzSet)
      cout << kz << ",)" << endl;
   else cout << "-)" << endl;

   // --- calculate memory consumption M(dx,dy,dz,nStates), empirical
   double memConsumption = 0.;
   if (Hamiltonian == 0)
      memConsumption = (dim(0) * dim(1) * dim(2))
         * ( 2.3125e-3 + 0.0625e-3 * nStates);
   if (Hamiltonian == 1)
      memConsumption = (dim(0) * dim(1) * dim(2))
         * ( 2.09375e-3 + 0.0625e-3 * nStates);
   cout << "\n\nApproximated memory consumption: " << round(memConsumption) << " MB\n\n";
}

void SxKdotP8x8::printEnergies () const
{
   
   cout.precision(5);
   const SxPWSet &waves = getWavesRef();
   PsiG comp, psiI;
   int i, iComp;
   double sum, sum0, sum1;
   cout << "k*p-Hamiltonian was called " << step/nStates << " times in calculation.\n";
   if (Hamiltonian == 0)
      cout << "State\tS up\t S dn\thh up\thh dn\tlh up\tlh dn\tso up\tso dn\n";
   else if (Hamiltonian == 1)
      cout << "State\tS up\t S dn\tP up\tZ up\tM up\tP dn\tZ dn\tM dn\n";
   sum0 = 0.;
   sum1 = 0.;
   for (i=0; i < nStates; i++)  {
      cout << i << " :\t";\
      for (iComp = 0; iComp < 8; iComp++)  {
  
         // TO BE DONE
         psiI = waves(i,0,0);
         
         comp = psiI.getComponent(iComp);
         
         sum = dot(comp, comp);
         if (iComp == 0)
            sum0 += sum;
         if (Hamiltonian == 0) {
            if (iComp == 1) sum1 += sum;
         } else if (Hamiltonian == 1)  {
            if (iComp == 4) sum1 += sum;
         }
         cout << (int(sum*1000)) / 10. << "\t";
      }
      cout << "\n";
   }
   cout << endl;
}

SxRho &SxKdotP8x8::getRho ()
{
   return *this;
}

void SxKdotP8x8::computeRho (const SxPsiSet &wavesIn, const SxFermi &fermi)
{  
   SX_CHECK (dynamic_cast<const SxPW *> (&wavesIn));
   const SxPW &wavesRef = *dynamic_cast<const SxPW *> (&wavesIn);
   SxRho::computeRho (fermi.focc, wavesRef);
}

PrecEnergy SxKdotP8x8::getEnergy (const SxPsiSet &psiSet,
                                  const SxFermi &fermi)
{
   SX_CHECK (dynamic_cast<const SxPWSet *>(&psiSet));   
   wavesPtr = static_cast<const SxPWSet *> (&psiSet);

   update (fermi);
   
   return eKin; 
}

void SxKdotP8x8::calculateStrain (
               SxArray<SxString> &inputFiles,
               SxVector<Double> &aHi,
               SxVector<Double> &aLo,
               SxVector<Double> &DHi,
               SxVector<Double> &DLo,
               SxVector<Double> &ZBParLo,
               SxVector<Double> &ZBParHi,
               double factor,
               bool Ham)
{
   cout << "--- calculate Strain contribution ---" << endl;
   if (Ham)
      cout << "zincblende structure" << endl;
   else
      cout << "wurtzite structure" << endl;
   cout << "nInputFiles: " << inputFiles.getSize() << endl;
   cout << "a(0): " << aHi(0) << " to " << aLo(0) << endl;
   cout << "D(0): " << DHi(0) << " to " << DLo(0) << endl;
   cout << "ZBParameter(0): " << ZBParHi(0) << " to " << ZBParLo(0) << endl;
   cout << "factor: " << factor << endl;
   /*
   int ispin, col, row, idx;
   sxarray<sxmeshr> vextr; // strain field components
   sxarray<psir> sigma;
   sxarray<sxmeshr> d; // d1-6;
   sxmeshr mix;
   double hi, lo;
   double minimum, maximum;
   
   cout << "calculating strain contribution...\n" << endl;
   
   // --- resize required meshr's
   h_strain.resize(nspin);
   sigma.resize(6);

   // --- sigma array:
   // xx     - 0
   // xy, yx - 1
   // xz, zx - 2
   // yy     - 3
   // yz, zy - 4
   // zz     - 5
   
   mix.resize(externalpotential.getsize());
   minimum = - externalpotential.minval();
   maximum = - externalpotential.maxval();
   for (idx = 0; idx < externalpotential.getsize(); idx++) 
      mix(idx) = mixing(minimum, 1., maximum, 0., externalpotential(idx));
   cout << "mixing from: " << mix.minval() << " to " << mix.maxval() << endl;
   
   // --- initialize h_strain with 0
   for (ispin = 0; ispin < nspin; ispin++)
      h_strain(ispin).resize(nspin);
   for (col = 0; col < nspin; col++)
      for (row = 0; row < nspin; row++)
         h_strain(col)(row) = zero;

   // --- resize strain field arrays 
   for (idx = 0; idx < 6; idx++)
      sigma(idx).resize(externalpotential.getsize());

   if (hamiltonian == 1)  { // wurtzite material
      d.resize(6);
      int i;
      for (i = 0; i < 6; i++)
         d(i).resize(externalpotential.getsize());
      sxmeshr a1, a2;// l1, l2, m1, m2, m3, n1, n2; 
      a1.resize(externalpotential.getsize());
      a2.resize(externalpotential.getsize());
      psir h, hc, k, kc;
      h.resize(externalpotential.getsize());
      k.resize(externalpotential.getsize());
      hc.resize(externalpotential.getsize());
      kc.resize(externalpotential.getsize());
      
      // --- read data from files
      
      for (idx = 0; idx < 6; idx++)  {
         sxio io (inputfiles(idx), sxio::binary_read_only);
         sxmatrix3<double> cell;
         sxvector3<int>    dim;
         vextr.resize (1);
         vextr = io.readmesh (&cell, &dim);
         sigma(idx)= vextr(0) * 1./factor;
         cout << "open " << inputfiles(idx) << endl;
         io.close();
         cout << "d[" << idx << "]: " << dhi(idx) << " to " << dlo(idx) << endl;
      }
      
      // --- determine a,l,m,n
     
      for (idx = 0; idx < externalpotential.getsize(); idx++)  {
         hi = ahi(0);
         lo = alo(0);
         a1(idx) =  -  (mix(idx) * (hi - lo) + lo);
         
         hi = ahi(1);
         lo = alo(1);
         a2(idx) =  -  (mix(idx) * (hi - lo) + lo);
         
         for (i = 0; i < 6; i++)
            d(i)(idx) = - (mix(idx) * (dhi(i) - dlo(i)) + dlo(i));
         
      }
      // --- print parameter
      cout << "strain parameter set: \n"
         << "a1 (perp): " << a1.minval() << " to " << a1.maxval() << endl
         << "a2 (par) : " << a2.minval() << " to " << a2.maxval() << endl << endl;
     
      // --- the hamiltonian itself
      // --- diagonal terms
      h_strain(0)(0) = a2 * (sigma(0) + sigma(3)) + a1 * sigma(5);
      h_strain(1)(1) = a2 * (sigma(0) + sigma(3)) + a1 * sigma(5);
      h_strain(2)(2) = (d(0) + d(2)) * sigma(5)
         + (d(1) + d(3)) * (sigma(0) + sigma(3));
      h_strain(3)(3) = d(0) * sigma(5)
         + d(1) * (sigma(0) + sigma(3));
      h_strain(4)(4) = (d(0) + d(2)) * sigma(5)
         + (d(1) + d(3)) * (sigma(0) + sigma(3));
      h_strain(5)(5) = (d(0) + d(2)) * sigma(5)
         + (d(1) + d(3)) * (sigma(0) + sigma(3));
      h_strain(6)(6) = d(0) * sigma(5)
         + d(1) * (sigma(0) + sigma(3));
      h_strain(7)(7) = (d(0) + d(2)) * sigma(5)
         + (d(1) + d(3)) * (sigma(0) + sigma(3));
     
      // --- upper nondiagonal terms
      h  = d(5) * (sigma(2) + i * sigma(4));
      hc = d(5) * (sigma(2) - i * sigma(4));
      k  = d(4) * (sigma(0) - sigma(3) + 2. * i * sigma(1));
      kc = d(4) * (sigma(0) - sigma(3) - 2. * i * sigma(1));

      h_strain(2)(3) = -hc;
      h_strain(2)(4) = -kc;
      h_strain(3)(2) = -h; 
      h_strain(3)(4) = hc;
      h_strain(4)(2) = -k;
      h_strain(4)(3) = h;

      h_strain(5)(6) = -hc;
      h_strain(5)(7) = -kc;
      h_strain(6)(5) = -h; 
      h_strain(6)(7) = hc;
      h_strain(7)(5) = -k;
      h_strain(7)(6) = h;

      // --- lower nondiagonal terms
     
   } else if (hamiltonian == 0)  { // zincblend material
      sxmeshr ac, av, b, d;
      ac.resize(externalpotential.getsize());
      av.resize(externalpotential.getsize());
      b.resize(externalpotential.getsize());
      d.resize(externalpotential.getsize());

      // --- read data from files
      for (idx = 0; idx < 6; idx++)  {
         sxio io (inputfiles(idx), sxio::binary_read_only);
         sxmatrix3<double> cell;
         sxvector3<int>    dim;
         vextr.resize (1);
         vextr = io.readmesh (&cell, &dim);
         sigma(idx)= vextr(0);
         cout << "open " << inputfiles(idx) << endl;
         io.close();
      }

      // --- assign material parameters
      for (idx = 0; idx < externalpotential.getsize(); idx++)  {
         hi = zbparhi(0);
         lo = zbparlo(0);
         ac(idx) =  (mix(idx) * (hi - lo) + lo);
         
         hi = zbparhi(1);
         lo = zbparlo(1);
         av(idx) =  (mix(idx) * (hi - lo) + lo);
         
         hi = zbparhi(2);
         lo = zbparlo(2);
         b(idx) =  (mix(idx) * (hi - lo) + lo);
         
         hi = zbparhi(3);
         lo = zbparlo(3);
         d(idx) =  (mix(idx) * (hi - lo) + lo);
      }

      cout << "zincblende strain parameter set: \n";
      cout << "ac: " << ac.maxval() << " to " << ac.minval() << endl;
      cout << "av: " << av.maxval() << " to " << av.minval() << endl;
      cout << "b: " << b.maxval() << " to " << b.minval() << endl;
      cout << "d: " << d.maxval() << " to " << d.minval() << endl;
      
      // --- the hamiltonian

      psir p, q, r, s, r_conj, s_conj, u, v, u_conj, v_conj;
      
      p = av * (sigma(0) + sigma(3) + sigma(5));
      q = b * (sigma(5) - .5 * (sigma(0) + sigma(3)));
      r = s32 * b * (sigma(0) - sigma(3)) - i * d * sigma(1);
      r_conj = s32 * b * (sigma(0) - sigma(3)) + i * d * sigma(1);
      s = -d * (sigma(2) - i * sigma(4));
      s_conj = -d * (sigma(2) + i * sigma(4));

      u = zero;
      u_conj = zero;
      v = zero;
      v_conj = zero;

      // --- diagonal elements
      h_strain(0)(0) = ac * (sigma(0) + sigma(3) + sigma(5));
      h_strain(1)(1) = ac * (sigma(0) + sigma(3) + sigma(5));
      h_strain(2)(2) = -p + q;
      h_strain(3)(3) = -p - q;
      h_strain(4)(4) = -p - q;
      h_strain(5)(5) = -p + q;
      h_strain(6)(6) = -av * (sigma(0) + sigma(3) + sigma(5));
      h_strain(7)(7) = -av * (sigma(0) + sigma(3) + sigma(5));
      
      // --- nonzero nondiagonal elements

      h_strain(0)(2) = -v_conj;          h_strain(2)(0) = -v;
                                         
      h_strain(0)(4) = -s3 * v;          h_strain(4)(0) = -s3 * v_conj;
                                         
      h_strain(0)(5) = s2 * u;           h_strain(5)(0) = s2 * u_conj;
                                         
      h_strain(0)(6) = u;                h_strain(6)(0) = u_conj;
                                         
      h_strain(0)(7) = -s2 * v_conj;     h_strain(7)(0) = -s2 * v;
                                         
      h_strain(1)(2) = s2 * u;           h_strain(2)(1) = s2 * u_conj;
                                         
      h_strain(1)(3) = s3 * v_conj;      h_strain(3)(1) = s3 * v;
                                         
      h_strain(1)(5) = v;                h_strain(5)(1) = v_conj;
                                         
      h_strain(1)(6) = -s2 * v;          h_strain(6)(1) = -s2 * v_conj;
                                         
      h_strain(1)(7) = -u;               h_strain(7)(1) = -u_conj;
                                         
      h_strain(2)(3) = -s_conj;          h_strain(3)(2) = -s;
                                         
      h_strain(2)(4) = r;                h_strain(4)(2) = r_conj;
                                         
      h_strain(2)(6) = s32 * s;          h_strain(6)(2) = s32 * s_conj;
                                         
      h_strain(2)(7) = -s2 * q;          h_strain(7)(2) = -s2 * q;
                                         
      h_strain(3)(5) = r;                h_strain(5)(3) = r_conj;
                                         
      h_strain(3)(6) = -s2 * r;          h_strain(6)(3) = -s2 * r_conj;
                                         
      h_strain(3)(7) = 1./s2 * s;        h_strain(7)(3) = 1./s2 * s_conj;
                                         
      h_strain(4)(5) = s_conj;           h_strain(5)(4) = s;
 
      h_strain(4)(6) = 1./s2 * s_conj;   h_strain(6)(4) = 1./s2 * s;

      h_strain(4)(7) = s2 * r_conj;      h_strain(7)(4) = s2 * r;

      h_strain(5)(6) = s2 * q;           h_strain(6)(5) = s2 * q;

      h_strain(5)(7) = s32 * s_conj;     h_strain(7)(5) = s32 * s;

   }

   // --- print hamiltonian
   double elem;
   cout << "strain hamiltonian: \n";
   char color[13];
   for (col = 0; col < nspin; col++)  {
      for (row = 0; row < nspin; row++)  {
               elem = h_strain(col)(row).sum().re;
               if ((col < 2) || (row < 2))  {
                  sprintf(color, "%c[%d;%d;%dm", 0x1b, 0, 31, 40);
                  sxprintf("%s", color);
               }  else  {
                  sprintf(color, "%c[%d;%d;%dm", 0x1b, 0, 33, 40);
                  sxprintf("%s", color);
               }
               if (col == row)  {
                  sprintf(color, "%c[%d;%d;%dm", 0x1b, 0, 35, 40);
                  sxprintf("%s", color);
               }
               sxprintf("%*.3e\t", 5, elem);
               sprintf(color, "%c[%d;%d;%dm", 0x1b, 0, 42, 40);
               sxprintf(color);
      }
      cout << endl;
   }
   */
}
      

void SxKdotP8x8::set (const SxPsiSet &psiSet, const SxFermi &fermi)
{
   SX_CHECK (dynamic_cast<const SxPWSet *> (&psiSet));
   wavesPtr = dynamic_cast<const SxPWSet *>(&psiSet);
   compute (fermi, true, true);
   fermiPtr = &fermi;
}

void SxKdotP8x8::readRho (const SxBinIO &io)
{
   SxRho::readRho (io);
}

void SxKdotP8x8::normalizeRho ()
{
   SxRho::normalizeRho ();
}

void SxKdotP8x8::writeRho (const SxString &filename) const
{
   SxRho::writeRho (filename);
}

SxDiracVec<TPrecCoeffG::TReal> 
SxKdotP8x8::preconditioner (const PsiG &psi,
                                 Preconditioner type) const
{
   if (type == 0) {/*get rid of compiler warnings*/}
   SxDiracVec<TPrecCoeffG::TReal> x, x2, x3, x4, n, K;
   x.resize(psi.getSize());
   x = D / ((D ^ psi.absSqr()).chop());
   x2 = x.sqr();
   x3 = x.cub();
   x4 = x2.sqr();
   n  = 27. + 18.*x + 12.*x2 + 8.*x3;
   K  = n / (n + 16.*x4);

   return K; 
}


