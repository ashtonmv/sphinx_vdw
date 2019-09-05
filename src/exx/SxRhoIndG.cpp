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

#include <SxRhoIndG.h>
#include <SxProjector.h>

SxRhoIndG::SxRhoIndG ()
   : GPtr    (NULL),
     gkPtr   (NULL),
     nk      (-1)
{
   // --- set uninitialized variables to nonsense values
   nv = nc = nG = -1;
   eX = 0.;
   fftFactorRtoG = 0.;
}

SxRhoIndG::SxRhoIndG (const SxGBasis &G, const SxGkBasis &gk, int nG_)
   : GPtr    (&G),
     gkPtr   (&gk),
     nk      (gkPtr->getNk())
{
   // --- set uninitialized variables to nonsense values
   nv = nc = -1;
   eX = 0.;

   // --- initialise R basis
   R = G.getRBasis ();

   // --- set number of G components used
   nG = nG_;
   SX_CHECK (nG % 2 == 1);    // An odd number of G components has to be used.

   // --- resize rhoInd
   rhoIndG.resize (nG-1);  // The (G = 0) component is not considered!
   rhoIndG.setBasis (GPtr);

   // --- fft prefactor
   //   ( To check the validity of this factor, use SxRhoIndG::checkFFTfactor! )
   fftFactorRtoG = 1. / R.fft3d.scaleFor;
}

void SxRhoIndG::checkFFTfactor (const SxPW &waves)
{
   SX_CHECK (GPtr);

   const SxGBasis    *GkPtr;  // Don't mix it up with *gkPtr!
   const SxGBasis    &G    = *GPtr;
   PsiG               psiG = waves (0,0,0);
   PrecCoeffG         g0;
   double             g0abs, eps = 1.e-10;

   GkPtr = &( (*gkPtr)(0) );
   psiG.setBasis (GkPtr);

   PsiR psiR  = ( R | psiG );
   PsiG rhoVV = ( G | (psiR * psiR.conj()) );
   g0    = rhoVV(0) * fftFactorRtoG;
   g0abs = sqrt (g0 * g0.conj());

   if (fabs(g0abs - 1.) < eps)  {
      cout << "SxRhoIndG::testFFTfactor --> (G=0)-component: " << g0
           << " ... ok" << endl;
   }  else  {
      cout << "SxRhoIndG::testFFTfactor --> (G=0)-component: " << g0
           << " ... failed" << endl;
      SX_QUIT;
   }
}

void
SxRhoIndG::compute(const SxPW &waves, const SxFermi &fermi, SxFockGk &FockOp)
{
   SX_CHECK  (waves.getNSpin() == 1, waves.getNSpin());
               // spin pol. not yet implemented
   SX_CHECK (waves.getNk() == nk, waves.getNk(), nk);
   SX_CHECK      (GPtr);

   int              ik;
   int              iG;
   int              iv, ic, vIdx, cIdx;
   PsiG             vPsiG, cPsiG, rhoVC, vFock;
   PsiR             vPsiR, cPsiR;
   PrecCoeffG       vcFock;
   PrecEps          vEps, cEps, vcEps;
   const SxGBasis  *GkPtr;  // Don't mix it up with *gkPtr!
   const SxGBasis  &G = *GPtr;
   PrecWeights      kWeight;

   // - - - - - - - - - - - - - - - - - - - - - - - - - -
   nv = fermi.getNValenceBands(0,0);
   nc = fermi.getNConductionBands(0,0);

   // TODO: The following lines should go to SxFermi.
   // --- check for equality of nv at all k-points
   for (ik = 0; ik < nk; ik++)  {
      if (fermi.getNValenceBands(0,ik) != nv)  {
         sxprintf ("Number of valence states is not the same for all k-points. ");
         sxprintf ("Sorry.\n");
         SX_QUIT;
      }
   }

   // --- find the smallest number of conduction bands
   for (ik = 0; ik < nk; ik++)  {
      if (fermi.getNConductionBands(0,ik) < nc)  {
         nc = fermi.getNConductionBands(0,ik);
      }
   }
   // - - - - - - - - - - - - - - - - - - - - - - - - - -

   sxprintf ("SxRhoIndG::compute ... using %d G vectors "
           "and %d conduction bands\n", nG, nc);
   cout << "SxRhoIndG::compute ... starting\n";  cout.flush();

   rhoIndG.set (0.);
   eX = 0.;

   for (ik = 0; ik < nk; ik++)  {
      GkPtr   = &( (*gkPtr)(ik) );
      kWeight = gkPtr->weights(ik);
      FockOp.compute (ik, waves, fermi);

      vPsiG.setBasis (GkPtr);
      cPsiG.setBasis (GkPtr);

      for (iv = 0; iv < nv; iv++)  {
         vIdx  = fermi.getValenceBandIdx (iv,0,ik);
         vEps  = fermi.eps (vIdx,0,ik);
         vPsiG = waves (vIdx,0,ik);
         vPsiR = ( R | vPsiG );
         vFock = FockOp * vPsiG;

         // --- exchange contribution to total energy
         eX += kWeight * (PrecEnergy)(vFock ^ vPsiG).chop();

         for (ic = 0; ic < nc; ic++)  {
            cIdx  = fermi.getConductionBandIdx (ic,0,ik);
            cEps  = fermi.eps (cIdx,0,ik);
            cPsiG = waves (cIdx,0,ik);
            cPsiR = ( R | cPsiG );

            rhoVC  = ( G | (vPsiR * cPsiR.conj()) );  // like in Abdullah's code
//          vcFock = ( vPsiG ^ (FockOp * cPsiG) ).chop();
            vcFock = (vFock ^ cPsiG).chop ();
            vcEps  = 2. * kWeight / (vEps - cEps);

            // --- compute E(G) + E(-G)^(*)
            for (iG = 0; iG < nG-1; iG += 2)  {
               rhoIndG(iG+0) += vcEps * (   vcFock * rhoVC(iG+2)
                                        + (vcFock * rhoVC(iG+1)).conj()  );
               rhoIndG(iG+1) += vcEps * (   vcFock * rhoVC(iG+1)
                                        + (vcFock * rhoVC(iG+2)).conj()  );
            }  // :iG
         }  // :ic
      }  // :iv
   }  // :ik

   // --- multiply with FFT factor
   // 
   //   ( Within the upper calculation the following scheme is applied:
   //     <G|psi_i> --> <R|psi_i> --> <G|psi_i*psi_j> =: rhoVC without
   //     considering any FFT prefactor. However, in the second trafo
   //     occurs a (pointwise) product of two wavefunctions. For one
   //     function, missing factors cancel due to the consecutive
   //     performance of backward and forward trafo. For the other
   //     wavefunction, however, one needs to care about the FFT factors. )
   rhoIndG *= fftFactorRtoG;

   // --- divide by omega
   rhoIndG *= R.cell.volume;  // ???????????? Why * and not / ?
   rhoIndG = rhoIndG.conj();    // Why????

   cout << "SxRhoIndG::compute ... done\n";  cout.flush();
}

SxMeshG SxRhoIndG::getRhoIndG () const
{
   SX_CHECK (rhoIndG.getSize() > 0, rhoIndG.getSize());
   return (rhoIndG);
}

SxMeshR SxRhoIndG::getRhoIndR () const
{
   SX_CHECK (rhoIndG.getSize() > 0, rhoIndG.getSize());
   return ( R | rhoIndG );
}

PrecEnergy SxRhoIndG::getXEnergy () const
{
   SX_CHECK (rhoIndG.getSize() > 0, rhoIndG.getSize());
   return eX;
}
