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
//Ref1: PRB (54), 11169-11186 (1996)
//Ref2: Comp. Mat. Sci. (6), 15-50 (1996)

#include <SxRhoMixerG.h>
#include <SxSymMatrix.h>

SxRhoMixerG::SxRhoMixerG (MixerType typeIn, double rhoMixingIn, int maxStepsIn, 
      double spinMix)
   : maxSteps (maxStepsIn),
     type (typeIn),
     renormModus (SxRhoMixerG::RenormOn),
     spinMixing (spinMix),
     normR (-1.),
     rhoNorm (-1),
     precondPtr (&identity),
     UP (SxXC::UP), DN (SxXC::DOWN)
{
   rhoMixing = identity.a = rhoMixingIn;
}


SxRhoMixerG::~SxRhoMixerG ()
{
   // empty
}


void SxRhoMixerG::setType (MixerType t)
{
   type = t;
}

void SxRhoMixerG::setNorm (double rhoNormIn, double dOmegaIn)
{
   SX_CHECK (rhoNormIn > 0., rhoNormIn);
   SX_CHECK (dOmegaIn  > 0., dOmegaIn);
   rhoNorm = rhoNormIn;
   dOmega  = dOmegaIn;
}

void SxRhoMixerG::setNormModus (NormHandling modus)
{
   renormModus = modus;
}


void SxRhoMixerG::setPreconditioner (const SxPreconditioner &preconditioner)
{
   precondPtr = &preconditioner;
}


void SxRhoMixerG::setMaxSteps (int m)
{
   maxSteps = m;
}


void SxRhoMixerG::addRhoIn (const SxVector<Double> &in)
{
   RhoR rhoR(1);
   rhoR(0).resize (in.getSize());
   for (int i=0; i < in.getSize(); ++i)  rhoR(0)(i) = in(i);
   addRhoIn (rhoR);
}

void SxRhoMixerG::addRhoIn (const RhoR &in)
{
   SX_CHECK (rhoNorm > 0. || renormModus == RenormOff,
             rhoNorm, renormModus);
               // Set NormHandling::RenormOff or call setRhoNorm first!

   ssize_t iSpin, nSpin = in.getSize();
   RhoR rhoR(nSpin);
   for (iSpin=0; iSpin < nSpin; iSpin++)  rhoR(iSpin).copy (in(iSpin));

   SX_CHECK (nSpin > 0);
   const SxRBasis *rBasis 
      = dynamic_cast<const SxRBasis *> (in(0).getBasisPtr ());
   SX_CHECK (rBasis);
   const SxGBasis &G = rBasis->getGBasis ();
   RhoG rhoG(nSpin); 
   for (iSpin = 0; iSpin < nSpin; iSpin++)  {
      rhoG(iSpin) = rhoR(iSpin).to (G);
   }
   
   rhoIn.append (rhoG);
}

void SxRhoMixerG::addRhoInG (const RhoG &rhoG)
{
   SX_CHECK (rhoNorm > 0. || renormModus == RenormOff,
             rhoNorm, renormModus);
               // Set NormHandling::RenormOff or call setRhoNorm first!

   RhoG rhoCopy(rhoG.getSize ());
   for (int iSpin = 0; iSpin < rhoG.getSize (); ++iSpin)
      rhoCopy(iSpin).copy (rhoG(iSpin));
   rhoIn.append (rhoCopy);
}




void SxRhoMixerG::addRhoOut (const SxVector<Double> &out)
{
   RhoR rhoR(1);
   rhoR(0).resize (out.getSize());
   for (int i=0; i < out.getSize(); ++i)  rhoR(0)(i) = out(i);
   addRhoOut (rhoR);
}


void SxRhoMixerG::addRhoOut (const RhoR &out)
{
   ssize_t iSpin, nSpin = out.getSize();
   ssize_t m            = rhoIn.getSize() - 1;
   RhoR rhoR(nSpin);
   for (iSpin=0; iSpin < nSpin; iSpin++)  rhoR(iSpin).copy (out(iSpin));
   
   SX_CHECK (nSpin > 0);
   const SxRBasis *rBasis 
      = dynamic_cast<const SxRBasis *> (out(0).getBasisPtr ());
   SX_CHECK (rBasis);
   const SxGBasis &G = rBasis->getGBasis ();
   RhoG rhoG(nSpin); 
   for (iSpin = 0; iSpin < nSpin; iSpin++)  {
      rhoG(iSpin) = rhoR(iSpin).to (G);
   }
   rhoOut << rhoG;

   // --- update residual vector
   R << ( rhoOut(m) - rhoIn(m) );  // ref2 (78), spinpol.
   
   // --- remove G=0 component
   for (iSpin = 0; iSpin < nSpin; iSpin++)  
      R(m)(iSpin)(0) = 0.;
   
   if (m > 0)  {                             
      dRhoIn << ( rhoIn(m) - rhoIn(m-1) );
      dR << ( R(m) - R(m-1) );            
      // --- remove G=0 component
      for (iSpin = 0; iSpin < nSpin; iSpin++)  {
         dRhoIn(m-1)(iSpin)(0) = 0.;
         dR(m-1)(iSpin)(0)     = 0.;
      }
   }

   SX_CHECK (rhoIn.getSize() == rhoOut.getSize(),
             rhoIn.getSize(),   rhoOut.getSize());
}


void SxRhoMixerG::addResidue (const RhoR &rIn)
{
   SX_CHECK (rIn.getSize() == 1,  rIn.getSize()); // spinpol not implemented
   SX_CHECK (rhoIn.getSize() > 0, rhoIn.getSize());
   ssize_t iSpin, nSpin = rhoIn(0).getSize();
   ssize_t m            = rhoIn.getSize() - 1;
   
   const SxRBasis *rBasis 
      = dynamic_cast<const SxRBasis *> (rIn(0).getBasisPtr ());
   SX_CHECK (rBasis);
   const SxGBasis &G = rBasis->getGBasis ();
   // --- change to G space
   RhoG rG(nSpin); 
   for (iSpin = 0; iSpin < nSpin; iSpin++)  {
      rG(iSpin) = rIn(iSpin).to (G);
   }
   R << rG;

   // --- remove G=0 component
   for (iSpin = 0; iSpin < nSpin; iSpin++)  R(m)(iSpin)(0) = 0.;

   // --- update of residual vector - in case you don't use ::addResidue,
   //     this is done in ::addRhoOut => never use both routines at time!!!
   if (m > 0)  {
      dRhoIn << ( rhoIn(m) - rhoIn(m-1) );
      dR     << ( R(m) - R(m-1) );
      // --- remove G=0 component
      for (iSpin = 0; iSpin < nSpin; iSpin++)  {
         dRhoIn(m-1)(iSpin)(0) = 0.;
         dR(m-1)(iSpin)(0)     = 0.;
      }
   }

//   // --- update of rhoOut
//   RhoG rhoOutTmp(nSpin);
//   for (iSpin=0; iSpin < nSpin; iSpin++)  
//      rhoOutTmp(iSpin).copy (rG(iSpin) + rhoIn(m)(iSpin));
//   rhoOut << rhoOutTmp;

//   ^^^ Das geht so einfach nicht, weil in ::addRhoOut, was man dann be-
//       nutzen müsste, das Residuum R nochmal einen drangehangen bekäme.

   SX_CHECK (rhoIn.getSize() == R.getSize(),
             rhoIn.getSize(),   R.getSize());
}

void SxRhoMixerG::addResidueG (const RhoG &rG)
{
   SX_CHECK (rG.getSize() == 1,  rG.getSize()); // spinpol not implemented
   SX_CHECK (rhoIn.getSize() > 0, rhoIn.getSize());
   //ssize_t iSpin, nSpin = rhoIn(0).getSize();
   ssize_t m            = rhoIn.getSize() - 1;
   
   R << rG;

   // --- remove G=0 component
   //for (iSpin = 0; iSpin < nSpin; iSpin++)  R(m)(iSpin)(0) = 0.;

   // --- update of residual vector - in case you don't use ::addResidue,
   //     this is done in ::addRhoOut => never use both routines at time!!!
   if (m > 0)  {
      dRhoIn << ( rhoIn(m) - rhoIn(m-1) );
      dR     << ( R(m) - R(m-1) );
      // --- remove G=0 component
      /*
      for (iSpin = 0; iSpin < nSpin; iSpin++)  {
         dRhoIn(m-1)(iSpin)(0) = 0.;
         dR(m-1)(iSpin)(0)     = 0.;
      }
      */
   }

   SX_CHECK (rhoIn.getSize() == R.getSize(),
             rhoIn.getSize(),   R.getSize());
}


RhoG SxRhoMixerG::getRhoIn (int i) const
{
   SX_CHECK (i >= 0 && i < rhoIn.getSize(), i, rhoIn.getSize());

   return rhoIn(i);
}




RhoR SxRhoMixerG::getMixedRho () const
{
   SX_CHECK   (rhoIn.getSize()  >= 1, rhoIn.getSize());
   SX_CHECK  (rhoOut.getSize() >= 1 || renormModus == RenormOff,
              rhoOut.getSize(), renormModus);
   SX_CHECK (rhoIn.getSize() == rhoOut.getSize() ||
             renormModus == RenormOff,
             rhoIn.getSize(), rhoOut.getSize(), renormModus);
   SX_CHECK  (rhoIn.getSize() == R.getSize(),
              rhoIn.getSize(),   R.getSize());

   ssize_t iSpin, nSpin = rhoIn(0).getSize();
   ssize_t m            = rhoIn.getSize() - 1;

   RhoG rhoOptG;
   if (type == Linear)  rhoOptG = getMixedRhoLinear ();
   else                 rhoOptG = getMixedRhoPulay ();
   RhoR rhoOpt(nSpin);
   const SxGBasis *gBasis 
      = dynamic_cast<const SxGBasis *> (rhoOptG(0).getBasisPtr ());
   SX_CHECK (gBasis);
   for (iSpin = 0; iSpin < nSpin; iSpin++)  {
      rhoOpt(iSpin) = rhoOptG(iSpin).to (gBasis->getRBasis ());
   }

   if (renormModus == SxRhoMixerG::RenormOn)  {
      // --- get current norm
      double rhoOptNorm = 0.;
      for (iSpin=0; iSpin < nSpin; iSpin++)  
         rhoOptNorm += rhoOpt(iSpin).sum() * dOmega;

      // --- renormalize rhoOpt
      for (iSpin=0; iSpin < nSpin; iSpin++)
         rhoOpt(iSpin) *= rhoNorm/rhoOptNorm;
   }

   // --- compute norm of residue in G-space !
   normR = sqrt ( double((R(m)(UP) ^ R(m)(UP)).chop()) );   // ref2 (79)
   if (nSpin > 2)
      normRPol = sqrt ( double((R(m)(DN) ^ R(m)(DN)).chop()) );   // ref2 (79)

   // --- clean up
   removeFirst ();

   return rhoOpt;
}


RhoG SxRhoMixerG::getMixedRhoG () const
{
   SX_CHECK   (rhoIn.getSize()  >= 1, rhoIn.getSize());
   SX_CHECK  (rhoOut.getSize() >= 1 || renormModus == RenormOff,
              rhoOut.getSize(), renormModus);
   SX_CHECK (rhoIn.getSize() == rhoOut.getSize() ||
             renormModus == RenormOff,
             rhoIn.getSize(), rhoOut.getSize(), renormModus);
   SX_CHECK  (rhoIn.getSize() == R.getSize(),
              rhoIn.getSize(),   R.getSize());

   ssize_t nSpin = rhoIn(0).getSize();
   ssize_t m            = rhoIn.getSize() - 1;

   RhoG rhoOptG;
   if (type == Linear)  rhoOptG = getMixedRhoLinear ();
   else                 rhoOptG = getMixedRhoPulay ();
 
   // --- compute norm of residue in G-space !
   normR = sqrt ( double((R(m)(UP) ^ R(m)(UP)).chop()) );   // ref2 (79)
   if (nSpin > 1)
      normRPol = sqrt ( double((R(m)(DN) ^ R(m)(DN)).chop()) );   // ref2 (79)

   // --- clean up
   removeFirst ();

   return rhoOptG;
}


RhoG SxRhoMixerG::getMixedRhoLinear () const
{
   ssize_t nSpin = rhoIn(0).getSize(), m = rhoIn.getSize()-1, iSpin;
   //bool spin = (nSpin > 1);
   RhoG rhoOpt (nSpin);

   double a = rhoMixing;
   const SxPreconditioner &K = *precondPtr;

   SxMeshR rhoMix, rhoSpin;

   if (precondPtr == &identity)  {  // TODO: ugly!!!
         
         for (iSpin = 0; iSpin < nSpin; iSpin++)  {
            rhoOpt(iSpin) = rhoIn(m)(iSpin) + a*(R(m)(iSpin));
//            rhoOpt(iSpin)(0) = spinMixing * rhoOut(m)(iSpin)(0) + 
//                               (1. - spinMixing) * rhoIn(m)(iSpin)(0);
         }
      }  else  {
         
         for (iSpin = 0; iSpin < nSpin; iSpin++)  {
            rhoOpt(iSpin) = rhoIn(m)(iSpin) + K*(R(m)(iSpin));
//            rhoOpt(iSpin)(0) = spinMixing * rhoOut(m)(iSpin)(0) + 
//                               (1. - spinMixing) * rhoIn(m)(iSpin)(0);
         }
   }
      
   return rhoOpt;
}


RhoG SxRhoMixerG::getMixedRhoPulay () const
{
   ssize_t i, j, iSpin, nSpin = rhoIn(0).getSize(), m = rhoIn.getSize()-1;
   bool spin = (nSpin > 1);
   RhoG rhoOpt (nSpin);
   Real8 a = rhoMixing;
   SxMeshR rhoSpin;

   SX_CHECK (rhoIn.getSize() == rhoOut.getSize() ||
             renormModus == RenormOff,
             rhoIn.getSize(), rhoOut.getSize(), renormModus);
   SX_CHECK  (rhoIn.getSize() == R.getSize(), rhoIn.getSize(), R.getSize());
   SX_CHECK  (rhoIn.getSize() == dRhoIn.getSize()+1,
              rhoIn.getSize(),   dRhoIn.getSize());

   const SxPreconditioner &K = *precondPtr;

   // --- 1st step: linear mixing
   if (m == 0)  {

      rhoOpt = getMixedRhoLinear ();

   }  else  {
   
      SX_CHECK (dR.getSize() == R.getSize()-1, dR.getSize(), R.getSize());

      // --- compute Pulay matrix
      SxSymMatrix<TPrecRhoR> A(m);
      SxVector<TPrecRhoR>    B(m), alpha;
      for (i=0; i < m; i++)  {
         for (j=i; j < m; j++)  {
            A(i,j) = (dR(i)(UP) ^ dR(j)(UP)).chop();
            if (spin) A(i,j) = A(i,j) + (dR(i)(DN) ^ dR(j)(DN)).chop();   
         }
         B(i) = (dR(i)(UP) ^ R(m)(UP)).chop();
         if (spin) B(i) = B(i) + (dR(i)(DN) ^ R(m)(DN)).chop();        
      }
      // TODO: expand is to be removed
      alpha = -(A.inverse().expand() ^ B);     // ref2, (90)

      //cout << "alpha = " << alpha << endl;
      //cout << "A = " << A.expand() << endl;
      //cout << "B = " << B << endl;

      // --- compute new optimal input charge density
      if (precondPtr == &identity)  {  // TODO: ugly!!!
         for (iSpin=0; iSpin < nSpin; iSpin++)  {
            rhoOpt(iSpin) = rhoIn(m)(iSpin) +  a*(R(m)(iSpin));
            for (i=0; i < m; i++)  {
               rhoOpt(iSpin) += alpha(i) * (dRhoIn(i)(iSpin) 
                                             +  a*(dR(i)(iSpin)));
            }
//            rhoOpt(iSpin)(0) = spinMixing * rhoOut(m)(iSpin)(0) + 
//                            (1. - spinMixing) * rhoIn(m)(iSpin)(0);
         }

      } else  {
         for (iSpin=0; iSpin < nSpin; iSpin++)  {
            rhoOpt(iSpin) = rhoIn(m)(iSpin) +  K*(R(m)(iSpin));
            for (i=0; i < m; i++)  {
               rhoOpt(iSpin) += alpha(i) * (dRhoIn(i)(iSpin) 
                                             +  K*(dR(i)(iSpin)));
            }
//            rhoOpt(iSpin)(0) = spinMixing * rhoOut(m)(iSpin)(0) + 
//                            (1. - spinMixing) * rhoIn(m)(iSpin)(0);
         }
      }
   }
   return rhoOpt;
}


void SxRhoMixerG::removeFirst () const
{
   if (rhoIn.getSize()  > maxSteps)    rhoIn.removeFirst ();
   if (rhoOut.getSize() > maxSteps)    rhoOut.removeFirst ();
   if (R.getSize()      > maxSteps)    R.removeFirst ();

   if (0 < dRhoIn.getSize() && dRhoIn.getSize() > maxSteps-1)  
      dRhoIn.removeFirst ();
   if (0 < dR.getSize()     && dR.getSize()     > maxSteps-1) 
      dR.removeFirst ();
}

// --------------------------------------------------------------------------
// Kerker preconditioning
// --------------------------------------------------------------------------
SxRhoMixerG::SxKerker::SxKerker (const SxGBasis &G, double a, double q0)
   : SxPreconditioner ()
{
   if (fabs(q0) < 1e-10 )  q0 = 1e-10;  // avoid singularity
   
   K = a * G.g2 / (G.g2 + q0*q0);  // ref2, (82)
   K(0) = 0.; // ensure that K*R[m] contains no norm
}

SxMeshG SxRhoMixerG::SxKerker::operator* (const SxMeshG &meshG) const
{
   return ( K * meshG );  // apply preconditioner in |G> space
}

