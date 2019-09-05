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

#include <SxXC.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <SxProjector.h>
#include <SxRadBasis.h>
#include <SxConstants.h>
#include <SxVector3.h>
#include <SxHamiltonian.h>
#include <SxFSAction.h>
#include <SxXCFunctional.h>
#include <SxLoopMPI.h>  // LoopMPI


//------------------------------------------------------------------------------
// Bibliography:
//
// ref1 - PRB 23, 5048 (1981)
// ref2 - PRB 54, 16533 (1996)
// ref3 - K. Burke's original code, in src/api/k.burkes.f
// ref4 - PRB 33, 8800 (1986)
// ref5 - PRB 45, 13244 (1992)
// ref6 - PRL 77, 3865 (1996)
// ref7 - "Electronic Structure of Solids", 11-20 (1991), 
//        ed. P.Ziesche, H.Eschrig (Akademie Verlag, Berlin)
// ref8 - http://www.sphinxlib.de/documentation/lectures/pbe_white_bird.{ps,pdf}
// ref9 - PRB 59, 10031 (1999)
// ref10 - B. Grabowski, diploma thesis
//------------------------------------------------------------------------------


SxXC::SxXC (const SxXC &in)
   : xFunctional(in.xFunctional),
     cFunctional(in.cFunctional),
     xcFunctional(in.xcFunctional),
     eX(0.), eC(0.), eXc(0.),
     eExchange(0.), eCorrelation (0.), eVxc(0.),
     vXcFile(in.vXcFile),
     nlcc (in.nlcc),
     phiCore (in.phiCore),
     nSpin (in.nSpin),
     calcX(in.calcX),
     calcC(in.calcC)
{
   vXc.resize (nSpin);
   vEx.resize (nSpin);
   vCor.resize (nSpin);
   for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
      if (in.vXc (iSpin).getSize () > 0) vXc (iSpin).copy (in.vXc(iSpin));
      if (in.vEx (iSpin).getSize () > 0) vEx (iSpin).copy (in.vEx(iSpin));
      if (in.vCor(iSpin).getSize () > 0) vCor(iSpin).copy (in.vCor(iSpin));
   }
   rhoCoreR = in.rhoCoreR;
   rhoXcR.resize (nSpin);
}

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
SxXC::SxXC (const SxPseudoPot &pot, const SxGBasis *gBasisPtr, int nSpin_)
   : calcX(true), calcC(true)
{
   nSpin = nSpin_;
   init ();
   if (pot.nlcc)  computePhiCore (pot, gBasisPtr);
}

SxXC::SxXC (int nSpin_)
   : calcX(true), calcC(true)
{
   nSpin = nSpin_;
   init ();
}

void SxXC::init ()
{
   nlcc = false;
   vXc.resize    (nSpin);
   vEx.resize    (nSpin);
   vCor.resize   (nSpin);
   rhoXcR.resize (nSpin);
}

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
SxXC::~SxXC () 
{
   // empty
}

void SxXC::read (const SxSymbolTable *hamiltonian, const SxCell &cell)
{
   xFunctional = cFunctional = xcFunctional = getXCFunctional (hamiltonian); 
   if (hamiltonian->containsGroup ("xcMesh"))  {
      setXcMesh (hamiltonian->getGroup("xcMesh"), cell);
      origXcBasis = xcBasisPtr;
   }
} 

void SxXC::setXcMesh (const SxSymbolTable *meshGrp, const SxCell &cell)
{
   try  {
      SxMesh3D mesh;
      if (meshGrp->contains ("mesh"))  {
         mesh = SxVector3<Int> (meshGrp->get("mesh")->toIntList ());
      } else if (meshGrp->contains ("eCut"))  {
         double eCut = meshGrp->get ("eCut")->toReal ();
         mesh = SxGBasis::getMeshSize (eCut, cell, 1.);
      } else if (meshGrp->contains ("meshAccuracy"))  {
         double accuracy = meshGrp->get ("meshAccuracy")->toReal ();
         double eCut     = SxGBasis::getECut (meshGrp->topLevel ());
         mesh = SxGBasis::getMeshSize (eCut, cell, accuracy);
      }
      if (origXcBasis && origXcBasis->fft3d.mesh == mesh)
         xcBasisPtr = origXcBasis;
      else if (!(xcBasisPtr && xcBasisPtr->fft3d.mesh == mesh)) {
         unsigned oldMode = 0;
         SxFFT::quickFFTPlanner (SxFFT::Estimate, &oldMode);
         xcBasisPtr = SxPtr<SxRBasis>::create (mesh(0), mesh(1), mesh(2), cell);
         SxFFT::restorePlannerMode (oldMode);
      } else
         {} // mesh has not changed
      cout << "xc will be computed on " << mesh << " mesh.\n";
   }  catch (const SxException &e)  {
      e.print ();
      SX_EXIT;
   }
}

SxXC::XCFunctional SxXC::getXCFunctional (const SxSymbolTable *top)
{
   XCFunctional xcFunctional = LDA;
   try  {
      const SxSymbolTable *hamiltonian;
      if (top->contains("xc"))
         hamiltonian = top;
      else
         hamiltonian = SxHamiltonian::getHamiltonianGroup (top->topLevel ());

      if (hamiltonian->contains("xc"))  {
         xcFunctional = (XCFunctional)hamiltonian->get("xc")->toInt();
      }
      if (isHybrid(xcFunctional))  {
         if (xcFunctional == SxXC::HSE06)  {
            if (hamiltonian->contains("omegaHSE"))
               SxXCFunctional::omegaHSE = hamiltonian->get("omegaHSE")
                                          ->toReal ();
            else
               SxXCFunctional::omegaHSE = 0.11;
         } else if (xcFunctional == SxXC::PBE0) {
            SxXCFunctional::omegaHSE = 0.;
         } else {
            cout << "Unknown hybrid functional" << xcFunctional << endl;
            SX_EXIT;
         }
         if (hamiltonian->contains ("alphaHybrid"))  {
            SxXCFunctional::alphaHybrid = hamiltonian->get ("alphaHybrid")
                                          ->toReal ();
            cout << "Setting alphaHybrid to " << SxXCFunctional::alphaHybrid
                 << endl;
         }
      }
   }  catch (const SxException &e)  {
      e.print ();
      SX_QUIT;
   }
   return xcFunctional;
}

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void SxXC::computePhiCore (const SxPseudoPot &pot, const SxGBasis *gBasisPtr)
{
   if (pot.nlcc)  {
      SX_CHECK (gBasisPtr);
      SX_CHECK (gBasisPtr->structPtr);
      double omega = gBasisPtr->structPtr->cell.volume;
      SX_CHECK (omega > 0.);
      int nSpecies = pot.getNSpecies ();
      nlcc = true;
      SxRadBasis radialBasis (pot.rad, pot.logDr);
      SxDiracVec<TReal8>    rad, radCore;
      Real8                 logDr;
      phiCore.resize(nSpecies);

      for (int is=0; is < nSpecies; is++)  {
         
         // --- is species core corrected?
         if (pot.radCore(is)(0).getSize() > 0)  {
            rad      = SxDiracVec<TReal8> (pot.radCore(is)(0));
            radCore  = SxDiracVec<TReal8> (pot.radCore(is)(1));
            logDr    = pot.logDr(is);  
            //TODO: radCore should have its own logDr!!!
         
            // --- Ylm(0): single atoms are spherical (s-like)
            //cout << "Core density projection should be done Dirac-like" << endl;
            // phiCore(is) = (*gBasisPtr) | radCore;
            phiCore(is) = 
               radialBasis.toPWBasis(rad, radCore, *gBasisPtr, 0, logDr)
               / sqrt(omega); // see core normalization in SxPseudoPot
         }  else  { // not core corrected
            phiCore(is) = SxDiracVec<TPrecPhi> ();
         }
      }
   }
}


//==============================================================================
// LDA routines
//==============================================================================
void SxXC::computeXC () 
{
   if (nlcc)  {

      const SxGBasis *pwBasis = NULL;
      int nSpecies = int(phiCore.getSize ());
      SxMeshG rhoCoreG;
    
      for (int is=0; is < nSpecies; is++)  {
         if (phiCore(is).getSize () > 0)  {
            if (!pwBasis)  {
               // --- get G basis and resize rhoCoreG
               pwBasis = dynamic_cast<const SxGBasis *>
                         (phiCore(is).getBasisPtr ());
               SX_CHECK (pwBasis);
               rhoCoreG.resize (pwBasis->ng);
               rhoCoreG.set (0.);
               rhoCoreG.setBasis(pwBasis);
            }
            rhoCoreG += pwBasis->structureFactors(is) * phiCore(is);
         }
      }
      const SxRBasis &R = pwBasis->getRBasis ();
      rhoCoreR = rhoCoreG.to (R);
    
      sxprintf ("NLCC: %g core electrons\n", rhoCoreR.sum() * R.dOmega);
      if (nSpin > 1)  rhoCoreR /= (PrecRhoR)nSpin;

      // --- for GGAs: rhoCoreR must not be zero!
      SxMeshR::Iterator it;
      for (it = rhoCoreR.begin(); it != rhoCoreR.end(); it++)
         if ( *it < 1e-12 )  *it = 1e-12;
   }
}


enum XCTimer { updateXC1, updateXC2, updateXC3 };
SX_REGISTER_TIMERS (XCTimer)
{
   regTimer (updateXC1, "updateXC1");
   regTimer (updateXC2, "updateXC2");
   regTimer (updateXC3, "updateXC3");
}


void SxXC::updateXC (const RhoR &rhoR, const SxPWSet *,
                     const SxFermi *)
{
   const SxRBasis *rBasisPtr 
      = dynamic_cast<const SxRBasis *>(rhoR(0).getBasisPtr());
   SX_CHECK (rBasisPtr);

   int iSpin;

   // --- add core charge density in case of NLCC
   for (iSpin=0; iSpin < nSpin; iSpin++)  {
      if (nlcc)  rhoXcR(iSpin) = rhoR(iSpin) + rhoCoreR;
      else       rhoXcR(iSpin) = rhoR(iSpin);
   }


   if (xcBasisPtr && (xcBasisPtr.getPtr () != rBasisPtr))  {
      //
      SX_START_TIMER (updateXC1);
      const SxGBasis &G = rBasisPtr->getGBasis ();
      // takes time depending on the mode of FFTW initialization
      {
         unsigned oldMode = 0;
         SxFFT::quickFFTPlanner (SxFFT::Estimate, &oldMode);
         xcBasisPtr->registerGBasis (G);
         SxFFT::restorePlannerMode (oldMode);
      }
      SX_STOP_TIMER (updateXC1);
      //
      SX_START_TIMER (updateXC2);
      for (iSpin=0; iSpin < nSpin; iSpin++)  {
         // Fourier interpolation
         rhoXcR(iSpin) = rhoXcR(iSpin).to(G).to(*xcBasisPtr);
         vXc(iSpin) = SxDiracVec<Double> (*xcBasisPtr);
      }
      SX_STOP_TIMER (updateXC2);
   }



   // --- XC functional
   SX_START_TIMER (updateXC3);
   switch (xcFunctional)  {
      case LDA        : if (nSpin == 1) computeLDA  ();
                        else            computeLSDA ();
                        break;
      case LDA_PW     : if (nSpin == 1) computeLDA_PW  ();
                        else            computeLSDA_PW ();
                        break;                  
      case PBE        : // khr -- improve parallelization
                        updateGGA_PBE_WB ();
                        break;
      case PBE_LDA    : {
                           SxArray<SxMeshR> tmpVc(nSpin);
                           PrecEnergy eCorTmp;
                           bool calcXo = calcX, calcCo = calcC;

                           // L(S)DA correlation
                           calcX = false;
                           if (nSpin == 1) computeLDA  (); 
                           else            computeLSDA ();
                    
                           for (iSpin=0; iSpin < nSpin; iSpin++) 
                              tmpVc(iSpin).copy(vXc(iSpin));
                           eCorTmp = eXc;
                           calcX = calcXo;
                        
                           // PBE exchange
                           calcC = false;
                           updateGGA_PBE_WB ();
                           calcC = calcCo;
                      
                           // add them together
                           for (iSpin=0; iSpin < nSpin; iSpin++) 
                              vXc(iSpin)+=tmpVc(iSpin);
                           eXc += eCorTmp;
                        }
                        break;
      case READ_VXC   : readVXC (*rBasisPtr);
                        break;
      case HSE06      :
      case PBE0       : updateGGA_WB ();
                        break;
      default         : sxprintf ("Unknown XC functional.\n");
                        SX_EXIT;
   }
   SX_STOP_TIMER (updateXC3);


   // --- back-transform to rBasis if xcBasis was used
   if (xcBasisPtr && (xcBasisPtr.getPtr () != rBasisPtr))  {
      const SxGBasis &G = rBasisPtr->getGBasis ();
      for (iSpin=0; iSpin < nSpin; iSpin++)  {
         // Fourier interpolation
         vXc(iSpin) = vXc(iSpin).to (G).to (*rBasisPtr); 
      }
   }

}

void SxXC::readVXC (const SxRBasis &R) 
{
   // reading extern exchange correlation potential for debugging
   cout << "Reading vXC from file \"" << vXcFile << "\"." << endl;

   VxcR vXcExt;
   try  {
      SxBinIO io (vXcFile, SxBinIO::BINARY_READ_ONLY);
      
      // --- read and check nSpin
      int nSpinExt = io.getDimension ("nMeshes");
      if (nSpinExt != nSpin && nSpin <= 0)
         throw SxException ("External vXC has wrong spin.", __FILE__, __LINE__);
      nSpin = nSpinExt;
      vXc.resize (nSpin);
      
      SxMatrix3<Double> cellExt;
      SxVector3<Int>    dimExt;

      // --- read in cell, mesh, and mesh dim
      vXcExt = io.readMesh (&cellExt, &dimExt);
      // check cell
      if ((R.cell - cellExt).absSqr ().sum() > 1e-10)
         throw SxException ("vXC cell is wrong.", __FILE__, __LINE__);
      // check mesh dim
      if ((R.fft3d.mesh - dimExt).normSqr () > 0)
         throw SxException ("vXC mesh dimension is wrong.", __FILE__, __LINE__);
   }  catch (const SxException &e)  {
      e.print ();
      cout << "Couldn't read xc potential, exiting ..." << endl;
      SX_QUIT;
   }

   for (int iSpin = 0; iSpin < nSpin; ++iSpin)  vXc(iSpin) = vXcExt(iSpin);
}  

//==============================================================================
// LDA - Perdew/Zunger '81
//==============================================================================
void SxXC::computeLDA ()
{
   const SxRBasis *rBasisPtr = dynamic_cast<const SxRBasis *>
                              (rhoXcR(0).getBasisPtr ());
   SX_CHECK (rBasisPtr);
   int n = rBasisPtr->getMeshSize();
   SX_CHECK (nSpin == 1, nSpin);

   // --- set up xc functional "machine"
   int mode = SxXCFunctional::ComputeEandV;
   if (calcX) mode |= SxXCFunctional::ComputeX;
   if (calcC) mode |= SxXCFunctional::ComputeC;
   SxXCFunctional xc (nSpin, mode);

   // --- get alias for density and potential
   const SxMeshR &rho = rhoXcR(0);
   SxMeshR       &vxc = vXc(0);
   if (vxc.getSize () != n) vxc = SxMeshR (rBasisPtr);

   // --- compute xc energy and potential
   double eXc_ = 0.;
#ifdef USE_OPENMP
#pragma omp parallel for firstprivate(xc) reduction(+:eXc_)
#endif
   for (int i = 0; i < n; ++i)  {
      xc.computeLDA (rho(i));
      vxc(i) = xc.vXc(0);
      eXc_  += xc.eXc;
   }
   eXc = eXc_ * rBasisPtr->dOmega;

   // --- setup eX/vEx and eC/vCor if needed
   if (!calcX)  {
      eCorrelation = eC = eXc;
      vCor(UP) = vxc;
   }
   if (!calcC)  {
      eExchange = eX = eXc;
      vEx(UP) = vxc;
   }
}

void SxXC::computeLSDA ()
{
   const SxRBasis *rBasisPtr = dynamic_cast<const SxRBasis *>
                              (rhoXcR(0).getBasisPtr ());
   int n = rBasisPtr->getMeshSize();
   SX_CHECK (rBasisPtr);
   SX_CHECK (nSpin == 2, nSpin);

   // --- set up xc functional "machine"
   int mode = SxXCFunctional::ComputeEandV;
   if (calcX) mode |= SxXCFunctional::ComputeX;
   if (calcC) mode |= SxXCFunctional::ComputeC;
   SxXCFunctional xc (nSpin, mode);

   // --- get alias for densities and potentials
   const SxMeshR &rhoUp   = rhoXcR(0), 
                 &rhoDown = rhoXcR(1);
   SxMeshR       &vxcUp   = vXc(0),
                 &vxcDown = vXc(1);
   if (vxcUp.getSize ()   != n) vxcUp   = SxMeshR (rBasisPtr);
   if (vxcDown.getSize () != n) vxcDown = SxMeshR (rBasisPtr);

   // --- compute xc energy and potentials
   double eXc_ = 0.; // openMP cannot make class variables private
#ifdef USE_OPENMP
#pragma omp parallel for firstprivate(xc) reduction(+:eXc_)
#endif
   for (int i = 0; i < n; ++i)  {
      xc.computeLSDA (rhoUp(i),rhoDown(i));
      vxcUp(i)   = xc.vXc(0);
      vxcDown(i) = xc.vXc(1);
      eXc_   += xc.eXc;
   }
   eXc = eXc_ * rBasisPtr->dOmega;

   // --- setup eX/vEx and eC/vCor if needed
   if (!calcX)  {
      eCorrelation = eC = eXc;
      vCor(UP)   = vxcUp;
      vCor(DOWN) = vxcDown;
   }
   if (!calcC)  {
      eExchange = eX = eXc;
      vEx(UP)   = vxcUp;
      vEx(DOWN) = vxcDown;
   }
}

//==============================================================================
// LDA - Perdev/Wang '91
//==============================================================================
void SxXC::computeLDA_PW ()
{
   const SxRBasis *rBasisPtr = dynamic_cast<const SxRBasis *>
                              (rhoXcR(0).getBasisPtr ());
   SX_CHECK (rBasisPtr);
   int n = rBasisPtr->getMeshSize();
   SX_CHECK (nSpin == 1, nSpin);

   // --- set up xc functional "machine"
   int mode = SxXCFunctional::ComputeEandV;
   if (calcX) mode |= SxXCFunctional::ComputeX;
   if (calcC) mode |= SxXCFunctional::ComputeC;
   SxXCFunctional xc (nSpin, mode);

   // --- get alias for density and potential
   const SxMeshR &rho = rhoXcR(0);
   SxMeshR       &vxc = vXc(0);
   if (vxc.getSize () != n) vxc = SxMeshR (rBasisPtr);

   // --- compute xc energy and potential
   double eXc_ = 0.; // openMP cannot make class variables private
#ifdef USE_OPENMP
#pragma omp parallel for firstprivate(xc) reduction(+:eXc_)
#endif
   for (int i = 0; i < n; ++i)  {
      double rhoHalf = 0.5 * rho(i);
      xc.computeLSDA_PW (rhoHalf,rhoHalf);
      vxc(i) = xc.vXc(0);
      eXc_   += xc.eXc;
   }
   eXc = eXc_ * rBasisPtr->dOmega;

   // --- setup eX/vEx and eC/vCor if needed
   if (!calcX)  {
      eCorrelation = eC = eXc;
      vCor(UP) = vxc;
   }
   if (!calcC)  {
      eExchange = eX = eXc;
      vEx(UP) = vxc;
   }
}

void SxXC::computeLSDA_PW ()
{
   const SxRBasis *rBasisPtr = dynamic_cast<const SxRBasis *>
                              (rhoXcR(0).getBasisPtr ());
   int n = rBasisPtr->getMeshSize();
   SX_CHECK (rBasisPtr);
   SX_CHECK (nSpin == 2, nSpin);

   // --- set up xc functional "machine"
   int mode = SxXCFunctional::ComputeEandV;
   if (calcX) mode |= SxXCFunctional::ComputeX;
   if (calcC) mode |= SxXCFunctional::ComputeC;
   SxXCFunctional xc (nSpin, mode);

   // --- get alias for densities and potentials
   const SxMeshR &rhoUp   = rhoXcR(0), 
                 &rhoDown = rhoXcR(1);
   SxMeshR       &vxcUp   = vXc(0),
                 &vxcDown = vXc(1);
   if (vxcUp.getSize ()   != n) vxcUp   = SxMeshR (rBasisPtr);
   if (vxcDown.getSize () != n) vxcDown = SxMeshR (rBasisPtr);

   // --- compute xc energy and potentials
   eXc = 0.;
   for (int i = 0; i < n; ++i)  {
      xc.computeLSDA_PW (rhoUp(i),rhoDown(i));
      vxcUp(i)   = xc.vXc(0);
      vxcDown(i) = xc.vXc(1);
      eXc   += xc.eXc;
   }
   eXc *= rBasisPtr->dOmega;

   // --- setup eX/vEx and eC/vCor if needed
   if (!calcX)  {
      eCorrelation = eC = eXc;
      vCor(UP)   = vxcUp;
      vCor(DOWN) = vxcDown;
   }
   if (!calcC)  {
      eExchange = eX = eXc;
      vEx(UP)   = vxcUp;
      vEx(DOWN) = vxcDown;
   }
}


//==============================================================================
// GGA utility routines
//==============================================================================
SxArray<SxMeshR> SxXC::gradient (SxMeshR &meshR)
{
   //SX_EXIT;
   SX_CHECK (dynamic_cast<const SxRBasis*>(meshR.getBasisPtr ()));
   const SxRBasis &R = *dynamic_cast<const SxRBasis*>(meshR.getBasisPtr ());
   const SxGBasis &gBasis = R.getGBasis ();
   SxArray<SxMeshR> gradR(3);  // xyz


// disabled temporarily for compatibility with SxParallelHierarchy

//#ifdef USE_LOOPMPI   // LoopMPI, #defined in <SxLoopMPI.h>
//   {
//      // khr:  Poor man's parallelization over three dimensions.
//      //       TODO: Implement parallelization on a lower level,
//      //             need to talk to CF first.
//      for (int d=0; d < 3; d++)
//      {
//         if (SxLoopMPI::myWork (d)) // LoopMPI
//         {
//            gradR(d) = ( (I * gBasis.gVec.colRef(d)) * meshR.to (gBasis) ).to (R);
//         }
//      }
//      //
//      for (int d=0; d < 3; d++)
//      {
//         int vecLen=0;
//         if (SxLoopMPI::myWork (d)) // LoopMPI
//         {
//            vecLen=(int)gradR(d).getSize();
//         }
//         vecLen=SxLoopMPI::bcast(vecLen, SxLoopMPI::whoseWork (d)); // LoopMPI
//         if (! SxLoopMPI::myWork (d)) // LoopMPI
//         {
//            gradR(d).resize(vecLen);
//         }
//         gradR(d).bcastMPI(SxLoopMPI::whoseWork (d)); // LoopMPI
//      }
//   }
//#else
//   {
      // previous, non parallel case
      for (int d=0; d < 3; d++)
      {
         gradR(d) = ( (I * gBasis.gVec.colRef(d)) * meshR.to (gBasis) ).to (R);
      }
//   }
//#endif

   return gradR;
}


enum PBEWBTimer { updatePBEWB1, updatePBEWB2, updatePBEWB3 };
SX_REGISTER_TIMERS (PBEWBTimer)
{
   regTimer (updatePBEWB1, "updatePBEWB1");
   regTimer (updatePBEWB2, "updatePBEWB2");
   regTimer (updatePBEWB3, "updatePBEWB3");
}


//==============================================================================
// GGA: PBE with White & Bird discretization and for spin polarized systems
//==============================================================================
void SxXC::updateGGA_PBE_WB ()
{

   // TODO: improve parallelization in gradient() to scale beyond 3 cores
   SX_START_TIMER(updatePBEWB1);
   const SxRBasis *rBasisPtr = dynamic_cast<const SxRBasis *>
                              (rhoXcR(0).getBasisPtr ());
   SX_CHECK (rBasisPtr);
   int meshSize = rBasisPtr->getMeshSize();
   // set exchange-correlation, exchange and correlation energy and potential
   eXc = 0.;
   // --- Begin of variable preparation ----------------------------------------
   // prepare meshes for exchange and correlation part
   SX_LOOP(iSpin)
      vXc(iSpin)  = SxMeshR(rBasisPtr);
   vEx.resize(0);
   vCor.resize(0);

   // --- compute gradient, gradRhoSpin needed for exchange part
   //                       gradRhoTotal needed for correlation part
   // NOTE: the gradients are used first for the gradients of the
   //       density and then overwritten on the fly by auxiliary terms for
   //       the White&Bird contribution. This ugly procedure saves a lot
   //       of memory
   SxArray<SxArray<SxMeshR> > gradients(nSpin);
   SX_LOOP(iSpin)
      gradients(iSpin) = gradient (rhoXcR(iSpin));
   SX_STOP_TIMER(updatePBEWB1);
   // --- End of variable preparation ------------------------------------------


   // TODO: improve parallelization
   SX_START_TIMER(updatePBEWB2);

   int xcFlags = SxXCFunctional::ComputeEandV;

   if (calcC) xcFlags |= SxXCFunctional::ComputeC;
   if (calcX) xcFlags |= SxXCFunctional::ComputeX;

   SxXCFunctional xc(nSpin, xcFlags);
   int up = 0, down = nSpin - 1; // nSpin = 1 => up=down

   {
#ifdef USE_LOOPMPI
      // -- calculate offsets
      int nelem = meshSize / SxLoopMPI::nr();
      int i0    = nelem * SxLoopMPI::me();
      int i1 = i0+nelem;
      if ( SxLoopMPI::me() == (SxLoopMPI::nr()-1))
         i1 = meshSize;
#else
      int i0 = 0, i1 = meshSize;
#endif

      // -- work in stripes on the arrays
      double myExc = 0.;
#ifdef USE_OPENMP
#pragma omp parallel for firstprivate(xc) reduction(+:myExc)
#endif
      for (int i=i0; i < i1; i++)
      {
         SxVector3<TPrecRhoR> gradRhoUp, gradRhoDown, gradRhoTl;
         // --- collect gradients
         for (int iDir = 0; iDir < 3; ++iDir)
         {
            gradRhoUp(iDir)   = gradients(up)(iDir)(i);
            gradRhoDown(iDir) = gradients(down)(iDir)(i);
         }
         if (nSpin == 1)
            gradRhoTl = gradRhoUp;
         else
            gradRhoTl = gradRhoUp + gradRhoDown;

         double gradRhoNormUp   = gradRhoUp.norm (),
                gradRhoNormDown = gradRhoDown.norm (),
                gradRhoNormTl   = gradRhoTl.norm ();

         // --- compute PBE functional
         xc.computePBE (rhoXcR(up)(i), gradRhoNormUp,
                        rhoXcR(down)(i), gradRhoNormDown,
                        gradRhoNormTl);

         myExc += xc.eXc;

         // derivative wrt grad rho (total)
         double fc_gradRho = (gradRhoNormTl > 0.)
                           ? (xc.eXc_gradRhoTl / gradRhoNormTl)
                           : 0.;

         for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
            // derivative with respect to rho
            vXc(iSpin)(i) = xc.vXc(iSpin); // 1st term in ref. 10, Eq. A.23b/c
            // derivative wrt grad rho(spin) (note: d|g|/d g_x = g_x / |g|)
            double g = (iSpin == 0) ? gradRhoNormUp : gradRhoNormDown;
            double fx_gradRho = (g > 0.) ? (xc.eXc_grad(iSpin) / g) : 0.;
            // compute R' terms in ref. 10, A.23b/c without FFT factors
            for (int iDir = 0; iDir < 3; ++iDir)  {
               gradients(iSpin)(iDir)(i)
                  = gradients(iSpin)(iDir)(i) * fx_gradRho // ref. 10, A.23b
                  + gradRhoTl(iDir) * fc_gradRho;          // ref. 10, A.23c
            }
         }
      }
      eXc = myExc;
   }
#ifdef USE_LOOPMPI
   {
      // -- sync operations, SxLoopMPI
      eXc = SxLoopMPI::sum(eXc);
      //
      for (int ip=0; ip<SxLoopMPI::nr(); ip++)
      {
         // -- calculate offsets
         int nelem=meshSize/SxLoopMPI::nr();
         int i0=nelem*ip;
         if ( ip == (SxLoopMPI::nr()-1))
         {
            nelem += meshSize%SxLoopMPI::nr();
         }
         // -- do MPI communication, TODO: further optimization!
         for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
            SxLoopMPI::bcast( &vXc(iSpin)(i0),        nelem, ip);
            for (int iDir = 0; iDir < 3; ++iDir)
               SxLoopMPI::bcast( &gradients(iSpin)(iDir)(i0), nelem, ip);
         }
      }

   }
#endif

   eXc *= rBasisPtr->dOmega;
   SX_STOP_TIMER(updatePBEWB2);

   SX_START_TIMER(updatePBEWB3);
   const SxRBasis &R = *rBasisPtr;
   const SxGBasis &G = R.getGBasis ();

   // --- add White & Bird contribution of gradients
   // => FFT [i.e. exp(iG(R-R'))] and iG parts of ref. 10, Eq. A.22, A.23
#ifdef USE_LOOPMPI  // #defined in <SxLoopMPI.h>
   {
      SxArray<PsiG> vXcTmp(nSpin);
      for (int iSpin = 0; iSpin < nSpin; ++iSpin) {
         vXcTmp(iSpin) = PsiG(G);
         vXcTmp(iSpin).set (0.);
         for (int iDir = 0; iDir < 3; ++iDir) {
            if (SxLoopMPI::myWork (iDir + 3 * iSpin)) // LoopMPI
            {
               SxMeshR &fxc_gradR = gradients(iSpin)(iDir);
               PsiG fxc_gradG = G | fxc_gradR;
               fxc_gradR.resize (0); // cut short a memory usage peak
#ifdef USE_OPENMP
#pragma omp parallel for if (G.ng > sxChunkSize)
#endif
               for (ssize_t ig = 0; ig < G.ng; ++ig)
                  vXcTmp(iSpin)(ig) += G.gVec(ig + iDir * G.ng) * fxc_gradG(ig);
            } else {
               gradients(iSpin)(iDir).resize (0); // cut memory usage peak
            }
         }
      }
      // do MPI synchronization after compute loop to MPI-parallelize over spin
      // possible alternative: bring vXc to G space, and FFT later to
      // the normal DFT mesh (which is usually 8x smaller)
      for (int iSpin = 0; iSpin < nSpin; ++iSpin)
      {
         SxLoopMPI::sum (vXcTmp(iSpin));
         vXc(iSpin) += (R | vXcTmp(iSpin)).imag ();// R|(-iG * x)=imag(R|G * x)
      }
   }
#else
   {
      for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
         PsiG fxc_G;
         for (int iDir = 0; iDir < 3; ++iDir)  {
            SxMeshR &fxc_gradR = gradients(iSpin)(iDir);
            PsiG fxc_gradG = G | fxc_gradR;
            fxc_gradR.resize (0); // cut short a memory usage peak
            if (iDir == 0)  {
               fxc_G = fxc_gradG;
#ifdef USE_OPENMP
#pragma omp parallel for if (G.ng > sxChunkSize)
#endif
               for (ssize_t ig = 0; ig < G.ng; ++ig)
                  fxc_G(ig) *= G.gVec(ig + iDir * G.ng);
            } else {
#ifdef USE_OPENMP
#pragma omp parallel for if (G.ng > sxChunkSize)
#endif
               for (ssize_t ig = 0; ig < G.ng; ++ig)
                  fxc_G(ig) += G.gVec(ig + iDir * G.ng) * fxc_gradG(ig);
            }
         }
         vXc(iSpin) += (R | fxc_G).imag ();// R|(-iG * x)=imag(R|G * x)
      }
   }
#endif // LoopMPI
   SX_STOP_TIMER(updatePBEWB3);

}




//==============================================================================
// generic GGA with White & Bird discretization and for spin polarized systems
//==============================================================================
void SxXC::updateGGA_WB ()
{
   const SxRBasis *rBasisPtr = dynamic_cast<const SxRBasis *>
                              (rhoXcR(0).getBasisPtr ());
   SX_CHECK (rBasisPtr);
   int iSpin, meshSize = rBasisPtr->getMeshSize();
   
   // set exchange-correlation, exchange and correlation energy and potential
   eXc = 0.;   
   
   // --- Begin of variable preparation ----------------------------------------
   // prepare meshes for exchange and correlation part
   for (iSpin=0; iSpin < nSpin; iSpin++)  {
      vXc(iSpin)  = SxMeshR(rBasisPtr);
   }
   vEx.resize(0);
   vCor.resize(0);
   SxArray<SxMeshR>   fx_gradRho(nSpin);
   SxMeshR            fc_gradRho(rBasisPtr);
   
   for (iSpin=0; iSpin < nSpin; iSpin++)
      fx_gradRho(iSpin) = SxMeshR(rBasisPtr); 
   
   // --- compute gradient, gradRhoSpin needed for exchange part
   //                       gradRhoTotal needed for correlation part
   SxArray<SxArray<SxMeshR> > gradRhoSpin(nSpin);
   for (iSpin=0; iSpin < nSpin; iSpin++)  {
      gradRhoSpin(iSpin) = gradient (rhoXcR(iSpin));
   }
   SxMeshR rhoTotal;
   SxArray<SxMeshR> gradRhoTotal(3);
   if (nSpin == 1)  {
      rhoTotal     = rhoXcR(0);
      gradRhoTotal = gradRhoSpin(0);
   } else  {
      // sum spin up and spin down
      rhoTotal = rhoXcR(0) + rhoXcR(1);
      for (int idir = 0; idir < 3; ++idir)  {
         gradRhoTotal(idir) = gradRhoSpin(0)(idir) + gradRhoSpin(1)(idir);
      }
   }

   // --- End of variable preparation ------------------------------------------
   

   SxVector3<TPrecRhoR> gradRhoUp, gradRhoDown;
   SxVector3<TPrecRhoR> gradRhoTl;

   int xcFlags = SxXCFunctional::ComputeEandV;
   if (calcC) xcFlags |= SxXCFunctional::ComputeC;
   if (calcX) xcFlags |= SxXCFunctional::ComputeX;

   SxXCFunctional xc(nSpin, xcFlags);
   int up = 0, down = nSpin - 1; // nSpin = 1 => up=down
   for (int i=0; i < meshSize; i++)  {

      // --- collect gradients
      for (int iDir = 0; iDir < 3; ++iDir)  {
         gradRhoUp(iDir)   = gradRhoSpin(up)(iDir)(i);
         gradRhoDown(iDir) = gradRhoSpin(down)(iDir)(i);
         gradRhoTl(iDir)   = gradRhoTotal(iDir)(i);
      }

      double gradRhoNormUp   = gradRhoUp.norm (),
             gradRhoNormDown = gradRhoDown.norm (),
             gradRhoNormTl   = gradRhoTl.norm ();

      // --- compute functional
      if (xcFunctional == PBE)  {
         xc.computePBE (rhoXcR(up)(i),   gradRhoNormUp,
                        rhoXcR(down)(i), gradRhoNormDown,
                        gradRhoNormTl);
      } else if (xcFunctional == PBE0)  {
         xc.computePBEHybrid (rhoXcR(up)(i),   gradRhoNormUp,
                              rhoXcR(down)(i), gradRhoNormDown,
                              gradRhoNormTl, SxXCFunctional::alphaHybrid);
      } else if (xcFunctional == HSE06)  {
         xc.computeHSE (rhoXcR(up)(i),   gradRhoNormUp,
                        rhoXcR(down)(i), gradRhoNormDown,
                        gradRhoNormTl);
      } else {
         SX_EXIT;
      }

      eXc += xc.eXc;
      for (iSpin = 0; iSpin < nSpin; ++iSpin)  {
         // derivative with respect to rho
         vXc(iSpin)(i) = xc.vXc(iSpin);
         // derivative wrt grad rho(spin) (note: d|g|/d g_x = g_x / |g|)
         double g = (iSpin == 0) ? gradRhoNormUp : gradRhoNormDown;
         fx_gradRho(iSpin)(i) = (g > 0.) ? (xc.eXc_grad(iSpin) / g) : 0.; 
      }
      // derivative wrt grad rho (total)
      fc_gradRho(i) = (gradRhoNormTl > 0.) ? xc.eXc_gradRhoTl / gradRhoNormTl 
                                           : 0.;
   }
   eXc *= rBasisPtr->dOmega;

   const SxRBasis &R = *rBasisPtr;
   const SxGBasis &G = R.getGBasis ();
   // --- add White & Bird contribution of gradients
   for (int iDir = 0; iDir < 3; ++iDir)  {
      for (iSpin = 0; iSpin < nSpin; ++iSpin)  {
         SxMeshR fxc_gradR = gradRhoSpin(iSpin)(iDir) * fx_gradRho(iSpin)
                           + gradRhoTotal      (iDir) * fc_gradRho;
          //vXc(iSpin) -= (R | (I*G.gVec.colRef(iDir) * (G|fxc_gradR))).real();
            vXc(iSpin) += (R | (G.gVec.colRef(iDir) * (G | fxc_gradR))).imag ();
      }
   }
}


void SxXC::enableExchange ()
{
   calcX = true;
}

void SxXC::disableExchange ()
{
   calcX = false;
}

void SxXC::enableCorrelation ()
{
   calcC = true;
}

void SxXC::disableCorrelation ()
{
   calcC = false;
}

bool SxXC::getExchangeStatus () const
{
   return calcX;
}

bool SxXC::getCorrelationStatus () const
{
   return calcC;
}


/*
void SxXC::testPBE ()
{
   Real8 H, H00, H01, H02, H03, H10, H11, H12, H13; 
   Real8 t, phi, phi0, phi1, phi0_3, phi1_3;
   SxArray<Real8> vC (2), dVc(2);
   Real8 vX, s, fX, rs, zeta;
   Real8 fxFac = pow(2.,1./3.);
   int i;
   SxString fileId;

   // ---- plot fig. 3 of ref2 and fig.1 of 6
   try { SxFSAction::rm ("pbe.fx-0.dat"); } catch (const SxException&) { }
   try { SxFSAction::rm ("pbe.fx-2.dat"); } catch (const SxException&) { }
   try { SxFSAction::rm ("pbe.fx-10.dat"); } catch (const SxException&) { }
   try { SxFSAction::rm ("pbe.fx-inf.dat"); } catch (const SxException&) { }
   SxList<double> rsVals; rsVals << 1e-10 << 2. << 10.;
   for (s=0; s < 3.5; s += 0.1)  {
      for (i=0; i < rsVals.getSize(); ++i)  {
         rs = rsVals(i);
         for (zeta=0.; zeta <= 1.001; zeta += 1.)  {
            fileId = SxString() + rs + "-" + zeta + ".dat";
            if      (fabs(zeta-0.) < 1e-10)  { exchangePBE (0.,s,0.,0.,&vX, &fX); } 
            else if (fabs(zeta-1.) < 1e-10)  { exchangePBE (0.,s/fxFac,0.,0.,&vX, &fX); 
               fX *= pow(2.,1./3.); }  
            else  { SX_EXIT; }
            (SxString() + s + "\t" + fX + "\n").appendToFile ("pbe.fx-"+fileId);
            phi = pow(1.+zeta,2./3.) + pow(1.-zeta,2./3.);
            t = s * 1.2277 / phi / sqrt(rs);
            correlationPBE (rs, 0., 0.5*phi, t, 0., 0., 0., &vC, &dVc, &H);
            (SxString() + s + "\t" + H + "\n").appendToFile ("pbe.fc-"+fileId);
            (SxString() + s + "\t" + (fX+H)+ "\n").appendToFile ("pbe.fxc-"+fileId);
         }
      }
   }

   // --- plot fig. 7,8 of ref2, i.e. zeta=0,1
   try { SxFSAction::rm ("pbe.H0-0.0.dat"); } catch (const SxException&) { }
   try { SxFSAction::rm ("pbe.H0-0.5.dat"); } catch (const SxException&) { }
   try { SxFSAction::rm ("pbe.H0-2.0.dat"); } catch (const SxException&) { }
   try { SxFSAction::rm ("pbe.H0-6.0.dat"); } catch (const SxException&) { }
   try { SxFSAction::rm ("pbe.H1-0.0.dat"); } catch (const SxException&) { }
   try { SxFSAction::rm ("pbe.H1-0.5.dat"); } catch (const SxException&) { }
   try { SxFSAction::rm ("pbe.H1-2.0.dat"); } catch (const SxException&) { }
   try { SxFSAction::rm ("pbe.H1-6.0.dat"); } catch (const SxException&) { }
   phi0 = 0.5 * ( pow (1. + 0., 2./3.) + pow (1. - 0., 2./3.) ); // ref2,(20)
   phi1 = 0.5 * ( pow (1. + 1., 2./3.) + pow (1. - 1., 2./3.) ); // ref2,(20)
   phi0_3 = phi0*phi0*phi0; phi1_3 = phi1*phi1*phi1;
   for (t=0.; t < 5.; t += 0.1)  {
      nSpin = 1;
      correlationPBE (1e-10, 0., phi0, t, 0., 0., 0., &vC, &dVc, &H00);
      correlationPBE (0.5,   0., phi0, t, 0., 0., 0., &vC, &dVc, &H01);
      correlationPBE (2.0,   0., phi0, t, 0., 0., 0., &vC, &dVc, &H02);
      correlationPBE (6.0,   0., phi0, t, 0., 0., 0., &vC, &dVc, &H03);
      nSpin = 2;
      correlationPBE (1e-10, 1., phi1, t, 0., 0., 0., &vC, &dVc, &H10);
      correlationPBE (0.5,   1., phi1, t, 0., 0., 0., &vC, &dVc, &H11);
      correlationPBE (2.0,   1., phi1, t, 0., 0., 0., &vC, &dVc, &H12);
      correlationPBE (6.0,   1., phi1, t, 0., 0., 0., &vC, &dVc, &H13);
      H00 /= phi0_3; H01 /= phi0_3; H02 /= phi0_3; H03 /= phi0_3;
      H10 /= phi1_3; H11 /= phi1_3; H12 /= phi1_3; H13 /= phi1_3;

      // --- zeta=0
      (SxString() + t + "\t" +H00 + "\n").appendToFile ("pbe.H0-0.0-dat");
      (SxString() + t + "\t" +H01 + "\n").appendToFile ("pbe.H0-0.5-dat");
      (SxString() + t + "\t" +H02 + "\n").appendToFile ("pbe.H0-2.0-dat");
      (SxString() + t + "\t" +H03 + "\n").appendToFile ("pbe.H0-6.0-dat");
      // --- zeta=1
      (SxString() + t + "\t" +H10 + "\n").appendToFile ("pbe.H1-0.0-dat");
      (SxString() + t + "\t" +H11 + "\n").appendToFile ("pbe.H1-0.5-dat");
      (SxString() + t + "\t" +H12 + "\n").appendToFile ("pbe.H1-2.0-dat");
      (SxString() + t + "\t" +H13 + "\n").appendToFile ("pbe.H1-6.0-dat");
   }
}
*/

void SxXC::computeXPotEnergy (const RhoR &rhoR)
{
   PrecEnergy eVEx = 0.;

   for (int iSpin = 0; iSpin < nSpin; iSpin++)
      eVEx += tr(rhoR(iSpin) * vEx(iSpin));

   sxprintf ("eVEx = %.15f\n", eVEx);
}

void SxXC::computeXCPotEnergy (const RhoR &rhoR)
{
   PrecEnergy eVEx, eVCor;

   eVxc = eVEx = eVCor = 0.;

   for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
      eVxc  += tr (rhoR(iSpin) * vXc(iSpin));
      eVEx  += tr (rhoR(iSpin) * vEx(iSpin));
      eVCor += tr (rhoR(iSpin) * vCor(iSpin));
   }

   sxprintf ("eVxc  = %.15f\n", eVxc);
   sxprintf ("eVEx  = %.15f\n", eVEx);
   sxprintf ("eVCor = %.15f\n", eVCor);
}

void SxXC::shiftXC (bool xc, bool x, bool c)
{
   SxMeshG            vG;

   for (int iSpin = 0; iSpin < nSpin; iSpin++)  {

      // --- xc potential
      if (xc)  {
         SX_CHECK(dynamic_cast<const SxRBasis*> (vXc(iSpin).getBasisPtr ()));
         const SxRBasis &R = *dynamic_cast<const SxRBasis*> 
                              (vXc(iSpin).getBasisPtr ());
         const SxGBasis &G = R.getGBasis ();
         vG         = vXc(iSpin).to (G);
//sxprintf ("xc potential: (G=0)-component: (%g,%g) --> shifted to 0\n",
//          vG(0).re, vG(0).im);
         vG(0) = 0.;
         vXc(iSpin) = vG.to (R);
      }

      // --- x potential
      if (x)  {
         SX_CHECK(dynamic_cast<const SxRBasis*> (vEx(iSpin).getBasisPtr ()));
         const SxRBasis &R = *dynamic_cast<const SxRBasis*> 
                              (vEx(iSpin).getBasisPtr ());
         const SxGBasis &G = R.getGBasis ();
         vG = vEx(iSpin).to (G);
         vG(0) = 0.;
         vEx(iSpin) = vG.to (R);
      }

      // --- c potential
      if (c)  {
         SX_CHECK(dynamic_cast<const SxRBasis*> (vCor(iSpin).getBasisPtr ()));
         const SxRBasis &R = *dynamic_cast<const SxRBasis*> 
                              (vCor(iSpin).getBasisPtr ());
         const SxGBasis &G = R.getGBasis ();
         vG = vCor(iSpin).to (G);
         vG(0) = 0.;
         vCor(iSpin) = vG.to (R);
      }
   }
}

bool SxXC::rereadTable (const SxSymbolTable *table)
{
   if (table && table->containsGroup ("xcMesh"))  {
      try {
         const SxSymbolTable *meshGrp =  table->getGroup ("xcMesh");
         SxCell cell;
         if (vXc.getSize () > 0 && vXc(0).getSize () > 0)
            cell = vXc(0).getBasis<SxRBasis> ().cell;
         else if (xcBasisPtr)
            cell = xcBasisPtr->cell;
         else
            cell = SxAtomicStructure (table->topLevel ()).cell;

         setXcMesh (meshGrp, cell);
      } catch (const SxException &e)  {
         e.print ();
         SX_EXIT;
      }
      return true;
   }
   if (xcBasisPtr != origXcBasis)  {
      cout << "| xc will be computed on ";
      if (origXcBasis)  {
         cout << origXcBasis->fft3d.mesh;
      } else {
         cout << "default";
         if (vXc.getSize () > 0 && vXc(0).getSize () > 0)
            cout << " (" << vXc(0).getBasis<SxRBasis> ().fft3d.mesh << ')';
      }
      cout << " mesh." << endl;
      xcBasisPtr = origXcBasis;
      return true;
   }
   return false;
}
