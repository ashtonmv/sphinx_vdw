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

#include <SxPAWRho.h>
#include <SxProjector.h>

SxPAWRho::SxPAWRho (const SxConstPtr<SxPAWPot> &pawPot)
   : potPtr (pawPot)
{
   SX_CHECK (potPtr);
}

/// Assign from a density
void SxPAWRho::operator= (const SxDensity &in)
{
   SX_CHECK (in.checkType<SxPAWRho> ());
   const SxPAWRho &inRef = in.getRef<SxPAWRho> ();
   pwRho = inRef.pwRho;
   Dij = inRef.Dij;
   blockAO = inRef.blockAO;
   potPtr = inRef.potPtr;
   //psCore = inRef.psCore;
}

void SxPAWRho::operator= (const SxPAWRho &in)
{
   pwRho = in.pwRho;
   Dij = in.Dij;
   blockAO = in.blockAO;
   potPtr = in.potPtr;
   //psCore = in.psCore;
}

///\name Vector-like operations
///@{
/// Add a density
void SxPAWRho::operator+= (const SxDensity &x)
{
   if (x.checkType<SxRho> ())  {
      // add a filtered, pure pseudo-density 
      cout << "Adding pseudo-rho to PAW rho" << endl;
      pwRho += x.getRef<SxRho> ();
      return;
   }
   SX_CHECK (x.checkType<SxPAWRho> ());
   const SxPAWRho &xRef = x.getRef<SxPAWRho> ();
   SX_CHECK (getNSpin () == xRef.getNSpin (),
             getNSpin (), xRef.getNSpin ());
   SX_CHECK (blockAO.getSize () == xRef.blockAO.getSize ()
             || xRef.blockAO.getSize () == 0,
             blockAO.getSize (), xRef.blockAO.getSize ());
   pwRho += xRef.pwRho;
   Dij += xRef.Dij;
   for (int iAO = 0; iAO < xRef.blockAO.getSize (); iAO++)
      *blockAO(iAO) += *xRef.blockAO(iAO);
   /*
   if (psCore.getSize () == 0 && xRef.psCore.getSize () > 0)
      psCore = xRef.psCore.getCopy ();
   */
}

/// Subtract a density
void SxPAWRho::operator-= (const SxDensity &x)
{
   SX_CHECK (x.checkType<SxPAWRho> ());
   const SxPAWRho &xRef = x.getRef<SxPAWRho> ();
   SX_CHECK (getNSpin () == xRef.getNSpin (),
             getNSpin (), xRef.getNSpin ());
   SX_CHECK (blockAO.getSize () == xRef.blockAO.getSize ()
             || xRef.blockAO.getSize () == 0,
             blockAO.getSize (), xRef.blockAO.getSize ());
   pwRho -= xRef.pwRho;
   Dij -= xRef.Dij;
   for (int iAO = 0; iAO < xRef.blockAO.getSize (); iAO++)
      *blockAO(iAO) -= *xRef.blockAO(iAO);
   /*
   if (psCore.getSize () == 0 && xRef.psCore.getSize () > 0)
      psCore = xRef.psCore.getCopy ();
   */
}

/// axpy-like operation
void SxPAWRho::plus_assign_ax (double a, const SxDensity &x)
{
   SX_CHECK (x.checkType<SxPAWRho> ());
   const SxPAWRho &xRef = x.getRef<SxPAWRho> ();
   SX_CHECK (getNSpin () == xRef.getNSpin (),
             getNSpin (), xRef.getNSpin ());
   SX_CHECK (blockAO.getSize () == xRef.blockAO.getSize (),
             blockAO.getSize (), xRef.blockAO.getSize ());
   pwRho.plus_assign_ax (a, xRef.pwRho);
   Dij.plus_assign_ax (a, xRef.Dij);
   for (int iAO = 0; iAO < blockAO.getSize (); iAO++)
      blockAO(iAO)->plus_assign_ax (a, *xRef.blockAO(iAO));

   /*
   if (psCore.getSize () == 0 && xRef.psCore.getSize () > 0)
      psCore = xRef.psCore.getCopy ();
   */
}

/// axpy-like operation
void SxPAWRho::plus_assign_aspin (double a, const SxDensity &x)
{
   SX_CHECK (x.checkType<SxPAWRho> ());
   const SxPAWRho &xRef = x.getRef<SxPAWRho> ();
   SX_CHECK (getNSpin () == 2, getNSpin ());
   SX_CHECK (xRef.getNSpin () == 1, xRef.getNSpin ());
   SX_CHECK (blockAO.getSize () == xRef.blockAO.getSize (),
             blockAO.getSize (), xRef.blockAO.getSize ());

   pwRho.plus_assign_aspin (a, xRef.pwRho);
   Dij.plus_assign_aspin (a, xRef.Dij);
   for (int iAO = 0; iAO < blockAO.getSize (); iAO++)
      blockAO(iAO)->plus_assign_aspin (a, *xRef.blockAO(iAO));

   /*
   if (psCore.getSize () == 0 && xRef.psCore.getSize () > 0)
      psCore = xRef.psCore.getCopy ();
   */
}

/// Scalar product
double SxPAWRho::operator| (const SxDensity &x) const
{
   SX_CHECK (x.checkType<SxPAWRho> ());
   SX_CHECK (potPtr);
   const SxPAWRho &xRef = x.getRef<SxPAWRho> ();
   SX_CHECK (getNSpin () == xRef.getNSpin (),
             getNSpin (), xRef.getNSpin ());

   // plane-wave part
   double res = (pwRho | xRef.pwRho);

   // get volume for volume averaging (prefactor convention)
   SX_CHECK (pwRho.rBasisPtr);
   double vol = pwRho.rBasisPtr->cell.volume;

   // --- sum 1-center part over spin and atoms
   SxRadialMesh rho1, rho2; 
   SX_LOOP2(iSpin, is)  {
      int nr = (int)potPtr->phiAE(is).nRows ();
      SxDiracVec<Double> r3 = potPtr->rad(is).cub ();
      SX_LOOP(ia) {
         SxDiracVec<Double> rho12(nr);
         rho12.set (0.);
         rho1 = potPtr->computeRho (Dij(iSpin,is,ia), (int)is, potPtr->phiAE(is));
         rho2 = potPtr->computeRho (xRef.Dij(iSpin,is,ia), (int)is, potPtr->phiAE(is));
         for (int lm = 0; lm < rho1.meshData.nCols (); ++lm)
            rho12 +=   rho1.meshData.colRef(lm)
                     * rho2.meshData.colRef(lm);
         rho1 = potPtr->computeRho (Dij(iSpin,is,ia), (int)is, potPtr->phiPS(is));
         rho2 = potPtr->computeRho (xRef.Dij(iSpin,is,ia), (int)is, potPtr->phiPS(is));
         for (int lm = 0; lm < rho1.meshData.nCols (); ++lm)
            rho12 -=   rho1.meshData.colRef(lm)
                     * rho2.meshData.colRef(lm);
         res += (rho12 * r3).integrate (potPtr->logDr(is)) / vol;
      }
   }

   // blockAO do not contribute to the scalar product/norm

   return res;
}

double SxPAWRho::normSqr () const
{
   SX_CHECK (potPtr);
   
   // plane-wave part
   double res = pwRho.normSqr ();
   
   // get volume for volume averaging (prefactor convention)
   SX_CHECK (pwRho.rBasisPtr);
   double vol = pwRho.rBasisPtr->cell.volume;

   // --- sum 1-center part over spin and atoms
   SxRadialMesh rho; 
   SX_LOOP2(iSpin, is)  {
      int nr = (int)potPtr->phiAE(is).nRows ();
      SxDiracVec<Double> r3 = potPtr->rad(is).cub ();
      SxDiracVec<Double> rho2(nr);
      rho2.set (0.);
      SX_LOOP (ia)  {
         rho = potPtr->computeRho (Dij(iSpin,is,ia),(int)is,potPtr->phiAE(is));
         for (int lm = 0; lm < rho.meshData.nCols (); ++lm)
            rho2 += rho.meshData.colRef(lm).sqr ();
         rho = potPtr->computeRho (Dij(iSpin,is,ia),(int)is,potPtr->phiPS(is));
         for (int lm = 0; lm < rho.meshData.nCols (); ++lm)
            rho2 -= rho.meshData.colRef(lm).sqr ();
      }
      res += (rho2 * r3).integrate (potPtr->logDr(is)) / vol;
   }

   // blockAO do not contribute to the scalar product/norm

   return fabs(res); // physically res >=0. But numerically?
}

SxDensity SxPAWRho::operator- (const SxDensity &x) const
{
   SX_CHECK (x.checkType<SxPAWRho> ());
   SX_CHECK (potPtr);
   const SxPAWRho &xRef = x.getRef<SxPAWRho> ();
   SX_CHECK (getNSpin () == xRef.getNSpin (),
             getNSpin (), xRef.getNSpin ());
   SX_CHECK (blockAO.getSize () == xRef.blockAO.getSize (),
             blockAO.getSize (), xRef.blockAO.getSize ());

   int nSpin = getNSpin ();
   SxPtr<SxPAWRho> res = SxPtr<SxPAWRho>::create (potPtr);
   res->pwRho.rBasisPtr = pwRho.rBasisPtr;
   res->pwRho.rhoR.resize (nSpin);
   // plane-wave part
   for (int iSpin = 0; iSpin < nSpin; ++iSpin)
      res->pwRho(iSpin) = pwRho(iSpin) - xRef.pwRho(iSpin);
   // Dij
   res->Dij.resize (Dij.atomInfo, nSpin);
   res->Dij.computeDiff (Dij, xRef.Dij);
   res->blockAO.resize (blockAO.getSize ());
   for (int iAO = 0; iAO < blockAO.getSize (); iAO++)  {
      res->blockAO(iAO) = SxPtr<SxBlockDensityMatrix>::create ();
      res->blockAO(iAO)->resize (blockAO(iAO)->getNSpin (),
                                 blockAO(iAO)->getNSite ());
      res->blockAO(iAO)->computeDiff (*blockAO(iAO), *xRef.blockAO(iAO));
   }
   /*
   if (psCore.getSize () == pwRho(0).getSize ())
      res->psCore = psCore.getCopy ();
   else if (xRef.psCore.getSize () == pwRho(0).getSize ())
      res->psCore = xRef.psCore.getCopy ();
   */
   return SxDensity(res);
}

SxDensity SxPAWRho::getCopy () const
{
   int nSpin = getNSpin ();
   SxPtr<SxPAWRho> res = SxPtr<SxPAWRho>::create (potPtr);
   res->pwRho.rBasisPtr = pwRho.rBasisPtr;
   res->pwRho.rhoR.resize (nSpin);
   res->Dij.resize (Dij.atomInfo, nSpin);
   for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
      // plane-wave part
      res->pwRho(iSpin).copy (pwRho(iSpin));
      // --- 1-center part
      res->Dij.copyRho (Dij);
   }
   res->blockAO.resize (blockAO.getSize ());
   for (int iAO = 0; iAO < blockAO.getSize (); iAO++)  {
      res->blockAO(iAO) = SxPtr<SxBlockDensityMatrix>::create ();
      res->blockAO(iAO)->copyRho (*blockAO(iAO));
   }
   /*
   if (psCore.getSize () == pwRho(0).getSize ())
      res->psCore = psCore.getCopy ();
   */
   return SxDensity(res);
}

SxDensity SxPAWRho::spin () const
{
   SX_CHECK (getNSpin () == 2, getNSpin ());
   SxPtr<SxPAWRho> res = SxPtr<SxPAWRho>::create (potPtr);
   res->pwRho.rBasisPtr = pwRho.rBasisPtr;
   res->pwRho.rhoR.resize (1);
   res->pwRho(0) = (pwRho(0) - pwRho(1));
   // --- 1-center part
   res->Dij.resize (Dij.atomInfo, 1);
   SX_LOOP2(is,ia)
      res->Dij(0,is,ia) = Dij(0,is,ia) - Dij(1,is,ia);
   /*
   if (psCore.getSize () == pwRho(0).getSize ())
      res->psCore = psCore.getCopy ();
   */
   res->blockAO.resize (blockAO.getSize ());
   for (int iAO = 0; iAO < blockAO.getSize (); iAO++)  {
      res->blockAO(iAO) = SxPtr<SxBlockDensityMatrix>::create ();
      res->blockAO(iAO)->resize (1, blockAO(iAO)->getNSite ());
      SX_LOOP(iSite)
         res->blockAO(iAO)->rho(0,iSite) = (*blockAO(iAO))(0,iSite)
                                         - (*blockAO(iAO))(1,iSite);
   }
   return SxDensity(res);
}

/// Check if this is a spin-polarized density
bool SxPAWRho::hasSpin () const
{
   return getNSpin () == 2;
}


/// Renormalize
void SxPAWRho::renormalize ()
{
   cout << "TODO: PAW renormalize" << endl;
}

void SxPAWRho::readRho (const SxString &file)
{
   SxBinIO io;
   try {
      io.open (file, SxBinIO::BINARY_READ_ONLY);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   pwRho.readRho (io);
   Dij.readRho (io);
   if (io.contains("extraRho"))  {
      int nExtra = io.getDimension ("nExtraRho");
      int nSite = io.getDimension ("nSite");
      SxVector<Int> aoBlockSize(nSite), nExtraSite(nExtra);
      io.read ("nExtraSite", &nExtraSite, nExtra);
      io.read ("extraRhoDim", &aoBlockSize, nSite);
      blockAO.resize (nExtra);
      int iTlSite = 0, offset = 0;
      SX_LOOP(iAO)  {
         blockAO(iAO) = SxPtr<SxBlockDensityMatrix>::create ();
         blockAO(iAO)->resize (getNSpin (), nExtraSite(iAO));
         SX_LOOP(iSite) {
            int N = aoBlockSize(iTlSite++);
            SX_LOOP(iSpin)  {
               (*blockAO(iAO))(iSpin,iSite).reformat (N, N);
            }
         }
         blockAO(iAO)->readRho (io, "extraRho", &offset);
      }
   }

   if (potPtr)  {
      if (Dij.getNSpecies () != potPtr->getNSpecies ())  {
         cout << "Unexpected number of species in '" << file << "'." << endl;
         cout << "Found " << Dij.getNSpecies ()
              << ", but expected " << potPtr->getNSpecies () << endl;
         SX_QUIT;
      }
      for (int is = 0; is < Dij.getNSpecies (); ++is)
      if (Dij(0,is,0).nCols () != potPtr->getNProj (is))  {
         cout << "Unexpected number of partials in '" << file
              << "' for species " << (is+1) << endl;
         cout << "Found " << Dij(0,is,0).nCols ()
              << ", but expected " << potPtr->getNProj (is) << endl;
         SX_QUIT;
      }

   }
}

void SxPAWRho::writeRho (const SxString &file) const
{
   SX_CHECK (pwRho.rBasisPtr);
   const SxRBasis &R  = *pwRho.rBasisPtr;
   const SxAtomicStructure *strPtr = R.getGBasis ().structPtr;
   try  {
      SxBinIO io(file, SxBinIO::BINARY_WRITE_ONLY);
      const SxVector3<Int> mesh = R.getMesh ();
      io.writeMesh (pwRho.rhoR, R.cell, mesh);  // TODO: ugly!!!
      if (strPtr) {
         SX_CHECK (Dij.getNSpecies () == strPtr->getNSpecies (),
                   Dij.getNSpecies (), strPtr->getNSpecies ());
#ifndef NDEBUG
         for (int is = 0; is < Dij.getNSpecies (); ++is)
            SX_CHECK (Dij.getNAtoms(is) == strPtr->getNAtoms(is),
                      Dij.getNAtoms(is), strPtr->getNAtoms(is));
#endif
         strPtr->write (io);
      }
      Dij.writeRho (io);
      int nSite = 0;
      SxVector<Int> aoBlockSize, nExtraSite;
      if (blockAO.getSize () > 0)  {
         nExtraSite.resize (blockAO.getSize ());
         SX_LOOP(iAO) {
            nSite += blockAO(iAO)->getNSite ();
            nExtraSite(iAO) = blockAO(iAO)->getNSite ();
         }
         io.addDimension ("nSite", nSite);
         aoBlockSize.resize (nSite);
         io.addDimension ("nExtraRho", (int)blockAO.getSize ());
         io.write ("nExtraSite", nExtraSite, "nExtraRho");
         int nElem = 0;
         for (int iAO = 0, iTlSite = 0; iAO < blockAO.getSize (); iAO++)  {
            for (int iSite = 0; iSite < blockAO(iAO)->getNSite (); iSite++, iTlSite++)  {
               SX_CHECK (   (*blockAO(iAO))(0, iSite).nRows ()
                         == (*blockAO(iAO))(0, iSite).nCols (),
                         (*blockAO(iAO))(0, iSite).nRows (),
                         (*blockAO(iAO))(0, iSite).nCols ());
               aoBlockSize(iTlSite) = (int)(*blockAO(iAO))(0, iSite).nRows ();
               nElem += aoBlockSize(iTlSite) * aoBlockSize(iTlSite);
            }
         }
         io.addDimension ("nCoeffExtraRho", nElem);
         io.write ("extraRhoDim", aoBlockSize, "nSite");
         io.addDoubleVar ("extraRho", "nMeshes", "nCoeffExtraRho");
      }

      io.setMode (SxBinIO::WRITE_DATA);
      io.writeMesh (pwRho.rhoR, R.cell, mesh);  // TODO: ugly!!!
      if (strPtr) strPtr->write (io);
      Dij.writeRho (io);
      if (blockAO.getSize () > 0)  {
         io.write ("extraRhoDim", aoBlockSize, "nSite");
         io.write ("nExtraSite", nExtraSite, "nExtraRho");
         int offset = 0;
         SX_LOOP(iAO) blockAO(iAO)->writeRho (io, "extraRho", &offset);
      }
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}

double SxPAWRho::getSpin () const
{
   SX_CHECK (pwRho.rBasisPtr);
   SX_CHECK (potPtr);
   const SxPAWPot &pawpot = *potPtr;
   double rhoSpin = 0.;
   int nSpin = getNSpin ();
	if (nSpin == 1) return 0.;
   for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
   
		if (iSpin==0) {
			// --- pw part
			rhoSpin += pwRho(iSpin).sum () * pwRho.rBasisPtr->dOmega;
		} else {
			rhoSpin -= pwRho(iSpin).sum () * pwRho.rBasisPtr->dOmega;
		}
      for (int is = 0; is < Dij.getNSpecies (); ++is) {
			for (int ia = 0; ia < Dij.getNAtoms (is); ++ia)  {
				for (int ipt = 0; ipt < pawpot.getNProjType (is); ++ipt)  {
				 	int offsetI = pawpot.offset(is)(ipt);
					for (int jpt = 0; jpt < pawpot.getNProjType (is); ++jpt)  {
						if (pawpot.lPhi(is)(ipt) == pawpot.lPhi(is)(jpt))  {
							int offsetJ = pawpot.offset(is)(jpt);
							int nm = 2 * pawpot.lPhi(is)(ipt) + 1;
							for (int m = 0; m < nm; ++m)  {
								if (iSpin == 0) {
									rhoSpin += Dij(iSpin,is,ia)(offsetI + m, offsetJ + m)
												 * pawpot.deltaS(is)(ipt,jpt);
								} else {
                           rhoSpin -= Dij(iSpin,is,ia)(offsetI + m, offsetJ + m)
                                     * pawpot.deltaS(is)(ipt,jpt);
                        }
                     }
                  }
               }
            }
         }
      }
   }
   return rhoSpin;
}

SxArray<double> SxPAWRho::getSpinMom (const SxAtomicStructure &str) const
{
     SX_CHECK (potPtr);
     const SxPAWPot &pawpot = *potPtr;
     int nAtom = str.getNAtoms();
     SxArray<double> SpinMom(nAtom);
     for (int iAtom = 0; iAtom < nAtom; iAtom++) {
        int is = str.getISpecies(iAtom);
        int ia = iAtom - str.atomInfo->offset(is);
        int npt = pawpot.getNProjType (is);
        SpinMom (iAtom) = 0;
        for (int ipt = 0; ipt < npt; ipt++) {
           int li = pawpot.lPhi(is)(ipt);
           int offsetI = pawpot.offset(is)(ipt) + li;
           for (int jpt = 0; jpt < npt; jpt++) {
              int lj = pawpot.lPhi(is)(jpt);
              int offsetJ = pawpot.offset(is)(jpt) + lj;
              if ( li == lj) {
                 double logdr = pawpot.logDr (is);
                 SxDiracVec<Double> rad = pawpot.rad(is);
                 double rPAW = pawpot.rc (is) * sqrt(-log(1e-4)); // calculation of radius of the PAW sphere
                 SxDiracVec<Double> rW (rad.getSize ());
                 for (int ir = 0; ir < rad.getSize (); ++ir)  {
                    if (rad(ir) < rPAW) rW(ir) = rad(ir) * rad(ir) * rad(ir);
                    else rW(ir) = 0.;
                 }

                 double integral = (rW * pawpot.phiAE(is).colRef(jpt) * pawpot.phiAE(is).colRef(ipt)).integrate (logdr);
                 for (int m=-li; m<=li; m++) {
                    SpinMom (iAtom)  += integral * (Dij(0,is,ia)(offsetI+m,offsetJ+m) - Dij(1,is,ia)(offsetI+m,offsetJ+m));
                 }
              }
           }
        }
     }
     return SpinMom;
}

double SxPAWRho::getNorm () const
{
   SX_CHECK (pwRho.rBasisPtr);
   SX_CHECK (potPtr);
   const SxPAWPot &pawpot = *potPtr;
   double rhoNorm = 0.;
   int nSpin = getNSpin ();
   for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
      // --- pw part
      rhoNorm += pwRho(iSpin).sum () * pwRho.rBasisPtr->dOmega;

      for (int is = 0; is < Dij.getNSpecies (); ++is) {
         for (int ia = 0; ia < Dij.getNAtoms (is); ++ia)  {
            for (int ipt = 0; ipt < pawpot.getNProjType (is); ++ipt)  {
               int offsetI = pawpot.offset(is)(ipt);
               for (int jpt = 0; jpt < pawpot.getNProjType (is); ++jpt)  {
                  if (pawpot.lPhi(is)(ipt) == pawpot.lPhi(is)(jpt))  {
                     int offsetJ = pawpot.offset(is)(jpt);
                     int nm = 2 * pawpot.lPhi(is)(ipt) + 1;
                     for (int m = 0; m < nm; ++m)  {
                        rhoNorm += Dij(iSpin,is,ia)(offsetI + m, offsetJ + m)
                                 * pawpot.deltaS(is)(ipt,jpt);
                     }
                  }
               }
            }
         }
      }
   }
   return rhoNorm;
}
////// Meine Baustelle
//   Die Funktion fuer die atomicChargeDensity anpassen

void SxPAWRho::atomicChargeDensity (const SxAtomicStructure    &str,
                                    const SxConstPtr<SxPAWPot> &potPtrIn,
                                    const SxRBasis             &R,
                                    SxArray<SxArray<double> >  &atomSpin)
{
   SX_CHECK (potPtrIn);
   potPtr = potPtrIn;
   const SxPAWPot &pot = *potPtr;
   if (!pwRho.rBasisPtr) pwRho.rBasisPtr = &R;
   
   // --- basis sets
   const SxGBasis &G = R.getGBasis ();                 // |G>

   int nSpin = 2; // we set up a spin-polarized density...
   int nSpecies = str.getNSpecies ();
   int iSpecies, l;
   int iSpin;

   // resize pw part
   pwRho.rhoR.resize (nSpin);
   // --- 1-center part
   Dij.resize (str.atomInfo, pot, nSpin);
   Dij.set (0.);

   // --- print-out
   cout << SX_SEPARATOR;
   cout << "|  Initial occupations\n";
   cout << SX_SEPARATOR;
   for (iSpecies=0; iSpecies < nSpecies; iSpecies++)  {
      cout << "| " << pot.elementName(iSpecies) << ":  ";
      for (int ipt = 0; ipt < pot.lPhi(iSpecies).getSize (); ++ipt)  {
         l = pot.lPhi(iSpecies)(ipt);
         switch (l)  {
            case 0 : sxprintf ("  s="); break;
            case 1 : sxprintf ("  p="); break;
            case 2 : sxprintf ("  d="); break;
            case 3 : sxprintf ("  f="); break;           
            default: sxprintf ("  %d=", l);
         }
         sxprintf ("%2.2f", pot.foccInit(iSpecies)(ipt)); 
      }
      sxprintf ("\n");
   }	
   cout << SX_SEPARATOR;

   SxDiracVec<Double> r3, rhoRad, psRho;
   Real8 focc, Y00 = SQRT_1_4PI;
   Real8 logDr=0., norm, normCorePS = 0.;
   SxMeshG rhoUpG (G), rhoDownG(G), rhoAtomG;
   rhoUpG.set (0.);
   rhoDownG.set (0.);
   for (iSpecies=0; iSpecies < str.getNSpecies (); iSpecies++)  {
      logDr   = pot.logDr(iSpecies);
      r3 = pot.rad(iSpecies).cub ();
      rhoRad = SxDiracVec<Double> (pot.getRadBasis ()(iSpecies));
      // rhoRad <<= pot.rhoCorePS (iSpecies);
      rhoRad.set (0.);
      normCorePS += sqrt(FOUR_PI) * str.getNAtoms(iSpecies)
                  * (r3 * pot.rhoCorePS(iSpecies)).integrate (logDr);
      bool hasInitialRho = pot.rhoInit(iSpecies).getSize () > 0;
      double nElecPerAtom = 0.;
      for (int ipt = 0; ipt < pot.lPhi(iSpecies).getSize (); ++ipt)  {
         focc    = pot.foccInit(iSpecies)(ipt);

         if (!hasInitialRho)  {
            psRho   = pot.phiPS(iSpecies).colRef(ipt).sqr ();

            // --- enforce normalization via pseudo-rho
            norm    = (r3 * psRho).integrate (logDr)
               + pot.deltaS(iSpecies)(ipt,ipt);
            if  (fabs(norm - 1.) > 1e-6 && fabs(focc) > 1e-12)  {
               cout << "WARNING: " << pot.elementName(iSpecies)
                  << " channel " << (ipt+1) << "(l=" << pot.lPhi(iSpecies)(ipt)
                  << ") had to be renormalized (norm=" << norm << ")." << endl;
            }
            // always renormalize
            focc /= norm;

            // radial density
            rhoRad += (Y00 * focc) * psRho;
         }

         // --- 1-center density matrix
         int nm = 2 * pot.lPhi(iSpecies)(ipt) + 1;
         double fAvg = focc / (nm * nSpin);
         for (iSpin = 0; iSpin < nSpin; ++iSpin)  {
            for (int ia = 0; ia < str.getNAtoms(iSpecies); ++ia)  {
               for (int im = 0; im < nm; ++im)  {
                  int ipl = pot.offset(iSpecies)(ipt) + im;
                  Dij(iSpin,iSpecies, ia)(ipl, ipl) = fAvg;
               }
            }
         }
         nElecPerAtom += focc;
      }

      PsiG rhoCoreG = (G | pot.rhoCorePS (iSpecies));

      // --- Fourier transform to G
      if (hasInitialRho)  {
         cout << pot.prettyName(iSpecies);
         cout << ": Taking initial atomic density from potential file" << endl;
         rhoRad = (G | pot.rhoInit(iSpecies));
         rhoRad -= rhoCoreG;      // original line   rhoRad -= rhoCoreG;
      } else {
         rhoRad.handle->auxData.is = iSpecies;
         rhoRad.handle->auxData.l  = 0;
         rhoRad.handle->auxData.m  = 0;
      }
      rhoAtomG = ( G | rhoRad); // total density
      for (int ia = 0; ia < str.getNAtoms(iSpecies); ++ia)  {
         double s = atomSpin(iSpecies)(ia) / nElecPerAtom;
         Dij(0,iSpecies,ia) *= 1. + s;
         Dij(1,iSpecies,ia) *= 1. - s;
         PsiG T = G.getPhaseFactors(iSpecies,ia); // exp(-iGtau)
         rhoUpG   += 0.5 * (1. + s) * T * rhoAtomG + 0.5 * T * rhoCoreG;
         rhoDownG += 0.5 * (1. - s) * T * rhoAtomG + 0.5 * T * rhoCoreG;
      }
   }
   
   pwRho(0) = R.symmetrize ( R | rhoUpG );
   pwRho(1) = R.symmetrize ( R | rhoDownG );

   double nTot = getNorm ();
   cout << "Atomic charge density:" << nTot << endl;
   cout << "Core charge density:  " << normCorePS << endl;
   cout << "Valence density:      " << (nTot - normCorePS) << endl;

   double spinTotal = 0.;
   SX_LOOP2(is,ia) spinTotal += atomSpin(is)(ia);
   cout << "Total spin from labels:  " << spinTotal << endl;
   cout << "Total spin from PAW rho: " << getSpin () << endl; 


}


////// Hier endet meine Baustelle








void SxPAWRho::atomicChargeDensity (const SxAtomicStructure    &str,
                                    const SxConstPtr<SxPAWPot> &potPtrIn,
                                    const SxRBasis             &R,
                                          int                  nSpin)
{
   SX_CHECK (potPtrIn);
   potPtr = potPtrIn;
   const SxPAWPot &pot = *potPtr;
   if (!pwRho.rBasisPtr) pwRho.rBasisPtr = &R;
   
   // --- basis sets
   const SxGBasis &G = R.getGBasis ();                 // |G>

   int nSpecies = str.getNSpecies ();

   // resize pw part
   pwRho.rhoR.resize (nSpin);
   // --- 1-center part
   Dij.resize (str.atomInfo, pot, nSpin);
   Dij.set (0.);

   // --- print-out
   cout << SX_SEPARATOR;
   cout << "|  Initial occupations\n";
   cout << SX_SEPARATOR;
   for (int iSpecies=0; iSpecies < nSpecies; iSpecies++)  {
      cout << "| " << pot.elementName(iSpecies) << ":  ";
      for (int ipt = 0; ipt < pot.lPhi(iSpecies).getSize (); ++ipt)  {
         int l = pot.lPhi(iSpecies)(ipt);
         switch (l)  {
            case 0 : sxprintf ("  s="); break;
            case 1 : sxprintf ("  p="); break;
            case 2 : sxprintf ("  d="); break;
            case 3 : sxprintf ("  f="); break;           
            default: sxprintf ("  %d=", l);
         }
         sxprintf ("%2.2f", pot.foccInit(iSpecies)(ipt)); 
      }
      sxprintf ("\n");
   }	
   cout << SX_SEPARATOR;

   SxDiracVec<Double> r3, rhoRad, psRho;
   Real8 focc, Y00 = SQRT_1_4PI;
   Real8 logDr=0., norm, normCorePS = 0.;
   SxMeshG rhoG (G), rhoAtomG;
   rhoG.set (0.);
   //for (iSpecies=0; iSpecies < str.getNSpecies (); iSpecies++)  {
   SX_LOOP(iSpecies)  {
      logDr   = pot.logDr(iSpecies);
      r3 = pot.rad(iSpecies).cub ();
      rhoRad = SxDiracVec<Double> (pot.getRadBasis ()((int)iSpecies));
      rhoRad <<= pot.rhoCorePS (iSpecies);
      normCorePS += sqrt(FOUR_PI) * str.getNAtoms((int)iSpecies)
                  * (r3 * pot.rhoCorePS(iSpecies)).integrate (logDr);
      bool hasInitialRho = pot.rhoInit(iSpecies).getSize () > 0;
      //for (int ipt = 0; ipt < pot.lPhi(iSpecies).getSize (); ++ipt)  {
      SX_LOOP(ipt) {
         focc    = pot.foccInit(iSpecies)(ipt);

         if (!hasInitialRho)  {
            psRho   = pot.phiPS(iSpecies).colRef(ipt).sqr ();

            // --- enforce normalization via pseudo-rho
            norm    = (r3 * psRho).integrate (logDr)
               + pot.deltaS(iSpecies)(ipt,ipt);
            if  (fabs(norm - 1.) > 1e-6 && fabs(focc) > 1e-12)  {
               cout << "WARNING: " << pot.elementName(iSpecies)
                  << " channel " << (ipt+1) << "(l=" << pot.lPhi(iSpecies)(ipt)
                  << ") had to be renormalized (norm=" << norm << ")." << endl;
            }
            // always renormalize
            focc /= norm;

            // radial density
            rhoRad += (Y00 * focc) * psRho;
         }

         // --- 1-center density matrix
         int nm = 2 * pot.lPhi(iSpecies)(ipt) + 1;
         if (pot.DijInit(iSpecies).getSize () > 0 && hasInitialRho)  {
         //if (false) {
            SX_LOOP3(iSpin,ia,jpt)  {
               if (pot.lPhi(iSpecies)(ipt) != pot.lPhi(iSpecies)(jpt)) continue;
               for (int im = 0; im < nm; ++im)  {
                  int ipl = pot.offset(iSpecies)(ipt) + im;
                  int jpl = pot.offset(iSpecies)(jpt) + im;
                  Dij(iSpin,iSpecies,ia)(ipl,jpl)
                     = pot.DijInit(iSpecies)(ipt,jpt) / nSpin;
               }
            }
         }  else {
            double fAvg = focc / (nm * nSpin);
            SX_LOOP3(iSpin, ia, im(nm))  {
               ssize_t ipl = pot.offset(iSpecies)(ipt) + im;
               Dij(iSpin,iSpecies, ia)(ipl, ipl) = fAvg;
            }
         }
      }
      // --- Fourier transform to G
      if (hasInitialRho)  {
         cout << pot.prettyName((int)iSpecies);
         cout << ": Taking initial atomic density from potential file" << endl;
         rhoRad = pot.rhoInit(iSpecies);
      } else {
         rhoRad.handle->auxData.is = (int)iSpecies;
         rhoRad.handle->auxData.l  = 0;
         rhoRad.handle->auxData.m  = 0;
      }
      rhoAtomG = ( G | rhoRad);

      rhoG += G.structureFactors(iSpecies) * rhoAtomG;
   }
   
   pwRho(0) = R.symmetrize ( R | rhoG );
   if (nSpin == 2)  {
      // duplicate spin channels
      pwRho(0) *= 0.5;
      pwRho(1).copy (pwRho(0));
   }

   double nTot = getNorm ();
   cout << "Atomic charge density:" << nTot << endl;
   cout << "Core charge density:  " << normCorePS << endl;
   cout << "Valence density:      " << (nTot - normCorePS) << endl;

}

void SxPAWRho::atomicChargeDensity (const SxAtomicStructure    &str,
                                    const SxConstPtr<SxPAWPot> &potPtrIn,
                                    const SxRBasis             &R,
                                    const SxArray2<SxVector<Double> > &focc)
{
   int nSpin = int(focc.getDim (0));
   SX_CHECK (focc.getDim (1) == str.getNAtoms (),
             focc.getDim (1), str.getNAtoms ());

   SX_CHECK (potPtrIn);
   potPtr = potPtrIn;
   const SxPAWPot &pot = *potPtr;
   if (!pwRho.rBasisPtr) pwRho.rBasisPtr = &R;
   
   // --- basis sets
   const SxGBasis &G = R.getGBasis ();                 // |G>

   int nSpecies = str.getNSpecies ();

   // resize pw part
   pwRho.rhoR.resize (nSpin);
   // --- 1-center part
   Dij.resize (str.atomInfo, pot, nSpin);
   Dij.set (0.);
   SxArray<SxArray<bool> > channelUsed (nSpecies);

   // --- print-out
   cout << SX_SEPARATOR;
   cout << "|  Initial occupations\n";
   cout << SX_SEPARATOR;
   for (int iSpecies=0, iTlAtom = 0; iSpecies < nSpecies; iSpecies++)  {
      cout << "| " << pot.elementName(iSpecies) << ":" << endl;
      channelUsed(iSpecies).resize (pot.lPhi(iSpecies).getSize ());
      channelUsed(iSpecies).set (false);
      for (int ia = 0; ia < str.getNAtoms (iSpecies); ++ia, ++iTlAtom)  {
         cout << "| atom " << (ia+1);

         int ja;
         for (ja = ia + 1; ja < str.getNAtoms (iSpecies); ++ja) {
            int jTlAtom = iTlAtom + ja - ia;
            if ( (focc(0, jTlAtom) - focc(0, iTlAtom)).normSqr () > 1e-18)
               break;
            if (nSpin == 1) continue;
            if ( (focc(1, jTlAtom) - focc(1, iTlAtom)).normSqr () > 1e-18)
               break;
         }
         ja--;

         if (ja > ia) cout << "..." << (ja + 1);
         cout << ": ";
         // --- show focc
         for (int ipt = 0; ipt < pot.lPhi(iSpecies).getSize (); ++ipt)  {
            int l = pot.lPhi(iSpecies)(ipt);
            switch (l)  {
               case 0 : sxprintf ("  s="); break;
               case 1 : sxprintf ("  p="); break;
               case 2 : sxprintf ("  d="); break;
               case 3 : sxprintf ("  f="); break;           
               default: sxprintf ("  (l=%d)=", l);
            }
            sxprintf ("%2.2f", focc(0, iTlAtom)(ipt));
            if (nSpin == 2)
               sxprintf ("/%2.2f", focc(1, iTlAtom)(ipt));

            for (int iSpin = 0; iSpin < nSpin; ++iSpin)
               if (fabs(focc(iSpin, iTlAtom)(ipt)) > 1e-9)
                  channelUsed(iSpecies)(ipt) = true;
         }
         sxprintf ("\n");

         iTlAtom += ja - ia;
         ia = ja;
      }
   }	
   cout << SX_SEPARATOR;

   Real8 Y00 = SQRT_1_4PI;
   Real8 normCorePS = 0.;
   SxMeshG rhoG (G), rhoAtomG;

   // --- pseudo core
   rhoG.set (0.);
   for (int iSpecies=0; iSpecies < str.getNSpecies (); iSpecies++)  {
      double logDr   = pot.logDr(iSpecies);
      SxDiracVec<Double> r3 = pot.rad(iSpecies).cub ();
      normCorePS += sqrt(FOUR_PI) * str.getNAtoms(iSpecies)
                  * (r3 * pot.rhoCorePS(iSpecies)).integrate (logDr);
      // --- Fourier transform to G
      rhoAtomG = ( G | pot.rhoCorePS (iSpecies));
      rhoG += G.structureFactors(iSpecies) * rhoAtomG;
   }
   pwRho(0) = R.symmetrize ( R | rhoG );
   if (nSpin == 2)  {
      // duplicate spin channels
      pwRho(0) *= 0.5;
      pwRho(1).copy (pwRho(0));
   }

   // --- valence occupation
   for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
      rhoG.set (0.);
      for (int iSpecies=0; iSpecies < str.getNSpecies (); iSpecies++)  {
         double logDr   = pot.logDr(iSpecies);
         SxDiracVec<Double> r3 = pot.rad(iSpecies).cub ();
         for (int ipt = 0; ipt < pot.lPhi(iSpecies).getSize (); ++ipt)  {
            
            // shortcut for unused channels
            if (!channelUsed(iSpecies)(ipt)) continue;

            int nm = 2 * pot.lPhi(iSpecies)(ipt) + 1;

            SxDiracVec<Double> rhoRad = pot.phiPS(iSpecies).colRef(ipt).sqr ();

            // --- enforce normalization via pseudo-rho
            double norm = (r3 * rhoRad).integrate (logDr)
                        + pot.deltaS(iSpecies)(ipt,ipt);
            if  (fabs(norm - 1.) > 1e-6)  {
               cout << "WARNING: " << pot.elementName(iSpecies)
                    << " channel " << (ipt+1) << "(l="
                    << pot.lPhi(iSpecies)(ipt) 
                    << ") had to be renormalized (norm=" << norm << ")." 
                    << endl;
            }
            // always renormalize
            rhoRad *= Y00 / norm;

            // --- Fourier transform to G
            rhoRad.handle->auxData.is = iSpecies;
            rhoRad.handle->auxData.ia = -1;
            rhoRad.handle->auxData.l  = 0;
            rhoRad.handle->auxData.m  = 0;
            rhoAtomG = ( G | rhoRad);

            // --- set up density with atom-specific focc
            for (int ia = 0; ia < str.getNAtoms(iSpecies); ++ia)  {
               double f = focc(iSpin, str.getIAtom(iSpecies,ia))(ipt);
               if (fabs(f) < 1e-12) continue;

               // --- pw density
               rhoG.plus_assign_ax (f, rhoAtomG 
                                       * G.getPhaseFactors (iSpecies, ia));

               // --- 1-center density matrix
               for (int im = 0; im < nm; ++im)  {
                  int ipl = pot.offset(iSpecies)(ipt) + im;
                  Dij(iSpin, iSpecies, ia)(ipl, ipl) = f / (norm * nm);
               }
            }
         }
      }
      pwRho(iSpin) += R.symmetrize ( R | rhoG );
   }

   double nTot = getNorm ();
   cout << "Atomic charge density:" << nTot << endl;
   cout << "Core charge density:  " << normCorePS << endl;
   cout << "Valence density:      " << (nTot - normCorePS) << endl;

}

void SxPAWRho::syncMPI ()
{
#ifdef USE_LOOPMPI
   pwRho.syncMPI ();
   Dij.syncMPI ();
#endif
}

