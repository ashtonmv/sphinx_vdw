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

#include <SxWannier.h>
#include <SxBasis.h>
#include <SxProjector.h>


#ifndef SX_STANDALONE


SxWannier::SxWannier ()
{
   meshDim = 0;
   dSpread = alpha = eCut = -1.;
   nk = nb = nSpin = iBottom = iTop = nBands = -1;
   repetition = 0;
   translation = 0;
   useGuess = false;
}

SxWannier::SxWannier (const SxSymbolTable *table)
{
   SxSymbolTable  *wavesGroup = NULL, *spreadGroup = NULL, *bandsGroup = NULL;
   SxSymbolTable  *graphics = NULL;

   try  {
      wavesGroup  = table->getGroup("waves");
      wavesfile   = wavesGroup->get("file")->toString();
      spreadGroup = table->getGroup("spread");
      maxSteps    = spreadGroup->get("maxSteps")->toInt();
      dSpread     = spreadGroup->get("dSpread")->toReal();
      alpha       = spreadGroup->get("alpha")->toReal();
      bandsGroup  = table->getGroup("bands");
      iBottom     = -1 + bandsGroup->get("nBottom")->toInt();  // starts at 0
      iTop        = -1 + bandsGroup->get("nTop")->toInt();     // starts at 0
      if (table->containsGroup("graphics"))  {
         graphics    = table->getGroup("graphics");
         repetition  = SxVector3<Int> (graphics->get("repetition")
                                               ->toIntList());
         translation = SxVector3<Int> (graphics->get("translation")
                                               ->toIntList());
      }

      // --- read waves file
      SxBinIO io (wavesfile, SxBinIO::BINARY_READ_ONLY);
      SxPtr<SxGkBasis> gkSetPtr = SxPtr<SxGkBasis>::create (io, true);
      gkSet = *gkSetPtr;
      gkSet.setupN123inv ();
      waves.read (io);
      waves.setGkBasisPtr (gkSetPtr);
      cell.read (io);
      io.read ("eCut", &eCut);
      io.read ("meshDim", &meshDim);
      io.close ();

      // --- get number of spin-channels, k-points, and bands
      nSpin = waves.getNSpin ();
      nk = waves.getNk ();
      nBands = iTop - (iBottom - 1);

      // --- perform initialisation, if required
      if (table->containsGroup("initialization"))  {
         initialGuess = SxGauss (waves.getGkBasisPtr(), cell, &*table);
         initOrbitals ();
         uInitPtr = &(initialGuess.projections);
         useGuess = true;
      }  else  {
         uInitPtr = &waves;
         useGuess = false;
      }
   }  catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }

   // --- initialise finite differences in k space
   initKDiffs ();
   nb = kDiffs.getNb ();  // number of k-points in the finite diff. formula

   // --- check input parameters iBottom, iTop
   if (iBottom < 0)  {
      sxprintf ("\n"
                "Error:  iBottom chosen too small. "
                "Must be 1 or larger. Sorry.\n");
      SX_QUIT;
   }

   if (iTop >= waves.getNStates())  {
      sxprintf ("\n"
                "Error:  iTop chosen too large. Must not be larger than number "
                "of waves in\n"
                "        wavesfile (%d). Sorry.\n",
                waves.getNStates());
      SX_QUIT;
   }

   // --- resize initial overlap tensor M
   int iSpin, ik, ib;
   MkbInit.resize (nk);

   for (ik = 0; ik < nk; ik++)  {
      MkbInit(ik).resize (nb);

      for (ib = 0; ib < nb; ib++)  {
         MkbInit(ik)(ib).resize (nSpin);

         for (iSpin = 0; iSpin < nSpin; iSpin++)  {
            MkbInit(ik)(ib)(iSpin).reformat (nBands,nBands);
         }
      }
   }

   // --- resize transformation tensors U and D
   Uk.resize (nSpin);
   D.resize (nSpin);

   for (iSpin = 0; iSpin < nSpin; iSpin++)  {
      Uk(iSpin).reformat (nk, nBands, nBands);
      D(iSpin).reformat (nk, nBands, nBands);
   }

   // --- initialise |G> basis - TODO: probably not necessary
   G.set (meshDim, cell, 4. * eCut);
   initialGuess.setGBasis (G);
}

void SxWannier::print () const
{
   // --- dump parameters from input file
   cout << SX_SEPARATOR;
   sxprintf ("| Input parameters\n");
   cout << SX_SEPARATOR;
   sxprintf ("| Bloch wavefunctions:\n");
   sxprintf ("|\n");
   sxprintf ("|   %s\n", wavesfile.ascii());
   sxprintf ("|\n");
   sxprintf ("| Physics:\n");
   sxprintf ("|\n");
   sxprintf ("|   nBottom:      %d (lowest band index considered, "
             "minimum = 1)\n", iBottom+1);
   sxprintf ("|   nTop:         %d (highest band index considered, "
             "minimum = 1)\n", iTop+1);
   sxprintf ("|   dSpread:      %g\n", dSpread);
   sxprintf ("|   alpha:        %g\n", alpha);
   sxprintf ("|   meshDim:      (%d,%d,%d)",
             meshDim(0), meshDim(1), meshDim(2));
   sxprintf ("  -->  size of real space mesh\n");
   sxprintf ("|   eCut:         %g\n", eCut);
   sxprintf ("|\n");
   sxprintf ("| Graphics:\n");
   sxprintf ("|\n");
   sxprintf ("|   repetition:   (%d,%d,%d)",
             repetition(0), repetition(1), repetition(2));
   sxprintf ("     -->  repeating factors for real space mesh\n");
   sxprintf ("|   translation:  (%d,%d,%d)",
             translation(0), translation(1), translation(2));
   sxprintf ("     -->  shifts Wannier functions from origin\n");

   sxprintf ("|\n");
   sxprintf ("| Initial guess:\n");
   sxprintf ("|\n");
   if (useGuess)  {
      initialGuess.print ();
   }  else  {
      sxprintf ("|   Bloch waves\n");
   }
}

//-----------------------------------------------------------------------------
//   initialization of orbitals using Gaussian bells
//-----------------------------------------------------------------------------
void SxWannier::initOrbitals ()
{
   cout << SX_SEPARATOR;
   cout << "| Initialise orbitals" << endl;

   SX_CHECK (iBottom >= 0 && iTop >= iBottom, iBottom, iTop);

   initialGuess.setRealSpaceMesh (SxGauss::FFT);
   initialGuess.computeGaussiansFFTMesh ();
   initialGuess.write ();
   initialGuess.projectOld (waves, iBottom, iTop);
   initialGuess.orthonormalizeProjections (Loewdin);
   
   // --- get initial unitary transformation
   if (nSpin > 1)  SX_QUIT;  // spin polarisation not yet implemented
   P.reformat (nk, nBands, nBands);
   SxArray<SxArray<SxDiracMat<TPrecCoeffG> > > p1, p2;
   p1 = initialGuess.getProjectionTrafo();
   p2 = initialGuess.getLoewdinTrafo();

   for (int ik = 0; ik < nk; ik++)  P.set (p1(ik)(0) ^ p2(ik)(0), ik);

   P.writeDirTensor3 ("P.sxb");
}

void SxWannier::initKDiffs ()
{
   kDiffs = SxKFinDiffs (cell, gkSet);
   kDiffs.checkKSur ();
}

void SxWannier::testScalarProduct ()
{
   PsiG       a, b;
#ifndef NDEBUG
   SxGBasis  &gaBas = gkSet(0);
#endif
   SxGBasis  &gbBas = gkSet(7);

   //a = PsiG (initialGuess.projections(0,0,0));
   //b = PsiG (initialGuess.projections(1,0,1));
   a = PsiG (waves(3,0,0));
   b = PsiG (waves(0,0,7));

   SX_CHECK (a.getBasisPtr() == &gaBas);
   SX_CHECK (b.getBasisPtr() == &gbBas);

   const SxGBasis *aBasPtr = dynamic_cast<const SxGBasis *> (a.getBasisPtr());
   SX_CHECK (aBasPtr);

   aBasPtr->n123invSetup ();

   SxComplex16 spr = gbBas.scalarProduct (*aBasPtr, a, b);
   cout << "spr = " << spr << endl;
}

void SxWannier::testKShift1 ()
{
   int                 ik = 2;
   SxGBasis           &gBasis = gkSet(ik);
   PsiG                f = PsiG (waves(0,0,ik));
   SxVector3<Double>   shiftVec;
   SxArray<Coord>      rMesh;
   SxRBasis            R;
   R.set (meshDim, cell);

   gBasis.registerBasis (R);

   shiftVec = cell.getReciprocalCell().relToCar (RelVec (0,-1,0));
   rMesh = initialGuess.computeMeshPoints (meshDim, RelVec (1,1,1), cell);

   cout << "shiftVec = " << shiftVec << endl;

   //PsiG fShifted = gBasis.shiftK (f, shiftVec, rMesh, &R);
   PsiG fShifted = gBasis.shiftK (f, shiftVec, rMesh);
   cout << "test 1:  fShifted = " << fShifted << endl;
   cout << endl;
   //cout << "f = " << f << endl;
}

void SxWannier::testKShift2 ()
{
   int              ik = 2;
   SxGBasis        &gBasis = gkSet(ik);
   PsiG             f = PsiG (waves(0,0,ik));
   SxVector3<Int>   deltaG (0,1,0);

   gBasis.n123invSetup ();

   PsiG fShifted = gBasis.shiftK (f, deltaG);
   cout << "test 2:  fShifted = " << fShifted << endl;

   // scalar test
   PsiG a;
   SxGBasis &gaBas = gkSet(5);

   a = PsiG (waves(1,0,5));
   SX_CHECK (a.getBasisPtr() == &gaBas);

   gaBas.n123invSetup ();

   SxComplex16 spr = gBasis.scalarProduct (gaBas, a, fShifted);
   cout << "test 2: spr = " << spr << endl;
}

// ----------------------------------------------------------------------------
//    calculate matrix elements < u_n,k | u_m,k+b >
// ----------------------------------------------------------------------------
void SxWannier::initOverlapMatrices ()
{
   cout << SX_SEPARATOR;
   cout << "| Calculating initial overlap matrices ... ";  cout.flush ();

   int           iSpin, ik, iq, ib, i, j;
   SxGBasis     *gkPtr, *gqPtr;
   RelVec        dGRel;  // TODO: Nachnormieren im Falle von rel. coord.
   PsiG          u1, u2;
   SxComplex16   spr;
   const SxPW   &uInit = *uInitPtr;

   for (ik = 0; ik < nk; ik++)  {
      kDiffs.computeSurroundingKPoints (ik);
      gkPtr = &( gkSet(ik) );

      for (ib = 0; ib < nb; ib++)  {
         iq    = kDiffs.kSurIdx(ib);
         dGRel = kDiffs.deltaGRel(ib);
         gqPtr     = &( gkSet(iq) );

         for (iSpin = 0; iSpin < nSpin; iSpin++)  {
            for (i = iBottom; i <= iTop; i++)  {
               u1 = uInit(i,iSpin,ik);

               for (j = iBottom; j <= iTop; j++)  {
                  u2  = gqPtr->shiftK (uInit(j,iSpin,iq), dGRel);
                  spr = gqPtr->scalarProduct (*gkPtr, u1, u2);
                  SX_CHECK (spr.abs() < 1.+1.e-4, spr.abs());
                  MkbInit(ik)(ib)(iSpin)(i,j) = spr;
               }  // :j
            }  // :i
         }  // :iSpin
      }  // :ib
   }  // :ik

   cout << "done" << endl;
   writeMkb (MkbInit, "Mkb-init.dat");
}

//-----------------------------------------------------------------------------
//   Wannier centers:  r(bar)_n := <0n|r|0n>
//-----------------------------------------------------------------------------
SxWannier::WFcenters
SxWannier::computeCenters (const SxWannier::WFTensor5 &Mkb) const
{
   // works out Eq. (31) in PRB B, 56, 12847 (1997)
   SX_CHECK (nBands == iTop - (iBottom - 1), nBands, iTop, iBottom);
   SX_CHECK   (nk > 0, nk);
   SX_CHECK  (nb == kDiffs.getNb(), nb, kDiffs.getNb());

   int                    ik, ib, iSpin, i;
   double                 wb;
   SxVector3<Double>      bVec;
   WFcenters              wrn(nSpin);
   SxDiracMat<Complex16>  M;

   for (iSpin = 0; iSpin < nSpin; iSpin++)  {
      wrn(iSpin).resize (nBands);

      for (i = 0; i < nBands; i++)  {
         wrn(iSpin)(i).set (0.);
      }

      for (ib = 0; ib < nb; ib++)  {
         bVec = kDiffs.bVecs(ib);
         wb   = kDiffs.bWeights(ib);

         for (ik = 0; ik < nk; ik++)  {
            M = Mkb(ik)(ib)(iSpin);  // reference

            for (i = 0; i < nBands; i++)  {
               wrn(iSpin)(i) -= wb * bVec
                                * (lnC ( M.row(i)(i) )).im
                                / (double)nk;
            }  // :i
         }  // :ik
      }  // :ib
   }  // :iSpin

   return wrn;
}

//----------------------------------------------------------------------------
//  spread functional
//----------------------------------------------------------------------------
SxWannier::WFspreads
SxWannier::computeSpreads (const SxWannier::WFTensor5 &Mkb,
                           const SxWannier::WFcenters &wrn) const
{
   // works out Eq. (7) in PRB B, 65, 035109 (2001)
   SX_CHECK (nBands == iTop - (iBottom - 1), nBands, iTop, iBottom);
   SX_CHECK   (nk > 0, nk);
   SX_CHECK  (nb == kDiffs.getNb(), nb, kDiffs.getNb());

   cout << "| Spread functional:\n|" << endl;

   int                    ik, ib, iSpin, i;
   WFspreads              omega(2);
   double                 omegaI, omegaD, omegaOD;
   double                 dSum, nSum, x, wb;
   double                 toAngstSquare = 1. / (1.8897261 * 1.8897261);
   SxComplex16            y;
   SxDiracMat<Complex16>  M;
   SxDiracVec<Complex16>  diagM;
   Coord                  bVec, rnVec;

   for (iSpin = 0; iSpin < nSpin; iSpin++)  {
      if (nSpin > 1)  {
         if (iSpin == SPIN_UP)  cout << "|\n|   Spin up:" << endl;
         else                   cout << "|\n|   Spin down:" << endl;
      }
      omega(iSpin) = 0.;
      omegaI = omegaD = omegaOD = 0.;

      for (ib = 0; ib < nb; ib++)  {
         bVec = kDiffs.bVecs(ib);
         wb   = kDiffs.bWeights(ib);

         for (ik = 0; ik < nk; ik++)  {
            M     = Mkb(ik)(ib)(iSpin);  // reference
            diagM = M.diag();
            dSum  = (diagM.conj() * diagM).sum();
            y     = (M.conj() * M).sum();

            omegaI  += wb * ((double)nBands - y.re);
            omegaOD += wb * (double)(y - dSum);

            nSum = 0.;
            for (i = 0; i < nBands; i++)  {
               rnVec  = wrn(iSpin)(i);
               x      = (lnC (diagM(i))).im + (bVec * rnVec).sum();
               nSum  += x * x;
            }  // :i

            omegaD += wb * nSum;
         }  // :ik
      }  // :ib

      omegaI       *= toAngstSquare / ((double)nk);
      omegaD       *= toAngstSquare / ((double)nk);
      omegaOD      *= toAngstSquare / ((double)nk);
      omega(iSpin)  = omegaI + omegaD + omegaOD;

      cout << "|   Omega  = " << omega(iSpin) << endl;
      cout << "|   OmegaI = " << omegaI
           << ",  OmegaD = " << omegaD
           << ",  OmegaOD = " << omegaOD << endl;

   }  // :iSpin

   return omega;
}

//-----------------------------------------------------------------------------
//     gradient of spread functional
//-----------------------------------------------------------------------------
SxWannier::WFgradients
SxWannier::computeGradients (const SxWannier::WFTensor5 &Mkb,
                             const SxWannier::WFcenters &wrn) const
{
   // works out Eq. (57) in PRB B, 56, 12847 (1997)
   typedef SxArray<SxArray<SxDiracMat<Complex16> > >   Tensor4;

   Tensor4      Rkb(nk);       // Eq. (45)
   Tensor4      Tkb(nk);       // Eq. (51)
   Coord        bVec, rnVec;
   double       wb;
   SxComplex16  Mij, Mjj, q, rTilde;
   int          ik, ib, iSpin, i, j;

   WFgradients            dWk(nk);
   SxDiracMat<Complex16>  dW(nBands,nBands), rkb, tkb, U;
   double                 w;

   for (iSpin = 0; iSpin < nSpin; iSpin++)  {
      for (ik = 0; ik < nk; ik++)  {

         // --- resize the tensors w.r.t. the differences b
         Rkb(ik).resize (nb);
         Tkb(ik).resize (nb);

         for (ib = 0; ib < nb; ib++)  {
            bVec = kDiffs.bVecs(ib);

            // --- reformat the tensors w.r.t. the band indeces
            Rkb(ik)(ib).reformat      (nBands, nBands);
            Tkb(ik)(ib).reformat      (nBands, nBands);

            for (j = 0; j < nBands; j++)  {
               Mjj = Mkb(ik)(ib)(iSpin)(j,j);

               // --- avoid division by 0 in 'rTilde'
               if (Mjj.abs() < 1.e-20)  Mjj = 1.e-20;

               rnVec = wrn(iSpin)(j);
               q     = (lnC (Mjj)).im + (bVec * rnVec).sum ();

               for (i = 0; i < nBands; i++)  {
                  Mij              = Mkb(ik)(ib)(iSpin)(i,j);
                  Rkb(ik)(ib)(i,j) = Mij * Mjj.conj();
                  rTilde           = Mij / Mjj;
                  Tkb(ik)(ib)(i,j) = rTilde * q;
               }  // :i
            }  // :j
         }  // :ib
      }  // :ik

      // --- calculate dW, defined in Eq. (57)
      for (ik = 0; ik < nk; ik++)  {
         dWk(ik).resize (nSpin);

         dW.set (0.);
         w = 0.;
         for (ib = 0; ib < nb; ib++)  {
            wb   = kDiffs.bWeights(ib);
            rkb  = Rkb(ik)(ib);
            tkb  = Tkb(ik)(ib);
            dW  += wb * (sAnti(rkb) - sSymm(tkb));  // Eq. (57)
            w   += wb;
         }  // :ib
         
         dWk(ik)(iSpin) = dW * alpha / w;
         antiCheck (dWk(ik)(iSpin));
      }  // :ik
   }  // :iSpin

   return dWk;
}

void SxWannier::testComputeFunctions ()
{
   WFcenters  wrn = computeCenters (MkbInit);
   printCenters (wrn);

   WFspreads  spreads = computeSpreads (MkbInit, wrn);

   WFgradients g = computeGradients (MkbInit, wrn);
   cout << "g(0)(0) = " << g(0)(0) << endl;
   cout << "g(1)(0) = " << g(1)(0) << endl;
   cout << "g(7)(0) = " << g(7)(0) << endl;
}

void SxWannier::printCenters (const SxWannier::WFcenters &wrn) const
{
   cout << "| Wannier centers:\n|\n";

   int iSpin, i;

   for (iSpin = 0; iSpin < nSpin; iSpin++)  {
      if (nSpin > 1)  {
         if (iSpin == SPIN_UP)  cout << "|   Spin up:" << endl;
         else                   cout << "|   Spin down:" << endl;
      }

      for (i = 0; i < nBands; i++)  {
         sxprintf ("|   (%10.7f, %10.7f, %10.7f)\n",
                   wrn(iSpin)(i)(0), wrn(iSpin)(i)(1), wrn(iSpin)(i)(2));
      }
   }
   cout << "|" << endl;
}

//-----------------------------------------------------------------------------
//   m i n i m i s a t i o n  -  c y c l e
//-----------------------------------------------------------------------------
void SxWannier::computeMLWFs ()
{
   int            iteration;
   int            iSpin, ik, iq, ib;
   SxString       spin;
   WFcenters      wrn = computeCenters (MkbInit);
   WFspreads      spread;
   WFspreads      spreadOld = computeSpreads (MkbInit, wrn);
   WFgradients    dWk;
   WFTensor5      Mkb = MkbInit;  // TODO: danger of referencing?
   SxString       filenameU, nameU = "U";
   SxString       filenameD, nameD = "D";

   for (iSpin = 0; iSpin < nSpin; iSpin++)  {
      cout << SX_SEPARATOR;
      if (nSpin == 1)  {
         sxprintf ("| Minimisation of spread functional\n");
      }  else  {
         if (iSpin == SPIN_UP)  spin = "\"up\"";
         else                   spin = "\"down\"";
         sxprintf ("| Minimisation of spread functional "
                   "for spin channel %s\n", spin.ascii());
      }

      // --- initialise trafos
      Uk(iSpin).setDiagonals (1.);

      // --- start iterations
      for (iteration = 0; iteration < maxSteps; iteration++)  {
         cout << SX_SEPARATOR;
         sxprintf ("| Iteration %d\n|\n", iteration);

         // --- compute centers and gradients
         wrn = computeCenters (Mkb);
         printCenters (wrn);
         dWk = computeGradients (Mkb, wrn);

         // --- update unitary trafo matrices U, see Eq. (60)
         for (ik = 0; ik < nk; ik++)  {
            //uniCheck ( Uk(iSpin)(ik) );
            SxDiracMat<Complex16> M (Uk(iSpin)(ik) ^ mExp (dWk(ik)(iSpin)));
            Uk(iSpin).set (M, ik);
            uniCheck( Uk(iSpin)(ik) );
         }

         // --- update overlap matrices M, see Eq. (61)
         for (ik = 0; ik < nk; ik++)  {
            kDiffs.computeSurroundingKPoints (ik);

            for (ib = 0; ib < nb; ib++)  {
               iq = kDiffs.kSurIdx(ib);
               Mkb(ik)(ib)(iSpin) = Uk(iSpin)(ik).adjoint()
                                    ^ ( MkbInit(ik)(ib)(iSpin)
                                        ^ Uk(iSpin)(iq) );
            }
         }

         // --- calculate new spread
         spread = computeSpreads (Mkb, wrn);

         // --- convergence?
         if (fabs(spread(iSpin) - spreadOld(iSpin)) < dSpread)  {
            cout << SX_SEPARATOR;
            sxprintf ("| Maximally localised Wannier functions: ");
            if (nSpin == 1)  {
               sxprintf ("convergence reached.\n");
            }  else  {
               sxprintf("\n"
                        "Convergence reached for spin channel %s.",
                        spin.ascii());
            }
            break;
         }

         // --- max. number of steps exceeded
         if (iteration >= maxSteps - 1)  {
            cout << SX_SEPARATOR;
            sxprintf ("| WARNING:  Maximally localised Wannier functions:\n"
                      "|           Maximum number of steps exceeded.");
            if (nSpin == 1)  {
               sxprintf ("\n");
            }  else  {
               sxprintf (" for spin channel %s.\n", spin.ascii());
            }
            sxprintf("|           Convergence not yet reached.\n");
            break;
         }

         spreadOld = spread;
      }  // :iteration

      // --- write the altogether trafo tensor D = P * U to output file
      if (nSpin == 1)  {
         D(0) = P * Uk(0);
         filenameU = nameU+".sxb";
         filenameD = nameD+".sxb";
      }  else  {
         D(iSpin) = P * Uk(iSpin);
         if (iSpin == SPIN_UP)  {
            filenameU = nameU+"-up.sxb";
            filenameD = nameD+"-up.sxb";
         }  else  {
            filenameU = nameU+"-down.sxb";
            filenameD = nameD+"-down.sxb";
         }
      }

      Uk(iSpin).writeDirTensor3 (filenameU);
      D(iSpin).writeDirTensor3 (filenameD);

   }  // :iSpin
}

//-----------------------------------------------------------------------------
//   Wannier functions in real space
//-----------------------------------------------------------------------------
SxWannier::WFwannier SxWannier::computeMLWFsRealSpace ()
{
   int             ik, iSpin, n, m;
   SxGBasis       *gkPtr;
   SxRBasis        q (meshDim, cell);
   double          kWeight;
   Coord           k, R = cell.relToCar (Coord (translation));
   int             ir, nr = repetition.product() * meshDim.product();
   WFwannier       wannierFunctions(nSpin);
   PsiR            colUk, phaseFactor(nr);
   PsiG            uN, uM;
   PsiR            uNR, uNRRepeated;
   SxArray<Coord>  rMesh = SxGauss::computeMeshPoints (meshDim,
                                                       repetition, cell);
   double          phase;
   SxComplex16     norm;
   const SxPW     &uInit = *uInitPtr;

   cout << SX_SEPARATOR;
   sxprintf ("| Calculation of Wannier functions in real space ... \n");
   cout.flush ();

   /*
   SxPW uInit;
   //SxDiracVec<Complex16> w(nr);

   SxString ukmn = "ukmn.sxb";
   Uk(0).readUkmn (ukmn);
   Uk(0).writeDirTensor3 ("readInU.sxb");
   try  {
      SxBinIO io ("uInit.sxb", SxBinIO::BINARY_READ_ONLY);
      uInit.read (io);
      io.close ();
   }  catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }
   for (ik = 0; ik < nk; ik++)  {
      for (int i = 0; i < nBands; i++)  {
         sxprintf ("Summe: ik = %d, i = %d --> s = (%g,%g)\n",
                 ik, i,
                 uInit(i,0,ik).sum().re,
                 uInit(i,0,ik).sum().im);
      }
   }
   */

   for (iSpin = 0; iSpin < nSpin; iSpin++)  {
      wannierFunctions(iSpin).resize (nBands);

      for (n = 0; n < nBands; n++)  {
         wannierFunctions(iSpin)(n).resize (nr);
         wannierFunctions(iSpin)(n).set (0.);
      }

      for (ik = 0; ik < nk; ik++)  {
         gkPtr   = &( gkSet(ik) );
         k       = kDiffs.kVec(ik);
         kWeight = kDiffs.weights(ik);
         cout << "ik = " << ik << " --> k = " << k << endl;

         // -------------------------------------------------------------
         //     phase: (r - R) k = kr - kR
         // -------------------------------------------------------------
         //     Explanation (see Eqs. (1), (3), (59)):
         //     The variable 'uN' refers to |u_nk> of Eq. (59).
         //     To calculate the WANNIER functions, Eq. (1) needs the
         //     BLOCH functions |psi_nk>, which are given by Eq. (2)
         //     with
         //           psi_nk (r) = exp{+ikr} u_nk (r)  (in real space).
         //
         //     The second, negative summand in 'phase' specifies,
         //     which WANNIER function
         //                             w_n (r-R) = w_nR (r) = < r | nR >
         //
         //     shall be computed (Eq. (1)). Therefore, R is specified in
         //     the input file by 'translation'.
         for (ir = 0; ir < nr; ir++)  {
            phase           = (rMesh(ir) - R) ^ k;
            phaseFactor(ir) = SxComplex16 (cos(phase), sin(phase));
         }

         for (n = 0; n < nBands; n++)  {
            colUk = Uk(iSpin)(ik).colRef(n);
            uN.resize (gkPtr->ng);

            uN.set (0.);
            for (m = 0; m < nBands; m++)  {
               // uM = < G | u_mk > in Eq. (59)
               uM  = uInit(m,iSpin,ik);
               // uN = < G | u_nk > in Eq. (59)
               uN += colUk(m) * uM;
            }

            // --- test norm of |u_nk>
            norm = (uN ^ uN).chop();
            if (fabs(norm.re-1.) > 1.e-4 || norm.im > 1.e-4)  {
               sxprintf ("| WARNING:  Norm of Wannier function differs from "
                         "one:  |u_nk|^2 = (%g,%g)\n", norm.re, norm.im);
            }

            // --- real space representation of |u_nk>
            uN.setBasis (gkPtr);
            uNR = ( q | uN );

            // --- < q | u_nk > is a periodic function in real space
            uNRRepeated = repeatFunction (uNR, meshDim, repetition);

            // --- The multiplication of < r | u_nk > with the phase factor
            //     exp{-ik(r-R)} performes
            //
            //        (1)  < r | u_nk > --> psi_nk, according to Eq. (3), and
            //        (2)  gives the factor exp{-ikR} of Eq. (1).
            //
            //     Multiplying by the weights and summation over the k-points
            //     performs the integration in Eq. (1), but - presumably -
            //     without the prefactor V / (2*pi)^3.
            wannierFunctions(iSpin)(n) += kWeight
                                          * (uNRRepeated * phaseFactor);
         }  // :n
      }  // :ik
   }  // :iSpin

   cout << "done." << endl;

   return wannierFunctions;
}

void SxWannier::writeMLWFs ()
{
   Coord               range1 = cell.col(0) * repetition(0);
   Coord               range2 = cell.col(1) * repetition(1);
   Coord               range3 = cell.col(2) * repetition(2);
   SxCell              range (range1, range2, range3);
   RelVec              dimOut = meshDim * repetition;
   int                 iSpin, i, ir, nr = dimOut.product();
   SxDiracVec<Double>  wfAbs(nr);
   double              x, y, imagPart;
   SxString            filename, name = "WF";

   WFwannier wf = computeMLWFsRealSpace ();

   for (iSpin = 0; iSpin < nSpin; iSpin++)  {
      for (i = 0; i < nBands; i++)  {

         // --- get absolute of Wannier functions
         imagPart = 0.;
         for (ir = 0; ir < nr; ir++)  {
            x = wf(iSpin)(i)(ir).re;
            y = wf(iSpin)(i)(ir).im;
            wfAbs(ir) = sqrt (x*x + y*y);
            if (fabs(y) > imagPart)  imagPart = y;
         }

         // --- dump warning, if large imaginary parts has been encountered
         if (fabs(imagPart) > 1.e-3)  {
            sxprintf ("| WARNING:  largest imaginary part in WF %d: ", i+1);
            sxprintf ("%.5g\n", imagPart);
         }

         // --- write functions to output file(s)
         if (nSpin == 1)  {
            filename = name+"-"+(i+1)+".sxb";
         }  else  {
            if (iSpin == SPIN_UP)  filename = name+"-"+(i+1)+"-up.sxb";
            else                   filename = name+"-"+(i+1)+"-down.sxb";
         }

         try  {
            SxBinIO io (filename, SxBinIO::BINARY_WRITE_ONLY);
            io.writeMesh (wfAbs, range, dimOut);
            io.setMode (SxBinIO::WRITE_DATA);
            io.writeMesh (wfAbs, range, dimOut);
            io.close ();
         }  catch (SxException e)  {
            e.print ();
            SX_QUIT;
         }
      }  // :i
   }  // :iSpin
}

//-----------------------------------------------------------------------------
//   logarithm of complex number (except for additive 2pi*i*n, n elem. N)
//-----------------------------------------------------------------------------
SxComplex16 SxWannier::lnC (SxComplex16 &arg) const
{
   double       x, y, r, phi;
   SxComplex16  z;

   x = arg.re;
   y = arg.im;

   if (x*x < 1.e-40)  {
      cout << SX_SEPARATOR;
      sxprintf ("| WARNING:  logarithm of pure imaginary number encountered!\n");
      SX_EXIT;
      z.re = 0;
      z.im = -1000.; //TODO must be some other value, hv been too lazy to think
      return z;
   }  else  {
      r    = sqrt ( x*x + y*y );
      phi  = atan (y / x);
      if (x < 0. && y > 0.)  phi += PI;
      if (x < 0. && y < 0.)  phi -= PI;
      z.re = log (r);
      z.im = phi;
      return z;
   }
}

//-----------------------------------------------------------------------------
//   exponential function of complex number
//-----------------------------------------------------------------------------
SxComplex16 SxWannier::expC (const SxComplex16 &z) const
{
   double       x, y, p;
   SxComplex16  h;

   x = z.re;
   y = z.im;

   p    = exp(x);
   h.re = p * cos(y);
   h.im = p * sin(y);

   return h;
}

//-----------------------------------------------------------------------------
//   exp( B )  where  B  is an anti-Hermitian matrix
//-----------------------------------------------------------------------------
SxDiracMat<Complex16> SxWannier::mExp (const SxDiracMat<Complex16> &B) const
{
   // See also R. Wuest, Hoehere Mathematik fuer Physiker, Bd. II, p. ????,
   // and Ilan Schnell, Diplomarbeit, Universitaet Bremen (1999), page 47.

   SX_CHECK (B.nRows() == B.nCols(), B.nRows(), B.nCols());

   SxDiracMat<Complex16>  H, C;
   SxDiracMat<Complex16>::Eigensystem eigSysH;
   SxDiracVec<Complex16>  eigValsH;
   SxDiracMat<Complex16>  eigVecsH, evInvH;

   // --- get dimension of matrix
   int dim = int(B.nRows());
   C.reformat (dim,dim);
   C.set (0.);

   // --- H := iB  ==>  Is B anti-Hermitian, then H is Hermitian.
   H = I * B;

   // --- diagonalize H
   eigSysH  = H.eigensystem (true);  // 'true' to sort the EWs
   eigValsH = eigSysH.vals;
   eigVecsH = eigSysH.vecs;
   evInvH   = eigVecsH.inverse();

   // --- create diagonal matrix, considering exp{B} = exp{-i*iB} = exp{-iH}
//   C = C.identity (-I * eigValsH, dim);
   for (int i = 0; i < dim; i++)  C(i,i) = expC (-I * eigValsH(i));

   return ( (eigVecsH ^ C) ^ evInvH );
}

//-----------------------------------------------------------------------------
//   unitary check
//-----------------------------------------------------------------------------
void SxWannier::uniCheck (const SxDiracMat<Complex16> &U) const
{
   // --- Checks if the product U*U^+ = E equals the identity,
   //     works for quadratic matrices only.

   SX_CHECK (U.nRows() == U.nCols(), U.nRows(), U.nCols());

   int                    i, j, d = int(U.nRows());
   SxDiracMat<Complex16>  eins(d,d);
   SxDiracMat<Complex16>  adj;
   SxDiracMat<Complex16>  p;
   double                 eiRe, eiIm, pRe, pIm, x, y;
   double                 epsilon = 1.e-2;

   eins = eins.identity();
   adj  = U.adjoint();
   p    = (U ^ adj);

   for (i = 0; i < d; i++)  {
      for (j = 0; j < d; j++)  {
         eiRe = eins.row(i)(j).re;
         eiIm = eins.row(i)(j).im;
         pRe  = p.row(i)(j).re;
         pIm  = p.row(i)(j).im;
         x    = sqrt ( pow ((eiRe - pRe), 2.) );
         y    = sqrt ( pow ((eiIm - pIm), 2.) );

         if (x > epsilon || y > epsilon)  {
            sxprintf ("| Check for unitarity failed. Sorry.\n");
            sxprintf ("\tepsilon = %g\n", epsilon);
            sxprintf ("\tRe{ U*U+_(%d,%d) } - Re{ E_(%d,%d) = %g\n", i,j,i,j,x);
            sxprintf ("\tIm{ U*U+_(%d,%d) } - Im{ E_(%d,%d) = %g\n", i,j,i,j,y);
            SX_QUIT;
         }
      }  // :j
   }  // :i
   //sxprintf ("  uniCheck: OK\n");
}

//-----------------------------------------------------------------------------
//   anti-Hermitian check 
//-----------------------------------------------------------------------------
void SxWannier::antiCheck (const SxDiracMat<Complex16> &W) const
{
   // --- Checks if a matrix is anti-Hermitian or not,
   //     works for quadratic matrices only.

   SX_CHECK (W.nRows() == W.nCols(), W.nRows(), W.nCols());

   int          i, j, d = int(W.nRows());
   SxComplex16  aIJ, aJI;
   double       aIJre, aIJim, aJIre, aJIim, x, y;
   double       epsQuit = 1.e-2, epsWarn = 1.e-5;

   for (i = 0; i < d; i++)  {
      for (j = 0; j < d; j++)  {
         aIJ   = W.row(i)(j);
         aJI   = W.row(j)(i);
         aIJre = aIJ.re;
         aIJim = aIJ.im;
         aJIre = aJI.re;
         aJIim = aJI.im;
         x     = fabs (aIJre + aJIre);
         y     = fabs (aIJim - aJIim);

         if (x > epsQuit || y > epsQuit)  {
            sxprintf ("| Check for anti-Hermitecity failed. Sorry.\n");
            sxprintf ("\tepsilon = %g\n", epsQuit);
            sxprintf ("\tRe{ a_(%d,%d) } + Re{ a_(%d,%d) } = %g\n", i,j,j,i,x);
            sxprintf ("\tIm{ a_(%d,%d) } - Im{ a_(%d,%d) } = %g\n", i,j,j,i,y);
            SX_QUIT;
         }

         if (x > epsWarn || y > epsWarn)  {
            cout << SX_SEPARATOR;
            sxprintf ("| WARNING: matrix not really anti-Hermitian. ");
            sxprintf ("Be careful!\n");
            sxprintf ("| epsilon = %g\n", epsWarn);
            sxprintf ("\tRe{ a_(%d,%d) } + Re{ a_(%d,%d) } = %g\n", i,j,j,i,x);
            sxprintf ("\tIm{ a_(%d,%d) } - Im{ a_(%d,%d) } = %g\n", i,j,j,i,y);
         }
      }  // :j
   }  // :i
   //sxprintf ("   antiCheck: OK\n");
}

//-----------------------------------------------------------------------------
//   operators  A[ B ] := (B - B ) / 2 ,    S[ B ] := (B + B ) / 2i
//-----------------------------------------------------------------------------
SxDiracMat<Complex16> SxWannier::sAnti (const SxDiracMat<Complex16> &B) const
{
   return (B - B.adjoint()) / 2.;
}

SxDiracMat<Complex16> SxWannier::sSymm (const SxDiracMat<Complex16> &B) const
{
   return (B + B.adjoint()) / (2.*I);
}

//-----------------------------------------------------------------------------
//   write WFTensor5 to output file (ASCII - for debugging)
//-----------------------------------------------------------------------------
void
SxWannier::writeMkb (const SxWannier::WFTensor5 &Mkb, const SxString &filename)
{
   SX_CHECK (Mkb.getSize() == nk, Mkb.getSize(), nk);
   SX_CHECK (Mkb(0).getSize() == nb, Mkb(0).getSize(), nb);
   SX_CHECK (Mkb(0)(0).getSize() == nSpin, Mkb(0)(0).getSize(), nSpin);
   SX_CHECK (Mkb(0)(0)(0).nRows() == nBands, Mkb(0)(0)(0).nRows(), nBands);
   SX_CHECK (Mkb(0)(0)(0).nCols() == nBands, Mkb(0)(0)(0).nCols(), nBands);

   FILE *fp = NULL;

   fp = fopen (filename.ascii(), "w");
   if (!fp)  {
      sxprintf ("Can't open file \"%s\". Sorry.\n", filename.ascii());
      SX_QUIT;
   }

   // --- create separator
   SxString sep = "";
   for (int w = 0; w < 160; w++)  sep += "-";
   sep += "\n";

   // --- print header line
   fprintf (fp, "\n");

   // --- print data
   int          ik, ib, iSpin, i, j;
   SxComplex16  c;

   for (iSpin = 0; iSpin < nSpin; iSpin++)  {
      for (ik = 0; ik < nk; ik++)  {
kDiffs.computeSurroundingKPoints (ik);

         for (ib = 0; ib < nb; ib++)  {
            fprintf (fp, "%s", sep.ascii());
            //fprintf (fp, "iSpin = %d,  M(ik = %d, ib = %d)\n", iSpin, ik, ib);
            fprintf (fp, "iSpin = %d,  M(ik = %d, ib = %d)  "
                         "-->  jk = %d  -->  ",
                         iSpin, ik, ib, kDiffs.kSurIdx(ib));
            fprintf (fp, "dG = {%d,%d,%d}\n",
                         kDiffs.deltaGRel(ib)(0),
                         kDiffs.deltaGRel(ib)(1),
                         kDiffs.deltaGRel(ib)(2));
            fprintf (fp, "%s", sep.ascii());

            for (i = 0; i < nBands; i++)  {
               for (j = 0; j < nBands; j++)  {
                  c = Mkb(ik)(ib)(iSpin)(i,j);
                  fprintf (fp, "|(%g, %g)| = %g\t\t", c.re, c.im, c.abs());
               }  // :j
               fprintf (fp, "\n\n");
            }  // :i
            fprintf (fp, "\n\n\n");
         }  // :ib
      }  // :ik
   }  //:iSpin

   fclose (fp);
}

//-----------------------------------------------------------------------------
//   repeatition of a function
//-----------------------------------------------------------------------------
SxDiracVec<Complex16>
SxWannier::repeatFunction (const SxDiracVec<Complex16> &function,
                           const RelVec &dim_, const RelVec &repeat)
{
   SX_CHECK (function.getSize() == dim_.product(),
             function.getSize(), dim_.product());
   SX_CHECK  (dim_(0) > 0, dim_(0));
   SX_CHECK  (dim_(1) > 0, dim_(1));
   SX_CHECK  (dim_(2) > 0, dim_(2));
   SX_CHECK  (repeat(0) > 0, repeat(0));
   SX_CHECK  (repeat(1) > 0, repeat(1));
   SX_CHECK  (repeat(2) > 0, repeat(2));

   if (repeat == RelVec (1,1,1))  return function;

   SxMesh3D               dim (dim_);
   SxMesh3D               dimOut = dim * repeat;
   SxDiracVec<Complex16>  functionOut( dimOut.product() );
   SxComplex16            val;
   int                    x, y, z, i, j, k, idx, xOff, yOff, zOff;

   // --- access points of the original mesh
   for (z = 0; z < dim(2); z++)  {
      for (y = 0; y < dim(1); y++)  {
         for (x = 0; x < dim(0); x++)  {

            idx = (int)dim.getMeshIdx (x, y, z, SxMesh3D::Positive);
            val = function (idx);

            // --- translation to neighboured unit cells
            for (k = 0; k < repeat(2); k++)  {
               zOff = z + k * dim(2);
               for (j = 0; j < repeat(1); j++)  {
                  yOff = y + j * dim(1);
                  for (i = 0; i < repeat(0); i++)  {
                     xOff = x + i * dim(0);

                     idx = (int)dimOut.getMeshIdx (xOff, yOff, zOff,
                                                   SxMesh3D::Positive);
                     functionOut(idx) = val;
                  }
               }
            }

         }
      }
   }

   return functionOut;
}


#else /* SX_STANDALONE */


#include <SxCLI.h>
#include <SPHInX.h>

int main (int argc, char **argv)
{
   // --- start timers

   SxCLI cli (argc, argv);

   SxString input = cli.option ("-i|-w|--input", "file",
      "input file for computation of Wannier functions\n\n"
      "Note: This is not the usual S/PHI/nX input file,\n"
      "      but a special one for the Wannier functions.")
      .toString ("wannier.sx");

   cli.authors = "Matthias Wahn, wahn@fhi-berlin.mpg.de";
   cli.version ("1.0");
   cli.finalize ();

   // --- dump headline
   cout << SX_SEPARATOR;
   cout << "|               C R E A T E   W A N N I E R   F U N C T I O N S"
        << endl;
   cout << SX_SEPARATOR;

   // --- parse input file
   SxParser parser;
   SxParser::Table table = parser.read (input);

   // --- instantiate Wannier function object
   SxWannier w (&*table);

   // --- dump input parameters
   w.print ();

   // --- initialise overlap matrices
   w.initOverlapMatrices ();

   // --- compute GWFs
   w.computeMLWFs ();

   // --- write a real space version of the GWFs to .sxb files
   w.writeMLWFs ();

   return 0;
}


#endif /* SX_STANDALONE */
