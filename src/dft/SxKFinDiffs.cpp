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

#include <SxKFinDiffs.h>

SxKFinDiffs::SxKFinDiffs () : SxKPoints ()
{
   foldingMP = 0;
   folding = 0.;
   nbVecs = -1;
}

SxKFinDiffs::SxKFinDiffs (const SxCell &cell_, const SxGkBasis &gkBasis)
   : SxKPoints (gkBasis)
{
   SX_CHECK (nk > 0, nk);
   SX_CHECK (cell_.volume > 0., cell_.volume);

   cell = cell_;  // TODO: Ask, if this is allowed!
   recCell = cell_.getReciprocalCell();

   // --- get folding from unsymmetrised MP mesh
   foldingMP = SxKFinDiffs::getFoldingMP (true);  // include Gamma point
   double f1 = double (foldingMP(0));
   double f2 = double (foldingMP(1));
   double f3 = double (foldingMP(2));
   folding = SxVector3<Double> (f1, f2, f3);
//   folding = SxVector3<Double> (foldingMP);  // TODO: works as well, probably

   // --- get type of the lattice structure
   lattice = SxCell::CellType (cell);
   lattice.print ();
   typedef SxCell::CellType LatType;
   LatType::Type type = lattice.getType ();

   switch (type)  {
      case LatType::FaceCenteredCubic : initFCC ();
                                        break;
      case LatType::Hexagonal         : initWurtzite ();
                                        break;
      default : sxprintf ("\n"
                          "Error:  This add-on can't deal yet with the %s\n"
                          "        lattice structure. Please choose an other "
                          "system or volunteer\n"
                          "        for the required implementations. Sorry.\n",
                          lattice.getTypeName().ascii());
                SX_QUIT;
                break;
   }

   computeKVecsRel ();
   //checkKIdcs ();

   cout << "" << endl;
}

void SxKFinDiffs::initFCC ()
{
   // --- check MP mesh -- also in non-debug cases
   if ( !(foldingMP(0) == foldingMP(1) && foldingMP(1) == foldingMP(2)) )  {
      sxprintf ("Error:  In case of %s, a uniform MP mesh (a x a x a)\n"
                "        has to be used. The detected mesh has the "
                "shape %d x %d x %d. Sorry.\n",
                lattice.getTypeName().ascii(),
                foldingMP(0), foldingMP(1), foldingMP(2));
      SX_QUIT;
   }

   // --- TODO: check cell vectors

   // --- set up b-vectors in Cartesian coordinates
   double f = double (foldingMP(0));
   double normFac = 2. * PI / (lattice.getA() * f);
   nbVecs = 8;
   bVecs.resize (nbVecs);
   bVecs(0) = SxVector3<Double> ( 1, 1, 1) * normFac;
   bVecs(1) = SxVector3<Double> (-1, 1, 1) * normFac;
   bVecs(2) = SxVector3<Double> ( 1,-1, 1) * normFac;
   bVecs(3) = SxVector3<Double> ( 1, 1,-1) * normFac;
   bVecs(4) = SxVector3<Double> (-1,-1, 1) * normFac;
   bVecs(5) = SxVector3<Double> (-1, 1,-1) * normFac;
   bVecs(6) = SxVector3<Double> ( 1,-1,-1) * normFac;
   bVecs(7) = SxVector3<Double> (-1,-1,-1) * normFac;

   // --- set up b-vectors in relative coordinates
   Coord  relVecDouble;
   bVecsRel.resize (nbVecs);
   for (int b = 0; b < nbVecs; b++)  {
      relVecDouble = f * recCell.carToRel (bVecs(b));
      bVecsRel(b)(0) = int ( lround(relVecDouble(0)) );
      bVecsRel(b)(1) = int ( lround(relVecDouble(1)) );
      bVecsRel(b)(2) = int ( lround(relVecDouble(2)) );

      sxprintf ("bRel( %d ) = ", b);  cout << bVecsRel(b) << endl;
   }

   // --- check relative b-vectors
   int bx, by, bz;
   for (int b = 0; b < nbVecs; b++)  {
      bx = bVecsRel(b)(0);
      by = bVecsRel(b)(1);
      bz = bVecsRel(b)(2);

      if (    (bx != -1 && bx != 0 && bx != 1)
           || (by != -1 && by != 0 && by != 1)
           || (bz != -1 && bz != 0 && bz != 1) )
      {
         sxprintf ("\n"
                   "Bug in program: Finite difference vectors (b-vectors) "
                   "not correctly\n"
                   "                calculated.\n"
                   "                Please report this error to "
                   "the S/PHI/nX team.\n"
                   "                Sorry for that.\n");
         SX_QUIT;
      }
   }

   // --- weights of the b-vectors
   double fccWeight = 1. / (nbVecs * normFac * normFac);  // = 3 / (nb * b^2)
   bWeights.resize (nbVecs);
   for (int b = 0; b < nbVecs; b++)  bWeights(b) = fccWeight;

   // --- resize arrays of surrounding k-points and G shift vectors
   kSurIdx.resize (nbVecs);
   deltaGCart.resize (nbVecs);
   deltaGRel.resize (nbVecs);
}

void SxKFinDiffs::initWurtzite ()
{
   // --- check MP mesh -- also in non-debug cases
   if (foldingMP(0) != foldingMP(1))  {
      sxprintf ("Error:  In case of a hexagonal structure, the MP mesh must "
                "be of the following\n"
                "        form: a x a x c\n"
                "        The detected mesh, however has the shape "
                "%d x %d x %d. Sorry.\n",
//                lattice.getTypeName().ascii(),
                foldingMP(0), foldingMP(1), foldingMP(2));
      SX_QUIT;
   }

   // --- TODO: check cell vectors

   // --- set up b-vectors in Cartesian coordinates
   double fxy = double (foldingMP(0));
   double fz  = double (foldingMP(2));
   double normFacXY = 4. * PI / SQRT3 / (lattice.getA() * fxy);
   double normFacZ  = 2. * PI / (lattice.getC() * fz);
   nbVecs = 8;
   bVecs.resize (nbVecs);
   bVecs(0) = SxVector3<Double> ( SQRT3/2.,  0.5,  0.0) * normFacXY;
   bVecs(1) = SxVector3<Double> (   0.0,     1.0,  0.0) * normFacXY;
   bVecs(2) = SxVector3<Double> (-SQRT3/2.,  0.5,  0.0) * normFacXY;
   bVecs(3) = SxVector3<Double> (-SQRT3/2., -0.5,  0.0) * normFacXY;
   bVecs(4) = SxVector3<Double> (   0.0,    -1.0,  0.0) * normFacXY;
   bVecs(5) = SxVector3<Double> ( SQRT3/2., -0.5,  0.0) * normFacXY;
   bVecs(6) = SxVector3<Double> (   0.0,     0.0, -1.0) * normFacZ;
   bVecs(7) = SxVector3<Double> (   0.0,     0.0,  1.0) * normFacZ;

   // --- set up b-vectors in relative coordinates
   Coord relVecDouble;
   bVecsRel.resize (nbVecs);
   for (int b = 0; b < nbVecs; b++)  {
      relVecDouble = folding * recCell.carToRel (bVecs(b));
      bVecsRel(b)(0) = int ( lround(relVecDouble(0)) );
      bVecsRel(b)(1) = int ( lround(relVecDouble(1)) );
      bVecsRel(b)(2) = int ( lround(relVecDouble(2)) );

      sxprintf ("bRel( %d ) = ", b);  cout << bVecsRel(b) << endl;
   }

   // --- check relative b-vectors
   int bx, by, bz;
   for (int b = 0; b < nbVecs; b++)  {
      bx = bVecsRel(b)(0);
      by = bVecsRel(b)(1);
      bz = bVecsRel(b)(2);

      if (    (bx != -1 && bx != 0 && bx != 1)
           || (by != -1 && by != 0 && by != 1)
           || (bz != -1 && bz != 0 && bz != 1) )
      {
         sxprintf ("\n"
                   "Bug in program: Finite difference vectors (b-vectors) "
                   "not correctly\n"
                   "                calculated.\n"
                   "                Please report this error to "
                   "the S/PHI/nX team.\n"
                   "                Sorry for that.\n");
         SX_QUIT;
      }
   }

   // --- weights of the b-vectors
   double hexWeightXY = 1. / (3. * normFacXY * normFacXY);
   double hexWeightZ  = 1. / (2. * normFacZ  * normFacZ);
   bWeights.resize (nbVecs);
   for (int b = 0; b < nbVecs; b++)
      bWeights(b) = (b < 6) ? hexWeightXY : hexWeightZ;

   // --- resize arrays of surrounding k-points and G shift vectors
   kSurIdx.resize (nbVecs);
   deltaGCart.resize (nbVecs);
   deltaGRel.resize (nbVecs);
}

int SxKFinDiffs::getNb () const
{
   return nbVecs;
}

void SxKFinDiffs::computeSurroundingKPoints (int ik)
{
   SX_CHECK (ik >= 0 && ik < nk, ik, nk);
   SX_CHECK (kSurIdx.getSize() == nbVecs, kSurIdx.getSize(), nbVecs);
   SX_CHECK (deltaGCart.getSize() == nbVecs, deltaGCart.getSize(), nbVecs);
   SX_CHECK (deltaGRel.getSize() == nbVecs, deltaGRel.getSize(), nbVecs);

   RelVec qRel, kRel = kVecsRel(ik), gRel;

   for (int b = 0; b < nbVecs; b++)  {
      gRel = 0;
      qRel = kRel + bVecsRel(b);
      /*
      if (b == 5)  {
      cout << "gRel:     " << gRel << endl;
      cout << "kRel:     " << kRel << endl;
      cout << "bVecsRel: " << bVecsRel(b) << endl;
      cout << "qRel:     " << qRel << endl;
      }
      */

      // --- if necessary, map q -> q'
      if (qRel(0) < 0)  {
         SX_CHECK (qRel(0) == -1, qRel(0));
         qRel(0) += foldingMP(0);
         gRel(0)  = -1;
      }  else if (qRel(0) >= foldingMP(0))  {
         SX_CHECK (qRel(0) == foldingMP(0), qRel(0), foldingMP(0));
         qRel(0) -= foldingMP(0);
         gRel(0)  = 1;
      }

      if (qRel(1) < 0)  {
         SX_CHECK (qRel(1) == -1, qRel(1));
         qRel(1) += foldingMP(1);
         gRel(1)  = -1;
      }  else if (qRel(1) >= foldingMP(1))  {
         SX_CHECK (qRel(1) == foldingMP(1), qRel(1), foldingMP(1));
         qRel(1) -= foldingMP(1);
         gRel(1)  = 1;
      }

      if (qRel(2) < 0)  {
         SX_CHECK (qRel(2) == -1, qRel(2));
         qRel(2) += foldingMP(2);
         gRel(2)  = -1;
      }  else if (qRel(2) >= foldingMP(2))  {
         SX_CHECK (qRel(2) == foldingMP(2), qRel(2), foldingMP(2));
         qRel(2) -= foldingMP(2);
         gRel(2)  = 1;
      }
      /*
      if (b == 5)  {
      cout << "qPrRel:   " << qRel << endl;
      cout << "gRel:     " << gRel << endl;
      }
      */

      // --- get index of q'
      kSurIdx(b) = getKPointIdx (qRel);

      // --- get G shift vector in Cartesian and relative coordinates
      //     TODO: invert minus signs above and drop them here ...
      deltaGCart(b) = -recCell.relToCar (Coord (gRel) );
      deltaGRel(b)  = -RelVec (gRel);
   }
}

int SxKFinDiffs::getKPointIdx (const RelVec &kRel)
{
   SX_CHECK (foldingMP.product() == nk, foldingMP.product(), nk);
   SX_CHECK (kRel(0) >= 0 && kRel(0) < foldingMP(0), kRel(0), foldingMP(0));
   SX_CHECK (kRel(1) >= 0 && kRel(1) < foldingMP(1), kRel(1), foldingMP(1));
   SX_CHECK (kRel(2) >= 0 && kRel(2) < foldingMP(2), kRel(2), foldingMP(2));

   return (kRel(0) * foldingMP(1) + kRel(1)) * foldingMP(2) + kRel(2);
}

void SxKFinDiffs::checkKIdcs ()
{
   RelVec kRel;
   int    kIdx;

   for (int ik = 0; ik < nk; ik++)  {
      cout << "idx = " << ik << "  -->  ";
      kRel = kVecsRel(ik);
      cout << "kRel = " << kRel << "  -->  ";
      kIdx = getKPointIdx (kRel);
      cout << "idx = " << kIdx << endl;
      if (ik != kIdx)  {
         sxprintf ("Error:  SxKFinDiffs::getKPointIdx failed in finding "
                   "correct index. Sorry.\n");
         SX_QUIT;
      }
   }
}

void SxKFinDiffs::checkKSur (int ik)
{
   SX_CHECK (ik >= 0 && ik < nk, ik, nk);

   int   qIdx;
   Coord k = kVec(ik), q, diff;  //cout << "k    =      " << k << endl;

   for (int b = 0/*5*/; b < nbVecs/*6*/; b++)  {
      qIdx  = kSurIdx(b);  //cout << "qIdx =      " << qIdx << endl;
      q     = kVec(qIdx);  //cout << "q    =      " << q << endl;
                           //cout << "qRel =      " << kVecsRel(b) << endl;
      //cout << "dGRel =  " << recCell.carToRel (deltaG(b)) << endl;
      q    -= deltaGCart(b);   //cout << "deltaG(b) = " << deltaG(b) << endl;
      diff  = q - k;

      cout << "b( " << b << " ) = " << bVecs(b) << "  -->  "
           << "diff = " << diff << endl;

      if ( (bVecs(b) - diff).norm() > 1.e-4 )  SX_QUIT;
   }
   sxprintf ("ik = %d -> q-points ok.\n", ik);
}

void SxKFinDiffs::checkKSur ()
{
   for (int ik = 0; ik < nk; ik++)  {
      computeSurroundingKPoints (ik);
      cout << "kSurIdx: " << kSurIdx << endl;
      cout << "deltaGCart:  " << deltaGCart << endl;
      checkKSur (ik);
   }
   sxprintf ("q-points for all k-points correct.\n");
}

SxVector3<Int> SxKFinDiffs::getFoldingMP (bool gamma)
{
   if (gamma == false)  SX_EXIT;  // not yet tested
   cout << SX_SEPARATOR;

   // --- offset in Cartesian coordinates
   Coord offset = kVec(0);

   // --- check for Gamma-point
   if (gamma == true)  {
      if (fabs(offset(0)) + fabs(offset(1)) + fabs(offset(2)) > 1.e-4)  {
         sxprintf ("\n"
                   "Error:  MP mesh doesn't contain Gamma-point. Sorry.\n");
         SX_QUIT;
      }
   }

   // --- get relative coordinates
   SxArray<Coord> kVecRel;
   kVecRel.resize (nk);
   for (int ik = 0; ik < nk; ik++)
      kVecRel(ik) = recCell.carToRel (kVec(ik) - offset);

   // --- get folding
   SxVector3<Int> i, fold(1);

   i(2) = 1;
   while (fabs(kVecRel(i(2))(2) - kVecRel(0)(2)) > 1.e-4)  {
      i(2)++;
      fold(2)++;
      SX_CHECK (i(2) <= nk, i(2), nk);
   }

   i(1) = i(2);
   while (fabs(kVecRel(i(1))(1) - kVecRel(0)(1)) > 1.e-4)  {
      i(1) += i(2);
      fold(1)++;
      SX_CHECK (i(1) <= nk, i(1), nk);
   }

   i(0) = i(1);
   while (i(0) < nk && fabs(kVecRel(i(0))(0) - kVecRel(0)(0)) > 1.e-4)  {
      i(0) += i(1);
      fold(0)++;
      SX_CHECK (i(0) <= nk, i(0), nk);
   }

   cout << "| MP folding: " << fold << endl;
   cout << "| offset:     " << offset * fold << endl;  // offset in input file

   // --- may happen in non-debug cases
   if (fold.product() != nk)  {
      sxprintf ("\n"
                "Error:  Can't find MP folding. Sorry. Maybe, one of the "
                "following\n"
                "        reasons applies:\n\n"
                "   (1)  The MP mesh contains more than one k-point per tile.\n"
                "   (2)  The mesh is too dense.\n"
                "   (3)  The k-points in the input file have been generated by "
                "a different\n"
                "        unit cell (likely).\n"
                "   (4)  You try to use a symmetry-reduced k-point mesh "
                "(most likely).\n");
      SX_QUIT;
   }

   return fold;
}

void SxKFinDiffs::computeKVecsRel ()
{
   Coord relVecDouble;
   SX_CHECK (folding.product() > 0.5, folding.product());
   kVecsRel.resize (nk);

   for (int ik = 0; ik < nk; ik++)  {
      relVecDouble = folding * recCell.carToRel (kVec(ik));
      kVecsRel(ik)(0) = int ( lround(relVecDouble(0)) );
      kVecsRel(ik)(1) = int ( lround(relVecDouble(1)) );
      kVecsRel(ik)(2) = int ( lround(relVecDouble(2)) );

      sxprintf ("kRel( %d ) = ", ik);  cout << kVecsRel(ik) << endl;
   }
}

