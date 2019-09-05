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

#include <SxSymFockGk.h>
#include <SxProjector.h>
#include <SxFockTimer.h>

SxSymFockGk::SxSymFockGk ()
   : nG       (-1),
     gkSetPtr (NULL),
     nk       (-1),
     nGkMax   (-1)
{
   lattice = None;
}

SxSymFockGk::SxSymFockGk (const SxGkBasis      &gkSet,
                                int             nG_,
                          const SxSymbolTable  *table)
   : cell  (SxAtomicStructure(table->topLevel()).cell),
     nG       (nG_),
     gkSetPtr (&gkSet),
     nk       (gkSetPtr->getNk()),
     nGkMax   (gkSetPtr->nGkMax)
{
   init (table->topLevel());
}

SxSymFockGk::SxSymFockGk (const SxGkBasis      &gkSet,
                                int             nG_,
                          const SxCell         &cell_,
                          const SxSymbolTable  *table)
   : cell  (cell_),
     nG       (nG_),
     gkSetPtr (&gkSet),
     nk       (gkSetPtr->getNk()),
     nGkMax   (gkSetPtr->nGkMax)
{
   init (table->topLevel());
}

void SxSymFockGk::init (const SxSymbolTable *top)
{
   // --- resize Fock operator
   Fock = SxDiracSymMat<TPrecCoeffG> (nG);

   // --- recognize structure type
   typedef SxCell::CellType CellType;
   cellType = CellType (cell);
   cellType.print ();
   CellType::Type type = cellType.getType ();

   switch (type)  {
      case CellType::SimpleCubic       : lattice = SC;
                                         break;
      case CellType::FaceCenteredCubic : lattice = FCC;
                                         break;
      case CellType::BodyCenteredCubic : lattice = BCC;
                                         break;
      case CellType::Hexagonal         : lattice = Wurtzite;
                                         break;
      case CellType::Tetragonal        : lattice = PrimitiveTetragonal;
                                         break;
      default : lattice = None;
                sxprintf ("Unknown lattice structure. "
                          "Sorry.\n");
                SX_QUIT;
                break;
   }
   

   // --- in case of wurtzite, create necessary tables
   if (lattice == Wurtzite)  createWurtziteTable ();

   // --- in case of primitive tetragonal structure, create necessary tables
   if (lattice == PrimitiveTetragonal)  {
      // --- check, that a = b
      if ( fabs(cellType.getA() - cellType.getB()) > 1.e-5 )  {
         sxprintf ("Little error: When using the exchange (Fock) operator "
                   "for a primitive\n"
                   "              tetragonal structure, the cell must be set "
                   "up so, that a = b.\n"
                   "It's currently not the case, but you can easily change "
                   "that in the input "
                   "              file. Try again!\n");
      }

      createPrimTetraTable ();
   }

   // --- get stars of k-point set
   const SxAtomicStructure str (&*top);
   SxKPoints kPoints (str.cell, &*top, true); // Use time reversal symmmetry!
   SxStars   stars (/*str.*/cell, kPoints.kVec, kPoints.weights, true);
   stars.compute ();
   rotations = stars.getRotations();
   qVecsSpecial = kPoints.kVec;
   weights = stars.weights;
cout << "rotations = " << rotations << endl;

   // --- prepare help functions for Gygi-Baldereschi's singularity trick
   setupSingTable (stars);
}

void SxSymFockGk::createWurtziteTable ()
{
   wurtziteTable = SxDiracMat<TPrecG> (15, 2, SxList<PrecG> ()
                   <<  0.   <<  0.
                   <<  0.25 <<  2.41198
                   <<  0.5  <<  4.6399
                   <<  1.   <<  8.1361
                   <<  1.5  << 10.502
                   <<  2.   << 12.4069
                   <<  3.   << 14.8682
                   <<  4.   << 16.4539
                   <<  5.   << 17.6677
                   <<  6.   << 18.7455
                   <<  7.   << 19.8621
                   <<  8.   << 20.4293
                   << 10.   << 21.7925
                   << 12.   << 22.548
                   << 24.   << 26.4598);
}

void SxSymFockGk::createPrimTetraTable ()
{
   // --- created by a Monte Carlo integration
   primTetraTable = SxDiracMat<TPrecG> (15, 2, SxList<PrecG> ()
                    <<  0.   <<  0.
                    <<  0.25 <<  2.52589
                    <<  0.5  <<  4.8609
                    <<  1.   <<  8.55935
                    <<  1.5  << 11.3161
                    <<  2.   << 13.3453
                    <<  3.   << 15.8271
                    <<  4.   << 17.7231
                    <<  5.   << 19.2718
                    <<  6.   << 20.3032
                    <<  7.   << 21.3782
                    <<  8.   << 22.4556
                    << 10.   << 23.681
                    << 12.   << 24.7692
                    << 24.   << 29.9053);
}

void SxSymFockGk::setupSingTable (const SxStars &stars)
{
   PrecG              a = 0., c = 0.; // lattice constants
   double             pMax = 0.; // max. radius of a sphere inside of (b1,b2,b3)
   double             lambda = 0.;   // width of a Gaussian
   double             eps = 1.e-4;
   int                ik, iqSpecial;
   SxVector3<TPrecG>  k,q;
   PrecG              sum, FItgrl;
   PrecWeights        qWeight;
   SxArray<PrecG>     ac;       // lattice constants in case of Wurtzite
   double             caRatio;  // ratio c/a

   singTable.resize (nk);

   // --- get lattice constants from cell
   switch (lattice)  {
      case SC       : a = cellType.getA ();
                      FItgrl = a * a * 9.9774204 / FOUR_PI_2;
                      break;

      case FCC      : a = cellType.getA ();
                      FItgrl = a * a * 4.423758 / FOUR_PI_2;
                      break;

      case BCC      : a = cellType.getA ();
                      FItgrl = a * a * 6.7172039 / FOUR_PI_2;
                      break;

      case Wurtzite : a = cellType.getA ();
                      c = cellType.getC ();
                      caRatio = c / a;
                      FItgrl  = interpolateWurtziteF (caRatio)
                                * a * a / FOUR_PI_2;
                      break;

      case PrimitiveTetragonal :
                      a = cellType.getA ();
                      c = cellType.getC ();
                      caRatio = c / a;
                      FItgrl = interpolatePrimTetraF (caRatio)
                               * a * a / FOUR_PI_2;
                      break;

      // --- experimental feature => results not at all reliable!!! FIXME
      case Cos      : pMax   = cell.getHeights().maxval() / 2.;
                      if (eps * pMax * pMax >= 1)
                         sxprintf ("Think about your unit cell!\n");
                      FItgrl = TWO_PI * pMax;
                      break;

      // --- experimental feature => results not at all reliable!!! FIXME
      case Gauss    : pMax   = cell.getHeights().maxval() / 2.;
                      if (eps * pMax * pMax >= 1)
                         sxprintf ("Think about your unit cell!\n");
                      lambda = pMax / sqrt (-log(eps * pMax * pMax));
                      FItgrl = TWO_PI * SQRT_PI * lambda;  // Gauss
                      break;

      default       : sxprintf ("No lattice type chosen. Sorry.\n");
                      SX_EXIT;
                      break;
   }

   // --- setup helper tables
   SxMatrix3<TPrecG>        S;
   SxVector3<TPrecG>        qRot;
   int                      iS, nS, iSymOp;
   int                      nqSpecial = (int)qVecsSpecial.getSize();
   SX_CHECK (cell.symGroupPtr);
   SxArray<SymMat>    symOp     = cell.symGroupPtr->getSymmorphic ();
   int                nSymOp    = (int)symOp.getSize();
cout << "nqSpecial = " << nqSpecial << endl;

   for (ik = 0; ik < nk; ik++)  {
      k  = gkSetPtr->getK(ik);

      sum = 0.;
      for (iqSpecial = 0; iqSpecial < nqSpecial; iqSpecial++)  {
cout << "iqSpecial = " << iqSpecial << endl;
         q       = qVecsSpecial(iqSpecial);
//cout << "q = " << q << endl;
         qWeight = stars.getWeight(iqSpecial);
         nS      = (int)stars.getStar(iqSpecial).getSize();
cout << "nS = " << nS << endl;

         // It seems to me that the following sum over the symmetries is not
         // necessary because of symmetry reasons. But I do not overview it
         // right now, and the following version it is not wrong at least ...
         for (iS = 0; iS < nS; iS++)  {
            iSymOp = rotations(iqSpecial)(iS);
            if (iSymOp < nSymOp)  S =         symOp(iSymOp);
            else  {               S = (-1.) * symOp(iSymOp - nSymOp);
cout << "   (time rev. symmetry used)" << endl;
            }
cout << "S = " << S << endl;
            qRot = S ^ q;
//cout << "qRot = " << qRot << endl;

            // remove singularity from sum ...
            if ( (k-qRot).normSqr() > 1.e-5 )  {
               switch (lattice)  {
                  case SC       : sum += qWeight * scF (k-qRot, a);
                                  break;

                  case FCC      : sum += qWeight * fccF (k-qRot, a);
//cout << "k = " << k << endl;
//cout << "a = " << a << endl;
//cout << "fccF (k-qRot, a) = " << fccF (k-qRot, a) << endl;
                                  break;

                  case BCC      : sum += qWeight * bccF (k-qRot, a);
                                  break;

                  case Wurtzite : sum += qWeight * wurtziteF (k-qRot, a, c);
                                  break;

                  case Cos      : sum += qWeight * cosF (k-qRot, pMax);
                                  break;

                  case Gauss    : sum += qWeight * gaussF (k-qRot, lambda);
                                  break;

                  default       : SX_EXIT;
                                  break;
               }  // switch
            } // if
         }  // :iS
      }  // :iqSpecial
      singTable(ik) = FItgrl - sum;
//cout << "FItgrl = " << FItgrl << endl;
//cout << "sum = " << sum << endl;
cout << "singTable(ik) = " << singTable(ik) << endl;
//SX_EXIT;
   }  // :ik
}

PrecG SxSymFockGk::scF (const SxVector3<TPrecG> &p, PrecG a)
{
   return 0.500 * a * a / ( 3.0 - cos(a*p(0)) - cos(a*p(1)) - cos(a*p(2)) );
}

PrecG SxSymFockGk::fccF (const SxVector3<TPrecG> &p, PrecG a)
{
   return 0.250 * a * a / ( 3.0 - cos(a*p(0)/2.) * cos(a*p(1)/2.) 
                                - cos(a*p(1)/2.) * cos(a*p(2)/2.) 
                                - cos(a*p(2)/2.) * cos(a*p(0)/2.) ) ;
}

PrecG SxSymFockGk::bccF (const SxVector3<TPrecG> &p, PrecG a)
{
   return 0.125 * a * a / ( 1.0 - cos(a*p(0)/2.)
                                * cos(a*p(1)/2.) * cos(a*p(2)/2.) );
}

PrecG SxSymFockGk::wurtziteF (const SxVector3<TPrecG> &p, PrecG a, PrecG c)
{
   double g = 3./2. * a * a / (c * c);
   return 0.750 * a * a / ( 3.0 + g - cos(a*p(0))
                            - 2. * cos(a*p(0)/2.) * cos(a*p(1)/2. * SQRT3)
                            - g * cos(c*p(2)) ) ;
}

PrecG SxSymFockGk::primTetraF (const SxVector3<TPrecG> &p, PrecG a, PrecG c)
{
   double g = (a * a) / (c * c);  // = (a/c)^2
   return 0.5 * a * a / (2.0 + g - cos(a*p(0)) - cos(a*p(1))
                                                   - g * cos(c*p(2)));
}

PrecG SxSymFockGk::cosF (const SxVector3<TPrecG> &p, double pMax)
{
   Coord q (p);

   cell.map (&q, SxCell::Origin);  // for periodicity: maps point into
                                       // the 1st reduced zone

   const PrecG q2 = q.normSqr();
   const PrecG qN = sqrt (q2);

   if (qN >= pMax)  return 0.;
   else             return 0.5 * (cos(PI * qN / pMax) + 1.) / q2;
}

PrecG SxSymFockGk::gaussF (const SxVector3<TPrecG> &p, double lambda)
{
   Coord q (p);

   cell.map (&q, SxCell::Origin);  // for periodicity: maps point into
                                       // the 1st reduced zone

   const PrecG q2 = q.normSqr();
   const PrecG w2 = lambda * lambda;

   return exp (-q2/w2) / q2;
}

PrecG SxSymFockGk::interpolateWurtziteF (const PrecG caRatio)
{
   SxDiracVec<TPrecG> x = wurtziteTable.colRef(0);
   SxDiracVec<TPrecG> y = wurtziteTable.colRef(1);

   //                      _____
   // --- c/a in units of V 2/3
   const double r = caRatio / sqrt (2./3.);

   // --- may happen in non debug cases
   if (r < 0. || r > 24.)  {
      cout << SX_SEPARATOR;
      sxprintf ("| Unfortunately, you are investigating a very strange system.\n"
              "| The ratio c/a of your Wurtzite structure is out of range:\n"
              "|                 _____                  _____\n"
              "|      0 < c/a / V 2/3  < 24, but c/a / V 2/3  = %g\n", r);
      sxprintf ("|\n"
              "| Sorry.\n");
      SX_QUIT;
   }

   // --- find sampling interval of r
   //     |-----|-----|-----|-----|-----> ... <-----|
   //     0   x(1)  x(2)  x(3)  x(4)      ...     x(14)
   int i = 1;
   while ( r > x(i) && i < x.getSize()-1 )  i++;

   // --- linear interpolation
   double s, t, res;
   s   = (x(i) - x(i-1)) / (y(i) - y(i-1));
   t   = y(i) - s * x(i);
   res = s * r + t;

   cout << SX_SEPARATOR;
   sxprintf ("|                        _____\n"
           "| For debugging:  c/a / V 2/3  = %g\n", r);
   sxprintf ("|                 I / (a/2pi)^2 = %g\n", res);
   sxprintf ("|                 uses sampling interval %d\n", i);

   return res;
}

PrecG SxSymFockGk::interpolatePrimTetraF (const PrecG caRatio)
{
   SxDiracVec<TPrecG> x = primTetraTable.colRef(0);
   SxDiracVec<TPrecG> y = primTetraTable.colRef(1);

   //                      _____
   // --- c/a in units of V 2/3  (for analogy with Wurtzite only)
   const double r = caRatio / sqrt (2./3.);

   // --- may happen in non-debug cases
   if (r < 0. || r > 24.)  {
      cout << SX_SEPARATOR;
      sxprintf ("| Unfortunately, you are investigating a very strange system.\n"
              "| The ratio c/a of your structure is out of range:\n"
              "|                 _____                  _____\n"
              "|      0 < c/a / V 2/3  < 24, but c/a / V 2/3  = %g\n", r);
      sxprintf ("|\n"
              "| Sorry.\n");
      SX_QUIT;
   }

   // --- find sampling interval of r
   //     |-----|-----|-----|-----|----->  ...  <-----|
   //     0   x(1)  x(2)  x(3)  x(4)       ...      x(14)
   int i = 1;
   while ( r > x(i) && i < x.getSize()-1 )  i++;

   // --- linear interpolation
   double s, t, res;
   s   = (x(i) - x(i-1)) / (y(i) - y(i-1));
   t   = y(i) - s * x(i);
   res = s * r + t;

   cout << SX_SEPARATOR;
   sxprintf ("|                        _____\n"
           "| For debugging:  c/a / V 2/3  = %g\n", r);
   sxprintf ("|                 I / (a/2pi)^2 = %g\n", res);
   sxprintf ("|                 uses sampling interval %d\n", i);

   return res;
}

void SxSymFockGk::compute (int ik, const SxPW &waves, const SxFermi &fermi)
{
   SX_CLOCK (Timer::FockCompute);
   cout << "SxSymFockGk::compute ... starting\n";  cout.flush();

   SX_CHECK  (waves.getNSpin() == 1, waves.getNSpin());
               // spin pol. not yet implemented
   SX_CHECK (0 <= ik && ik < nk, ik, nk);
   SX_CHECK (waves.getNk() == nk, waves.getNk(), nk);

   int                 iG, iGp, iG2, iGSum, iGpSum;
   int                 iv, vIdxk, vIdxq, ng;
   int                 iqSpecial, nqSpecial = (int)qVecsSpecial.getSize();
cout << "nqSpecial = " << nqSpecial << endl;
   int                 iS, nS, iSymOp,
                       nSymOp = cell.symGroupPtr->getNSymmorphic ();
   bool                timeRevSym;
   PsiG                psiGk, psiGq;
   PrecCoeffG          x, y;
   SymMat              T;
   SxMatrix3<Int>      symOpRel;
   SxVector3<Int>      gRotRel, gpRotRel;
   SxVector3<Double>   kRot;
   const SxGBasis     *gkPtr, *gqPtr;  // Don't mix it up with *gkSetPtr!
   const int           nv = fermi.getNValenceBands(0,0);
   double              qGSk;

   // Note, that SxCell is constructed column-wisely whereas
   // SxMatrix3<T> is constructed row-wisely currently.
   const SxCell        recCell = cell.getReciprocalCell().transpose();

   // TODO: Dump warning, if nv is not the same for all k-points.

   Fock.set (0.);
   gkPtr = &( (*gkSetPtr)(ik) );

   for (iv = 0; iv < nv; iv++)  {
      for (iqSpecial = 0; iqSpecial < nqSpecial; iqSpecial++)  {
cout << "iqSpecial = " << iqSpecial << endl;
cout << "rotations(iqSpecial) = " << rotations(iqSpecial) << endl;
         vIdxq  = fermi.getValenceBandIdx (iv,0,iqSpecial);
cout << "vIdxq = " << vIdxq << endl;
         psiGq  = waves (vIdxq,0,iqSpecial);
//cout << "|psiGq|^2 = " << (psiGq ^ psiGq).chop() << endl;
         gqPtr  = &( (*gkSetPtr)(iqSpecial) );
         ng     = gqPtr->ng;
cout << "ng = " << ng << endl;
         nS     = (int)rotations(iqSpecial).getSize();
cout << "nS = " << nS << endl;

         for (iS = 0; iS < nS; iS++)  {
cout << "iS = " << iS << " (von " << nS << ")" << endl;
            iSymOp = rotations(iqSpecial)(iS);
cout << "iSymOp = " << iSymOp << endl;
            // T := S^-1
            if ((timeRevSym = (iSymOp >= nSymOp)))  {
               T = cell.symGroupPtr->getSymmorphic(iSymOp-nSymOp).transpose();
               T *= -1.;
cout << "   (time rev. symmetry used)" << endl;
            } else {
               T = cell.symGroupPtr->getSymmorphic(iSymOp).transpose();
            }
cout << "T = " << T << endl;
            kRot = T ^ (gkSetPtr->getK(ik));          // T k = S^-1 k

            // get T in relative coordinates
            symOpRel = recCell.carToRel (T);

            for (iG = 0; iG < nG; iG++)  {
               gRotRel = gkPtr->getGRot (symOpRel, iG);

               for (iGp = iG; iGp < nG; iGp++)  {
                  gpRotRel = gkPtr->getGRot (symOpRel, iGp);

                  for (iG2 = 0; iG2 < ng; iG2++)  {
                     iGSum  = gqPtr->getIdxGSum (iG2, gRotRel);
                     iGpSum = gqPtr->getIdxGSum (iG2, gpRotRel);

                     if (iGSum >= 0 && iGpSum >= 0)  {
                        x = psiGq(iGSum);
                        y = psiGq(iGpSum);

                        // get |q_special - S^-1 ^ k + G2|^2
                        qGSk = (gqPtr->getG(iG2) - kRot).normSqr();

                        if ( !(gkSetPtr->equal (0.,qGSk)) )  {
                           if (!timeRevSym)
                              Fock(iG,iGp) += x * y.conj() / qGSk
                                                * weights(iqSpecial);
                           else
                              Fock(iG,iGp) += x.conj() * y / qGSk
                                                * weights(iqSpecial);
//                              * gkSetPtr->weights(iqSpecial);
//cout << "iG = " << iG << endl;
//cout << "iGp = " << iGp << endl;
/*
if (iG == 0 && iGp == 1)  {
cout << "iG2 = " << iG2 << endl;
cout << "Fock(iG,iGp) = " << Fock(iG,iGp) << endl;
}
*/
                        }  // if
                     }  // if
                  }  // :iG2
               }  // :iGp
            }  // :iG
         }  // :iSymOp
      }  // :iqSpecial
cout << "Fock(0,1) = " << Fock(0,1) << endl;
//cout << "Fock = " << Fock.expand().colRef(0) << endl;

      // --- Gygi-Baldereschi singularity handling
      vIdxk = fermi.getValenceBandIdx (iv,0,ik);
cout << "vIdxk = " << vIdxk << endl;
cout << "ik = " << ik << endl;
      psiGk = waves(vIdxk,0,ik);

      for (iG = 0; iG < nG; iG++)  {
//singTable(ik) = 5.65237;
         x = singTable(ik) * psiGk(iG);

         for (iGp = iG; iGp < nG; iGp++)  {
            y = psiGk(iGp);
            Fock(iG,iGp) += x * y.conj();
/*
if (iG == 0 && iGp == 0)  {
cout << "iG = " << iG << " --> x = " << x << endl;
cout << "iGp = " << iGp << " --> y = " << y << endl;
cout << "singsingTable(ik) = " << singTable(ik) << endl;
cout << "psiGk(iG) = " << psiGk(iG) << endl;
cout << "Fock(0,0) = " << Fock(0,0) << endl;
}
*/
         }  // :iGp
      }  // :iG
   }  // :iv

   // --- give it the correct prefactor
//   Fock *= -FOUR_PI / cell.volume;  // *= not yet implented for triangular
                                      // matrices. Therefore the factor is
                                      // handled in ::operator*, see below.
cout << "Fock(0,1) = " << Fock(0,1) * (-FOUR_PI / cell.volume) << endl;
cout << "Fock = " << Fock.expand().colRef(1) * (-FOUR_PI / cell.volume) << endl;

   cout << "SxSymFockGk::compute ... done\n" << endl;  cout.flush();
}

void SxSymFockGk::compute2 (int ik, const SxPW &waves, const SxFermi &fermi)
{
   SX_CLOCK (Timer::FockCompute);
   cout << "SxSymFockGk::compute2 ... starting\n";  cout.flush();
cout << "Muh!" << endl;

   SX_CHECK  (waves.getNSpin() == 1, waves.getNSpin());
               // spin pol. not yet implemented
   SX_CHECK (0 <= ik && ik < nk, ik, nk);
   SX_CHECK (waves.getNk() == nk, waves.getNk(), nk);

   int                 iG, iGp, iG2, iGSum, iGpSum;
   int                 iv, vIdxk, vIdxq, ng;
   int                 iqSpecial, nqSpecial = (int)qVecsSpecial.getSize();
   int                 iS, nS, iSymOp,
                       nSymOp = cell.symGroupPtr->getNSymmorphic ();
   bool                timeRevSym;
   PsiG                psiGk, psiGq;
   PrecCoeffG          x, y;
   SymMat              T;
   SxMatrix3<Int>      symOpRel;
   SxVector3<Int>      gRotRel, gpRotRel;
   SxVector3<Double>   kRot;
   const SxGBasis     *gkPtr, *gqPtr;  // Don't mix it up with *gkSetPtr!
   const int           nv = fermi.getNValenceBands(0,0);
   double              qGSk;

   // Note, that SxCell is constructed column-wisely whereas
   // SxMatrix3<T> is constructed row-wisely currently.
   const SxCell        recCell = cell.getReciprocalCell().transpose();

   // TODO: Dump warning, if nv is not the same for all k-points.

   Fock.set (0.);
   gkPtr = &( (*gkSetPtr)(ik) );

      for (iqSpecial = 0; iqSpecial < nqSpecial; iqSpecial++)  {
cout << "iqSpecial = " << iqSpecial << endl;
         gqPtr  = &( (*gkSetPtr)(iqSpecial) );
         ng     = gqPtr->ng;
         nS     = (int)rotations(iqSpecial).getSize();

         for (iS = 0; iS < nS; iS++)  {
cout << "iS = " << iS << endl;
            iSymOp = rotations(iqSpecial)(iS);
            // T := S^-1
            if ((timeRevSym = (iSymOp >= nSymOp)))  {
               T = cell.symGroupPtr->getSymmorphic(iSymOp - nSymOp).transpose();
               T *= -1.;
            } else {
               T = cell.symGroupPtr->getSymmorphic(iSymOp).transpose();
            }
            kRot = T ^ (gkSetPtr->getK(ik));          // T k = S^-1 k

            // get T in relative coordinates
            symOpRel = recCell.carToRel (T);

            for (iG = 0; iG < nG; iG++)  {
               gRotRel = gkPtr->getGRot (symOpRel, iG);

               for (iGp = iG; iGp < nG; iGp++)  {
                  gpRotRel = gkPtr->getGRot (symOpRel, iGp);

                  for (iG2 = 0; iG2 < ng; iG2++)  {
                     iGSum  = gqPtr->getIdxGSum (iG2, gRotRel);
                     iGpSum = gqPtr->getIdxGSum (iG2, gpRotRel);

                        // get |q_special - S^-1 ^ k + G2|^2
                        qGSk = (gqPtr->getG(iG2) - kRot).normSqr();
                        if ( !(gkSetPtr->equal (0.,qGSk)) )  {

                     if (iGSum >= 0 && iGpSum >= 0)  {
   for (iv = 0; iv < nv; iv++)  {
         vIdxq  = fermi.getValenceBandIdx (iv,0,iqSpecial);
         psiGq  = waves (vIdxq,0,iqSpecial);
                        x = psiGq(iGSum);
                        y = psiGq(iGpSum);

                           if (!timeRevSym)
                              Fock(iG,iGp) += x * y.conj() / qGSk
                                                * weights(iqSpecial);
                           else
                              Fock(iG,iGp) += x.conj() * y / qGSk
                                                * weights(iqSpecial);
                        }  // if
   }  // :iv
                     }  // if
                  }  // :iG2
               }  // :iGp
            }  // :iG
         }  // :iSymOp
      }  // :iqSpecial

      for (iv = 0; iv < nv; iv++)  {
      // --- Gygi-Baldereschi singularity handling
      vIdxk = fermi.getValenceBandIdx (iv,0,ik);
      psiGk = waves(vIdxk,0,ik);

      for (iG = 0; iG < nG; iG++)  {
         x = singTable(ik) * psiGk(iG);

         for (iGp = iG; iGp < nG; iGp++)  {
            y = psiGk(iGp);
            Fock(iG,iGp) += x * y.conj();
         }  // :iGp
      }  // :iG
   }  // :iv

   // --- give it the correct prefactor
//   Fock *= -FOUR_PI / cell.volume;  // *= not yet implented for triangular
                                      // matrices. Therefore the factor is
                                      // handled in ::operator*, see below.

   cout << "SxSymFockGk::compute2 ... done\n" << endl;  cout.flush();
}

PsiG SxSymFockGk::operator* (const PsiG &psiG)
{
   SX_CLOCK (Timer::FockApply);
   int nGin = (int)psiG.getSize();

   SX_CHECK (nGin >= nG, nGin, nG);

   PsiG psiGcopy(nGin);  // is not temporary, unlike to psiG
                         // SxIdx can't handle temporary objects.
   psiGcopy <<= psiG;
   PsiG out(nGin, 0.);

   // --- multiplication
//   out( SxIdx(0,nG-1) ) = Fock ^ psiGcopy( SxIdx(0,nG-1) );
//                             ^^^^^
//   Kann einer mal 'ne Matrixmultiplikation fuer Dreiecksmatrizen 'reinhauen?
//   ;-) ;-)
//
//   work-around:
   out( SxIdx(0,nG-1) ) = Fock.expand() ^ psiGcopy( SxIdx(0,nG-1) );

   // --- prefactor, which has been omitted in ::compute (...)
   out *= -FOUR_PI / cell.volume;

   // --- basis pointer
   out.setBasis (psiG.getBasisPtr());

   return out;
}
