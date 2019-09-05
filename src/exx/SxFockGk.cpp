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

#include <SxFockGk.h>
#include <SxRBasis.h>
#include <SxFockTimer.h>

SxFockGk::SxFockGk ()
   : nG      (-1),
     gkPtr   (NULL),
     nk      (-1),
     nGkMax  (-1)
{
   lattice = None;
}

/*
SxFockGk::SxFockGk (const SxGkBasis &gk_, int nG_)
   : cellPtr (&(gk_.rBasis->cell)),
     nG      (nG_),
     gkPtr   (&gk_),
     nk      (gkPtr->getNk()),
     nGkMax  (gkPtr->nGkMax)
{
   init ();
}
*/

SxFockGk::SxFockGk (const SxGkBasis &gk_, int nG_, const SxCell &cell_)
   : cell   (cell_),
     nG     (nG_),
     gkPtr  (&gk_),
     nk     (gkPtr->getNk()),
     nGkMax (gkPtr->nGkMax)
{
   init ();
}

void SxFockGk::init ()
{
   // --- setup weight_q/|k-q+G1|^2
   setupAuxTables ();

   // --- resize Fock operator
   Fock.reformat (nG, nG);

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
      default : lattice = None;
                sxprintf ("Unknown lattice structure. "
                          "Sorry.\n");
                SX_QUIT;
                break;
   }
   

   // --- in case of wurtzite, create necessary tables
   if (lattice == Wurtzite)  createWurtziteTable ();

   // --- prepare help functions for Gygi-Baldereschi's singularity trick
   setupSingTable ();
}

void SxFockGk::createWurtziteTable ()
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

void SxFockGk::setupAuxTables ()
{
   // --- setup weight_q/|k-q-G|^2
   int                ik, iq, ig, ng, idx;
   const SxGBasis    *Gq;
   SxVector3<TPrecG>  k, q, GplusQ;
   PrecG              kmGq;
   PrecWeights        qWeight;

   kqG.resize (nk * nk * nGkMax);
   kqG = 0.;  // Must not be -1., since the expression appears in a sum later!

   for (ik = 0; ik < nk; ik++)  {
      k = gkPtr->getK(ik);

      for (iq = 0; iq < nk; iq++)  {
         //q       = gkPtr->getK(iq);
         qWeight = gkPtr->weights(iq);
         Gq      = &( (*gkPtr)(iq) );
         ng      = Gq->ng;

         for (ig = 0; ig < ng; ig++)  {
            GplusQ = Gq->getG(ig);
            idx    = (ik * nk + iq) * nGkMax + ig;
            kmGq   = (k - GplusQ).normSqr();
            if ( !(gkPtr->equal (0.,kmGq)) )  kqG(idx) = qWeight / kmGq;
         }
      }
   }
}

void SxFockGk::setupSingTable ()
{
   PrecG              a = 0., c = 0.; // lattice constants
   double             pMax = 0.; // max. radius of a sphere inside of (b1,b2,b3)
   double             lambda = 0.;    // width of a Gaussian
   double             eps = 1.e-4;
   int                ik, iq;
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
   for (ik = 0; ik < nk; ik++)  {
      k = gkPtr->getK(ik);

      sum = 0.;
      for (iq = 0; iq < nk; iq++)  {
         q       = gkPtr->getK(iq);
         qWeight = gkPtr->weights(iq);

         if (ik != iq)  {  // to remove the singularity from the sum ...
            switch (lattice)  {
               case SC       : sum += qWeight * scF (k-q, a);
                               break;

               case FCC      : sum += qWeight * fccF (k-q, a);
                               break;

               case BCC      : sum += qWeight * bccF (k-q, a);
                               break;

               case Wurtzite : sum += qWeight * wurtziteF (k-q, a, c);
                               break;

               case Cos      : sum += qWeight * cosF (k-q, pMax);
                               break;

               case Gauss    : sum += qWeight * gaussF (k-q, lambda);
                               break;

               default       : SX_EXIT;
                               break;
            }  // switch
         }  // if
      }  // :iq
      singTable(ik) = FItgrl - sum;
cout << "singTable(ik) = " << singTable(ik) << endl;
   }  // :ik
}

PrecG SxFockGk::scF (const SxVector3<TPrecG> &p, PrecG a)
{
   return 0.500 * a * a / ( 3.0 - cos(a*p(0)) - cos(a*p(1)) - cos(a*p(2)) );
}

PrecG SxFockGk::fccF (const SxVector3<TPrecG> &p, PrecG a)
{
   return 0.250 * a * a / ( 3.0 - cos(a*p(0)/2.) * cos(a*p(1)/2.) 
                                - cos(a*p(1)/2.) * cos(a*p(2)/2.) 
                                - cos(a*p(2)/2.) * cos(a*p(0)/2.) ) ;
}

PrecG SxFockGk::bccF (const SxVector3<TPrecG> &p, PrecG a)
{
   return 0.125 * a * a / ( 1.0 - cos(a*p(0)/2.)
                                * cos(a*p(1)/2.) * cos(a*p(2)/2.) );
}

PrecG SxFockGk::wurtziteF (const SxVector3<TPrecG> &p, PrecG a, PrecG c)
{
   double g = 3./2. * a * a / (c * c);
   return 0.750 * a * a / ( 3.0 + g - cos(a*p(0))
                            - 2. * cos(a*p(0)/2.) * cos(a*p(1)/2. * SQRT3)
                            - g * cos(c*p(2)) ) ;
}

PrecG SxFockGk::cosF (const SxVector3<TPrecG> &p, double pMax)
{
   Coord q (p);

   cell.map (&q, SxCell::Origin);  // for periodicity: maps point into
                                       // the 1st reduced zone

   const PrecG q2 = q.normSqr();
   const PrecG qN = sqrt (q2);

   if (qN >= pMax)  return 0.;
   else             return 0.5 * (cos(PI * qN / pMax) + 1.) / q2;
}

PrecG SxFockGk::gaussF (const SxVector3<TPrecG> &p, double lambda)
{
   Coord q (p);

   cell.map (&q, SxCell::Origin);  // for periodicity: maps point into
                                       // the 1st reduced zone

   const PrecG q2 = q.normSqr();
   const PrecG w2 = lambda * lambda;

   return exp (-q2/w2) / q2;
}

PrecG SxFockGk::interpolateWurtziteF (const PrecG caRatio)
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
   //     |-----|-----|-----|-----|----->  ···  <-----|
   //     0   x(1)  x(2)  x(3)  x(4)       ···      x(14)
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
   sxprintf ("|                 I / (a/2pi)² = %g\n", res);
   sxprintf ("|                 uses sampling interval %d\n", i);

   return res;
}
      
      
      
            
void SxFockGk::compute (int ik, const SxPW &waves, const SxFermi &fermi)
{
   SX_CLOCK (Timer::FockCompute);
   cout << "SxFockGk::compute ... starting\n";  cout.flush();
   
   SX_CHECK  (waves.getNSpin() == 1, waves.getNSpin());
               // spin pol. not yet implemented
   SX_CHECK (0 <= ik && ik < nk, ik, nk);
   SX_CHECK (waves.getNk() == nk, waves.getNk(), nk);

   int              iG, iGp, iG1, iGSum, iGpSum;
   int              iv, iq, vIdxk, vIdxq, idx, ng, iL;
   PsiG             psiGk, psiGq;
   PrecCoeffG       x, y;
   SxVector<Int>    L(nG);          // temporary lookup table
   const SxGBasis  *GkPtr, *GqPtr;  // Don't mix it up with *gkPtr!
   const int        nv = fermi.getNValenceBands(0,0);

   // TODO: Dump warning, if nv is not the same for all k-points.

   Fock.set (0.);
   GkPtr = &( (*gkPtr)(ik) );

   for (iv = 0; iv < nv; iv++)  {
      for (iq = 0; iq < nk; iq++)  {
         vIdxq = fermi.getValenceBandIdx (iv,0,iq);
         psiGq = waves (vIdxq,0,iq);
         GqPtr = &( (*gkPtr)(iq) );
         ng    = GqPtr->ng;

         for (iG1 = 0; iG1 < ng; iG1++)  {
            idx = (ik * nk + iq) * nGkMax + iG1;

            for (iL = 0; iL < nG; iL++)  {
               L(iL) = gkPtr->getIdxGSum (iL, GkPtr, iG1, GqPtr);
            }  // :iL

            for (iGp = 0; iGp < nG; iGp++)  {
               iGpSum = L(iGp);

               if (iGpSum >= 0)  {
                  x = kqG(idx) * psiGq(iGpSum);

                  for (iG = iGp; iG < nG; iG++)  {
                     iGSum = L(iG);

                     if (iGSum >= 0)  {
                        y = psiGq(iGSum);
                        Fock(iG,iGp) += x.conj() * y;
                     }
                  }  // :iG
               }
            }  // :iGp
         }  // :iG1
      }  // :iq

      // --- Gygi-Baldereschi singularity handling
      vIdxk = fermi.getValenceBandIdx (iv,0,ik);
      psiGk = waves (vIdxk,0,ik);

      for (iGp = 0; iGp < nG; iGp++)  {
         x = singTable(ik) * psiGk(iGp);

         for (iG = iGp; iG < nG; iG++)  {
            y = psiGk(iG);
            Fock(iG,iGp) += x.conj() * y;
            // --- fill lower triangle
            Fock(iGp,iG) = Fock(iG,iGp).conj();
         }  // :iG
      }  // :iGp
   }  // :iv

   // --- give it the correct prefactor
   Fock *= -FOUR_PI / cell.volume;

   cout << "SxFockGk::compute ... done\n" << endl;  cout.flush();
}

void SxFockGk::computeSymFock (int ik, const SxPW &waves, const SxFermi &fermi)
{
   SX_EXIT;
}

PsiG SxFockGk::operator* (const PsiG &psiG)
{
   SX_CLOCK (Timer::FockApply);
   int nGin = (int)psiG.getSize();

   SX_CHECK (nGin >= nG, nGin, nG);
   
   PsiG psiGcopy(nGin);  // is not temporary, unlike to psiG
                         // SxIdx can't handle temporary objects.
   psiGcopy <<= psiG;
   PsiG out(nGin, 0.);

   // --- multiplication
   out( SxIdx(0,nG-1) ) = Fock ^ psiGcopy( SxIdx(0,nG-1) );

   // --- basis pointer
   out.setBasis (psiG.getBasisPtr());

   return out;
}
