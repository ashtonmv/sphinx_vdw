// ---------------------------------------------------------------------------
//                                          
//      The ab-initio based multiscale library 
//     
//                  S / P H I / n X
//                  http://www.sphinxlib.de
//     
//      Contact:    Sixten Boeck, boeck@mpie.de
//                  Algorithm Design and Modeling Group
//                  Computational Materials Design
//                  Max-Planck-Institute for Iron Research
//                  40237 Duesseldorf, Germany
//     
//      Authors:    see src/AUTHORS
//     
// ---------------------------------------------------------------------------

#include <SxCLI.h>
#include <SxRadBasis.h>
#include <SxRadGBasis.h>
#include <SxProjector.h>
#include <SxGBasis.h>
#include <SxProjector.h>
#include <SxFermi.h>
#include <SxFFT1d.h>
#include <SxDifferential.h>
#include <SxSymMatrix.h>
#include <SxSimpleParser.h>
#include <SxDefectAlignUtil.h>

using namespace SxDefectAlignUtil;

class SxIntegrator
{
   protected:
      double val;
      double last[4];
      ssize_t n;
   public:
      SxIntegrator () : val(0.), n(0L) { }

      inline void operator+= (double x)
      {
         val += x;
         if (n & 1) val += x;
         if (!n) val -= 0.5 * x;
         n++;
         last[n & 3] = x;
      }

      inline void operator-= (double x)
      {
         operator+= (-x);
      }

      double integral (double dx = 1.) const
      { 
         SX_CHECK (n>=3, n);
         double res = (2. * val - last[n & 3]) /3.;
         if ((n & 1) == 0)  {
            res -= ( 15. * last[n & 3]      // 4/3 - 1/3 - 3/8
                    - 1. * last[(n+1) & 3]  // 1/3 - 3/8
                    + 5. * last[(n+2) & 3]  // 4/3 - 9/8
                    -11. * last[(n+3) & 3]) // 2/3 - 9/8
                 / 24.;
         } else {
         }
         return res * dx;
      }

      void reset ()  {
         n = 0; val = 0.;
      }
};


class ChargeDistribution
{
   public:
      class Gaussian {
         public:
            double posZ;
            double betaZ, betaPara;
            double Q;
            Gaussian (const SxSymbolTable *table);
            Gaussian () {
               posZ = betaZ = betaPara = Q = sqrt(-1.);
            }
            double getVal (double z, double kPara) const
            {
               double dz = (z - posZ) / betaZ;
               return Q / (sqrt(TWO_PI) * betaZ)
                      * exp(-0.5 * (sqr(dz) + sqr(kPara * betaPara)));
            } 
      };
      class Slab {
         public:
            double fromZ, toZ;
            double Q;
            Slab (const SxSymbolTable *table);
            Slab () { fromZ = toZ = Q = sqrt(-1.); }
      };
      SxArray<Gaussian> charges;
      SxArray<Slab> slabs;
      SxVector<Double> getQ(double k, ssize_t N, double delta) const;
      SxVector<Double> getQ(double k, ssize_t N, double delta, int shiftZ) const;
      ChargeDistribution (const SxSymbolTable* table);
};

ChargeDistribution::ChargeDistribution (const SxSymbolTable *table)
{
   SxStack<Gaussian> stack;
   SxStack<Slab> slabStack;
   SYMBOLPARSE(table) {
      FOREACH_SYMBOLGROUP("charge")  {
         stack << Gaussian(SYMBOLGROUP_TABLE);
      }
      FOREACH_SYMBOLGROUP("Qslab")  {
         slabStack << Slab(SYMBOLGROUP_TABLE);
      }
   } 
   charges = stack;
   slabs = slabStack;
}

SxVector<Double> ChargeDistribution::getQ(double k, ssize_t N, double delta) const
{
   SxVector<Double> res(N);
   res.set (0.);
   if (charges.getSize () > 0)
      SX_LOOP(iq)  {
         double betaZ = charges(iq).betaZ;
         int nBroad = int(10. * betaZ / delta);
         ssize_t z0 = ssize_t((charges(iq).posZ) / delta);
         for (ssize_t iz = z0 - nBroad; iz <= z0 + nBroad; ++iz)  {
            ssize_t izMod = iz % N; if (izMod < 0) izMod += N;
            res(izMod) += charges(iq).getVal(double(iz) * delta, k);
         }
      }
   if (fabs(k) < 1e-12 && slabs.getSize () > 0)  {
      SX_LOOP(is)  {
         int from = int(slabs(is).fromZ / delta);
         int to = int(slabs(is).toZ / delta);
         if (to == from) to++; // sheet charge layer 
         double val = slabs(is).Q / ((to - from) * delta);
         for (ssize_t iz = from; iz < to; ++iz)  {
            ssize_t izMod = iz % N; if (izMod < 0) izMod += N;
            res(izMod) += val;
         }
      }
   }
   return res;
}

SxVector<Double> ChargeDistribution::getQ(double k, ssize_t N, double delta,
                                          int shiftZ) const
{
   SxVector<Double> res(N);
   res.set (0.);
   for (int iq = 0; iq < charges.getSize (); iq++)  {
      double betaZ = charges(iq).betaZ;
      int nBroad = int(10. * betaZ / delta);
      ssize_t z0 = ssize_t((charges(iq).posZ) / delta);
      for (ssize_t iz = z0 - nBroad; iz <= z0 + nBroad; ++iz)  {
         ssize_t izMod = iz + shiftZ;
         if (izMod >= 0 && izMod < N)
            res(izMod) += charges(iq).getVal(double(iz) * delta, k);
      }
   }
   if (fabs(k) < 1e-12 && slabs.getSize () > 0)  {
      SX_LOOP(is)  {
         int from = int(slabs(is).fromZ / delta);
         int to = int(slabs(is).toZ / delta);
         if (to == from) to++; // sheet charge layer 
         double val = slabs(is).Q / ((to - from) * delta);
         for (ssize_t iz = from; iz < to; ++iz)  {
            ssize_t izMod = iz + shiftZ;
            if (izMod >= 0 && izMod < N)
               res(izMod) += val;
         }
      }
   }
   return res;
}

ChargeDistribution::Gaussian::Gaussian (const SxSymbolTable *table)
{
   SYMBOLPARSE(table) {
      posZ     = SYMBOLGET("posZ");
      betaZ    = SYMBOLGET("betaZ")    || SYMBOLGET("beta") || 1.;
      betaPara = SYMBOLGET("betaPara") || SYMBOLGET("beta") || 1.;
      Q        = SYMBOLGET("Q");
   }
}

ChargeDistribution::Slab::Slab (const SxSymbolTable *table)
{
   SYMBOLPARSE(table) {
      fromZ = SYMBOLGET("fromZ");
      toZ   = SYMBOLGET("toZ");
      Q     = SYMBOLGET("Q");
   }
}


// delta iteration for layered system
SxMatrix<Double> setupLayer(const SxVector<Double> &eps,
                            bool pbc = false)
{
   ssize_t N = eps.getSize ();
   SxMatrix<Double> DQ2(2*N, 2*N);
   DQ2.set (0.);
   SX_LOOP(i)  
   {
      // [eps(z) - eps(z+delta)] / [eps(z) + eps(z+delta)]
      double dezp = (i < N - 1) ? ((eps(i) - eps(i+1))/(eps(i) + eps(i+1)))
                                : 0.;
      // [eps(z) - eps(z-delta)] / [eps(z) + eps(z-delta)]
      double dezn = (i > 0)     ? ((eps(i) - eps(i-1))/(eps(i) + eps(i-1)))
                                : 0.;
      DQ2(0+i,0+i)   = -1.;
      DQ2(0+i,N+i)   = dezp;
      if (i < N - 1)
         DQ2(0+i,1+i)   = 1. + dezp;
      if (pbc && i == N-1)  {
         dezp = (eps(i) - eps(0))/(eps(i) + eps(0));
         DQ2(0+i,N+i)   = dezp;
         DQ2(0+i,0)   = 1. + dezp;
      }

      DQ2(N+i,N+i)   = -1.;
      DQ2(N+i,0+i)   = dezn;
      if (i > 0)
         DQ2(N+i,N+i-1) = 1. + dezn;
      if (pbc && i == 0)  {
         dezn = (eps(i) - eps(N-1))/(eps(i) + eps(N-1));
         DQ2(N+i,0+i)   = dezn;
         DQ2(N+i,N+N-1) = 1. + dezn;
      }
   }

   return DQ2;
}

/* @brief return potential V(z,k) * k via the image charge method
   \f[
     V(z,k) = \frac{2\pi}{k \epsilon(z)} 
              \left(q(z,0) + \sum_{n=1}^{\infty} exp(-k \Delta n)
              (q(z,n) + q(z,-n))\right)
   \f]
   where q(z,n) is the image charge for layer z situated at z+n*delta.
   These image charges are generated iteratively from n->n+1
              
   @param q     charge distribution
   @param dk    \f$\Delta k\f$
   @param pbc   periodic boundary condition?
   */
SxVector<Double> 
imageChargeMethod (const SxVector<Double> &q,
                   const SxVector<Double> &eps,
                   double dk,
                   bool  pbc,
                   bool ddk = false)
{
   SX_CHECK (q.getSize () == eps.getSize (), q.getSize (), eps.getSize ());
   ssize_t N = eps.getSize ();
   double expdk = exp(-dk);
   SxVector<Double> work(6 * N);
   double absQ0 = 0.;
   // 0/1   (q+/q-) for either source or dest
   // 2/3   (q+/q-) for either dest or source
   // 4     deps [eps(i)-eps(i+1)]/[eps(i) + eps(i+1)]
   // 5     resulting potential
   // this setup optimizes memory access, since the iterations work
   // local in the work array

   // --- initialize
   SX_LOOP(i)  {
      work(0+6*i) = work(1+6*i) = work(5+6*i) = q(i);
      absQ0 += 2. * fabs(q(i));
      if (i < N-1) work(4+6*i) = (eps(i) - eps(i+1))/(eps(i) + eps(i+1));
   }
   if (absQ0 < 1.) absQ0 = 1.;
   if (ddk) for (ssize_t i = 0; i < N; ++i) work(5+6*i)=0.;
   if (pbc)
      work(6*N-2) = (eps(N-1) - eps(0))/(eps(N-1) + eps(0));
   else
      work(6*N-2) = 0.;

   int source = 0, dest = 2;
   double expdkn = 1.; // exp(-dk * n)
   ssize_t itmax = 10*1000*1000;
   for (ssize_t n = 1; n < itmax; ++n)  {
      double absQ = 0.;
      expdkn *= expdk;

      for (ssize_t i = 0; i < N; ++i)  {
         double qp, qn; // new image charges at +n/-n

         // charges from interface i/i+1
         if (i == N-1)  {
            if (pbc)  {
               double dezp = work(4+6*i);
               qp =  dezp       * work(6*i +source+1)
                  + (1. + dezp) * work(6*0 +source+0);
            } else {
               qp = 0.;
            }
         } else {
            double dezp = work(4+6*i);
            qp =  dezp       * work(6*i    +source+1)
               + (1. + dezp) * work(6*(i+1)+source+0);
         }
         work(6*i+dest+0) = qp;

         // charges from interface i/i-1
         if (i == 0)  {
            if (pbc)  {
               double dezn = -work(6*(N-1)+4);
               qn =  dezn       * work(6*i     +source+0)
                  + (1. + dezn) * work(6*(N-1) +source+1);
            } else {
               qn = 0.;
            }
         } else {
            double dezn = -work(4+6*(i-1));
            qn =  dezn       * work(6*i    +source+0)
               + (1. + dezn) * work(6*(i-1)+source+1);
         }
         work(6*i + dest + 1) = qn;

         // get potential (except prefactors)
         if (ddk)
            work(6*i + 5) -= double(n) * expdkn * (qp + qn);
         else
            work(6*i + 5) += expdkn * (qp + qn);
         SX_CHECK_NUM(qp);
         SX_CHECK_NUM(qn);
         // 1-norm of image charges
         absQ += fabs(qp) + fabs(qn);
      }

      // --- break-off criteria for n-loop
      // convergence
      if (absQ * expdkn < 1e-12 * absQ0) break;
      // divergence
      if (absQ > 1e10 * absQ0 || n+1 >= itmax) {
         cout << "Image charge method does not converge within"
              << itmax << " steps." << endl;
         cout << "|Q|_1 = " <<  absQ * expdkn << endl;
         SxMatrix<Double> M = setupLayer (eps,pbc);
         SX_LOOP(ii) M(ii,ii) += 1.;
         cout << "Eigenvalues of iteration matrix:" << endl;
         cout << M.eigenvalues () << endl;
         cout << M.eigenvalues ().absSqr () << endl;

         for (ssize_t i = 0; i < N; ++i)
            cout << "% " << work(6*i+dest+0) << ' ' << work(6*i+dest+1) << endl;
         SX_EXIT; 
      }

      // switch source and dest
      source = dest;
      dest = 2 - source;

      /*
      if (pbc)  {
         double S = 0.;
         for (ssize_t i = 0; i < N; ++i)
            S += work(6*i+dest+0) + work(6*i+dest+1);
         cout << "@ " << absQ << ' ' << S << endl;
      }
      */
   }

   SxVector<Double> V(N);
   SX_LOOP(i)  {
      V(i) = TWO_PI * work(6*i+5) / eps(i);
   }
   return V;

}

SxVector<Double> Vk0 (const SxVector<Double> &q,
                      const SxVector<Double> &eps,
                      double delta,
                      SxVector<Double> *VF = NULL)
{
   SX_CHECK (q.getSize () == eps.getSize (), q.getSize (), eps.getSize ());
   ssize_t N = eps.getSize ();
   SxVector<Double> V(N); // the potential to be computed
   if (VF) VF->resize (N);
   V(0) = 0.; // arbitrary
   double VE = 0., avgVE = VE; // potential for 0 charge, unit field
   double Q = 0.; // total charge, proportional to eps(z) dV/dz
   // starting Q is arbitrary, will be corrected below
   double avgV = 0.;
   // --- integrate 
   //       d/dz eps(z) d/dz V(z) = q(z)
   //    layerwise: charges sit in the center of each layer
   //    potential is piecewise linear in [0,1/2] and [1/2,1]
   for (ssize_t i = 1; i < N; ++i)  {
      V(i) = V(i-1) + 0.5 * Q  * delta * (1./eps(i-1) + 1./eps(i));
      VE  +=          0.5 * 1. * delta * (1./eps(i-1) + 1./eps(i));
      avgV += V(i);
      avgVE += VE;
      Q+= q(i) * delta;
   }
   // potential drop across the unit cell
   double dV = V(N-1) + 0.5 * Q  * delta * (1./eps(N-1) + 1./eps(0));
   VE +=                0.5 * 1. * delta * (1./eps(N-1) + 1./eps(0));

   // --- now correct for potential drop and set average potential to 0
   dV /= VE;
   avgV -= dV * avgVE;
   avgV /= double(N);
   V(0) -= avgV;
   VE = 0.;
   if (VF) (*VF)(0) = 0.;
   for (ssize_t i = 1; i < N; ++i)  {
      VE +=           0.5 * 1. * delta * (1./eps(i-1) + 1./eps(i));
      V(i) -= dV * VE + avgV;
      if (VF) (*VF)(i) = VE;
   }

   V *= -FOUR_PI; // prefactor
   return V;
}

SxVector<Double> Vk0cut (const SxVector<Double> &q,
                         const SxVector<Double> &eps,
                         double delta,
                         int cut,
                         SxVector<Double> *VF = NULL)
{
   SX_CHECK (q.getSize () == eps.getSize (), q.getSize (), eps.getSize ());
   ssize_t N = eps.getSize ();
   SxVector<Double> V(N); // the potential to be computed
   if (VF) VF->resize (N);

   V(cut) = 0.; // arbitrary
   double VE = 0.; // potential for 0 charge, unit field
   if (VF) (*VF)(cut) = 0.;
   double Q = 0.; // total charge, proportional to eps(z) dV/dz
   // starting Q is arbitrary, will be corrected below
   // --- integrate 
   //       d/dz eps(z) d/dz V(z) = q(z)
   //    layerwise: charges sit in the center of each layer
   //    potential is piecewise linear in [0,1/2] and [1/2,1]

   // integrate from cut to top
   for (ssize_t i = cut + 1; i < N; ++i)  {
      V(i) = V(i-1) + 0.5 * Q  * delta * (1./eps(i-1) + 1./eps(i));
      VE  +=          0.5 * 1. * delta * (1./eps(i-1) + 1./eps(i));
      Q+= q(i) * delta;
      if (VF) (*VF)(i) = VE;
   }
   // transition N-1 => 0
   if (cut > 0)  {
      V(0) = V(N-1) + 0.5 * Q  * delta * (1./eps(N-1) + 1./eps(0));
      VE  +=          0.5 * 1. * delta * (1./eps(N-1) + 1./eps(0));
      Q+= q(0) * delta;
      if (VF) (*VF)(0) = VE;
   }
   // integrate from 0 to cut
   for (ssize_t i = 1; i < cut; ++i)  {
      V(i) = V(i-1) + 0.5 * Q  * delta * (1./eps(i-1) + 1./eps(i));
      VE  +=          0.5 * 1. * delta * (1./eps(i-1) + 1./eps(i));
      Q+= q(i) * delta;
      if (VF) (*VF)(i) = VE;
   }

   V *= -FOUR_PI; // prefactor
   return V;
}


void addDrop (SxVector<Double> *V0,
              const SxVector<Double> &VF,
              const SxVector<Double> &eps,
              double delta,
              int cutZ,
              double dropV /* *area */)
{
   ssize_t N = VF.getSize ();
   double VE = VF(N-1) - VF(0)
             + (0.5/eps(N-1) + 0.5/eps(0)) * delta;
   (*V0) -= (dropV / VE) * VF;
   for (int i = cutZ; i < N; i++) (*V0)(i) += dropV;
}
             

int main (int argc, char **argv)
{
   // --- command line parsing
   SxCLI cli(argc, argv);

   SxString inputFile = cli.option("-i|--input","input file",
                                   "input file").toString ("system.sx");

   ssize_t N = cli.option ("-n", "integer", "number of points")
           .toInt (-1);
   cli.last ().defaultValue = "default: auto";
   ssize_t Np = cli.option ("-p", "integer", "points per z-mesh point")
           .toInt (1);
   double ecut = cli.option ("--ecut", "cutoff", "cutoff energy")
                    .toDouble ();
   double dk = cli.option ("--dk", "dk", "dk for lateral integration")
                    .toDouble (1e-3);
   bool oldBroadening = cli.option ("--old-broadening",
         "use dielectric constant broadening for epsilon rather than "
         "inverse epsilon. This feature is for backward compatibility. The "
         "new default broadening on the inverse better preserves dielectric "
         "interface positions.").toBool ();

   int shiftZGroup = cli.newGroup ("shift");
   double shiftQ = cli.option ("--shift", "shift", "shift all charges in z-direction")
                    .toDouble ();
   double onlyProfile = cli.option ("-q|--onlyProfile", "only k=0 periodic profile for alignment")
                    .toBool ();

   cli.setGroup(cli.generalGroup);
   double shiftV = cli.option ("-C|--potentialShift","shift",
                               "shift periodic potential (in eV)")
                   .toDouble (0.) / HA2EV;

   int vGroup = cli.newGroup ("read potentials");
   SxString refPotFile = cli.option ("--vref", "potential file",
                                     "reference potential")
                         .toString ();
   SxString defPotFile = cli.option ("--vdef", "potential file",
                                     "defect potential")
                         .toString ();
   bool defIsEV = !cli.option ("--defInHartree", 
                               "defect potential file in Hartree").toBool ();
   bool refIsEV = !cli.option ("--refInHartree", 
                               "reference pot file in Hartree").toBool ();
   double avgWidth = cli.option ("--average", "length", "local average (bohr)")
                     .toDouble (0., 0.);
   FileType fileType = getFileType (cli);
   //bool printOrigPot
   //   = cli.option ("--printPot", "print defect and reference potential")
   //     .toBool ();
   //SxString structDefFile
   //   = cli.option ("--structDef", "SPHInX structure file", 
   //                 "evaluate potentials at atomic coordinates given in the SPHInX file").toString ("");
   //SxString structRefFile
   //   = cli.option ("--structRef", "SPHInX structure file", 
   //                 "evaluate potentials at atomic coordinates given in the SPHInX file").toString ("");
   //double betaAvg = cli.option ("--atomAverage","length",
   //                             "Gaussian broadening for atomic-sphere averages").toDouble (1.5, 0.);
   //bool atomField
   //   = cli.option ("--field", "compute long-range fields at atoms")
   //     .toBool ();

   cli.finalize ();
   
   initSPHInXMath ();

   // --- read cell
   bool printFullVLine = false;
   SxCell defCell, refCell;
   SxMesh3D defMesh, refMesh;
   // vDef and vRef in Hartree
   SxMeshR defPot, refPot;
   SxAtomicStructure structRef, structDef;
   if (cli.groupAvailable (vGroup))  {
      printFullVLine = true;
      defPot = getPot (defCell, defMesh, defPotFile, fileType, &structDef);
      if (defIsEV) defPot /= HA2EV;
      cout << "cell defect = " << defCell << endl;
      refPot = getPot (refCell, refMesh, refPotFile, fileType, &structRef);
      if (refIsEV) refPot /= HA2EV;
      cout << "cell bulk = " << refCell << endl;
#ifdef SXDA_STRUCTURE_PRINTSX
      // --- print structure in sx format (for debugging)
      {
         FILE *fp = sxfopen ("defectStruct.sx", "w");
         structDef.fprint (fp);
         fclose (fp);
         fp = sxfopen ("refStruct.sx", "w");
         structRef.fprint (fp);
         fclose (fp);
      }
#endif
   }

   SxParser parser;
#ifdef SX_NOINSTALLPATH
   parser.setValidation (false);
#endif
   SxParser::Table table = parser.read (inputFile, "std/imagecharge.std");
   SxCell cell;
   if (table->containsGroup ("structure")) {
      cell = SxCell (&*table);
   } else if (defCell.volume > 0.) {
      cell = defCell;
      cout << "Periodic supercell taken from " << defPotFile << endl;
      cout << cell << endl;
   } else {
      cout << "Periodic supercell must be specified in input file '" 
           << inputFile << "'." << endl;
      SX_QUIT;
   }
   
   // --- check cell orientation
   if (   fabs(cell.basis(2) ^ Coord(1,1,0)) > 1e-6
       || fabs(cell.basis(0) ^ Coord(0,0,1)) > 1e-6
       || fabs(cell.basis(1) ^ Coord(0,0,1)) > 1e-6
       || cell(2,2) < 0.)
   {
      cout << "The supercell's a3 axis must be parallel to 0,0,1!" << endl;
      cout << "The supercell's a1 and a2 axes must be orthogonal to 0,0,1!"
           << endl;
      cout << "a3(z) must be positive." << endl;
      cout << "Cell: " << cell << endl;
      SX_QUIT;
   }

   cout << "Cutoff = " << ecut << " Ry" << endl;
   if (Np > 1) 
      cout << "Points per FFT mesh: " << Np << endl;
   SxMesh3D mesh = SxGBasis::getMeshSize (ecut, cell);
   if (N > 0)
      mesh(2) = int(N * Np);
   else
      N = mesh(2) * Np;
   cout << "N = " << N << endl;
   double delta = cell.basis(2).norm ()/double(N);
   double area = cell.volume / cell.basis(2).norm ();

   ChargeDistribution charges(&*table);
   if (cli.groupAvailable(shiftZGroup))  {
      SX_LOOP(iq) charges.charges(iq).posZ += shiftQ;
      cout << "Shifting charge by " << shiftQ << endl;
   }
   if (fabs(shiftV) > 1e-10)
      cout << "Shifting potential by " << shiftV * HA2EV << " eV" << endl;

   SxPtr<ChargeDistribution> background;
   int bgCutZ = -1;
   double bgDropV = 0.;
   if (table->containsGroup ("background"))  {
      const SxSymbolTable* bgGroup = table->getGroup ("background");
      background = SxPtr<ChargeDistribution>::create (bgGroup);
      SYMBOLPARSE (bgGroup) {
         double cut = SYMBOLGET("cut") || 0.;
         bgCutZ = int(lround(cut/delta) % N);
         if (bgCutZ < 0) bgCutZ += (int)N;
         bgDropV = SYMBOLGET("dropV") || 0.;
         bgDropV /= HA2EV;
      }
   }

   SxCell recCell = cell.getReciprocalCell ();

   // --- set up eps
   SxVector<Double> eps(N);
   eps.set (1.);
   int isoBottom = -1, isoTop = -1, cutZ = -1, electrodeZ = 0;
   double zField = 0., dropV = 0.;
   double shiftZ0 = 0.;
   SYMBOLPARSE(&*table)  {
      FOREACH_SYMBOLGROUP("slab")  {
         // --- get dielectric slab parameters
         double epsilon = SYMBOLGET("epsilon");
         double fromZ = SYMBOLGET("fromZ");
         double toZ = SYMBOLGET("toZ");
         double broadening = SYMBOLGET("broadening") || 0.7;
         bool oldBroaden2 = SYMBOLGET("oldBroadening").toBool ();

         // --- make it positive
         while (fromZ < 10. * broadening) fromZ += cell(2,2);
         while (toZ < fromZ) toZ += cell(2,2);
         if (oldBroadening || oldBroaden2)  {
            cout << "Using old broadening for slab epsilon="
                 << epsilon << " at " << fromZ << "..." << toZ << endl;
         }

         // --- set the slab epsilon
         for (int iz = min(int((fromZ -10. * broadening)  / delta), 0);
              iz <= int((toZ + 10. * broadening) / delta);
              ++iz)
         {
            double z = iz * delta;
            int izMod = iz % int(N);
            if (oldBroadening || oldBroaden2)  {
               eps(izMod) += (epsilon - eps(izMod)) * 0.25
                           * (1. + erf ((z - fromZ) / broadening))
                           * (1. + erf((toZ - z)   / broadening));
            } else {
               double epsInvOrig = 1./eps(izMod);
               eps(izMod) = 1./(epsInvOrig
                               + (1./epsilon - epsInvOrig) * 0.25
                                 * (1. + erf ((z - fromZ) / broadening))
                                 * (1. + erf((toZ - z)   / broadening)));
            }
         }
      }
      SYMBOLGROUP("isolated") {
         double fromZ = SYMBOLGET("fromZ");
         //if (fromZ < 0.) {
         //   cout << "fromZ for isolated region must be > 0" << endl;
         //   SX_QUIT;
         //}
         double toZ = SYMBOLGET("toZ");
         isoBottom = int(lround(floor(fromZ/delta)));
         isoTop    = int(lround(ceil(toZ/delta)));
         while (isoBottom < 0)  {
            isoBottom += (int)N;
            isoTop += (int)N;
            shiftZ0 += cell(2,2);
         }

      }
      SYMBOLGROUP("periodic") {
         if (HAVE_SYMBOL("cut"))  {
            double cut = SYMBOLGET("cut");
            cutZ = int(lround(cut/delta) % N);
            if (cutZ < 0) cutZ += (int)N;
            zField = (SYMBOLGET("zField") || 0.) / HA2EV; 
            double electrode = SYMBOLGET("zElectrode") || 0.;
            electrodeZ = int(lround(electrode/delta));
         }
         dropV = SYMBOLGET("dropV") || 0.;
         dropV /= HA2EV;
      }
   }
   
   if (fabs(shiftZ0) > 1e-12)  {
      SX_LOOP(iq) charges.charges(iq).posZ += shiftZ0;
      SX_LOOP(is) {
         charges.slabs(is).fromZ += shiftZ0;
         charges.slabs(is).toZ += shiftZ0;
      }
   }
   
   // --- periodic case
   cout << "--- Periodic" << endl;
   SxMesh3D mesh2D(mesh(0), mesh(1), 1);
   int meshSize = mesh2D.product ();
   double periodic = 0.;
   {
      // --- k = 0
      SxVector<Double> qk = charges.getQ(0, N, delta);
      double Q = qk.sum () * delta;
      cout << "Q=" << Q << endl;
      SxVector<Double> qkN = qk - qk.sum () / double(N);
      SxVector<Double> VF, V0;
      if (cutZ < 0 || fabs(dropV) > 1e-6)  {
         V0 = Vk0(qkN, eps, delta, &VF);
         if (cutZ >= 0) addDrop (&V0, VF, eps, delta, cutZ, dropV * area);
         V0 -= V0.sum () / double(N);
      } else  {
         V0 = Vk0cut(qk, eps, delta, cutZ, &VF);
         V0 += VF * zField;
         double VR = V0(cutZ-1)/area
                   - (zField - FOUR_PI * Q/area)*(cutZ-1 - electrodeZ) * delta;
         V0 -= VR * area;
         qkN <<= qk;
      }
      V0 += shiftV * area;
      /*
      cout << "Left field dependence:   " << dot(VF, qk) * delta << endl;
      cout << "Right field dependence:   " << dot(VF, qk) * delta - VF(N-1) * Q << endl;
      if (cutZ < 0)  {
         if (fabs(eps(0) - eps(N-1)) > 1e-3)  {
            double ifPos = (VF(N-1) - delta * double(N-1) / eps(N-1))
                         / (1./eps(0) - 1./eps(N-1));
            cout << "Interface: " << ifPos << endl;
            double effPos = dot(VF, qk) * delta / Q * eps(0);
            cout << "Effective position (wrt interface): "
                 << ((effPos < ifPos)
                  ? (effPos - ifPos)
                  : (effPos - ifPos) * eps(0)/eps(N-1))
                 << endl;
         } else {
            cout << "Dielectric position shift: " 
                 << (VF(N-1) * eps(0) - delta * double(N-1)) << endl;
            cout << "Effective position: " << dot(VF, qk) * delta / Q * eps(0)
                 << endl;
         }
      }
      */
      periodic += dot(qkN, V0) * delta;
      if (background)  {
         // --- handle background potential
         SxVector<Double> bg = background->getQ(0,N,delta);
         //cout << "bg: " << bg.sum () * delta << endl;
         bg -= bg.sum () / double(N);
         SxVector<Double> Vbg = Vk0(bg, eps, delta, &VF);
         addDrop (&Vbg, VF, eps, delta, bgCutZ, bgDropV * area);
         Vbg -= Vbg.sum () / double(N);
         // contribution to potential
         V0 += Vbg;
         // contribution to energy
         periodic += dot(qk, Vbg) * 2. * delta; // 2: compensate double-counting 0.5 below

         // polarization contribution
         double epsInvAvg = 0.;
         SX_LOOP(i) epsInvAvg += 1./eps(i) ;
         epsInvAvg /= double(N);
         double ePol = -0.5 / FOUR_PI * area / (delta * double(N))
                     * (1./epsInvAvg - 1./eps(bgCutZ)) * bgDropV * bgDropV;
         if (fabs(eps(bgCutZ) - 1.) < 1e-3 && fabs(bgDropV) > 1e-3)  {
            cout << "Slab polarization energy: " << (ePol * HA2EV) << " eV"
                 << endl;
            periodic += ePol * /* compensate factors below: */ 2. * area;
            cout << "Note: Static (field-free case) slab dipole is not considered"
                 << endl;
            cout << "      Effective background field: "
                 << bgDropV / (delta * double(N) * epsInvAvg)
                 << " Hartree/bohr" << endl;
         } 

      }
      if (FILE *fp = fopen("system-profile.dat","w"))  {
         for (ssize_t iz = 0; iz < N; iz += Np) {
            fprintf (fp, "%.6f %.10f %.10f %.10f\n",
                     double(iz) * delta, eps(iz), qk(iz), V0(iz) * HA2EV/area);
         }
         fclose (fp);
      }
      if (printFullVLine)  {
         FILE *fp = sxfopen ("vline-eV.dat", "w");
         SxVector<Double> vRef, vDef, vLR;
         vRef = readLine (refCell, refMesh, refPot, 2, N, defCell,refPotFile);
         vDef = readLine (defCell, defMesh, defPot, 2, N, defCell,defPotFile);
         vLR = V0;
         if (avgWidth > 1e-16)  {
            vRef = average (vRef, avgWidth / delta);
            vDef = average (vDef, avgWidth / delta);
            vLR  = average (vLR, avgWidth / delta);
         }
         for (ssize_t iz = 0; iz < N; iz++) {
            double vModel = vLR(iz) * HA2EV/area;
            double vDFT = (vDef(iz) - vRef(iz)) * HA2EV;
            fprintf (fp, "%.6f\t%.4f\t%.4f\t%.4f\n", double(iz) * delta,
                     vModel, vDFT, vDFT - vModel);
         }
         fclose (fp);
      }
      if (onlyProfile) return 0;
      for (int i = 1; i < meshSize; ++i)  {
         double g2 = recCell.relToCar(mesh2D.getMeshVec(i, SxMesh3D::Origin)).normSqr ();
         SxVector<Double> qkp = charges.getQ(sqrt(g2), N, delta);
         SxVector<Double> Vk = imageChargeMethod(qkp,eps,sqrt(g2)*delta,true,false);
         periodic += dot (qkp, Vk) * sqr(delta) / sqrt(g2);
      }
      periodic /= area;
      periodic *= 0.5;
      periodic += shiftV * Q;
      if (cutZ >= 0 && fabs(dropV) <= 1e-6)  {
         double VL = V0(cutZ)/area 
                   + zField * double(N - cutZ + electrodeZ) * delta;
         double VR = V0(cutZ-1)/area
                   - (zField - FOUR_PI * Q/area)*(cutZ-1 - electrodeZ) * delta;
         double QL = zField / FOUR_PI * area;
         double QR = - (QL + Q);
         cout << "Electrode @ z=" << electrodeZ * delta << endl;
         cout << "QL = " << QL << " VL=" << VL * HA2EV << " eV" << endl;
         cout << "QR = " << QR << " VR=" << VR * HA2EV << " eV" << endl;
         periodic += 0.5 * (QL * VL + QR * VR);
      }
   }

   /*
   double k = 1e-3;
   SxVector<Double> qk = charges.getQ(0, N);
   qk -= qk.sum () / N;
   SX_LOOP(iz)  {
      cout << iz << ' ' << qk(iz) << endl;
   }
   cout << endl;
   SxVector<Double> V0 = Vk0(qk,eps,delta);
   cout << V0.sum () << endl;
   cout << V0.normSqr () << endl;
   //while (k > 1e-6)   {
      SxVector<Double> Vk = imageChargeMethod(qk,eps,k*delta,true,true) 
                          * delta * delta;
      Vk -= Vk.sum () / N;
      cout << Vk.sum () << endl;
      cout << Vk.normSqr () << endl;
      double Ek = dot (qk, imageChargeMethod(qk,eps,k*delta,true,false)),
             Edeka = dot (qk, imageChargeMethod(qk,eps,k*delta,true,true))
                   * delta;
      cout << k << ' ' 
           << Ek
           << ' '
           << Edeka
           << ' ' << Ek - Edeka * k
           << endl;
      k *= 0.5;
      SX_LOOP(i)  {
         cout << "@ " << i << " " << V0(i) << " " << Vk(i) << endl;
      }
      cout << "@ &" << endl;
   //}
   */
   cout.precision (12);
   SxIntegrator itg; 

   // shift Z=0 to isolated case's bottom
   double shiftIso = 0.;
   cout << "--- Isolated" << endl;
   if (isoBottom  >= 0)  {
      while (isoTop <= isoBottom) isoTop += (int)N;
      N = isoTop - isoBottom;
      SxVector<Double> epsNew(N);
      SX_LOOP(iz) epsNew(iz) = eps( (iz + isoBottom) % eps.getSize ());
      shiftIso = -isoBottom * delta;
      eps = epsNew;
      cout << "Isolated from " << isoBottom * delta -shiftZ0 << " to " 
           << isoTop * delta -shiftZ0 << " (" << N << " points)" << endl;

      // --- k = 0
      SxVector<Double> qk = charges.getQ(0, N, delta, -isoBottom);
      double Q = qk.sum () * delta;
      cout << "Q=" << Q << endl;
      SxVector<Double> VF;
      Vk0(qk,eps,delta, &VF);
      //cout << "Left D-field dependence:   " << dot(VF, qk) * delta << endl;
      //cout << "Right D-field dependence:   " << dot(VF, qk) * delta - VF(N-1) * Q << endl;
      if (fabs(eps(0) - eps(N-1)) > 1e-3)  {
         double ifPos = (VF(N-1) - delta * double(N-1) / eps(N-1))
                      / (1./eps(0) - 1./eps(N-1));
         double effPos = dot(VF, qk) * delta / Q * eps(0);
         double effPosIf = effPos - ifPos;
         ifPos -= shiftIso + shiftZ0; // change from within isolated to absolute
         cout << "Interface: " << ifPos << endl;
         cout << "zL = " << (isoBottom * delta - shiftZ0) << " eps = " << eps(0) << endl;
         cout << "E-field dependence left:  zeff = " << (ifPos + effPosIf)
              << endl;
         cout << "zR = " << (isoTop * delta - shiftZ0) << " eps = " << eps(N-1) << endl;
         cout << "E-field dependence right: zeff = " 
              << (ifPos + effPosIf * eps(N-1)/eps(0)) << endl;
      } else {
         cout << "No dielectric interface detected." << endl;
         double dielShift = VF(N-1) * eps(0) - delta * double(N-1);
         cout << "Dielectric position shift: " 
              << dielShift << endl;
         double zeffL = dot(VF, qk) * delta / Q * eps(0) - shiftIso - shiftZ0;
         cout << "zL = " << (isoBottom * delta - shiftZ0) << " eps = " << eps(0) << endl;
         cout << "E-field dependence left:  zeff = " << zeffL << endl;
         cout << "zR = " << (isoTop * delta - shiftZ0) << " eps = " << eps(N-1) << endl;
         cout << "E-field dependence right: zeff = " << (zeffL - dielShift) 
              << endl;
      }
      // SX_LOOP(iz) cout << (delta * iz - shiftIso) << ' ' << eps(iz) << ' ' << qk(iz) << endl;
   } else {
      isoBottom = 0;
   }

   for (double k = 0; k*k < ecut; k+=dk)  {
      SxVector<Double> qk = charges.getQ(k, N, delta, -isoBottom);
      double v = dot (qk, imageChargeMethod(qk,eps,k*delta,false));
      itg += v;
   }
   double iso = itg.integral (dk);
   iso *= delta * delta; // for qk continuous -> discrete; for z-integration
   iso *= 0.5 / TWO_PI;
   cout << "---" << endl;
   cout << "isolated energy = " << (iso * HA2EV) << " eV" << endl;

   cout << "periodic energy = " << (periodic * HA2EV) << " eV" << endl;
   cout << "iso - periodic energy = " << (iso - periodic) * HA2EV << " eV" << endl;

}
