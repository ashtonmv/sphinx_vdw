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

#include <SxCLI.h>
#include <SxFermi.h>
#include <SxSpectrum.h>
#include <SxPartialWaveBasis.h>
#include <SxAtomicStructure.h>
#include <SxPAWSet.h>
#include <SxProjector.h>
#include <SxYlm.h>
#include <SxDifferential.h>

void getIntegrals0123 (const SxDiracVec<Double> &rad,
                       double logdr,
                       double rPAW,
                       const SxDiracVec<Double> &phiAE,
                       const SxDiracVec<Double> &psiGrad,
                       const SxDiracVec<Double> &psiCore,
                       double q,
                       int    lamda,
                       int    m,
                       double &in1,
                       double &in2,
                       double &in3,
                       double &in0)
{
   SxDiracVec<Double> r3 (rad.getSize ()), 
                      r2 (rad.getSize ()),
                      r0 (rad.getSize ());
   for (int ir = 0; ir<rad.getSize(); ++ir)  {
      if (rad(ir) < rPAW)  {
         r2(ir) = rad(ir) * rad(ir) 
                * SxYlm::jsb(lamda, q * rad(ir));
         r3(ir) = rad(ir) * r2(ir);
         r0(ir) = rad(ir) * rad(ir) * rad(ir);
      }
      else r3(ir) = 0., r2(ir) = 0., r0(ir) = 0.;
   }
   in1 = (r3 * psiCore * phiAE).integrate (logdr);
   in3 = (r2 * psiCore * phiAE).integrate (logdr);
   in2 = (r2 * psiGrad * phiAE).integrate (logdr)
         - m * in3;
   if (lamda == 0) 
      in0 = (r0*psiCore*phiAE).integrate (logdr);
}

int main (int argc, char **argv)
{
   initSPHInXMath ();
   
   SxCLI cli(argc, argv);

   cli.preUsageMessage 
      = SxString(" Add-on calculates energy loss near edge structure for PAW\n"
                 "  Gaussian broadening applied to each state").wrap ();
   cli.authors = "Siyuan Zhang s.zhang@mpie.de";
   
   SxString wavesFile
            = cli.option ("-w|--waves","file","input waves file")
             .toString   ("waves.sxb");
   
   SxString inputFile 
      = cli.option("--input","file","read k-vectors from this input file")
        .toString ("input.sx");
   
   SxString outFile
            = cli.option ("-o","file","prefix of the output file;"
                          ".dat is always added to the end of the output file")
             .toString   ("elnes");

   SxString coreFile = cli.option ("-c|--core","file","Abinit core waves file")
                      .toString   ("corewf.abinit");
   
   cli.newGroup ("single atom");
   int atom = cli.option ("-a|--atom", 
                          "number from the input file (starting from 0)")
             .toInt(0);

   int multiAtom = cli.newGroup ("multiple atoms");
   int iSpecies = cli.option ("--iSpecies", "species id", 
                              "which species (starting from 1)").toInt () - 1;
   SxArray<int> atomList = cli.option ("-A", "atom list", 
          "list of atoms  within species (comma-separated, ranges with "
          "<start>-<end>, starting from 1").toIdxList ();
   cli.setGroup (cli.generalGroup);

   int n = cli.option ("-n", "quantum number n of the core wave")
          .toInt(1);

   int l = cli.option ("-l", "quantum number l of the core wave")
          .toInt(0);

   double energy = cli.option ("-e|--energy", "energy [keV] of beam electrons")
                  .toDouble(200.) * 1000. / HA2EV;
   if (energy <= 0)  {
      cout << "Beam energy must be positive, set to the default(200)!" << endl;
      energy = 200000. / HA2EV;
   }

   double eloss = cli.option ("--el", "energy loss [eV] of the edge onset")
                 .toDouble(400.) / HA2EV;
   if (eloss <= 0)  {
      cout << "Energy loss must be positive, set to the default(400)!" << endl;
      eloss = 400. / HA2EV;
   }

   double alpha = cli.option ("--alpha", "convergence angle in mrad")
                 .toDouble(0.);
   if (alpha < 0)  {
      cout << "alpha should NOT be negative, set to the default(0.)!" << endl;
      alpha = 0.;
   }

   double beta = cli.option ("--beta", "collection angle in mrad")
                 .toDouble(10.);
   if (beta < 0)  {
      cout << "beta should NOT be negative, set to the default(10.)!" << endl;
      beta = 10.;
   }

   int nq = cli.option ("--nq", "radial sampling of a round detector")
           .toInt(10);
   if (nq < 1)  {
      cout << "nq must be positive, set to the default(10)!" << endl;
      nq = 10;
   }

   int nth = cli.option ("--nth", "azimuthal sampling of a round detector")
            .toInt(10);
   if (nth < 1)  {
      cout << "nth must be positive, set to the default(10)!" << endl;
      nth = 10;
   }

   int lamdaMax = cli.option ("--lamda", "maximum l for lamda expansion")
             .toInt(3);
   if (lamdaMax < 1)  {
      cout << "lamda must be positive, set to the default(3)!" << endl;
      lamdaMax = 3;
   }

   int orth = cli.option ("--orth", "orthogonalise PAW partial waves: "
                          "default no; 1: orthogonalisation")
             .toInt(0);
   if ((orth < 0) || (orth > 1))  {
      cout << "orth out of range, set to the default!" << endl;
      orth = 0;
   }

   int mono = cli.option ("--mono", "final state selection rule: "
                          "default no selection; 0: exclude monopole term")
             .toInt(-1);
   if ((mono < -1) || (mono > 0))  {
      cout << "mono out of range, set to the default!" << endl;
      mono = -1;
   }

   int pole = cli.option ("--pole", "final state selection rule: "
                          "default(-1) no selection; specify nth pole")
             .toInt(-1);
   if ((pole < -1) || (pole > lamdaMax))  {
      cout << "pole out of range, set to the default!" << endl;
      pole = -1;
   }

   int rayl = cli.option ("--rayl", "rayleigh expansion of exp(iqr): "
                          "default(-1) up to lamdaMax; specify nth order")
             .toInt(-1);
   if ((rayl < -1) || (rayl > lamdaMax))  {
      cout << "rayl out of range, set to the default!" << endl;
      rayl = -1;
   }

   int wop = cli.option ("--wop", "1: write complete outputs; default: 0")
            .toInt(0);

   SxVector3<Double> euler (0.,0.,0.);
   if (cli.option ("--euler", "degree", "between q,r coordinates").exists ())  {
      euler = SxVector3<Double> (cli.last ().toList3());
   }
   cli.last ().required (false);
   cli.last ().defaultValue = "default: 0,0,0";

   int rel = cli.option ("--rel", "0: non-relativistic; default: relativistic")
             .toInt(1);

   double broadening = cli.option ("-b|--broad", "energy [eV]", 
                                   "Gaussian weight factor, eV")
                       .toDouble(0.1) / HA2EV;

   double shift = cli.option ("-s|--shift", "energy [eV]", 
          "shift borders of energy interval from (energyMin, energyMax) "
          "to (energyMin-shift, energyMax+shift), eV")
                  .toDouble(2. * broadening * HA2EV) / HA2EV;
   cli.last ().defaultValue = "default: 2 * broadening";

   SxList<double> range = cli.option("--range","range","Emin:Emax (eV)")
                          .toDoubleList ();
   if (range.getSize () != 2 && range.getSize () != 0 && !cli.error)  {
      cout << "Illegal range." << endl;
      cli.setError ();
   }
   
   int nPerPeak = cli.option ("--fine", "integer", 
                              "energy resolution parameter, approximately "
                              "number of points used to display 1/2 Gauss peak")
                  .toInt(15);
   
   cli.finalize ();

   // --- constants
   const double C = 137.036;
   double gamma = 1. + energy / C / C;
   double v0 = C * sqrt ((1. + 1. / gamma) * (1. - 1. / gamma));
   double k0 = gamma * v0;
   double thE = eloss / gamma / v0 / v0;
   cout << "energy = " << energy << " Hartree" << endl 
        << "energy loss = " << eloss << " Hartree" << endl
        << "gamma = " << gamma << endl
        << "v0 = " << v0 << endl
        << "k0 = " << k0 << endl
        << "thetaE = " << thE << " rad" << endl;
   SxComplex16 fac2 = 0., fac3 = 0.;
   if (rel) fac2 = I * sqrt (FOUR_PI / 3.) * v0 / C / C;
   cout << "fac2 = sqrt(4pi/3)I*v_0/c^2 = " << fac2 << endl;

   double qmax = k0 * (alpha + beta) / 1000.;
   double qmin = k0 * thE / 100.;
   if (qmax <= qmin)  {
      cout << "Collection angle is too small to justify multiple sampling: "
           << "Set nq back to 1!" << endl;
      nq = 1;
   }
   if (nq == 1) nth = 1;

   // --- Euler angles between q and r coordinates
   double th1 = euler(0) * PI / 180., 
          th2 = euler(1) * PI / 180., 
          th3 = euler(2) * PI / 180.;
   cout << "theta1 = " << th1 << endl
        << "theta2 = " << th2 << endl
        << "theta3 = " << th3 << endl;
   SxMatrix3<Double> rot (cos(th1)*cos(th3)-cos(th2)*sin(th1)*sin(th3), 
                         -cos(th1)*sin(th3)-cos(th2)*cos(th3)*sin(th1), 
                          sin(th1)*sin(th2),
                          cos(th3)*sin(th1)+cos(th1)*cos(th2)*sin(th3), 
                          cos(th1)*cos(th2)*cos(th3)-sin(th1)*sin(th3), 
                         -cos(th1)*sin(th2),
                          sin(th2)*sin(th3), 
                          cos(th3)*sin(th2), 
                          cos(th2));
   cout << "rot = " << rot << endl;
   
   // --- wigner3j table
   SxYlm Y;
   int lamdaMax2 = (lamdaMax + 1) * (lamdaMax + 1);
   int lamdaMax4 = lamdaMax2 * lamdaMax2;
   int lamdaMax6 = lamdaMax2 * lamdaMax4;
   SxDiracVec<Complex16> cg(lamdaMax6);
   int i = 0;
   for (int l1 = 0; l1 <= lamdaMax; ++l1)  {
      for (int m1 = -l1; m1 <= l1; ++m1)  {
         for (int l2 = 0; l2 <= lamdaMax; ++l2)  {
            for (int m2 = -l2; m2 <= l2; ++m2)  {
               for (int l3 = 0; l3 <= lamdaMax; ++l3)  {
                  for (int m3 = -l3; m3 <= l3; ++m3)  {
                     //wigner3j(i) = Y.wigner3j(l1, l2, l3, m1, m2, m3);
                     cg(i) = sqrt ((2. * l1 + 1.) * (2. * l2 + 1.) * 
                                   (2. * l3 + 1.) / FOUR_PI)
                           * Y.wigner3j(l1, l2, l3, m1, m2, m3) 
                           * Y.wigner3j(l1, l2, l3, 0, 0, 0);
                     if ((m1 % 2) == 1) cg(i) *= -1;
                     ++i;
                  }
               }
            }
         }
      }
   }
   
   // --- read input file
   SxParser parser;
   SxParser::Table table = parser.read(inputFile);
   SxPtr<SxPAWPot> pawPotPtr =SxPtr<SxPAWPot>::create (&*table);
   SxPAWPot &pawPot = *pawPotPtr; 
  
   // --- read waves file (not waves yet)
   int nk = -1;
   SxFermi fermi;
   SxVector<Double> weights;
   SxAtomicStructure structure;
   SxPtr<SxGkBasis> gkBasisPtr;
   SxPtr<SxPartialWaveBasis> pBasisPtr;
   SxBinIO io;
   try  {
      io.open (wavesFile, SxBinIO::BINARY_READ_ONLY);
      nk = io.getDimension ("nk");
      weights.resize (nk);
      io.read ("kWeights", &weights, nk);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   fermi.read (io);

   structure.read (io);

   // --- check the atom specification
   if (cli.groupAvailable (multiAtom)) {
      if (iSpecies < 0 || iSpecies >= structure.getNSpecies ())  {
         cout << "Invalid species " << (iSpecies + 1) << "!" << endl;
         cout << "Structure has " << structure.getNSpecies () << " species." 
              << endl;
         SX_QUIT;
      }
      // --- check atom number
      SX_LOOP(ja)  {
         if (atomList(ja) >= structure.getNAtoms(iSpecies) || atomList(ja)<0) {
            cout << "Invalid atom id " << atomList(ja) + 1 << "!" << endl;
            cout << "Species " << (iSpecies + 1) << " ("
                 << pawPot.prettyName (iSpecies) << ") has only "
                 << structure.getNAtoms(iSpecies) << " atoms." << endl;
            SX_QUIT;
         }
      }
   } else {
      // --- set up atom list with single atom
      if (atom >= structure.getNAtoms() || atom < 0) {
         cout << "Invalid atom id " << atom << "!" << endl;
         cout << "Structure has " << structure.getNAtoms () << " atoms." 
              << endl;
         SX_QUIT;
      }
      iSpecies = structure.getISpecies (atom);
      atomList.resize (1);
      atomList(0) = atom - structure.atomInfo->offset(iSpecies);
   }

   // --- read core waves
   SxDiracVec<Double> rad = pawPot.rad(iSpecies);
   double logdr = pawPot.logDr (iSpecies);
   double rPAW = pawPot.rc (iSpecies) * sqrt(-log(1e-4));
   cout << "Species " << iSpecies << " has rPAW=" << rPAW 
        << " and logdr=" << logdr << endl;
   pawPot.readCoreAbinit (coreFile, iSpecies);
   int nCore = (int)pawPot.lCore(iSpecies).getSize ();
   cout << "Number of core waves: " << nCore << endl;
   SxDiracVec<Double> psiCore, psiGrad;
   for (int iCore = 0; iCore < nCore; ++iCore)  {
      int lCore = pawPot.lCore(iSpecies)(iCore);
      cout << "wave" << iCore << " has l = " << lCore << endl;
      if (l == lCore)  {
         psiCore = pawPot.psiCoreAE(iSpecies).colRef(iCore+n-l-1);
         cout << "Core wave: wave" << iCore+n-l-1 << endl;
         SxDifferential radialGrad (4);
         psiGrad = radialGrad.apply(psiCore) / logdr;
         if (wop)  {
            SxString file = outFile + ".psiCore" + SxString(iCore);
            writePlot (file, rad, psiCore);
            writePlot (file + "Grad", rad, psiGrad);
         }
         iCore = nCore; // quit the loop
      }
   }

   int npt = pawPot.getNProjType (iSpecies);
   int lMax = pawPot.lMax (iSpecies);

   // --- orthogonalise All-Electron partial waves
   for (int ipt = 0; ipt < npt; ++ipt)  {
      SxDiracVec<Double> psii = pawPot.phiAE(iSpecies).colRef(ipt);
      int lf = pawPot.lPhi(iSpecies)(ipt);
      // --- orthogonality between partial and core waves
      for (int iCore = 0; iCore < nCore; ++iCore)  {
         SxDiracVec<Double> 
         corei = pawPot.psiCoreAE(iSpecies).colRef(iCore);
         int li = pawPot.lCore(iSpecies)(iCore);
         if (lf == li)  {
            SxDiracVec<Double> r3 (rad.getSize ());
            for (int ir = 0; ir<rad.getSize(); ++ir)  {
               if (rad(ir) < rPAW)  {
                  r3(ir) = rad(ir) * rad(ir) * rad(ir);
               }
               else r3(ir) = 0.;
            }
            double integ = (r3*psii*corei).integrate(logdr);
            cout << "integral r3 valence" << ipt << " and core" 
                 << iCore << ": " << integ << endl;
            if (orth)  {
               pawPot.phiAE(iSpecies).colRef(ipt) -= integ * corei;
               double inte = (r3*psii*corei).integrate(logdr);
               cout << "integral r3 orth valence" << ipt 
                    << " and core" << iCore << ": " << inte << endl;
            }
         }
         // --- orthogonality between core waves
         if (ipt == npt -1)  {
            for (int jCore = iCore; jCore < nCore; ++jCore)  {
               SxDiracVec<Double>
               corej =  pawPot.psiCoreAE(iSpecies).colRef(jCore);
               int lj = pawPot.lCore(iSpecies)(jCore);
               if (li == lj)  {
                  double 
                  integ = (rad.cub()*corei*corej).integrate(logdr);
                  cout << "integral r3 core" << iCore << " and core"
                       << jCore << ": " << integ << endl;
               }
            }
         }
      }
   }

   gkBasisPtr = SxPtr<SxGkBasis>::create (io, false);
   gkBasisPtr->changeTau (structure);

   // create partial wave basis
   pBasisPtr = SxPtr<SxPartialWaveBasis>::create (pawPotPtr, structure);
   pBasisPtr->createProjBasis (*gkBasisPtr);

   // construct PAW waves container
   int nSpin = fermi.getNSpin ();

   SxArray2<SxDiracVec<Complex16> > projPsis;
   projPsis.reformat (nSpin, nk);

   // TODO: read waves by k-point (needs ReadOnDemand for SxPAWSet...)
   SxPAWSet waves (gkBasisPtr, pBasisPtr, fermi.getNStates (), nSpin);
   const SxPartialWaveBasis &pBasis = *pBasisPtr;
   
   waves.read (io, SxPWSet::KeepGkBasis);
   io.close ();
   SX_LOOP2 (iSpin, ik)  {
      projPsis(iSpin, ik) = pBasis | waves(iSpin, ik);
   }
  
   // --- get minimum and maximum
   // TODO: eloss +
   double eMin = fermi.eps(0,0,0),
          eMax = fermi.eps(0,0,0);
   if (range.getSize () == 0)  {
      double eMinK,eMaxK;
      for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
         for (int ik = 0; ik < nk; ++ik)  {
            eMinK = fermi.eps(iSpin,ik).minval ();
            if (eMinK < eMin) eMin = eMinK;
            eMaxK = fermi.eps(iSpin,ik).maxval ();
            if (eMaxK > eMax) eMax = eMaxK;
         }
      }
   } else {
      eMin = range(0) / HA2EV;
      eMax = range(1) / HA2EV;
      shift = 0.;
   }

   // --- get the start indices in the <pi|psi> list
   SxArray<int> ipStart(atomList.getSize ());
   {
      int ip0 = 0;
      for (int jSpecies = 0; jSpecies < iSpecies; ++jSpecies)
         ip0 += structure.getNAtoms (jSpecies) * pawPot.getNProj (jSpecies);
      SX_LOOP(ia)  {
         ipStart(ia) = ip0 + atomList(ia) * pawPot.getNProj (iSpecies);
         if (wop) cout << "ipBasis=" << ipStart(ia) << endl;
      }
   }

   // --- precompute Clebsch-Gordan coupling coefficients factors
   int ncg = (lMax + 1) * (lMax + 1) * lamdaMax2 * (2 * l + 1);
   SxDiracVec<Complex16> cg1 (ncg), cg2 (ncg), cg3 (ncg);
   cg1.set (0.), cg2.set (0.), cg3.set (0.);
   SxDiracVec<Double> in0 (npt);

   for (int m = -l; m <= l; ++m)  {
      for (int ipt = 0; ipt < npt; ++ipt)  {
         int lf = pawPot.lPhi(iSpecies)(ipt);
         for (int lamda = 0; lamda<=lamdaMax; ++lamda)  {
            for (int mu = -lamda; mu <= lamda; ++mu)  {
               int lamdamu = Y.combineLm (lamda, mu);
               for (int mf = -lf; mf <= lf; ++mf)  {

                  int lmf = Y.combineLm (lf, mf);
                  int icg = lmf + sqr(lMax+1) * (lamdamu + lamdaMax2 * (l + m));

                  double include = 1.;
                  int ilmf = Y.combineLm (lf, -mf) * lamdaMax4;
                  int lm = Y.combineLm (l, m);
                  if ((l < 1) || (abs(l - 1) < abs(m))) include = 0.;

                  cg1(icg) = cg(ilmf + lamdamu * lamdaMax2 + lm);

   cg2(icg) =   cg(Y.combineLm (l+1, -m) * lamdaMax4 + lm * lamdaMax2 + 2)
              * cg(ilmf + lamdamu * lamdaMax2 + Y.combineLm (l + 1, m))
            + include 
              * cg(Y.combineLm (l-1, -m) * lamdaMax4 + lm * lamdaMax2 + 2)
              * cg(ilmf + lamdamu * lamdaMax2 + Y.combineLm (l - 1, m));

   cg3(icg) =   cg(Y.combineLm (l+1, -m) * lamdaMax4 + (lm + 1) * lamdaMax2 + 1)
              * cg(ilmf + lamdamu * lamdaMax2 + Y.combineLm (l + 1, m))
            + include 
              * cg(Y.combineLm (l-1, -m) * lamdaMax4 + (lm + 1) * lamdaMax2 + 1)
              * cg(ilmf + lamdamu * lamdaMax2 + Y.combineLm (l - 1, m));

                 if (l < abs(m + 1)) cg3(icg) = 0.;
              }
            }
         }
      }
   }

   // --- calculate spectrum
   SxSpectrum elnes;
   elnes.nPerPeak = nPerPeak;

   for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
      // --- set up new spectrum
      elnes.set (eMin - shift, eMax + shift, broadening, nq * nth);
      
      // --- add peaks
      for (int q1 = 0; q1 < nq; ++q1)  {
         double qr = 0.;
         if (nq > 1) qr = exp ((q1 / (nq - 1.)) * log (qmax) 
                        + ((nq - 1. - q1) / (nq - 1.)) * log (qmin));
         if (wop) cout << "q1 = " << q1 << "; qr = " << qr << endl;

         // --- integration
         int nin = npt * (lamdaMax + 1);
         SxDiracVec<Double> in1 (nin), in2 (nin), in3 (nin);
         in1.set (0.), in2.set (0.), in3.set (0.);
         in0.set (0.);
         double q = sqrt(sqr(qr) + sqr(k0*thE));
         for (int ipt = 0; ipt < npt; ++ipt)  {
            for (int lamda = 0; lamda<=lamdaMax; ++lamda)  {
               int iin = ipt * (lamdaMax + 1) + lamda;
               getIntegrals0123(rad, logdr, rPAW,
                                pawPot.phiAE(iSpecies).colRef(ipt), psiGrad,
                                psiCore, q, lamda, -l,
                                in1(iin), in2(iin), in3(iin), in0(ipt));
            }
         }

         for (int q2 = 0; q2 < nth ; ++q2)  {
            double qth = double(q2) / nth * TWO_PI;
            SxVector3<Double> qvec (qr*cos(qth), qr*sin(qth), k0*thE);
            qvec = rot ^ qvec;
            //double q = qvec.norm();
            if (q > 1e-12 && fabs(q - qvec.norm ()) > 1e-12*q)
               sxprintf ("error q = %e\n", q - qvec.norm ());
            double fac = 4.*gamma*gamma / (q*q - eloss*eloss/C/C);
            if (wop && (q2 == 0)) 
               cout << "fac = 4gamma^2/(q^2-El^2/c^2) = " << fac << endl;
            if (wop) cout << "qvec = " << qvec << endl;

            SxVector<Double> normSlm(lamdaMax2), Slm(lamdaMax2);
            SxVector<Complex16> Ylm(lamdaMax2);
            normSlm.set (0.), Slm.set (0.), Ylm.set (0.);
            Y.getYlmArray(lamdaMax, qvec, &normSlm);
            for (int il = 0; il <= lamdaMax; ++il)  {
               for (int im = -il; im <= il; ++im)  {
                  int ilm = Y.combineLm (il, im);
                  Slm(ilm) = Y.getYlmNormFactor(il, im) * normSlm(ilm);
                  if (im > 0)  {
                     Ylm(ilm) = (Slm(ilm) + I*Slm(ilm - 2 * im))/SQRT2;
                     Ylm(ilm - 2 * im) 
                   = (Slm(ilm) - I*Slm(ilm - 2 * im))/SQRT2;
                     if ((im % 2) == 1) Ylm(ilm - 2 * im) *= -1;
                  }
                  if (im == 0) Ylm(ilm) = Slm(ilm);
               }
            }

            bool firstq = (q1 == 0) && (q2 == 0);

            for (int m = -l; m <= l; ++m)  {
               if (rel) fac3 = -fac2 * sqrt(2.*(l - m)*(l + m + 1.));
               if (firstq) 
                  cout << "m = " << m << endl << 
                  "fac3 = -fac2*sqrt(2(l-m)(l+m+1)) = " << fac3 << endl;

               int nStates = fermi.getNStates ();
#ifdef USE_OPENMP
#pragma omp parallel for collapse(2)
#endif
               for (int ik = 0; ik < nk; ++ik)  {
                  for (int in = 0; in < nStates; ++in)  {
                     double knweight = weights(ik) * (2. / double(nSpin) 
                                            - fermi.focc(in, iSpin, ik));

                     bool firstkn = (ik == 0) && (in == 0);
                     if (!firstkn && fabs(knweight) < 1e-14) continue;

                     for (int ia = 0; ia < atomList.getSize (); ++ia)  {
                        int ipBasis = ipStart(ia);
                        // --- calculate projected weight
                        SxComplex16 element = 0.;

                        for (int ipt = 0; ipt < npt; ++ipt)  {
                           int lf = pawPot.lPhi(iSpecies)(ipt);
                           bool smono = (mono != 0) || (lf != l);
                           bool spole = (pole==-1) || (pole==abs(lf-l));
                           SxComplex16 iL = 1.;
                           for (int lamda = 0; lamda<=lamdaMax; ++lamda, iL *= I)  {
                              int iin = ipt * (lamdaMax + 1) + lamda;
                              bool srayl = (rayl==-1) || (rayl==lamda);

                              // --- selection rule
                              if (srayl && spole && smono)  {
                                 for (int mu = -lamda; mu <= lamda; ++mu)  {
                                    int lamdamu = Y.combineLm (lamda, mu);
                                    SxComplex16 fac1 = FOUR_PI * Ylm(lamdamu).conj () * iL;
                                    for (int mf = -lf; mf <= lf; ++mf)  {

   int lmf = Y.combineLm (lf, mf);
   int icg = lmf + (lMax + 1) * (lMax + 1) * (lamdamu + lamdaMax2 * (l + m));


   element += fac * fac1 * projPsis(iSpin,ik)(ipBasis + lf + mf, in) *
   (cg1 (icg) * in1 (iin) + fac2 * cg2 (icg) * in2 (iin) 
   + fac3 * cg3 (icg) * in3 (iin)); 

                                    } //mf
                                 } //mu
                              }
                           } //lamda
                           ipBasis += 2 * lf + 1;
                           if (firstq && (m == -l) && firstkn) 
                              cout << "ipBasis = " << ipBasis << endl;
                        } //ipt
                        double peakWeight =  knweight * element.absSqr();
                        if (qr == 0) 
                           elnes.addPeak (fermi.eps(in,iSpin,ik),
                                          peakWeight, 0);
                        else
                           elnes.addPeak (fermi.eps(in,iSpin,ik), 
                                          qr * qr / nth * peakWeight,
                                          q1 * nth + q2);
                     } // ia
                  } //in
               } //ik
               if (wop)  {
                  SxString file = outFile + ".integ.qr" + SxString(q1);
                  writePlot (file, in1, in2, in3);
               }
            } //m
         } //q2 (theta)
      } //q1 (r)
      if (wop) writePlot (outFile + ".cg", cg1, cg2, cg3);
      writePlot (outFile + ".in0", in0, in0);
      
      // compute spectrum
      elnes.compute ();
      // change from Hartree to eV 
      elnes.spectra  /= (iSpin == 0 ? 1. : -1.) * HA2EV;
      elnes.energies *= HA2EV;

      ssize_t nE = elnes.energies.getSize ();
      cout << "nE = " << nE << endl;
      SxVector<Double> eSpec (nE);
      eSpec.set (0.);

      // ---compute total spectrum
      if (nq > 1)  {
         double logdq = (log(qmax) - log(qmin)) / (nq - 1.);
         cout << "logdq = " << logdq << endl;
         for (int q1 = 0; q1 < nq; ++q1)  {
            for (int q2 = 0; q2 < nth; ++q2)  {
               eSpec += logdq * elnes.spectra.colRef (q1 * nth + q2);
            }
         }
      }

      // ---output
      try {
         SxString name = outFile;
         if (nSpin == 2)
            name += (SxList<SxString> () << ".up" << ".down")(iSpin);
         io.open (name + ".dat", SxBinIO::ASCII_WRITE_ONLY);
         elnes.fprint (io.fp);
         io.close ();

         if (nq > 1)  {
            io.open (name + ".eSpec", SxBinIO::ASCII_WRITE_ONLY);
            for (int iE = 0; iE < nE; ++iE)  {
               sxfprintf(io.fp, "%f", elnes.energies(iE));
               sxfprintf(io.fp, "\t%f", eSpec(iE));
               sxfprintf(io.fp, "\n");
            }
            io.close ();
         }
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
   } //iSpin

}
