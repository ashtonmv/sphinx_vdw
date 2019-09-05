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
#include <SxEPRHyper.h>
#include <SxElemDB.h>
#include <SxProjector.h>
#ifndef SX_STANDALONE

//  Ref 1: Phys. Rev. B 47, 4244 (1993)
//  Ref 2: Phys. Rev. B 62, 6158 (2000)

void SxEPRHyper::readHfData (const SxSymbolTable *table, bool gyrRatioTable)
{
   SxElemDB elemDB;

   try {
      const SxSymbolTable *potGroup = SxSpeciesData::findTable(table);
      const SxSymbolTable *species = potGroup->getGroup("species");
      int is, nSpecies = species->getNItems ("species");
      gyromagneticRatio.resize (nSpecies);
      rhoNucleus.resize (nSpecies);
      nucCharge.resize(nSpecies);

      for (is = 0; is < nSpecies; species = species->nextSibling ("species"), is++) {

         SxString elName = species->get("element")->toString ();

         if ( species->contains ("gyromagneticRatio") )  {

            gyromagneticRatio(is) = species->get ("gyromagneticRatio")->toReal ();
            cout << "Using gyromagneticRatio from species { }" << endl;

         } else if ( (gyrRatioTable) && elemDB.getGyrRatio (elName) != 0 )  {

            gyromagneticRatio(is) = elemDB.getGyrRatio (elName);
            cout << "Using gyromagneticRatio from elements.sx" << endl;

         } 

         else gyromagneticRatio(is) = 0.;
            
         rhoNucleus(is) = (species->contains ("rhoNucleus"))
                               ? species->get ("rhoNucleus")->toReal ()
                               : 1.;
         nucCharge(is) = elemDB.getAtomicNumber(elName);
         
      }
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}


// =========================
// Single-projector approach
// =========================

SxEPRHyper::SxRadPsi SxEPRHyper::readFort39 (FILE *fp)
{
   SX_CHECK (fp);
   // --- read until first '#'
   while (fgetc(fp) != '#') { 
      if (feof(fp))  {
         cout << "Unexpected end of file while reading fort.39-file" << endl;
         SX_QUIT;
      }
   }
   // --- read until end of line
   while (fgetc(fp) != '\n' && !feof(fp)) ;
   if (feof(fp))  {
      cout << "Unexpected end of file while reading fort.39-file" << endl;
      SX_QUIT;
   }

   // --- read rad and u
   SxStack<double> rad, u;
   double rPoint, uPoint;
   while (fscanf(fp, "%lf%lf", &rPoint, &uPoint) == 2)  {
      rad << rPoint;
      u << uPoint;
   }

   // --- setup return container
   SxRadPsi res;
   res.rad = rad;
   res.u = u;
   cout << "nPoints = " << res.rad.getSize () << endl;
   return res;
}


void SxEPRHyper::readAePsData (const SxArray<SxString> &psiFiles,
                               bool nonRelIso)
//TODO: clean-up function (parts might be separated singleProj/PAW)
{
   int nSpecies = int(gyromagneticRatio.getSize ());
   SX_CHECK (nSpecies > 0, nSpecies);
   SX_CHECK (nSpecies >= psiFiles.getSize (), nSpecies,
             psiFiles.getSize ());
   // --- print <r-3> matrix element differences
   rm3diff.resize (nSpecies);
   sEnhancement.resize (nSpecies);
   deltaRhoPSAE.resize (nSpecies);

   //TODO: SxArray(is) <SxMatrix>: <psiAE|Delta|psiAE>: l+1 x l+1 - Matrix

   int iFile = 0;
   
   for (int is = 0; is < nSpecies; ++is)  {
      if (fabs (gyromagneticRatio(is)) < 1e-10)  {

         cout << "Skipping species " << (is+1) 
              << " since gyromagneticRatio=0" << endl;
         sEnhancement(is) = 1.;
         rm3diff(is) = 0.;
      
      } else {
      
         // --- read file
         if (iFile >= psiFiles.getSize ())  {
            cout << "Missing ae/ps file for species " << (is+1) << endl;
            SX_QUIT;
         }
         FILE *fp = fopen (psiFiles(iFile).ascii (), "r");
         if (!fp)  {
            cout << "Cannot open " << psiFiles(iFile) << endl;
            SX_QUIT;
      
         }
      
         SxRadPsi ps, ae;
         // pseudo s
         ps = readFort39 (fp);
         // all-electron s
         ae = readFort39 (fp);

         int nPS = int(ps.rad.getSize ()), nAE = int(ae.rad.getSize ());
         SX_CHECK (nPS > 1, nPS);
         SX_CHECK (nAE > 1, nAE);
         double logDrAE = log(ae.rad(nAE-1)/ae.rad(0)) / (nAE-1);
         double logDrPS = log(ps.rad(nPS-1)/ps.rad(0)) / (nPS-1);

         // linear extrapolation for pseudo-waves
         // CPC 119 67 (1999), p74
         double sPs = (  ps.rad(1) * ps.u(0)/ps.rad(0) 
                      - ps.rad(0) * ps.u(1)/ps.rad(1))
                      / (ps.rad(1) - ps.rad(0));

         double s0;

         if ( nonRelIso == false && nucCharge(is) == 0 )  {

            cout << "+-----------------------------------------------+\n"
                 << "| WARNING: nuclear charge (nucCharge) not set.  |\n"
                 << "+-----------------------------------------------+"
                 << endl;
            cout << "Switching to non-relativistic calculation scheme." 
                 << endl;
            nonRelIso = true;

         }

         if (nonRelIso)  {
            
            // exponential extrapolation for ae-waves
            double f0 = ae.u(0)/ae.rad(0), f1 = ae.u(1)/ae.rad(1);
            double sAe = f0 * pow(f0/f1, -ae.rad(0) / (ae.rad(0) - ae.rad(1)));

            cout << "Assuming non-relativistic pseudopotential" << endl;            
            cout << "rho^AE(0)/rho^PS(0) [original] : " 
                 << sqr(ae.u(0)/ps.u(0)) << endl;
            cout << "rho^AE(0)/rho^PS(0) [extrapolated]: " 
                 << sqr(sAe/sPs) << endl;
         
            deltaRhoPSAE(is) = (1./FOUR_PI) *
                        (sqr (ae.u(0)/ae.rad(0)) - sqr(ps.u(0)/ps.rad(0)));

            sEnhancement(is) = sqr(sAe/sPs);
            s0 = sqr (sAe)/ FOUR_PI; // 4pi = (Y_00)^2, s0=rho^AE(0)
         
         } else  {

            // --- relativistic correction: Phys. Rev. B 35 3271 (1987)
            // a: analytical integral:int [Co*r^(2lambda-2)*delta_T,{r,0,\inf}]
            double rTh = nucCharge(is) * (1./137) * (1./137);
            double lambda = sqrt (1. - (nucCharge(is)/137.)
                                     * (nucCharge(is)/137.));
            double Co = ae.u(0) / pow (ae.rad(0),lambda);         
            double a = pow (2., 3-2*lambda) * Co*Co*PI*pow (rTh, 2*lambda-2) 
                       * (lambda-1) * 1./sin(2.*PI*lambda);

            SxDiracVec<Double> intFac(nAE);

            for (int i = 0; i < nAE; i++)
               intFac(i) = (sqr (ae.u(i)/(ae.rad(i))) 
                            - Co*Co*pow (ae.rad(i), 2*lambda-2))
                            * (rTh/2.) * ae.rad(i) 
                            / ((ae.rad(i)+(rTh/2.))*(ae.rad(i)+(rTh/2.)));
               
            double aeFac = a + intFac.integrate(logDrAE);

            cout << "Assuming scalar-relativistic pseudopotential" << endl;
            cout << "rho^AE(0)/rho^PS(0) [original] : " 
                 <<  sqr(ae.u(0)/ps.u(0)) << endl;
            cout << "rho^AE(0)/rho^PS(0) [corrected] : " 
                 << aeFac * (1. / sqr (sPs)) << endl;

            deltaRhoPSAE(is) = (1./FOUR_PI) *
                               (aeFac - sqr(ps.u(0)/ps.rad(0)));

            sEnhancement(is) = aeFac * (1. / sqr (sPs));
            cout << "enhancement: " << sEnhancement(is) << endl;
            s0 = aeFac/ FOUR_PI; // 4pi = (Y_00)^2, s0=rho^AE(0)

         }

         // prefactor for hyperfine interaction. ref 1, eq. (1), (6), (7)
         // mu0/4pi g_e mu_e /a_0^3= 12.531338 T(esla)
         // mu0  ... permeability of vacuum
         // g_e  ... electron gyromagnetic ratio = 2.0023193043617
         // mu_e ... Bohr magneton = hbar e/ 2 m_e 
         // a_0  ... Bohr radius (a_0^{-3} is unit of spin density)
         // nuclear gyromagneticRatio in MHz/T
         double prefacA = 12.531338 * gyromagneticRatio(is);
         cout << "As-free = " << (8. * PI / 3.) * prefacA * s0 
              << " MHz" << endl;

         bool anisotropic = psPot.pseudoPsi(is).getSize () >= 2;

         if (anisotropic) {

            // pseudo p 
            ps = readFort39 (fp);
            // all-electron p
            ae = readFort39 (fp);
            fclose (fp);

            cout << "Using only s- and p-like contributions..." << endl;       

            // --- compute difference in <r-3> matrix elements
            double rm3ae = ((ae.u/ae.rad).sqr ()).integrate (logDrAE); 
            double rm3ps = ((ps.u/ps.rad).sqr ()).integrate (logDrPS); 
            rm3diff(is) = rm3ae - rm3ps;

            cout << "Integral( ae.u * ae.r ): " << 
                 (ae.u.sqr () * ae.rad).integrate (logDrAE) << endl;
            cout << "Integral( ps.u * ps.r ): " <<
                 (ps.u.sqr () * ps.rad).integrate (logDrPS) << endl;

            cout << "Ap-free = " << prefacA * 0.4 * rm3ae << " MHz" << endl;

         }

         iFile++;
         
      }
   }
}

void SxEPRHyper::compute (const SxAtomicStructure &structure, 
                          const SxPtr<SxMuPW> &muPtr, const SxPW &waves, 
                          const SxFermi &fermi)
{
   SX_CHECK (dynamic_cast<const SxGBasis *>(spinDensityG.getBasisPtr ()));
   const SxGBasis &gBasis 
      = *dynamic_cast<const SxGBasis *>(spinDensityG.getBasisPtr ());

   radBasisPtr = SxPtr<SxRadBasis>::create (psPot.rad, psPot.logDr);

   const SxRadBasis &r = *radBasisPtr;
   int nSpecies = structure.getNSpecies ();

   // ref 1, Sec. II C.2 text before eq. (15)
   // the - is i^2 from the expansion of the pw in Bessel function & d-wave
   double YdPrefactor = -sqrt (FOUR_PI / 5.); 
   for (int is = 0; is < nSpecies; ++is)  {
      // skip species without nuclear spin
      if (fabs(gyromagneticRatio(is)) < 1e-10) continue;
      bool anisotropic = psPot.pseudoPsi(is).getSize () >= 2;

      // prefactor for hyperfine interaction. ref 1, eq. (1), (6), (7)
      // mu0/4pi g_e mu_e /a_0^3= 12.531338 T(esla)
      // mu0  ... permeability of vacuum
      // g_e  ... electron gyromagnetic ratio = 2.0023193043617
      // mu_e ... Bohr magneton = hbar e/ 2 m_e 
      // a_0  ... Bohr radius (a_0^{-3} is unit of spin density)
      // nuclear gyromagneticRatio in MHz/T
      double prefacA = 12.531338 * gyromagneticRatio(is);
      cout << "gyromagnetic ratio for species " << (is+1) << ": "
           << gyromagneticRatio(is) << " MHz/T" << endl;

      SxDiracMat<Double> anisoProjG(gBasis.ng, 5);
      SxDiracMat<Double> anisoProjCutG(gBasis.ng, 5);
      double rm3Rad = 0.;

      if (anisotropic)   {
            SxDiracVec<Double> anisoProj, anisoProjCut;
            SxDiracVec<Double> phi = psPot.pseudoPsi(is)(1);
            double logDr = psPot.logDr(is);

            // --- r^-3 matrix elements (with cutoff)
            // denominator of ref. 1 eq. (23)
            rm3Rad = (phi.sqr () 
                            * fCut (r.radFunc(is), rcut)).integrate (logDr);

            // --- full r-3 projector
            anisoProj = YdPrefactor / r.radFunc(is).cub ();
            anisoProj.handle->auxData.is = is;
            anisoProj.handle->auxData.l = 2;
            anisoProj.setBasis (&r);

            // --- cut-off r-3 projector
            anisoProjCut = fCut (r.radFunc(is), rcut) * anisoProj;
            anisoProjCut.handle->auxData.is = is;
            anisoProjCut.handle->auxData.l = 2;
            anisoProjCut.setBasis (&r);

            // --- Fourier transformation to G-space
            for (int m = -2; m <= 2; ++m)  {
               anisoProj.handle->auxData.m = m;
               anisoProjG.colRef (m+2) <<= (gBasis | anisoProj);
            }
            for (int m = -2; m <= 2; ++m)  {
               anisoProjCut.handle->auxData.m = m;
               anisoProjCutG.colRef (m+2) <<= (gBasis | anisoProjCut);
            }
      }

      //TODO: SxTimer (getPhaseFactors might be expensive -> store them)
      for (int ia = 0; ia < structure.getNAtoms (is); ++ia)  {

         double spinAtNucleus = (spinDensityG.conj () 
                                * gBasis.getPhaseFactors (is, ia)).sum ().re
                                / sqrt (structure.cell.volume);

         // SxTimer
         // ref 1, eq. (28)
         double deltaRho = 0.;
         double aFermi = 0.;

         if (muPtr.getPtr()) {

          SX_CHECK(muPtr);
          const SxMuPW &mu = *muPtr;

          SxVector<Double> weights = waves.getGkBasis ().weights;
   
          int nOrb = mu.getNStates ();
          int nk = waves.getNk ();
          int muPos = -1;

          cout  << "+-------------------------------------------------+\n"
                << "| WARNING: Symmetry issues are not fully resolved |\n"
                << "+-------------------------------------------------+" 
				    << endl;

            // --- determine i (could be improved by making a list of iAtom)
            for (int i = 0; i < nOrb; i++) {

               if ( (mu.psiOrb(i).iSpecies == is) && (mu.psiOrb(i).iAtom == ia)
                     && (mu.psiOrb(i).l == 0) ) {
                  
                     muPos = i;
                     break;
                  }
         
            }

            for (int ik = 0; ik < nk; ik++) {

               for (int iSpin = 0; iSpin <= 1; iSpin++) {

                  PsiG muI = mu (muPos,iSpin,ik);
                  
                  int nStates = fermi.getNStates (ik);

                  for (int i = 0; i < nStates; i++) {

                     double coeffAbsSqr = 
                                        dot (waves(i,iSpin,ik),muI).absSqr ();
                     double fOcc=fermi.focc (i,iSpin,ik);

                     if ( fOcc > 1e-6 ) {

                        deltaRho += fOcc 
                                    * weights (ik) 
                                    * coeffAbsSqr
                                    * deltaRhoPSAE (is)
                                    * ((iSpin==0) ? 1 : -1);

                     }
                  
                  } // nStates
               
               } // iSpin
          
            } // ik
         
            aFermi = (spinAtNucleus + deltaRho);

         } // aIsoProj
         else { 

            // rescale density with rho(0) at high cutoff            
            spinAtNucleus *= rhoNucleus(is);
            aFermi = (spinAtNucleus * sEnhancement(is));

         }

         cout << "is = " << (is+1) << "; ia = " << (ia+1) << endl;
         cout << "spinAtNucleus: " << spinAtNucleus << endl;
         cout << "Fermi contact interaction: "
              << (8. * PI / 3.) *  prefacA * aFermi << " MHz" << endl;
         
         if (anisotropic)  {

            SxMatrix3<Double> tensorA;
            tensorA.set (0.);
            
            for (int m = -2; m <= 2; ++m)  {
               SxComplex16 spinYdrm3Ps = (spinDensityG.conj ()
                                          * gBasis.getPhaseFactors (is, ia)
                                          * anisoProjG.colRef(2+m)).sum ();

               // ref 1, eq. (23)
               SxComplex16 kappa = (spinDensityG.conj ()
                                    * gBasis.getPhaseFactors (is, ia)
                                    * anisoProjCutG.colRef(2+m)).sum ()
                                    / (0.4 * rm3Rad);

               // ref 1, eq. (25)
               double spinYdrm3 = spinYdrm3Ps + kappa * 0.4 * rm3diff(is);

               // --- translate Y_d to cartesian coordinates
               switch (m)  {
                  case -2 : tensorA(0,1) = tensorA(1,0) 
                                         = sqrt(3.)*spinYdrm3;
                            break;
                  case -1 : tensorA(1,2) = tensorA(2,1) 
                                         = -sqrt(3.)*spinYdrm3;
                            break;
                  case  0 : tensorA(0,0) -= spinYdrm3;
                            tensorA(1,1) -= spinYdrm3;
                            tensorA(2,2) += 2. * spinYdrm3;
                            break;
                  case  1 : tensorA(0,2) = tensorA(2,0) 
                                         = -sqrt(3.)*spinYdrm3;
                            break;
                  case  2 : tensorA(0,0) += sqrt(3.)*spinYdrm3;
                            tensorA(1,1) -= sqrt(3.)*spinYdrm3;
             
                }
            
               }

            tensorA *= prefacA;
 
            SxMatrix<Double>::Eigensystem eigA 
                               = SxMatrix<Double>(tensorA).eigensystem ();

            // --- sort eigenvectors by abs()
            SxDiracVec<Double> eigAbs(3), eigAbsSorted(3);

            for (int i = 0; i < 3; i++) eigAbs(i) = fabs (eigA.vals(i).re);

            SxArray<long int> sortIdx = eigAbs.getSortIdx();

            for (int i = 2; i >= 0 ; i--) {
               
               eigAbsSorted(i) = eigA.vals(sortIdx(i)).re;

               cout << "main value = " << eigA.vals(sortIdx(i)).re << " MHz, "
                    << "axis = ["<< eigA.vecs(sortIdx(i),0)  << ", "
                    << eigA.vecs(sortIdx(i),1) << ", "
                    << eigA.vecs(sortIdx(i),2) << "]" << endl; 
            }

            cout << "Uniaxiality: " <<  
            (1/3.)*(eigAbsSorted(2)-(0.5*(eigAbsSorted(1)+eigAbsSorted(0)))) 
                 << " MHz" << endl;
            cout << "Rhombicity: " 
                 << 0.5*(fabs(eigAbsSorted(1)-fabs(eigAbsSorted(0))))  
                 << " MHz" << endl;

         }
         
         cout << endl;

      } // iAtom

   } // iSpecies

}   
#else

#include <SxCLI.h>
#include <SxPAWRho.h>
int main(int argc, char **argv)
{
   initSPHInXMath ();

   SxCLI cli(argc, argv);
   cli.authors = "G. Pfanner, C. Freysoldt";

   cli.preUsageMessage 
      = "This add-on computes EPR hyperfine coupling constants.";

   SxString inFile =
      cli.option ("-i|--input","file","SPHInX input file")
      .toString ("input.sx");

   SxString relStrFile =
      cli.option ("--relaxedStr","file","SPHInX structure file with relaxed structure")
      .toString (""); 

   SxString rhoFile
      = cli.option ("-r|--rho","file","spin-polarized density file")
        .toString ("rho.sxb");

   SxString wavesFileName
      = cli.option ("-w|--waves","SPHI/n/X waves file")
        .toString ("waves.sxb");

   SxArray<SxString> psiFiles
      = cli.option ("--psifile", "file", 
                    "fort.39 file (one option per species)").toList ();

   SxFFT::plannerCLI (cli);

   // compute projections
   bool aIsoProj = cli.option ("--proj","Compute aIso from projections")
                   .toBool ();

   // bool
   bool singleProj = 
      cli.option ("--sProj","Use single-projector implementation (PRB 47, 4244 (1993))")
                 .toBool ();

   // look-up gyrRatio from table
   bool gyrRatioTable = 
      cli.option ("--gyrRat","Use gyromagnetic Ratio from elements.sx")
                 .toBool ();

   // relativistic correction
   bool nonRelIso = 
      cli.option ("--nonrel","Do not apply relativistic correction to aIso")
                   .toBool ();

   // factor by which log-grid is interpolated 
   int refine = cli.option ("--refine","int","radial grid refinement")
                .toInt (1,1);

   // minimum for interpolation grid
   double r0New = cli.option ("--r0", "double", "radial grid min")
                  .toDouble (-1.);

   //rCut-radius (rCut of l=0 pseudopotential)
   double rcut = cli.option ("--rcut", "cutoff", "cutoff radius")
                 .toDouble (1.);

   cli.finalize ();

   // --- read input file
   // TODO: load atomic-structure from waves if parameter fails
   //       structure.read (io);
   SxParser parser;
   SxParser::Table table = parser.read (inFile);
   SxAtomicStructure structure;
   if (relStrFile.getSize () == 0)  {

      structure = SxAtomicStructure (&*table);
      cout << "+-------------------------------------------------+\n"
           << "| WARNING: reading atomic structure from input.sx |\n" 
           << "+-------------------------------------------------+" << endl;
   } else {
      SxParser parser2;
      SxParser::Table table2 = parser2.read (relStrFile, "std/structure.std");
      structure = SxAtomicStructure (&*table2);
   }

   // --- read density file
   SxRBasis R;
   SxVector3<Int> mesh;
   SxBinIO io;
   
   try  {
         io.open (rhoFile, SxBinIO::BINARY_READ_ONLY);
         io.read ("dim", &mesh);
         R.set (mesh, SxCell(io));
   } catch (SxException e)  {
   e.print ();
   SX_EXIT;
   }

  // if (singleProj) //TODO: only define rho if singleProj 
     SxRho rho(io, &R); // const SxRho &rho=in.getRef<SxRho> ()?
  // }
  // else {
       SxPAWRho in;

      in.pwRho.rBasisPtr = &R;
      in.readRho (rhoFile);
      in.potPtr = SxPtr<SxPAWPot>::create (&*table);
      
   //}

   io.close ();

   int nSpin = int(rho.rhoR.getSize ());
   if (nSpin != 2)  {
         cout << "Density file contains " << nSpin 
          	  << " spin channel(s) instead of 2." << endl;
         SX_QUIT;
   }


   //setup plane-wave density
   double gCut = SxGBasis::getGCut(SxGBasis::getECut(&*table));
   SxGBasis gBasis(mesh, structure, gCut);

   SxEPRHyper eprHyper;
   eprHyper.readHfData (&*table, gyrRatioTable);

   if (singleProj) {

      SxPtr<SxMuPW> muPtr;

      // TODO: use rcut for each species from table
      //       might be distilled from SxGkBasis
         
      SxPseudoPot psPot(&*table);

      SxPW waves;
      SxFermi fermi;
      SxPtr<SxGkBasis> GkPtr;

      // setup for explicit computation of projections (single-projector approach)
      // TODO: outsource as a routine
      if (aIsoProj)  {

         // --- read waves.sxb
         try {
            SxBinIO wavesFile (wavesFileName, SxBinIO::BINARY_READ_ONLY);
            waves.read (wavesFile);
            GkPtr = SxPtr<SxGkBasis>::create (wavesFile);
            GkPtr->changeTau (structure);
            waves.setGkBasisPtr (GkPtr);
            fermi.read (wavesFile);
            wavesFile.close ();
         }
         catch (SxException e)  {
            e.print ();
            SX_QUIT;
         }

         // basis set and atomic orbitals
         SxConstPtr<SxRadBasis> rPtr = 
         SxConstPtr<SxRadBasis>::create(psPot.rad, psPot.logDr);  
         const SxRadBasis r = *rPtr;   
  
         SxArray< SxArray<SxDiracVec<TReal8> > > psProj = psPot.pseudoPsi;
  
         for (int iSpecies = 0; iSpecies < psProj.getSize (); iSpecies++) {

            for (int j = 0; j < psProj(iSpecies).getSize (); j++) {
      
               psProj(iSpecies)(j) = 
                  fCut (psPot.rad(iSpecies),rcut) * psProj(iSpecies)(j);
               double Norm = 1./ ((psProj(iSpecies)(j).conj () 
                             * psPot.pseudoPsi(iSpecies)(j)
                             * psPot.rad(iSpecies).cub ()
                             ).integrate(psPot.logDr(iSpecies)));
               psProj(iSpecies)(j) *= Norm;
         
            }

         }

         muPtr = SxPtr<SxMuPW>::create (psProj,rPtr,psPot.pseudoFocc,GkPtr);

         /*
         SxAtomicOrbitals psPotWaves(psPot.pseudoPsi,r);
         SxAtomicOrbitals psProjWaves(psProj,r);

         cout << "(psWave|gBasis)(gBasis|psProj):" 
         << ((psPotWaves(0,0,0,0,0)|gBasis)*(gBasis|psProjWaves(0,0,0,0,0)))
           .sum() 
         << endl;

         */


      }

      for (int iSpecies = 0; iSpecies < structure.getNSpecies (); ++iSpecies) {

         // --- interpolate phi on finer mesh
         cout << "Refining mesh for species " << (iSpecies + 1) 
              << " by factor " << refine << "..." << endl;
         double logDrOrig = psPot.logDr(iSpecies);
         double r0 = psPot.rad(iSpecies)(0);

         // --- define finer mesh
         int nrOrig = int(psPot.rad(iSpecies).getSize ());
         double logDr = logDrOrig / refine;
         int nExtra = (r0New < 0.) ? 0 : int(log (r0/r0New) / logDr);
         int nFine = (nrOrig-1)*refine + 1 + nExtra;
         SxDiracVec<Double> rad(nFine);

         for (int i = 0; i < nFine; ++i) 
            rad(i) = r0 * exp ((i-nExtra) * logDr);
      
         // interpolate
         for (int l = 0; l < psPot.pseudoPsi(iSpecies).getSize (); ++l)  {
            // update rad,psi,logDr in psPot
            psPot.pseudoPsi(iSpecies)(l) 
               = interpolateRad (psPot.pseudoPsi(iSpecies)(l),
                                 psPot.rad(iSpecies)(0),
                                 psPot.logDr(iSpecies),
                                 rad);

         }

         psPot.rad(iSpecies) = rad;
         psPot.logDr(iSpecies) = logDr;
      }
   
      // --- perform EPR hyperfine calculations
      eprHyper.spinDensityG = gBasis | (rho(0)-rho(1));
      eprHyper.rcut = rcut;
      eprHyper.psPot = psPot;
      cout << SX_SEPARATOR;
         
      eprHyper.readAePsData (psiFiles, nonRelIso);

      if (aIsoProj)  eprHyper.compute (structure, muPtr, waves, fermi);
      else eprHyper.compute (structure);
   }

   if (!singleProj)   {

      SxPAWRho pawRhoSpin = in.spin().getRef<SxPAWRho> ();
      pawRhoSpin.pwRho.rhoR(0).setBasis (&R);
      pawRhoSpin.potPtr = in.potPtr;
   
      eprHyper.spinDensityG = gBasis | pawRhoSpin.pwRho(0);
      
      int nSpecies = pawRhoSpin.Dij.getNSpecies ();

      for (int is = 0; is < nSpecies; ++is)  {

         double prefacA = 12.531338 * eprHyper.gyromagneticRatio(is); //TODO: unify with other occurrences

         int nAtoms = int(pawRhoSpin.Dij.getNAtoms (is));

         for (int ia = 0; ia < nAtoms; ++ia)  {
            
            SxRadialMesh rhoRadPS = pawRhoSpin.potPtr->computeRhoPS (pawRhoSpin.Dij(0,is,ia),is,0);
            SxRadialMesh rhoRadAE = pawRhoSpin.potPtr->computeRhoAE (pawRhoSpin.Dij(0,is,ia),is,0);
            
            double spinAtNucleusPW = (eprHyper.spinDensityG.conj () 
                                     * gBasis.getPhaseFactors (is, ia)).sum ().re
                                     / sqrt (structure.cell.volume);

            /*
            SxDiracVec<Double> rMesh(nr);
            rMesh(0) = r0;
            for (int i=1; i < nr; i++) rMesh(i) = r0 * exp(dex*i);
            */
            const SxDiracVec<Double> rad  = pawRhoSpin.potPtr->getRadBasis ().radFunc(is);
            double spinCorrection = rhoRadAE(0,0,0) - rhoRadPS(0,0,0)
                                  // linear extrapolation to nucleus
                                  - (rhoRadAE(1,0,0) - rhoRadAE(0,0,0) - rhoRadPS(1,0,0) + rhoRadPS(0,0,0))
                                    /(rad(1) - rad(0))*rad(0);
                          

            //double spinAtNucleus = spinAtNucleusPW + (1/sqrt(FOUR_PI)) * (rhoRadAE(0,0,0) - rhoRadPS(0,0,0));
            double spinAtNucleus = spinAtNucleusPW + spinCorrection/sqrt(FOUR_PI);
            double aFermi = spinAtNucleus; //TODO: scalar-relativistic correction
            
            cout << "planeWave(0): " << spinAtNucleusPW << endl;
            cout << "rhoRadAE(0,0,0): " << rhoRadAE(0,0,0)/sqrt(FOUR_PI) << endl;
            cout << "rhoRadPS(0,0,0): " << rhoRadPS(0,0,0)/sqrt(FOUR_PI) << endl;

            cout << (8. * PI / 3.) *  prefacA * (spinAtNucleusPW - rhoRadPS(0,0,0)/sqrt(FOUR_PI)) << "should be zero" << endl;
            
            cout << "is = " << (is+1) << "; ia = " << (ia+1) << endl;
            cout << "spinAtNucleus: " << spinAtNucleus << endl;
            cout << "Fermi contact interaction: "
                 << (8. * PI / 3.) *  prefacA * aFermi << " MHz" << endl;

         }


      }
   }

   return 0;
}
#endif
