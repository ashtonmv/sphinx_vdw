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
#include <SxVector.h>
#include <SxCell.h>
#include <SxFFT1d.h>
#include <SxGBasis.h>
#include <SxPreconditioner.h>
#include <SxRBasis.h>
#include <SxProjector.h>
#include <SxNeighbors.h>
#include <SxTimer.h>
#include <SxDefectAlignUtil.h>

using namespace SxDefectAlignUtil;

double rhoRec (double g2, double beta2, double gamma2, double x)
{
   return         x / sqr(1. + gamma2 * g2)
          + (1. - x) * exp( - 0.25 * beta2 * g2);
}

double rhoRecLimit0 (double beta2, double gamma2, double x)
{
   // rhoRec(G->0) -> 1 + rhoRecLimit0 * G^2
   return       x  * (-2. * gamma2)
          + (1.-x) * (-0.25 * beta2);
}

void computeStruct (const SxAtomicStructure & defStr,
                    const SxAtomicStructure & refStr,
                    double betaAvg,
                    const SxMeshR &vDef, 
                    const SxMeshR &vRef, 
                    const SxMesh3D &defMesh, 
                    const SxMesh3D &refMesh, 
                    const Coord &pos,
                    double eCut,
                    double q, double dielecConstant,
                    const SxMatrix3<Double> epsTensor,
                    double beta2, double gamma2, double expNorm,
                    bool field)
{
   SX_CHECK(vDef.getBasisPtr ());
   SX_CHECK(vRef.getBasisPtr ());
   double betaAvg2 = betaAvg * betaAvg;
   
   // construct GBasis
   SxGBasis defG(defMesh, defStr, eCut);
   SxGBasis refG(refMesh, refStr, eCut);

   // evaluate potential at atom positions (with Gaussian) for defect cell
   PsiG vDefG = defG | vDef;
   PsiG vRefG = refG | vRef;
   PsiG gaussDef = exp((-0.5*betaAvg2) * defG.g2);
   gaussDef /= sqrt(defStr.cell.volume);
   PsiG gaussRef = exp((-0.5*betaAvg2) * refG.g2);
   gaussRef /= sqrt(refStr.cell.volume);
   PsiG vLrDef(defG);
   vLrDef(0) = FOUR_PI * q / dielecConstant 
          * rhoRecLimit0 (beta2, gamma2, expNorm);
   for (int ig = 1; ig < defG.ng; ++ig) {
      Coord g = defG.getG (ig);
      vLrDef(ig) = FOUR_PI / (g ^ epsTensor ^ g)
              * q * rhoRec(defG.g2(ig), beta2, gamma2, expNorm)
              * exp ( -I * (g ^ pos));
   }
   if (epsTensor == SxMatrix3<Double> (1,0,0,0,1,0,0,0,1))  {
      vLrDef(0) *= dielecConstant;
      vLrDef /= dielecConstant;
   }
   vLrDef /= sqrt(defStr.cell.volume);

   ofstream outFile("vAtoms.dat");
   SxGrid grid (refStr, 10);
   for (int is = 0; is < refStr.getNSpecies (); ++is)  {
      for (int defIa = 0; defIa < defStr.getNAtoms(is); ++defIa)  {
         // map to defect cell
         Coord atomPos = defStr(is,defIa);
         int refIdx = refStr.find(atomPos, grid);
         int refIa = -1;
         if (refIdx >= 0 && refStr.getISpecies (refIdx, &refIa) == is)  {
            PsiG defPhase = gaussDef * defG.getPhaseFactors (is, defIa);
            PsiG refPhase = gaussRef * refG.getPhaseFactors (is, refIa);
            double vDefAtom = dot(defPhase, vDefG).re;
            double vLrAtom = dot(defPhase, vLrDef).re;
            double vRefAtom = dot(refPhase, vRefG).re;
               outFile << defStr.cell.getMapped (defStr(is,defIa) - pos, 
                     SxCell::WignerSeitz).norm ()
               << '\t' << vLrAtom * HA2EV
               << '\t' << (vDefAtom - vRefAtom) * HA2EV
               << '\t' << (vDefAtom - vRefAtom - vLrAtom) * HA2EV;
            if (field)  {
               for (int idir = 0; idir < 3; ++idir)
                  outFile << '\t' << 
                     (-dot(defG.getPhaseFactors (is, defIa), vLrDef * defG.gVec.colRef(idir)).im 
                      / sqrt(defStr.cell.volume)); 
            }
            SX_LOOP(iDir) outFile << '\t' << atomPos(iDir);
            outFile << endl;
         }
      }
      outFile << endl;
   }
}

int main (int argc, char** argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   SxCLI cli(argc, argv);
   cli.authors = "C. Freysoldt";
   cli.version ("2.2");

   double eCut, q, vAlign, dielecConstant, avgWidth;
   SxMatrix3<Double> epsTensor;
   Coord pos(0.,0.,0.);
   SxString refPotFile, defPotFile;

   int sxInput = cli.newGroup ("sxdefectalign input file");
   SxString inFile = cli.option ("-i|--input", "file",
                                 "sx-style sxdefectalign input file")
                     .toString ();
   
   int clInput = cli.newGroup ("cell from command line");
   cli.excludeGroup(sxInput);
   Coord col1 (cli.option ("--axis1","vector", "a1 (units: bohr)")
               .toList3 ()),
         col2 (cli.option ("--axis2","vector", "a2 (units: bohr)")
               .toList3 ()),
         col3 (cli.option ("--axis3","vector", "a3 (units: bohr)")
               .toList3 ());
   
   cli.newGroup ("parameters from command line");
   cli.excludeGroup(sxInput);
   eCut = cli.option ("--ecut", "energy (Ry)", "cutoff energy in Ry")
                 .toDouble ();

   double step = cli.option ("--gstep", "step size", "|G| step size")
                 .toDouble (1e-4, 1e-10);

   double beta = cli.option ("--beta", "length", 
                             "Gaussian decay exp(-r^2/beta^2)")
                 .toDouble (1.);
   
   int expGroup = cli.newGroup ("exponential model");
   cli.excludeGroup(sxInput);
   double gamma = cli.option ("--gamma", "length", 
                             "exponential decay exp(-r/gamma)")
                  .toDouble ();
   double expNorm = cli.option ("--expnorm", "<0..1>", 
                             "relative norm of exponential part")
                    .toDouble ();
   bool printRho
      = cli.option ("--printRho", "print model rho to rhoModel.dat")
        .toBool ();
   
   cli.newGroup ("CSRB screening");
   cli.excludeGroup(sxInput);
   double q2tf = cli.option ("--csrb", "(screening vector)^2", "square of "
         "Thomas-Fermi screening vector for CSRB screening").toDouble (0.); 
   
   cli.newGroup("model charge specifications");
   cli.excludeGroup(sxInput);
   q = cli.option ("-q|--charge", "excess electrons", "number of "
         "additional excess electrons of the defect: a defect charge of -1 " 
         "corresponds to 1 excess electrons, while a defect charge of +1 " 
         "corresponds to -1 excess electrons")
              .toDouble (1.);
   cli.option ("--pos|--center", "vector", 
                   "defect center position");
   cli.last ().optional = true;
   if (cli.last ().exists ())  {
      pos = Coord(cli.last ().toList3 ());
   }

   bool relative = cli.option ("--relative", "pos is in relative coordinates")
                   .toBool ();

   bool tensor = cli.option ("--tensor", "tensor",
   "dielectric tensor. Specify as\n"
   "eps_xx,eps_yy,eps_zz (3 values) for a diagonal tensor\n"
   "eps_xx,eps_yy,eps_zz,eps_yz,eps_xz,eps_xy (6 values) or ...\n"
   "xx,xy,xz,xy,yy,yz,xz,yz,zz                (9 values) for a full tensor")
                 .exists ();
   cli.last ().optional = true;
   if (tensor)  {
      SxList<double> epsVals = cli.last ().toDoubleList ();
      epsTensor.set (0.);
      switch (epsVals.getSize ())  {
         case 3:
            for (int i = 0; i < 3; ++i) epsTensor(i,i) = epsVals(i);
            break;
         case 6:
            for (int i = 0; i < 3; ++i) epsTensor(i,i) = epsVals(i);
            epsTensor(1,2) = epsTensor(2,1) = epsVals(3);
            epsTensor(0,2) = epsTensor(2,0) = epsVals(4);
            epsTensor(1,0) = epsTensor(0,1) = epsVals(5);
            break;
         case 9:
            epsTensor = SxMatrix3<Double>(epsVals(0), epsVals(1), epsVals(2),
                                          epsVals(3), epsVals(4), epsVals(5),
                                          epsVals(6), epsVals(7), epsVals(8));
            break;
         default:
            cout << "Failed to interpret --tensor with " << epsVals.getSize ()
                 << "values - must be 3, 6, or 9" << endl;
            cli.setError ();
      }
   }

   cli.option ("-e|--eps", "constant", "dielectric constant");
   if (tensor && cli.last ().exists ())  {
      cout << "Both dielectric constant and tensor are given!" << endl;
      cout << "Dielectric constant will be ignored" << endl;
   }
   dielecConstant = cli.last ().toDouble (1.);


   int vGroup = cli.newGroup ("read potentials");
   cli.excludeGroup(sxInput);
   refPotFile = cli.option ("--vref", "potential file",
                                 "reference potential")
                         .toString ();
   defPotFile = cli.option ("--vdef", "potential file",
                                 "defect potential")
                         .toString ();
   bool defIsEV = !cli.option ("--defInHartree", 
                               "defect potential file in Hartree").toBool ();
   bool refIsEV = !cli.option ("--refInHartree", 
                               "reference pot file in Hartree").toBool ();
   FileType fileType = getFileType (cli);
   bool printOrigPot
      = cli.option ("--printPot", "print defect and reference potential")
        .toBool ();
   SxString structDefFile
      = cli.option ("--structDef", "SPHInX structure file", 
                    "evaluate potentials at atomic coordinates given in the SPHInX file").toString ("");
   SxString structRefFile
      = cli.option ("--structRef", "SPHInX structure file", 
                    "evaluate potentials at atomic coordinates given in the SPHInX file").toString ("");
   double betaAvg = cli.option ("--atomAverage","length",
                                "Gaussian broadening for atomic-sphere averages").toDouble (1.5, 0.);
   bool atomField
      = cli.option ("--field", "compute long-range fields at atoms")
        .toBool ();

   cli.setGroup (cli.generalGroup);
   avgWidth = cli.option ("--average", "length", "local average (bohr)")
                     .toDouble (0., 0.);
   vAlign = cli.option ("-C|--align", "align", 
                               "potential alignment constant (eV)")
                   .toDouble (0.) / HA2EV;
   cli.finalize ();

   if (tensor)  {
      if ((epsTensor - epsTensor.transpose ()).absSqr ().sum () > 1e-12)  {
         cout << "Dielectric tensor must be symmetric" << endl;
         cout << "Tensor is " << epsTensor << endl;
         SX_QUIT;
      }
      if (q2tf > 1e-12)  {
         cout << "CSRB screening cannot be used with dielectric tensors"
              << endl;
         SX_QUIT;
      }
      SxSphereGrid grid(SxSphereGrid::Grid_110);
      double epsInv = 0.;
      SX_LOOP(i) {
         Coord xyz = grid.getXyz((int)i);
         epsInv += grid.weights(i) / (xyz ^ epsTensor ^ xyz);
      }
     dielecConstant = 1./epsInv;
   } else {
      epsTensor = SxMatrix3<Double> (1, 0, 0, 0, 1, 0, 0, 0, 1);
   }

   if (!cli.groupAvailable (expGroup))  {
     expNorm = 0.;
     gamma = 1.;
   }
   double beta2 = beta*beta, gamma2 = gamma * gamma;


   // --- read cell
   bool printFullVLine = false;
   SxCell defCell, refCell;
   SxMesh3D defMesh, refMesh;
   // vDef and vRef in Hartree
   SxMeshR defPot, refPot;
   SxAtomicStructure structRef, structDef;

   if (cli.groupAvailable (sxInput)) {
      SxParser parser;
      SxConstPtr<SxSymbolTable> table = parser.read (inFile);
      SxSymbolTable *mainGroup = table->getGroup("DefectAlign");
      
      // --- cell
      SxSymbolTable *subGroup = mainGroup->getGroup("cellGroup");
      if (subGroup->contains("fromVElStat")) 
         cout << "Cells are taken from electrostatic potentials." << endl;
      else if (subGroup->contains("file"))  {
         SxString cellFile = subGroup->get("file")->toString();
         SxParser cellParser;
         SxParser::Table cellTable = cellParser.read (cellFile,"std/structure.std");
         try {
            defCell = SxCell(&*cellTable);
            cout << "cell defect = " << defCell << endl;
         } catch (SxException e)  {
            e.print ();
            SX_EXIT;
         }
      } else if (subGroup->contains("cell"))  {
         defCell = CellMat(subGroup->get("cell")->toList ()).transpose ();
         cout << "cell defect = " << defCell << endl;
      } else {
         cout << "No valid option in cell group!" << endl;
         SX_QUIT;
      }
      
      // --- electrostatic potentials
      SxSymbolTable *vTable = mainGroup->getGroup("vElStat");
      if (vTable->contains("fileType"))  {
         SxString fileMode = vTable->get("fileType")->toString();
         if (fileMode == "VASP") fileType = VASP_LOCPOT;
         else if (fileMode == "SPHINX") fileType = sxb;
         else {
            cout << "Unknown filetype for electrostatic potential." << endl;
            SX_QUIT;
         }
      } else fileType = sxb;
      subGroup = vTable->getGroup("vDefect");
      defPotFile = subGroup->get("file")->toString();
      defPot = getPot (defCell, defMesh, defPotFile, fileType, &structDef);
      defPot /= HA2EV;
      cout << "cell defect = " << defCell << endl;
      subGroup = vTable->getGroup("vReference");
      refPotFile = subGroup->get("file")->toString();
      refPot = getPot (refCell, refMesh, refPotFile, fileType, &structRef);
      if (refIsEV) refPot /= HA2EV;
      cout << "cell bulk = " << refCell << endl;
      printFullVLine = true;

      // --- parameters
      subGroup = mainGroup->getGroup("parameters");
      eCut = subGroup->get("eCut")->toReal ();

      // --- model charge
      subGroup = mainGroup->getGroup("modelCharge");
      if (subGroup->contains("dielecConstant")) 
         dielecConstant = subGroup->get("dielecConstant")->toReal ();
      else 
         dielecConstant = 1.0;

      q = 0;
      pos = Coord(0.,0.,0);
      double weight = 0.0;
      for (SxSymbolTable *gaussian = subGroup->getGroup ("gaussian");
           gaussian != NULL;
           gaussian = gaussian->nextSibling ("gaussian"))
      {
         double thisQ = gaussian->get("electrons")->toReal ();
         pos += fabs(thisQ) * Coord(gaussian->get("pos")->toList ());
         weight += fabs(thisQ);
         q += thisQ;
      }
      pos /= weight;
   } else  {
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
      } else if (inFile.getSize () > 0)  {
         SxParser parser;
         SxParser::Table table = parser.read (inFile,"std/structure.std");
         try {
            defCell = SxCell(&*table);
         } catch (SxException e)  {
            e.print ();
            SX_EXIT;
         }
      } else if (cli.groupAvailable (clInput))  {
         defCell = SxCell (col1, col2, col3);
      } else  {
         cout << "No cell provided, specify electrostatic potentials or SPHinX structure/inputfile, or"
            " use command line option!" << endl;
         SX_QUIT;
      }
   }

   SxFFT::quickFFTPlanner ();

   SxPtr<SxRBasis> RDef, RRef;
   if (defMesh.product () != 0)  {
      RDef = RDef.create (defMesh, defCell);
      defPot.setBasis(&*RDef);
   }
   if (refMesh.product () != 0)  {
      RRef = RRef.create (refMesh, refCell);
      refPot.setBasis(&*RRef);
   }

   cout << "Excess Electrons = " << q << endl;
   cout << "Located at " << pos << endl; 

   SxCell recCell = defCell.getReciprocalCell ();

   if (relative) pos = defCell.relToCar (pos);

   // --- setup atomic structures
   cout << "Atomic structure ";
   if (structDefFile.getSize () > 0 && structRefFile.getSize () > 0)  {
      SxParser parser;
      structRef = SxAtomicStructure (&*parser.read(structRefFile, "std/structure.std"));
      structDef = SxAtomicStructure (&*parser.read(structDefFile, "std/structure.std"));
   } else if (structDef.getNAtoms () == 0 && structRef.getNAtoms () == 0) {
      cout << "not ";
   }
   cout << "specified." << endl;

   // --- compute energies
   double eIso = 1.;

   // --- isolated
   double lastVal = 0.;
   if (q2tf < 1e-10)  {
      for (double g = step; g*g < eCut; g += 2. * step)  {
         eIso += 4. * sqr(rhoRec(sqr(g), beta2, gamma2, expNorm));
         eIso += 2. * (lastVal = sqr(rhoRec(sqr(g+step), beta2, gamma2, expNorm)));
      }
      eIso -= lastVal;
      eIso *= sqr(q) * step / (3. * PI) ;
   } else {
      for (double g = step; g*g < eCut; g += 2. * step)  {
         eIso += 4. * sqr(rhoRec(sqr(g), beta2, gamma2, expNorm))
              * SxPreconditioner::invCSRB(sqr(g),q2tf,dielecConstant);
         eIso += 2. * (lastVal = sqr(rhoRec(sqr(g+step), beta2, gamma2, expNorm)))
              * SxPreconditioner::invCSRB(sqr(g+step),q2tf,dielecConstant);
      }
      eIso -= lastVal;
      eIso *= sqr(q) * step / (3. * PI) ;
   }
   if (tensor) eIso /= dielecConstant;

   // --- periodic
   double ePeriodic = 0.;
   SxMesh3D mesh = SxGBasis::getMeshSize (eCut, defCell);
   if (defMesh.product () == 0) defMesh = mesh;
   int meshSize = mesh.product (), ng = 0;
   for (int i = 1; i < meshSize; ++i)  {
      Coord g = recCell.relToCar(mesh.getMeshVec(i, SxMesh3D::Origin));
      double g2 = g.normSqr ();
      if (g2 < eCut)  {
         if (q2tf < 1e-10)
            ePeriodic += sqr(rhoRec(g2, beta2, gamma2, expNorm))
                       / (g ^ epsTensor ^ g);
         else
            ePeriodic += sqr(rhoRec(g2, beta2, gamma2, expNorm)) / g2
                       * SxPreconditioner::invCSRB(sqr(g2),q2tf,dielecConstant);
         ng++;
      }
   }
   cout << "ng=" << ng << endl;
   ePeriodic *= sqr(q) * TWO_PI / defCell.volume;
   if (q2tf < 1e-10)  {
      ePeriodic += sqr(q) * FOUR_PI / defCell.volume
                   * rhoRecLimit0 (beta2, gamma2, expNorm)
                   / (tensor ? dielecConstant : 1.);
   } else {
      ePeriodic += sqr(q) * FOUR_PI / defCell.volume
                   * rhoRecLimit0 (beta2, gamma2, expNorm)
                 / dielecConstant;

      eIso *= dielecConstant;
      ePeriodic *= dielecConstant;
   }


   if (structDef.getNAtoms () > 0 && structRef.getNAtoms () > 0)  {
      structRef.epsEqual = 0.5;
      computeStruct (structDef, structRef, betaAvg, defPot, refPot, defMesh,
                     refMesh, pos, eCut, q, dielecConstant, epsTensor, beta2, 
                     gamma2, expNorm, atomField);
   }

   // --- compute potential
   for (int idir = 0; idir < 3; idir++)  {
      int nz = defMesh(idir);
      double dz = defCell.basis(idir).norm () / nz;
      SxFFT1d fft1d (SxFFT::Forward, nz);
      double dg = recCell.basis(idir).norm ();
      double epsEff;
      if (tensor)
        epsEff = recCell.basis(idir) ^ epsTensor ^ recCell.basis(idir)
               / (dg * dg);
      else
        epsEff = dielecConstant;
      Coord posRel = defCell.carToRel (pos);
      SxVector<Complex16> vInG(nz), vInR(nz);
      vInG(0) = 0.;
      vInG(0) = FOUR_PI * q / dielecConstant 
              * rhoRecLimit0 (beta2, gamma2, expNorm);
      if (idir == 0)  {
         cout << "V average: " <<  (vInG(0).re/defCell.volume * HA2EV) << " eV"
              << endl;
      }
      for (int z = 1; z < nz; ++z)  {
         double g = ((2 * z < nz) ?  z : (z - nz)) * dg;
         vInG(z) = FOUR_PI / ( epsEff * sqr(g) )
            * q * rhoRec(sqr(g), beta2, gamma2, expNorm)
            * exp ( - TWO_PI * I * g/dg * posRel(idir));
         if (q2tf > 1e-10)
            vInG(z) *= dielecConstant
               * SxPreconditioner::invCSRB(sqr(g),q2tf,dielecConstant);
      }
      if (nz % 2 == 0) vInG(nz/2) = 0.;
      fft1d.fftForward (nz, vInG.elements, vInR.elements);
      vInR /= defCell.volume;
      vInR += vAlign;


      if (printRho)  {
         SxVector<Complex16> rhoInG(nz), rhoInR(nz);
         for (int z = 0; z < nz; ++z)  {
            double g = ((2 * z < nz) ?  z : (z - nz)) * dg;
            rhoInG(z) = rhoRec(sqr(g), beta2, gamma2, expNorm)
                      * exp ( - TWO_PI * I * g/dg * posRel(idir));
         }
         if (nz % 2 == 0) rhoInG(nz/2) = 0.;
         fft1d.fftForward (nz, rhoInG.elements, rhoInR.elements);
         rhoInR /= defCell.volume;
         FILE *fp = fopen ("rhoModel.dat","w");
         if (!fp)  {
            cout << "Cannot open rhoModel.dat for writing" << endl;
            SX_EXIT;
         }
         for (int z = 0; z < nz; ++z)  {
            fprintf(fp, "%f\t%.12g\n", z * dz, rhoInR(z).re);
         }
         fclose (fp);

      }

      if (avgWidth > 1e-16)  {
         cout << "Averaging (" << (avgWidth / dz) << " points)" << endl;
         vInR = average (vInR.real (), avgWidth / dz);
      }
   
      SxString filename = "vline-eV-a"+SxString(idir)+".dat";
      FILE *fp = fopen (filename.ascii(),"w");
      if (!fp)  {
         cout << "Cannot open vline.dat for writing" << endl;
         SX_EXIT;
      }
      
      for (int z = 0; z < nz; ++z)  {
         fprintf(fp, "%f\t%f\n", z * dz, vInR(z).re  * HA2EV);
      }
   
      // --- read potential files
      // at this point only vRef has to be read in
      if (printFullVLine)  {
         SxVector<Double> vRef, vDef;
         vRef = readLine (refCell, refMesh, refPot, idir, nz, defCell,refPotFile);
         vDef = readLine (defCell, defMesh, defPot, idir, nz, defCell,defPotFile);
         if (avgWidth > 1e-16)  {
            vRef = average (vRef, avgWidth / dz);
            vDef = average (vDef, avgWidth / dz);
         }

         fprintf(fp, "&\n");
         for (int z = 0; z < nz; ++z)  {
            fprintf(fp, "%f\t%f\t%f\n", z * dz, 
                  (vDef(z) - vRef(z)) * HA2EV,
                  (vDef(z) - vRef(z) - vInR(z).re) * HA2EV);
         }
         if (printOrigPot)  {
            fprintf(fp, "&\n");
            for (int z = 0; z < nz; ++z)  {
               fprintf(fp, "%f\t%f\t%f\n", z * dz, 
                     vDef(z) * HA2EV, vRef(z) * HA2EV);
            }
         }
      }

      fclose (fp);
   }

   cout << "vAlign=" << (vAlign *HA2EV) << " eV" << endl;

   // --- report energies
   cout << SX_SEPARATOR;
   cout << "=== Intermediate results (" << (tensor ? "" : "un")
        << "screened) ===" << endl;
   cout << "Isolated energy       : " << eIso << endl;
   cout << "Periodic energy       : " << ePeriodic << endl;
   cout << "Difference (Hartree)  : " << ePeriodic - eIso << endl;
   cout << "Difference (eV)       : " << (ePeriodic - eIso) * HA2EV << endl;

   cout << SX_SEPARATOR;
   if (tensor)  {
      cout << "Calculation performed with epsilon = " << epsTensor << endl
           << "Spherical harmonic average         = ";

   } else
      cout << "Calculation performed with epsilon = ";
   cout << dielecConstant << endl;

   cout << SX_SEPARATOR;
   cout << "Defect correction (eV): ";
   if (tensor)
      cout << ((eIso - ePeriodic)  - q * vAlign) * HA2EV;
   else
      cout << ((eIso - ePeriodic) / dielecConstant - q * vAlign) * HA2EV;
   cout << " (incl. screening & alignment)" << endl;
}
