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

#include <SxHubbardMO.h>
#include <SxRotation.h>
#include <SxPAWOverlap.h>
#include <SxSimpleParser.h>
#include <SxProjector.h>
#include <SxNeighbors.h>
#include <SxHubbardU.h>
#include <SxTimer.h>
#include <SxFermi.h>

// Reference 1: SPHInX Hubbard U implementation notes
namespace Timer {
   enum HubbardMOTimers { MOCompute, AOMOForce };
}

SX_REGISTER_TIMERS(Timer::HubbardMOTimers)
{
   using namespace Timer;
   regTimer (MOCompute, "Hub. MO compute");
   regTimer (AOMOForce, "Hub. AO-MO force");
}

SxHubbardMO::SxHubbardMO (int siteOffsetIn)
   : nSite(0), siteOffset(siteOffsetIn), distFrom (-1.), dDist (0.),
     nAoPerSite (0), nInterpolate (100), nRad(200), verbose(false)
{
   //empty
}

SxDiracVec<Double>
SxHubbardMO::truncateShapeG (const SxDiracVec<Double> &shape,
                           double rCut,
                           double width) const
{
   const SxRadGBasis *radG
      = dynamic_cast<const SxRadGBasis*>(shape.getBasisPtr ());
   SX_CHECK (radG);
   SxRadRBasis R(0., rCut + 6. * width, nRad);
   const SxDiracVec<Double> &r = R.getRadRFunc ();

   SxDiracVec<Double> shapeR = R | shape;
   if (verbose) {
     SxBinIO(prefix + "shape.dat", SxBinIO::ASCII_WRITE_ONLY)
        .writeXYPlot (r, shapeR);
   }

   // ref. 1, Eq. (19)
   SX_LOOP(i)
      shapeR(i) *=  0.5 * erfc((r(i) - rCut)/width);

   if (verbose) {
      SxBinIO(prefix + "proj.dat", SxBinIO::ASCII_WRITE_ONLY)
         .writeXYPlot (r, shapeR);
   }

   return (*radG) | shapeR;
}


void SxHubbardMO::setupBox (double rSize, double eCut)
{
   SxCell &setupCell = setupStructure.cell;
   setupCell = CellMat (rSize, 0., 0., 0., rSize, 0., 0., 0., rSize);
   if (verbose)
      cout << "Setup cell for Hubbard MO: " << setupCell << endl;

   unsigned oldFFTMode;
   SxFFT::quickFFTPlanner (SxFFT::Estimate, &oldFFTMode);
   SxMesh3D setupMesh = SxGBasis::getMeshSize (eCut, setupCell);
   Coord k = setupCell.getReciprocalCell ().relToCar (Coord(0.25, 0.25, 0.25));
   setupG.set (setupMesh, setupCell, eCut, k);
   SxFFT::restorePlannerMode (oldFFTMode);

   setupG.structPtr = &setupStructure;
}

void SxHubbardMO::setupAO (const SxDiracVec<Double> &ao,
                           double rCut, double truncWidth,
                           const SxPtr<SxPAWPot> &pawPotPtr,
                           double eCut)
{
   SX_CHECK (shapeProj.getSize () == 0);
   // radial G basis
   double dg = 0.003;
   double gMax = sqrt(eCut)*1.1;
   int ngRad = int(gMax/dg) + 1;
   radGPtr = SxPtr<SxRadGBasis>::create (0., gMax, ngRad, SxRadGBasis::Linear);
   SxRadGBasis &radG = *radGPtr;

   // get shape in our radial G basis
   SxDiracVec<Double> shape = radG | ao;
   shapeAO = radG.toSpline (shape);

   // --- get truncated shape (without PAW normalization)
   shapeProj = truncateShapeG (shapeAO, rCut, truncWidth);

   if (verbose)  {
      SxBinIO(prefix + "shapeTruncG.dat", SxBinIO::ASCII_WRITE_ONLY)
         .writeXYPlot (radG.getRadGFunc (), shapeProj);
   }

   // --- add atoms with correct species to the setup structure
   int iSpecies = ao.handle->auxData.is;
   if (setupStructure.getNAtoms () > 0)  {
      setupG.changeTau (setupStructure);
      setupG.setupRealYlm (pawPotPtr->lMax (iSpecies));
   }

   l = ao.handle->auxData.l;
   mMO = abs(ao.handle->auxData.m);

   // --- setup PAW projectors in radG and setupG
   int npt = pawPotPtr->getNProjType (iSpecies);
   int npl = pawPotPtr->getNProj (iSpecies);
   // all projectors (types) in radial G space
   SxDiracVec<Double> allProjRadG;
   allProjRadG.reformat (radG.getNElements (), npt);
   // all projectors in setupG space
   PsiG projSpec;
   if (setupStructure.getNAtoms () > 0)  {
      projSpec.reformat (setupG.getNElements (), npl);
      projSpec.setBasis (setupG);
   }
   // <p|shapeProj> for single-atom normalization of Hubbard projector
   SxDiracVec<Double> pTrunc(npt);
   for (int ipt = 0, ipl = 0; ipt < npt; ipt++) {
      // get PAW projector
      int lProj = pawPotPtr->lPhi(iSpecies)(ipt);
      SxDiracVec<Double> projRadG = pawPotPtr->pPS(iSpecies).colRef (ipt);
      projRadG.handle->auxData.l = lProj;
      projRadG.handle->auxData.ia = -1; // no phase factor
      // project to radial G space
      projRadG = radG | projRadG;
      allProjRadG.colRef (ipt) <<= projRadG;
      // --- project to setupG space
      if (setupStructure.getNAtoms () > 0)  {
         for (int m = -lProj; m <= lProj; m++, ipl++)  {
            projRadG.handle->auxData.m = m;
            projSpec.colRef (ipl) <<= setupG | projRadG;
         }
      }

      // ref. 1, Eq. (20)
      // --- compute <p|shapeProj>
      if (l == lProj)
         pTrunc(ipt) = tr(projRadG * shapeProj);
      else
         pTrunc(ipt) = 0.;
   }

   // --- compute S|shapeProj> with S for single atom
   // ref. 1, Eq. (20)
   shapeProj += allProjRadG ^ pawPotPtr->deltaS(iSpecies) ^ pTrunc;

   if (verbose) {
      SxRadRBasis R(0., rCut + 6. * truncWidth, nRad);
      const SxDiracVec<Double> &r = R.getRadRFunc ();

      SxDiracVec<Double> shapeR = R | shapeProj;
      SxBinIO(prefix + "projS.dat", SxBinIO::ASCII_WRITE_ONLY)
         .writeXYPlot (r, shapeR);
      SxBinIO(prefix + "shapeTruncSG.dat", SxBinIO::ASCII_WRITE_ONLY)
         .writeXYPlot (radG.getRadGFunc (), shapeProj);
   }

   if (setupStructure.getNAtoms () > 0)  {
      // --- setup PAW overlap operator for single molecule
      SxPtr<SxAOBasis> proj = SxPtr<SxAOBasis>::create (projSpec, iSpecies);
      // --- setup refOrbMap (needed in SxPAWOverlap)
      proj->refOrbMap.resize (pawPotPtr->getNSpecies ());
      proj->refOrbMap(iSpecies).resize (npl);
      for (int ipt = 0, ipl = 0; ipt < npt; ipt++) {
         int lProj = pawPotPtr->lPhi(iSpecies)(ipt);
         for (int m = -lProj; m <= lProj; m++, ipl++)  {
            proj->refOrbMap(iSpecies)(ipl) = SxAOBasis::AoIndex (ipt, lProj, m);
         }
      }
      pBasis = SxPtr<SxPartialWaveBasis>::create (pawPotPtr, setupStructure);
      pBasis->projectors = proj;
   } else {
      // --- normalize projector
      double normAO = tr (shape * shape);
      // PAW norm corrections
      SxDiracVec<Double> pAO(npt);
      SX_LOOP(ipt)  {
         // get PAW projector
         int lProj = pawPotPtr->lPhi(iSpecies)(ipt);

         // --- compute <p|AO>
         if (l == lProj)
            pAO(ipt) = tr(allProjRadG.colRef (ipt) * shape);
         else
            pAO(ipt) = 0.;
      }
      normAO += dot(pAO, pawPotPtr->deltaS(iSpecies) ^ pAO);
      shapeProj /= tr (shapeProj * shape) / sqrt(normAO);
   }

   // spline the AO projector
   shapeProj = radG.toSpline (shapeProj);

}

SxVector<Double> SxHubbardMO::getRot (const Coord &axis) const
{
   SX_CHECK (axis.normSqr () > 1e-12);
   SxMatrix<Double> rotL; // rotation matrix from z-axis to molecular axis
                          // for the spherical harmonics at our L
   Coord molecularAxis = axis / axis.norm ();
   Coord rotAxis = molecularAxis.x (Coord(0.,0.,1.));
   double rotAxisNorm = rotAxis.norm ();
   SymMat rotZ;
   if (rotAxisNorm < 1e-8)  {
      rotZ = SymMat(1., 0., 0., 0., 1., 0., 0., 0., 1.);
   } else {
      rotZ = SxRotation (rotAxis / rotAxisNorm, acos(molecularAxis(2)));
   }
   SxYlm::SxClebschTable clebsch
      = SxYlm::getClebschGordan (max(l,1), max(l,1), 1, SxYlm::ComplexYlm);
   rotL = SxYlm::computeYlmRotMatrices(l, rotZ, clebsch)(l);

   SxVector<Double> res;
   res.reformat (2, 2 * l + 1);
   for (int M = 0; M < 2*l+1; M++)  {
      res(0,M) = rotL(M, l - mMO);
      res(1,M) = rotL(M, l + mMO);
   }
   return res;
}


void SxHubbardMO::setupNormalization (double distFromIn,
                                      double distTo)
{
   distFrom = distFromIn;
   SX_CHECK (setupStructure.getNAtoms () == 2, setupStructure.getNAtoms ());
   SX_CHECK (dynamic_cast<const SxRadGBasis*> (shapeAO.getBasisPtr ()));
   SX_CHECK (dynamic_cast<const SxRadGBasis*> (shapeProj.getBasisPtr ()));
   SX_CHECK (fabs(fabs(sign) - 1.) < 1e-12, sign);

   SxOverlap S (SxPtr<SxPAWOverlap>::create (pBasis, pBasis->getPotPtr ()));

   SxVector<Double> distData(nInterpolate);
   dDist = (distTo-distFrom)/double(nInterpolate - 1);
   setupStructure.ref(0) = Coord(0., 0., 0.);
   PsiG aoG = setupG | shapeAO;
   PsiG projG = setupG| shapeProj;
   SX_LOOP(i)  {
      double dist = distFrom +  int(i) * dDist;
      // ref 1, Eq. (24)
      PsiG phase = setupG.getPhaseFactors(Coord(0.,0.,dist)) + sign;
      PsiG mo = aoG * phase;
      PsiG moProj = projG * phase;
      setupStructure.ref(1) = Coord(0., 0., dist);
      setupG.changeTau (setupStructure);
      // ref. 1, Eq. (28)
      distData(i) = (mo | S | mo).re / dot(moProj, mo).absSqr ();
      //distData(i) = (mo | S | mo).re;
      //distData(i) = dot(moProj, mo).absSqr ();
   }
   pNorm.compute (distData, true);
}

bool SxHubbardMO::validateDist (double dist) const
{
   SX_CHECK (dDist > 0., dDist);
   if (dist < distFrom)  {
      cout << "Warning: atomic distance d=" << dist
           << " < distMin=" << distFrom << endl;
      return false;
   }
   int pNormSize = (int)pNorm.getCoeff()(0).getSize ();
   double distMax = distFrom + dDist * (pNormSize - 1);
   if (dist > distMax) {
      cout << "Warning: atomic distance d=" << dist
           << " > distMax=" << distMax << endl;
      return false;
   }
   return true;
}

void SxHubbardMO::finalize ()
{
   // --- destroy SxGBasis in place
   setupG.~SxGBasis (); // destruct
   new (&setupG) SxGBasis (); // reconstruct an empty object

   setupStructure.resize (0);

   // destroy vectors
   shapeAO = SxDiracVec<Double> ();
   shapeProj = SxDiracVec<Double> ();
   // destroy basis sets
   radGPtr = SxPtr<SxRadGBasis> ();
   pBasis = SxPtr<SxPartialWaveBasis> ();

}

static void checkIot (int iot, int nTypes, const SxString& what, int is)
{
   if (iot >= nTypes)  {
      cout << "In HubbardU.MO:" << endl;
      cout << "Invalid iot=" << iot << endl;
      cout << what << " has only " << nTypes
           << "orbital types for species " << (is+1) << endl;
      cout << "Counting iot starts at 0." << endl;
      SX_QUIT;
   }
}

void SxHubbardMO::read (const SxSymbolTable *table,
                        const SxAtomicStructure &structure,
                        const SxPtr<SxPAWPot> &potPtr)
{
   SxAtomicOrbitals aoQuamol;
   SxDiracVec<Double> ao;
   double rCut = 0., truncWidth = 0.7;
   double distTo = -1.;
   double rSize = -1.;
   double eCut = SxGBasis::getECut (table->topLevel ());
   int signIn = 0;
   int is = findMolecules (table, structure, *potPtr);
   SYMBOLPARSE(table)  {
      SYMBOLGROUP ("MO") {
         verbose = SYMBOLGET("verbose").toBool ();
         SYMBOLGROUP ("orbital")  {
            bool fromPot = SYMBOLGET("fromPotential").toBool ();
            int iot = SYMBOLGET("iot");
            if (fromPot)  {
               checkIot (iot, potPtr->getNProjType (is), "Potential",
                         is);
               ao = potPtr->phiPS(is).colRef (iot);
               ao.handle->auxData.l = potPtr->lPhi(is)(iot);
            } else {
               SxString quamolFile = SYMBOLGET("file");
               aoQuamol.read (quamolFile);
               int isFile = SYMBOLGET("is") || is;
               checkIot (iot, aoQuamol.getNOrbTypes (), quamolFile, isFile);
               ao.handle->auxData.l = aoQuamol.getL (is, iot);
               ao = aoQuamol(isFile, iot);
            }
         }
         mMO = SYMBOLGET("mMO");
         if (mMO > ao.handle->auxData.l)  {
            cout << "Invalid mMO=" << mMO << "from atomic orbital with l="
                 << ao.handle->auxData.l << endl;
            SX_QUIT;
         }
         rCut = SYMBOLGET("rCut");
         truncWidth = SYMBOLGET("cutWidth") || 0.7;
         distTo = SYMBOLGET("maxDist");
         distFrom = SYMBOLGET("minDist") || 0.5 * distTo;
         nInterpolate = SYMBOLGET("nInterpolate") || 100;
         nRad = SYMBOLGET("nRadGrid") || 200;
         signIn = SYMBOLGET("sign");
         rSize = SYMBOLGET("setupBoxSize") || 40.;
         if (abs(signIn) != 1)  {
            cout << "Invalid MO sign=" << signIn << endl;
            cout << "Only +1 or -1 allowed" << endl;
            SX_QUIT;
         }
      } else {
         SX_EXIT;
      }
   }
   sign = signIn;
   ao.handle->auxData.is = is;
   ao.handle->auxData.m = abs(mMO);
   setupBox (rSize, eCut);
   for (int iSpecies = 0; iSpecies < structure.getNSpecies (); iSpecies++)
      setupStructure.newSpecies ();
   // two atoms in setup structure
   SX_CHECK (setupStructure.getNAtoms () == 0, setupStructure.getNAtoms ());
   setupStructure.addAtom (is, Coord(0.,0.,0.));
   setupStructure.addAtom (is, Coord(0.,0.,1.));
   setupStructure.endCreation ();

   setupAO (ao, rCut, truncWidth, potPtr, eCut);
   nAoPerSite = 2 * (2 * l + 1);
   setupNormalization (distFrom, distTo);
}

void SxHubbardMO::readAO (const SxSymbolTable *table,
                          const SxAtomicStructure &structure,
                          const SxPtr<SxPAWPot> &potPtr)
{
   SxAtomicOrbitals aoQuamol;
   SxDiracVec<Double> ao;
   double rCut = 0., truncWidth = 0.7;
   double eCut = SxGBasis::getECut (table->topLevel ());
   int is = findAtoms (table, structure, *potPtr);
   SYMBOLPARSE(table)  {
      SYMBOLGROUP ("AO") {
         verbose = SYMBOLGET("verbose").toBool ();
         SYMBOLGROUP ("orbital")  {
            bool fromPot = SYMBOLGET("fromPotential").toBool ();
            int iot = SYMBOLGET("iot");
            if (fromPot)  {
               checkIot (iot, potPtr->getNProjType (is), "Potential",
                         is);
               ao = potPtr->phiPS(is).colRef (iot);
               ao.handle->auxData.l = potPtr->lPhi(is)(iot);
            } else {
               SxString quamolFile = SYMBOLGET("file");
               aoQuamol.read (quamolFile);
               int isFile = SYMBOLGET("is") || is;
               checkIot (iot, aoQuamol.getNOrbTypes (), quamolFile, isFile);
               ao.handle->auxData.l = aoQuamol.getL (is, iot);
               ao = aoQuamol(isFile, iot);
            }
         }
         rCut = SYMBOLGET("rCut");
         truncWidth = SYMBOLGET("cutWidth") || 0.7;
         nRad = SYMBOLGET("nRadGrid") || 200;
      } else {
         SX_EXIT;
      }
   }
   ao.handle->auxData.is = is;
   ao.handle->auxData.m = 0;
   setupAO (ao, rCut, truncWidth, potPtr, eCut);
   nAoPerSite = 2 * l + 1;
}


static void inputProblem (const SxSymbolTable *table)
{
   cout << "Hubbard U in "
        << table->parserFilename
        << " line " << table->parserLineNumber
        << ":";
}

int SxHubbardMO::findMolecules (const SxSymbolTable* moGroup,
                                const SxAtomicStructure &structure,
                                const SxSpeciesData &speciesInfo)
{
   SxGrid grid(structure,10);
   SxNeighbors nn;
   int iSpecies = -1;
   SxString label;
   const SxArray<SxString> *labels = NULL;
   double maxDist;
   SYMBOLPARSE(moGroup)  {
      SxString element = SYMBOLGET("element") || "";
      iSpecies << SYMBOLGET("species");
      iSpecies--;
      label = SYMBOLGET("label") || "";
      if (element.getSize () > 0)  {
         iSpecies = speciesInfo.find (element);
      }
      if (label.getSize () > 0)  {
         if (!structure.hasLabels ())  {
            inputProblem (SYMBOLGROUP_TABLE);
            cout << "No labels defined in structure" << endl;
            SX_QUIT;
         }
         labels = &structure.getLabels ();
      }
      maxDist = SYMBOLGET("maxDist");
   }
   SxStack<int> atomId;
   SxArray<bool> atomDone(structure.getNAtoms ());
   SxStack<Coord> deltaR;
   atomDone.set (false);
   for (int iTlAtom = (iSpecies >= 0) ? structure.getIAtom (iSpecies, 0)
                                      : 0;
        iTlAtom < structure.getNAtoms ();
        iTlAtom++)
   {
      if (atomDone(iTlAtom)) continue;
      if (labels && (*labels)(iTlAtom) == label)  {
         int is = structure.getISpecies (iTlAtom);
         if (iSpecies >= 0 && iSpecies != is)  {
            cout << "Species mismatch: species for labeled atom "
                 << (iTlAtom+1) << " is " << (is+1)
                 << "(" << speciesInfo.chemName(is) << "), but should be "
                 << (iSpecies + 1) << "(" << speciesInfo.chemName(iSpecies)
                 << ")." << endl;
            SX_QUIT;
         }
         iSpecies = is;
      }
      atomDone(iTlAtom)=true;
      if (iSpecies < 0) continue; // we rely on labels, but haven't found any
                                  // atom yet
      nn.compute (grid, structure, structure.getAtom (iTlAtom), maxDist,
                  SxNeighbors::StoreIdx | SxNeighbors::StoreAbs);
      if (nn.getSize () == 0) continue;
      int jTl = -1;
      bool haveProblem = false;
      SX_LOOP(in)  {
         int jTlAtom = nn.idx(in);
         if (labels && (*labels)(jTlAtom) != label) continue;
         if (structure.getISpecies (jTlAtom) != iSpecies) continue;
         if (jTl > -1)  {
            cout << "WARNING: atom " << (iTlAtom+1) << " @ "
                 << structure.constRef (iTlAtom)
                 << " has multiple neighbors within maxDist=" << maxDist << ':'
                 << endl;
            cout << "atom " << (jTl + 1) << " @ " << structure.getAtom (jTl)
                 << endl;
            jTl=-1;
            haveProblem = true;
         }
         if (haveProblem)  {
            // report atom
            cout << "atom " << (jTlAtom + 1) << " @ "
                 << structure.getAtom (jTlAtom) << endl;
         } else {
            jTl = jTlAtom;
            deltaR << Coord(0.,0.,0.);
            deltaR << nn.absPositions(in) - structure.getAtom (jTl);
         }
      }
      if (haveProblem) {
         deltaR.pop ();
         deltaR.pop ();
         continue;
      }
      if (jTl < 0) continue;
      cout << "Found " << speciesInfo.chemName(iSpecies)
           << "2 molecule with atoms "
           << (iTlAtom+1) << " and " << (jTl+1) << '.' << endl;
      atomId << iTlAtom;
      atomId << jTl;
      atomDone(jTl) = true;
   }
   atomIdx = atomId;
   nSite = (int)atomIdx.getSize () / 2;
   latShift = deltaR;
   return iSpecies;
}

int SxHubbardMO::findAtoms(const SxSymbolTable* aoGroup,
                           const SxAtomicStructure &structure,
                           const SxSpeciesData &speciesInfo)
{
   int iSpecies = -1;
   SxString label;
   const SxArray<SxString> *labels = NULL;
   SYMBOLPARSE(aoGroup)  {
      SxString element = SYMBOLGET("element") || "";
      iSpecies << SYMBOLGET("species");
      iSpecies--;
      label = SYMBOLGET("label") || "";
      if (element.getSize () > 0)  {
         iSpecies = speciesInfo.find (element);
      }
      if (label.getSize () > 0)  {
         if (!structure.hasLabels ())  {
            inputProblem (SYMBOLGROUP_TABLE);
            cout << "No labels defined in structure" << endl;
            SX_QUIT;
         }
         labels = &structure.getLabels ();
      }
      if (!labels)  {
         // --- all atoms of this species
         if (iSpecies < 0)  {
            inputProblem (SYMBOLGROUP_TABLE);
            cout << "no species or labels defined" << endl;
            SX_QUIT;
         }
         // --- generate list of atoms
         nSite = (int)structure.getNAtoms (iSpecies);
         atomIdx.resize (nSite);
         SX_LOOP (ia)
            atomIdx (ia) = structure.getIAtom (iSpecies, ia);
         return iSpecies;
      }
   }

   // --- labels
   SxStack<int> atomId;
   for (int iTlAtom = (iSpecies >= 0) ? structure.getIAtom (iSpecies, 0)
                                      : 0;
        iTlAtom < structure.getNAtoms ();
        iTlAtom++)
   {
      if ((*labels)(iTlAtom) != label) continue;
      int is = structure.getISpecies (iTlAtom);
      if (iSpecies >= 0 && iSpecies != is)  {
         cout << "Species mismatch: species for labeled atom "
              << (iTlAtom+1) << " is " << (is+1)
              << "(" << speciesInfo.chemName(is) << "), but should be "
              << (iSpecies + 1) << "(" << speciesInfo.chemName(iSpecies)
              << ")." << endl;
         SX_QUIT;
      }
      iSpecies = is;
      atomId << iTlAtom;
   }
   atomIdx = atomId;
   nSite = (int)atomIdx.getSize ();
   return iSpecies;
}

void SxHubbardMO::setupProjGk (const SxGkBasis &gk)
{
   aoProj = SxPtr<SxAOBasis>::create ();
   shapeProj.handle->auxData.m = NONE_M;
   // --- generate reference orbitals
   aoProj->refOrbitals.resize (gk.getNk ());
   if (latShift.getSize ())
      aoProj->extraPhase.resize (gk.getNk ());
   int nm = 2 * l + 1;
   int iSpecies = shapeProj.handle->auxData.is;
   const SxAtomicStructure *strPtr = NULL;
   for (int ik = 0; ik < gk.getNk (); ++ik) {
      if (!SxLoopMPI::myWork(ik)) continue;
      if (!strPtr) strPtr = gk(ik).structPtr;
      SX_CHECK (strPtr);
      aoProj->refOrbitals(ik).resize (strPtr->getNSpecies ());
      SxDiracVec<Complex16> &refOrbs = aoProj->refOrbitals(ik)(iSpecies);
      refOrbs.reformat (gk(ik).ng, nm);
      refOrbs.setBasis (gk(ik));
      // project from radial G to |G+k> (without Ylm)
      SxDiracVec<Double> shapeProjG = gk(ik) | shapeProj;
      for (int m = -l; m <= l; ++m)  {
         // --- multiply with Ylm
         SxDiracVec<Double> Ylm = gk(ik).getYlm (l,m);
         double norm = SxYlm::getYlmNormFactor(l,m);
         for (int ig = 0; ig < gk(ik).ng; ++ig)
            refOrbs(ig, l + m) = norm * Ylm(ig) * shapeProjG(ig);
      }
      if (latShift.getSize ())  {
         // --- phase correction for molecules across boundaries
         aoProj->extraPhase(ik).resize (atomIdx.getSize () * nm);
         for (ssize_t iAtom = 0, iPh = 0; iAtom < atomIdx.getSize (); ++iAtom)  {
            SxComplex16 phase = SxComplex16::phase(-gk.getK(ik) ^ latShift(iAtom));
            for (int im = 0; im < nm; ++im)
            aoProj->extraPhase(ik)(iPh++) = phase;
         }
         VALIDATE_VECTOR(aoProj->extraPhase(ik));
      }
   }
   if (!strPtr) { // this MPI task has no k-points...
      aoProj = SxPtr<SxAOBasis> (); // remove AO basis
      return;
   }
   aoProj->cacheRefOrb = SxAOBasis::CacheAll;
   // --- generate orbital map
   aoProj->orbitalMap.resize (nm * atomIdx.getSize ());
   int offset = strPtr->atomInfo->offset (iSpecies);
   for (int io = 0, ia = 0; ia < atomIdx.getSize (); ++ia)  {
      int iAtom = atomIdx(ia) - offset;
      SX_CHECK (strPtr->getIAtom(iSpecies, iAtom) == atomIdx(ia),
                iSpecies, iAtom, atomIdx(ia));
      for (int m = 0; m < nm; ++m, ++io)  {
         aoProj->orbitalMap(io) = SxAOBasis::OrbitalIndex(iSpecies, iAtom, m);
      }
   }
}

void SxHubbardMO::compute (SxHubbardU *hubbardU,
                           const SxBlockDensityMatrix& Pij,
                           const SxAtomicStructure &structure)
{
   if (latShift.getSize ())
      computeMO (hubbardU, Pij, structure);
   else
      computeAO (hubbardU, Pij, structure);
}


void SxHubbardMO::computeMO (SxHubbardU *hubbardU,
                           const SxBlockDensityMatrix& Pij,
                           const SxAtomicStructure &structure)
{
   SX_CLOCK(Timer::MOCompute);
   trafoForce.resize ((int)atomIdx.getSize ());
   trafoForce.set (0.,0.,0.);
   if (!trafoForce.atomInfo
       || !trafoForce.atomInfo->isChild (structure.atomInfo))
   {
      // --- our selection of atoms is a subset of the total structure
      //     => create a corresponding parentMap to allow for implicit
      //        remapping from our forces to the global forces
      SxPtr<SxAtomInfo> info = SxAtomInfo::derive (structure.atomInfo);
      info->parentMap.resize (atomIdx.getSize ());
      info->resize (structure.getNSpecies ());
      info->nAtoms = 0;
      SX_LOOP(ia)  {
         info->parentMap(ia) = (int)atomIdx(ia);
         info->nAtoms(structure.getISpecies (atomIdx(ia)))++;
      }
      info->setupOffset ();
      trafoForce.replaceInfo (info);
   }
   int nm = 2 * l + 1;
   double fFull = Pij.getNSpin () == 1 ? 2. : 1.;
   double fInv = 1. / fFull;
   // rescale energy to single spin
   hubbardU->energy *= fInv;
   hubbardU->eDoubleCounting *= fInv;
   SX_LOOP2(iSpin, iSite)  {
      Coord axis = structure(atomIdx(2*iSite+1))
                 - structure(atomIdx(2*iSite  ));
      axis += latShift(2*iSite + 1) - latShift(2*iSite);
      double dist = axis.norm ();
      if (!validateDist (dist))  { SX_EXIT; }
      SxVector<Double> orbToSite = getRot (axis);

      // d orbToSite / d axis(i)
      SxArray<SxVector<Double> > orbToSiteDeriv(3);
      double dx = 3e-6 * dist;
      for (int iDir = 0; iDir < 3; iDir++)  {
         SX_CLOCK (Timer::AOMOForce);
         Coord dir(0,0,0);
         dir(iDir) = dx;
         orbToSiteDeriv(iDir) = (getRot(axis+dir) - getRot(axis-dir))
                              / (2. * dx);
      }

      // contract from AO to MO space
      SxMatrix<Double> molD(nm,nm);
      molD.set (0.);
      const SxMatrix<Double> &D = Pij(iSpin, iSite);
      SX_CHECK (D.nCols () == 2 * nm, D.nCols (), 2 * nm);
      SX_CHECK (D.nRows () == 2 * nm, D.nRows (), 2 * nm);
      // sum over atoms (tau, tau') ref. 1, Eq. (33) using Eq. (26)
      for (int m1 = 0; m1 < nm; m1++)  {
         for (int m2 = 0; m2 < nm; m2++) {
            molD(m1,m2) =         D(m1     , m2) + D(m1 + nm, nm + m2)
                        + sign * (D(m1 + nm, m2) + D(m1     , nm + m2));
         }
      }
      //{
      //   SxMatrix<Complex16>::Eigensystem eig = SxMatrix<Complex16> (molD).eigensystem ();
      //   cout << "molD eigenspace:" << endl;
      //   for (int mm = 0; mm < nm; ++mm)
      //      cout << eig.vals(mm).re << ": " << eig.vecs.colRef(mm) << endl;
      //}
      // ref. 1, Eq. (33) rotational part
      SxVector<Double> hubbardP = orbToSite ^ molD ^ orbToSite.transpose ();
      // normalize
      double normDeriv = getPNormDeriv (dist);
      double projNorm = getPNorm (dist);
      // ref. 1, Eq. (33) normalization; also renormalize for nSpin=1 <=> f=2
      hubbardP *= projNorm * fInv;

      // compute HubbardU stuff
      SxVector<Double> Hij
         = hubbardU->computeIncr((int)iSite + siteOffset, hubbardP);

      SX_START_TIMER(Timer::AOMOForce);
      // transformation force
      SX_LOOP(iDir)  {
         double f;
         // --- ref. 1, Eq. (41a)
         f  = dot(orbToSite ^ molD ^ orbToSiteDeriv(iDir).transpose (), Hij);
         f += dot(orbToSiteDeriv(iDir) ^ molD ^ orbToSite.transpose (), Hij);
         f *= projNorm * fFull;
         trafoForce.ref(2*iSite+1)(iDir) += f;
         trafoForce.ref(2*iSite  )(iDir) -= f;
      }
      // normalization force along the molecular axis
      // ref. 1, Eq. (44)
      double fDist = (Hij * hubbardP).sum () * normDeriv / projNorm * fFull;
      Coord fNorm = axis * (fDist/dist);
      trafoForce.ref(2*iSite)   -= fNorm;
      trafoForce.ref(2*iSite+1) += fNorm;
      SX_STOP_TIMER (Timer::AOMOForce);

      // Hamiltonian in AO space
      // ref. 1, Eq. (37)
      Hij *= projNorm;
      Hij = orbToSite.transpose () ^ Hij ^ orbToSite;
      SxMatrix<Double> &H = hamProj(iSpin, iSite);
      SX_CHECK (H.nCols () == 2 * nm, H.nCols (), 2 * nm);
      SX_CHECK (H.nRows () == 2 * nm, H.nRows (), 2 * nm);
      for (int m1 = 0; m1 < nm; m1++)  {
         for (int m2 = 0; m2 < nm; m2++) {
            H(m1   ,m2   ) =        Hij(m1,m2);
            H(m1+nm,m2   ) = sign * Hij(m1,m2);
            H(m1   ,m2+nm) = sign * Hij(m1,m2);
            H(m1+nm,m2+nm) =        Hij(m1,m2);
         }
      }
   }
   hubbardU->energy *= fFull;
   hubbardU->eDoubleCounting *= fFull;
}

void SxHubbardMO::computeAO (SxHubbardU *hubbardU,
                             const SxBlockDensityMatrix& Pij,
                             const SxAtomicStructure &structure)
{
   SX_CLOCK(Timer::MOCompute);
   trafoForce.resize ((int)atomIdx.getSize ());
   if (!trafoForce.atomInfo
       || !trafoForce.atomInfo->isChild (structure.atomInfo))
   {
      // --- our selection of atoms is a subset of the total structure
      //     => create a corresponding parentMap to allow for implicit
      //        remapping from our forces to the global forces
      SxPtr<SxAtomInfo> info = SxAtomInfo::derive (structure.atomInfo);
      info->parentMap.resize (atomIdx.getSize ());
      info->resize (structure.getNSpecies ());
      info->nAtoms = 0;
      SX_LOOP(ia)  {
         info->parentMap(ia) = (int)atomIdx(ia);
         info->nAtoms(structure.getISpecies (atomIdx(ia)))++;
      }
      info->setupOffset ();
      trafoForce.replaceInfo (info);
      trafoForce.set (0.,0.,0.); // never changes
   }
   int nm = 2 * l + 1;
   double fFull = Pij.getNSpin () == 1 ? 2. : 1.;
   double fInv = 1. / fFull;
   // rescale energy to single spin
   hubbardU->energy *= fInv;
   hubbardU->eDoubleCounting *= fInv;
   SX_LOOP2(iSpin, iSite)  {
      const SxMatrix<Double> &D = Pij(iSpin, iSite);
      SX_CHECK (D.nCols () == nm, D.nCols (), nm);
      SX_CHECK (D.nRows () == nm, D.nRows (), nm);
      if (verbose) {
         SxMatrix<Complex16>::Eigensystem eig = SxMatrix<Complex16> (D).eigensystem ();
         cout << "occ eigenspace (site " << (iSite+1) << "):" << endl;
         for (int mm = 0; mm < nm; ++mm)
            cout << eig.vals(mm).re << ": " << eig.vecs.colRef(mm) << endl;
      }
      // compute HubbardU stuff from renormalized D (for nSpin=1 <=> f=2)
      hamProj(iSpin, iSite)
         = hubbardU->computeIncr((int)iSite + siteOffset, fInv * D);
   }
   hubbardU->energy *= fFull;
   hubbardU->eDoubleCounting *= fFull;
}

void SxHubbardMO::addToRho (SxBlockDensityMatrix *Pij,
                            const PsiG &waves,
                            double weight,
                            const SxDiracVec<Double> &focc) const
{
   SX_CHECK (Pij);
   SxDiracVec<Complex16> aoPsi = (*aoProj) | waves;
   int iSpin = waves.handle->auxData.iSpin;
   SX_CHECK (iSpin >=0 && iSpin < Pij->getNSpin (),
             iSpin, Pij->getNSpin ());
   for (ssize_t iState = 0; iState < waves.nCols (); ++iState)  {
      int offset = 0;
      // ref. 1, Eq. (31)
      SX_LOOP(iSite)  {
         SX_LOOP2(io,jo)  {
            (*Pij)(iSpin,iSite)(io,jo) += (  aoPsi(offset + io, iState).conj ()
                                           * aoPsi(offset + jo, iState)).re
                                           * focc(iState) * weight;
         }
         offset += nAoPerSite;
      }
      SX_CHECK (offset == aoPsi.nRows (), offset, aoPsi.nRows ());
   }
}

void SxHubbardMO::symmetrize (SxBlockDensityMatrix *Pij,
                              const SxAtomicStructure &structure,
                              const SxYlmRotGroup &ylmRot) const
{
   SX_CHECK (Pij->getNSite () == getNSite (),
             Pij->getNSite (), getNSite ());
   // --- TODO
   // symmetrization function in SxBlockDensityMatrix
   // sym (iSite, L1, offset1, L2, offset2, ylmRot, mapSite, mapSpin)
   // => symmetrize subblock 2L1+1 x 2L2+1 (at offset offset1,offset2)
   //    across all equivalent sites, rotating the sites and potentially
   //    the spin according to the maps

   // map iTlAtom -> atomIdx
   SxArray<int> mapAtoms(structure.getNAtoms ());
   mapAtoms.set (-1);
   SX_LOOP(i) mapAtoms(atomIdx(i))= (int)i;

   SX_CHECK (structure.cell.symGroupPtr);
   const SxSymGroup &syms = *structure.cell.symGroupPtr;
   int nSym = syms.getNSymmorphic ();
   int nAtomsPerMolecule = (latShift.getSize () > 0) ? 2 : 1;
   SxVector<Int> symAtomIdx;
   symAtomIdx.reformat(nAtomsPerMolecule,nSym);
   SxGrid grid (structure, 10);

   int nm = 2 * l + 1;
   for (int iSite = 0, idx = 0; iSite < getNSite ();
        ++iSite, idx += nAtomsPerMolecule)
   {

      // --- find equivalent atoms
      for (int iSym = 0; iSym < nSym; ++iSym)  {
         SymMat S = syms.getSymmorphic (iSym);
         for (int ia = 0; ia < nAtomsPerMolecule; ++ia)  {
            Coord pos = structure.getAtom(atomIdx(idx+ia));
            int jTl = structure.find (S ^ pos, grid);
            symAtomIdx(ia, iSym) = mapAtoms(jTl);
         }
         // the two atoms must belong to the same molecule (=site)
         SX_CHECK ((nAtomsPerMolecule == 1) ||
                   (symAtomIdx(0, iSym) / 2 == symAtomIdx(1, iSym) / 2),
                   symAtomIdx(0, iSym), symAtomIdx(1, iSym));
      }
      // check if symmetrization has been performed already
      if (symAtomIdx.minval () < idx) continue;

      SxMatrix<Double> sym (nm, nm), rotRho(nm,nm);
      for (int iSpin = 0; iSpin < Pij->getNSpin (); ++iSpin)  {
         for (int ia = 0; ia < nAtomsPerMolecule; ++ia)  {
            for (int ja = 0; ja < nAtomsPerMolecule; ++ja)  {
               sym.set (0.);
               // --- construct symmetrized molDensMat for this site
               for (int iSym = 0; iSym < nSym; ++iSym)  {
                  int rotSite = symAtomIdx(ia, iSym) / nAtomsPerMolecule;
                  int offset, offset2;
                  if (nAtomsPerMolecule == 1)  {
                     offset = offset2 = 0;
                  } else {
                     offset  = (symAtomIdx(ia, iSym) & 1) * nm;
                     offset2 = (symAtomIdx(ja, iSym) & 1) * nm;
                  }
                  // --- collect rotRho
                  const SxMatrix<Double> &rotD = (*Pij)(iSpin,rotSite);
                  const SxMatrix<Double> &Dl  = ylmRot(iSym)(l);
                  for (int m = 0; m < nm; ++m)
                     for (int m2 = 0; m2 < nm; ++m2)
                        rotRho(m,m2) = rotD(offset + m, offset2 + m2);
                  // --- symmetrize
                  // ref. 1, Eq. (45)
                  sym += Dl.transpose () ^ rotRho ^ Dl;
               }
               sym /= double(nSym);
               // --- distribute symmetrized matrix
               for (int iSym = 0; iSym < nSym; ++iSym)  {
                  // rotate with symmetry
                  const SxMatrix<Double> &Dl  = ylmRot(iSym)(l);
                  // ref. 1, Eq. (46)
                  rotRho = Dl ^ sym ^ Dl.transpose ();

                  int rotSite = symAtomIdx(ia, iSym) / nAtomsPerMolecule;
                  int offset, offset2;
                  if (nAtomsPerMolecule == 1)  {
                     offset = offset2 = 0;
                  } else {
                     offset  = (symAtomIdx(ia, iSym) & 1) * nm;
                     offset2 = (symAtomIdx(ja, iSym) & 1) * nm;
                  }
                  // write out
                  SxMatrix<Double> &rotD = (*Pij)(iSpin,rotSite);
                  for (int m = 0; m < nm; ++m)
                     for (int m2 = 0; m2 < nm; ++m2)
                        rotD(offset + m, offset2 + m2) = rotRho(m,m2);
               }
            }
         }
      }
   }

}

SxBlockDensityMatrix
SxHubbardMO::computeGradPij (const SxFermi &fermi,
                                   const SxDiracVec<Complex16> &P,
                                   const SxArray<SxDiracMat<Complex16> > &gradP) const
{
   // note: for gradients, we use the "spin" index of SxBlockDensityMatrix
   // for the direction
   SxBlockDensityMatrix res;
   res.resize (3, getNSite ());
   int nAO = getNAOperSite ();
   SX_LOOP2(iDir, iSite) {
      res(iDir, iSite).reformat (nAO, nAO);
      res(iDir, iSite).set (0.);
   }
   SX_CHECK (P.handle);
   int ik = P.handle->auxData.ik;
   int iSpin = P.handle->auxData.iSpin;
   SX_LOOP(iState)  {
      int offset = 0;
      double weight = fermi.kpPtr->weights(ik) * fermi.focc(iState, iSpin, ik);
      SX_LOOP(iSite)  {
         SX_LOOP(iDir)  {
            SxMatrix<Double> &D = res(iDir, iSite);
            // cf. definition of Q matrix in ref. 1, Eq. (39)
            SX_LOOP2(iAO, jAO)
               D(iAO, jAO) -= weight * (gradP(iDir)(offset+iAO,iState).conj()
                                        *         P(offset+jAO,iState)     ).re;
         }
         offset += nAO;
      }
      SX_CHECK (offset == P.nRows (), offset, P.nRows ());
   }
   return res;
}


void SxHubbardMO::symmetrizeGradPij (SxBlockDensityMatrix *Pij,
                                     const SxAtomicStructure &structure,
                                     const SxYlmRotGroup &ylmRot) const
{
   SX_CHECK (Pij->getNSite () == getNSite (),
             Pij->getNSite (), getNSite ());
   // note: for gradients, we use the "spin" index of SxBlockDensityMatrix
   // for the direction
   SX_CHECK(Pij->getNSpin () == 3, Pij->getNSpin ());
   // map iTlAtom -> atomIdx
   SxArray<int> mapAtoms(structure.getNAtoms ());
   mapAtoms.set (-1);
   SX_LOOP(i) mapAtoms(atomIdx(i))= (int)i;

   SX_CHECK (structure.cell.symGroupPtr);
   const SxSymGroup &syms = *structure.cell.symGroupPtr;
   int nSym = syms.getNSymmorphic ();
   int nAtomsPerMolecule = (latShift.getSize () > 0) ? 2 : 1;
   SxVector<Int> symAtomIdx;
   symAtomIdx.reformat(nAtomsPerMolecule,nSym);
   SxGrid grid (structure, 10);

   int nm = 2 * l + 1;
   for (int iSite = 0, idx = 0; iSite < getNSite ();
        ++iSite, idx += nAtomsPerMolecule)
   {

      // --- find equivalent atoms
      for (int iSym = 0; iSym < nSym; ++iSym)  {
         const SymMat &S = syms.getSymmorphic (iSym);
         for (int ia = 0; ia < nAtomsPerMolecule; ++ia)  {
            Coord pos = structure.getAtom(atomIdx(idx+ia));
            int jTl = structure.find (S ^ pos, grid);
            symAtomIdx(ia, iSym) = mapAtoms(jTl);
         }
         // the two atoms must belong to the same molecule (=site)
         SX_CHECK ((nAtomsPerMolecule == 1) ||
                   (symAtomIdx(0, iSym) / 2 == symAtomIdx(1, iSym) / 2),
                   symAtomIdx(0, iSym), symAtomIdx(1, iSym));
      }
      // check if symmetrization has been performed already
      if (symAtomIdx.minval () < idx) continue;

      SxMatrix<Double> rotRho(nm,nm);
      SxArray<SxMatrix<Double> > sym(3);
      for (int ia = 0; ia < nAtomsPerMolecule; ++ia)  {
         for (int ja = 0; ja < nAtomsPerMolecule; ++ja)  {
            SX_LOOP(iDir)  {
               sym(iDir).reformat (nm, nm);
               sym(iDir).set (0.);
            }
            for (int iSym = 0; iSym < nSym; ++iSym)  {
               int rotSite = symAtomIdx(ia, iSym) / nAtomsPerMolecule;
               int offset, offset2;
               if (nAtomsPerMolecule == 1)  {
                  offset = offset2 = 0;
               } else {
                  offset  = (symAtomIdx(ia, iSym) & 1) * nm;
                  offset2 = (symAtomIdx(ja, iSym) & 1) * nm;
               }
               const SymMat &S = syms.getSymmorphic (iSym);
               // --- collect rotRho
               const SxMatrix<Double> &Dl  = ylmRot(iSym)(l);
               for (int iDir = 0; iDir < 3; ++iDir)  {
                  const SxMatrix<Double> &rotD = (*Pij)(iDir,rotSite);
                  for (int m = 0; m < nm; ++m)
                     for (int m2 = 0; m2 < nm; ++m2)
                        rotRho(m,m2) = rotD(offset + m, offset2 + m2);
                  // ref. 1, Eq. (47), sum over kk'
                  rotRho = Dl.transpose () ^ rotRho ^ Dl;
                  // ref. 1, Eq. (47), sum over beta
                  for (int jDir = 0; jDir < 3; ++jDir)
                     sym(jDir).plus_assign_ax (S(iDir,jDir), rotRho);
               }
            }
            SX_LOOP(iDir) sym(iDir) /= double(nSym);
            // --- distribute symmetrized matrix
            for (int iSym = 0; iSym < nSym; ++iSym)  {
               // rotate with symmetry
               const SymMat &S = syms.getSymmorphic (iSym);
               const SxMatrix<Double> &Dl  = ylmRot(iSym)(l);
               SX_LOOP(iDir)  {
                  rotRho.set (0.);
                  // ref. 1, Eq. (48), sum over beta
                  SX_LOOP(jDir)
                     rotRho.plus_assign_ax (S(iDir,jDir), sym(jDir));

                  // ref. 1, Eq. (48), sum over kk'
                  rotRho = Dl ^ rotRho ^ Dl.transpose ();

                  int rotSite = symAtomIdx(ia, iSym) / nAtomsPerMolecule;
                  int offset, offset2;
                  if (nAtomsPerMolecule == 1)  {
                     offset = offset2 = 0;
                  } else {
                     offset  = (symAtomIdx(ia, iSym) & 1) * nm;
                     offset2 = (symAtomIdx(ja, iSym) & 1) * nm;
                  }
                  // write out
                  SxMatrix<Double> &rotD = (*Pij)(iDir,rotSite);
                  for (int m = 0; m < nm; ++m)
                     for (int m2 = 0; m2 < nm; ++m2)
                        rotD(offset + m, offset2 + m2) = rotRho(m,m2);
               }
            }
         }
      }
   }
}

SxAtomicStructure
SxHubbardMO::getForce (const SxBlockDensityMatrix &gradPij,
                       ssize_t iSpin) const
{
   SxAtomicStructure f2 = trafoForce.getNewStr ();
   f2.set (0., 0., 0.);
   int nm = 2 * l + 1;
   int nAtomsPerMolecule = (latShift.getSize () > 0) ? 2 : 1;
   SX_LOOP2(iDir, iSite)  {
      for (int ia = 0; ia < nAtomsPerMolecule; ++ia)  {
         int iTl = int(nAtomsPerMolecule * iSite + ia);
         for (int lm = 0; lm < nm; ++lm)  {
            int iAO = ia * nm + lm;
            for (int jAO = 0; jAO < nAtomsPerMolecule * nm; ++jAO)  {
               // ref. 1, Eq. (40)
               f2.ref(iTl)(iDir) += gradPij(iDir,iSite)(iAO,jAO)
                                  * hamProj(iSpin,iSite)(iAO, jAO);
               f2.ref(iTl)(iDir) += gradPij(iDir,iSite)(iAO,jAO)
                                  * hamProj(iSpin,iSite)(jAO, iAO);
            }
         }
      }
   }
   return f2;
}
