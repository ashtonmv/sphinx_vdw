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
#include <SxDipoleCorrZ.h>
#include <SxFFT1d.h>
#include <SxFileIO.h>

bool SxDipoleCorrZ::checkCell (const SxCell &cell)
{
   Coord z(0., 0., 1.);
   if (fabs(cell.basis(0) ^ z) > 1e-6
       || fabs(cell.basis(1) ^ z) > 1e-6
       || fabs(cell.basis(2).x(z).norm ()) > 1e-6)
   {
      cout << "Dipole correction implemented only for cells where\n"
              "the third axis is the z-axis, and the other two are\n"
              "orthogonal to it" << endl;
      return false;
   }
   return true;
}

SxDipoleCorrZ::SxDipoleCorrZ (const SxMeshR &rhoRTotal,
                              const SxMeshR &rhoHartreeR,
                              double chargeIn)
: fixedDipoleZ (false), cutZ (0.), dipole (0.), extraField (0.),
  charge(chargeIn)
{
   update (rhoRTotal, rhoHartreeR);
}

void SxDipoleCorrZ::update (const SxMeshR &rhoRTotal,
                            const SxMeshR &rhoHartreeR)
{
   // --- get cell and mesh dimensions and factors
   ssize_t meshSize = rhoRTotal.getSize ();
   const SxRBasis *rPtr = dynamic_cast<const SxRBasis *> 
                          (rhoRTotal.getBasisPtr ());
   SX_CHECK (rPtr);
   SxMesh3D mesh = rPtr->fft3d.mesh;
   const SxCell &cell = rPtr->cell;
   
   int z, nz = mesh(2);
   double dZ = cell(2,2) / nz;
   
   // --- check for correct cell
   if (! SxDipoleCorrZ::checkCell (cell)) {
      SX_QUIT;
   }
   
   // --- find density minimum (from total electron density)
   SxMeshR avgRho(mesh(2));
   avgRho.set (0.);
   SxMeshR::Iterator it = rhoRTotal.begin ();
   for (ssize_t ir = 0; ir < meshSize; ++ir, ++it)
      avgRho(mesh.getMeshVec (ir, SxMesh3D::Positive)(2)) += *it;
   avgRho /= (mesh(0) * mesh(1));

   SxVector<Double> zVals (mesh(2));
   for (int i = 0; i < avgRho.getSize (); i++)  zVals (i) = i * dZ;

   if (fixedDipoleZ) {
      dipolePos = cell.getMapped (Coord(0.,0.,dipolePos), SxCell::Positive)(0);
      rhoMinZ = int(dipolePos / dZ);
   } else  {
      rhoMinZ = getDipolePos(avgRho, dZ, 0.01 * avgRho.sum () * dZ);
      dipolePos = rhoMinZ * dZ;
      //fixedDipoleZ = true;
   }

   rhoMin = avgRho(rhoMinZ);

   // --- determine the absolute cutZ position
   {
      double c = fabs(cell(2,2));
      double cutZnow = dipolePos;
      if (cutZnow - cutZ > 0.8 * c) cutZnow -= c;
      if (cutZnow - cutZ < -0.8 * c) cutZnow += c;
      if (fabs(cutZnow - cutZ) > 0.3 * c && cutZ != 0.)  {
         cout << "WARNING: big shift in cut from " << cutZ
              << " to " << cutZnow << endl;
         cout << "         The total energy may jump!" << endl;
      }
      if (cutZnow == 0.) cutZnow = 1e-17;
      cutZ = cutZnow;
      cout << "cut=" << cutZ << endl;
   }

   // --- get dipole (from total charge density incl. nuclei)
   double oldDipole = dipole;
   double dipole_ = 0.;
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:dipole_)
#endif
   for (ssize_t ir = 0; ir < meshSize; ++ir)  {
      z = mesh.getMeshVec(ir, SxMesh3D::Positive)(2);
      if (z < rhoMinZ)  {
         dipole_ += z * rhoHartreeR(ir);
      } else {
         dipole_ += (z - nz) * rhoHartreeR(ir);
      }
   }
   dipole = dipole_ * dZ* rPtr->dOmega;
   dipole -= charge * dipolePos;

   if (SxLoopMPI::me()==0 && oldDipole!=0. && fabs(oldDipole - dipole) > 0.2)  {
      cout << "WARNING: BIG DIPOLE CHANGE: PRINT avgRhos !!!!!!!!" << endl;
      SxMeshR avgRhoHartree(mesh(2));
      avgRhoHartree.set (0.);
      for (ssize_t ir = 0; ir < meshSize; ++ir)  {
         z = mesh.getMeshVec(ir, SxMesh3D::Positive)(2);
         avgRhoHartree(z) += rhoHartreeR(ir);
      }
      try {
         SxFileIO avgRhoTotFile, avgRhoHartreeFile;
         avgRhoTotFile.open ("avgRhoTotal.dat","a");
         avgRhoHartreeFile.open ("avgRhoHartree.dat","a");
         for (int i = 0; i < zVals.getSize (); i++)  {
            SxString line = SxString(zVals(i)) + "\t" + SxString(avgRho(i)) + "\n";
            avgRhoTotFile.write (line);
            line = SxString(zVals(i)) + "\t" + SxString(avgRhoHartree(i)) + "\n";
            avgRhoHartreeFile.write (line);
         }
         SxString xmgBreak = "&\n";
         avgRhoTotFile.write (xmgBreak);
         avgRhoHartreeFile.write (xmgBreak);
      } catch (const SxException &e)  {
         e.print ();
         cout << "error ignored..." << endl;
      }
   }

   cout << SX_SEPARATOR;
   cout << "| Dipole correction" << endl;
   cout << "| dipole position @ z = " << (rhoMinZ * dZ)
         << ": density = " << rhoMin << endl;
   cout << "| dipole: " << dipole << endl;
   // Do we have a charged system ?
   if (fabs(charge) > 1e-10) {
      cout << "| q=" << charge << endl;
   }
   cout << SX_SEPARATOR;

}

int SxDipoleCorrZ::getDipolePos (const SxMeshR &avgRho, double dZ, double accRhoMax)
{

   unsigned int dimZ = int(avgRho.getSize ());
   int result;
   double accRho = avgRho.minval(&result);

   SxStack<int> vacRegionStack;
   vacRegionStack << result;

   int left = result - 1;
   int right = result + 1;
   while ((accRho < accRhoMax) && (vacRegionStack.getSize () < dimZ))  {
      double rhoLeft = avgRho((dimZ + left) % dimZ);
      double rhoRight = avgRho((dimZ + right) % dimZ);
      if (rhoLeft <= rhoRight)  {
         vacRegionStack << left;
         accRho += rhoLeft * dZ;
         left--;
      } else {
         vacRegionStack << right;
         accRho += rhoRight * dZ;
         right++;
      }
   }

   SxVector<Int> vacRegion (vacRegionStack);
   result = (dimZ + vacRegion.sum () / int (vacRegion.getSize ())) % dimZ;

   return result;

}

void SxDipoleCorrZ::correctPotential (SxMeshR *vHartreeRPtr) const
{
   SX_CHECK (vHartreeRPtr);
   const SxRBasis *rPtr = dynamic_cast<const SxRBasis *> 
                          (vHartreeRPtr->getBasisPtr ());
   SX_CHECK (rPtr);
   
   const SxCell &cell = rPtr->cell;
   const SxMesh3D mesh = rPtr->fft3d.mesh;
   int nz = mesh(2);
   
   // field multiplied with dZ
   double dZ = cell(2,2) / nz;
   double field = FOUR_PI * dipole / cell.volume * dZ;
   double bg = 0.;
   if (fabs(charge) > 1e-10)  {
      bg = -TWO_PI * charge / cell.volume * dZ * dZ;
      field += bg * (nz - 2 * rhoMinZ);
   }
   field += extraField * dZ;
   // make average of field to zero
   double dV = (rhoMinZ - nz + 0.5 * (nz - 1)) * field;
   if (fabs(charge) > 1e-10)  {
      dV = 0.;
      // average of bg potential
      dV -= bg * (rhoMinZ * (rhoMinZ - nz) + nz * nz / 3.);
      // average of field potential
      dV -= (rhoMinZ - 0.5 * nz) * field;
      // align such that zAlign (effective position of charge) is at 0
      dV -= TWO_PI * dipole * dipole/(cell.volume * charge);
      double fieldFromCharge = -FOUR_PI * charge/cell.volume * cell(2,2);
      cout << "field from charge " << fieldFromCharge << endl;
      cout << "zAlign=" << getZalign () << endl;
   }
   
   SxMeshR &vHartreeR = *vHartreeRPtr;
   // --- add saw-tooth potential
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
   for (ssize_t ir = 0; ir < vHartreeR.getSize (); ++ir)  {
      int z = mesh.getMeshVec(ir, SxMesh3D::Positive)(2);
      if (z < rhoMinZ)
         vHartreeR(ir) += z * field + dV + bg * z * z;
      else
         vHartreeR(ir) += (z - nz) * field + dV + bg * (z - nz) * (z - nz);
   }
}

void SxDipoleCorrZ::correctForces(SxAtomicStructure *forcePtr,
                                  const SxAtomicStructure &str,
                                  const SxSpeciesData &specData) const
{
   SX_CHECK (forcePtr);
   SxAtomicStructure &fHartree = *forcePtr;
   if (fabs(charge) < 1e-10)  {
      // --- neutral systems: field force is constant
      double field = FOUR_PI * dipole / fHartree.cell.volume;
      field += extraField;
      for (int is=0; is < fHartree.getNSpecies (); is++)  {
         double fField = field * specData.valenceCharge(is);
         for (int ia=0; ia < fHartree.getNAtoms(is); ia++)  {
            fHartree.ref(is,ia)(2) += fField;
         }
      }
   } else {
      // --- charged systems
      double field = FOUR_PI * dipole / fHartree.cell.volume;
      field += extraField;
      double bg = -FOUR_PI * charge / str.cell.volume;
      double c = str.cell(2,2);
      field += bg * (0.5 * c - dipolePos);
      SX_LOOP2(is,ia)  {
         double z = str.cell.getMapped (str(is,ia), SxCell::Positive)(2);
         double fField = field + bg * z;
         if (fabs(z) >= fabs(dipolePos))
            fField -= bg * c;
         fHartree.ref(is,ia)(2) += fField * specData.valenceCharge (is);
      }
   }
}

double SxDipoleCorrZ::getVRZero (const SxMeshR &vCorrected) const
{
   const SxRBasis &R = vCorrected.getBasis<SxRBasis> ();
   const SxMesh3D &mesh = R.fft3d.mesh;
   int z = rhoMinZ - 1;
   if (z < 0) z += R.fft3d.mesh(2);
   // --- find the average potential at rhoMinZ - 1
   double vAvg = 0.;
   for (int i = 0; i < R.fft3d.mesh(0); i++)  {
      for (int j = 0; j < R.fft3d.mesh(1); j++)  {
         vAvg += vCorrected(mesh.getMeshIdx(i,j,z, SxMesh3D::Positive));
      }
   }
   vAvg /= R.fft3d.mesh(0) * R.fft3d.mesh(1);

   // now determine the field (to the right of the slab = left of cut)
   double A = R.cell.volume / R.cell(2,2);
   double fieldRight = -FOUR_PI * charge / A
                     + extraField;
   double dZ = R.cell(2,2) / double (R.fft3d.mesh(2));
   return vAvg - (cutZ - dZ) * fieldRight;
}

double SxDipoleCorrZ::getVRZero (const PsiG &vCorrected) const
{
   const SxGBasis &G = vCorrected.getBasis<SxGBasis> ();
   SX_CHECK (G.structPtr);
   const SxCell &cell = G.structPtr->cell;
   const SxMesh3D &mesh = G.fft3d(0).mesh;
   double piZ = TWO_PI * (rhoMinZ - 1) / double(mesh(2));
   // --- find the average potential at rhoMinZ - 1
   double vAvg = 0.;
   for (int ig = 0; ig < G.n123(0).getSize (); ++ig)  {
      SxVector3<Int> idx = mesh.getMeshVec (G.n123(0)(ig), SxMesh3D::Positive);
      if (idx(0) == 0 && idx(1) == 0)  {
         vAvg += (vCorrected(ig) * exp(-I * piZ * double(idx(2)))).re;
         cout << vCorrected(ig) << ' ' << idx(2) << endl;
      }
   }
   vAvg /= sqrt(cell.volume);

   // now determine the field (to the right of the slab = left of cut)
   double c = cell(2,2);
   double A = cell.volume / c;
   double fieldRight = -FOUR_PI * charge / A
                     + extraField;
   return vAvg - cutZ * fieldRight;
}

double SxDipoleCorrZ::getElectrodeEnergy (const SxCell &cell) const
{
   double c = cell(2,2);
   double A = cell.volume / c;
   double QL = extraField /(-FOUR_PI) * A;
   //double QR = -(QL+charge);
   double zCR = dipole/charge;
   double VL = -FOUR_PI * dipole/A - cutZ * (FOUR_PI * charge / A);
   cout << "VL=" << VL * HA2EV << endl;
   double VC = (-FOUR_PI * charge / A + extraField) * (cutZ + zCR);
   cout << "VC=" << VC * HA2EV << endl;
   // note: right-hand potential defines the energy zero
   // note: central potential is part of the Hartree energy
   return 0.5 * (QL * VL);
}

