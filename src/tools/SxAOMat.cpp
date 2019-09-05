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
#include <SxPW.h>
#include <SxBinIO.h>
#include <SxProjector.h>
#include <SxAOBasis.h>
#include <SxKGrid.h>
#include <SxPseudoPot.h>
#include <SxFermi.h>
#include <SxGrid.h>
#include <SxHamSolver.h>
#include <SxPAWOverlap.h>

enum AoMatTimer { RtoK, kToR, Ham, Startup, Diag, CalcBS};
SX_REGISTER_TIMERS(AoMatTimer)
{
   regTimer (RtoK, "R->k");
   regTimer (kToR, "k->R");
   regTimer (Ham, "Ham");
   regTimer (Startup, "startup");
   regTimer (Diag, "diag");
   regTimer (CalcBS, "calcBS");
}
typedef Double SMatType;
typedef SxNArray<double, 3> SxClebschTable;

SxArray<int>  createTBInteraction (const SxAOBasis &aoBasis, int iOrb)
{
   SxArray<int> result;
   int is = aoBasis.orbitalMap(iOrb).is;
   int ioAtom = aoBasis.orbitalMap(iOrb).io;
   int m = aoBasis.refOrbMap(is)(ioAtom).m;
   for (int jOrb = 0; jOrb < aoBasis.getNOrb(); jOrb++)  {
      int js = aoBasis.orbitalMap(jOrb).is;
      int joAtom = aoBasis.orbitalMap(jOrb).io;
      int jm = aoBasis.refOrbMap(js)(joAtom).m;
      if ((m == 0) && (jm ==  0)) result.append(jOrb);
      if ((m != 0) && (m != -jm)) result.append(jOrb);
   }
   return result;
} 

template <class T>
SxDiracVec<T> nice (const SxDiracVec<T> &in)
{
   SxDiracVec<Complex16> res;
   res.copy (in);
   for (int i = 0; i < res.getSize (); ++i)  {
      if (fabs(res(i).re) < 1e-12) res(i).re = 0.;
      if (fabs(res(i).im) < 1e-12) res(i).im = 0.;
   }
   return res;
}


/** \brief Rotate and remove k-phase

  @note This attaches an i-phase to each orbital with odd l
  to compensate for the i^l phase convention of atomic orbitals.
  */
SxDiracMat<Complex16> 
rotateAoMatrix (const SxDiracMat<Complex16> &S,
                const SymMat &sym,
                const SxAOBasis &aoBasis,
                const SxAtomicStructure &structure,
                const SxGrid &grid,
                const SxArray<SxMatrix<Double> > &ylmRot,
                const Coord &k)
{
   Coord rotK = sym ^ k;
   int nOrb = aoBasis.getNOrb ();

   SxDiracMat<Complex16> rotS(nOrb, nOrb);

   // --- rotate the rows
   for (int iOrb = 0; iOrb < nOrb; /* done inside */)  {
      int is = aoBasis.orbitalMap(iOrb).is;
      int ia = aoBasis.orbitalMap(iOrb).ia;
      int nOrbLocal = int(aoBasis.refOrbMap(is).getSize ());
      int ioAtom = aoBasis.orbitalMap(iOrb).io;
      int l = aoBasis.refOrbMap(is)(ioAtom).l;
      int nl = 2 * l + 1;
      SX_CHECK (aoBasis.refOrbMap(is)(ioAtom).m == -l,
                aoBasis.refOrbMap(is)(ioAtom).m, l);

      // --- collect the rows
      SxDiracMat<Complex16> cols, rowsT;
      cols = S(SxIdx(iOrb * nOrb, (iOrb + nl) * nOrb - 1));
      cols.reshape (nOrb, nl);
      rowsT = cols.conj ();

      // rotate the rows
      if (l > 0)
         rowsT = rowsT ^ toVector(ylmRot(l)).transpose ();
      if (l & 1)
         rowsT *= -I; // remove imaginary phase

      // --- find equivalent atom
      Coord rotPos = sym ^ structure.getAtom(is,ia);
      int iaRot = structure.find (rotPos, grid)
                - structure.atomInfo->offset(is);
      int iOrbRot = iOrb + nOrbLocal * (iaRot - ia);
      
      // --- phase shift
      double phaseShift = 0.;
      phaseShift -= rotK ^ structure(is, iaRot);
      phaseShift += k ^ structure.getAtom(is, ia);
      //cout << exp(I * phaseShift) << endl;
      rowsT *= exp(-I * phaseShift);

      // place rotated rows into cols at rotated atom 
      rotS(SxIdx(iOrbRot * nOrb, (iOrbRot + nl) * nOrb - 1)) << rowsT;

      iOrb += nl;
   }
   // exchange rows and columns
   SxDiracMat<Complex16> rotS2 = rotS.transpose ();
   // --- rotate the columns
   for (int iOrb = 0; iOrb < nOrb; /* done inside */)  {
      int is = aoBasis.orbitalMap(iOrb).is;
      int ia = aoBasis.orbitalMap(iOrb).ia;
      int nOrbLocal = int(aoBasis.refOrbMap(is).getSize ());
      int ioAtom = aoBasis.orbitalMap(iOrb).io;
      int l = aoBasis.refOrbMap(is)(ioAtom).l;
      int nl = 2 * l + 1;
      SX_CHECK (aoBasis.refOrbMap(is)(ioAtom).m == -l,
                aoBasis.refOrbMap(is)(ioAtom).m, l);

      // --- collect the columns
      SxDiracMat<Complex16> cols, colsRot;
      cols = rotS2(SxIdx(iOrb * nOrb, (iOrb + nl) * nOrb - 1));
      cols.reshape (nOrb, nl);

      // rotate the columns
      if (l > 0)
         colsRot = cols ^ toVector(ylmRot(l)).transpose ();
      else
         colsRot = cols * 1.;
      if (l & 1)
         colsRot *= I; // remove imaginary phase

      // --- find equivalent atom
      Coord rotPos = sym ^ structure.getAtom(is,ia);
      int iaRot = structure.find (rotPos, grid)
                - structure.atomInfo->offset(is);
      int iOrbRot = iOrb + nOrbLocal * (iaRot - ia);
      
      // --- phase shift
      double phaseShift = 0.;
      phaseShift -= rotK ^ structure(is, iaRot);
      phaseShift += k ^ structure.getAtom(is, ia);
      //cout << exp(I * phaseShift) << endl;
      colsRot *= exp(I * phaseShift);

      // place rotated columns at rotated atom
      rotS(SxIdx(iOrbRot * nOrb, (iOrbRot + nl) * nOrb - 1)) << colsRot;

      iOrb += nl;
   }
   return rotS;
}


SxDiracMat<Complex16>
aoMatrixRtoK (const SxArray<SxDiracMat<SMatType> > AinR,
              const SxAtomicStructure &structure,
              const SxMesh3D &meshR,
              const SxAOBasis &aoBasis,
              const Coord &k
             )
{
   SX_CLOCK (RtoK);
   ssize_t nR = meshR.getSize ();
   int nOrb = aoBasis.getNOrb ();
   SX_CHECK (nR == AinR.getSize (), nR, AinR.getSize ());
   SxCell superCell(structure.cell.basis(0) * meshR(0),
                    structure.cell.basis(1) * meshR(1),   
                    structure.cell.basis(2) * meshR(2));
   superCell = superCell.getRegularCell ();
   SxMesh3D tm(3,3,3);

   SxDiracMat<Complex16> AinK(nOrb, nOrb);
   AinK.set (0.);
   for (int iR = 0; iR < nR; ++iR)  {
      SxVector3<Int> relR = meshR.getMeshVec (iR, SxMesh3D::Origin);
      Coord vecR = structure.cell.relToCar (relR);
      const SxDiracMat<SMatType> &AR = AinR(iR);

      int offsetJ = 0;
      for (int js = 0; js < structure.getNSpecies (); ++js)  {
         int noJ = int(aoBasis.refOrbMap(js).getSize ());
         for (int ja = 0; ja < structure.getNAtoms (js); ++ja)  {

            int offsetI = 0;
            for (int is = 0; is < structure.getNSpecies (); ++is)  {
               int noI = int(aoBasis.refOrbMap(is).getSize ());
               for (int ia = 0; ia < structure.getNAtoms (is); ++ia)  {

                  //Coord d = structure(is,ia) - structure(js,ja) - vecR;
                  Coord d = structure(is,ia) - structure(js,ja) - vecR;
                  superCell.map (&d, SxCell::WignerSeitz);
                  SxComplex16 phase = exp(-I * (k ^ d));

                  // --- handle d's at boundary of WS cell
                  double d2 = d.normSqr ();
                  int nb = 1;           
                  for (int itm = 1; itm < 27; ++itm)  {
                     SxVector3<Int> tRel = tm.getMeshVec(itm, SxMesh3D::Origin);
                     Coord dShift = d + superCell.relToCar (tRel);
                     if (fabs(dShift.normSqr () - d2) < 1e-6 * d2)  {
                        nb++;
                        phase += exp(-I * (k ^ dShift));
                     }
                  }
                  phase /= nb;

                  phase *= exp (I*(k^(structure(is,ia) - structure(js,ja))));

                  for (int jo = offsetJ; jo < offsetJ + noJ; ++jo)  {
                     for (int io = offsetI; io < offsetI + noI; ++io)  {
                        AinK(io, jo) += phase * AR(io, jo);
                     }
                  }
                  offsetI += noI;
               }
            }
            offsetJ += noJ;
         }
      }
   }
   return AinK;
}

// --- Calculate matrix to rotate (dist) to (0,0,1)
SxMatrix3<Double> getRotationMatrix (const Coord &dist)
{
   // REF: wikipedia, rotation matrix
   // alpha is angle (x, ProjXY(dist)(x))
   // beta is angle (z, dist(z))
   //
   //      | cos(beta) 0 -sin(beta)|   | cos(alpha)   sin(alpha) 0 |
   //rot = |    0      1     0     | * | -sin(alpha)  cos(alpha) 0 |
   //      | sin(beta) 0  cos(beta)|   |     0            0      1 | 
   
   SxMatrix3<Double> result;
   Coord d = 1.0 * dist;
   if (d.norm() < 1e-8 || (fabs(d(0)) + fabs(d(1)) < 1e-8))  { 
      // identical atom or d = alpha * z -> no Rotation
      result(0,0) = 1.0, result(0,1) = 0.0, result(0,2) = 0.0;
      result(1,0) = 0.0, result(1,1) = 1.0, result(1,2) = 0.0;
      result(2,0) = 0.0, result(2,1) = 0.0, result(2,2) = 1.0;
      if(d(2) < 0) result(2,2) = -1.;
   } else {
      d.normalize();
      Coord dProjXY;
      dProjXY(0) = d(0), dProjXY(1) = d(1), dProjXY(2) = 0.0;
      dProjXY.normalize();
      // Rotate all to 0,0,1
      double sign = dProjXY(1) < 0 ? -1.0 : 1.0;
      double ca = dProjXY(0), sa = sign * sqrt(1-ca*ca);
      double cb = d(2), sb = sqrt(1-cb*cb);

      SxMatrix3<Double> rotToX, rotToZ;
      rotToX(0,0) =  ca;  rotToX(0,1) =  sa; rotToX(0,2) = 0.0;
      rotToX(1,0) = -sa;  rotToX(1,1) =  ca; rotToX(1,2) = 0.0;
      rotToX(2,0) = 0.0;  rotToX(2,1) = 0.0; rotToX(2,2) = 1.0;
     
      rotToZ(0,0) =  cb;  rotToZ(0,1) = 0.0; rotToZ(0,2) = -sb;
      rotToZ(1,0) = 0.0;  rotToZ(1,1) = 1.0; rotToZ(1,2) = 0.0;
      rotToZ(2,0) =  sb;  rotToZ(2,1) = 0.0; rotToZ(2,2) =  cb;

      result = rotToZ ^ rotToX;
   }

   // --- For DEBUG
   if (d.norm () > 1e-8)  {
      d.normalize ();
      if ( ((result ^ d)- Coord(0.,0.,1.)).norm () > 1e-6)  {
         cout << "WRONG ROTATION!" << endl;
         cout << "Deviationnorm is " << ((result ^ d)- Coord(0.,0.,1.)).norm () << endl;
         cout << "Determinant of rotation matrix is: " << result.determinant() << endl;
         d.print(); cout << " goes to "; (result ^ d).print (); cout << endl;
         SX_EXIT;
      }
   }

   return result;
}

void switchOffInteraction (const SxArray<SxArray<int> > &blockInteraction, 
                           SxConstPtr<SxAOBasis> aoPtr, 
                           SxDiracMat<SMatType> &SinR, 
                           SxDiracMat<SMatType> &HinR,
                           const int iSpecies, const int jSpecies,
                           const int offsetI, const int offsetJ)
{
   int noI = int(aoPtr->refOrbMap(iSpecies).getSize ());
   int noJ = int(aoPtr->refOrbMap(jSpecies).getSize ());
   for (int iBlock = 0; iBlock < blockInteraction.getSize (); iBlock++)  {
      if ( (iSpecies == blockInteraction(iBlock)(0) && jSpecies == blockInteraction(iBlock)(2))
        || (iSpecies == blockInteraction(iBlock)(2) && jSpecies == blockInteraction(iBlock)(0)))  {
         for (int io = offsetI; io < offsetI + noI;)  {
            int iL = aoPtr->refOrbMap(iSpecies)(io-offsetI).l;
            int nOrbsI = 2*iL+1;
            for (int jo = offsetJ; jo < offsetJ + noJ;)  {
               int jL = aoPtr->refOrbMap(jSpecies)(jo-offsetJ).l;
               int nOrbsJ = 2*jL+1;
                                 
               if ( (iSpecies == blockInteraction(iBlock)(0) 
                        && iL == blockInteraction(iBlock)(1) 
                        && jL == blockInteraction(iBlock)(3))
                  ||(iSpecies == blockInteraction(iBlock)(2) 
                        && iL == blockInteraction(iBlock)(3) 
                        && jL == blockInteraction(iBlock)(1)))  {

                  SxDiracMat<SMatType> block;
                  // Switchoff Overlap S 
                  block = SinR.getBlock(io, jo, nOrbsI, nOrbsJ);
                  block.set(0.0);
                  SinR.setBlock(block, io, jo);
                  // Switchoff Hamiltoncontribution H
                  block = HinR.getBlock(io, jo, nOrbsI, nOrbsJ);
                  block.set(0.0);
                  HinR.setBlock(block, io, jo);

               }

               jo += nOrbsJ;
            } // jo
            io += nOrbsI;
         } // io
      } // if iSpecies && jSpecies
   } // blockInteraction
}

void writeMatElems (SxBinIO &io, SxArray<SxDiracMat<SMatType> > &SinR, SxArray<SxDiracMat<SMatType> > &HinR, 
                    SxAtomicStructure &structure, SxArray<Coord> &vecR, SxAOBasis &aoBasis, 
                    SxMesh3D &meshR, int iSpecies)
{
   SxCell superCell(structure.cell.basis(0) * meshR(0),
                    structure.cell.basis(1) * meshR(1),   
                    structure.cell.basis(2) * meshR(2));

   //calculate number of neighbors and number of all orbitals
   int nR = int(vecR.getSize());
   int nNeighbors = 0;
   int nOrbs = int(aoBasis.refOrbMap(iSpecies).getSize ());

   // Each atom of species has nR * nTlAtoms neighbors
   // -> nNeigbours = nAtoms(iSpecies) * nR * nTlAtoms
   nNeighbors = structure.getNAtoms(iSpecies) * nR * structure.getNAtoms ();

   // --- Calculate total number of orbitals to exchange with
   SxVector<Int> nNeighborOrbs (nNeighbors);
   int iNeighbor = 0;
   for (int ia = 0; ia < structure.getNAtoms (iSpecies); ia++)
      for (int iR = 0; iR < nR; iR++) 
         for (int js = 0; js < structure.getNSpecies (); js++)  
            for (int ja = 0; ja < structure.getNAtoms (js); ja++, iNeighbor++)  
               nNeighborOrbs(iNeighbor) = int(aoBasis.refOrbMap(js).getSize ());
   int nNOrbitals = nNeighborOrbs.sum();

   SxMatrix<Double> dist (nNeighbors,3);
   SxVector<Int> neighborSpecies (nNeighbors);
   SxDiracMat<SMatType> SElems (nNOrbitals, nOrbs), HElems(nNOrbitals, nOrbs);
   SxMatrix<Double> RElems(3 * nNeighbors,3);

   ssize_t offsetI = 0;
   for (int is = 0; is < iSpecies; is++)  {
      offsetI += structure.getNAtoms(is) 
                 * aoBasis.refOrbMap(is).getSize ();
   }
   
   ssize_t noI = aoBasis.refOrbMap(iSpecies).getSize ();
   iNeighbor = 0;
   int writeOffset = 0;
   for (int ia = 0; ia < structure.getNAtoms (iSpecies); ia++, offsetI += noI)  {
      for (int iR = 0; iR < vecR.getSize(); iR++)  {
         int offsetJ = 0;
         for (int js = 0; js < structure.getNSpecies (); js++)  {
            int noJ = int(aoBasis.refOrbMap(js).getSize ());
            for (int ja = 0; ja < structure.getNAtoms (js); 
                     ja++, offsetJ += noJ, writeOffset += noJ, iNeighbor++)  {
               Coord d = structure(js,ja) + vecR(iR) - structure(iSpecies,ia);
               superCell.map (&d, SxCell::WignerSeitz);
               dist(iNeighbor,0) = d(0);
               dist(iNeighbor,1) = d(1);
               dist(iNeighbor,2) = d(2);
               neighborSpecies(iNeighbor) = js;

               SxMatrix<Double> RotMat = getRotationMatrix(d);
               RElems.setBlock(RotMat, 3 * iNeighbor, 0);

               SxDiracMat<SMatType> block = SinR(iR).getBlock(offsetI, offsetJ, noI, noJ);
               SElems.setBlock(block.transpose (), writeOffset, 0);

               block = HinR(iR).getBlock(offsetI, offsetJ, noI, noJ);
               HElems.setBlock(block.transpose (), writeOffset, 0);
            }
         }
      }
   }

   try {
      // write variables
      io.addDimension("nOrbs", nOrbs);
      io.addDimension("nNeighbors", nNeighbors);
      io.addDimension("nCols", nNOrbitals);
      io.addDimension("RCols", 3 * nNeighbors);
      io.write("neighborSpecies", neighborSpecies, "nNeighbors");
      io.write("nNeighborOrbs", nNeighborOrbs, "nNeighbors");

      io.write("distances", dist, "nNeighbors", "coordDim");
      io.write("R-Matrices", RElems, "RCols", "coordDim");
      io.write("S-Matrices", SElems, "nCols", "nOrbs");
      io.write("H-Matrices", HElems, "nCols", "nOrbs");


   } catch (SxException e)   {
      e.print ();
      SX_EXIT;
   }
}

int main (int argc, char **argv)
{

   SX_START_TIMER (Startup);
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   SxCLI cli (argc,argv);
   cli.authors = "C. Freysoldt and B. Lange";

   cli.newGroup("Scissor operator");

   bool scissor = cli.option ("--scissor", "apply scissor shift").toBool();

   double scissorU = cli.option ("-u", "double value", "energetic shift for band eigenenergies")
                     .toDouble();

   double scissorRangeLower = cli.option ("--rLow", "double value", "lowest energy value to apply scissor shift")
                              .toDouble();
   double scissorRangeUpper = cli.option ("--rUp", "double value", "highest energy value to apply scissor shift")
                              .toDouble();
   double dCutoff           = cli.option ("--cut", "double value", "maximal distance to apply scissor shift correction in real space")
                              .toDouble();

   cli.newGroup("General group");

   double neighborCut   = cli.option ("--nCut", "double value", "maximal Neighbor distance")
                              .toDouble(-1.0);
   
   SxString inFile = cli.option ("-i|--input", "file", "S/PHI/nX input file")
                     .toString ("input.sx");

   SxString kFile = cli.option ("-k|--kFile", "file", "S/PHI/nX k-points file")
                     .toString ("");

   SxString aoFile = cli.option ("--basis", "file", "netcdf ao file")
                     .toString ("");
   
   SxString interactionFile = cli.option ("--interaction", "file", "interactionfile").toString ("");

   SxString rhoFile = 
      cli.option ("-r|--rho","file","density file")
      .toString ("rho.sxb");

   SxMesh3D meshR
      = SxVector3<Int> (cli.option ("-m", "vector", "R-mesh").toIntList3 ());

   bool printR = cli.option ("--printR","print R-space matrices").toBool ();
   bool printS = cli.option ("--printS","print S eigenvalues").toBool ();
   bool write = cli.option ("--write", "write Rotation matrices").toBool();

   cli.finalize ();


   SxParser parser;
   SxParser::Table table = parser.read (inFile);

   // setup HamSolver
   SxHamSolver hamSolver(rhoFile, table);
   
   // setup GkBasis
   SxKPoints kPoints (hamSolver.structure.cell,
                      SxList<SxVector3<TPrecG> > () << Coord (0,0,0),
                      SxList<double> () << 1.0, meshR);
   double eCut = SxGBasis::getECut(&*table);
   SxPtr<SxGkBasis> gkBasisPtr = SxPtr<SxGkBasis>::create(kPoints, hamSolver.structure, 
                                                          hamSolver.mesh, eCut, true); 
   SxGkBasis &gkBasis = *gkBasisPtr;
   
   // setup Hamiltonian
   SxAtomicStructure &structure = hamSolver.structure;
   hamSolver.setupHam(rhoFile, table, gkBasisPtr);
   SxHamiltonian &ham = *hamSolver.hamPtr;

   // done, now start calculation
   const SxSymGroup &syms = *structure.cell.symGroupPtr;
   int nSym = syms.getNSymmorphic ();

   bool PAW = false;
   if (table->containsGroup("PAWHamiltonian")) PAW = true;

   // Overlap
   SxPtr<SxOverlapBase> SPtr;
   SxArray<SxArray<SxDiracVec<Double> > > pseudoPsi;
   SxConstPtr<SxRadBasis> radPtr;
   int lmax = 0;
   if (PAW)  {
      SxPtr<SxPAWPot> pawPotPtr = SxPtr<SxPAWPot> (hamSolver.potPtr);
      SxPtr<SxPartialWaveBasis> pBasis = SxPtr<SxPartialWaveBasis>::create (pawPotPtr, structure);
      pBasis->createProjBasis (gkBasis);
      SPtr = SxPtr<SxPAWOverlap>::create (pBasis, pawPotPtr);
      pseudoPsi = pawPotPtr->getPhiPS ();
      radPtr = radPtr.create(pawPotPtr->rad, pawPotPtr->logDr);
      lmax = pawPotPtr->getLMax();
   } else {
      SxPtr<SxPseudoPot> pseudoPotPtr = SxPtr<SxPseudoPot> (hamSolver.potPtr);
      SPtr = SxPtr<SxPWOverlap>::create ();
      pseudoPsi = pseudoPotPtr->getPseudoPsi ();
      radPtr = radPtr.create(pseudoPotPtr->rad, pseudoPotPtr->logDr);
      lmax = pseudoPotPtr->getLMax();
   }

   cout << "Setting up ao basis" << endl;
   SxPtr<SxAOBasis> aoPtr;

   if (aoFile.getSize () == 0)  {
      aoPtr = aoPtr.create (gkBasis, *radPtr, pseudoPsi);
   } else {
      try {
         SxBinIO io (aoFile, SxBinIO::BINARY_READ_ONLY);
         SxAtomicOrbitals orbs(io);
         radPtr = orbs.getRadBasisPtr ();
         for (int is = 0; is < orbs.getNSpecies (); ++is)  {
            for (int iot = 0; iot < orbs.getNOrbTypes (is); ++iot)  {
               cout << "Orbital is = " << is << ", iot = " << iot 
                    << ": l = " << orbs.getL (is, iot) << endl;
               // --- orthonormalize l-channels
               if (!PAW) {
                  for (int jot = 0; jot < iot; ++jot)  {
                     if (orbs.getL (is, iot) == orbs.getL(is, jot))  {
                        cout <<  tr(orbs(is, iot) * orbs(is, jot)) << endl;
                        orbs(is, iot) -= orbs(is, jot)
                           * tr(orbs(is, iot) * orbs(is, jot));
                     }
                  }
                  orbs(is, iot) /= sqrt (tr (orbs(is, iot).sqr ()));
               }
            }
         }
         lmax = orbs.getLMax ();
         aoPtr = aoPtr.create (gkBasis, orbs);
      } catch (SxException e) {
         e.print ();
         SX_EXIT;
      }
   }

   SxArray<SxArray<int> > blockInteraction; 
   if (interactionFile.getSize () != 0)  {
      SxParser blockerParser;
      SxParser::Table blockerTable = parser.read (interactionFile);
      try   {
         // read Interactionblocker
         const SxSymbolTable *blocker;
         blocker = blockerTable->getGroup ("blockInteraction");
         int nBlockers = blocker->getNItems ("blockInteraction");
         blockInteraction.resize(nBlockers);

         int iBlocker = 0;
         for (blocker = blockerTable->getGroup ("blockInteraction");
              blocker != NULL;
              blocker = blocker->nextSibling ("blockInteraction"),iBlocker++)  {
            blockInteraction(iBlocker).resize(4);
            blockInteraction(iBlocker)(0) = blocker->get("is1")->toInt ();
            blockInteraction(iBlocker)(1) = blocker->get("l1")->toInt ();
            blockInteraction(iBlocker)(2) = blocker->get("is2")->toInt ();
            blockInteraction(iBlocker)(3) = blocker->get("l2")->toInt ();
         }
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
   }

   SxAOBasis &aoBasis = *aoPtr;
   aoBasis.setOverlapCaching (SxAOBasis::Recompute);
   //aoBasis.setInvOverlapCaching (SxAOBasis::CacheCurrentK);

   SxYlmRotGroup ylmRot 
      = SxYlm::computeYlmRotMatrices (syms.getSymmorphic (), lmax);

   SxAOBasis::TPsi psiAO, psi;
   int nk = gkBasis.getNk ();
   int nOrb = aoBasis.getNOrb ();
   ssize_t nR = meshR.getSize ();

   // Talk a little
   sxprintf("Basis has a total number of %i orbitals.\n",nOrb);
   sxprintf("Your choosen R mesh forms %i unit cells.\n",int(nR));

   SxArray<Coord> vecR(nR);
   SxArray<SxDiracMat<SMatType> > SinR(nR);
   SxArray<SxDiracMat<SMatType> > HinR(nR);
   SxArray<SxDiracMat<SMatType> > UinR(nR);
   SxClebschTable complexClebschTable;
   SxArray<SxMatrix<Double> > Dmm, Dmm2;
   if (write) complexClebschTable = SxYlm::getClebschGordan (lmax, lmax, lmax, SxYlm::ComplexYlm);

   for (int iR = 0; iR < nR; ++iR)  {
      SxVector3<Int> relR = meshR.getMeshVec (iR, SxMesh3D::Origin);
      vecR(iR) = structure.cell.relToCar (relR);
      SinR(iR).reformat (nOrb, nOrb);
      SinR(iR).set (0.);
      HinR(iR).reformat (nOrb, nOrb);
      UinR(iR).reformat (nOrb, nOrb);
   }
   SX_STOP_TIMER (Startup);

   SxVector<Int> symAtomIdx(nSym);
   SxGrid grid (structure, 10);
   int nSpin = hamSolver.nSpin;
   for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
   for (int ik = 0; ik < nk; ++ik)  {
      cout << "Computing S...(" << ik+1 << " of " << nk << ")" << endl;
      for (int iR = 0; iR < nR; ++iR)  {
         HinR(iR).set (0.);
         UinR(iR).set (0.);
      }

      SxDiracMat<Complex16> Sk(nOrb, nOrb);
      SxDiracMat<Complex16> Hk(nOrb, nOrb);
      SxDiracMat<Complex16> U(nOrb, nOrb);
      U.set(0.0);
      {
         SX_CLOCK (Ham);
         for (int io = 0; io < nOrb; ++io)  {
            PsiG nu = aoBasis.getAOinG (ik, io);
            nu.handle->auxData.i     = io;
            nu.handle->auxData.iSpin = iSpin;
            nu.handle->auxData.ik    = ik;
            PsiG Hnu = ham | nu;
            Hk.colRef (io) << (aoBasis | Hnu);
            if (iSpin == 0 || scissor)  {
               PsiG Snu = SPtr->apply(nu);
               Sk.colRef (io) << (aoBasis | Snu);
            }
         }
      }
      
      if (printS)  { 
         if (iSpin == 0)  {
            cout << "S: " << Sk.eigenvalues ().real () << endl;
            cout << "S=" << nice(Sk) << endl;
         }
         cout << "H: " << Hk.eigenvalues ().real () << endl;
         cout << "H=" << nice(Hk) << endl;
      }

      // perform Scissorshift
      if (scissor)  {
         cout << "Perform scissorshift with U = " << scissorU 
              << " eV in Range " << scissorRangeLower 
              << " < epsilon < " << scissorRangeUpper 
              << "." << endl;
         SxDiracMat<Complex16> L;
         L = Sk.choleskyDecomposition (UpperRight).inverse ();
         SxDiracMat<Complex16>::Eigensystem eig = (SxDiracMat<Complex16>(L.adjoint () ^ Hk ^ L)).eigensystem();
         cout << "E: " << HA2EV * eig.vals.real() << endl;
         for(int iRow = 0; iRow < nOrb; iRow++)  {
            if ((eig.vals(iRow).re * HA2EV >= scissorRangeLower) && 
                (eig.vals(iRow).re * HA2EV<= scissorRangeUpper)) 
               U(iRow, iRow) = scissorU/ HA2EV;
         }
         // Rotate in space
         U = (eig.vecs ^ U ^ eig.vecs.adjoint ());
         // Cholesky Rotation
         U = (Sk ^ L ^ U ^ L.adjoint() ^ Sk);
         // apply Scissor Shift
         Hk += U;
      }

      for (int iSym = 0; iSym < nSym; ++iSym)  {
         SymMat sym = syms.getSymmorphic(iSym);
         Coord rotK = sym ^ gkBasis.getK (ik);

         // --- k->R transform of S
         if (iSpin == 0)  {
            SxDiracMat<Complex16> rotS;
            rotS = rotateAoMatrix (Sk, sym, aoBasis, structure, grid,
                                   ylmRot(iSym), gkBasis.getK (ik));

            SX_START_TIMER (kToR);
            for (int iR = 0; iR < nR; ++iR)  {
               SxComplex16 phase = exp(-I * (rotK ^ vecR(iR)));
               double w = gkBasis.weights(ik) / nSym;
               SinR(iR).plus_assign_ax (w, (phase * rotS).real ());
               //SinR(iR).plus_assign_ax (w, (phase * rotS));
            }
            SX_STOP_TIMER (kToR);
         }

         // --- k->R transform of H and U
         SxDiracMat<Complex16> rotH;
         SxDiracMat<Complex16> rotU;
         rotH = rotateAoMatrix (Hk, sym, aoBasis, structure, grid,
                                ylmRot(iSym), gkBasis.getK (ik));
         rotU = rotateAoMatrix (U,  sym, aoBasis, structure, grid,
                                ylmRot(iSym), gkBasis.getK (ik));
         SX_START_TIMER (kToR);
         for (int iR = 0; iR < nR; ++iR)  {
            SxComplex16 phase = exp(-I * (rotK ^ vecR(iR)));
            double w = gkBasis.weights(ik) / nSym;
            HinR(iR).plus_assign_ax (w, (phase * rotH).real ());
            UinR(iR).plus_assign_ax (w, (phase * rotU).real ());
         }
         SX_STOP_TIMER (kToR);
      }
   }

   if (printR || scissor || write || blockInteraction.getSize () > 0)  {
      SxCell superCell(structure.cell.basis(0) * meshR(0),
                       structure.cell.basis(1) * meshR(1),   
                       structure.cell.basis(2) * meshR(2));
      int offsetI = 0;
      for (int is = 0; is < structure.getNSpecies (); ++is)  {
         int noI = int(aoBasis.refOrbMap(is).getSize ());
         for (int ia = 0; ia < structure.getNAtoms (is); ++ia)  {
            for (int iR = 0; iR < nR; ++iR)  {
               int offsetJ = 0;
               for (int js = 0; js < structure.getNSpecies (); ++js)  {
                  int noJ = int(aoBasis.refOrbMap(js).getSize ());
                  for (int ja = 0; ja < structure.getNAtoms (js); ++ja)  {
                     Coord d = structure(js,ja) + vecR(iR) - structure(is,ia);
                     superCell.map (&d, SxCell::WignerSeitz);
                     cout << SX_SEPARATOR;
                     cout << structure(is, ia) << "<=>"
                        << (structure(is,ia) + d) << "; d=" << d.norm ()
                        << endl << endl;
                     cout << "coords = ";
                     SxVector<Double>(d).print();

                     // Switchoff Interaction
                     if ( (d.norm () > 1e-6) && (blockInteraction.getSize () > 0) )
                        switchOffInteraction (blockInteraction, aoPtr, SinR(iR), HinR(iR),
                                              is, js, offsetI, offsetJ);

                     if ( (d.norm() > neighborCut) && (neighborCut > 0.0) )  {
                        SxDiracMat<SMatType> block;
                        block = SinR(iR).getBlock(offsetI, offsetJ, noI, noJ);
                        block.set(0.0);
                        SinR(iR).setBlock(block, offsetI, offsetJ);
                        block = HinR(iR).getBlock(offsetI, offsetJ, noI, noJ);
                        block.set(0.0);
                        HinR(iR).setBlock(block, offsetI, offsetJ);
                     }

                     if (scissor && d.norm() < dCutoff) {
                        cout << "Apply Correction!" << endl;
                        SxDiracMat<SMatType> block = 
                           HinR(iR).getBlock(offsetI, offsetJ, noI, noJ);
                        block += UinR(iR).getBlock(offsetI,offsetJ,noI,noJ);
                        HinR(iR).setBlock(block,offsetI,offsetJ);
                     }
                     
                     // --- Matrix rotation in bonding direction
                     if (write)  {
                        SxMatrix3<Double> rotMat = getRotationMatrix(d);
                        cout << "R = [";
                        sxprintf("[%.6f, %.6f, %.6f],\n",rotMat(0,0),rotMat(0,1),rotMat(0,2));
                        sxprintf("[%.6f, %.6f, %.6f],\n",rotMat(1,0),rotMat(1,1),rotMat(1,2));
                        sxprintf("[%.6f, %.6f, %.6f]];\n",rotMat(2,0),rotMat(2,1),rotMat(2,2));
                        Dmm = SxYlm::computeYlmRotMatrices (lmax, rotMat,
                                                            complexClebschTable);
                        
                        cout << "Before rotation" << endl;
                        cout << "S = [";
                        SxDiracMat<SMatType> block = SinR(iR).getBlock(offsetI, offsetJ, noI, noJ);
                        for (int iRow = 0; iRow < block.nRows(); iRow++)  {
                           cout << "[";
                           for (int iCol = 0; iCol < block.nCols(); iCol++)  {
                              sxprintf("%.8f",block(iRow,iCol));
                              if (iCol < block.nCols() - 1 ) sxprintf(", ");
                              else sxprintf("]");
                           }
                           if (iRow < block.nRows() - 1 ) sxprintf(",\n");
                              else sxprintf("];\n");
                        }
                        cout << "H = [";
                        block = HinR(iR).getBlock(offsetI, offsetJ, noI, noJ);
                        for (int iRow = 0; iRow < block.nRows(); iRow++)  {
                           cout << "[";
                           for (int iCol = 0; iCol < block.nCols(); iCol++)  {
                              sxprintf("%.8f",block(iRow,iCol));
                              if (iCol < block.nCols() - 1 ) sxprintf(", ");
                              else sxprintf("]");
                           }
                           if (iRow < block.nRows() - 1 ) sxprintf(",\n");
                              else sxprintf("];\n");
                        }

                        //Perform rotation
                        for (int io = offsetI; io < offsetI + noI;)  {
                           int iL = aoPtr->refOrbMap(is)(io-offsetI).l;
                           int nOrbsI = 2*iL+1;
                           SxDiracMat<Double> rotLeft = toVector(Dmm(iL));
                           for (int jo = offsetJ; jo < offsetJ + noJ;)  {
                              int jL = aoPtr->refOrbMap(js)(jo-offsetJ).l;
                              int nOrbsJ = 2*jL+1;
                              SxDiracMat<Double> rotRight = toVector(Dmm(jL).transpose ());
                              // Rotate S 
                              block = SinR(iR).getBlock(io, jo, nOrbsI, nOrbsJ);
                              block = (rotLeft ^ block ^ rotRight);
                              SinR(iR).setBlock(block, io, jo);
                              //Rotate H
                              block = HinR(iR).getBlock(io, jo, nOrbsI, nOrbsJ);
                              block = (rotLeft ^ block ^ rotRight);
                              HinR(iR).setBlock(block, io, jo);
                              jo += nOrbsJ;
                           }
                           io += nOrbsI;
                        }
                     }

                     // --- Print out matrices
                     if (write) {
                        cout << "After rotation" << endl;
                        cout << "S = [";
                        SxDiracMat<SMatType> block = SinR(iR).getBlock(offsetI, offsetJ, noI, noJ);
                        for (int iRow = 0; iRow < block.nRows(); iRow++)  {
                           cout << "[";
                           for (int iCol = 0; iCol < block.nCols(); iCol++)  {
                              sxprintf("%.8f",block(iRow,iCol));
                              if (iCol < block.nCols() - 1 ) sxprintf(", ");
                              else sxprintf("]");
                           }
                           if (iRow < block.nRows() - 1 ) sxprintf(",\n");
                              else sxprintf("];\n");
                        }
                        cout << "H = [";
                        block = HinR(iR).getBlock(offsetI, offsetJ, noI, noJ);
                        for (int iRow = 0; iRow < block.nRows(); iRow++)  {
                           cout << "[";
                           for (int iCol = 0; iCol < block.nCols(); iCol++)  {
                              sxprintf("%.8f",block(iRow,iCol));
                              if (iCol < block.nCols() - 1 ) sxprintf(", ");
                              else sxprintf("]");
                           }
                           if (iRow < block.nRows() - 1 ) sxprintf(",\n");
                              else sxprintf("];\n");
                        }
                     } else  {
                        cout << "S:" << endl;
                        for (int io = offsetI; io < offsetI + noI; ++io)  {
                           for (int jo = offsetJ; jo < offsetJ + noJ; ++jo)  {
                              sxprintf ("% .4f ", SinR(iR)(io, jo));
                           }
                           cout << endl;
                        }
                        cout << "H:" << endl;
                        for (int io = offsetI; io < offsetI + noI; ++io)  {
                           for (int jo = offsetJ; jo < offsetJ + noJ; ++jo)  {
                              sxprintf ("% .4f ", HinR(iR)(io, jo));
                           }
                           cout << endl;
                        }
                     }
                     offsetJ += noJ;
                  } // ja
               } // js
            } // iR
            offsetI += noI;
         } // ia
      }// is

      if (write)  {
         for (int iSpecies = 0; iSpecies < structure.getNSpecies (); iSpecies++)  {
            SxString fileName = "matElem-" + SxString(iSpecies) + ".sxb";
            SxBinIO io(fileName, SxBinIO::BINARY_WRITE_ONLY);
            io.setMode (SxBinIO::WRITE_HEADER);
            writeMatElems(io, SinR, HinR, structure, vecR, aoBasis, meshR, iSpecies);
            io.setMode (SxBinIO::WRITE_DATA);
            writeMatElems(io, SinR, HinR, structure, vecR, aoBasis, meshR, iSpecies);
            io.close();
         }
      }

      if (write)  {
         SxString SName = "SElemAbs.dat";
         SxString HName = "HElemAbs.dat";
         FILE *Sfp = fopen (SName.ascii (), "w");
         FILE *Hfp = fopen (HName.ascii (), "w");
         for (int iR = 0; iR < nR; iR++)  {
            for (int iOrb=0; iOrb < nOrb; iOrb++)  {
               int is = aoBasis.orbitalMap(iOrb).is;
               int ia = aoBasis.orbitalMap(iOrb).ia;
               for (int jOrb=iOrb; jOrb < nOrb; jOrb++)  {
                  int js = aoBasis.orbitalMap(jOrb).is;
                  int ja = aoBasis.orbitalMap(jOrb).ia;
                  double dist = (structure(js,ja) + vecR(iR) 
                        - structure(is,ia)).norm ();
                  sxfprintf(Sfp,"%.6f",dist);
                  sxfprintf(Hfp,"%.6f",dist);
                  sxfprintf(Sfp,"\t%e\n",fabs(SinR(iR)(iOrb,jOrb)));
                  sxfprintf(Hfp,"\t%e\n",fabs(HinR(iR)(iOrb,jOrb)));
               }
            }
         }
         fclose(Sfp);
         fclose(Hfp);
      }
   }

   //Grep contributions
   /*
   if (write)  {
      SxCell superCell(structure.cell.basis(0) * meshR(0),
                       structure.cell.basis(1) * meshR(1),   
                       structure.cell.basis(2) * meshR(2));
      SxArray<SxArray<SxArray<SxVector<Double> > > > toWrite(structure.getNSpecies ());
      for (int is = 0; is < structure.getNSpecies(); is ++)  
         toWrite(is).resize(structure.getNSpecies());
      // Choose atom pair
      for (int iOrb = 0; iOrb < aoBasis.getNOrb(); iOrb++)  {
         SxArray<int> nonZero = createTBInteraction (aoBasis,iOrb);
         int length = nonZero.getSize() + 1;
         SxVector<Double> data(length);
         for (int iR = 0; iR < nR; ++iR)  {
            Coord d = structure(js,ja) + vecR(iR) - structure(is,ia);
            superCell.map (&d, SxCell::WignerSeitz);
            data(0) = d.norm();
            for (int j = 0; j < nonZero.getSize(); j++)  {
               int jOrb = nonZero(iOrb)(j);
               data (j+1) = HinR(iR)(iOrb, jOrb);
            }
         }
      }

      for (int is = 0; is < structure.getNSpecies(); is ++)  { 
         for (int js = 0; js < structure.getNSpecies(); js ++)  {
            FILE *TBFilePtr;
            SxString TBFile = "TB-" + SxString(is) + "-" + SxString(js) + ".dat";
            if ( (TBFilePtr = fopen(TBFile.ascii(), "w")) == NULL)  {
               cout << "Cannot open " << TBFile.ascii() << endl;
               SX_EXIT;
            }
            for (int i = 0; i < toWrite(is)(js).getSize(); i++)  {
               for (int j = 0; j < toWrite(is)(js)(i).getSize(); j++)  {
                  fprintf(TBFilePtr, "%lf\t", toWrite(is)(js)(i)(j));
               }
               fprintf(TBFilePtr, "\n");
            }
         }
      }
   }
   */
   
   
   /*
   for (int ik = 0; ik < kPoints.getNk (); ++ik)  {
      Coord kVec = kPoints.getK (ik);
      SxDiracMat<Complex16> Hk, Sk, L;
      Sk = aoMatrixRtoK (SinR, structure, meshR, aoBasis, kVec);
      //cout << nice(Sk) << endl;
      //cout << "S: " << Sk.eigenvalues ().real () << endl;
      Hk = aoMatrixRtoK (HinR, structure, meshR, aoBasis, kVec);
      //cout << "H: " << Hk.eigenvalues ().real () << endl;
      L = Sk.choleskyDecomposition (UpperRight).inverse ();
      cout << "E = "
           << SxDiracMat<Complex16> (L.adjoint () ^ Hk ^ L)
              .eigenvalues ().real () * HA2EV
           << endl;
   }
   */
   if ((kFile.getSize () > 0) && !write)  {
      table = parser.read (kFile, "std/basis.std");
      SxKPoints bs (structure.cell, &*table);
      ofstream epsFile;
      if (nSpin == 1)  {
         epsFile.open ("eps-ao.dat");
      } else {
         if (iSpin == 0)
            epsFile.open ("eps-ao.up.dat");
         else if (iSpin == 1)
            epsFile.open ("eps-ao.down.dat");
         else
            SX_EXIT;
      }
      SX_START_TIMER (CalcBS);
      for (int ik = 0; ik < bs.getNk (); ++ik)  {
         Coord kVec = bs.getK (ik);
         SxDiracMat<Complex16> Hk, Sk, L;
         Sk = aoMatrixRtoK (SinR, structure, meshR, aoBasis, kVec);
         Hk = aoMatrixRtoK (HinR, structure, meshR, aoBasis, kVec);
         int nr = int(Sk.nRows ());
         if (printS)  {
            cout << SX_SEPARATOR;
            cout << "ik=" << ik << endl;
            bs.getK(ik).print();
            SxDiracMat<Complex16>::Eigensystem seig = Sk.eigensystem ();
            cout << "S: " << endl; 
            seig.vals.real ().print(true);
            cout << SX_SEPARATOR;
            if (seig.vals (0).re < 1e-12)  {
               int i = 0;
               while (seig.vals(i).re < 1e-12) i++;
               SxDiracVec<Complex16> nonSing
                  = seig.vecs (SxIdx(nr * i, nr*nr -1));
               nonSing.reshape (nr, nr-i);
               Sk = nonSing.adjoint () ^ Sk ^ nonSing;
               Hk = nonSing.adjoint () ^ Hk ^ nonSing;
            }
         }
         SX_START_TIMER (Diag);
         L = Sk.choleskyDecomposition (UpperRight).inverse ();
         Hk = L.adjoint () ^ Hk ^ L;
         SxDiracVec<Double> eps = SxDiracMat<Complex16> (Hk).eigenvalues ().real ();
         SX_STOP_TIMER (Diag);
         epsFile << (ik+1);
         for (int i = 0; i < eps.getSize (); ++i)
           epsFile << '\t' << eps(i) * HA2EV; 
         // fill in dummy eigenvalues from singular values
         for (int i = int(eps.getSize ()); i < nr; ++i)
           epsFile << '\t' << 1e8;
         epsFile << endl;
      }
      SX_STOP_TIMER (CalcBS);
      epsFile.close ();
   }
   }
   printTiming ();
   return 0;
}
