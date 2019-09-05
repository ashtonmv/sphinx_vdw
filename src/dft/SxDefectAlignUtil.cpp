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

#include <SxDefectAlignUtil.h>
#include <SxLocpot.h>
#include <SxFileParser.h>
#include <SxFFT1d.h>
#include <SxGrid.h>
#ifdef SXDA_STRUCTURE_PRINTSX
#include <SxElemDB.h>
#endif

namespace SxDefectAlignUtil {

SxMeshR getPot (SxCell &cell, SxMesh3D &mesh, SxString &fileName,
                FileType fileType, SxAtomicStructure *structure)
{
   SX_CHECK(structure);
   SxMeshR result;

   if (fileType == sxb)  {
      try {
         SxBinIO io(fileName, SxBinIO::BINARY_READ_ONLY);
         result = io.readMesh (&cell, &mesh)(0);
         io.close ();
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
   } else if (fileType == VASP_LOCPOT)  {
      try  {
         SxLocpot locpot(fileName);
         locpot.read ();
         cell = locpot.getCell ();
         mesh = locpot.getMesh ();
         result = locpot.getPotential ();
         *structure = locpot.getStructure ();
      } catch (SxException e) {
         e.print ();
         SX_EXIT;
      }
   } else if (fileType == Socorro)  {
      SxFileParser fp(fileName);
      fp.read ("LATTICE VECTORS");
      fp >> cell(0,0) >> cell(1,0) >> cell(2,0)
         >> cell(0,1) >> cell(1,1) >> cell(2,1)
         >> cell(0,2) >> cell(1,2) >> cell(2,2);
      cell.setup ();
      fp.read ("ATOMS");
      int na = fp.getInt ();
      structure->cell = cell;
      for (int ia = 0; ia < na; ++ia)  {
         Coord posRel;
         fp >> posRel(0) >> posRel(1) >> posRel(2);
         structure->addAtom (cell.relToCar (posRel));
      }
      structure->endCreation ();
      fp.read ("MESH SIZES");
      fp >> mesh(0) >> mesh(1) >> mesh(2);
      result.resize (mesh.getSize ());
      fp.read ("POTENTIAL");
      for (int i3 = 0; i3 < mesh(2); i3++)
         for (int i2 = 0; i2 < mesh(1); i2++)
            for (int i1 = 0; i1 < mesh(0); i1++)
               fp >> result(mesh.getMeshIdx(i1,i2,i3,SxMesh3D::Positive));
      result *= 0.5 * HA2EV; // Rydberg to eV
   } else if (fileType == QuantumEspresso)  {
      SxFileParser fp(fileName);
      fp.verbose = true;

      // read title lines
      fp.nextLine (2);
      // line 3: number of atoms, coordinate system offset
      fp.topic ("mesh+cell");
      Coord offset;
      int na = fp.getInt ();
      fp >> offset(0) >> offset(1) >> offset(2);
      if (offset.normSqr () > 1e-12)  {
         cout << "ERROR: quantum espresso mesh file '" << fileName 
              << "' contains a non-zero mesh offset in line 2: " << offset
              << endl;
         cout << "I do not know how to deal with the situation..." << endl;
         cout << "Please contact the developer (freysoldt@mpie.de)" << endl;
         SX_QUIT;
      }
      fp.nextLine ();

      // line 4-6: mesh and unit cell
      for (int i = 0; i < 3; ++i)  {
         Coord meshDir;
         fp >> mesh(i);
         for (int j = 0; j < 3; ++j)
            cell(j,i) = fp.getDouble () * mesh(i);
         fp.nextLine ();
      }
      //cell /= A2B;
      cell.setup ();
      structure->cell = cell;
      cout << "mesh=" << mesh << endl;

      // --- line 7..6+N
      fp.topic ("atoms");
      SxList<int> speciesMap;
      for (int ia = 0; ia < na; ++ia)  {
         Coord posAbs;
         int Z = fp.getInt ();
         double Zr = fp.getDouble ();
         if (fabs(Z - Zr) > 1e-6)  {
           cout << "Suspicious quantum espresso format ";
           fp.where ();
           cout << endl << "Integer atomic number is " << Z << ", but real is " << Zr << endl;
         }
         int is = (int)speciesMap.findPos (Z);
         if (is == -1) {
            cout << "New species Z=" << Z << endl;
            structure->newSpecies ();
            is = (int)speciesMap.getSize ();
            speciesMap << Z;
         }
         fp >> posAbs(0) >> posAbs(1) >> posAbs(2);
         //posAbs /= A2B;
         structure->addAtom (is, posAbs);
         fp.nextLine ();
      }
      structure->endCreation ();
#ifdef SXDA_STRUCTURE_PRINTSX
      {
         SxElemDB db;
         SxArray<SxString> elem(speciesMap.getSize ());
         for (ssize_t is = 0; is < speciesMap.getSize (); ++is)  {
            elem(is) = db.getChemSymbol (speciesMap(is));
         }
         structure->atomInfo->meta.update (SxAtomicStructure::Elements, elem);
      }
#endif

      // --- line 7+N ...
      fp.topic ("potential");
      result.resize (mesh.getSize ());
      for (int i1 = 0; i1 < mesh(0); i1++)
         for (int i2 = 0; i2 < mesh(1); i2++)
            for (int i3 = 0; i3 < mesh(2); i3++)
               fp >> result(mesh.getMeshIdx(i1,i2,i3,SxMesh3D::Positive));
      result *= 0.5 * HA2EV; // Rydberg to eV
   } else {
      cout << "Unsupported potential type! SxDefectAlign quits here!" << endl;
      SX_QUIT;
   }
   cell.setup ();

   return result;
}

SxVector<Double> readLine (const SxCell &potCell,
                           const SxMesh3D &potMesh,
                           const SxMeshR &potData,
                           int idir,
                           ssize_t n,
                           const SxCell &cell,
                           const SxString &file)
{

   // --- average over two dimensions
   SxVector<Double> res(potMesh(idir));
   res.set (0.);
   for (int i = 0; i < potData.getSize (); ++i)  {
      res(potMesh.getMeshVec(i, SxMesh3D::Positive)(idir)) += potData(i);
   }
   res /= potMesh.product () / potMesh(idir);
   
   // --- determine multiplicity
   Coord relBasis = potCell.carToRel (cell.basis(idir));
   for (int d = 0; d < 3; ++d)  {
      if (idir != d && fabs(relBasis(d)) > 1e-6)  {
        cout << "Incompatible cell in file " << file << endl;
        cout << "a" << (idir+1) 
             << " in " << file << "'s system:" << relBasis << endl;
        SX_EXIT;
      }
   }
   int mul = (int)lround(relBasis(idir));
   if (fabs(mul - relBasis(idir)) > 1e-6 || mul <= 0)  {
        cout << "Incompatible cell in file " << file << endl;
        cout << "a" << (idir+1) 
             << " in " << file << "'s system: " << relBasis << endl;
        SX_EXIT;
   }
   if (mul > 1) cout << "multiplicity = " << mul << endl;
   if (n != potMesh(idir) || mul > 1)  {
      // --- Fourier interpolation
      SxVector<Complex16> resC(res);
      SxVector<Complex16> oldInG(potMesh(idir)), newInG(n);
      
      {
         SxFFT1d oldFFT(SxFFT::Reverse, potMesh(idir));
         oldFFT.fftReverse (potMesh(idir), resC.elements, oldInG.elements);
      }
      
      newInG.set (0.);
      newInG(0) = oldInG(0);
      int iLimit = min(potMesh(idir), int(n/mul));
      for (int i = 1; 2 * i < iLimit; ++i)  {
         newInG(mul * i) = oldInG(i);
         newInG(n-mul*i) = oldInG(potMesh(idir)-i);
      }

      resC.resize (n);
      {
         SxFFT1d newFFT(SxFFT::Forward, (int)n);
         newFFT.fftForward (int(n), newInG.elements, resC.elements);
      }
      res = resC.real ();
   }
   return res;
}

SxVector<Double> average (const SxVector<Double> &x, double w)
{
   int nAvg = int(floor(0.5 * w)), n = int(x.getSize ()), idx;
   double fw = w - (2*nAvg-1), sum;
   SxVector<Double> res(x.getSize ());
   for (int i = 0; i < n; ++i)  {
      // lower boundary
      idx = (i - nAvg) % n;
      if (idx < 0) idx += n;
      sum = 0.5 * x(idx);
      // upper boundary
      idx = (i + nAvg) % n;
      sum += 0.5 * x(idx);
      sum *= fw;

      // integer average
      for (int j = 1-nAvg; j < nAvg; ++j)  {
         idx = (i + j) % n;
         if (idx < 0) idx += n;
         sum += x(idx);
      }
      res(i) = sum / w;
   }
   return res;
}

int getMappedAtom (const Coord &pos,
                   SxAtomicStructure &structure,
                   const int iSpecies)
{
   int result = -1;
   // grid for neighbor search
   SxGrid grid (structure, 10);
   
   double epsEqualOrig = structure.epsEqual;
   structure.epsEqual = 0.5;
   int idx = structure.find(pos, grid);

   if (idx > -1)  {
      int species = structure.getISpecies(idx, &result);
      if (species != iSpecies) result = -1;
   }
   structure.epsEqual = epsEqualOrig;

   return result;
}

FileType getFileType (SxCLI &cli)
{
   FileType fileType = sxb;
   if (cli.option ("--vasp", "potentials are VASP LOCPOT files").toBool ())  {
      fileType = VASP_LOCPOT;
   }
   if (cli.option ("--socorro", "potentials are socorro files").toBool ())  {
      fileType = Socorro;
   }
   if (cli.option ("--qe", "potentials are quantum espresso files").toBool ())  {
      fileType = QuantumEspresso;
   }
   return fileType;
}
   
}


