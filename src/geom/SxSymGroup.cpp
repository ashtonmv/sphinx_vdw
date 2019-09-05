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

#include <SxSymGroup.h>
#include <SxCell.h>
#include <SxRotation.h>

SxSymGroup::SxSymGroup(const SxArray<SymMat> &symmorphic, 
                       const SxMatrix3<TPrecTauR> &cell,
                       const SxMesh3D &meshIn)
 : mesh(meshIn),
   nSymmorphic(int(symmorphic.getSize ()))
{
   primCell = cell; // may be empty
   SX_CHECK (fabs(cell.determinant ()) > 1e-12 || mesh.product () <=1);
   primOps.resize (nSymmorphic);
   for (int i = 0; i < nSymmorphic; ++i)  {
      primOps(i).rot = symmorphic(i);
      primOps(i).trans.set(0.);
   }
   check ();
}

SxSymGroup::SxSymGroup(const SxCell          &primCellIn,
                       const SxArray<SxSymOp> &primOpsIn, 
                       const SxCell          &cell)
 : primOps(primOpsIn)
{
   int nSym = int(primOps.getSize ());
   // --- symmorphic symmetries must be first nSymmorphic syms
   nSymmorphic = 0;
   for (int iSym = 0; iSym < nSym; ++iSym)  {
      if (primOps(iSym).isSymmorphic ())  {
         if (nSymmorphic < iSym)  {
            // bring symmorphic symmetry to front
            SxSymOp aux = primOps(iSym);
            primOps(iSym) = primOps(nSymmorphic);
            primOps(nSymmorphic) = aux;
         }
         nSymmorphic++;
      }
   }
   // set up internal periodicity
   SxCell genCell(primCellIn);
   if (cell.volume > 0.)
      mesh = genCell.computeGen (cell);
   else
     mesh = SxMesh3D (1,1,1);;
   primCell = genCell;
   check ();
}

SxArray2<int>
SxSymGroup::getSymMulTable (const SxArray<SymMat> &syms, double epsSym)
{
   int nSym = int(syms.getSize ());
   if (nSym == 0) return SxArray2<int> ();

   SxArray2<int> table(nSym, nSym);
   // --- Checking: are the symmorphic symmetries a full group?
   SymMat result;
   bool found;
   SxArray<bool> isFound (nSym);
   for (int i = 0; i < nSym; i++)  {
      isFound.set (false);
      for (int j = 0; j < nSym; j++)  {
         result =syms(i) ^ syms(j);
         // find result
         found = false;
         for (int k = 0; k < nSym; k++)  {
            if (isFound(k)) continue; // each result can be obtained only once
            if ((syms(k) - result).absSqr().sum() > epsSym)
               continue;
            isFound(k) = true;
            table(i,j) = k;
            found = true;
            break;
         }
         if (!found)  {
            cout << "Symmetry inconsistency error: symmetries are no group\n";
            for (int k = 0; k < nSym; ++k)  {
                cout << "Symmetry " << (k+1) << ": " 
                     << syms(k) << " = " << SxRotation(syms(k)).getName () << endl;
            }
            cout << "epsSym = " << epsSym << endl;
            cout << "Offending multiplication:" << endl;
            cout << SxRotation(syms(i)).getName () << " ^ "
                 << SxRotation(syms(j)).getName () << " = "
                 << SxRotation(result).getName () << endl;
            cout << "cartesian coordinates:" << endl;
            cout << syms(i) << " ^ ";
            cout << syms(j) << " = ";
            cout << result << endl;
            /*
            if (cell.volume > 0.)  {
               cout << "relative coordinates:" << endl;
               cout << cell.carToRel(getSymmorphic(i)) << " ^ ";
               cout << cell.carToRel(getSymmorphic(j)) << " = ";
               cout << cell.carToRel(result) << endl;
            }
            */
            SX_EXIT;
         }
      }
   }
   return table;
}

void SxSymGroup::check () const
{
   int nSym = getNSymmorphic ();
   SxCell cell;
   if (primCell.determinant () > 1e-12 && mesh.normSqr () > 0)  {
      cell = SxCell(primCell.col(0) * mesh(0),
                    primCell.col(1) * mesh(1),
                    primCell.col(2) * mesh(2));
   }
   // --- Checking: are symmorphic symmetries compatible with lattice?
   if (cell.volume > 0.)  {
      for (int i = 0; i < nSym; i++)  {
         for (int iDir = 0; iDir < 3; ++iDir)  {
            Coord rotBasis = getSymmorphic(i) ^ cell.basis(iDir);
            cell.map(&rotBasis, SxCell::Origin);
            if (rotBasis.norm () > cell.epsSym)  {
               cout << "Symmetry group inconsistency." << endl;
               cout << "Symmetry " << getSymmorphic(i)
                    << " is incompatible with lattice."
                    << endl;
               cout << "Basis vector " << (iDir+1) << " differs from lattice by "
                    << rotBasis << endl;
               SX_EXIT;
            }
         }
      }
   }
   // --- Checking: are the symmorphic symmetries a full group?
   getSymMulTable (getSymmorphic (), cell.epsSym);

   if (primOps.getSize () > nSymmorphic)  {
      cout << "TODO: check for non-symmorphic symmetry consistency" << endl;
   }
}

SxSymOp SxSymGroup::operator()(int iSym) const
{
   SX_CHECK (iSym >= 0 && iSym < getSize (), iSym, getSize ());
   int idx = int(iSym / primOps.getSize ());
   ssize_t iOp = iSym % primOps.getSize ();
   SxSymOp res(primOps(iOp));
   res.trans += primCell ^ mesh.getMeshVec (idx,SxMesh3D::Positive);
   return res;
}

SxArray<SxSymOp> SxSymGroup::operator/ (const SxArray<SxSymOp> &kernel) const 
{
   ssize_t iK, nK = kernel.getSize (),
       nSym = getSize ();
   if (nK == 0)  {
      cout << "Error in SxSymGroup::operator/: empty kernel provided" << endl;
      SX_EXIT;
   }  else if (nK == nSym)  {
      cout << "Warning from SxSymGroup::operator/: all syms in kernel?" << endl;
   }
   SxList<SxSymOp> res; 
   SxArray<SxSymOp> symDone (nK);
   for (iK = 0; iK < nK; iK++)  { symDone(iK) = kernel(iK); }
   cout << "SxSymGroup::" << getSize () << "nSym / nKern" << nK << endl;

   SxSymOp symId;
   int iDone, nDone = -1;
   CellMat primInv = primCell.inverse ();
   SxVector3<Double> relPrim;
   SxVector3<Int> nullVec (0, 0, 0), relPrimInt;
   for (int iSym = 0; iSym < nSym; iSym++)  {
      cout << iSym;
      nDone = int(symDone.getSize ());
      for (iDone = 0; iDone < nDone; iDone++)  {//List.contains not usable: periodicity
         symId = (operator() (iSym)) ^ symDone(iDone).inverse ();
         relPrim = primInv ^ symId.trans;
         relPrimInt(0) = int(lround(relPrim(0)));
         relPrimInt(1) = int(lround(relPrim(1)));
         relPrimInt(2) = int(lround(relPrim(2)));
         if (symId.rot.trace () > 3 -1e-2 
         && (relPrim - relPrimInt).normSqr () < 1e-3 
         && (relPrimInt % mesh) == nullVec)
         { cout << "found "; break; }  
      }
      if (iDone == nDone)  {//symOp not found
         cout << (operator() (iSym)) << "added " << endl;
         res << (operator() (iSym));
         symDone.resize (nDone + nK, true);
         for (iK = 0; iK < nK; iK++)  //list is slow
         {  symDone(nDone + iK) = (operator() (iSym)) ^ kernel(iK); }
         if (nDone + nK > nSym)  {
            cout << "Warning from SxSymGroup::operator/: too many symmetries generated!" << endl
                 << "Kernel has redundant symms, or symms not recognised?" << (nDone + nK) << endl;
         }  
      }
   }
   cout << endl;
   if (nDone + nK < nSym)  {
      cout << "Warning from SxSymGroup::operator/: too few symmetries generated!" << endl
           << "Original set contains doubles?" << (nDone + nK) << endl;
   }
   return res; 
}

SxArray<SymMat> SxSymGroup::getSymmorphic () const
{
   SxArray<SymMat> res(nSymmorphic);
   for (int i = 0; i < nSymmorphic; ++i)
      res(i) = primOps(i).rot;
   return res;
}

void SxSymGroup::read(const SxBinIO &io, const SxCell &cell)
{
   primCell = cell;
   mesh = SxMesh3D (1,1,1);
   SxArray<SymMat> ops;
   try {
      nSymmorphic = io.getDimension ("nSymmetries");
      ops.resize (nSymmorphic);
      io.read ("symOps", &ops);
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
   primOps.resize (nSymmorphic);
   for (int i = 0; i < nSymmorphic; ++i)  {
      primOps(i).rot = ops(i);
      primOps(i).trans.set(0.);
   }
   check ();
}

void SxSymGroup::write(SxBinIO &io) const
{
   if (nSymmorphic == 0) return;
   try  {
      io.addDimension ("xyz", 3);
      io.addDimension ("nSymmetries", nSymmorphic);
      io.write("symOps", getSymmorphic (), "nSymmetries", "xyz");
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}

void SxSymGroup::print () const
{
   cout << SX_SEPARATOR;
   cout << "| Symmetry group (" << getSize () << " symmetries)" << endl
        << SX_SEPARATOR;
   if (mesh.product () != 0 && fabs(primCell.determinant ()) > 1e-12)  {
      cout << "| Pure translations: "
         << mesh(0) << 'x' << mesh(1) << 'x' << mesh(2) << endl;
      cout << "| with basis " << SxCell(primCell) << endl;
      cout << SX_SEPARATOR;
   }
   cout << "| Primitive symmetries: " << primOps.getSize () 
        << " (" << nSymmorphic << " symmorphic)" << endl;
   for (int i = 0; i < primOps.getSize (); ++i)
      cout << "| " << primOps(i) 
                   << " " << SxRotation(primOps(i).rot).getName () << endl;
   cout << SX_SEPARATOR;
   cout.flush ();
}

void SxSymGroup::fprintsx (FILE *output, bool nonsymmorphic) const
{
   int nSym = nonsymmorphic ? getNPrimitive () : getNSymmorphic ();
   sxfprintf(output, "   symmetry {\n");
   if (mesh.product () > 1 && fabs(primCell.determinant ()) > 1e-12)  {
      sxfprintf (output,
                 "      primCell = [[% 14.8f, % 14.8f, % 14.8f]\n"
                 "                  [% 14.8f, % 14.8f, % 14.8f]\n"
                 "                  [% 14.8f, % 14.8f, % 14.8f]];\n",
                 primCell(0,0), primCell(1,0), primCell(2,0),
                 primCell(0,1), primCell(1,1), primCell(2,1),
                 primCell(0,2), primCell(1,2), primCell(2,2));
      sxfprintf (output, "      superMesh = [%d, %d, %d];\n",
                 mesh(0), mesh(1), mesh(2));
   }
   for (int iSym = 0; iSym < nSym; ++iSym)  {
      if (iSym == getNSymmorphic ())
         sxfprintf(output, "\n      // --- nonsymmorphic primitives\n");
      SxSymOp fullS = (*this)(iSym);
      SymMat &S = fullS.rot;
      sxfprintf(output, "      // Symmetry Number %i, %s\n"
                        "      operator { S = [[%f,%f,%f],\n"
                        "                      [%f,%f,%f],\n"
                        "                      [%f,%f,%f]];",
                        iSym+1, SxRotation::getName(S).ascii(),
                        S(0,0), S(0,1), S(0,2), S(1,0), S(1,1), S(1,2),
                        S(2,0), S(2,1), S(2,2));
      if (!fullS.isSymmorphic ())
         sxfprintf(output, "\n                 shift = [%.8f, %.8f, %.8f];",
                           fullS.trans(0), fullS.trans(1), fullS.trans(2) );
      sxfprintf(output, " }\n");

   }
   sxfprintf(output, "   }\n");
}

void SxSymGroup::read (const SxSymbolTable *table)
{
   SX_CHECK (table);
   if (table->containsGroup ("symmetry"))
      table = table->getGroup ("symmetry");
   SxStack<SxSymOp> syms;
   try {
      SxSymbolTable *sym;
      if (table->containsGroup ("operator"))  {
         for ( sym  = table->getGroup ("operator");
               sym != NULL;
               sym  = sym->nextSibling ("operator"))
         {
            SymMat S(sym->get("S")->toList());
            Coord trans(0.,0.,0.);
            if (sym->contains("shift"))
               trans = Coord(sym->get("shift")->toList ());
            syms << SxSymOp (S, trans);
         }
      }
      if (table->contains ("superMesh"))  {
         mesh = SxVector3<Int> (table->get("superMesh")->toIntList ());
         primCell = CellMat (table->get ("primCell")->toList ()).transpose ();
      } else {
         primCell = SxCell ();
         mesh.set (1);
      }
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   if (syms.getSize () == 0)  {
      // no symmetries ? => identity
      syms << SxSymOp (SymMat(1,0,0,0,1,0,0,0,1));
   }
   primOps = syms;
   int nSym = int(primOps.getSize ());
   // --- symmorphic symmetries must be first nSymmorphic syms
   nSymmorphic = 0;
   for (int iSym = 0; iSym < nSym; ++iSym)  {
      if (primOps(iSym).isSymmorphic ())  {
         if (nSymmorphic < iSym)  {
            // bring symmorphic symmetry to front
            SxSymOp aux = primOps(iSym);
            primOps(iSym) = primOps(nSymmorphic);
            primOps(nSymmorphic) = aux;
         }
         nSymmorphic++;
      }
   }
}

