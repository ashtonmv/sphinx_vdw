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

#include <SxCell.h>
#include <SxMathLib.h>
#include <SxSymGroup.h>
#include <SxRotation.h>


// ------- constructor
SxCell::SxCell () 
   : CellMat (0,0,0,0,0,0,0,0,0),
     inv (0,0,0,0,0,0,0,0,0),
     volume (0.),
     epsRegular (EPS_REGULAR_DEFAULT),
     epsSym (EPS_SYM_DEFAULT),
     reg (0,0,0,0,0,0,0,0,0),
     regInv (0,0,0,0,0,0,0,0,0)
{ 
}

SxCell::SxCell (const CellMat& inCell) 
   : CellMat (inCell),
     epsRegular (EPS_REGULAR_DEFAULT),
     epsSym (EPS_SYM_DEFAULT)
{ 
   setup ();
}

SxCell::SxCell (const CellMat& inCell, const SxPtr<SxSymGroup> &symGrp) 
   : CellMat (inCell),
     epsRegular (EPS_REGULAR_DEFAULT),
     epsSym (EPS_SYM_DEFAULT)
{ 
   symGroupPtr = symGrp; // may be empty
   setup ();
}

SxCell::SxCell (const SxCell &cell) 
   : CellMat (cell),
     epsRegular (cell.epsRegular),
     epsSym (cell.epsSym)
{
   symGroupPtr = cell.symGroupPtr; // may be empty
   setup ();
}

SxCell::SxCell (const Coord& a1, const Coord& a2, const Coord &a3)
   : epsRegular (EPS_REGULAR_DEFAULT),
     epsSym (EPS_SYM_DEFAULT)
{
   set (a1, a2, a3);
}

SxCell::SxCell (const SxList<Coord> &list)
   : epsRegular (EPS_REGULAR_DEFAULT),
     epsSym (EPS_SYM_DEFAULT)
{
   set (list(0), list(1), list(2));
}

SxCell::SxCell (const SxBinIO &io)
   : epsRegular (EPS_REGULAR_DEFAULT),
     epsSym (EPS_SYM_DEFAULT)
{
   read (io);
}

SxCell::SxCell (const SxSymbolTable *table) 
   : epsRegular (EPS_REGULAR_DEFAULT),
     epsSym (EPS_SYM_DEFAULT)
{ 
   SX_CHECK (table);
   const SxSymbolTable *str;
   const SxString strName = "structure";
   try {
      // get structure group
      str = (table->name == strName) ? table : table->getGroup (strName);

      // set cell
      // Basis vectors are entered as rows, but are cols internally.
      // That's why we transpose.
      CellMat::operator= (CellMat(str->get("cell")->toList ()).transpose ());

      if (table->contains ("epsSym"))
         epsSym = table->get ("epsSym")->toReal ();

      // --- structure.symmetry
      if (str->containsGroup ("symmetry"))  {
         symGroupPtr = SxPtr<SxSymGroup>::create ();
         symGroupPtr->read (table);
      }
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   setup ();
}

SxCell::~SxCell ()
{
   /* empty */
}



SxCell& SxCell::operator= (const SxCell &in)
{
   CellMat::operator= (in);
   volume      = in.volume;
   inv         = in.inv;
   reg         = in.reg;
   regInv      = in.regInv;
   symGroupPtr = in.symGroupPtr;
   epsSym      = in.epsSym;
   return *this;
}


void SxCell::set (const CellMat &inCell)
{
   CellMat::operator= (inCell);
   setup ();
}

void SxCell::set (const Coord &a1, const Coord &a2, const Coord &a3)
{
   *((CellMat*)this) = CellMat(a1(0), a2(0), a3(0),
                               a1(1), a2(1), a3(1),
                               a1(2), a2(2), a3(2));
   setup ();
}

void SxCell::write (SxBinIO &io) const
{
   try {
      io.addDimension ("xyz", 3);

      // write cell (don't try to declare it twice, though)
      if (!io.contains("cell") || (io.ncMode == SxBinIO::WRITE_DATA))
         io.write ("cell", static_cast<const SxMatrix3<Double>&>(*this), "xyz");

      // --- symmetries
      if (symGroupPtr) symGroupPtr->write (io);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}

void SxCell::read (const SxBinIO &io)
{
   try {
      io.read ("cell", static_cast<CellMat*>(this));
      if (io.contains("symOps"))  {
         // --- read symmetries
         symGroupPtr = SxPtr<SxSymGroup>::create ();
         symGroupPtr->read (io, *this);
      } else {
         // no symmetries available
         symGroupPtr = SxPtr<SxSymGroup> ();
      }
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   setup ();
}

int SxCell::setFromVectors (const SxList<Coord> &points)
{
   const SxList<Coord> *thePoints = &points;
   SxPtr<SxList<Coord> > updatedPoints, newPoints;
   SxList<Coord>::ConstIterator ptIt;
   int nDim;
   Coord mappedPoint;
   for (int count = 0; count < 100; ++count)  {

      if (count > 50)  {
         // --- try vector reduction
         SxList<Coord>::Iterator a, b;
         double aSqr, abSqr;
         for (a = updatedPoints->begin (); a != updatedPoints->end (); ++a)  {
            aSqr = (*a).normSqr ();
            for (b = updatedPoints->begin (); b != updatedPoints->end (); ++b) {
               if (a == b) continue; // no self-reduction!
               // factor 0.99 ignores negligible improvements/numerical noise
               if ((abSqr = (*a + *b).normSqr ()) < 0.99 * aSqr)  {
                  if (abSqr > 1e-10)  {
                     *a += *b;
                     aSqr = abSqr;
                  }
               } else if ((abSqr = (*a - *b).normSqr ()) < 0.99 * aSqr)  {
                  if (abSqr > 1e-10)  {
                     *a -= *b;
                     aSqr = abSqr;
                  }
               }
            }
         }
      }
      nDim = setFromGoodVectors (*thePoints);
      
      if (nDim < 2) return nDim; // no extra treatment here
      if (nDim == 2) setCol(2, col(0).x(col(1)).normalize ());
      setup ();
      
      // --- check that all points are lattice points
      for (ptIt = thePoints->begin (); ptIt != thePoints->end (); ++ptIt)  {
         mappedPoint = getMapped(*ptIt, Origin);
         if (mappedPoint.normSqr () > 1e-4)  {
            // add non-lattice points for next cycle
            if (!newPoints) newPoints = newPoints.create ();
            (*newPoints) << mappedPoint;
         }
      }
      if (!newPoints)   {
         // everything is fine
         if (nDim == 2) setCol(2, Coord(0.,0.,0.));
         return nDim;
      }
      // add current basis vectors
      for (int d = 0; d < nDim; ++d) newPoints->append(col(d));
      
      // --- switch to new list of points
      updatedPoints = newPoints;
      thePoints = updatedPoints.getPtr ();
      newPoints = SxPtr<SxList<Coord> > ();
   }
   // --- crash if we can't find lattice after 100 cycles
   cout << "points:\n" << points << endl;
   cout << "guessed cell:" << (*this) << endl;
   cout << "Internal error in SxCell::setFromVectors, SxCell.cpp:" << __LINE__
        << ".\nPoints seem to be not a lattice." << endl;
   SX_EXIT;
   return -1; // avoid warning
}


int SxCell::setFromGoodVectors (const SxList<Coord> &points)
{
   ssize_t n = points.getSize ();
   CellMat::set(0.);

   // find 1st basis vector: shortest non-zero vector
   double maxLength = -1.,l;
   SxList<Coord>::ConstIterator it;
   int i;
   Coord a1, a2, a3;
   for (i = 0, it=points.begin (); i < n; ++i, ++it)  {
      l = (*it).normSqr ();
      if (l > 1e-7 && (maxLength < 0. || l < maxLength))  {
         maxLength = l;
         a1 = *it;
      }
   }
   if (maxLength < 0.) return 0;
   setCol(0, a1);

   // find 2nd basis vector: shortest "most orthogonal"
   double scp, minScp = -1.;
   maxLength = -1.;
   Coord vec;
   for (i = 0, it = points.begin (); i < n; ++i, ++it)  {
      vec = *it;
      vec -= round ((vec ^ a1)/a1.normSqr ()) * a1; // most orthogonal
      if ( ((vec).x(a1)).normSqr () < 1e-7) continue; // colinear
      scp = fabs(a1 ^ vec);
      if (scp < 1e-7)  {
         minScp = 0.;
         if (maxLength < 0. || vec.normSqr () < maxLength)  {
            maxLength = vec.normSqr ();
            a2 = vec;
         }
      } else if (minScp < 0. || scp < minScp)  {
         minScp = scp;
         a2 = vec;
      }
   }

   if (a2.normSqr () < 1e-7) return 1;
   setCol(1, a2);
   
   // find 3rd basis vector: shortest "most near"
   Coord normal = a1.x(a2);
   normal.normalize ();
   maxLength = minScp = -1.;
   for (i = 0, it = points.begin (); i < n; ++i, ++it)  {
      scp = fabs(normal ^ (*it) );
      if (scp < 1e-7) continue; // in-plane
      if ( (scp < minScp + 1e-7) || (minScp < 0.) )  {
         minScp = scp;
         if (maxLength < 0. || (*it).normSqr () < maxLength)  {
            maxLength = (*it).normSqr ();
            a3 = *it;
         }
      }
   }
   if (a3.normSqr () < 1e-7)  {
      return 2;
   }
   setCol (2, a3);
   *this = getRegularCell ();
   return 3;
      
}

void SxCell::setup ()
{
  volume = fabs (this->determinant ());
  SX_CHECK(volume > 1e-10, volume); 
  if (volume < 1e-10)  {
    cout << "Invalid cell: " << endl;
    cout << (*this) << endl; 
    cout << "Cell has no volume!" << endl;
    SX_QUIT;
  }
  inv = CellMat::inverse ();

  // --- Wigner-Seitz mapper
  reg = getRegularLattice (*this, epsRegular);
  regInv = reg.inverse ();
}

PrecTauR SxCell::shortestDist (const Coord &cart1, const Coord &cart2) const
{
   Coord d = cart1 - cart2;
   
   // Map into Wigner Seitz Cell
   map(&d,SxCell::WignerSeitz);

   return d.norm ();
}

// --------------- regular cells
CellMat SxCell::getRegularLattice (const CellMat &cell, double eps)
{
   SX_CHECK (fabs(cell.determinant ()) > 1e-10);
   int i, j;
   PrecTauR proj;
   bool done;
   SxVector3<TPrecTauR> basis[3], redVec;
   basis[0] = cell.col(0);
   basis[1] = cell.col(1);
   basis[2] = cell.col(2);
   int count = 10000; // break do loop (should never happen)
   do  {
      done = true;
      if (count < 20)  { // this is going to fail - provide debug output
         cout << "---" << endl
              << "a1 = " << basis[0] << endl
              << "a2 = " << basis[1] << endl
              << "a3 = " << basis[2] << endl;
      }
      for (j= 0; j < 3 && done; j++)  {
         for (i = 0; i < 5; i++)  {
            if (i==j) continue;
            if (i < 3)
               redVec = basis[i];
            else if ( fabs(basis[(j+1)%3] ^ basis[(j+2)%3])
                     /(basis[(j+1)%3].norm () *  basis[(j+2)%3].norm ()) < 0.1)
               break; // don't try reduction by sum or difference of the
                      // other two vectors if they are largely orthogonal
            else if (i == 3)
               redVec = basis[(j+1)%3] + basis[(j+2)%3]; // sum of two others
            else
               redVec = basis[(j+1)%3] - basis[(j+2)%3]; // diff of two others

            proj = (redVec ^ basis[j]) / redVec.normSqr ();
            if (count < 20)  // this is going to fail - provide debug output
               cout << j << " " << i << " " << proj << endl;
            if (fabs (proj - 0.5) < eps) {
               proj = 1.0; // == round(0.5);
               // this might be going to fail, but the cell is OK
               // within the eps error bar => don't try to improve further
               if (count <= 20) continue;
            }
            if (fabs (proj + 0.5) < eps) proj = 0.; // -0.5 is fine
            proj = round(proj);
            if (fabs(proj) > 1e-12)  {
               done = false;
               basis[j] -= proj * redVec;
               break;
            }
         }
      }
   } while (!done && --count);
   if (!count)  {
      cout << "Infinite loop error in SxCell::getRegularLattice." << endl;
      cout << "Inform the developers, attach input file" << endl;
      SX_EXIT;
   }
   return CellMat(basis[0](0), basis[1](0), basis[2](0),
                  basis[0](1), basis[1](1), basis[2](1),
                  basis[0](2), basis[1](2), basis[2](2));
}

/// get heights between parallel planes
Coord SxCell::getHeights () const
{
   // lattice vectors
   const Coord a1 = basis(0);
   const Coord a2 = basis(1);
   const Coord a3 = basis(2);

   // normal vectors on the planes
   Coord  n1 = (a2.x (a3)).normalize();
   Coord  n2 = (a3.x (a1)).normalize();
   Coord  n3 = (a1.x (a2)).normalize();

   // heights
   double h1 = (n1 ^ a1);
   double h2 = (n2 ^ a2);
   double h3 = (n3 ^ a3);

   return Coord (fabs(h1), fabs(h2), fabs(h3));
}

Coord SxCell::getBoundingBox (double rCut) const
{
   SX_CHECK (volume > 0., volume);
   Coord area;
   for (int dim = 0; dim < 3; ++dim)
      area(dim) = basis((dim+1) % 3).x (basis((dim+2) % 3)).norm ();

   // circumscribing box
   return (rCut/volume) * area;
}


SxList<SymMat> SxCell::getLatticeSymmetries () const
{
   return (getRegularCell ().getLatticeSymmetries (true));
}

bool SxCell::isLatticeSymmetry (const SxMatrix3<Double> &sym) const
{
   SX_CHECK(volume > 0.);
   for (int i = 0; i < 3; ++i)  {
      Coord rotBasisRel = inv ^ (sym ^ basis(i));
      for (int j = 0; j < 3; ++j)
         if (fabs(rotBasisRel(j) - round(rotBasisRel(j))) > epsSym)
            return false;
   }
   return true;
}

SxList<SymMat> SxCell::getLatticeSymmetries (bool regular) const
{
   /* Algorithm: the symmetry matrix in relative coordinates has
      only -1, 0, or 1 as elements if the cell is regular. So we
      try these 3^9 = 19683 matrices. Because the matrix must
      be unitary, many of these matrices can be excluded 
      immediately (leaving 6960 matrices).
   */
   SX_CHECK (regular);
   SxMatrix3<Int> sym;
   SymMat symCart;
   SxList<SymMat> symmetries;
   int i,j;
   bool isUnitary;
   // first col
   for (sym(0,0) = -1; sym(0,0) <= 1; sym(0,0)++)  {
      for (sym(0,1) = -1; sym(0,1) <= 1; sym(0,1)++)  {
         for (sym(0,2) = -1; sym(0,2) <= 1; sym(0,2)++)  {
            // continue if 1st col is 0
            if (sym(0,0) == 0 && sym(0,1) == 0 && sym(0,2) == 0) continue;
            
   // back-shift for readability: second col
   for (sym(1,0) = -1; sym(1,0) <= 1; sym(1,0)++)  {
      for (sym(1,1) = -1; sym(1,1) <= 1; sym(1,1)++)  {
         for (sym(1,2) = -1; sym(1,2) <= 1; sym(1,2)++)  {
            // continue if 2nd col is 0
            if (sym(1,0) == 0 && sym(1,1) == 0 && sym(1,2) == 0) continue;
            
   // back-shift for readability: third col
   for (sym(2,0) = -1; sym(2,0) <= 1; sym(2,0)++)  {
      // continue if 1st row is 0
      if (sym(0,0) == 0 && sym(1,0) == 0 && sym(2,0) == 0) continue;
      for (sym(2,1) = -1; sym(2,1) <= 1; sym(2,1)++)  {
         // continue if 2nd row is 0
         if (sym(0,1) == 0 && sym(1,1) == 0 && sym(2,1) == 0) continue;
         for (sym(2,2) = -1; sym(2,2) <= 1; sym(2,2)++)  {

   // back-shift for readability: matrix is set up
   if (sym(2,2) == 0)  {
      // continue if 3rd col is 0
      if (sym(2,0) == 0 && sym(2,1) == 0) continue;
      // continue if 3rd row is 0
      if (sym(0,2) == 0 && sym(1,2) == 0) continue;
   }
   
   // matrices must have a det of 1 in rel. coordinates
   if (abs(sym.determinant ()) != 1) continue;
   
   symCart = relToCar (SymMat(sym)); // from now on floating point
   
   // determinant must be 1 if symCart is unitary
   if (fabs(fabs(symCart.determinant ()) - 1.) > epsSym) continue;

   // matrix must be unitary in cartesian coords
   isUnitary = true;
   for (i = 0; i < 3; i++)  {
      for (j = 0; j < 3; j++)  {
         if (fabs( (symCart.row(i) ^ symCart.row(j)) 
                     - ((i == j) ? 1. : 0.)            ) > epsSym )  {
            isUnitary = false;
            break; // break j loop
         }
      }
      if (!isUnitary) break; // break i loop
   }
   if (!isUnitary) continue; 

   {
      // Loewdin orthogonalisation (improves upon numerical noise)
      SxMatrix<Complex16> S2 = SxMatrix3<Complex16>(symCart.transpose () ^ symCart), U;
      SxMatrix<Complex16>::Eigensystem eig = S2.eigensystem ();
      SxMatrix<Complex16> diag(3,3);
      diag.set (0.);
      SX_LOOP (d) diag(d,d) = 1./sqrt(eig.vals(d).re);
      U = eig.vecs ^ diag ^ eig.vecs.adjoint ();
      SxMatrix3<Double> U3(U.real ());
      symCart = symCart ^ U3;
   }
   
   // cleanup numerical noise
   for (i = 0; i < 3; ++i)  {
      for (j = 0; j < 3; ++j)  {
         if (fabs(symCart(i)(j)) < 1e-14) 
            symCart(i)(j) = 0.;
         else if (fabs(symCart(i)(j) - 1.) < 1e-14)
            symCart(i)(j) = 1.;
         else if (fabs(symCart(i)(j) + 1.) < 1e-14)
            symCart(i)(j) = -1.;
      }
   }
   SX_CHECK (isLatticeSymmetry (symCart));
         
   // Yeah! The cartesian transformation is unitary and maps the lattice
   // onto itself (the matrix in relative coordinates is integer)
   symmetries << symCart;
            
         } // 2,2
      } // 2,1
   } // 2,0
         } // 1,2
      } // 1,1
   } // 1,0
         } // 0,2
      } // 0,1 
   } // 0,0

   int nSym = int(symmetries.getSize ());

   // --- sorting
   SxArray<SymMat> allSym = symmetries;
   SxArray<SxSymType> type(nSym);
   for (int iSym = 0; iSym < nSym; ++iSym)  {
      type(iSym) = SxRotation(allSym(iSym)).getType ();
   }
   SxArray<ssize_t> order = type.getSortIdx ();
   symmetries.resize (0);
   for (int iSym = 0; iSym < nSym; ++iSym)  {
      symmetries << allSym(order(iSym));
   }

   // --- Checking: are the symmetries a full group?
   SxSymGroup::getSymMulTable (symmetries, epsSym);
   
   return symmetries;
}

SxCell SxCell::getReciprocalCell () const
{
   return SxCell (TWO_PI * inv.transpose (), symGroupPtr);
}
   
int SxCell::moduleDiv3 (const SxVector3<Int> &v, int d)
{
   SX_CHECK (v(0) != 0 || v(1) != 0 || v(2) !=0);
   SxVector3<Int> nv = v;
   // find smallest factor n that makes (n * v) % d == (0,0,0)
   for (int n = 1; n <= d; ++n, nv += v, nv %= d)
      if (nv(0) % d == 0 && nv(1) % d == 0 && nv(2) % d == 0) return n;
   SX_EXIT; return -1;
}

SxVector3<Int> SxCell::computeGen (const SxCell &superCell)
{
   SX_CHECK (volume > 0.);
   // supercell in relative coordinates
   SxMatrix3<Int> relS;
   for (int d=0; d < 3; ++d) 
      relS.setCol(d, round(carToRel(superCell.basis(d))));
   SX_CHECK( ((*this)^relS) == superCell);
   
   int det = abs(relS.determinant ()); 
   SX_CHECK(det != 0); // number of prim cells in supercell
   SX_CHECK(fabs(superCell.volume - det * volume) < 1e-3,
            superCell.volume, det * volume);
   
   // intMap maps relative coordinates of the primitive cell to those of the
   // supercell multiplied by the number of primitive cells. This makes all
   // numbers integer. 
   SxMatrix3<Int> intMap = det * SxMatrix3<Double>(relS).inverse ();
   //cout << "DS^-1 = " << intMap << ", D=" << det << endl;
   //cout << (intMap ^ S) << endl;

   // these will become the new lattice vectors (relative coords)
   SxVector3<Int> gv0(1,0,0), gv1(0,1,0), gv2(0,0,1);

   // The idea is to determine how many times we can go along the different
   // directions before being equivalent to a supercell lattice vector.
   // This number is called rank (for lack of better name) and is
   // determined by moduleDiv3.
   // The rank is then the maximum multiplicity along that direction.

   // r0,r1,r2 will become the mesh (
   int r0 = moduleDiv3(intMap ^ gv0, det),
       r1 = det + 1, // will be minimized
       r2 = det + 1, // will be minimized
       newRank;
   // --- start with basis vector that has largest rank
   if ((newRank = moduleDiv3(intMap ^ gv1, det)) > r0)  {
      r0 = newRank;
      gv1 = gv0;
      gv0 = SxVector3<Int>(0,1,0);
   }
   if ((newRank = moduleDiv3(intMap ^ gv2, det)) > r0)  {
      r0 = newRank;
      gv2 = gv0;
      gv0 = SxVector3<Int>(0,0,1);
   }
   // cout << "r0 = " << r0 << endl;
   
   // --- minimize rank(gv1 + n gv0)
   int nOpt = 0;
   for (int i = 0; i < r0; ++i)  {
      if ((newRank = moduleDiv3(intMap ^ (gv1 + i * gv0), det)) < r1)  {
         nOpt = i;
         r1 = newRank;
         if (r1 == 1) break; // can't become less
      }
   }
   gv1 += nOpt * gv0;
   //cout << "gv1 = " << gv1 << endl << "r1 = " << r1 << endl;

   // --- minimize rank(gv2 + m gv1 + n gv0)
   nOpt = 0; int mOpt = 0;
   for (int i = 0; i < r0; ++i)  {
      for (int j = 0; j < r1; ++j)  {
         if ((newRank = moduleDiv3(intMap^(gv2 + i * gv0 + j * gv1), det)) < r2)
         {
            nOpt = i;
            mOpt = j;
            r2 = newRank;
            // TODO: trust that no further bugs are found and enable shortcuts
            //if (r0 * r1 * r2 == det) break; // done
         }
      }
      //if (r0 * r1 * r2 == det) break; // done
   }
   gv2 += nOpt * gv0 + mOpt * gv1;
   // cout << "gv2 = " << gv2 << endl << "r2 = " << r2 << endl;

   SX_CHECK (r0 * r1 * r2 == det); 

   // reset cell
   set ((*this) ^ gv0, (*this) ^ gv1, (*this) ^ gv2);
   return SxVector3<Int>(r0, r1, r2);
}

ostream &operator<< (ostream &out, const SxCell &cell)
{
   if (cell.volume < 0.)  {
      out << "cell:UNINITIALIZED";
   } else {
      out << "[a1=" << cell.basis(0) << ",a2=" << cell.basis(1) << ",a3=" 
          << cell.basis(2) << "]";
   }
   return out;
}


// -------- Cell Type ---------

SxCell::CellType::CellType ()
{
   type = Unknown;
   a = b = c = -1.;
   bcCos = acCos = abCos = -2.;  // cos can't be -2
}

SxCell::CellType::CellType (const SxCell &cell)
{
   SX_CHECK (cell.volume > 0., cell.volume);

   // --- cell vectors
   const Coord p1 = cell.basis(0);
   const Coord p2 = cell.basis(1);
   const Coord p3 = cell.basis(2);

   // --- norms of the cell vectors
   double p1Norm = p1.norm ();
   double p2Norm = p2.norm ();
   double p3Norm = p3.norm ();

   // --- angles between the cell vectors
   double cosAlpha = (p2 ^ p3) / (p2Norm * p3Norm);
   double cosBeta  = (p3 ^ p1) / (p3Norm * p1Norm);
   double cosGamma = (p1 ^ p2) / (p1Norm * p2Norm);

   // --- find equality between cell vectors
   bool eq12 = (fabs (p1Norm - p2Norm) < 1.e-5);
   bool eq13 = (fabs (p1Norm - p3Norm) < 1.e-5);
   bool eq23 = (fabs (p2Norm - p3Norm) < 1.e-5);

   // --- find special angles
   bool alpha60  = (fabs (cosAlpha - 0.5)   < 1.e-5);  // cos 60 = 1/2
   bool beta60   = (fabs (cosBeta  - 0.5)   < 1.e-5);
   bool gamma60  = (fabs (cosGamma - 0.5)   < 1.e-5);
   bool alpha90  = (fabs (cosAlpha)         < 1.e-5);  // cos 90 = 0
   bool beta90   = (fabs (cosBeta)          < 1.e-5);
   bool gamma90  = (fabs (cosGamma)         < 1.e-5);
   bool alpha109 = (fabs (cosAlpha + 1./3.) < 1.e-5);  // cos 109.47 = -1/3
   bool beta109  = (fabs (cosBeta  + 1./3.) < 1.e-5);
   bool gamma109 = (fabs (cosGamma + 1./3.) < 1.e-5);
//   bool alpha120 = (fabs (cosAlpha + 0.5)   < 1.e-5);  // cos 120 = -1/2
//   bool beta120  = (fabs (cosBeta  + 0.5)   < 1.e-5);
   bool gamma120 = (fabs (cosGamma + 0.5)   < 1.e-5);


   // --- find lattice type
 
   // orthogonal
   if (alpha90 && beta90 && gamma90)  {
      // sc
      if (eq12 && eq13 && eq23)  {
         type = SimpleCubic;
         typeName = "simple cubic (sc)";
         a = b = c = p1Norm;
         bcCos = acCos = abCos = 0.;
      }  else

      // tetragonal
      if (eq12 || eq13 || eq23)  {
         type = Tetragonal;
         typeName = "tetragonal";
         a = p1Norm;
         b = p2Norm;
         c = p3Norm;
         bcCos = acCos = abCos = 0.;
      }  else

      // orthorhombic
      if ( !(eq12 || eq13 || eq23) )  {
         type = Orthorhombic;
         typeName = "orthorhombic";
         a = p1Norm;
         b = p2Norm;
         c = p3Norm;
         bcCos = acCos = abCos = 0.;
      }  else SX_EXIT;
   }  else

   // hexagonal (cell vectors must be given in a special order)
   if (eq12 && alpha90 && beta90 && gamma120)  {
      type = Hexagonal;
      typeName = "hexagonal";
      a = b = p1Norm;
      c = p3Norm;
      bcCos = acCos = 0.;
      abCos = -0.5;
   }  else

   // cubic
   if (eq12 && eq23 && eq13)  {
      // fcc
      if (alpha60 && beta60 && gamma60)  {
         type = FaceCenteredCubic;
         typeName = "face centered cubic (fcc)";
         a = cbrt (4. * cell.volume);  // V = a^3 / 4
         b = c = -1.;
         bcCos = acCos = abCos = 0.5;
      }  else

      // bcc
      if (alpha109 && beta109 && gamma109)  {
         type = BodyCenteredCubic;
         typeName = "body centered cubic (bcc)";
         a = cbrt (2. * cell.volume);  // V = a^3 / 2
         b = c = -1.;
         bcCos = acCos = abCos = -1./3.;
      }
   }  else

   // unknown lattice type
   {
      type = Unknown;
      typeName = "unknown";
      a = b = c = -1.;
      bcCos = acCos = abCos = -2.;
   }
}

SxCell::CellType::Type SxCell::CellType::getType () const
{
   return type;
}

SxString SxCell::CellType::getTypeName () const
{
   return typeName;
}

double SxCell::CellType::getA () const
{
   SX_CHECK (a > 0., a);
   return a;
}

double SxCell::CellType::getB () const
{
   SX_CHECK (b > 0., b);
   return b;
}

double SxCell::CellType::getC () const
{
   SX_CHECK (c > 0., c);
   return c;
}

double SxCell::CellType::getCosAB () const
{
   SX_CHECK (abCos > -2., abCos);
   return abCos;
}

double SxCell::CellType::getCosBC () const
{
   SX_CHECK (bcCos > -2., bcCos);
   return bcCos;
}

double SxCell::CellType::getCosAC () const
{
   SX_CHECK (acCos > -2., acCos);
   return acCos;
}

void SxCell::CellType::print () const
{
   cout << SX_SEPARATOR;
   cout << "| Lattice structure: " << typeName << endl;
   cout << "|" << endl;
   cout << "|   aLat = " << a << endl;
   if (b > 0.)  cout << "|   bLat = " << b << endl;
   if (c > 0.)  cout << "|   cLat = " << c << endl;
   cout << "|" << endl;
   cout << "| Angles between the unit cell vectors a1, a2, a3 in degrees:"
        << endl;
   cout << "|" << endl;
   cout << "|   <)(a1,a2) = " << acos (abCos) * RAD2DEG << " deg" << endl;
   cout << "|   <)(a2,a3) = " << acos (bcCos) * RAD2DEG << " deg" << endl;
   cout << "|   <)(a3,a1) = " << acos (acCos) * RAD2DEG << " deg" << endl;
   cout << SX_SEPARATOR;
}


SxMatrix3<Int> SxCell::getConventional () const
{
   SxArray<SymMat> latSym = getLatticeSymmetries ();
   if (latSym.getSize () == 48)  {

   } else {
      cout << "WARNING: in getConventional: failed to detect cubic cell type."
           << endl;
      cout << "Cell is " << (*this) << endl;
      cout << "with " << latSym.getSize () << " symmetries." << endl;
      SX_QUIT;
   }
   SxList<Coord> aCub; 
   for (int iSym = 0; iSym < latSym.getSize (); ++iSym)  {
      SxSymType type = SxRotation::getType(latSym(iSym));
      if (type.axisCount == 4)  {
         bool newAxis = true;
         Coord symAxis = type.opCoord;
         for (int i = 0; i < aCub.getSize (); ++i)  {
            if (fabs (symAxis ^ aCub(i)) > 1e-10 ) {
               newAxis = false;
            }
         }
         if (newAxis) aCub << symAxis;
      }
   }
   SX_CHECK (aCub.getSize () == 3, aCub.getSize ());
   if (fabs(aCub(0)(0)) < fabs(aCub(1)(0)))  {
      Coord ax = aCub(1);
      aCub(1) = aCub(0);
      aCub(0) = ax;
   }
   if (fabs(aCub(0)(0)) < fabs(aCub(2)(0)))  {
      Coord ax = aCub(2);
      aCub(2) = aCub(0);
      aCub(0) = ax;
   }
   if (fabs(aCub(1)(1)) < fabs(aCub(2)(1)))  {
      Coord ay = aCub(2);
      aCub(2) = aCub(1);
      aCub(1) = ay;
   }
   for (int i = 0; i < 3; ++i)  {
      if (aCub(i)(i) < 0.) aCub(i) = -aCub(i);
   }
   for (int mul = 1; mul <= 4; mul *=2)  {
      double aLat = cbrt (mul * volume);
      SxMatrix3<Int> res;
      bool ok = true;
      for (int i = 0; i < 3; ++i)  {
         Coord aRed = carToRel (aLat * aCub(i));
         if ((aRed - round(aRed)).normSqr () < 1e-10)  {
            for (int j = 0; j < 3; ++j) res(j,i) = int(lround(aRed(j)));
         } else {
            ok = false;
         }
      }
      if (ok) return res;
   }
   cout << "Failed to deduce conventional cubic cell" << endl;
   cout << "Cell = " << (*this) << endl;
   SX_EXIT;
}

SxMatrix3<Double> SxCell::getStandardRot () const
{
   CellType type(*this);
   CellMat newCell;
   double a = type.getA (), b = type.getB (), c = type.getC ();
   switch (type.getType ())  {
      case CellType::SimpleCubic:
         newCell = CellMat (a,0,0, 0,a,0, 0,0,a); break;
      case CellType::Tetragonal:
      case CellType::Orthorhombic:
         newCell = CellMat (a,0,0, 0,b,0, 0,0,c); break;
      case CellType::FaceCenteredCubic:
         newCell = 0.5 * CellMat (0,a,a, a,0,a, a,a,0); break;
      case CellType::BodyCenteredCubic:
         newCell = 0.5 * CellMat (-a,a,a, a,-a,a, a,a,-a); break;
      case CellType::Hexagonal:
         newCell = CellMat (a, -0.5      *a,0,
                            0, sqrt(0.75)*a,0,
                            0,      0,       c); break;
      case CellType::Unknown:
         return SxMatrix3<Double> (1.,0.,0.,0.,1.,0.,0.,0.,1.);
      default:
         SX_EXIT;
   }
   return newCell ^ inv;
}
