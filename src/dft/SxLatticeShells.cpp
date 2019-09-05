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

#include <SxLatticeShells.h>
#include <SxSortedList.h>

SxLatticeShells::SxLatticeShells ()
{
   // empty
}

SxLatticeShells::SxLatticeShells (const SxCell &cell_, int nShells_,
                                  const Coord &referenceVec_)
{
   cell          = cell_;
   nShells       = nShells_;
   referenceVec  = referenceVec_;

   shells.resize (nShells+1);  // i.e., for 0 <= iShell <= nShells
}

void SxLatticeShells::find ()
{
   cout << SX_SEPARATOR;
   cout << "| Finding shells of neighboured lattice vectors ... ";
   cout.flush ();

   int    x, y, z, n = nShells;
   Coord  R, xyz;
   double length2;

   SxSortedList<double> lengthes2;
   SxList<Coord >       vectors;
   int                  pos;

   SxList<Coord>        allVecsList;
   SxList<int>          allShellsList;

   // --- using inversion symmetry of the lattice
   for (x = 0; x <= n; x++)  {
      for (y = -n; y <= n; y++)  {
         for (z = -n; z <= n; z++)  {
            xyz = Coord ((double)x, (double)y, (double)z);
            R = cell ^ xyz;
            length2 = R.normSqr();
            pos = (int)lengthes2.append (length2);
            vectors.insert (pos, R);
            if (x > 0)  {
               pos = (int)lengthes2.append (length2);
               vectors.insert (pos, -R);
            }
         }
      }
   }

   SX_CHECK (vectors.getSize() == 8*n*n*n + 12*n*n + 6*n + 1,
             vectors.getSize(), 8*n*n*n + 12*n*n + 6*n + 1);
   SX_CHECK (vectors.getSize() == lengthes2.getSize(),
             vectors.getSize(), lengthes2.getSize());

   // --- create arrays
   double eps = 1.e-4;
   double h = 0.;
   for (int i = 0, iShell = 0; i < lengthes2.getSize(); i++)  {
      R = vectors(i) + referenceVec;

      // --- distinguish the shells
      if (fabs(lengthes2(i) - h) < eps)  {
         shells(iShell).append (R);
      }  else  {
         h = lengthes2(i);
         iShell++;
         if (iShell > nShells)  break;
         shells(iShell).append (R);
      }

      allVecsList.append (R);
      allShellsList.append(iShell);
   }
   
   // --- creating arrays from lists
   allNeighbourVecs = allVecsList;
   allShells = allShellsList;

   // --- dump results
   cout << "done." << endl;  cout.flush ();
   cout << "|" << endl;

   sxprintf ("| %d shells of neighboured lattice vectors considered.\n",
             nShells);
   cout  <<  "| reference vector: " << referenceVec << endl;
   cout  <<  "|" << endl;

   for (int iShell = 0; iShell <= nShells; iShell++)
      sxprintf ("| shell %d --> %2d neighboured lattice vectors found.\n",
                iShell, (int)shells(iShell).getSize());
   sxprintf ("|   total --> %2d neighboured lattice vectors found.\n",
                (int)allNeighbourVecs.getSize());

   for (int iShell = 0; iShell <= nShells; iShell++)  {
      cout << SX_SEPARATOR;
      sxprintf ("| neighbour vectors in shell %d:\n", iShell);
      cout << SX_SEPARATOR;
      printShell (shells(iShell));
      cout << endl;
   }
   cout << SX_SEPARATOR;
   cout << "| all neighbour vectors:" << endl;
   cout << SX_SEPARATOR;
   printShell (allNeighbourVecs);
   cout << endl;
}

SxArray<Coord> SxLatticeShells::getShell (int i) const
{
   SX_CHECK (i >= 0 && i <= nShells, i, nShells);

   SxArray<Coord> shell = shells(i);  // SxList --> SxArray
   return shell;
}

SxArray<Coord> SxLatticeShells::getAllNeighbVecs () const
{
   SX_CHECK (allNeighbourVecs.getSize() > 0, allNeighbourVecs.getSize());
   return allNeighbourVecs;
}

Coord SxLatticeShells::getNeighbVec (int idx) const
{
   SX_CHECK (idx >= 0 && idx <= allNeighbourVecs.getSize(),
             idx, allNeighbourVecs.getSize());

   return allNeighbourVecs(idx);
}

SxArray<int> SxLatticeShells::getAllShells () const
{
   SX_CHECK (allShells.getSize() > 0, allShells.getSize());
   return allShells;
}

int SxLatticeShells::getShellOfIdx (int idx) const
{
   SX_CHECK (idx >= 0 && idx <= allShells.getSize(),
             idx, allShells.getSize());

   return allShells(idx);
}

int SxLatticeShells::getNAllVecs () const
{
   return (int)allNeighbourVecs.getSize ();
}

void SxLatticeShells::printShell (const SxList<Coord> &list) const
{
   SxArray<Coord> array = list;
   printShell (array);
}

void SxLatticeShells::printShell (const SxArray<Coord> &array) const
{
   int     i, j, n = (int)array.getSize ();
   double  a;
   Coord   vec;

   sxprintf ("{");
   for (i = 0; i < n; i++)  {
      vec = array(i);

      a = vec(0);
      sxprintf ("(%3.2f", a);
      for (j = 1; j < 3; j++)  {
         a = vec(j);
         sxprintf (", %3.2f", a);
      }
      if (i < n-1)  { sxprintf ("), "); }  else  { sxprintf (")"); }
   }
   sxprintf ("}");
}
