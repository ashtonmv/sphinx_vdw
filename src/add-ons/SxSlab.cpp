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
#include <SxConfig.h>
#include <SxSlab.h>
#include <SxSortedList.h>
#include <SxSymFinder.h>
#include <SxAtomicStructure.h>
#include <SxNeighbors.h>
#include <SxSpeciesData.h>

#ifndef SX_STANDALONE

#define DEBUG_PRINT_VAR(var) \
do  {\
cout << "--- DEBUG (" << __FILE__ << ":" << __LINE__ << ")" << endl; \
cout << #var << " = " << (var) << endl;\
} while (0)
#define DEBUG_PRINT(message) \
do  {\
cout << "--- DEBUG (" << __FILE__ << ":" << __LINE__ << ")" << endl; \
cout << #message << endl; \
} while (0)
#define DEBUG_STATEMENT(code) \
do  {\
   cout << "--- DEBUG (" << __FILE__ << ":" << __LINE__ << ")" << endl; \
   code; \
} while (0)

// ----------------------- SxSlab --------------------------------
SxSlab::SxSlab ()
   : hkl (0),
     axisRotation (1., 0., 0., 0., 1., 0., 0., 0., 1.)
{
   // empty
}


void SxSlab::setBulk (const SxAtomicStructure &strIn)
{
   bulkStructure = strIn;
   hkl.set (0); // reset plane
}


void SxSlab::setPlane (const RelVec& miller)
{
   SX_CHECK (miller.absSqr().sum() > 0);
   hkl = miller;
}

void SxSlab::setPlane (int h, int k, int l)
{
   SX_CHECK (h != 0 || k != 0 || l != 0);
   hkl(0) = h; hkl(1) = k; hkl(2) = l;
}

int SxSlab::getGcd (int x, int y)
{
   x = abs(x); y = abs(y);
   if (x > y) return getGcd (y, x);
   if (x == 0) return y;
   int mod = y % x;
   if (mod == 0) return x;
   return getGcd (mod, x); 
}

bool SxSlab::shortenVector (SxVector3<Int>* vPtr)
{
   SX_CHECK (vPtr);
   SxVector3<Int> &v = *vPtr;
   SX_CHECK (v.absSqr().sum() != 0);
   int gcd = getGcd(abs(v(1)),abs(v(2)));
   gcd = getGcd( abs(v(0)), gcd);
   for (int i = 0; i < 3; i++) v(i) /= gcd;
   return (gcd != 1);
}

void SxSlab::integerGramSchmidt(RelVec *v1Ptr, 
                                RelVec *v2Ptr,
                                const CellMat &cell)
{
   SX_CHECK (v1Ptr);
   SX_CHECK (v2Ptr);
   RelVec &v1 = *v1Ptr;
   RelVec &v2 = *v2Ptr;
   bool change;
   int overlap;
   CellMat metric = cell.transpose () ^ cell;
   int count = 0; // hexagonal lattice would loop forever
   do  {
      change = false;
      // --- make 0 orthogonal to 1
      overlap = int ( 2. * (v2 ^ (metric ^ v1) / (cell ^ v2).absSqr ().sum ()));
      overlap = overlap / 2 + overlap % 2; // round
      change |= (overlap != 0);
      v1 -= v2 * overlap;
      change |= shortenVector(v1Ptr);

      // --- make 1 orthogonal to 0
      overlap = int ( 2. * (v2 ^ (metric ^ v1) / (cell ^ v1).absSqr ().sum ()));
      overlap = overlap / 2 + overlap % 2; // round
      change |= (overlap != 0);
      v2 -= v1 * overlap;
      change |= shortenVector(v2Ptr);
   } while (change && count++ < 400/* just a number */);
}

/// Find a cell where the first two basis vectors are orthogonal to hkl
/// and call findElementaryCell afterwards
void SxSlab::findCell ()
{
   SxList<RelVec> bList;
   // three vectors orthogonal to (h,k,l)
   bList << RelVec (0, -hkl(2), hkl(1));
   bList << RelVec (-hkl(2), 0, hkl(0));
   bList << RelVec (-hkl(1), hkl(0), 0);

   // --- get 2 in-plane vectors
   
   for (int i = 0, idx = 0; i < 3; i++)
      // --- remove zero vectors, shorten others
      if (bList(idx).absSqr ().sum () == 0)
         bList.remove(idx);
      else  {
         shortenVector(&bList(idx));
         idx++;
      }
   // --- avoid colinear vectors 0 and 1
   while (bList.getSize () > 1 && bList(0).x(bList(1)).absSqr ().sum () == 0) 
      bList.remove(1);
   if (bList.getSize () < 2)  {
      cout << "Can't find linear independant vectors in net plane..." << endl;
      SX_EXIT;
   }
   // replace bList with list of 0 and 1 (resize() doesn't work due to bug)
   bList = (SxList<RelVec> () << bList(0) << bList(1));
   // make angle between 0 and 1 as large as possible
   integerGramSchmidt(&bList(0), &bList(1), bulkStructure.cell);

   // --- get non-inplane vector b3
   int idx = 0;
   while (hkl(idx) == 0) idx++; // because b3 ^ normal = normal(idx) != 0
   RelVec b3 (0); // (0,0,0)
   b3(idx) = 1; // set one element to 1, i.e. (1,0,0) (0,1,0) or (0,0,1)
   bList << b3;

   SxMatrix3<Int> relCell;
   relCell = SxMatrix3<Int> (bList(0), bList(1), bList(2)).transpose ();
   // axis orientation such that determinant is positive
   if (relCell.determinant () < 0)
      relCell = SxMatrix3<Int> (bList(1), bList(0), bList(2)).transpose ();

   /*
   DEBUG_PRINT (relative cell before reduction);
   DEBUG_STATEMENT (relCell.print (););
   */

   findElementaryCell (relCell);
}

/// Store shortest vector unless 0.
void SxSlab::storeShortest(const SxVector3<Double> &newVec,
                                 SxVector3<Double> *minVecPtr)
{
   SX_CHECK (minVecPtr);
   double length = minVecPtr->absSqr ().sum ();
   if (length < 1e-8 || newVec.absSqr ().sum () < length) 
      *minVecPtr = newVec;
}

/// Get elementary cell
void SxSlab::findElementaryCell (const SxMatrix3<Int> &relCell)
{
   // lattice points in cell
   SxList<RelVec> latticePoints = getLatticePointsInCell (relCell);

   SxList<RelVec>::Iterator relPointIt;
   SxList<SxVector3<Double> > netplanePoints;
   SxList<RelVec> relNetplanePoints;
   SxVector3<Double> cartCoord, c1(0.), c2, c3(0.);
   int dist, mindist = -1; // distance from netplane 0
      
   // --- loop over lattice points in order to find
   // * c1 in the netplane through origin
   // * netplane through origin (for c2 later)
   // * c3 to the closest netplane
   for (relPointIt = latticePoints.begin ();
        relPointIt != latticePoints.end ();
        relPointIt++)
   {
      if ((*relPointIt).absSqr ().sum () == 0) continue; // ignore (0,0,0)
      cartCoord = bulkStructure.cell.relToCar(*relPointIt);
      dist = abs(*relPointIt ^ hkl);
      if (dist == 0)  {
         // --- c1 candidates
         netplanePoints << cartCoord; 
         relNetplanePoints << *relPointIt;
         storeShortest (cartCoord, &c1);
      } else {
         // --- c3 candidates
         if (dist < mindist || mindist < 0)  {
            // get closest non-zero netplane
            mindist = dist;
            c3.set(0.);
         }
         if (dist == mindist) storeShortest(cartCoord, &c3);
      }
   } // loop over lattice points

   /*
   DEBUG_PRINT_VAR (netplanePoints.getSize ());
   DEBUG_PRINT_VAR (relNetplanePoints);
   */

   // --- find c2
   SxList<SxVector3<Double> >::Iterator netVecIt;
   int relShift = 0;
   double areasq, minArea = -1.;
   double cosine = 0., shiftedCosine, minCosine = -1.;
   bool equal;
   for (netVecIt = netplanePoints.begin (),
        relPointIt = relNetplanePoints.begin ();
        netVecIt != netplanePoints.end ();
        netVecIt++, relPointIt++)
   {
      areasq = ((*netVecIt).x (c1)).absSqr ().sum (); // square of cell area
      if (areasq < 1e-8) continue; // colinear with c1
      equal = fabs (areasq - minArea) < 1e-8;
      if ( minArea < 0. || (areasq < minArea && !equal))  {
         // minimal area
         minArea = areasq;
         minCosine = -1.;
         equal = true;
      }
      // continue unless we are in closest row
      if (!equal) continue; 

      cosine = *netVecIt ^ c1;
      relShift = 0;
      // --- look for shorter vectors in neighbour cells
      shiftedCosine = (*netVecIt - c1) ^ c1;
      equal = fabs(fabs(shiftedCosine) - fabs(cosine)) < 1e-8;
      if (   (!equal && fabs(shiftedCosine) < fabs(cosine))
          || ( equal && (shiftedCosine < cosine)))  {
         cosine = shiftedCosine;
         relShift = -1;
      }
      shiftedCosine = (*netVecIt + c1) ^ c1;
      equal = fabs(fabs(shiftedCosine) - fabs(cosine)) < 1e-8;
      if (   (!equal && fabs(shiftedCosine) < fabs(cosine))
          || ( equal && (shiftedCosine < cosine)))  {
          cosine = shiftedCosine;
          relShift = +1;
      }

      equal = fabs(fabs(cosine) - minCosine) < 1e-8;
      if (   minCosine < 0.                         // first one in this row
            || (fabs(cosine) < minCosine && !equal) // smallest one
            || (equal && cosine < 0.))              // larger angle
      {
         minCosine = fabs(cosine);
         c2    = *netVecIt   + double(relShift) * c1;
      }
   } // loop over netplane points
   newCell = CellMat (c1, c2, c3).transpose ();
   if (newCell.determinant () < 0.)
      newCell = CellMat (c2, c1, c3).transpose ();
}

void SxSlab::rotateCell ()
{
   // rotate cell such that c1 || x and z || normal
   SxVector3<Double> xAxis, yAxis, zAxis;
   // zAxis  = newCell.col(0).x(newCell.col(1));
   zAxis = bulkStructure.cell.inv.transpose () ^ hkl;
   zAxis.normalize ();
   xAxis  = newCell.col(0); xAxis.normalize ();
   yAxis  = xAxis.x (zAxis);
   axisRotation = CellMat (xAxis, yAxis, zAxis).transpose ().inverse (); 
   newCell =  axisRotation ^ newCell; 
   if (newCell(2,2) < 0.)  {
      // we want the z component of the 3rd axis to be positive
      newCell = SxCell(newCell.col(1), newCell.col(0), -newCell.col(2));
   }
   if (newCell.determinant () < 0.)  {
      // --- positive determinant
      newCell = SxCell(newCell.col(0), -newCell.col(1), newCell.col(2));
   }
   // DEBUG_PRINT_VAR (newCell);
}

void SxSlab::getSlab (int nLayer, int nVacLayer, bool zOrthogonal)
{
   SxVector3<Double> relCoord;
   
   // get coordinates for n-layered, rotated bulk system
   slabStructure = axisRotation ^ bulkStructure;
   slabStructure.cell = newCell;
   slabStructure %= newCell; // important in order to get a compact slab
   slabStructure = slabStructure.repeat (1, 1, nLayer);
   SxPtr<SxAtomInfo> aInfo = SxAtomInfo::derive (bulkStructure.atomInfo);
   aInfo->parentMap.resize (slabStructure.getNAtoms ());
   for (int ia = 0; ia < slabStructure.getNAtoms (); ++ia)
      aInfo->parentMap(ia) = ia / nLayer;
   aInfo->nAtoms = slabStructure.atomInfo->nAtoms;
   aInfo->setupOffset ();
   aInfo->meta = bulkStructure.atomInfo->meta;
   slabStructure.replaceInfo (aInfo);
   
   // --- change z direction of cell
   if (zOrthogonal)  {
      slabStructure.cell(0,2) = 0.; 
      slabStructure.cell(1,2) = 0.; 
      slabStructure.cell(2,2) = newCell(2,2) * (nLayer + nVacLayer);
      slabStructure.cell.setup ();
   } else  {
      Coord c3 = newCell.col(2) * double(nLayer + nVacLayer);
      // --- minimize in-plane component of c3
      double c3z = c3(2);
      c3(2) = 0.;
      newCell.map (&c3, SxCell::Origin);
      c3(2) = c3z;
      // --- set up new cell
      slabStructure.cell = SxCell(newCell.col (0), newCell.col(1), c3);
   }
   
   // --- move atoms into new cell and sort them by z, y, x
   slabStructure %= slabStructure.cell;
   slabStructure.sort (true);
}

void SxSlab::getSaturatedSlab (int nLayer, int nVacLayer,
                               bool zOrthogonal,
                               Direction dir,
                               double bondThreshold,
                               double newBondLength,
                               double identicalThreshold)
{
   SX_CHECK (bondThreshold > 0., bondThreshold);
   SX_CHECK (newBondLength > 0., newBondLength);
   SX_CHECK (identicalThreshold >= 0., identicalThreshold);

   SxList<Coord> saturationAtoms;

   // get the slab
   getSlab (nLayer, nVacLayer, zOrthogonal);

   SxGrid grid(bulkStructure, 3), slabGrid(slabStructure, 3);
   SxNeighbors nn;
   SxVector<Int> map
      = slabStructure.atomInfo->getIdxMap (bulkStructure.atomInfo);
   Coord vacCenter = slabStructure.cell.basis(2)
                   * (nLayer + 0.5 * nVacLayer) / (nLayer + nVacLayer);
   for (int iaSlab = 0; iaSlab < slabStructure.getNAtoms (); ++iaSlab)  {
      int ia = map(iaSlab);
      nn.compute (grid, bulkStructure, bulkStructure.getAtom (ia),
                  bondThreshold, SxNeighbors::StoreAbs | SxNeighbors::StoreRel);
      for (int in = 0; in < nn.getSize (); ++in)  {
         // bottom layer saturation
         Coord nnPos = (axisRotation ^ nn.relPositions(in))
                     + slabStructure(iaSlab);
         if (slabStructure.find (nnPos, slabGrid) == -1)
         {
            double relPos = slabStructure.cell.carToRel (nnPos - vacCenter)(2);
            while (relPos < -0.5) relPos += 1.;

            if (dir == Both
                || (dir == Downward && relPos > 0.)
                || (dir == Upward && relPos < 0.))
            {
               Coord satAtom = slabStructure.getAtom (iaSlab)
                             + (axisRotation ^ nn.relPositions(in))
                               * (newBondLength / nn.relPositions(in).norm ());
               saturationAtoms << slabStructure.cell.getMapped(satAtom);
            }
         }
      }
   }

   // --- throw out double saturation atoms
   SxList<Coord>::Iterator atom1, atom2;
   int idx1, idx2, nIdentical;
   Coord shift, delta;
   SxList<int> toBeRemoved;
   double identicalSquare = identicalThreshold * identicalThreshold;
   for (atom1 = saturationAtoms.begin (), idx1 = 0;
        idx1 < saturationAtoms.getSize ();
        ++atom1, idx1++)  {
      SX_CHECK(toBeRemoved.getSize () == 0, toBeRemoved.getSize ());
      nIdentical = 1;
      shift.set(0.);
      for (atom2 = atom1, ++atom2, idx2 = idx1 + 1;
            idx2 < saturationAtoms.getSize ();
            ++atom2, idx2++)  {
         delta = *atom2 - *atom1;
         // map difference vector into -0.5 to +0.5 range (relative)
         bulkStructure.cell.map(&delta, SxCell::Origin);
         if (delta.normSqr () <= identicalSquare)  {
            shift += delta;
            // propably iterators don't like a removal of the element they
            // refer to, so we delay the removal
            toBeRemoved << idx2; 
            nIdentical++; 
         }
      }
      *atom1 += shift / double (nIdentical);
      while (toBeRemoved.getSize () > 0)  {
         // very important: work from end, or we would shift the indices
         saturationAtoms.remove(toBeRemoved.last ());
         toBeRemoved.removeLast ();
      }
      // workaround for iterator problems
      if (idx1 == 0) 
         atom1 = saturationAtoms.begin ();
      else
         (--atom1)++; // go back and forth in case the next element was removed
   }
   // --- Add the saturation atoms
   slabStructure.startCreation ();
   slabStructure.newSpecies ();

   for (atom1 = saturationAtoms.begin ();
        atom1 != saturationAtoms.end ();
        ++atom1)  {
      slabStructure.addAtom (slabStructure.cell.getMapped(*atom1));
   }
   slabStructure.endCreation ();
   slabStructure.sort ();
}

SxMatrix3<Int> SxSlab::intInverse (const SxMatrix3<Int> m, int *factor)
{
   SX_CHECK (factor);
   SxMatrix3<Int> inv(m(1).x (m(2)), m(2).x (m(0)), m(0).x (m(1)));
   int i, j, gcd = 0;
   for (i = 0; i < 3; i++)  {
      for (j = 0; j < 3; j++)  {
         if (inv(i,j) != 0)  {
            if (gcd != 0)  {
               gcd = getGcd (inv(i,j), gcd);
            }  else  {
               gcd = inv(i,j);
            }
         }
      }
   }
   SX_CHECK (gcd > 0, gcd);
   for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++)
         inv(i,j) /= gcd;
   *factor = (m.determinant () / gcd);
   return inv;
}

/// fill cell with lattice points
SxList<RelVec> SxSlab::getLatticePointsInCell (const SxMatrix3<Int> &relCell)
{
   RelVec minComp(0), maxComp(0), x, point;
   int i, j, boundary, inside;
   SxList<RelVec> latticePoints;
   for (i = 0; i < 3; i++)    // directions 
      for (j = 0; j < 3; j++) // vectors bj
      {
         if (minComp(i) > relCell(i,j))
            minComp(i) = relCell(i,j);
         if (maxComp(i) < relCell(i,j))
            maxComp(i) = relCell(i,j);
      }
   // determinant of relCell must be positive, because the relative
   // coordinates are expected to be between 0 and boundary, and it is
   // boundary = determinant / abs(gcd of inverse of matrix elements)
   SX_CHECK (relCell.determinant () > 0, relCell.determinant ());
   // integer part of inverse, i.e. boundary * relCell.inverse ();
   SxMatrix3<Int> inv = intInverse(relCell, &boundary);
   for (x(0) = minComp(0); x(0) <= maxComp(0); x(0)++)
      for (x(1) = minComp(1); x(1) <= maxComp(1); x(1)++)
         for (x(2) = minComp(2); x(2) <= maxComp(2); x(2)++)  {
            // get point in new relative coordinates (multiplied by boundary)
            point = inv ^ x;
            // check that it is inside boundaries
            inside = true;
            for (i = 0; i < 3; i++)
               inside &= (point(i) <= boundary) && (point(i) >= 0);
            if (inside)
               latticePoints << x; // store in old relative coordinates
         }
   return latticePoints;
}         

SxAtomicStructure
SxSlab::compute (Direction dir,
                      int nLayer,
                      int nVacLayer,
                      bool zOrthogonal,
                      double threshold,
                      double bondLength,
                      double minDist,
                      const SxString & chemSymbol)
{
   findCell ();
   rotateCell ();
   
   // --- get new atomic coordinates
   if (dir == None)  {
      getSlab (nLayer, nVacLayer, zOrthogonal);
   } else {
      getSaturatedSlab (nLayer, nVacLayer, zOrthogonal, dir, threshold,
                             bondLength, minDist);
      const SxArray<SxString> &oldNames = bulkStructure.getElements ();
      int nSpecies = bulkStructure.getNSpecies ();
      SxArray<SxString> chemName(nSpecies + 1);
      for (int is = 0; is < nSpecies; ++is) chemName(is) = oldNames(is);
      chemName(nSpecies) = chemSymbol;
      slabStructure.atomInfo->meta.attach (SxAtomicStructure::Elements,
                                           chemName);
   }

   SxAtomicStructure as (slabStructure, SxAtomicStructure::Copy);
   // move atoms to high-symmetry position
   as += SxSymFinder(as).getHighSymShift ();
   // --- try to put slab in the (vertical) center of the cell
   if (zOrthogonal)  {
      if (as.constRef(0)(2) < 0.) as += 0.5 * as.cell.basis(2);
      if (as.constRef(0)(2) > as.cell(2,2)) as -= 0.5 * as.cell.basis(2);
   }
   // move atoms back into cell
   as %= as.cell;

   return as;

}
#else  /* SX_STANDALONE */

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   SxCLI::CliArg *opt;
   cli.preUsageMessage =
      "This add-on serves to set up slab geometries from a known bulk geometry."
      " In order to do so, an adapted bulk cell is calculated where two basis "
      "vectors lie in the desired plane (given by the Miller indices). This "
      "cell is multiplied by the number of layers (-l option). If additional "
      "vacuum layers are specified (-v option), empty space is inserted.\n"
      "There's also a (still experimental) automatic saturation code.";
   SxString inFile = cli.option ("-i", "file", "S/PHI/nX input for bulk")
                     .toString("input.sx");
   SxString outFile 
      = cli.option ("-o", "file name", 
                    "write structure to output file instead of screen "
                    "(or STDOUT)")
        .toString ("");

   int nLayer      = cli.option ("-l", "number", "number of layers")
                     .toInt (1,1);
   int nVacLayer   = cli.option ("-v", "number", "number of vacuum layers")
                     .toInt (0,0);
   SxList<SxString> hklList
      = cli.option ("--hkl", "plane", "Miller indices. Enter in one, e.g. "
                    "'--hkl 1-10' or give one option for each idx, e.g. "
                    "'--hkl 1 --hkl -1 --hkl 0'. Note that these indices "
                    "refer to the given cell, even for fcc and alike, where "
                    "(h k l) traditionally refer to the cubic (super)cell.")
        .required ().toList (); 
   if (hklList.getSize () == 1)  {
      // --- "hkl" -> ("h","k","l")
      int i, idx = 0;
      for (i = 0; i < 3; i++, idx++)  {
         if (hklList(0).getSize () == idx)  {
            cout << "Invalid hkl expression '" << hklList(0);
            cout << "': missing " << SxString("hkl")(i) << "." << endl;
            cli.setError ();
            break;
         }
         if (hklList(0)(idx) == '-')  {
            if (++idx != hklList(0).getSize ())
               hklList << hklList(0).subString (idx - 1, idx);
            else  {
               cout << "hkl error: no number after '-' sign" << endl;
               cli.setError ();
               break;
            } 
         } else  {
            hklList << hklList(0).subString (idx, idx);
         }
      }
      if (idx < hklList(0).getSize () && hklList.getSize () == 4)  {
         cout << "Invalid hkl expression '" << hklList(0);
         cout << "': trailing characters after reading" << endl;
         cout << "h = " << hklList(1) << endl;
         cout << "k = " << hklList(2) << endl;
         cout << "l = " << hklList(3) << endl;
         cli.setError ();
      }
      hklList.remove (0);
   }

   SxList<double> glide = cli.option ("--glide", "2 vector", 
       "2D glide vector (relative). The c lattice vector will be shifted by "
       "the indicated fractions of the surface vectors.")
                          .toDoubleList ();
   if (glide.getSize () > 0 && glide.getSize () != 2)
   {
      cout << "Invalid glide vector with " << glide.getSize () << " elements";
      cli.setError ();
   }

   bool conventionalCell
      = cli.option ("--conventional", "hkl refer to conventional cubic cell")
        .toBool ();

   bool zOrthogonal 
      = cli.option ("--orthogonal", "set lateral components of c3 axis to 0")
        .toBool ();

   opt = &cli.option 
      ("--saturate", "code",
       "request automatic saturation: the code is\n"
       "'{UDB}:<maxDist>:<bondLength>[:<symbol>]'\n"
       "The first letter tells which surface to saturate: (U)pward, (D)ownward "
       "or (B)oth surfaces.\n"
       "<maxDist> (Bohr) is the threshold for the bulk bond distance. All bonds"
       " shorter than this threshold and broken at the surface are "
       "saturated.\n"
       "<bondLength> (Bohr) is the bonding length to the saturation atom, the "
       "direction of the bond is that of the broken bond.\n"
       "<symbol> is an optional chemical symbol, the default is H (hydrogen).");
   opt->required (false);
   SxSlab::Direction dir = SxSlab::None;
   double threshold = 0., bondLength = 0.;
   SxString chemSymbol;
   if (!cli.error && opt->exists (true))  {
      SxString satString = opt->toString ();
      SxList<SxString> tokens = satString.tokenize (':');
      if (tokens.getSize () < 3 || tokens.getSize () > 4)  {
         cout << "Can't understand saturation code: too little or too many ':'";
         cout << endl;
         cli.setError ();
      } else {
         if (tokens(0).toUpper () == "U") dir = SxSlab::Upward;
         if (tokens(0).toUpper () == "D") dir = SxSlab::Downward;
         if (tokens(0).toUpper () == "B") dir = SxSlab::Both;
         if (dir == SxSlab::None)  {
            cout << "Error: don't understand side code '" << tokens(0);
            cout << "'. Must be U, D or B." << endl;
            cli.setError ();
         }
         threshold = tokens(1).toDouble ();
         bondLength = tokens(2).toDouble ();
         chemSymbol = (tokens.getSize () > 3) ? tokens(3) : SxString("H");
      }
   }
   opt = &cli.option ("--mindist", "number", 
         "Only for saturations: minimal distance between saturation atoms."
         "If atoms are closer, their coordinates will be averaged. This may"
         "change the bond distance and the bond direction. You'll have to "
         "play around with <bondLength> to get the desired result - sorry, I'm "
         "lazy.");
   if (opt->exists (true) && dir == SxSlab::None)  {
      cout << "--mindist given without (valid) --saturation option" << endl;
      cli.setError ();
   }
   double minDist = opt->toDouble (0.01, 0.);

   cli.finalize ();

   // parse hkl-List
   if (hklList.getSize () != 3)  {
      cout << "Corrupt hkl-list. Check command line options." << endl;
      SX_EXIT;
   }
   RelVec miller (hklList(0).toInt (), 
                  hklList(1).toInt (), 
                  hklList(2).toInt ());
   if (miller.absSqr().sum() == 0)  {
      cout << "Invalid hkl = (0,0,0)." << endl;
      SX_EXIT;
   }

   SxSlab slab;

   // --- read input
   SxParser parser;
   SxParser::Table table = parser.read (inFile, "std/structure.std");
   SxAtomicStructure str (&*table);
   str.readElements (&*table);

   if (conventionalCell)  {
      cout << "conventional cell: " 
           << (str.cell ^ str.cell.getConventional () ) << endl;
      SxMatrix3<Double> cTrans = str.cell.getConventional ();
      miller =  (cTrans.inverse ().transpose () * cTrans.determinant ())
             ^ miller;
      SxSlab::shortenVector (&miller);

      cout << "primitive cell: ";
   }

   // str.print (); // TODO!
   cout << "{h, k, l} = " << miller << endl;
   cout << "Surface normal (cartesian coordinates) :";
   cout << (str.cell.inv.transpose () ^ miller) << endl;
   slab.setBulk  (str);
   slab.setPlane (miller);
  
   SxAtomicStructure newSlab;
   newSlab = slab.compute (dir, nLayer, nVacLayer, zOrthogonal, threshold,
                           bondLength, minDist, chemSymbol);

   if (glide.getSize () == 2)  {
      newSlab.cell = SxCell(newSlab.cell.basis(0),
                            newSlab.cell.basis(1),
                            newSlab.cell.basis(2) 
                            + glide(0) * newSlab.cell.basis(0)
                            + glide(1) * newSlab.cell.basis(1));
   }

   FILE *output;
   if (outFile.getSize () == 0)
      output = stdout;
   else
      output = fopen (outFile.ascii(),"w");

   newSlab.fprint (output);
   if (output != stdout) fclose (output);
   
   sxprintf ("Slab successfully generated.\n");
   
   return 0;

}

#endif /* SX_STANDALONE */
