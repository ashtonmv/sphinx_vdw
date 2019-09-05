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

#include <SxUtil.h>
#include <SxParser.h>
#include <SxAtomicStructure.h>
#include <SxPseudoPot.h>
#include <SxGkBasis.h>
#include <SxGBasis.h>
#include <SxRBasis.h>
#include <SxRadBasis.h>
#include <SxAtomicOrbitals.h>
#include <SxQuantumNumbers.h>
#include <SxProjector.h>
#include <SxConstants.h>

//= std/sxatomorb.std
//+ set validation;
//+ include <std/structure.std>; 
//+ include <std/basis.std>; 

//= sxatomorb.sx
//+ format sxatomorb;
//+ include <parameters.sx>;
//+ pseudoPot {
//+    species {
//+       include <species/ga-lda-ham.sx>;
//+    }
//+ }
//+ structure {
//+    cell = 10*[[1,0,0],[0,1,0],[0,0,1]];
//+    species {
//+       atom { coords = [0,0,0]; }
//+    }
//+ }
//+ basis {
//+    eCut = 10;
//+    kPoint { coords = [0,0,0]; weight=1; }
//+ }


/**
  \example atomorb.cpp

  This example demonstrates the Dirac projectors in SFHIngX to visualize
  atomic orbitals.
  The output are SFHIngX Binary Files (*.sxb) which can be vizualized
  with e.g. pxviewer (http://www.phinax.de).

  It computes
  \f[
     \rho_i(R) = |\langle R | \mu_{i_s,i_a,n,l,m}\rangle|^2
               = |\sum_{G} \langle R | G \rangle 
                          \langle G | \mu_{i_s,i_a,n,l,m} \rangle|^2
  \f]
  with \f$ \mu \f$ being the
  atomic orbitals of a given species \f$i_s\f$, atom \f$i_a\f$,
  principle quantum number \em n, angular quantum number \em l, and
  magnetic projection \em m. 
  To visualize the result it is squared in realspace \f$ | R \rangle \f$.

  
  Examples:\par
  - px and pz \image html  px.png 
              \image latex px.png
              \image html  pz.png
              \image latex pz.png
  - dz2:      \image html  dz2.png
              \image latex dz2.png
        

  @author Sixten Boeck
  */
int main ()
{
   // --- we need an atomic structure (incl. pseudopotentials)
   SxParser parser;
   SxParser::Table table = parser.read ("sxatomorb.sx");
   SxAtomicStructure str (&*table);
   SxPseudoPot psPot (&*table);

   // --- basis sets
   SxVector3<Int> mesh (SxGBasis::getMesh(&*table));
   double gCut (SxGBasis::getECut (&*table));

   SxRBasis   R  (mesh, str.cell);                                // |R>
   SxGBasis   G  (mesh, str, gCut);                               // |G>
   SxGkBasis  Gk (G, &*table);                                    // |G+k>
   SxPtr<SxRadBasis> rPtr;
   rPtr = rPtr.create (psPot.rad, psPot.logDr);                   // |r>

   // -- vectors
   SxAtomicOrbitals  orbitals (psPot.pseudoPsi, rPtr);    // <r|>

   int s = SxQuantumNumbers::s;
   int p = SxQuantumNumbers::p;
   int d = SxQuantumNumbers::d;

   // --- write s, p, and d orbitals in R space to disc
   //     use pxviewer for 3d-visualization, e.g. "pxviewer -sxb dx2y2.sxb"
   R.writeMesh3d ("s.sxb",   SUM(G,(R|G)*(G|orbitals(0,0,0,s, 0))).absSqr());
//   R.writeMesh3d ("px.sxb",   SUM(G,(R|G)*(G|orbitals(0,0,0,p,-1))).absSqr());
//   R.writeMesh3d ("py.sxb",   SUM(G,(R|G)*(G|orbitals(0,0,0,p, 0))).absSqr());
//   R.writeMesh3d ("pz.sxb",   SUM(G,(R|G)*(G|orbitals(0,0,0,p, 1))).absSqr());
//   R.writeMesh3d ("dz2.sxb",   SUM(G,(R|G)*(G|orbitals(0,0,0,d, 2))).absSqr());
}
