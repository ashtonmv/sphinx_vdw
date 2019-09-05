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

#ifndef _SX_TYPES_H_
#define _SX_TYPES_H_

#include <SxMatrix3.h>
#include <SxVector.h>
#include <SxVector.h>
#include <SxDirac.h>
#include <SxArray.h>
#include <SxMatrix.h>
#include <SxPrecision.h>
#include <SxBundle3.h>

// --- Fermi distribution
typedef SxBundle3<TPrecEps>                      Eps;
typedef SxBundle3<Int>                           EpsSort;
typedef SxBundle3<TPrecFocc>                     Focc;
// --- wave functions
typedef SxDiracVec<TPrecCoeffG>                  PsiG;
typedef SxDiracVec<TPrecCoeffG>                  PsiGI;
typedef SxDiracVec<TPrecCoeffR>                  PsiR;
// --- potentials, charge density
typedef SxDiracVec<TPrecRhoG>                    SxMeshG;
typedef SxDiracVec<TPrecRhoR>                    SxMeshR;
typedef SxArray<SxMeshR>                         VxcR;    // :iSpin
typedef SxArray<SxMeshG>                         RhoG;    // :iSpin
typedef SxArray<SxMeshR>                         RhoR;    // :iSpin
// --- forces
typedef SxList<SxList<SxVector3<Double> > >      Tau;     // :is,:ia,:{xyz}


#endif /* _SX_TYPES_H_ */
