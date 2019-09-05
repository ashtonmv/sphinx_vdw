// ---------------------------------------------------------------------------
//
//           The general purpose cross platform C/C++ framework
//
//                       S x A c c e l e r a t e
//
//           Home:       https://www.sxlib.de
//           License:    Apache 2
//           Authors:    see src/AUTHORS
//
// ---------------------------------------------------------------------------

#ifndef _SX_PRECISION_H_
#define _SX_PRECISION_H_

#include <SxComplex.h>
#include <SxTypeMapper.h>

// --- precision of low level types
typedef int                Int4;
typedef long int           Int8;
typedef float              Real4;
typedef double             Real8;
typedef SxComplex8         Cmplx8;
typedef SxComplex16        Cmplx16; 

// --- precision of type mappers
typedef SxTypeMapper<Real4, Real4>             TReal4;
typedef SxTypeMapper<Real8, Real8>             TReal8;
typedef SxTypeMapper<Cmplx8, Real4,  Cmplx8>   TCmplx8;
typedef SxTypeMapper<Cmplx16, Real8, Cmplx16>  TCmplx16;

// --- precision of several physical variables
//     beware of mixing different precisions. mixing may cause
//     slow down the code due to many type cast commands!
typedef TReal8             TPrecG;
typedef Real8              PrecG;
typedef Int                TPrecFFTIdx;
typedef Int4               PrecFFTIdx;
                        
//#define WAVES_SINGLE_PREC  1
#ifdef  WAVES_SINGLE_PREC
   typedef TCmplx8         TPrecCoeffG;
   typedef Cmplx8          PrecCoeffG;
   typedef TCmplx8         TPrecPotNl;
   typedef Cmplx8          PrecPotNl;
   typedef TCmplx16        TPrecCoeffR;
   typedef Cmplx16         PrecCoeffR;
#else  /* WAVES_SINGLE_PREC */
   typedef TCmplx16        TPrecCoeffG;
   typedef Cmplx16         PrecCoeffG;
   typedef TCmplx16        TPrecPotNl;
   typedef Cmplx16         PrecPotNl;
   typedef TCmplx16        TPrecCoeffR;
   typedef Cmplx16         PrecCoeffR;
#endif /* WAVES_SINGLE_PREC */
   
typedef TCmplx16        TPrecPhase;
typedef Cmplx16         PrecPhase;
typedef TReal8          TPrecPhi;
typedef Real8           PrecPhi;
 
typedef TCmplx16           TPrecPhaseG;
typedef Cmplx16            PrecPhaseG;
typedef TReal8             TPrecPhiG;
typedef Real8              PrecPhiG;
typedef TCmplx16           TPrecRhoG;
typedef Cmplx16            PrecRhoG;
typedef TReal8             TPrecRhoR;
typedef Real8              PrecRhoR;
typedef TCmplx16           TPrecEffPotG;
typedef Cmplx16            PrecEffPotG;   
typedef TReal8             TPrecEffPotR;
typedef Real8              PrecEffPotR;
                        
typedef TReal8             TPrecEps;
typedef Real8              PrecEps;
typedef TReal8             TPrecWeights;
typedef Real8              PrecWeights;
typedef TReal8             TPrecFocc;
typedef Real8              PrecFocc;
typedef TReal8             TPrecEnergy;
typedef Real8              PrecEnergy;    

typedef TReal8             TPrecTauR;
typedef Real8              PrecTauR;
                        
typedef TReal8             TPrecForcesR;
typedef Real8              PrecForcesR;   
typedef TCmplx16           TPrecForcesG;
typedef Cmplx16            PrecForcesG;   

// --- types of Basis elements, see SxBasis
typedef TPrecCoeffR  TRBasisType;
typedef TPrecCoeffG  TGBasisType;
typedef TPrecCoeffG  TGkBasisType;
typedef Double       TRadBasisType;
typedef Double       TRadRBasisType;
typedef Double       TRadGBasisType;
typedef TPrecCoeffR  TRRBasisType;
typedef Complex16    TAOBasisType;
typedef Complex16    TWannierBasisType;



#ifdef FFT_SINGLE_PREC
   typedef float         SX_FFT_REAL;
   typedef SxComplex8    SX_FFT_COMPLEX;
   typedef Complex8      SX_T_FFT_COMPLEX;
#else
   typedef double        SX_FFT_REAL;
   typedef SxComplex16   SX_FFT_COMPLEX;
   typedef Complex16     SX_T_FFT_COMPLEX;
#endif /* FFT_SINGLE_PREC */


#endif  // _SX_PRECISION_H_

