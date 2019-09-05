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
#ifndef _SX_BLASLIB_H_
#define _SX_BLASLIB_H_

#include <SxComplex.h>
#include <SxUtil.h>
#include <SxConfig.h>
#include <SxMath.h>

enum UPLO   { UpperRight, LowerLeft };
enum EIGCMD { All, ValuesOnly, VectorsOnly, OptSize };

//------------------------------------------------------------------------------
// norm of vectors
//------------------------------------------------------------------------------
SX_EXPORT_MATH float  norm2 (const float *vec, int n);
SX_EXPORT_MATH double norm2 (const double *vec, int n);
SX_EXPORT_MATH float  norm2 (const SxComplex8 *vec, int n);
SX_EXPORT_MATH double norm2 (const SxComplex16 *vec, int n);

//------------------------------------------------------------------------------
// scale vectors
//------------------------------------------------------------------------------
SX_EXPORT_MATH void scale (float *vec,  const float alpha, int n);
SX_EXPORT_MATH void scale (double *vec, const double alpha, int n);
SX_EXPORT_MATH void scale (SxComplex8 *vec,  const SxComplex8 &alpha, int n);
SX_EXPORT_MATH void scale (SxComplex16 *vec, const SxComplex16 &alpha, int n);


//------------------------------------------------------------------------------
// Y += a*X
//------------------------------------------------------------------------------
SX_EXPORT_MATH void axpy (float *yOut, const float &alpha, const float *xIn,
                          int n);
SX_EXPORT_MATH void axpy (double *yOut, const double &alpha, const double *xIn,
                          int n);
SX_EXPORT_MATH void axpy (SxComplex8 *yOut, const SxComplex8 &alpha, 
                          const SxComplex8  *xIn, int n);
SX_EXPORT_MATH void axpy (SxComplex16 *yOut, const SxComplex16 &alpha, 
                          const SxComplex16 *xIn, int n);


//------------------------------------------------------------------------------
// scalar product of vectors
//------------------------------------------------------------------------------
SX_EXPORT_MATH float  scalarProduct (const float *aVec, const float *bVec, 
                                     int n);
SX_EXPORT_MATH double scalarProduct (const double *aVec, const double *bVec, 
                                     int n);
SX_EXPORT_MATH SxComplex8 scalarProduct (const SxComplex8 *aVec, 
                                         const SxComplex8 *bVec, int n);
SX_EXPORT_MATH SxComplex16 scalarProduct (const SxComplex16 *aVec, 
                                          const SxComplex16 *bVec, int n);

//------------------------------------------------------------------------------
// general matrix-matrix multiplication
//------------------------------------------------------------------------------
SX_EXPORT_MATH void matmult (float *resMat, const float *aMat, 
                             const float *bMat, int aMatRows, int aMatCols, 
                             int bMatCols);
SX_EXPORT_MATH void matmult (double *resMat, const double *aMat, 
                             const double *bMat, int aMatRows, int aMatCols, 
                             int bMatCols);
SX_EXPORT_MATH void matmult (SxComplex8 *resMat, const SxComplex8 *aMat,
                             const SxComplex8 *bMat, int aMatRows, 
                             int aMatCols, int bMatCols);
SX_EXPORT_MATH void matmult (SxComplex16 *resMat, const SxComplex16 *aMat, 
                             const SxComplex16 *bMat, int aMatRows, 
                             int aMatCols, int bMatCols);

//------------------------------------------------------------------------------
// overlap matrices
//------------------------------------------------------------------------------
// computes top(A)^T* top(B)
// top(X) restricts matrix to sumSize many rows (from top)
SX_EXPORT_MATH void matovlp (float *resMat, 
                             const float *aMat, const float *bMat, 
                             int aMatRows, int aMatCols, int bMatRows,
                             int bMatCols, int sumSize);
SX_EXPORT_MATH void matovlp (double *resMat, 
                             const double *aMat, const double *bMat, 
                             int aMatRows, int aMatCols, int bMatRows,
                             int bMatCols, int sumSize);
SX_EXPORT_MATH void matovlp (SxComplex8 *resMat, 
                             const SxComplex8 *aMat, const SxComplex8 *bMat, 
                             int aMatRows, int aMatCols, int bMatRows,
                             int bMatCols, int sumSize);
SX_EXPORT_MATH void matovlp (SxComplex16 *resMat, 
                             const SxComplex16 *aMat, const SxComplex16 *bMat, 
                             int aMatRows, int aMatCols, int bMatRows,
                             int bMatCols, int sumSize);
//------------------------------------------------------------------------------
// Matrix decompositions
//------------------------------------------------------------------------------
SX_EXPORT_MATH void cholesky (float *resMat, enum UPLO,   // modifies inMat!!!
                              float *inMat, int n);
SX_EXPORT_MATH void cholesky (double *resMat, enum UPLO, 
                              double *inMat, int n);
SX_EXPORT_MATH void cholesky (SxComplex8 *resMat, enum UPLO, 
                              SxComplex8 *inMat, int n);
SX_EXPORT_MATH void cholesky (SxComplex16 *resMat, enum UPLO, 
                              SxComplex16 *inMat, int n);

SX_EXPORT_MATH void singularValueDecomp (float *mat, int nRows, int nCols,
                                         float *vals,
                                         float *left,
                                         float *right, // V^H
                                         bool zeroSpace);
SX_EXPORT_MATH void singularValueDecomp (double *mat, int nRows, int nCols,
                                         double *vals,
                                         double *left,
                                         double *right, // V^H
                                         bool zeroSpace);
SX_EXPORT_MATH void singularValueDecomp (SxComplex8 *mat, int nRows, int nCols,
                                         float *vals,
                                         SxComplex8 *left,
                                         SxComplex8 *right, // V^H
                                         bool zeroSpace);
SX_EXPORT_MATH void singularValueDecomp (SxComplex16 *mat, int nRows, int nCols,
                                         double *vals,
                                         SxComplex16 *left,
                                         SxComplex16 *right, // V^H
                                         bool zeroSpace);

//------------------------------------------------------------------------------
// matrix inversion
//------------------------------------------------------------------------------
SX_EXPORT_MATH void matInverse (float       *mat, int nRows, int nCols);
SX_EXPORT_MATH void matInverse (double      *mat, int nRows, int nCols);
SX_EXPORT_MATH void matInverse (SxComplex8  *mat, int nRows, int nCols);
SX_EXPORT_MATH void matInverse (SxComplex16 *mat, int nRows, int nCols);
SX_EXPORT_MATH void matInverseTri (float       *mat, int nRows, enum UPLO);
SX_EXPORT_MATH void matInverseTri (double      *mat, int nRows, enum UPLO);
SX_EXPORT_MATH void matInverseTri (SxComplex8  *mat, int nRows, enum UPLO);
SX_EXPORT_MATH void matInverseTri (SxComplex16 *mat, int nRows, enum UPLO);

//------------------------------------------------------------------------------
// Linear equation solver (least sqare based)
//------------------------------------------------------------------------------
SX_EXPORT_MATH void solveLinEq (float *mat, int nRows, int nCols, 
                                float *b,   int bCols);
SX_EXPORT_MATH void solveLinEq (double *mat, int nRows, int nCols, 
                                double *b,   int bCols);
SX_EXPORT_MATH void solveLinEq (SxComplex8 *mat, int nRows, int nCols, 
                                SxComplex8 *b,   int bCols);
SX_EXPORT_MATH void solveLinEq (SxComplex16 *mat, int nRows, int nCols, 
                                SxComplex16 *b,   int bCols);

//------------------------------------------------------------------------------
// eigensolver
//------------------------------------------------------------------------------
// modifies inMat!!!
SX_EXPORT_MATH int matEigensolver (SxComplex8  *eigVals, float *eigVecs, 
                                   float *inMat, int n, EIGCMD cmd=All, 
                                   int size=0);
SX_EXPORT_MATH int matEigensolver (SxComplex16 *eigVals, double *eigVecs, 
                                   double *inMat, int n, EIGCMD cmd=All, 
                                   int size=0);
SX_EXPORT_MATH int matEigensolver (SxComplex8  *eigVals, SxComplex8  *eigVecs,
                                   SxComplex8  *inMat, int n, EIGCMD cmd=All,
                                   int size=0);
SX_EXPORT_MATH int matEigensolver (SxComplex16 *eigVals, SxComplex16 *eigVecs,
                                   SxComplex16 *inMat, int n, EIGCMD cmd=All,
                                   int size=0);
// modified inMat!!!
SX_EXPORT_MATH void matEigensolverTri (float *eigVals, float *eigVecs,
                                       float *inMat,
                                       int n, enum UPLO, EIGCMD cmd=All);
SX_EXPORT_MATH void matEigensolverTri (double *eigVals, double *eigVecs,
                                       double *inMat,
                                       int n, enum UPLO, EIGCMD cmd=All);
SX_EXPORT_MATH void matEigensolverTri (float *eigVals, SxComplex8 *eigVecs,
                                       SxComplex8 *inMat,
                                       int n, enum UPLO, EIGCMD cmd=All);
SX_EXPORT_MATH void matEigensolverTri (double *eigVals, SxComplex16 *eigVecs,
                                       SxComplex16 *inMat,
                                       int n, enum UPLO, EIGCMD cmd=All);


//------------------------------------------------------------------------------
// Some missing pototypes
//------------------------------------------------------------------------------
extern "C" {
//   zgetri (int, void *, int, int *, void *, int, int &);
}


#endif /* _SX_BLASLIB_H_ */
