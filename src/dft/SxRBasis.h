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

#ifndef _SX_R_BASIS_H_
#define _SX_R_BASIS_H_

#include <SxString.h>
#include <SxBasis.h>
#include <SxCell.h>
#include <SxSymGroup.h>
#include <SxFFT3d.h>
#include <SxBinIO.h>
#include <SxTypes.h>
#include <SxDFT.h>


/** 
  \author Sixten Boeck
  */
class SxGBasis;

/** \brief Realspace basis

  \b SxRBasis = S/PHI/nX Realspace Basis

  \ingroup group_dft
  \author  Sixten Boeck
  */
class SX_EXPORT_DFT SxRBasis : public SxBasis
{
   public:

      typedef TRBasisType               TBasisType;
      typedef SxDiracVec<TPrecCoeffR>   TPsi;
      typedef SxDiracVec<TPrecRhoR>     TRho;

      /// Return basis type
      virtual SxString getType () const { return "|R>"; }

      SxRBasis ();
      SxRBasis (int nx, int ny, int nz, const SxCell &);
      SxRBasis (const SxVector3<Int> &, const SxCell &);

      /// Initialize basis
      void set (const SxVector3<Int> &, const SxCell &);
      /// Initialize basis
      void set (int nx, int ny, int nz, const SxCell &);

      virtual ~SxRBasis ();

      template<class T>
      SxDiracVec<T> symmetrize (const SxDiracVec<T> &) const;
      void writeMesh3d (const SxString &file, const TRho &mesh3d) const;
      RhoR readMesh3d  (const SxString &file) const;

      int getMeshSize () const { return fft3d.meshSize; }
      SxVector3<Int> getMesh () const { return fft3d.mesh; }

      virtual void print () const;
      void setNComp ( int comp);
      virtual int getNComp () const;

      /** \brief \f$ \Delta \Omega\f$ used for integration

          Sometimes realspace entities (like the charge density SxRho) have
          to be integrated. Since the realspace is represented on a regular
          grid (FFT grid) such integrations can numerically easily evaluated.
          If X(R) is an entity living on the realspace grid its integral 
          is
          \f[
             \int X(R) \Omega \Rightarrow \sum_i X(R_i) \Delta \Omega
          \f] and
          \f[ 
             \Delta \Omega = \frac{\Omega}{n_\mathrm{FFT}}
          \f]
          The number of realspace (FFT) mesh points are \f$ n_\mathrm{FFT}.\f$
          
       */
      Real8             dOmega;

      /** \brief number of components */
      int               nComp;

      /// Return number of mesh points
      virtual ssize_t getNElements () const {
         return fft3d.meshSize * nComp;
      }
//   protected:

      /// \brief The FFT object used to transform objects between |G> and |R> 
      mutable SxFFT3d   fft3d;

      /// The unit cell
      SxCell            cell;

   protected:
      /// The dual G-basis
      mutable const SxGBasis* gBasisPtr;

   public:

      /** Register the dual G basis */
      inline void registerGBasis (const SxGBasis &gBasis) const;

      /** Get the dual G basis */
      const SxGBasis& getGBasis () const
      {
         SX_CHECK(gBasisPtr);
         return *gBasisPtr;
      }

      /** Deregister a basis */
      virtual void deregister (const SxBasis *basis) const;

      REGISTER_PROJECTOR (SxRBasis, SxRBasis, identity);
      REGISTER_PROJECTOR (SxRBasis, SxGBasis, toGSpace);

      virtual Real8 tr (const SxDiracVec<Double> &) const;
      virtual Real8 tr (const SxDiracVec<Complex16> &) const;

      SxDiracVec<TBasisType>  identity (const SxRBasis *,
                                        const SxDiracVec<TBasisType> &) const;

      SxDiracVec<TPrecCoeffG> toGSpace (const SxGBasis *, 
                                        const SxDiracVec<TBasisType> &) const;
};

namespace Timer {
   enum RGBasisTimer {
      GBasis,
      Sym,
      GtoR, GtoR_FFT,
      RtoG, RtoG_FFT
   };
}

SX_REGISTER_TIMERS (Timer::RGBasisTimer)
{
   using namespace Timer;
   regTimer (GBasis,     "|G> Setup");
   regTimer (Sym,        "Symmetrization");
   regTimer (GtoR,       "G->R routine");
   regTimer (GtoR_FFT,   "G->R FFT call");
   regTimer (RtoG,       "R->G routine");
   regTimer (RtoG_FFT,   "R->G FFT call");
}

//------------------------------------------------------------------------------
// Symetrization of a 3d realspace mesh
//------------------------------------------------------------------------------
template<class T>
SxDiracVec<T> SxRBasis::symmetrize (const SxDiracVec<T> &meshIn) const
{
   SX_CHECK (meshIn.getBasisPtr() == this);
   SX_CLOCK (Timer::Sym);

   // --- do we have to symmetrize mesh???
   if (!cell.symGroupPtr || cell.symGroupPtr->getNSymmorphic () <= 1)  {
      return meshIn;
   }

   const SxSymGroup &S = *cell.symGroupPtr;
   int nOp = S.getNSymmorphic ();

   // --- Note, that SxCell is constructed column-wisely whereas
   //     SxMatrix3<T> is constructed row-wisely currently.
   SxCell meshCell (cell(0) / fft3d.mesh(0),
                    cell(1) / fft3d.mesh(1),
                    cell(2) / fft3d.mesh(2));//.transpose();

   SxArray<SxMatrix3<Int> >  symOpRel (nOp);
   for (int iOp=0; iOp < nOp; iOp++)
      symOpRel(iOp) = meshCell.carToRel (S.getSymmorphic(iOp));

   SxVector3<Int> mesh = fft3d.mesh;
   SxDiracVec<T> meshOut (meshIn.getSize());
   meshOut.setBasis (meshIn.getBasisPtr());


   double nOpInv = 1. / double(nOp);
#ifdef USE_OPENMP
#pragma omp parallel
#endif
   {
      SxArray<ssize_t> idxRot(nOp);
      const int nj = mesh(0), nk = mesh(1), nl = mesh(2);
#ifdef USE_OPENMP
#pragma omp for collapse(2)
#endif
      for (int j=0; j < nj; j++)  {
         for (int k=0; k < nk; k++)  {
            for (int l=0; l < nl; l++)  {
               SxVector3<Int> i(j,k,l);
               ssize_t iOut = fft3d.mesh.getMeshIdx(i, SxMesh3D::Positive);

               // --- collect all indices
               bool done = false;
               for (int iOp = 0; iOp < nOp; iOp++)  {
                  SxVector3<Int> rot = symOpRel(iOp) ^ i;
                  ssize_t idx = fft3d.mesh.getMeshIdx(rot, SxMesh3D::Unknown);
                  // only lowest index does the computation
                  if (idx < iOut) { done = true; break; }
                  idxRot(iOp) = idx;
               }
               // only lowest index does the computation
               if (!done)  {
                  // --- compute sum of all equivalent elements ...
                  typename T::Type res = 0.;
                  for (int iOp = 0; iOp < nOp; iOp++)
                     res += meshIn(idxRot(iOp));
                  res *= nOpInv;
                  // --- distribute result
                  for (int iOp = 0; iOp < nOp; iOp++)
                     meshOut(idxRot(iOp)) = res;
               }
            }
         }
      }
   }

   return meshOut;
}

#include <SxGBasis.h>

void SxRBasis::registerGBasis (const SxGBasis &gBasis) const
{
   registerBasis (gBasis);
   gBasisPtr = &gBasis;
}
#endif /* _SX_R_BASIS_H_ */
