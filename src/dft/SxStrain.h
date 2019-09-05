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
#ifndef _SX_STRAIN_H_
#define _SX_STRAIN_H_

#include <SxDFT.h>
#include <SxRBasis.h>
#include <SxBinIO.h>
#include <SxVector3.h>
#include <SxRho.h>
#include <SxHamiltonian.h>
#include <SxFermi.h>
#include <SxTypes.h>

/** \brief SxStrain

  \b SxStrain = Continuum theory strain field calculation
  from \ref Pryor97

  \author Oliver Marquardt, marquardt@mpie.de
  */

class SX_EXPORT_DFT SxStrain : public SxHamiltonian
{
   public:

      SxStrain (); // { }
      /// Constructor
      SxStrain (const SxPWSet &,
                const SxRBasis* rPtr);

      /// Destructor
      virtual ~SxStrain ();

      static const SxSymbolTable * getHamiltonianGroup (const SxSymbolTable *);
      PsiG operator* (const PsiG &) const;
      /// \brief Apply strain "Hamiltonian"
      // hook in for SxHamiltonian
      virtual PsiG apply (const PsiG &psi) const
      {
         return SxStrain::operator* (psi);
      }

   protected:
      const SxRBasis *rBasisPtr;
      SxRBasis R1;
      SxGBasis G1;
      SxArray<double> cP;
      SxVector3<Double> wgt;
      SxString inFile;

   public:
      virtual PrecEnergy getETrial();
      virtual void printEnergies () const;

      virtual void compute (const SxFermi &, bool, bool) {}
      virtual void writePotentials () const;
      /** \brief Write a R-mesh as ASCII file "<name>.dat"
          @param name name of the mesh
          @param data content of the mesh
        */
      void writeMeshASCII (const SxString &name, const PsiR &data) const;
      /** \brief Write a R-mesh as sxb file "<name>.sxb"
          @param name name of the mesh
          @param data content of the mesh
        */
      void writeMeshSxb (const SxString &name, const SxMeshR &data) const;

      virtual void set (const SxPsiSet &, const SxFermi &) {}
      virtual void read (const SxSymbolTable *); // --- read input parameter
      virtual PrecEnergy getEnergy (const SxPsiSet &,
                                    const SxFermi &);
      virtual double getEnergy () const { return 0.; }
   protected:
      bool first, moreIO, containsVacuum, readExternalStrains;
      SxArray<PsiR> firstDerRComp (const PsiR &) const;
      SxArray<PsiR> gIntegrate (const PsiG &);
      SxArray<PsiG> gVecs, gG2;
      PsiG scaledG2;
      /** \brief Initialize the gVecs (after wgt is known)

           @param gBasis3comp 3-component G basis for u
        */
      void initGvec (const SxGBasis &gBasis3comp);
      double energy, dV, fMax, fMin;
      double interpolate(double, double, double, double);
      SxMesh3D mesh3d;
      SxArray<PsiR> P, P2, P2a;
      /** \brief Gradient is calculated as follows:

        in Zincblend structure:
        \f[\frac{\delta F}{\delta u_i} = C_{11} \frac{\partial^2 u_i}{\partial r_i^2} +
        C_{12} \left(
          \frac{\partial^2 u_{i+1}}{\partial r_{i+1}\partial r_i}
        + \frac{\partial^2 u_{i+2}}{\partial r_{i+2}\partial r_i}\right)\f]

        \f[ + \frac{1}{2}C_{44}\left(\frac{\partial^2 u_x}{r_{i+1}^2} +
        2 \left[\frac{\partial^2 u_{i+1}}{\partial r_{i+1}\partial r_i} +
        \frac{\partial^2 u_{i+2}}{\partial r_{i+2}\partial r_i}\right] +
        \frac{\partial^2 u_{i+1}}{\partial r_{i+1}\partial r_i}\right)\f]

        \f[ -\frac{\partial \alpha e_0}{\partial r_i}\f]
        where the \f$C_{ij}\f$'s are the elastic constants and \f$u_i\f$ is the
        displacement vector. The strain field components
        \f$\sigma_{ij}\f$ are calculated as follows:
        \f[\sigma_{ij} = \frac{1}{2}\left(\frac{\partial u_i}{\partial r_j}
        + \frac{\partial u_j}{\partial u_i}\right)\f]
        in Wurtzite structure:
        \f[\frac{\delta F}{\delta u_x} = C_{11}\frac{\partial ^2 u_x}{\partial x^2}
        + C_{12}\frac{\partial^2 u_y}{\partial x\partial y}
        + C_{13}\frac{\partial^2 u_z}{\partial x\partial z}\f]

        \f[ + \frac{1}{2}C_{44}\left(\frac{\partial^2 u_x}{\partial z^2}
        + \frac{\partial u_z}{\partial z\partial x}\right)
        + \left(\frac{1}{2}C_{11} - C_{12}\right)
        \frac{\sigma_{xx}^{(0)}}{\partial x}
        -\frac{1}{2}C_{13}\frac{\partial \sigma_{zz}^{(0)}}{\partial x}\f]

        similar for \f$\frac{\delta F}{\delta u_y}\f$ and
        \f[\frac{\delta F}{\delta u_z} = C_{33}\frac{\partial ^2 u_x}{\partial x^2}
        + C_{13}\left(\frac{\partial^2 u_y}{\partial x\partial y}
        + \frac{\partial^2 u_z}{\partial x\partial z}\right)\f]

        \f[ + \frac{1}{2}C_{44}\left(\frac{\partial^2 u_x}{\partial z\partial x}
        + \frac{\partial^2 u_y}{\partial z\partial y}
        + \frac{\partial^2 u_z}{\partial x^2}
        + \frac{\partial^2 u_z}{\partial y^2}\right)
        -\frac{1}{2}C_{33}\frac{\sigma_{zz}^{(0)}}{\partial z}
        -C_{13}\frac{\sigma_{xx}^{(0)}}{\partial x} \f]

        with
        \f[\sigma_{xx}^{(0)} = \frac{a_{dot} - a_{matrix}}{a_{matrix}}\f]
        \f[\sigma_{zz}^{(0)} = \frac{c_{dot} - c_{matrix}}{c_{matrix}}\f]
        The piezoelectric potential is calculated from the polarization:
        \f[\varepsilon_0\varepsilon_r\Delta V_p(\mathbf{r})
        = - \vec{\nabla}\mathbf{P}(\mathbf{r})\f]
        \f[\varepsilon_0\varepsilon_r \mathbf{G}^2V_p(\mathbf{r})
        = i\mathbf{G}\mathbf{P}(\mathbf{G})\f]
        \f[V_p(\mathbf{G})=
        \frac{i\mathbf{G}}{\varepsilon_0\varepsilon_r\mathbf{G}^2}\mathbf{P}(\mathbf{G})\f]
        The polarization vector \f$\mathbf{P}\f$ is:
        \f[\mathbf{P} = \left(\begin{array}{c}
        e_{14}\sigma_{yz}\\
        e_{14}\sigma_{xz}\\
        e_{14}\sigma_{xy}
        \end{array}\right)\f]
        in Zincblend systems and
         \f[\mathbf{P} = \left(\begin{array}{c}
        e_{15}\sigma_{xz}\\
        e_{15}\sigma_{yz}\\
        e_{31}(\sigma_{xx} + \sigma_{yy}) + e33\sigma_{zz} + P_{spont}
        \end{array}\right)\f]
        in Wurtzite systems.
        */
      PsiR gradient (const PsiG &) const;
      /** \brief Calculate u in real space and decompose into components

           \note u in real space is rescaled!
        */
      SxArray<PsiR> uInR (const PsiG u) const;
      /// Return elastic constants (ijkl format)
      SxMeshR elConst (int, int, int, int) const;
      /// Check if Cijkl is non-zero
      bool nonZeroElement (int, int, int, int) const;

      enum Preconditioner { Payne, Arias };
   public:
      virtual SxDiracVec<TPrecCoeffG::TReal> preconditioner (const PsiG &psi) const {
         return preconditioner (psi, Payne);
      }
   protected:
      /// diagonal of gradient for preconditioner
      PsiG Diag;

      SxDiracVec<TPrecCoeffG::TReal>
      preconditioner (const PsiG &, Preconditioner type) const;

      mutable int step;
      int nMat, nParam, rSize, vacMat;
	   SxMatrix<Double> matParam; // all parameters of all materials
      SxArray<SxMeshR> materials, eIJ;
      SxMeshR parameters (int iParam);
      SxList<SxString> paramNames; // names of involved parameters
      PsiR one;
      mutable PsiG uStore;

      int firstOutsideBracket(SxString, char);
      SxArray<PsiR> parsePolarisation ();
      PsiR parseElement (SxString);
      double parseFunction (SxString);
      SxString replaceUnknown(SxString);

      void calculateSigma (const SxArray<PsiR> &uR);

      void calculateSigma (const PsiG &, const SxArray<PsiR> &);
      double getFreeEnergy();
      double fTot;
      SxMeshR strain(int i, int j);
      // --- global parameters
      SxArray<SxMeshR> oldC; // OLD CODE, elastic 'tensor'=vector
      SxArray<SxArray<SxMeshR> > C; // NEW CODE, elastic tensor
      SxArray<SxMeshR> latConst; // NEW CODE, lattice constants
      SxArray<SxArray<bool> > nonzeroC; // NEW CODE, =false if elastic constant cIJ = 0
      bool generalized, generalizedLat;  // To be removed in following versions.
      SxMeshR c11S, c33S, c44S, c66S, c12S, c13S, c15S;
      int crystalStruct;
      bool secondOrderPiezo, AParams, rotated, c15given;
      //double mixing(double, double, double, double, double);
      SxList <SxString> extChargeFiles;
      SxMeshR sigma0X, sigma0Y, sigma0Z, eXX, eYY, eZZ, eXY, eXZ, eYZ,
      pSpont, e14, e15, e31, e33, vP, epsilon, extChg, extPot,
      B114, B124, B156, A1, A2, Tr;
      SxArray<PsiR> sigma03;
      PsiR sigmaXX, sigmaXY, sigmaXZ, sigmaYY, sigmaYZ, sigmaZZ, zero,
           sigmaXXRot, sigmaYYRot, sigmaZZRot, sigmaXYRot, sigmaXZRot,
           sigmaYZRot;
};
#endif /* _SX_STRAIN_H_ */
