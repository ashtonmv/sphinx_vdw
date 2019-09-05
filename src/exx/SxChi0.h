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

#ifndef _SX_CHI_0_H_
#define _SX_CHI_0_H_

#include <SxExx.h>
#include <SxDirac.h>
#include <SxPW.h>
#include <SxFermi.h>
#include <SxRBasis.h>
#include <SxGBasis.h>
#include <SxGkBasis.h>

/** \brief Linear response \f$\chi_{0}\f$ in G space used in classical EXX
           calculations

    \b SxChi0 = S/PHI/nX linear response function

    ....

    \author Matthias Wahn, wahn@fhi-berin.mpg.de */
class SX_EXPORT_EXX SxChi0
{
   public:

      ///@name constructors and destructor
      //@{
      /** contructor */
      SxChi0 ();

      /** constructor */
      SxChi0 (const SxGBasis &G, const SxGkBasis &gk, int nG_);

      /** constructor */
      SxChi0 (SxGBasis &G, SxGkBasis &gk, int nG_, const SxCell &cell);

      /** constructor */
      SxChi0 (SxGBasis &G, SxGkBasis &gk, const SxCell &cell);

      /** destructor */
      ~SxChi0 ()  { /* empty */ }
      //@}

      /** computes the linear response function \f$\chi_{0}\f$ in
          reciprocal space

          \note Since the \f$(\mathbf{G}=\mathbf{0})\f$-components are
                cut away in \f$\chi_{0}\f$, the result is an (nG-1) x (nG-1)
                matrix.
       */
      void compute (const SxPW &waves, const SxFermi &fermi);
      void computeOld (const SxPW &waves, const SxFermi &fermi);

      /** computes the eigenvalues of \f$\chi_{0}\f$ */
      void computeEigVals ();

      /** computes eigenvalues and eigenvectors of \f$\chi_{0}\f$ */
      void computeEigSys ();

      /** \return the eigenvalues of \f$\chi_{0}\f$ */
      PsiG getEigVals ();

      /** \return the eigenvectors of \f$\chi_{0}\f$ */
      SxDiracMat<TPrecCoeffG> getEigVecs ();

      /** \return \f$ n_{\mathbf G}\f$ is number of G vector used for
                  \f$\chi_{0}\f$. Since the (G=0)-component is cut,
                  \f$\chi_{0}\f$ is stored as a \f$ (n_{\mathbf G}-1)\f$
                  square matrix. */
      int getNG () const;

      /** \return the matrix \f$\chi_{0}\f$ in reciprocal space */
      SxDiracMat<TPrecCoeffG> getChi0 ();

      /** \return the inverse of the linear response function in the
                  \f$(G\not=0)\f$ subspace of reciprocal space */
      SxDiracMat<TPrecCoeffG> getInverse ();

      /** \return the diagonal elements
                 \f$\chi_{0}(\mathbf{G},\mathbf{G})\f$ */
      SxDiracVec<TPrecCoeffG> getDiag ();

      /** \return the /f$(\mathbf{G}_{1},\mathbf{G}_{2})/f$-th element of
                  \f$\chi_{0}\f$ */
      PrecCoeffG operator() (const int iG, const int jG);

      /** \return the i-th eigenvalue */
      PrecCoeffG getEigVal (const int iG);

      /** \return the i-th eigenvector in G space */
      PsiG getEigVecG (const int iG);

      /** \return the i-th eigenvector in r space */
      PsiG getEigVecR (const int iG);

      /** multiplies \f$\chi_{0}\f$ as matrix with a vector */
      PsiG operator* (const PsiG &v);

      ///@name I/O
      //@{
      /** writes \f$\chi_{0}\f$ and its eigensystem to output file */
      void write (const SxString &filename) const;

      /** reads in \f$\chi_{0}\f$ and its eigensystem */
      void read (const SxString &filename);

      /** reads in the eigensystem only */
      void readEigSys (const SxString &filename) const;
      //@}

      ///@name Debugging
      //@{
      /** \brief test FFT prefactor

        When porting the code to other systems, there might occur
        some inconsistency with the FFT prefactors. This function is
        meant to have an easy check on this.

        The computation of \f$ \chi_{0} \f$ requires expressions of the
        form
        \f[ \langle v\mathbf{k} | \mathrm{e}^{-\mathrm{i}\mathbf{G}\cdot
                    \mathbf{r}} | c\mathbf{k} \rangle\; ,
        \f]
        where \f$ |v\mathbf{k}\rangle \f$ and \f$ |v\mathbf{k}\rangle \f$
        shall denote valence and conduction band states, respectively,
        and the operator is given in its real space form.

        This is conveniently performed using the Fourier trafo:
        \f[ \langle v\mathbf{k} | \mathrm{e}^{-\mathrm{i}\mathbf{G}\cdot
                    \mathbf{r}} | c\mathbf{k} \rangle
           = \int\limits_{\Omega} \mathrm{e}^{-\mathrm{i}\mathbf{G}\cdot
                   \mathbf{r}} \, \psi_{v\mathbf{k}}^{\ast}(\mathbf{r})\,
                   \psi_{c\mathbf{k}}(\mathbf{r})\, \mathrm{d}\mathbf{r}\; ,
        \f]
        with \f$\Omega\f$ being the volume of the unit cell.

        This test now just performs
        \f[ \langle v\mathbf{k} | \mathrm{e}^{-\mathrm{i}\mathbf{G}\cdot
                    \mathbf{r}} | v\mathbf{k} \rangle
        \f]
        for the lowest valence band at one \f$\mathbf{k}\f$-point. Note,
        that also the ket is a valence band here. If all prefactors are
        correct, the \f$(\mathbf{G}=\mathbf{0})\f$-component must be one,
        due to normalization.
       */
      void checkFFTfactor (const SxPW &waves);
      //@}

   protected:

      /** the matrix \f$\chi_{0}\f$ in reciprocal space */
      SxDiracMat<TPrecCoeffG>    Chi0;

      /** number of G components used in \f$\chi_{0}\f$ */
      int                        nG;

      /** pointer to basis in reciprocal space */
      const SxGBasis            *GPtr;

      /** pointer to the |G+k> basises */
      const SxGkBasis           *gkPtr;

      /** real space basis */
      SxRBasis                   R;

      /** number of k-points */
      const int                  nk;

      /** number of valence bands */
      int                        nv;

      /** number of conduction bands */
      int                        nc;

      /** eigenvalues of \f$\chi_{0}\f$ */
      SxDiracVec<TPrecCoeffG>    eigVals;

      /** eigenvectors of \f$\chi_{0}\f$ */
      SxDiracMat<TPrecCoeffG>    eigVecs;

      /** FFT prefactor used for the transformation of a function with the
          dimension of a charge density

          \note Check this factor using SxChi0::checkFFTfactor.
       */
      SX_FFT_REAL                fftFactorRtoG;

      /** size of G mesh */
      int                        sizeG;
};

#endif /* _SX_CHI_0_H_ */
