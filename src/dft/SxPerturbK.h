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

#ifndef _SX_PERTURB_K_H_
#define _SX_PERTURB_K_H_

#include <SxGkBasis.h>
#include <SxTimer.h>
#include <SxPseudoPot.h>
#include <SxQuantumNumbers.h>
#include <SxDFT.h>

/** \brief \f$k\cdot p\f$ perturbation theory

    \b SxClass = S/PHI/nX \f$k\cdot p\f$ perturbation theory class

    This class provides the perturbation operator

    \f[ \frac{\partial}{\partial\mathbf k} \hat H_{\mathbf k} \f]

    for the plane-wave band structure Hamiltonian \f$ \hat H_{\mathbf k}\f$.

    Note that it is a vector-like operator.

    The band-structure Hamiltonian works on lattice-periodic functions
    \f$u_{n\mathbf k}\f$ and corresponds to the Hamiltonian for the 
    corresponding Bloch wave 
    \f[
    \psi_{n\mathbf k} = u_{n\mathbf k} \cdot e^{i\mathbf k \cdot \mathbf r}
    \f]

    It reads
    \f[
    \hat H_{\mathbf k} = -\frac 12(\nabla + \mathbf k)^2 + \hat v_{eff}
    + \sum_\phi |G><G+k|\phi>\frac{1}{E_{KB}}<\phi|G+k><G|
    \f]

    The last term is the non-local pseudopotential in the Kleinman-Bylander 
    form, with 
    \f[
    \phi(R,\Omega) = v_l(R) \cdot \mu_l(R) \cdot Y_{lm} (\Omega)
    \f]

    The derivative with respect to k has no contributions from the effective
    potential term, because it doesn't depend on k. The kinetic part yields
    \f[
    \frac{\partial}{\partial\mathbf k} \hat T_{\mathbf k}
    = - \mathbf k \cdot (\nabla + \mathbf k) = \mathbf k \cdot \hat{\mathbf p}
    \f]
    This expression has given the name to the whole story.

    The non-local contribution is very often neglected, but there is
    no good reason to do so generally. It reads

    \f[
    \frac{\partial}{\partial\mathbf k} \hat V^{nl}_{\mathbf k}
    = \sum_\phi |G><G+k|-i\mathbf r \phi>\frac{1}{E_{KB}}<\phi|G+k><G|
              \qquad + \qquad
                |G><G+k|\phi>\frac{1}{E_{KB}}<-i\mathbf r\phi|G+k><G|
    \f]

    This expression can be evaluated by noting that 
    \f$ i \mathbf r \f$ is the product of a real spherical harmonic
    \f$X_{lm}\f$
    with l=1 and a trivial radial function 
    \f$ f(R) = \sqrt\frac{4\pi}{3}\cdot R\f$. Thus
    \f[
    |ir\phi> = \sum_{m=-1}^1 R \cdot X_{1m}(\Omega)  
                             \cdot v_l(R) \mu_l (R) \cdot X_{lm'}(\Omega)
             = \sum_{m=-1}^1 R \cdot v_l(R) \mu_l (R)  
               \cdot X_{1m}(\Omega)  \cdot X_{lm'}(\Omega)
    = \phi_{\mathbf r}
    \f]

    By recoupling the spherical harmonics with the help of Clebsch-Gordan
    coefficients, i.e.
    \f[
    X_{lm} X_{l'm'} = \sum_{LM} <l l' m m' | L M> X_{LM}
    \f]
    the product of spherical harmonics can be resolved, so that the 
    projection onto the plane-waves is straight-forward.
    Thus, a closed expression is obtained, very similar to the original 
    Kleinman-Bylander form:
    \f[
    \frac{\partial}{\partial\mathbf k} \hat V^{nl}_{\mathbf k}
    = - \sum_\phi |G><G+k|\phi_{\mathbf r}>\frac{1}{E_{KB}}<\phi|G+k><G|
                \qquad + \qquad 
                |G><G+k|\phi>\frac{1}{E_{KB}}<\phi_{\mathbf r}|G+k><G|
    \f]

    Note that the factorization of the operator is maintained, although
    the number of projectors is quadrupled (one \f$\phi\f$, three 
    \f$\phi_{\mathbf r}\f$).

    \author Christoph Freysoldt, freyso@fhi-berlin.mpg.de */
class SX_EXPORT_DFT SxPerturbK
{
   public:
      /// Constructor
      SxPerturbK ();
      /// Destructor
      ~SxPerturbK () {}

   protected:
      /// Known G-bases
      SxArray<const SxGBasis *> gBases;
      /// Non-local pseudopotential projectors ik:(ig:iOrb)
      SxArray<PsiGI> phiNl;
      /// Non-local pseudopotential derivative projectors ik:idir:(ig:iOrb)
      SxArray<SxArray<PsiGI> > phiNlR;
      /// Kleinman-Bylander energies is:l:n
      SxArray<SxArray<SxArray<double> > > eKB; 
      /** Quantum numbers (is,ia,n,l,m) for all projectors 
        \note i is used for the reference orbital index
        */
      SxArray<SxQuantumNumbers> orbInfo;

      /// Include contributions from non-local pseudopotential?
      bool nonLocal;

      /// Internal timer
      mutable SxTimer timer;

   public:

      /// Set non-local pseudopotential
      void set (const SxPseudoPot &psPot, const SxGkBasis &gk,
                const SxAtomicStructure &structure);

      /// Print internal timer
      void printTimer ();

      /** \brief <psiL|p|psiR> matrix elements
        
          If the #set routine has not been called, this returns the
          kinetic contribution only.
          
          If #set has been called, also the non-local pseudopotential in
          its Kleinman-Bylander form contributes to the matrix element.

          \param psiL  wavefunctions, nL columns
          \param psiR  wavefunctions, nR columns
          \return a (nL * nR) x 3 matrix. Each column is a
                  (nL x nR) matrix
        */
      SxDiracMat<TPrecCoeffG> 
      getMatrixElements (const SxGBasis::TPsi &psiL, 
                         const SxGBasis::TPsi &psiR) const;

      /** \brief Block size for projector matrix block algorithms */
      int blockSize;

};

#endif /* _SX_PERTURB_K_H_ */
