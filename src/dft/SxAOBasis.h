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

#ifndef _SX_AO_BASIS_H_
#define _SX_AO_BASIS_H_

#include <SxPrecision.h>
#include <SxBasis.h>
#include <SxGBasis.h>
#include <SxDirac.h>
#include <SxArray.h>
#include <SxGkBasis.h>
#include <SxPtr.h>
#include <SxDFT.h>
#include <SxOverlap.h>
#include <SxPWOverlap.h>
#include <SxAtomicOrbitals.h>



/** This class describes an abstract atomic orbital basis. The form
    of the atomic orbitals is not specified.

    \b SxAOBasis = S/PHI/nX Atomic Orbital Basis

    \ingroup group_dft
    \author  Christoph Freysoldt
    */
class SX_EXPORT_DFT SxAOBasis : public SxBasis
{
   public:

      typedef TAOBasisType              TBasisType;
      typedef SxDiracVec<TBasisType>    TPsi;
      
      /// Return basis type
      virtual SxString getType () const { return "|mu>"; }
      
      /** \brief Caching behaviour for overlap matrices and their inverse
        */
      enum Caching {
         /// undefined behaviour
         Unknown, 
         /// recompute on demand
         Recompute, 
         /// cache one matrix
         CacheCurrentK, 
         /// cache all matrices
         CacheAll };
      
      /// Empty constructor
      SxAOBasis ();

      /** \brief Constructor from pseudowaves
          \sa set
        */
      SxAOBasis (const SxGkBasis &gk, 
                 const SxRadBasis &rad,
                 const SxArray<SxArray<SxDiracVec<TReal8> > > &psi,
                 const SxConstPtr<SxOverlapBase> SPtrIN = SxPtr<SxPWOverlap>::create ());

      /** \brief Constructor from pseudowaves
         @param gk       G+k basis
         @param psi      radial orbitals (:iSpecies)
         @param l-values l-values for these
        */
      SxAOBasis (const SxGkBasis &gk, 
                 const SxArray<SxDiracMat<Double> > &psi,
                 const SxArray<SxArray<int> >       &lPsi);

      /** \brief Constructor from SxAtomicOrbitals
        \sa set
        */
      SxAOBasis (const SxGkBasis &gk, 
                 const SxAtomicOrbitals &ao,
                 const SxConstPtr<SxOverlapBase> SPtrIN = SxPtr<SxPWOverlap>::create ());

      /** \brief Set to pseudowaves
          \param gk  G+k basis
          \param rad radial basis
          \param psi radial wavefunctions is:l:ir
        */
      void set (const SxGkBasis &gk, 
                const SxRadBasis &rad,
                const SxArray<SxArray<SxDiracVec<TReal8> > > &psi,
                const SxConstPtr<SxOverlapBase> SPtrIN = SxPtr<SxPWOverlap>::create ());
      /** \brief Single-species single-k "light" constructor */
      SxAOBasis (const SxDiracVec<Complex16> &YlmGl, int iSpecies);

   protected:
      /// Initialize timer
      void init ();
      /// Compute reference orbitals
      void computeRefOrbitals(const SxGkBasis &gk,
                              const SxRadBasis &rad,
                              const SxArray<SxArray<SxDiracVec<TReal8> > > &psiRad);
   public:
      /// Caching behaviour for reference orbitals
      enum Caching cacheRefOrb;
   protected:
      /// k-point that is cached for reference orbitals
      mutable int refOrbCachedK;
   public:
      /// Pointer to Overlap Operator
      SxConstPtr<SxOverlapBase> SPtr;


      /// Destructor
      virtual ~SxAOBasis ();

      REGISTER_PROJECTOR (SxAOBasis, SxAOBasis, identity);
      REGISTER_PROJECTOR (SxAOBasis, SxGBasis, toPWBasis);

      /** \brief Identity projection, does not do anything.
        \note This allows to write
        \code
psiAO = (ao | psi);
        \endcode
        even if psi is in an AO basis already.
        */
      SxDiracVec<TBasisType> identity  (const SxAOBasis *,
                                        const SxDiracVec<TBasisType> &) const;
      /**
        \brief Project to plane-wave basis

       */
      SxDiracVec<TGBasisType> toPWBasis (const SxGBasis *,
                                         const SxDiracVec<TBasisType> &) const;

      /**
        \brief Project from plane-wave basis

        \note The projection routine is located here, SxGBasis::toAO just 
              calls this function.

       */
      SxAOBasis::TPsi fromPWBasis (const SxDiracVec<TBasisType> &psiG) const;

      /**
        \brief Project from plane-wave basis

       Projection from SxGBasis or SxPAWBasis
       */
      SxAOBasis::TPsi fromPWBasis (const PsiG &psiS, int ik) const;

      /** \brief Derivative of orbital projections with respect to atomic
                 position
        */
      SxArray<SxDiracMat<Complex16> > gradProject (const PsiG &psi) const;

      /** \brief Derivative of orbital projections with respect to atomic
                 position
          @param psiS  wavefunction with this basis' overlap operator applied
                      (if applicable)
          @param ik   Gk index of reference orbitals to use

          @note This routine is the interface for projections from a SxPAWBasis
                in the absence of an overlap operator in the AO basis. In this
                case, we simply ignore the PAW projection coefficients behind
                the G+k coefficients.

                The routine does NOT check that the PAW basis matches the
                G+k basis of the reference orbitals.

        */
      SxArray<SxDiracMat<Complex16> >
      gradProject (const PsiG &psiS, int ik) const;

      /** \brief The reference orbitals in a |G+k>-Basis ik:is:(ig,iOrb)
          
          For each species, the orbitals are stored once and shifted to
          the appropiate atomic position whenever needed. This reduces
          the necessary amount of memory.
          The iOrb index is a condensed index (n,l,m) and can be decrypted
          with refOrbMap
        */
      SxArray<SxArray<SxDiracMat<TGBasisType> > > refOrbitals;

      /** Additional phase per orbital */
      SxArray<SxDiracVec<Complex16> > extraPhase;

      /**
        \brief This is a container class for atomic orbital indices n,l,m
        */
      class AoIndex  {
         public:
            int n, l, m;
            AoIndex () : n(-1), l(-1), m(-1) {}
            AoIndex (int n_, int l_, int m_) : n(n_), l(l_), m(m_) {}
            ~AoIndex () { /* empty */ }
            // The standard copy and assignment operators work fine
            // because there's no memory management here.
            // i.e. these two are automatically defined:
            // AoIndex (const AoIndex);
            // operator= (const &AoIndex);
            
            bool operator== (const AoIndex &in) const
            {
               return (in.n == n && in.l == l && in.m == m);
            }
            bool operator!= (const AoIndex &in) const
            {
               return (in.n != n || in.l != l || in.m != m);
            }
            /** \brief Ordering comparison (hierarchy 1) n 2) l 3) m).

              \note
              (n1, l1, m1) < (n2, l2, m2) means
              -# n1 != n2 => n1 < n2
              -# n1 == n2 && l1 != l2 => l1 < l2
              -# n1 == n2 && l1 == l2 => m1 < m2
              */
            bool operator< (const AoIndex &in) const
            {
               if (n != in.n) return (n < in.n);
               if (l != in.l) return (l < in.l);
               return (m < in.m);
            }
            /** \brief Ordering comparison (hierarchy 1) n 2) l 3) m).

              \note
              (n1, l1, m1) > (n2, l2, m2) means
              -# n1 != n2 => n1 > n2
              -# n1 == n2 && l1 != l2 => l1 > l2
              -# n1 == n2 && l1 == l2 => m1 > m2
              */
            bool operator> (const AoIndex &in) const
            {
               if (n != in.n) return (n > in.n);
               if (l != in.l) return (l > in.l);
               return (m > in.m);
            }
      };
/** \brief Maps iRefOrb to (n,l,m) for each species (is:iRefOrb)

          \example
          \code
SxAOBasis::AOIndex nlm = aoBasis.refOrbMap(is)(iRefOrb);
int n = nlm.n;
int l = nlm.l;
int m = nlm.m;
\endcode
        */
      SxArray<SxArray<AoIndex> > refOrbMap;

      /**
        \brief This is a container class for atomic orbital indices is,ia,io
        */
      class OrbitalIndex  {
         public:
            int is, ia, io;
            OrbitalIndex () : is(-1), ia(-1), io(-1) {} 
            OrbitalIndex (int is_, int ia_, int io_) 
              : is(is_), ia(ia_), io(io_) {}
            ~OrbitalIndex () {/* empty */}
            // The standard copy and assignment operators work fine
            // because there's no memory management here.
            // i.e. these two are automatically defined:
            // OrbitalIndex (const OrbitalIndex);
            // operator= (const &OrbitalIndex);
            
            bool operator== (const OrbitalIndex &in) const
            {
               return (in.is == is && in.ia == ia && in.io == io);
            }
            bool operator!= (const OrbitalIndex &in) const
            {
               return (in.is != is || in.ia != ia || in.io != io);
            }
            /** \brief Ordering comparison (hierarchy 1) is 2) ia 3) io).

              \note
              (is, ia, io) < (is2, ja, jo) means
              -# is != js => is < js
              -# is == js && ia != ja => ia < ja
              -# is == js && ia == ja => io < jo
              */
            bool operator< (const OrbitalIndex &in) const
            {
               if (is != in.is) return (is < in.is);
               if (ia != in.ia) return (ia < in.ia);
               return (io < in.io);
            }
            /** \brief Ordering comparison (hierarchy 1) is 2) ia 3) io).

              \note
              (is, ia, io) > (js, ja, jo) means
              -# is != js => is > js
              -# is == js && ia != ja => ia > ja
              -# is == js && ia == ja => io > jo
              */
            bool operator> (const OrbitalIndex &in) const
            {
               if (is != in.is) return (is > in.is);
               if (ia != in.ia) return (ia > in.ia);
               return (io > in.io);
            }
      };

/** \brief Maps iOrb to (is,ia,io) for each orbital (mu)

          \example
          \code
int mu;
SxAOBasis::OrbitalIndex sao = aoBasis.orbitalMap(mu);
SxAOBasis::AoIndex nlm = aoBasis.refOrbMap(sao.is)(sao.io);
int is = sao.is;
int ia = sao.ia;
int n = nlm.n;
int l = nlm.l;
int m = nlm.m;
\endcode
*/
      SxArray<OrbitalIndex> orbitalMap; 

      int getNSpecies () const { return int(refOrbMap.getSize ()); }
      int getNOrb () const { return int(orbitalMap.getSize ()); }
      /// Return number of atomic orbitals
      virtual ssize_t getNElements () const  { return getNOrb (); }
      int getLMax () const;

      /// Blocksize for blocked projections
      int blockSize;
/** \brief Return all orbitals in |G+k> basis
        */
      SxDiracMat<TGBasisType> getAOinG (int ik) const;
      /** \brief Return a single orbital in |G+k> basis
        */
      SxDiracVec<TGBasisType> getAOinG(int ik, int iOrb) const;

   protected:
      //\name Overlap 
      //\{
      /** \brief The overlap matrices ik:(mu,nu)

          This stores the matrices
          \f[ S^k_{\mu,\nu} = \int dr^3 (\chi^k_{\mu})^*(r) \chi^k_nu(r) \f]
        
        */
      mutable SxArray<SxDiracMat<TGBasisType> > overlap;
      /// Caching behaviour for overlap matrices
      enum Caching cacheOverlap;
      /// k-point that is cached for overlap
      mutable int overlapCachedK;

      /** \brief The overlap matrices ik:(mu,nu)

          This stores the matrices
          \f[ S^k_{\mu,\nu} = \int dr^3 (\chi^k_{\mu})^*(r) \chi^k_nu(r) \f]
        
        */
      mutable SxArray<SxDiracMat<TGBasisType> > invOverlap;
      /// Caching behaviour for inverse overlap matrices
      enum Caching cacheInverse;
      /// k-point that is cached for inverse overlap
      mutable int invOverlapCachedK;
      
   public:
      /// Get overlap matrix for k-point ik
      SxDiracMat<TAOBasisType> getOverlap (int ik) const;
      /// Get inverse overlap matrix for k-point ik
      SxDiracMat<TAOBasisType> getInverseOverlap (int ik) const;
      /**  \brief Compute the overlap matrix
        \note You may consider caching the overlap matrices instead
              of computing them whenever needed.
        \sa setOverlapCaching, getOverlap 
        */
       SxDiracMat<TAOBasisType> calculateOverlap (int ik) const;


      /** \brief Set caching behaviour for overlap matrices
          @param mode the mode
          @param nk number of k-points to cache if mode is CacheAll
                 and k-dependent reference orbitals are not cached.
        */
      void setOverlapCaching (enum Caching mode, int nk = -1);
      /** \brief Set caching behaviour for inverse overlap matrices
          @param mode the mode
          @param nk number of k-points to cache if mode is CacheAll
                 and k-dependent reference orbitals are not cached.
       */
      void setInvOverlapCaching (enum Caching mode, int nk = -1);
      //\}

};

#endif /* _SX_AO_BASIS_H_ */

