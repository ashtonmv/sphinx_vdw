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


#ifndef _SX_PW_ATOM_ORBS_H_
#define _SX_PW_ATOM_ORBS_H_

#include <SxArray.h>
#include <SxRadBasis.h>
#include <SxPWSet.h>
#include <SxQuantumNumbers.h>
#include <SxAtomicOrbitals.h>
#include <SxGkBasis.h>
#include <SxDirac.h>
#include <SxDFT.h>

/** \brief Project atomic orbitals onto plane-waves \f$
        \langle \bf{G+k} | \Psi_{i,\sigma,\bf{k}} \rangle
        =
        \sum_r 
           \langle \bf{G+k} | r                   \rangle 
           \langle        r | \mu_{i_s,i_a,n,l,m} \rangle.
    \f$

    \b SxMuPW = S/PHI/nX \f$\mu \f$ projected onto plane-waves.

    This class allows atomic orbitals to be used with a plane-wave Hamilton
    operator class ::SxPWHamiltonian. It represents the atomic orbitals 
    \f[
       | \mu_{i_s,i_a,n,l,m} \rangle
    \f]
    with \f$i_s\f$ being the index of the species, \f$i_a\f$ the index of
    the atom (of that species), and \em n, \em l, and \em m being the quantum
    numbers.
    The orbitals initially are given on a radial grid
    \f[
       \langle r | \mu_{i_s,i_a,n,l,m}\rangle
    \f]
    and will be projected on demand onto the plane-wave grid
    \f[
        \langle \bf{G+k} | \Psi_{i,\sigma,\bf{k}} \rangle
        =
        \sum_r 
           \langle \bf{G+k} | r                   \rangle 
           \langle        r | \mu_{i_s,i_a,n,l,m} \rangle.
    \f]
    The plane-wave wavefunctions are represented by a state or band index \em i,
    the spin quantum number \f$\sigma\f$ and the \b k point \em k.
    
    \sa      \ref page_dirac
    \ingroup group_dirac
    \author  Sixten Boeck
    */
class SX_EXPORT_DFT SxMuPW : public SxPWSet
{
   public:
      /** \brief Get reference orbital index corresponding to iOrb
          TODO: make lookup table like SxPWHamiltonian::phiOrbNl
       */
      int getRefIdx (int iOrb) const
      {
         return phiIdx(iOrb);
      }
      
   public:
      /** \brief Constructor. Initialize all/no components

          This constructor does nothing but initializing some members.
          With the input parameter ::initOrbitals it is possible to control
          whether the orbital container is initialized, too.
          Switching off this flag is useful if some of the species
          or angular momentums should be excluded, e.g. disable the
          local component of the non-local pseudopotential in
          SxPWHamiltonian::computeVnl.
          \param in           The sampled \f$\langle r | \Phi_{i_s,l}\f$,
                              <b>Storage order:</b> :iSpecies,:l
          \param radBasis     Reference to \f$|r\rangle\f$
          \param radBasis     Reference to \f$|\bf{G+k}\rangle\f$
          \param initOrbitals if "true": call SxMuPW::addOrbitals for 
                              all components */
      SxMuPW (const SxArray<SxArray<SxRadBasis::TPsi> > &,
              SxConstPtr<SxRadBasis>,
              SxPtr<SxGkBasis>,
              bool initOrbitals = true);

      SxMuPW (const SxAtomicOrbitals &,
              SxPtr<SxGkBasis>,
              bool initOrbitals = true);
      /** \brief Constructor. Specify which components are to be Initialized.

          This constructor adds the orbitals automatically according to
          the provided \em occMu array
          \param in           The sampled 
                              \f$\langle r | \Phi_{i_s,l} \rangle \f$,
                              <b>Storage order:</b> :iSpecies,:l
          \param radBasis     Reference to \f$|r\rangle\f$
          \param occMu        Enables/disables certain \f$(i_s,l)\f$
                              components,
                              <b>Storage order:</b> :iSpecies,:l
          \param radBasis     Reference to \f$|\bf{G+k}\rangle\f$
       */
      SxMuPW (const SxArray<SxArray<SxRadBasis::TPsi> > &,
              SxConstPtr<SxRadBasis>,
              const SxArray<SxArray<bool> > &,
              SxPtr<SxGkBasis>);

      SxMuPW (const SxAtomicOrbitals &,
              const SxArray<SxArray<bool> > &,
              SxPtr<SxGkBasis>);

      /** \brief Destructor */
      virtual ~SxMuPW ();

      
      /** \brief Initialize orbitals manually

        This function can be used if a certain storage order of orbitals is
        required or if some orbitals should be ignored (e.g. 
        SxPWHamiltonian::computeVnl).
        This function can be used only if the constructor had been called
        with the \em initOrbitals parameter being false.
        When all orbitals are added the function ::finalize has to be
        called before the object can be used for further computaions.

        \sa     SxMuPW::finalize 
        */
      void addOrbital (int iSpecies, int l);

      /** \brief Optimize the internal storage of the orbitals.

        This function has to be called only if the orbitals were added
        manually (::addOrbital). 
        Internally the list structures are remapped to faster arrays which
        are used in further function calls.

        \sa addOrbital */
      void finalize ();
      
      /** \brief computes \f$ \langle \bf{G+k} | \mu_i\rangle\f$
         
          Projects an atomic orbital \f$|\mu_i\rangle\f$ onto plane-waves.
          The index i corresponds to the band index or to the n-tupel
          \f$(i_s,i_a,n,l,m)\f$, resp. 

          \todo Caching for mu(i,iSpin,ik) -> <Psi|H|Psi> */
      virtual const SxGBasis::TPsi operator() (int i, int iSpin, int ik) const;

      /** \brief computes \f$ \langle \bf{G+k} | \mu_i\rangle\f$
         
          Projects an atomic orbital \f$|\mu_i\rangle\f$ onto plane-waves.
          The index i corresponds to the band index or to the n-tupel
          \f$(i_s,i_a,n,l,m)\f$, resp. 

          \todo Caching for mu(i,iSpin,ik) -> <Psi|H|Psi> */
      virtual SxGBasis::TPsi operator() (int i, int iSpin, int ik)
      {
         return (*(const SxMuPW *)this)(i, iSpin, ik);
      }

      /** \brief Returns the number of orbitals. 
       
          This function returns the number of orbitals
          \f$ |\mu_{i_s,i_a,n,l,m}\rangle\f$ which can be projected
          on the plane-waves.  */
      virtual int getNStates (int =0) const;

      /** \brief Returns the number of spin channels

          The current implementation of SxMuPW does not deal with
          spin at all. This function returns always 1 (hardcoded).
          \todo spin polarization of SxMuPW */
      virtual int getNSpin () const;
      
      virtual int getNk () const;

      /** \brief Set the GkBasis */
      virtual void setGkBasisPtr (SxPtr<SxGkBasis>);

// protected:

      /** \brief Atomic orbitals \f$ 
          \langle r | \mu_{i_s,i_a,n,l,m} \rangle 
          \f$ */
      SxAtomicOrbitals mu;  

      /** \brief non-translated projected orbitals

          This list contains the projected orbitals 
          \f[
             \langle \bf{G+k} | \Phi_{i_s,n,l,m} \rangle
             = 
             \sum_r \langle \bf{G+k} | r \rangle
                    \langle r | \mu_{i_s,n,l,m} \rangle
          \f]
          Note, that they are not yet translated to the atomic positions.
          This translation is done in the ::operator() function on the fly */
      SxArray<SxArray<SxDiracVec<TPrecCoeffG> > >  muPhi; //:ik,:iOrb

      /** \brief Pointer to the \f$ r \rangle \f$ basis. */
      const SxRadBasis *radBasisPtr;

      /** \brief convert \f$i \rightarrow (i_s,n,l,m)\f$

          Lookup table to convert the band index i to the orbital
          index i'. The orbital index i' can be also refered as
          \f$i_s,n,l,m\f$. Note, that the atomic index \f$i_a\f$ is
          not taken into account in the orbital index.

          \sa   psiOrb */
      SxArray<int>               phiIdx;

      /** \brief convert \f$ i \rightarrow (i_s,i_a,n,l,m)\f$

           Lookup table to convert the band index i to the orbital
           quantum numbers \f$i_s,i_a,n,l,m\f$. Note, that the atomic
           index \f$i_a\f$ has been taken into account.

           \sa phiIdx */
      SxArray<SxQuantumNumbers>  psiOrb;

      void registerMemoryObservers ();

   private:

      bool isFinalized;
      SxList<SxQuantumNumbers> psiOrbList;
      SxList<int>              phiIdxList;
      SxArray<SxList<SxDiracVec<TPrecCoeffG> > >  muPhiList;

};

#endif /* _SX_PW_ATOM_ORBS_H_ */
