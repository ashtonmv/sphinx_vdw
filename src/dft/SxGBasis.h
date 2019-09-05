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

#ifndef _SX_G_BASIS_H_
#define _SX_G_BASIS_H_

#include <SxBasis.h>
#include <SxMatrix.h>
#include <SxList.h>
#include <SxArray.h>
#include <SxPrecision.h>
#include <SxDirac.h>
#include <SxBinIO.h>
#include <SxFFT3d.h>
#include <SxFFT2d1d.h>
#include <SxTypes.h>
#include <SxAtomicStructure.h>
#include <SxSymbolTable.h>
#include <SxYlm.h>
#include <SxRadBasis.h>
#include <SxDFT.h>

class SxRBasis;

/** \brief | \b G> basis

  \b SxGBasis = S/PHI/nX \f$ | \mathbf{G} \rangle \f$ basis

  This class represents the |G> basis. In particulary it defines how the
  projections from |G> to the realspace |R> can be performed by applying
  fast Fourier transformation (FFT).

  \ingroup group_dft
  \ingroup group_dirac
  \author  Sixten Boeck
  */
class SX_EXPORT_DFT SxGBasis : public SxBasis
{
   public:

      typedef TGBasisType                   TBasisType;
      typedef SxDiracVec<TPrecCoeffG>       TPsi;
      typedef SxDiracVec<TPrecCoeffG>       TPsiGI;
      typedef SxDiracVec<TPrecRhoG>         TRho;

      /** \brief list of all \b G vectors, {:ig,:{xyz}}

          Contains the set of |G> vectors. They will be computed in the
          ::compute routine as
          \f[
             \mathbf{G} = \mathbf{B}^t
                           \left( \begin{matrix}
                              i \\
                              j \\
                              k
                           \end{matrix} \right)
                        + \Delta \mathbf{G}
          \f] with
          \f[
             B = 2\pi \mathbf{A}^{-1}
          \f]
          Note, that the \b B matrix contains the reciprocal unit vectors
          in columns.
          \sa g2
       */
      SxDiracVec<TPrecG>                    gVec;              // :ig,:3
      /** \brief list of \f$ |\mathbf{G}|^2 \f$, {:ig}

           This array contains the absolute squared \b G vector values, i.e.,
           \f[
              |\mathbf{G}|^2 = \langle \mathbf{G} | \mathbf{G} \rangle
           \f] */
      SxDiracVec<TPrecG>                    g2;
      /** \brief structure factors, {:is,:ig}

           This array contains the structure factors for all species. They
           read
           \f[
              S_{i_s}(\mathbf{G}) = \sum_{i_a(i_s)}
                                     e^{-i \mathbf{G} \cdot \tau_{i_s,i_a})}
           \f]
           They will be computed/updated in the routine ::changeTau.
           \sa phaseFactors
       */
      SxArray<SxDiracVec<TPrecPhase> >      structureFactors;  // :is,:ig
      /** \brief Get phase factors for one atom

        \param is species id
        \param ia atom id

        This returns the phase factors for one atom. They read
        \f[
           P_{i_s,i_a}(\mathbf{G}) = e^{-i \mathbf{G} \cdot \tau_{i_s,i_a})}
        \f]
        Depending on #memMode the phase factors are either
        computed in changeTau, or they are calculated on demand.
        \sa phaseFactors
       */
      SxDiracVec<TPrecPhase> getPhaseFactors(int is, int ia) const;

      /** \brief Get phase factors for one atom and one G vector

        \param is species id
        \param ia atom id
        \param ig which G vector

        Depending on #memMode the phase factors are either
        computed in changeTau, or they are calculated on demand.
        \sa phaseFactors
       */
      SxComplex16 getPhaseFactors(int is, int ia, ssize_t ig) const;

      /** If set to true, 3D phase factors are stored. Otherwise, only
          partial 1D phase factors are stored, and the 3D ones are computed
          on demand.
          \brief Whether to save the 3D phase factors
          \sa phaseFactors
      */
      enum MemoryUsage { SaveTime, SaveMemory } memMode;
      /// This cleans the phase factors
      void cleanPhaseFactors () const;
   protected:
      /** \brief phase factors, {:is,:ia,:ig}

           This array contains the phase factors for all atoms. They
           read
           \f[
              P_{i_s,i_a}(\mathbf{G}) = e^{-i \mathbf{G} \cdot \tau_{i_s,i_a})}
           \f]
           They will be computed/updated in the routine ::changeTau.
           \sa structureFactors

           If memMode is set to 'SaveMemory' phaseFactors contains
           the 1dimensional phase factors in a packed manner, namely
           - 0     ... nx-1         for 1st reciprocal lattice vector.
           - nx    ... nx+ny-1      for 2nd reciprocal lattice vector.
           - nx+ny ... nx+ny+nz-1   for 3rd reciprocal lattice vector.

           They are NOT computed in changeTau, but on demand by
           getPhaseFactors. The memory can be freed by cleanPhaseFactors.

           The calculation on demand makes use of the factorization
           \f[ e^{-i\vec G \cdot \vec r}
               = e^{-iG_1 \cdot r_1}\cdot e^{-iG_2 \cdot r_2}
                 \cdot e^{-iG_3 \cdot r_3}
           \f]
           and only the 1-dimensional factors e^{-i G_i \cdot r_i} are
           stored permanently.
       */
      mutable SxArray<SxArray<SxDiracVec<TPrecPhase> > > phaseFactors;
      /** \brief Packed version of fft3d.mesh.getMeshVec(n123)

          This contains the packed version of fft3d.mesh.getMeshVec(n123), i.e.
          packedGrel(3*ig + dim) = fft3d(0).mesh.getMeshVec(n123(0)(ig))(dim)
          with dim = {0,1,2}.
        */
      mutable SxArray<ssize_t> packedGrel;
      /// If necessary, set up 1D phFac(is)(ia) and packedGrel
      void setupPhase1D (int is, int ia) const;
      /// Set up packed 1D phase factors for a given coordinate
      SxDiracVec<Complex16> setupPhase1D (const Coord &tau) const;
      /// Setup mapping from ig to packed 1D phase factors
      void setupPackedG () const;
   public:
      /** The index array to build up the fft mesh :iFFT:ig.

        Entities living on the |G> space are saved in a compressed way.
        Before performing the FFT they have to be inflated to the entire FFT
        box. In case of transformations from |R> to |G> the data in the
        full FFT box will be compressed again.
        This mapping is done using the n123 index array. It will be initialized
        in the ::compute routine and is being kept fixed during the entire
        calculation.
        Note: for each registered R-basis, there is one n123

        \sa     compute
        \sa     n123inv
        \author Sixten Boeck
       */
      mutable SxArray<SxVector<TPrecFFTIdx> > n123;

      /** \brief The inverse of the first n123 index array.

          Input: fft mesh index. Output: index in this G (or G+k) basis array.
          If the requested fft mesh idx is not contained in this G (or G+k)
          basis array, a -1 is being returned. In other words: -1 is being
          returned, if and only if G or G+k lies outside the energy cutoff
          radius.

          Before this table is ready to use, you have to run ::n123invSetup()
          first.

          \sa     n123
          \sa     n123invSetup()
          \author Matthias Wahn
          */
      mutable SxDiracVec<Int>                n123inv;

      /** \brief Number of \b G vectors */
      int                                    ng;

      /** \brief Number of components */
      int nComp;

      /** \brief setup the real Ylm (not normalized) */
      void setupRealYlm (int lmax) const;

      /** \brief Get the real Ylm (not normalized) */
      SxDiracVec<Double> getYlm(int l, int m) const;

   protected:
      /// The dual R-basis
      mutable const SxRBasis* rBasisPtr;

      /// The cashed Ylm
      mutable SxDiracMat<Double> realYlm;

      /// Pointer to 2D+1D for convolutions etc.
      SxPtr<SxFFT2d1d> mixedFFT;

   public:
      /// Has this G-basis a mixed 2+1 D FFT?
      /// Map mesh
      SxVector<TPrecFFTIdx> n231;

      bool hasMixedFFT () const { return mixedFFT; }

      /// Has this G-basis a mixed 2+1 D FFT for a given mesh?
      bool hasMixedFFT (const SxMesh3D &mesh) const {
         return mixedFFT && mixedFFT->realMesh == mesh;
      }

      /// Get ptr to mixed 2+1 D FFT
      const SxPtr<SxFFT2d1d>& getMixedFFT () const { return mixedFFT; }

      /** Register the dual R basis */
      void registerRBasis (const SxRBasis &rBasis) const
      {
         registerBasis (rBasis);
         rBasisPtr = &rBasis;
      }

      /** Get the dual R basis */
      const SxRBasis& getRBasis () const
      {
         SX_CHECK(rBasisPtr);
         return *rBasisPtr;
      }

      /// Setup the mixed FFT for fast convolutions etc.
      void setupMixedFFT (const SxRBasis &rBasis);

      /// Convolute
      PsiGI convolute (const PsiGI &allPsi,
                       const SxDiracVec<Double> &V) const;

      /// 2+1D FFT algorithm for density calculation
      void addToRho (const PsiGI &allPsi,
                     const SxDiracVec<Double> &focc,
                     SxDiracVec<Double> *rho) const;

      /// Returns number of G vectors
      virtual ssize_t getNElements () const { return ng * nComp; }
      /** \brief FFT transformator

          This object is used to transform entities from \mathbf G space
          to realspace SxRBasis */
      mutable SxArray<SxFFT3d>              fft3d;
      /** \brief Atomic structure

          The atomic coordinates have to be known in order to compute both
          the phase and the structure factors. */
      const SxAtomicStructure               *structPtr;



      SxGBasis (MemoryUsage memModeIn = SaveTime);
      /** Construct GBasis from existing one with new number of components */
      SxGBasis (const SxGBasis &, int);
      /** Construct a G basis with a offset vector, e.g. for setting
          up the |G+k> basis or |G+k+q> basis.
          \param dG   translation vector */
      SxGBasis (const SxVector3<Int> &mesh, const SxAtomicStructure &,
                PrecEnergy gCut,
                const SxVector3<TPrecG> &dG = SxVector3<TPrecG>(0.,0.,0.),
                bool useTimeRevSymIn=true,
                MemoryUsage memModeIn = SaveMemory);
      void set (const SxVector3<Int> &mesh, const SxAtomicStructure &,
                PrecEnergy gCut,
                const SxVector3<TPrecG> &dG = SxVector3<TPrecG>(0.,0.,0.),
                bool useTimeRevSymIn=true,
                MemoryUsage memModeIn = SaveMemory);
 
      /** Destructor */
      virtual ~SxGBasis ();

      /// Initialize G basis
      void set (const SxVector3<Int> &mesh,
                const SxCell &cell,
                PrecEnergy gCut,
                const SxVector3<TPrecG> &dG = SxVector3<TPrecG>(0.,0.,0.),
                bool useTimeRevSymIn=true,
                MemoryUsage memModeIn = SaveTime);

      SxGBasis& operator=(const SxGBasis &) {SX_EXIT;}
      /** \brief Computes |G>, |G|^2, and n123

        This routine prepares the G Basis by setting up the FFT mapping
        indices.
        It creates all those |G> vectors
        \f[
           \mathbf{G} = i \mathbf{a}_1 + j \mathbf{a}_2 + k \mathbf{a}_3
                      + \Delta \mathbf{G}
        \f]
        with \f$ \Delta \mathbf{G}\f$ is just any translation vector (usually
        zero, for |G+k> spaces it's the \b k vector, SxGkBasis). The indices
        \em i, \em j, and \em k vary from \f$ -\frac{n_{xyz}}{2}\f$ to
        \f$ \frac{n_{xyz}}{2} \f$. Only |G> vectors inside the cutoff
        sphere
        \f[
           |\mathbf{G}|^2 \leq g_{\mathrm{cut}}
        \f]
        are included.
        In order to safe memory the |G> vectors are not stored in the full
        FFT box. Only the sphere is considered because everything outside the
        cutoff sphere would vanish anyway. Hence, all entities given in
        |G> space are saved in a compressed fashion (see also ::n123) in
        a pseudo-one-dimensional array.
        This array samples all grid points inside the cutoff sphere
        like
           \image html  fft-1.png
           \image latex fft-1.png
        A further (implicit) transformation is required because the FFT grid
        ranges from 0 to mesh instead of -mesh/2 to +mesh/2. A simple
        translation by +mesh/2 would introduce a wrong frequency doubling.
        The implemented way for this "translation" is as demstrated here:
        to be continued...
           \image html  fft-2.png
           \image latex fft-2.png
        The original mesh is centered arround the origin with the quadrants
        1, 2, 3, and 4. The quadrands are mapped by the modulo function to the
        FFT input mesh 1', 2', 3', and 4'. This scheme works due to the
        sampling theorem.
        \author Sixten Boeck */
      void  compute ();

      /// Set number of components
      void setNComp (int comp);
      /// Get number of components
      virtual int getNComp () const;

      // make also general registration visible in present scope
      using SxBasis::registerBasis;

      /** Register an additional R basis */
      void registerBasis (const SxRBasis &) const;

      /** Register basis of unknown type */
      virtual void registerUnknown(const SxBasis &basis) const;

      /** Special code for deregistering an R basis */
      virtual void deregister (const SxBasis *basis) const;

      /// Change the FFT mesh size
      void replaceMesh (const SxFFT3d &newFFT);

      /// Get basis type
      virtual SxString getType () const  { return "|G>"; }
   protected:
      /// Set up the G-vectors
      void  update ();
      /// Set up the G-vectors if structPtr not available
      void  update (const SxCell &cell);
   public:
      /** \brief Updates structure and phase factors

          This routine has to be called whenever the atomic positions
          have changed. */
      void  changeTau (const SxAtomicStructure &tauList);

      /**
        \brief Calculate phase factors, is:(ig:ia)
        \param coord the structure
        \return an array (nSpecies) of matrices (ng x nAtoms(is))
               containing phase factors \f$e^{-iGr}\f$
        */
      SxArray<SxDiracMat<TPrecPhase> >
      getPhaseFactors (const SxAtomicStructure &structure) const;

      /** \brief Calculate phase factor for one coordinate
        \param tau coordinate
        \return phase factors \f$e^{-iG\tau}\f$
      */
      SxDiracVec<Complex16> getPhaseFactors(const Coord &tau) const
      {
         return composePhase(setupPhase1D(tau));
      }

      /** Compose 3D phase factors from packed 1D phase factors
        */
      SxDiracVec<Complex16>
      composePhase (const SxDiracVec<Complex16> &phFac) const;

      /** Multiply with 3D phase factors from packed 1D phase factors
          @param phFac   packed 1D phase factors
          @param resPtr  vector to be multiplied
        */
      void applyComposedPhase (const SxDiracVec<Complex16> &phFac,
                               SxDiracVec<Complex16> *resPtr) const;

      SxAtomicStructure get1DPackedVecs () const;

      /** \brief extract \c eCut from the input file */
      static Real8 getECut (const SxSymbolTable *);
      /** \brief extract the FFT mesh data \c mesh from the input file.

          The returned mesh is commensurable. */
      static SxVector3<Int> getMesh (const SxSymbolTable *);

      /** \brief extract the FFT mesh data \c mesh from the input file.

          The returned mesh is commensurable with the symmetries
          provided by the cell. */
      static SxVector3<Int> getMesh (const SxSymbolTable *,
                                     const SxCell &);

      /// Calculate the optimal mesh for a certain cut-off energy
      static SxVector3<Int> getMeshSize (Real8 eCut,
                                         const SxMatrix3<TPrecTauR> &aMat,
                                         Real8 meshAccuracy = 1.0);
      /** \brief Get a commensurable FFT mesh size

          This function computes FFT mesh size by calling
          SxGBasis::getMeshSize according to a given energy-cutoff and the
          given super cell dimensions. Before returning it varifies that
          the mesh is commensurable (SxGBasis::isCommensurableMesh).
          If this test fails it tries to find alternative FFT mesh which
          reflects the given cell geometry. */
      static SxVector3<Int>
         getCommensurableMesh (Real8 eCut, const SxCell &,
                               double meshAccuracy=1.0);
      /** Find symmetry commensurable mesh
        @param cell    unit cell with symmetries
        @param meshIn  starting point for mesh search
        */
      static SxVector3<Int> getCommensurableMesh (const SxCell &cell,
                                                  SxVector3<Int> meshIn);

      /** \brief Checks whether a mesh is commensurable

          SxGBasis::isCommensurableMesh checks if a given mesh reflects
          the exsiting cell geometry and symmetry operations. */
      static bool isCommensurableMesh (const SxVector3<Int> &,
                                       const SxCell &,
                                       bool verbose=false);

      /** Returns a copy of the specified g-vector. In time crucial
          part get the entire vector rather than using this function. */
      SxVector3<TPrecG> getG (int) const;

      SxVector3<TPrecG> getK () const;

      REGISTER_PROJECTOR (SxGBasis, SxGBasis, identity);
      REGISTER_PROJECTOR (SxGBasis, SxRBasis, toRealSpace);
      REGISTER_PROJECTOR (SxGBasis, SxRadBasis, toRadialSpace);
      REGISTER_PROJECTOR (SxGBasis, SxAOBasis, toAO);
      REGISTER_PROJECTOR (SxGBasis, SxPartialWaveBasis, toPartials);
      REGISTER_PROJECTOR (SxGBasis, SxPAWBasis, toPAW);
      REGISTER_PROJECTOR (SxGBasis, SxBasis, toAny);

      /** \brief Identity projector from \b G to \b G

          This is the identity projector
          \f$ \langle \mathbf{G} | \mathbf{G} \rangle \f$
          \sa \ref page_dirac */
      SxDiracVec<TBasisType> identity (const SxGBasis *,
                                       const SxDiracVec<TBasisType> &) const;
      /** \brief Projector from \b G space to realspace

          This is the projector
          \f$ \langle \mathbf{R} | \mathbf{G} \rangle \f$
          \sa \ref page_dirac */
      SxDiracVec<TBasisType> toRealSpace ( const SxRBasis *,
                                           const SxDiracVec<TGBasisType> &
                                         ) const;
      /** \brief Projector from \b G space to radial Space

          This is the projector
          \f$ \langle \mathbf{r_{is,ia,n,l,m}} | \mathbf{G} \rangle \f$
          \sa /ref page_dirac */
      SxDiracVec<TBasisType> toRadialSpace ( const SxRadBasis *,
                                             const SxDiracVec<TGBasisType> &
                                           ) const;
      /** \brief Projector from \b G space to atomic orbitals

          This is the projector
          \f$ \langle \mathbf{\mu} | \mathbf{G} \rangle \f$
          \sa \ref page_dirac */
      SxDiracVec<TAOBasisType>
      toAO (const SxAOBasis *aoBasis_,
            const SxDiracVec<TGBasisType> &psiG) const;

      /** \brief Projector (from \b G space to partial waves

          This is the projector
          \f$ \langle \mathbf{p} | \mathbf{G} \rangle \f$
          \sa \ref page_dirac */
      SxDiracVec<Complex16> toPartials (const SxPartialWaveBasis *,
                                        const SxDiracVec<TGBasisType> &) const;

      /** \brief Projector from \b G space to full PAW space

          This is the projector
          \f$ \langle G, \mathbf{p} | \mathbf{G} \rangle \f$
          \sa \ref page_dirac */
      SxDiracVec<Complex16> toPAW (const SxPAWBasis *,
                                   const SxDiracVec<TGBasisType> &) const;

      /** \brief Projector from \b G space to some anonymous basis

          This is the versatile interface: we try to dynamically determine
          the basis to project to.

          \sa \ref page_dirac */
      SxDiracVec<Complex16> toAny (const SxBasis *,
                                   const SxDiracVec<Complex16> &) const;

      /** \brief Laplacian operator in \b G space

          The Laplacian in \b G is
          \f[
             L \Psi(\mathbf{G}) = |\mathbf{G}|^2 |\Psi(\mathbf{G})|^2
          \f]
          \sa \ref page_dirac */
      virtual Real8 laplacian (const void *) const;

      /**
        \author Abdullah Al-Sharif
        \deprecated see list of add-ons
        */
      SxDiracVec<TPrecCoeffR>       GtoR111   (const SxDiracVec<TPrecCoeffG> &,
                                               int nr, int shift=0) const;

      /**
        \author Abdullah Al-Sharif
        \deprecated see list of add-ons

        Used for debugging in the gsEXX project.
        */
      SxVector<TPrecCoeffR>
         GtoR111 (const SxDiracVec<TPrecCoeffG> &,
                  int nr, int shift, const SxCell &cell) const;


      /** \brief Sets up the fft -> G(+k) basis index array. */
      void n123invSetup () const;


      /** \brief Yields the index of the component of a sum of two indeces
                 given this |G+k> basis.

          Be \f$ i_{1} \f$ the index of \f$ G_{1} \f$ and \f$ i_{2} \f$ the
          index of \f$ G_{2} \f$ in this \f$ |\mathbf{G}\rangle \f$ or
          \f$ |\mathbf{G}+\mathbf{k}\rangle \f$ basis set. The output of this
          function is then the index of the
          \f$ (\mathbf{G_{1}}+\mathbf{G_{2}}) \f$ component.

          \note   Before you can use this function, the n123inv table must
                  be initialised. Therefore use SxGBasis::n123invSetup ().

          \sa     SxGBasis::getIdxGSum (idxGp, *gpPtr, idxGq, *gqPtr),
                  SxGBasis::getIdxGDiff (idxG1, idxG2),
                  SxGBasis::getIdxGDiff (idxGp, *gpPtr, idxGq, *gqPtr)

          \author Matthias Wahn
          */
      inline int getIdxGSum (int idxG1, int idxG2) const
      {
         SX_CHECK (n123inv.getSize() > 0, n123inv.getSize());

         // get idx in 3-dim. fft mesh
         const PrecFFTIdx fftIdx1 = n123(0)(idxG1);
         const PrecFFTIdx fftIdx2 = n123(0)(idxG2);

         // get mesh vector with coordinates -mesh(i)/2 < vec(i) <= mesh(i)/2
         const SxVector3<Int> idxVec1 
            = fft3d(0).mesh.getMeshVec (fftIdx1, SxMesh3D::Origin);
         const SxVector3<Int> idxVec2 
            = fft3d(0).mesh.getMeshVec (fftIdx2, SxMesh3D::Origin);

         // add index vectors
         const SxVector3<Int> idxVecRes = idxVec1 + idxVec2;

         // since the range of this vector has doubled, cut it pursuant to
         // the fft mesh
         const SxVector3<Int> mesh = fft3d(0).mesh;
         int i;
         for (i = 0; i < 3; i++)
            if (2 * idxVecRes(i) <= -mesh(i) || 2 * idxVecRes(i) > mesh(i))
               return -1;

         // if vector is in fft range, get its index in this G (or G+k) basis
         const ssize_t fftIdxRes 
            = fft3d(0).mesh.getMeshIdx (idxVecRes, SxMesh3D::Origin);
         return n123inv(fftIdxRes);
      }

      /** \brief Yields the index of the component belonging to the sum
                 of two indeces, the latter may belong to arbitrary
                 \f$ | \mathbf{G} + \mathbf{p},\mathbf{q} \rangle \f$
                 basis sets.

          Be \f$ i_{1} \f$ the index of \f$ G_{1} \f$ in the
          \f$ |\mathbf{G}+\mathbf{p}\rangle \f$ basis set and
          \f$ i_{2} \f$ the index of \f$ G_{2} \f$ in the
          \f$ |\mathbf{G}+\mathbf{q}\rangle \f$ basis set. The output of this
          function is then the index of the
          \f$ (\mathbf{G_{1}}+\mathbf{G_{2}}) \f$ component in this
          \f$ | \mathbf{G} + \mathbf{k} \rangle \f$ or
          \f$ | \mathbf{G} \rangle \f$ basis set.

          @param   idxGp index of a coefficient in the |G+p> basis
          @param  *gpPtr pointer to the |G+p> basis
          @param   idxGq index of a coefficient in the |G+q> basis
          @param  *gqPtr pointer to the |G+q> basis

          @return  the index of the coefficient of the sum
                   \f$\mathbf{G_{1}} + \mathbf{G_{2}}\f$ of the vectors
                   denoted by the two indeces, given in this |G+k> basis,
                   i.e., the basis of the object where you put the dot and
                   the name of this function behind.

          \note    Don't mix up *gpPtr and *gqPtr with pointers to objects of
                   the SxGkBasis class! *gpPtr and *gqPtr here denote objects
                   of the SxGBasis class.

          \note    Before you can use this function, the n123inv table must
                   be initialised. Therefore use SxGBasis::n123invSetup ().

          \sa      SxGBasis::getIdxGSum (idxG1, idxG2),
                   SxGBasis::getIdxGDiff (idxG1, idxG2),
                   SxGBasis::getIdxGDiff (idxGp, *gpPtr, idxGq, *gqPtr)

          \author  Matthias Wahn
       */
      inline int getIdxGSum (int idxGp, const SxGBasis *gpPtr,
                             int idxGq, const SxGBasis *gqPtr) const
      {
         SX_CHECK     (gpPtr);
         SX_CHECK     (gqPtr);
         SX_CHECK (n123inv.getSize() > 0, n123inv.getSize());

         // get idx in 3-dim. mesh
         const PrecFFTIdx fftIdxGp = gpPtr->n123(0)(idxGp);
         const PrecFFTIdx fftIdxGq = gqPtr->n123(0)(idxGq);

         // get mesh vector with coordinates -mesh(i)/2 < vec(i) <= mesh(i)/2
         const RelVec idxVecGp
            = gpPtr->fft3d(0).mesh.getMeshVec (fftIdxGp, SxMesh3D::Origin);
         const RelVec idxVecGq 
            = gqPtr->fft3d(0).mesh.getMeshVec (fftIdxGq, SxMesh3D::Origin);

         // add index vectors
         const RelVec idxVecRes = idxVecGp + idxVecGq;

         // since the range of this vector may exceed the range of the fft
         // mesh belonging to this |G+k> basis, cut it pursuent to the mesh
         const RelVec mesh = fft3d(0).mesh;
         int i;
         for (i = 0; i < 3; i++)
            if (2 * idxVecRes(i) <= -mesh(i) || 2 * idxVecRes(i) > mesh(i))
               return -1;

         // if vector is in fft range, get its index in the |G+k> basis
         const ssize_t fftIdxRes
            = fft3d(0).mesh.getMeshIdx (idxVecRes, SxMesh3D::Origin);
         return n123inv(fftIdxRes);
      }

      /** \brief Yields the index of the component of a sum of a
                 \f$ \mathbf{G} \f$-vector belonging to a certain index and
                 another \f$ \mathbf{G} \f$-vector given in relative
                 coordinates.

                 Be \f$ i_{\mathbf{G}_1} \f$ the index of \f$ \mathbf{G}_1 \f$
                 in this \f$ |\mathbf{G}(+\mathbf{k})\rangle \f$ basis set.
                 Be \f$ \mathbf{G}_2 \f$ a further vector given in relative
                 coordinates \f$ \mathbf{G}_2 = \mathbf{G}_2^{\rm rel}\,B \f$,
                 with \f$ B \f$ denoting the Matrix of the reciprocal basis
                 vectors (line-wise).
                 The function's output is then the index
                 \f$ i_{{\mathbf{G}_1}+{\mathbf{G}_2}} \f$ of the sum of the
                 vectors. The index belongs to this
                 \f$ |\mathbf{G}(+\mathbf{k})\rangle \f$ basis set.

          \note  Before you can use this function, the n123inv table must
                 be initialised. Therefore use SxGBasis::n123invSetup ().

          \sa    SxGBasis::getIdxGSum (idxG1, idxG2),
                 SxGBasis::getIdxGSum (idxGp, *gpPtr, idxGq, *gqPtr),

          \author Matthias Wahn
       */
      inline int getIdxGSum (int idxG, RelVec &relVec) const
      {
         SX_CHECK (n123inv.getSize() > 0, n123inv.getSize());

         // get idx of 3-dim. fft mesh
         const PrecFFTIdx fftIdx = n123(0)(idxG);

         // get mesh vector with coordinates -mesh(i)/2 < vec(i) <= mesh(i)/2
         const RelVec idxVec
            = fft3d(0).mesh.getMeshVec (fftIdx, SxMesh3D::Origin);

         // add vectors in relative coordinates
         const RelVec idxVecRes = idxVec + relVec;

         // check, if result is still in range
         const SxVector3<Int> mesh = fft3d(0).mesh;
         int i;
         for (i = 0; i < 3; i++)
            if (2 * idxVecRes(i) <= -mesh(i) || 2 * idxVecRes(i) > mesh(i))
               return -1;

         // if vector is in fft range, get its index in this G (or G+k) basis
         const ssize_t fftIdxRes 
            = fft3d(0).mesh.getMeshIdx (idxVecRes, SxMesh3D::Origin);
         return n123inv(fftIdxRes);
      }

      /** \brief Yields the index of the component of the difference of two
                 indeces given this |G+k> basis.

          Be \f$ i_{1} \f$ the index of \f$ G_{1} \f$ and \f$ i_{2} \f$ the
          index of \f$ G_{2} \f$ in this \f$ |\mathbf{G}\rangle \f$ or
          \f$ |\mathbf{G}+\mathbf{k}\rangle \f$ basis set. The output of this
          function is then the index of the
          \f$ (\mathbf{G_{1}}-\mathbf{G_{2}}) \f$ component.

          \note   Before you can use this function, the n123inv table must
                  be initialised. Therefore use SxGBasis::n123invSetup ().

          \sa     SxGBasis::getIdxGSum (idxG1, idxG2),
                  SxGBasis::getIdxGSum (idxGp, *gpPtr, idxGq, *gqPtr),
                  SxGBasis::getIdxGDiff (idxGp, *gpPtr, idxGq, *gqPtr)

          \author Matthias Wahn, C. Freysoldt (low-level)
          */
      inline int getIdxGDiff (int idxG1, int idxG2) const
      {
         SX_CHECK (n123inv.getSize() > 0, n123inv.getSize());

         // get idx in 3-dim. fft mesh
         const PrecFFTIdx fftIdx1 = n123(0)(idxG1);
         const PrecFFTIdx fftIdx2 = n123(0)(idxG2);

         // get mesh vector with coordinates -mesh(i)/2 < vec(i) <= mesh(i)/2
         const SxVector3<Int> idxVec1
            = fft3d(0).mesh.getMeshVec (fftIdx1, SxMesh3D::Origin);
         const SxVector3<Int> idxVec2
            = fft3d(0).mesh.getMeshVec (fftIdx2, SxMesh3D::Origin);

         /*
         // substract index vectors
         const SxVector3<Int> idxVecRes = idxVec1 - idxVec2;

         // for debugging only -- if this vector is out of the range of
         // the fft mesh, something got wrong
         const SxVector3<Int> &mesh = fft3d(0).mesh;
         int i;
         for (i = 0; i < 3; i++)
            if (2 * idxVecRes(i) <= -mesh(i) || 2 * idxVecRes(i) > mesh(i))
               return -1;

         // if vector is in fft range, get its index in this G (or G+k) basis
         const int fftIdxRes 
            = fft3d(0).mesh.getMeshIdx (idxVecRes, SxMesh3D::Origin);
         return n123inv(fftIdxRes);
         */

         // --- roughly 30% faster low-level implementation
         int m2 = fft3d(0).mesh(2);
         int a1 = fftIdx1 / m2;
         int a2 = fftIdx2 / m2, ix2;
         int dx, dy, dz;

         // z-coord
         {
            dz = fftIdx1 - a1 * m2;
            if (dz + dz > m2) dz -= m2;
            int iz2 = fftIdx2 - a2 * m2;
            dz -= iz2;
            if (iz2 + iz2 > m2) dz += m2;
            if (dz < 0)  {
               if (dz + dz <= -m2) return -1;
               dz += m2;
            } else {
               if (dz + dz > m2) return -1;
            }
         }

         // y-coord (plus 1st step of x-coord)
         int m1 = fft3d(0).mesh(1);
         dx = a1 / m1;
         dy = a1 - dx * m1;
         if (dy + dy > m1) dy -= m1;
         ix2 = a2 / m1;
         {
            int iy2 = a2 - ix2 * m1;
            if (iy2 + iy2 > m1)
               dy += m1 - iy2;
            else
               dy -= iy2;
         }
         if (dy < 0)  {
            if (dy + dy <= -m1) return -1;
            dy += m1;
         } else {
            if (dy + dy > m1) return -1;
         }

         // x-coord
         int m0 = fft3d(0).mesh(0);
         if (dx+dx > m0) dx -= m0;
         dx -= ix2;
         if (ix2 + ix2 > m0) dx += m0;
         if (dx < 0)  {
            if (dx + dx <= -m0) return -1;
            dx += m0;
         } else {
            if (dx + dx > m0) return -1;
         }

         // if vector is in fft range, get its index in this G (or G+k) basis
         ssize_t fftIdxRes 
            = fft3d(0).mesh.getMeshIdx (dx, dy, dz, SxMesh3D::Origin);
         return n123inv(fftIdxRes);
      }

      /** \brief Yields the index of the component belonging to the difference
                 of two indeces, the latter may belong to arbitrary
                 \f$ | \mathbf{G} + \mathbf{p},\mathbf{q} \rangle \f$
                 basis sets.

          Be \f$ i_{1} \f$ the index of \f$ G_{1} \f$ in the
          \f$ |\mathbf{G}+\mathbf{p}\rangle \f$ basis set and
          \f$ i_{2} \f$ the index of \f$ G_{2} \f$ in the
          \f$ |\mathbf{G}+\mathbf{q}\rangle \f$ basis set. The output of this
          function is then the index of the
          \f$ (\mathbf{G_{1}}-\mathbf{G_{2}}) \f$ component in this
          \f$ | \mathbf{G} + \mathbf{k} \rangle \f$ or
          \f$ | \mathbf{G} \rangle \f$ basis set.

          @param   idxGp index of a coefficient in the |G+p> basis
          @param  *gpPtr pointer to the |G+p> basis
          @param   idxGq index of a coefficient in the |G+q> basis
          @param  *gqPtr pointer to the |G+q> basis

          @return  the index of the coefficient of the difference
                   \f$\mathbf{G_{1}} - \mathbf{G_{2}}\f$ of the vectors
                   denoted by the two indeces, given in this |G+k> basis,
                   i.e., the basis of the object where you put the dot and
                   the name of this function behind.

          \note    Don't mix up *gpPtr and *gqPtr with pointers to objects of
                   the SxGkBasis class! *gpPtr and *gqPtr here denote objects
                   of the SxGBasis class.

          \note    Before you can use this function, the n123inv tables must
                   be initialised. Therefore use SxGBasis::n123invSetup ().

          \sa      SxGBasis::getIdxGSum (idxG1, idxG2)
                   SxGBasis::getIdxGSum (idxGp, *gpPtr, idxGq, *gqPtr),
                   SxGBasis::getIdxGDiff (idxG1, idxG2)

          \author  Matthias Wahn
       */
      inline int getIdxGDiff (int idxGp, const SxGBasis *gpPtr,
                              int idxGq, const SxGBasis *gqPtr) const
      {
         SX_CHECK     (gpPtr);
         SX_CHECK     (gqPtr);
         SX_CHECK (n123inv.getSize() > 0, n123inv.getSize());

         // get idx in 3-dim. mesh
         const PrecFFTIdx fftIdxGp = gpPtr->n123(0)(idxGp);
         const PrecFFTIdx fftIdxGq = gqPtr->n123(0)(idxGq);

         // get mesh vector with coordinates -mesh(i)/2 < vec(i) <= mesh(i)/2
         const RelVec idxVecGp
            = gpPtr->fft3d(0).mesh.getMeshVec (fftIdxGp, SxMesh3D::Origin);
         const RelVec idxVecGq
            = gqPtr->fft3d(0).mesh.getMeshVec (fftIdxGq, SxMesh3D::Origin);

         // substract index vectors
         const RelVec idxVecRes = idxVecGp - idxVecGq;

         // since the range of this vector may exceed the range of the fft
         // mesh belonging to this |G+k> basis, cut it pursuent to the mesh
         const RelVec mesh = fft3d(0).mesh;
         int i;
         for (i = 0; i < 3; i++)
            if (2 * idxVecRes(i) <= -mesh(i) || 2 * idxVecRes(i) > mesh(i))
               return -1;

         // if vector is in fft range, get its index in the |G+k> basis
         const ssize_t fftIdxRes 
            = fft3d(0).mesh.getMeshIdx (idxVecRes, SxMesh3D::Origin);
         return n123inv(fftIdxRes);
      }

      inline RelVec getGRot (SxMatrix3<Int> &symOpRel, int idxG) const
      {
         // get idx in 3-dim. mesh
         PrecFFTIdx fftIdxG = n123(0)(idxG);

         // get mesh vector with coordinates -mesh(i)/2 < vec(i) <= mesh(i)/2
         RelVec idxVecG = fft3d(0).mesh.getMeshVec (fftIdxG, SxMesh3D::Origin);

         // rotate index vector
         RelVec idxVecRes = symOpRel ^ idxVecG;

         return idxVecRes;
      }


      /** \brief  Yields the index of \f$ \mathbf{G} \f$ in another
                  \f$ | \mathbf{G} + \mathbf{q} \rangle \f$ basis.

          You know: the index \f$ i_{\mathbf k} \f$ of a wavefunction
          coefficient \f$ c_{n\mathbf{k}}(\mathbf{G}) \f$ belonging to the
          vector \f$ \mathbf{G} \f$ in this
          \f$|\mathbf{G}(+\mathbf{k})\rangle\f$ basis.

          You want to know: the index \f$ i_{\mathbf q} \f$ of the
          wavefunction coefficient \f$ c_{m\mathbf{q}}(\mathbf{G}) \f$
          belonging to the same \f$ \mathbf{G} \f$ but this time in the
          \f$|\mathbf{G}+\mathbf{q}\rangle\f$ basis.

          @param  idxG   index of \f$\mathbf{G}\f$ in this |G+k> basis
          @param  gqPtr  pointer to the |G+q> basis

          @return index of \f$\mathbf{G}\f$ in the |G+q> basis

          \note   Before you can use this function, the n123inv table of
                  the \f$ |\mathbf{G}(+\mathbf{q})\rangle \f$ basis must
                  be initialised. Therefore use SxGBasis::n123invSetup ().

          \author Matthias Wahn
       */
      inline int conveyIdx (int idxG, const SxGBasis *gqPtr) const
      {
         SX_CHECK     (gqPtr);
         SX_CHECK (gqPtr->n123inv.getSize() > 0, n123inv.getSize());

         // get the fft mesh index, and then the index in the |G+q> basis
         return gqPtr->n123inv( n123(0)(idxG) );
      }

      /** \brief the scalar product \f$ \langle a | b \rangle \f$ of two
                 vectors/functions

          computes the scalar product \f$ \langle a | b \rangle \f$ of two
          functions/vectors \f$ | a \rangle \f$ and \f$ | b \rangle \f$
          given in the \f$ (\mathbf{G} + \mathbf{q}) \f$ and in the
          \f$ (\mathbf{G} + \mathbf{k}) \f$ space representation, respectively,
          thereby \f$ \mathbf{k} \f$ referring to
          this \f$ |\mathbf{G}(+\mathbf{k})\rangle \f$ basis.

          @param  aBasis  the \f$ (\mathbf{G} + \mathbf{q}) \f$
                          basis, on which \f$ | a \rangle \f$ "lives"
          @param  a       1st argument of the scalar product
          @param  b       2nd argument of the scalar product

          @return \f$ \langle a | b \rangle \f$

          \note   Before you can use this function, the n123inv table of
                  the basis of \f$ |a\rangle \f$ must be initialised.
                  Therefore use SxGBasis::n123invSetup ().

          \author Matthias Wahn
       */
      SxComplex16 scalarProduct (const SxGBasis &aBasis,
                                 const SxDiracVec<Complex16> &a,
                                 const SxDiracVec<Complex16> &b);

      /** \brief shifts a function \f$ f_{\mathbf{k}} \f$ given on this
                 this \f$ |\mathbf{G}+\mathbf{k}\rangle \f$ basis set by a
                 reciprocal lattice vector \f$ \mathbf{G} \f$ given in
                 Cartesian coordinates.

                 @param f           a periodic function \f$ f_{\mathbf{k}} \f$
                                    in its \f$ (\mathbf{G} + \mathbf{k}) \f$
                                    space representation
                 @param deltaGCart  a reciprocal lattice vector in Cartesian
                                    coordinates
                 @param rMesh       the mesh points of the real space in
                                    Cartesian coordinates

                 @return            \f$ f_{\mathbf{k}+\mathbf{G}} \f$ in
                                    \f$ (\mathbf{G} + \mathbf{k}) \f$
                                    space representation

                 \sa                SxGBasis::shiftK (f, deltaGRel)

                 \author            Matthias Wahn
       */
      SxDiracVec<Complex16> shiftK (const SxDiracVec<Complex16> &f,
                                    const SxVector3<TPrecG>     &deltaGCart,
                                    const SxArray<Coord>        &rMesh);

      /** \brief shifts a function \f$ f_{\mathbf{k}} \f$ given on this
                 this \f$ |\mathbf{G}+\mathbf{k}\rangle \f$ basis set by a
                 reciprocal lattice vector \f$ \mathbf{G} \f$ given in
                 relative coordinates.

                 @param f           a periodic function \f$ f_{\mathbf{k}} \f$
                                    in its \f$ (\mathbf{G} + \mathbf{k}) \f$
                                    space representation
                 @param deltaGRel   a reciprocal lattice vector in relative
                                    coordinates, i.e., \f$ \mathbf{G}^{\rm rel}
                                    = \mathbf{G}\,B^{-1} \f$,
                                    where \f$ B \f$ denotes the Matrix of the
                                    reciprocal basis vectors (line-wise).

                 @return            \f$ f_{\mathbf{k}+\mathbf{G}} \f$ in
                                    \f$ (\mathbf{G} + \mathbf{k}) \f$
                                    space representation

                 \sa                SxGBasis::shiftK (f, deltaGCart, rMesh)

                 \author            Matthias Wahn
       */
      SxDiracVec<Complex16>
         shiftK (const SxDiracVec<Complex16> &f, RelVec &deltaGRel);

      /** \brief Transfer (lattice-periodic part of) waves to this G-basis
          @param in waves in some G-basis
          @return waves in this G basis

       */
      SxDiracVec<TPrecCoeffG>
      transferWaves (const SxDiracVec<TPrecCoeffG> &in) const;

      /** \brief Computes the cut-off radius of the \b G basis.

          It reads \f$4 e_{\rm cut}\f$. */
      static PrecEnergy getGCut (PrecEnergy eCut);
      static PrecEnergy getGCut (const SxSymbolTable *);
      /** \deprecated: not required anymore*/
      int   getMaxNG ();

      /** \todo saving of G basis is not yet implemented */
      void write (SxBinIO &);

      /** Read a single G-basis from netCDF file
        @param io     netcdf file
        @param ngIn   number of G vectors
        @param offset G-vector data offset in file
        @param mesh   FFT mesh (FFTs will not be initialized if 0)

        The G-vector data is stored contiguously for all G-bases contained
        in the netcdf file. By specifying the offset and the number of G
        vectors, this data can properly be read in without knowing the
        "official" k-index.
        */
      void read (const SxBinIO &io, int ngIn, int offset, 
                 const SxVector3<Int> &mesh = SxVector3<Int> (0,0,0));

      /** \brief rearrange vectors according to the FFT index order

          When wavefunction files are transferred between computer platforms
          they FFT index order might differ. Hence, when they are read in
          the coefficients have to be resorted according to the current
          FFT index order. This routine returns the resorted vector */
      PsiG mapToFFT (const PsiG &, const SxVector<TPrecFFTIdx> &,
                     const SxMesh3D &) const;


      /** \brief print debug information */
      void print () const;


   protected:

      /** \brief radius of the cutoff sphere */
      PrecEnergy        gCut;
      /** \brief optional translation vector */
      SxVector3<TPrecG> dG;

      /** \brief apply time reversal symmetry G=-G */
      bool             useTimeRevSym;  // G = -G

      /** \brief Define which variables are observed by the memory tracker */
      void registerMemoryObservers ();

   public:
      /** \brief Auxiliary function for memory consumption

          \note This is necessary because some memory consumers are
                protected and hence not accessible outside this class.
          \todo Think about a friend function as alternative.
        */
      size_t getNBytes () const;
};



// --------------------------------------------------------------------------

template<>
inline size_t getNBytes<SxGBasis> (const SxGBasis &G)
{
   return G.getNBytes () + sizeof(SxGBasis);
}

inline size_t SxGBasis::getNBytes () const
{
   return   ::getNBytes (gVec)
          + ::getNBytes (g2)
          + ::getNBytes (structureFactors)
          + ::getNBytes (phaseFactors)
          + ::getNBytes (packedGrel)
//        + ::getNBytes (fft3d)
          + ::getNBytes (n123)
          + ::getNBytes (n123inv);
}

#endif /* _SX_G_BASIS_H_ */
