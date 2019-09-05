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

#ifndef _SX_GK_BASIS_H_
#define _SX_GK_BASIS_H_

#include <SxPrecision.h>
#include <SxArray.h>
#include <SxList.h>
#include <SxMatrix.h>
#include <SxDirac.h>
#include <SxVector.h>
#include <SxTypes.h>
#include <SxGBasis.h>
#include <SxDFT.h>
#include <SxSymbolTable.h>
#include <SxKPoints.h>

class SX_EXPORT_DFT SxGkBasis : public SxKPoints
{
   public:
     
      ///:ik,:iSpin,{ig x i}
      typedef SxGBasis::TPsi    TPsi;
      typedef SxGBasis::TPsiGI  TPsiGI;
      typedef SxGBasis::TRho    TRho;

      SxArray<SxPtr<SxGBasis> >              gBasisList;
      const SxGBasis                        *gBasis;
      PrecEnergy                             gkCut;

      /** Determines the index of a given \f$ \bf{G+k} \f$ 
          vector in the gBasis after subtracting
          the \f$ \bf{k} \f$ vector. 
          \sa lookup
      */
      SxMatrix<Int> GkInG;        // G+k --> G

      /** Determines the index of the vector \f$ \bf{G+G'} \f$ 
         in a subset of the gBasis determined by SxGkBasis::nGMax
          \sa lookup
      */
      SxMatrix<Int> GplusGinG;    // G+G'--> G
      
     
      /** Determines the index of a given \f$ \bf{G} \f$ 
          vector in the gkBasis after adding the proper
          \f$ \bf{k} \f$ vector. 
          \sa lookup
      */
      SxMatrix<Int> GinGk; 	    // G   --> G+k

      /** Determines the index of the vector \f$ \bf{G-G'} \f$ 
          in the gBasis.
          \sa lookup
      */
      SxMatrix<Int> GminusGinG;   // G-G'--> G

      /** This is the size of the largest SxGkBasis. 
          \sa SxGkBasis::nGkMin */ 
      int nGkMax; 

      /** This is the size of the smallest SxGkBasis. 
          \sa SxGkBasis::nGkMax */ 
      int nGkMin;

      /** This is the size of the subset of the SxGBasis that contains
          all \f$ \bf{G} \f$ vectors used to generate all
          SxGkBasis  */
      int nGMax;

      /** setup of the n123inv tables, currently mainly used for the exact
          exchange formalism */
      void setupN123inv ();

      /** setup nGkMin and nGkMax

          Usually those two parameters are setup automatically by calling
          the constructor of the SxGkBasis class. But when initializing it
          by reading from file, these parameters are omitted.

          \TODO write out and read in nGkMin/Max */
      void setupNGkMinMax ();

      /** \brief  Yields the index of G in another |G+q> basis.

          You know: the index \f$ i_{\mathbf k} \f$ of a wavefunction
          coefficient \f$ c_{n\mathbf{k}}(\mathbf{G}) \f$ belonging to the
          vector \f$ \mathbf{G} \f$ in the \f$|\mathbf{G}+\mathbf{k}\rangle\f$
          basis.

          You want to know: the index \f$ i_{\mathbf q} \f$ of the
          wavefunction coefficient \f$ c_{m\mathbf{q}}(\mathbf{G}) \f$
          belonging to the same \f$ \mathbf{G} \f$ but this time in the
          \f$|\mathbf{G}+\mathbf{q}\rangle\f$ basis.

          @param  index of \f$\mathbf{G}\f$ in the |G+k> basis
          @param  pointer to the |G+k> basis
          @param  pointer to the |G+q> basis

          @return index of \f$\mathbf{G}\f$ in the |G+q> basis

          \note   Don't mix *gkPtr up with a pointer to an object of
                  this SxGkBasis class! *gkPtr here denotes an object of
                  SxGBasis.
          */
      inline int transmitIdx (int gkBasIdxG, const SxGBasis *gkPtr,
                                             const SxGBasis *gqPtr) const
      {
         // get fft mesh index, and then the index in |G+q> basis
         return gqPtr->n123inv( gkPtr->n123(0)(gkBasIdxG) );
      }

      /** \brief  Yields the index of the component of sum of two indeces,
                  the latter belonging to two different |G+k> basis sets.

          @param  index of a coefficient in the |G+k> basis
          @param  pointer to the |G+k> basis
          @param  index of a coefficient in the |G+q> basis
          @param  pointer to the |G+q> basis

          @return the index of the coefficient of the sum of the
                  \f$\mathbf{G}\f$ vectors denoted by the to indeces,
                  given in the |G+q> basis.

          \note   The result is given in the |G+q> basis, i.e., where
                  the 2nd basis pointer points to.

          \note   Don't mix *gkPtr up with a pointer to an object of
                  this SxGkBasis class! *gkPtr here denotes an object of
                  SxGBasis.

          \sa     SxGBasis::getIdxGSum (idxG1, idxG2)
          */
      inline int getIdxGSum (int idxGk, const SxGBasis *gkPtr,
                             int idxGq, const SxGBasis *gqPtr) const
      {
         SX_CHECK (gkPtr);
         SX_CHECK (gqPtr);

         // get idx in 3-dim. mesh
         const PrecFFTIdx fftIdxGk = gkPtr->n123(0)(idxGk);
         const PrecFFTIdx fftIdxGq = gqPtr->n123(0)(idxGq);

         // get mesh vector with coordinates -mesh(i)/2 < vec(i) <= mesh(i)/2
         const RelVec idxVecGk
            = gkPtr->fft3d(0).mesh.getMeshVec (fftIdxGk, SxMesh3D::Origin);
         const RelVec idxVecGq
            = gqPtr->fft3d(0).mesh.getMeshVec (fftIdxGq, SxMesh3D::Origin);

         // add index vectors
         const RelVec idxVecRes = idxVecGk + idxVecGq;

         // since the range of this vector exceeds the range of the fft mesh,
         // cut it pursuent to the fft mesh
         const RelVec mesh = gqPtr->fft3d(0).mesh;
         int i;
         for (i = 0; i < 3; i++)
            if (2 * idxVecRes(i) <= -mesh(i) || 2 * idxVecRes(i) > mesh(i))
               return -1;

         // if vector is in fft range, get its index in the |G+q> basis
         const ssize_t fftIdxRes 
            = gqPtr->fft3d(0).mesh.getMeshIdx (idxVecRes, SxMesh3D::Origin);
         return gqPtr->n123inv(fftIdxRes);
      }

      /// Dummy function
      virtual ssize_t getNElements () const 
      {
         SX_EXIT;
         return -1;
      }

      void setNComp(int);
      
      SxDiracVec<Int> GIndx,  GIndp,  nStars; 
      SxMatrix<Int> MGIndx, MnStars;     //:ig,:jg
      SxDiracVec<Int> MGIndp;
      int nGIndp; 

      void symmetrize_g (SxDiracVec<TPrecCoeffG> &vector); 
      void symmetrize_g (SxMatrix<TPrecCoeffG> &matrix); 
      
      SxGkBasis ();
      SxGkBasis (const  SxBinIO &in, bool initFFT = true, bool saveMemory = false) {read (in, initFFT, saveMemory);};
      SxGkBasis (const  SxGBasis &gBasis,
                 const  SxSymbolTable *table);
      /** \brief Setup Gk basis from k-points, G-Basis and cutoff energy.
          \param kp         kPoints object
          \param G_         G-basis (what do we need that for?)
          \param cutoff     cutoff energy
          \param saveMemory if true, do not store phase factors
        */
      SxGkBasis (const SxKPoints &kp,
                 const SxGBasis &G_,
                 double cutoff,
                 bool saveMemory = false);
      SxGkBasis (const SxKPoints &kp,
                 const SxAtomicStructure &structure,
                 const SxMesh3D &mesh,
                 double cutoff,
                 bool saveMemory = false);

      SxGkBasis (const SxGkBasis &x) : SxKPoints (x) { SX_EXIT; }
      void operator= (const SxGkBasis &) { SX_EXIT; }

      virtual ~SxGkBasis ();
      void init (bool);

      SxGBasis &operator() (int);
      const SxGBasis &operator() (int) const;

      SxGBasis &operator() (SxAutoLoop &ik)
      {
         ik.setLimit (getNk ());
         return operator() (int(ik.i));
      }

      const SxGBasis &operator() (SxAutoLoop &ik) const
      {
         ik.setLimit (getNk ());
         return operator() (int(ik.i));
      }

      void changeTau (const SxAtomicStructure &tauList);

      /** \brief Return the structure from the first available G+k-Basis */
      const SxAtomicStructure &getTau () const
      {
         if (gBasis && gBasis->structPtr) return *gBasis->structPtr;
         for (int ik = 0; ik < gBasisList.getSize (); ++ik)
            if (gBasisList(ik) && gBasisList(ik)->structPtr)
               return *gBasisList(ik)->structPtr;
         cout << "Missing structure in GkBasis!" << endl;
         SX_EXIT;
      }

      /** Returns \em true if both \f$ \bf{G} \f$ vectors are equal 
          up to a threshold SxGkBasis::DLT  */
      bool equal(const SxVector3<TPrecG> &, const SxVector3<TPrecG> &) const;
      
      /** Returns \em true if both real numbers are equal 
          up to a threshold SxGkBasis::DLT  */
      bool equal(const PrecG &, const PrecG &) const;

      /**
         This function implements a binary search tree
         to search for the position of the vector \f$ G \f$ 
         in the basis pointed to by gB. 

         In exact exchange implementation we encountered the need to 
         lookup the vector \f$ \bf{(G+k + G'+k)} \f$ in the set of 
         \f$ \bf{(G+q)} \f$ 's. The direct application of this 
         lookup function, although possible, would dramatically slow
         down the code. If the indeces where to be tabulated this  
         would involve a huge matrix of the size (ng x ng x nk x nk)
         to span the indeces over \f$ \bf{G, G', k,}\f$ 
         and \f$ \bf{q} \f$.
         To optimize this lookup, we created three tables:
         - SxGkBasis::GkInG, 
         - SxGkBasis::GplusGinG, and 
         - SxGkBasis::GinGk.

         The first maps \f$ \bf{G+k} \mapsto \bf{G} \f$, the
         second maps \f$ \bf{G+G'} \mapsto \bf{G} \f$ and the third
         maps \f$ \bf{G} \f$ back to \f$ \bf{G+k} \f$.
         
         The procedure to use these tables is quit simple:
         lookup (G+k + G'+k) in G+q is found via

         \code
           int iG, jGp;    //  indeces of G+k and G'+k
           int ik, iq;  
           int lookup;
           lookup = GinGk(GplusGinG( GkinG(iG,ik), GkinG(iGp,ik) ), iq);
         \endcode

         The most general case of course would be
         lookup (G+k + G'+k') in G+q:

         \code
           int iG, jGp;    //  indeces of G+k and G'+k
           int ik, iq;
           int ikp;        //  use k' instead of k
           int lookup;
           lookup = GinGk(GplusGinG( GkinG(iG,ik), GkinG(iGp,ikp) ), iq);
         \endcode

         In the construction of the Hamiltonian and semi-local matrices
         we encountered the need to lookup (G+k - G'+k) in the SxGBasis.
         A fourth table was constructed:  SxGkBasis::GminusGinG.
         This table locates the vector \f$ G'' = (G+k - G'+k) \f$ 
         in the gBasis. To eliminate the k-dependence
         SxGkBasis::GkInG should be used, therefore this lookup is 
         performed via:
         \code
         int iG, iGp, ik;
         int lookup = GminusGinG( GkinG(iG, ik), GkinG(iGp,ik) )
         \endcode

         \param   G  an SxVector3<TPrecG> generic G vector.
         \param   gB a pointer to an SxGbasis.
         \returns index of \b G. 
         \author  Abdallah Qteish
         \author  Abdullah Al-Sharif
      */   
      int lookup (const SxVector3<TPrecG> &G, const SxGBasis *gB) const;

      /** 
        This is only a shortcut of the sentence
      \code  
        GinGk(GplusGinG( GkInG(igk,ik), GkInG(igq,iq) ), iq );
      \endcode  
      */
      inline int lookup (int ik, int igk, int iq, int igq); 

      const SxList<int> ngPerK () const;
      const SxList<int> ngPerK (PrecEnergy eCut) const;

      void write (const SxBinIO &) const;

      /** \brief Initialize from a file 
        @param initFFT whether to initialize the G->R FFT's
       */
      void read (const SxBinIO &, bool initFFT = true, bool saveMemory = false);

      /** Printout */
      void print () const;

      // get smallest and largest |G+k|^2
      double getMinGk () const;
      double getMaxGk () const;

      SxArray<ssize_t> getKSortIdx () const;

      inline int getTotalDataPoints () const
      {
         int result = 0;
         for(int ik = 0; ik < nk; ik++) result += (int)(*this)(ik).g2.getSize();
         return result;
      }


   protected:
      void setupBasisList (const SxAtomicStructure &structure, const SxMesh3D &mesh, bool);

      //void setupLookupTables ();

      void registerMemoryObservers ();
      
      
};

#endif /* _SX_GK_BASIS_H_ */
