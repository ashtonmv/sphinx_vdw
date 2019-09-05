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

#include <SxMatrix3.h>

#ifndef SXMAT
#  error "Undefined variable SXMAT"
#endif
#ifndef SXMATBASIS
#  error "Undefined variable SXMATBASIS"
#endif



// According to lookup rules in derived template classes base elements have 
// to make visible explictly by using 'this->...'.

/**
  @ingroup group_num
  */
template<class T>
class SXMAT : public SXMATBASIS<T>
{
   public:
      /**@name Construction and Destruction */
      //@{
      inline SXMAT();
      inline SXMAT (const SXMATBASIS<Int> &s);
      inline SXMAT (const SXMATBASIS<Float> &s);
      inline SXMAT (const SXMATBASIS<Double> &s);
      inline SXMAT (const SXMATBASIS<Complex8> &s);
      inline SXMAT (const SXMATBASIS<Complex16> &s);
      inline SXMAT (const SxMatrix3<T> &m);
      inline SXMAT (const SxList<typename T::Type> &l);
      inline SXMAT (ssize_t r, ssize_t c, SxList<typename T::Type> &l);
      inline SXMAT (ssize_t nr, ssize_t nc);
      inline SXMAT (void *in, ssize_t nElem);
      inline explicit SXMAT (ssize_t nElem);
      inline virtual ~SXMAT ();
      //@}

      bool isHermitian () const;


      /** Controls the behaviour of the eigensolver routine. */
      enum EigCmd {
         /// Eigensolver will compute both eigenvalues and eigenvectors
         ALL=-1, 
         /// Eigensolver will compute eigenvalues.
         EIGENVALS_ONLY=-2, 
         /// Eigensolver will compute eigenvectors
         EIGENVECS_ONLY=-3, 
         /// Eigensolver will compute optimal work size only.
         EIGENSIZE=-4
      };

      typedef SXMATBASIS<typename T::TCmplx>  Eigenvalues;
      typedef SXMATBASIS<T>                   Eigenvectors;
      class Eigensystem
      {
         public:
            Eigenvalues  vals;
            Eigenvectors vecs;

            Eigensystem () { }
            Eigensystem (const Eigensystem &in)  {
               this->vals = in.vals;
               this->vecs = in.vecs;
            }
            Eigensystem &operator= (const Eigensystem &in)  {
               if ( &in == this )  return *this;
               vals    = in.vals;
               vecs    = in.vecs;
               return *this;
            }
            void print (bool vectorForm=false) const  {
               SX_CHECK (vals.getSize() > 0 && vecs.getSize() > 0,
                         vals.getSize(),       vecs.getSize());
               //   SX_CHECK (vals.getSize() == vecs.getSize(),
               //               vals.getSize(),   vecs.getSize());

               printf ("Eigenvalues:\n");
               vals.print (vectorForm);
               printf ("Eigenvectors:\n");
               vecs.print (vectorForm);
            }
            
      };

      typename SXMAT<T>::Eigenvalues  eigenvalues  (ssize_t size=0, 
                                                    bool sorting=true);
      typename SXMAT<T>::Eigenvectors eigenvectors (ssize_t size=0, 
                                                    bool sorting=true);
      typename SXMAT<T>::Eigensystem  eigensystem  (ssize_t size=0, 
                                                    bool sorting=true) const;
      typename SXMAT<T>::Eigensystem  eigensystem  (bool sorting) const;
      typename SXMAT<T>::Eigensystem  eigensystem  (EigCmd, ssize_t, 
                                                    bool sorting=true) const;


      /* Returns the matrix
          \f[
             U^\frac{1}{2} = 
                 \left( \begin{matrix}
                   v_1 \\ v_2 \\ \vdots \\ v_n
                 \end{matrix} \right)
                 \left( \begin{matrix}
                   \varepsilon_1^{-\frac{1}{2}} & 0 & \cdots & 0  \\
                   0 & \varepsilon_2^{-\frac{1}{2}} & \cdots & 0  \\
                   \vdots & \vdots                        & \ddots & \vdots \\
                   0 & 0   &  \cdots & \varepsilon_n^{-\frac{1}{2}} 
                 \end{matrix} \right)
                 \left( \begin{matrix}
                   v_1 & v_2 & \cdots & v_n
                 \end{matrix} \right)
          \f] with v being the eigenvectors of current matrix M
          \f[
             M v = \varepsilon v.
          \f]

          \par Usage:
          Here is an example for orthogonalization using U matrix:
          \verbatim
             ...
             SXMAT<...> U, S;
             S = psi.adjoint() ^ psi;          // overlap matrix
             U = S.getU ();
             X = ...                           // e.g. steepest descent vector
             X = (U ^ X.adjoint()).adjoint();  // rotate vector X
             ...
          \endverbatim
          
          \brief Returns matrix \f$ U^{\frac{1}{2}}\f$.
          \see   SxPsiSet::orthogonalize
          \param psi vector \f$ \Psi \f$.
          \return \f$ U^{\frac{1}{2}}\f$
          */
      inline SXMATBASIS<T> getU ();

      
};



//------------------------------------------------------------------------------
// Constructors / Destructor
//------------------------------------------------------------------------------
template<class T> SXMAT<T>::SXMAT () : SXMATBASIS<T> () { }
template<class T> SXMAT<T>::SXMAT (const SXMATBASIS<Int> &s) 
   : SXMATBASIS<T> (s) { }
template<class T> SXMAT<T>::SXMAT (const SXMATBASIS<Float> &s) 
   : SXMATBASIS<T> (s) { }
template<class T> SXMAT<T>::SXMAT (const SXMATBASIS<Double> &s) 
   : SXMATBASIS<T> (s) { }
template<class T> SXMAT<T>::SXMAT (const SXMATBASIS<Complex8> &s) 
   : SXMATBASIS<T> (s) { }
template<class T> SXMAT<T>::SXMAT (const SXMATBASIS<Complex16> &s) 
   : SXMATBASIS<T> (s) { }
template<class T> SXMAT<T>::SXMAT (const SxMatrix3<T> &m) 
   : SXMATBASIS<T> (9) 
{ 
   this->elements[0]=m(0,0); this->elements[1]=m(1,0); this->elements[2]=m(2,0);
   this->elements[3]=m(0,1); this->elements[4]=m(1,1); this->elements[5]=m(2,1);
   this->elements[6]=m(0,2); this->elements[7]=m(1,2); this->elements[8]=m(2,2);
   this->reshape (3,3);
}
template<class T> SXMAT<T>::SXMAT (const SxList<typename T::Type> &l) 
   : SXMATBASIS<T> (l) { }
template<class T> SXMAT<T>::SXMAT (ssize_t r, ssize_t c, 
                                         SxList<typename T::Type> &l) 
   : SXMATBASIS<T> (l) 
{ 
      SXMAT<T> tmp; tmp.copy (*this); tmp.reshape(c, r);
      *this = tmp.transpose(); 
}
template<class T> SXMAT<T>::SXMAT (void *in, ssize_t nElem)
   : SXMATBASIS<T> (in, nElem) { }
template<class T> SXMAT<T>::SXMAT (ssize_t nr, ssize_t nc)
   : SXMATBASIS<T> (nr * nc) { this->reshape (nr, nc); }
template<class T> SXMAT<T>::SXMAT (ssize_t nElem)
   : SXMATBASIS<T> (nElem) { }
template<class T> SXMAT<T>::~SXMAT () { }



template<class T> 
bool SXMAT<T>::isHermitian () const
{
   ssize_t i, j, n = this->nRows();
   if (n != this->nCols())  return false;

   for (i=0; i < n; i++)  {
      // diag must be real
      if ( fabs(typename T::Cmplx ((*this)(i,i)).im) > 1e-12 )
         return false;  
      // off-diag elements:  M_ij = M~_ij
      for (j=0; j < i; j++)  {
         const typename T::Cmplx &m_ij = (*this)(i,j);
         const typename T::Cmplx &m_ji = (*this)(j,i);
         if (   fabs(m_ij.re - m_ji.re) > 1e-12
             || fabs(m_ij.im + m_ji.im) > 1e-12)
            return false;
      }
   }
   return true;

}



//------------------------------------------------------------------------------
// Eigensolver
//------------------------------------------------------------------------------
template<class T>
typename SXMAT<T>::Eigenvalues SXMAT<T>::eigenvalues (ssize_t size,bool sorting)
{
   return eigensystem (EIGENVALS_ONLY, size, sorting).vals;
}


template<class T>
typename SXMAT<T>::Eigenvectors SXMAT<T>::eigenvectors (ssize_t size,
                                                        bool    sorting)
{
   return eigensystem (EIGENVECS_ONLY, size, sorting).vecs;
}


template<class T>
typename SXMAT<T>::Eigensystem 
SXMAT<T>::eigensystem (bool sorting) const
{
   return eigensystem (ALL, 0, sorting);
}


template<class T>
typename SXMAT<T>::Eigensystem 
SXMAT<T>::eigensystem (ssize_t size, bool sorting) const
{
   return eigensystem (ALL, size, sorting);
}


template<class T>
typename SXMAT<T>::Eigensystem 
SXMAT<T>::eigensystem (EigCmd, ssize_t, bool sorting) const
{
   SX_CHECK (this->nRows() == this->nCols(), this->nRows(), this->nCols());

   Eigensystem eig;
   eig.vals = Eigenvalues  (this->nRows());
   eig.vecs = Eigenvectors (this->nRows() * this->nCols());
   eig.vecs.reshape (this->nRows(), this->nCols());

   // create work vector (original vector will be destroyed by LAPACK!)
   SXMATBASIS<T> inpVec; inpVec.copy (*this);

   matEigensolver (eig.vals.elements, eig.vecs.elements, inpVec.elements,
                   static_cast<int>(this->nRows()), All);

   // sort by ascending (real parts of) eigenvalues
   if (sorting)  {
      ssize_t n = this->nRows();
      SxArray<ssize_t> sortIdx (n); 
      sortIdx = (SXMATBASIS<typename T::TReal> (eig.vals)).getSortIdx();
      Eigenvectors tmpVecs; tmpVecs.copy (eig.vecs);
      for (ssize_t i=0; i < n; i++)  
         eig.vecs.colRef(i) <<= tmpVecs.colRef(sortIdx(i));
      eig.vals.sortByIdx (sortIdx);	

   }

   // --- Gram-Schmidt orthogonalization of degenerate eigenvectors
   //     for complex Hermitean matrices (real matrices have strange
   //     convention for eigenvectors)
   if (sizeof(typename T::Cmplx) == sizeof(typename T::Type))  {
      bool checkedHermitean = false;
      for (ssize_t i = 0; i < this->nRows(); ++i)  {
         if (fabs(typename T::Cmplx(eig.vals(i)).im) > 1e-12) break; // not Hermitean

         bool renormalize = false;
         for (ssize_t j = 0; j < i; ++j)  {
            typename T::Cmplx d = eig.vals(i) - eig.vals(j);
            if (d.absSqr () < 1e-16)  {
               if (!checkedHermitean)  {
                  for (ssize_t k = i; k < this->nRows (); ++k)
                     if (fabs(typename T::Cmplx(eig.vals(i)).im) > 1e-12)
                        return eig; // not Hermitean
                  if (! this->isHermitian ())
                     return eig; // not Hermitean
                  checkedHermitean = true;
               }
               typename T::Type scp = dot(eig.vecs.colRef(j), 
                                          eig.vecs.colRef(i));
               eig.vecs.colRef(i).plus_assign_ax (-scp, eig.vecs.colRef(j));
               renormalize = true;
            }
         }
         if (renormalize) eig.vecs.colRef(i).normalize ();
      }
   }

   return eig;
}



template<class T>
SXMATBASIS<T> SXMAT<T>::getU ()
{
   SXMAT<TPrecCoeffG> S, Ieps;
   SXMATBASIS<TPrecCoeffG> U;

   //printf (">>> %d,%d\n", nRows(), nCols());
   Eigensystem eig;
   SXMAT<T> V = *this;
   S     = V.adjoint() ^ V;                      // S being the overlap
   eig   = S.eigensystem(true);                  // eigenproblem: Sv = ev
   Ieps  = Ieps.identity(1. / sqrt(eig.vals));   // diagonal: 1/sqrt(eps)
   U     = eig.vecs ^ Ieps ^ eig.vecs.adjoint(); // U^

   return U;
}

#ifdef USE_LOOPMPI
template<>
inline void SxLoopMPI::sum (SXMAT<Int> &inout)
{
   // actual implementation is in SxLoopMPI
   SxLoopMPI::sum (inout.elements, inout.elements, inout.getSize ());
}

template<>
inline void SxLoopMPI::sum (SXMAT<Double> &inout)
{
   // actual implementation is in SxLoopMPI
   SxLoopMPI::sum (inout.elements, inout.elements, inout.getSize ());
}

template<>
inline void SxLoopMPI::sum (SXMAT<Complex16> &inout)
{
   // actual implementation is in SxLoopMPI
   SxLoopMPI::sum ((double*)inout.elements, (double*)inout.elements,
                   inout.getSize () * 2);
}

template <>
inline void SxLoopMPI::bcast (SXMAT<Double> &inout, int source)
{
   SxLoopMPI::bcast(inout.elements, inout.getSize(), source);
}

template <>
inline void SxLoopMPI::bcast (SXMAT<Complex16> &inout, int source)
{
   SxLoopMPI::bcast((double*)inout.elements, 2 * inout.getSize(), source);
}
#else
template<> inline void SxLoopMPI::sum (SXMAT<Int> &) { }
template<> inline void SxLoopMPI::sum (SXMAT<Double> &) { }
template<> inline void SxLoopMPI::sum (SXMAT<Complex16> &) { }
template<> inline void SxLoopMPI::bcast (SXMAT<Double> &, int) { }
template<> inline void SxLoopMPI::bcast (SXMAT<Complex16> &, int) { }
#endif


