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

#ifndef SXSYMMAT
#  error "Undefined variable SXSYMMAT"
#endif
#ifndef SXMATBASIS
#  error "Undefined variable SXMATBASIS"
#endif


// According to lookup rules in derived template classes base elements have 
// to be made visible explictly by using 'this->...'.

/**
  @ingroup group_num
  */
template<class T>
class SXSYMMAT : public SXMATBASIS<T>
{
   public:
      /**@name Construction and Destruction */
      //@{
      inline SXSYMMAT();
      inline SXSYMMAT (const SXMATBASIS<Float> &s);
      inline SXSYMMAT (const SXMATBASIS<Double> &s);
      inline SXSYMMAT (const SXMATBASIS<Complex8> &s);
      inline SXSYMMAT (const SXMATBASIS<Complex16> &s);
      inline SXSYMMAT (SxList<typename T::Type> &l);
      inline SXSYMMAT (void *in, ssize_t nElem);
      inline explicit SXSYMMAT (ssize_t nElem);
      inline virtual ~SXSYMMAT ();

      //@}

      /** Controls the behaviour of the eigensolver routine. */
      enum EigCmd {
         /// Eigensolver will compute both eigenvalues and eigenvectors
         ALL            = -1, 
         /// Eigensolver will compute eigenvalues.
         EIGENVALS_ONLY = -2, 
         /// Eigensolver will compute eigenvectors
         EIGENVECS_ONLY = -3, 
         /// Eigensolver will compute optimal work size only.
         EIGENSIZE      = -4
      };

      typedef SXMATBASIS<typename T::TReal>  Eigenvalues;
      typedef SXMATBASIS<T>                  Eigenvectors;
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

      typename SXSYMMAT<T>::Eigenvalues  eigenvalues  ();
      typename SXSYMMAT<T>::Eigenvectors eigenvectors ();
      /** \brief Get eigensystem.
          \param sorting If set to true (default), the eigenvalues
                         and eigenvectors are sorted, with increasing
                         eigenvalues. The effort should be negligible
                         compared to finding the eigenvectors/values.
       */
      typename SXSYMMAT<T>::Eigensystem  
      eigensystem  (bool sorting = true) const;
      /** \brief Get eigensystem using specified algorithm
          \param EigCmd  specifies eigensolver
          \param sorting If set to true (default), the eigenvalues
                         and eigenvectors are sorted, with increasing
                         eigenvalues. The effort should be negligible
                         compared to finding the eigenvectors/values.
       */
      typename SXSYMMAT<T>::Eigensystem  
      eigensystem  (EigCmd, bool sorting = true) const;


};



//------------------------------------------------------------------------------
// Constructors / Destructor
//------------------------------------------------------------------------------
template<class T> SXSYMMAT<T>::SXSYMMAT () 
   : SXMATBASIS<T> () { this->handle->parameters 
                          |= (IS_TRIANGULAR | UPPER_RIGHT); }
template<class T> SXSYMMAT<T>::SXSYMMAT (const SXMATBASIS<Float> &s) 
   : SXMATBASIS<T> (s) { }
template<class T> SXSYMMAT<T>::SXSYMMAT (const SXMATBASIS<Double> &s) 
   : SXMATBASIS<T> (s) { }
template<class T> SXSYMMAT<T>::SXSYMMAT (const SXMATBASIS<Complex8> &s) 
   : SXMATBASIS<T> (s) { }
template<class T> SXSYMMAT<T>::SXSYMMAT (const SXMATBASIS<Complex16> &s) 
   : SXMATBASIS<T> (s) { }
template<class T> SXSYMMAT<T>::SXSYMMAT (SxList<typename T::Type> &l) 
   : SXMATBASIS<T> (l) 
{
   this->handle->parameters |= (IS_TRIANGULAR | UPPER_RIGHT); 
   float x = (float)(l.getSize());
   ssize_t n = (ssize_t)(0.5 * sqrt(1. + 8.*x) - 0.5);  // from: x=n(n+1)/2
   this->reshape (n, n); 
}
template<class T> SXSYMMAT<T>::SXSYMMAT (void *in, ssize_t nElem)
   : SXMATBASIS<T> (in, nElem) 
{ 
   SX_EXIT; // not implemented yet
}
template<class T> SXSYMMAT<T>::SXSYMMAT (ssize_t nElem)
   : SXMATBASIS<T> (nElem+(nElem-1)*nElem/2) 
{ 
   this->handle->parameters |= (IS_TRIANGULAR | UPPER_RIGHT); 
   this->reshape (nElem, nElem); 
}
template<class T> SXSYMMAT<T>::~SXSYMMAT () { }




//------------------------------------------------------------------------------
// Eigensolver
//------------------------------------------------------------------------------
template<class T>
typename SXSYMMAT<T>::Eigenvalues 
SXSYMMAT<T>::eigenvalues ()
{
   return eigensystem (EIGENVALS_ONLY).vals;
}


template<class T>
typename SXSYMMAT<T>::Eigenvectors 
SXSYMMAT<T>::eigenvectors ()
{
   return eigensystem (EIGENVECS_ONLY).vecs;
}


template<class T>
typename SXSYMMAT<T>::Eigensystem 
SXSYMMAT<T>::eigensystem (bool sorting) const
{
   return eigensystem (ALL, sorting);
}


template<class T>
typename SXSYMMAT<T>::Eigensystem 
SXSYMMAT<T>::eigensystem (EigCmd cmd, bool sorting) const
{
   SX_CHECK (this->nRows() == this->nCols(), this->nRows(), this->nCols());
   SX_CHECK (cmd == ALL);  // current interface supports only 'All'

   Eigensystem eig;
   eig.vals = Eigenvalues  (this->nRows());
   eig.vecs = Eigenvectors (this->nRows() * this->nCols());
   eig.vecs.reshape (this->nRows(), this->nCols());

   // create work vector (original vector will be destroyed by LAPACK!)
   SXMATBASIS<T> inpVec; inpVec.copy (*this);

   SX_CHECK ( this->handle->parameters & UPPER_RIGHT );

   // --- change UpperRight -> LowerLeft and conjugate eigenvectors
   matEigensolverTri (eig.vals.elements, eig.vecs.elements, inpVec.elements,
                      static_cast<int>(this->nRows()), LowerLeft, All);
   if (sorting)  {
      ssize_t n = this->nRows ();
      SxArray<ssize_t> sortIdx;
      sortIdx = (SXMATBASIS<typename T::TReal> (eig.vals)).getSortIdx();
      // sort eigenvectors
      Eigenvectors tmpVecs = eig.vecs;
      eig.vecs = Eigenvectors (n * n); eig.vecs.reshape (n, n);
      for (ssize_t i = 0; i < n; ++i)
         eig.vecs.colRef(i) <<= tmpVecs.colRef(sortIdx(i)).conj ();
      // sort eigenvalues
      eig.vals.sortByIdx (sortIdx);

   } else {
      eig.vecs = eig.vecs.conj();
   }

   return eig;
}

template<class T>
std::ostream& operator<< (std::ostream &s, const SXSYMMAT<T> &in)
{
   SX_CHECK     (in.handle);
   SX_CHECK (in.handle->size >= 0, in.handle->size);
   SX_CHECK     (in.handle->parameters & IS_TRIANGULAR);

   ssize_t r, c;
   int width = s.width();

   if (in.handle->size == 0)  {            // empty

      s << "empty";

   }  else  {
      // --- matrix print out
      ssize_t n = in.nRows ();
      if (in.handle->parameters & UPPER_RIGHT)  {
         for (r=0; r < n; r++)  {
            s << '{';
            for (c=r; c < n - 1; c++)  {
               s << setw(width) << in(r,c) << ", ";
            }
            s << setw(width) << in(r,n-1) << '}' << endl;
         }
      }  else  { // lower left
         for (r=0; r < n; r++)  {
            s << '{';
            for (c=0; c < r; c++)  {
               s << setw(width) << in(r,c) << ", ";
            }
            s << setw(width) << in(r,r) << '}' << endl;
         }
      }
   }
   
   return s;
}


