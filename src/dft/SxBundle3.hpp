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

template<class T>
SxBundle3<T>::SxBundle3 (int nStates, int nSpin, int nk)
{
   bundle.resize (nk);
   int ik, iSpin;
   for (ik=0; ik < nk; ik++)  {
      bundle(ik).resize (nSpin);
      for (iSpin=0; iSpin < nSpin; iSpin++)  {
         bundle(ik)(iSpin).resize(nStates);
         bundle(ik)(iSpin).setBasis(NULL);
      }
   }
}


template<class T>
SxBundle3<T>::SxBundle3 (const SxBundle3<T> &in)
{
   int nStates, iSpin, nSpin, ik, nk;

   nk = int(in.bundle.getSize());
   bundle.resize (nk);

   for (ik=0; ik < nk; ik++)  {
      nSpin = int(in.bundle(ik).getSize());
      bundle(ik).resize (nSpin);

      for (iSpin=0; iSpin < nSpin; iSpin++)  {
         nStates = (int)in.bundle(ik)(iSpin).getSize();
         bundle(ik)(iSpin).resize (nStates);
         bundle(ik)(iSpin).copy ( in.bundle(ik)(iSpin) );
      }
   }
}



//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
template<class T>
SxBundle3<T>::~SxBundle3 ()
{

}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
template<class T>
SxBundle3<T> &SxBundle3<T>::operator= (const SxBundle3<T> &in)
{
   if ( &in == this )  return *this;

   int iSpin, nSpin, ik, nk;

   nk = int(in.bundle.getSize());
   bundle.resize (nk);

   for (ik=0; ik < nk; ik++)  {
      nSpin = int(in.bundle(ik).getSize());
      bundle(ik).resize (nSpin);

      for (iSpin=0; iSpin < nSpin; iSpin++) 
         bundle(ik)(iSpin).copy ( in.bundle(ik)(iSpin) );

   }
   return *this;
}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
template<class T>
SxDiracVec<T> &SxBundle3<T>::operator() (ssize_t iSpin, ssize_t ik)
{
   return bundle(ik)(iSpin);
}

template<class T>
const SxDiracVec<T> &SxBundle3<T>::operator() (ssize_t iSpin, ssize_t ik) const
{
   return bundle(ik)(iSpin);
}




template<class T>
typename T::Type &
SxBundle3<T>::operator() (ssize_t i, ssize_t iSpin, ssize_t ik)  
{
   return bundle(ik)(iSpin)(i);
}

template<class T>
const typename T::Type &
SxBundle3<T>::operator() (ssize_t i, ssize_t iSpin, ssize_t ik) const
{
   return bundle(ik)(iSpin)(i);
}




template<class T>
void SxBundle3<T>::set (const typename T::Type &s)
{
   ssize_t ik, nk, iSpin, nSpin;

   nk = bundle.getSize();

   for (ik=0; ik < nk; ik++)  {
      nSpin = bundle(ik).getSize();
      for (iSpin=0; iSpin < nSpin; iSpin++)  bundle(ik)(iSpin).set (s);
   }
}


template<class T>
int SxBundle3<T>::getNStates (int ik) const
{
   return (int)bundle(ik)(0).getSize();
}

template<class T>
int SxBundle3<T>::getNSpin (int ik) const
{
   return (int)bundle(ik).getSize();
}


template<class T>
void SxBundle3<T>::write (SxBinIO &io) const
{
   int iSpin, nSpin = bundle(0).getSize();
   int ik,    nk    = bundle.getSize();
#  ifdef _USE_H5   
      io.createGroup (groupName);
      for (ik=0; ik < nk; ik++)  {
         io.createGroup (groupName+"/k-"+ik);
         for (iSpin=0; iSpin < nSpin; iSpin++)  {
            io.writeData (groupName+"/k-"+ik+"/spin-"+iSpin, in(iSpin,ik));
         }
      }
#  else
      SX_CHECK (io.isOpen);
      fwrite (&nk,      sizeof(nk),      1, io.fp);
      fwrite (&nSpin,   sizeof(nSpin),   1, io.fp);
      int n;
      for (ik=0; ik < nk; ik++)  {
         n = getNStates(ik);
         fwrite (&n, sizeof(n), 1, io.fp); // nStates(ik)
         for (iSpin=0; iSpin < nSpin; iSpin++)  {
            fwrite (bundle(ik)(iSpin).elements, sizeof (typename T::Type), n, io.fp );
         }
      }
#  endif /* USE_H5 */      
}


template<class T>
void SxBundle3<T>::read (SxBinIO &io)
{
#  ifdef _USE_H5
      SX_EXIT;
#  else   
      SX_CHECK (io.isOpen);
      SX_CHECK (io.mode == SxBinIO::BINARY_READ_ONLY);
      int ik, nk, iSpin, nSpin, n;

      fread (&nk,      sizeof (nk),      1, io.fp);
      fread (&nSpin,   sizeof (nSpin),   1, io.fp);
      bundle.resize (nk);
      
      for (ik=0; ik < nk; ik++)  {
         bundle(ik).resize (nSpin);
         fread (&n, sizeof (n), 1, io.fp); // nStates (ik)
         
         for (iSpin=0; iSpin < nSpin; iSpin++)  {
            
            bundle(ik)(iSpin).resize (n);
            fread (bundle(ik)(iSpin).elements, sizeof (typename T::Type), n, io.fp );
         }
      }
#   endif /* _USE_H5 */       
}

template<class T>
int SxBundle3<T>::getNk () const
{
   return int(bundle.getSize ());
}

template <class T>
void SxBundle3<T>::synMPI ()  // LoopMPI
{
   int ik, nk, iSpin, nSpin;

   nk = int(bundle.getSize());

   SX_MPI_LEVEL("waves-k");
   SX_MPI_SOURCE(CurrentLevel, TaskGroupMaster);
   SX_MPI_TARGET(LevelSiblings, TaskGroupAll);
   for (ik=0; ik < nk; ik++)  {
      nSpin = int(bundle(ik).getSize());
      for (iSpin=0; iSpin < nSpin; iSpin++) {
         SxLoopMPI::bcast (bundle(ik)(iSpin), SxLoopMPI::whoseWork(ik));
      }
   }
}


//------------------------------------------------------------------------------
template<class T>
SxBinIO &operator<< (SxBinIO &io, const SxBundle3<T> &bundle)
{
   bundle.write (io);
   return io;
}



template<class T>
SxBinIO &operator>> (SxBinIO &io, SxBundle3<T> &bundle)
{
   bundle.read (io);
   return io;
}
