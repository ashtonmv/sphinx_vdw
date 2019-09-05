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

template<class T>
SxTensor3<T>::SxTensor3 ()
{
   nk = nn = nm = -1;
}
   
template<class T>
SxTensor3<T>::SxTensor3 (const int nk_, const int nn_, const int nm_)
{
   nk = nk_;  nn = nn_;  nm = nm_;

   tensor.resize (nk);

   int k;
   for (k = 0; k < nk; k++)  {
      tensor(k) = SxMatrix<T> (nn,nm);
//      tensor(k).reformat (nn,nm);
   }
}

template<class T>
SxTensor3<T>::SxTensor3 (const SxTensor3<T> &t)
{
   nk = t.nk, nn = t.nn, nm = t.nm;

   SX_CHECK (t.tensor.getSize () == nk, t.tensor.getSize (), nk);

   tensor.resize (nk);

   int k;
   for (k = 0; k < nk; k++)  {
      tensor(k) = SxMatrix<T> (nn,nm);
//      tensor(k).reformat (nn,nm);
      tensor(k).copy ( t(k) );
   }
}

template<class T>
void SxTensor3<T>::reformat (const int nk_, const int nn_, const int nm_)
{
   nk = nk_;  nn = nn_;  nm = nm_;

   tensor.resize (nk);

   int k;
   for (k = 0; k < nk; k++)  {
      tensor(k) = SxMatrix<T> (nn,nm);
//      tensor(k).reformat (nn,nm);
   }
}

template<class T>
void SxTensor3<T>::set (const SxMatrix<T> &mat, const int &k)
{
   SX_CHECK (k >= 0 && k < nk, k, nk);

//   tensor(k).copy (mat);
   tensor(k).set (mat);
}

template<class T>
void SxTensor3<T>::set (const typename T::Type &elem,
                        const int &k, const int &n, const int &m)
{
   SX_CHECK (k >= 0 && k < nk, k, nk);
   SX_CHECK (n >= 0 && n < nn, n, nn);
   SX_CHECK (m >= 0 && m < nm, m, nm);

   tensor(k)(n,m) = elem;
}

template<class T>
void SxTensor3<T>::setAllElem (const typename T::Type &in)
{
   SX_CHECK (tensor.getSize () == nk, tensor.getSize (), nk);

   int k;
   for (k = 0; k < nk; k++)  {
      tensor(k).set (in);
   }
}

template<class T>
void SxTensor3<T>::setDiagonals (const typename T::Type &in)
{
   SX_CHECK (tensor.getSize () == nk, tensor.getSize (), nk);
   SX_CHECK (nn == nm, nn, nm);  // otherwise not quadratic

   SxVector<T> diag(nn,in);

   int k;
   for (k = 0; k < nk; k++)  {
      tensor(k) = tensor(k).identity (diag);
   }
}

template<class T>
SxMatrix<T> &SxTensor3<T>::operator() (const int &k)
{
   SX_CHECK (k >= 0 && k < nk, k, nk);
   return tensor(k);
}

template<class T>
const SxMatrix<T> &SxTensor3<T>::operator() (const int &k) const
{
   SX_CHECK (k >= 0 && k < nk, k, nk);
   return tensor(k);
}

template<class T>
typename T::Type &SxTensor3<T>::operator() (const int &k,
                                            const int &n, const int &m)
{
   SX_CHECK (k >= 0 && k < nk, k, nk);
   SX_CHECK (n >= 0 && n < nn, n, nn);
   SX_CHECK (m >= 0 && m < nm, m, nm);
   return tensor(k)(n,m);
}

template<class T>
const typename T::Type &SxTensor3<T>::operator() (const int &k,
                                                  const int &n,
                                                  const int &m) const
{
   SX_CHECK (k >= 0 && k < nk, k, nk);
   SX_CHECK (n >= 0 && n < nn, n, nn);
   SX_CHECK (m >= 0 && m < nm, m, nm);
   return tensor(k)(n,m);
}

template<class T>
SxTensor3<T> SxTensor3<T>::operator= (const SxTensor3<T> &in)
{
   if (this == &in)  return *this;

   reformat (in.nk, in.nn, in.nm);

//   for (int ik = 0; ik < nk; ik++)
//      tensor(ik) = SxMatrix<T> (in(ik));  // references the matrices
   for (int ik = 0; ik < nk; ik++)
      tensor(ik).set (in(ik));  // without referencing

   return *this;
}

template<class T>
void SxTensor3<T>::print ()
{
   print (nk);
}

template<class T>
void SxTensor3<T>::print (const int &z)
{
   SX_CHECK (z > 0 && z <= nk, z, nk);

   cout << SX_SEPARATOR;
   if (z < nk)  {
      printf ("SxTensor3  ---  %d of %d matrices printed:\n\n", z, nk);
   }

   int k;
   for (k = 0; k < z; k++)  {
      printf ("k = %d:\n\n", k);
      cout << tensor(k) << endl;
   }

   cout << SX_SEPARATOR;
}

template<class T>
void SxTensor3<T>::print (const SxString &s)
{
   print (s, nk);
}

template<class T>
void SxTensor3<T>::print (const SxString &s, const int &z)
{
   SX_CHECK (z > 0 && z <= nk, z, nk);

   cout << SX_SEPARATOR;

   printf ("SxTensor3:  %s  ---  %d of %d matrices printed:\n\n",
           s.ascii (), z, nk);

   int k;
   for (k = 0; k < z; k++)  {
      printf ("%s --> k = %d:\n\n", s.ascii (), k);
      cout << tensor(k) << endl;
      if (k != z-1)  {
         cout << "- - - - - - - - - - - - - - - - - - - - - - -" << endl;
      }  else  {
         cout << SX_SEPARATOR;
      }
   }
}

template<class T>
void SxTensor3<T>::writeTensor3 (const SxString &file) const
{
   try  {
      SxBinIO io (file, SxBinIO::BINARY_WRITE_ONLY);
      writeTensor3 (io);
      io.setMode (SxBinIO::WRITE_DATA);
      writeTensor3 (io);
      io.close ();
   }  catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }
}

template<class T>
void SxTensor3<T>::writeTensor3 (const SxBinIO &io) const
{
   int iElem, nElem = nk * nm * nn;
   int ik, im, in;
   SxVector<T> kmnElems(nElem);

   // --- encode data: tensor --> one very long vector
   iElem = 0;
   for (ik = 0; ik < nk; ik++)  {
      for (im = 0; im < nm; im++)  {
         for (in = 0; in < nn; in++)  {
            kmnElems(iElem) = tensor(ik)(im,in);
            iElem++;
         }
      }
   }

   SX_CHECK (iElem == nElem, iElem, nElem);

   try  {
      // --- write header, create dimensions
      io.addDimension ("nk", nk);
      io.addDimension ("nm", nm);
      io.addDimension ("nn", nn);
      io.addDimension ("nElem", nElem);

      // --- write data
      io.write ("kmnElems", kmnElems, "nElem");
   }  catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }
}

template<class T>
void SxTensor3<T>::readTensor3 (const SxString &file)
{
   try  {
      SxBinIO io (file, SxBinIO::BINARY_READ_ONLY);
      readTensor3 (io);
      io.close ();
   }  catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }
}

template<class T>
void SxTensor3<T>::readTensor3 (const SxBinIO &io)
{
   int iElem, nElem;
   int ik, im, in, nkIn, nmIn, nnIn;
   SxVector<T> kmnElems;

   try  {
      nkIn  = io.getDimension ("nk");
      nmIn  = io.getDimension ("nm");
      nnIn  = io.getDimension ("nn");
      nElem = io.getDimension ("nElem");

      kmnElems.resize (nElem);

      io.read ("kmnElems", &kmnElems, nElem);
   }  catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }

   if (nkIn * nmIn * nnIn != nElem)  {
      sxprintf ("| Error:  file \"%s\" corrupt, nk x nm x nn = %d, but "
                "nElem = %d\n",
                io.filename.ascii(), nkIn * nmIn * nnIn, nElem);
      SX_QUIT;
   }

   reformat (nkIn, nmIn, nnIn);

   // --- decode data
   iElem = 0;
   for (ik = 0; ik < nk; ik++)  {
      for (im = 0; im < nm; im++)  {
         for (in = 0; in < nn; in++)  {
            tensor(ik)(im,in) = kmnElems(iElem);
            iElem++;
         }
      }
   }

   SX_CHECK (iElem == nElem, iElem, nElem);
}

template<class T>
void SxTensor3<T>::readUkmn (const SxString &file)
{
   int iElem, nElem;
   int ik, im, in, nkIn, dimIn;
   SxVector<T> kmnElems;

   try  {
      SxBinIO io (file, SxBinIO::BINARY_READ_ONLY);

      nkIn  = io.getDimension ("nk");
      dimIn = io.getDimension ("dim");
      nElem = io.getDimension ("nElem");

      kmnElems.resize (nElem);

      io.read ("kmnElems", &kmnElems, nElem);

      if (nkIn * dimIn * dimIn != nElem)  {
         sxprintf ("| Error:  file \"%s\" corrupt, nk x dim x dim = %d, but "
                   "nElem = %d\n",
                   io.filename.ascii(), nkIn * dimIn * dimIn, nElem);
         SX_QUIT;
      }

      io.close ();
   }  catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }

   reformat (nkIn, dimIn, dimIn);

   // --- decode data
   iElem = 0;
   for (ik = 0; ik < nk; ik++)  {
      for (im = 0; im < dimIn; im++)  {
         for (in = 0; in < dimIn; in++)  {
            tensor(ik)(im,in) = kmnElems(iElem);
            iElem++;
         }
      }
   }

   SX_CHECK (iElem == nElem, iElem, nElem);
}
