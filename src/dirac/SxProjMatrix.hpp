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
SxDiracVec<T>
SxProjMatrix<T>::getProjectedSingle (int iStart, int iEnd,
                                     const SxDiracVec<T> &projBlock,
                                     const SxDiracVec<T> &right) const
{
   SX_CHECK (iEnd - iStart + 1 == projBlock.nCols (),
             iEnd, iStart, projBlock.nCols ());
   // return (projBlock.adjoint () ^ right); // XPress
   return (right.conj () ^ projBlock).conj ()
          .reshape (projBlock.nCols ());
}

template<class T>
SxDiracVec<T>
SxProjMatrix<T>::SaveProjections::getProjectedSingle (int iStart, int iEnd,
                                     const SxDiracVec<T> &projBlock,
                                     const SxDiracVec<T> &right) const
{
   SxDiracVec<T> res;  
   if (cacheUsed)  {
      res.copy (savedProjections(SxIdx(iStart, iEnd)));
      res.reshape (projBlock.nCols ());
   } else {
      res = SxProjMatrix<T>::getProjectedSingle (iStart, iEnd, projBlock,
                                                 right);
      if (useCache)
         savedProjections(SxIdx(iStart, iEnd)) = res;
      return res;
   }
   return res;
}

template<class T>
typename T::Type
SxProjMatrix<T>::getProjectedSingle (int,
                                     const SxDiracVec<T> &projector,
                                     const SxDiracVec<T> &right) const
{
   return dot (projector, right);
}


template<class T>
typename T::Type
SxProjMatrix<T>::SaveProjections::getProjectedSingle (
      int iProj,
      const SxDiracVec<T> &projector,
      const SxDiracVec<T> &right) const
{
   if (cacheUsed)  {
      return savedProjections(iProj);
   } else {
      // perform projection
      typename T::Type proj = dot (projector, right);
      if (useCache) savedProjections(iProj) = proj;
      return proj;
   }
}
      
template <class T>
int
SxProjMatrix<T>::SaveProjections::prepare (const SxDiracVec<T> &right) const
{
   int nR = (int)right.nCols ();
   if (nR == 0 && right.getSize () == nElements) nR = 1;
   if (useCache)  {
      savedProjections = SxDiracVec<T> (nProj * nR);
      savedProjections.reshape (nProj, nR);
      cacheUsed = false;
   } else if (cacheUsed) {
      nR = (int)savedProjections.nCols ();
   }
   return nR;
}

template <class T>
void SxProjMatrix<T>::SaveProjections::finalize () const
{
   if (useCache) cacheUsed = true;
}

template<class T>
void SxProjMatrix<T>::applyAndAddSingle   (      SxDiracVec<T> *left,
                                           const SxDiracVec<T> &right) const
{
   SX_CHECK (left);
   SX_CHECK (left->getSize () == nElements, left->getSize (), nElements);
   
   // hook in for derived classes
   prepare (right);

   int ip, iBlock = 0;
   SxDiracVec<T> projBlock, projected, projCol;
   projBlock.reformat (nElements, blocksize);

   for (ip=0; ip < nProj; ++ip)  {
      
      // collect projectors
      projCol = projBlock.colRef(iBlock++);
      getProjector(ip, &projCol);

      // --- when block is full
      if (iBlock == projBlock.nCols () /* == current blocksize */)  {
         if (iBlock > 1)  {
            projected = getProjectedSingle (ip-iBlock+1, ip, projBlock, right);
            
            if (useFactors ())  {
               // multiply with inner factor
               typename SxDiracVec<T>::Iterator pIt = projected.begin ();
               for (int i = ip - iBlock + 1; i <= ip; ++i)
                  *pIt++ *= getFactor (i);
            }

            // add to left side
            *left += projBlock ^ projected;
         } else {
            typename T::Type proj = getProjectedSingle (ip, projBlock, right);
            // multiply with inner factor
            proj *= getFactor (ip);

            // add to left side
            left->plus_assign_ax (proj, projBlock);
         }
         
         iBlock = 0;
         // reformat next block
         if (blocksize + ip >= nProj && ip + 1 < nProj)
            projBlock.reformat (nElements, nProj - ip - 1);
      }
   }
   // hook in for derived classes
   finalize ();
}

template<class T>
SxDiracVec<T>
SxProjMatrix<T>::getProjectedMultiple (int iStart, int iEnd,
                                     const SxDiracVec<T> &projBlock,
                                     const SxDiracVec<T> &right) const
{
   SX_CHECK (iEnd - iStart + 1 == projBlock.nCols (),
             iEnd, iStart, projBlock.nCols ());
   //return (projBlock.adjoint () ^ right);
   // note: "right" may have more rows than we use
   return projBlock.overlap (right, projBlock.nRows ());
}

template<class T>
SxDiracVec<T>
SxProjMatrix<T>::SaveProjections::getProjectedMultiple (int iStart, int iEnd,
                                     const SxDiracVec<T> &projBlock,
                                     const SxDiracVec<T> &right) const
{
   SxDiracVec<T> res;  
   if (cacheUsed)  {
      int nCols = (int)savedProjections.nCols ();
      res.reformat (iEnd-iStart+1, nCols);
      int offset;
      SxIdx idx;
      for (int ic = 0; ic < nCols; ++ic)  {
         offset = nProj * ic;
         idx.start = offset + iStart;
         idx.end   = offset + iEnd;
         res.colRef(ic) <<= savedProjections (idx);
      }
   } else {
      res = SxProjMatrix<T>::getProjectedMultiple (iStart, iEnd, projBlock,
                                                 right);
      if (useCache)  {
         int nCols = (int)res.nCols ();
         int offset;
         SxIdx idx;
         for (int ic = 0; ic < nCols; ++ic)  {
            offset = nProj * ic;
            idx.start = offset + iStart;
            idx.end   = offset + iEnd;
            savedProjections(idx) = res.colRef(ic);
         }
      }
      return res;
   }
   return res;
}

template<class T>
SxDiracVec<T>
SxProjMatrix<T>::getProjectedMultiple (int,
                                     const SxDiracVec<T> &projector,
                                     const SxDiracVec<T> &right) const
{
   //return projector.conj ().reshape (1, nElements) ^ right;
   // note: "right" may have more rows than we use
   SxDiracVec<T> res = right.overlap (projector, projector.nRows ());
   // form conjugate in place
   for (int i = 0; i < res.getSize (); i++)
      res(i).im = -res(i).im;
   return res;
}


template<class T>
SxDiracVec<T>
SxProjMatrix<T>::SaveProjections::getProjectedMultiple (
      int iProj,
      const SxDiracVec<T> &projector,
      const SxDiracVec<T> &right) const
{
   SxDiracVec<T> res;  
   if (cacheUsed)  {
      // get a single row
      res = savedProjections.row (iProj);
   } else {
      // perform projection
      res = SxProjMatrix<T>::getProjectedMultiple(iProj, projector, right);
      if (useCache)   {
         int nCols = (int)right.nCols ();
         // set a single row
         typename SxDiracVec<T>::Iterator spIt, resIt;
         spIt = savedProjections.begin (); spIt += iProj;
         resIt = res.begin ();
         for (int ic = 0; ic < nCols; ++ic, spIt += nProj, ++resIt)
            *spIt = *resIt;
      }
   }
   return res;
}
      
template<class T>
void SxProjMatrix<T>::applyAndAddMultiple (      SxDiracVec<T> *left,
                                           const SxDiracVec<T> &right) const
{
   SX_CHECK (left);
   
   int nCols = prepare (right);

   SX_CHECK (nCols > 1);
   SX_CHECK (left->nRows () == nElements, left->nRows (), nElements);
   SX_CHECK (left->nCols () == nCols, left->nCols (), nCols);

   int ip, iBlock = 0;
   SxDiracVec<T> projBlock, projected, projCol;
   projBlock.reformat (nElements, blocksize);

   for (ip=0; ip < nProj; ++ip)  {
      
      // collect projectors
      projCol = projBlock.colRef(iBlock++);
      getProjector(ip, &projCol);

      // --- when block is full
      if (iBlock == projBlock.nCols () /* == current blocksize */)  {
         if (iBlock > 1)  {
            projected = getProjectedMultiple (ip - iBlock + 1, ip, projBlock,
                                              right);
            
            // get inner factor
            if (useFactors ())  {
               SxDiracVec<T> diag(iBlock);
               typename SxDiracVec<T>::Iterator diagIt = diag.begin ();
               for (int i = ip - iBlock + 1; i <= ip; ++i)
                  *diagIt++ = getFactor (i);
               // multiply with inner factor
               for (int ic = 0; ic < nCols; ++ic)
                  projected.colRef(ic) *= diag;
            }

            // add to left side
            *left += projBlock ^ projected;
         } else {
            projected = getProjectedMultiple (ip, projBlock, right);
            
            // multiply with inner factor
            typename T::Type factor = getFactor (ip);

            // add to left side
            typename SxDiracVec<T>::Iterator pIt = projected.begin ();
            for (int ic = 0; ic < nCols; ++ic, ++pIt)
               left->colRef(ic).plus_assign_ax(factor * *pIt,
                                               projBlock);
         }
         
         iBlock = 0;
         // reformat next block
         if (blocksize + ip >= nProj && ip + 1 < nProj)
            projBlock.reformat (nElements, nProj - ip - 1);
      }
   }
   // hook in for derived classes
   finalize ();
}

template<class T>
void SxProjMatrix<T>::projectionKernel (SxDiracVec<T> *res,
                                        const SxDiracVec<T> &left,
                                        const SxDiracVec<T> &right) const
{
   int nL = (int)left.nCols ();
   if (nL == 0 && left.getSize () == nElements) nL = 1;
   SX_CHECK (nL == 0 || res, nL);
   int nR = prepare (right);

   if (res)  {
      if (res->getSize () > 0 && res->nRows () == nL && res->nCols () == nR)  {
         res->set (0.);
      } else {
         *res = SxDiracVec<T> ();
      }
   }

   int ip, iBlock = 0;
   SxDiracVec<T> projBlock, projL, projR, projCol;
   projBlock.reformat (nElements, blocksize);

   for (ip=0; ip < nProj; ++ip)  {
      
      // collect projectors
      projCol = projBlock.colRef(iBlock++);
      getProjector(ip, &projCol);

      // --- when block is full
      if (iBlock == projBlock.nCols () /* == current blocksize */)  {
         if (iBlock > 1)  {
            // --- perform projection
            projR = getProjectedMultiple (ip-iBlock+1,ip,projBlock,right);
            projR.reshape (iBlock, nR);

            if (nL > 0)  {
               
               if (&left == &right)
                  projL = projR;
               else
                  //projL = projBlock.adjoint () ^ left;
                  projL = projBlock.overlap (left);
               // may not be the case for nL == 1
               projL.reshape (iBlock, nL);
               
               // --- calculate contribution to matrix

               if (useFactors ())  {
                  // get inner factor
                  SxDiracVec<T> diag(iBlock);
                  typename SxDiracVec<T>::Iterator diagIt = diag.begin ();
                  for (int i = ip - iBlock + 1; i <= ip; ++i)
                     *diagIt++ = getFactor (i);
                  // multiply with inner factor
                  for (int ic = 0; ic < nR; ++ic)
                     projR.colRef(ic) *= diag;
               }

               // get contribution of this block
               SxDiracVec<T> LR = projL.overlap (projR);

               // add to result
               if (res->getSize () == 0)
                  *res = LR;
               else
                  *res += LR;
            }
         } else {
            // --- perform projection
            projR = getProjectedMultiple (ip, projBlock, right);
            projR.reshape (1, nR);
            if (nL > 0)  {
               if (&left == &right) 
                  projL = projR.adjoint ();
               else
                  //projL = (projBlock.adjoint () ^ left).adjoint ();
                  projL = (projBlock.overlap (left)).adjoint ();
               projL.reshape (nL, 1);
               
               // --- calculate contribution to matrix
               
               // get inner factor
               typename T::Type factor = getFactor (ip);

               if (res->getSize () == 0)  {
                  res->reformat (nL, nR);
                  res->set (0.);
               }

               // add contribution of this block to result
               if (nL == 1)  {
                  res->plus_assign_ax(projL(0) * factor, projR);
               } else if (nR == 1)  {
                  res->plus_assign_ax(projR(0) * factor, projL);
               } else  {
                  // result += factor * (projL ^ projR);
                  res->plus_assign_ax (factor, projL ^ projR);
               }
            }
         }
         
         iBlock = 0;
         // reformat next block
         if (blocksize + ip >= nProj && ip + 1 < nProj)
            projBlock.reformat (nElements, nProj - ip - 1);
      }
   }
   finalize ();
}


template <class T>
void SxProjMatrix<T>::getProjector(int i, SxDiracVec<T> *target) const
{
   SX_CHECK (target);
   (*target) <<= getProjector (i);
}


template <class T>
void SxProjMatrix<T>::SaveProjections::clearCache ()
{ 
   savedProjections = SxDiracVec<T> (); 
   cacheUsed = false;
}


template <class T>
SxDiracVec<T> 
SxProjMatrix<T>::SaveProjections::getProjection(const SxDiracVec<T> &right)
{
   SX_CHECK (right.getSize () > 0);
   SX_CHECK(right.nRows () == nElements, right.nRows (), nElements);
   useCache = true;
   clearCache ();
   projectionKernel (NULL, SxDiracVec<T> (), right);
   SX_CHECK (cacheUsed);
   return savedProjections;
}

template <class T>
SxDiracVec<T>
SxProjMatrix<T>::SaveProjections::getProjectionFromExtended(const SxDiracVec<T> &right)
{
   SX_CHECK (right.getSize () > 0);
   SX_CHECK(right.nRows () >= nElements, right.nRows (), nElements);
   useCache = true;
   clearCache ();
   projectionKernel (NULL, SxDiracVec<T> (), right);
   SX_CHECK (cacheUsed);
   return savedProjections;
}

template<class T>
SxDiracVec<T>
SxProjMatrix<T>::SaveProjections::gradient (const SxDiracVec<T> &proj)
{
   SX_CHECK (proj.getSize () > 0 || cacheUsed);
   if (proj.getSize () > 0)  {
      SX_CHECK (proj.nRows () == nProj, proj.nRows (), nProj);
      savedProjections = proj;
      cacheUsed = true;
   }
   useCache = false;
   int nCols = (int)savedProjections.nCols ();
   if (nCols == 0) nCols = 1;
   SxDiracMat<T> res (nElements, nCols);
   res.set (0.);

   if (nCols == 1)
      applyAndAddSingle (&res, SxDiracVec<TPrecCoeffG> ());
   else
      applyAndAddMultiple (&res, SxDiracVec<TPrecCoeffG> ());

   return res;
}
