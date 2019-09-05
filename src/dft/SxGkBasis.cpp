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

#include <SxUtil.h>
#include <SxGkBasis.h>
#include <SxConstants.h>
#include <SxError.h>
//#include <main.h>
#include <stdio.h>
#include <SxRBasis.h>

#include <SxLoopMPI.h>


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
SxGkBasis::SxGkBasis () 
   : SxKPoints (),
     gBasis(NULL)
{
   nGkMin = nGkMax = nGMax = nGIndp = -1;
   gkCut = -1.;
   registerMemoryObservers ();
}




SxGkBasis::SxGkBasis (const SxGBasis &G_, 
                      const SxSymbolTable *table)
   : SxKPoints (G_.structPtr->cell, table)
{
   SX_CHECK (table);
   gBasis = &G_;
   
   // --- basis
   SxSymbolTable *basis;
   bool saveMemory = false;
   try {
      basis = table->getGroup ("basis");
      gkCut = basis->get("eCut")->toReal();
      
      if (basis->contains("saveMemory"))
         saveMemory = basis->get("saveMemory")->toAttribute ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   init (saveMemory);
   registerMemoryObservers ();
}

SxGkBasis::SxGkBasis (const SxKPoints &kp,
                      const SxGBasis &G_,
                      double cutoff,
                      bool saveMemory)
   : SxKPoints (kp),
     gBasis(&G_),
     gkCut(cutoff)
{
   init (saveMemory);
   registerMemoryObservers ();
}

SxGkBasis::SxGkBasis (const SxKPoints &kp,
                      const SxAtomicStructure &structure,
                      const SxMesh3D &mesh,
                      double cutoff,
                      bool saveMemory)
   : SxKPoints (kp),
     gBasis(NULL),
     gkCut(cutoff)
{
   setupBasisList (structure, mesh, saveMemory);
   print();
   registerMemoryObservers ();
}

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
SxGkBasis::~SxGkBasis ()
{
   // empty
}

void SxGkBasis::init (bool saveMemory)
{  
   cout << SX_SEPARATOR;
   cout << "| |G+k> Basis\n";
   cout << SX_SEPARATOR;

   SX_CHECK(gBasis);
   SX_CHECK(gBasis->structPtr);
   setupBasisList (*gBasis->structPtr, gBasis->fft3d(0).mesh, saveMemory);
   print ();
}


void SxGkBasis::changeTau (const SxAtomicStructure &tauList)
{
   SxList<SxGBasis *>::Iterator it;
   SX_MPI_LEVEL("waves-k");
   for (int ik=0; ik < nk; ik++)  {
      if (SxLoopMPI::myWork(ik))
         gBasisList(ik)->changeTau (tauList);
   }
}

SxGBasis &SxGkBasis::operator() (int ik)
{
   SX_CHECK (ik >= 0 && ik < nk, ik, nk);
   return *gBasisList(ik);
}

const SxGBasis &SxGkBasis::operator() (int ik) const
{
   SX_CHECK (ik >= 0 && ik < nk, ik, nk);
   return *gBasisList(ik);
}


const SxList<int> SxGkBasis::ngPerK () const
{
   SxList<int> ngList;
   for (int ik=0; ik < nk; ik++)  {
      ngList << gBasisList(ik)->ng;
   }
   return ngList;
}


const SxList<int> SxGkBasis::ngPerK (PrecEnergy eCut) const
{
   int ik, igk, ngk, n;
   SxPtr<SxGBasis> g;
   
   SxList<int> nPerK;
   for (ik=0; ik < nk; ik++)  {
      g = gBasisList(ik);
      ngk = g->ng;
      for (igk=0, n=0; igk < ngk; igk++)  
         if (g->g2(igk) < eCut)  n++;
      nPerK.append (n);
   }
   return nPerK;
}




void SxGkBasis::setupBasisList (const SxAtomicStructure &structure, const SxMesh3D &mesh, bool saveMemory)
{
   gBasisList.resize (nk);
   for (int ik=0; ik < nk; ik++)  {
      SxGBasis::MemoryUsage memMode = saveMemory ? SxGBasis::SaveMemory
                                                 : SxGBasis::SaveTime;
      SX_MPI_LEVEL("waves-k");
      if (SxLoopMPI::myWork(ik)) {
         gBasisList(ik) = SxPtr<SxGBasis>::create (mesh,
                                                   structure,
                                                   gkCut,
                                                   getK (ik), 
                                                   false,
                                                   memMode);
      } else {
         gBasisList(ik) = SxPtr<SxGBasis>::create ();
      }
   }
#ifdef USE_LOOPMPI
   {
      // get the number of g vectors for all k's
      SX_MPI_LEVEL("waves-k");
      for (int ik=0; ik < nk; ik++)  {
         int &ng = gBasisList(ik)->ng;
         ng = SxLoopMPI::bcast (ng, SxLoopMPI::whoseWork(ik));
      }
   }
#endif

   // --- get nGkMin, nGkMax
   nGkMin = gBasisList(0)->ng;
   nGkMax = 0;
   for (int ik = 0; ik < nk; ik++)  {
      nGkMax = maximum ( nGkMax , (gBasisList)(ik)->ng );
      nGkMin = minimum ( nGkMin , (gBasisList)(ik)->ng );
   }
   // sxprintf ("| Size of smallest |G+k> basis: nGkMin = %d\n", nGkMin);
   // sxprintf ("| Size of largest |G+k> basis:  nGkMax = %d\n", nGkMax);
   // cout << SX_SEPARATOR;
}

void SxGkBasis::setupN123inv ()
{
   SX_NO_MPI;
   // --- setup n123inv tables
   int ik;
   for (ik = 0; ik < nk; ik++)  {
      gBasisList(ik)->n123invSetup();
      sxprintf ("| ik = %d: n123inv --> setup done\n", ik);
   }
}

void SxGkBasis::setupNGkMinMax ()
{
   SX_NO_MPI;
   // --- get nGkMin, nGkMax
   nGkMin = gBasisList(0)->ng;
   nGkMax = 0;
   for (int ik = 0; ik < nk; ik++)  {
      nGkMax = maximum ( nGkMax , (gBasisList)(ik)->ng );
      nGkMin = minimum ( nGkMin , (gBasisList)(ik)->ng );
   }
   sxprintf ("| Size of smallest |G+k> basis: nGkMin = %d\n", nGkMin);
   sxprintf ("| Size of largest |G+k> basis:  nGkMax = %d\n", nGkMax);
   cout << SX_SEPARATOR;
}

/*
void SxGkBasis::setupLookupTables ()
{
   int ng; 
   int ik, ig, igp, L;
   SxPtr<SxGBasis> Gk;
   SxVector3<TPrecG> k, G;
   
   cout.flush();
   nGkMin = gBasis->ng; 
   nGkMax = 0;
   for (ik = 0; ik < nk; ik++) {
      nGkMax = maximum ( nGkMax , (gBasisList)(ik)->ng );
      nGkMin = minimum ( nGkMin , (gBasisList)(ik)->ng );
   }
   cout << "| Basis size:                         " 
        << "min: " << nGkMin << ", max: " << nGkMax << endl;
   cout << "| Constructing lookup tables:\n";
   
   cout.flush();
   nGMax = 0;
   GkInG.resize (nk * nGkMax); GkInG.reshape (nGkMax, nk);
   for (ik = 0; ik < nk; ik++) {
      k  = kVec(ik);
      Gk = (gBasisList)(ik);
      ng = Gk->ng;
      
      for (ig = 0; ig < ng; ig++) {
         G = Gk->getG(ig);
         L = lookup(G-k, gBasis);
         GkInG(ig,ik) = L ;
         if (L >= nGMax) nGMax = L+1;
      }
   }
   cout << "|                       Required " << nGMax << " G vectors.\n";
   
   cout.flush();
   GplusGinG.resize (nGMax * nGMax); GplusGinG.reshape (nGMax, nGMax);
   GminusGinG.resize (nGMax * nGMax); GminusGinG.reshape (nGMax, nGMax);

   for (ig = 0; ig < nGMax; ig++) 
      for (igp = 0; igp < nGMax; igp++) {
         L = lookup (gBasis->getG(ig) + gBasis->getG(igp), gBasis);
         if ( L >= nGMax || L == -1 )  GplusGinG(ig,igp) = nGMax;
         else  GplusGinG(ig,igp) = L;
         GminusGinG(ig,igp) = lookup (gBasis->getG(ig) - gBasis->getG(igp), gBasis);
      }

   cout.flush();
   GinGk.resize (nk * (nGMax+1)); GinGk.reshape ((nGMax+1), nk);
   GinGk.set(-1);
   for (ik = 0; ik < nk; ik++) {
      ng = gBasisList(ik)->ng;
      for (ig = 0; ig < ng; ig++) {
         GinGk(GkInG(ig,ik), ik) = ig;
      }
   }

   SxVector3<TPrecG> Gp;
   int iOp, nOp = geom->symOp.getSize(); 
   int currentIdx;
   int jg, i, j, ij;
   ng = gBasis->ng;  
   GIndx.resize(ng);   GIndx.set(-1); 
   GIndp.resize(ng);   nGIndp = ng;  
   nStars.resize(ng);  nStars.set(1); 
   MGIndx.reformat (ng,ng); MGIndx.set (-1);
   MnStars.reformat(ng,ng); MnStars.set(1); 
   MGIndp.resize(ng*ng);
   for (ig = 0;  ig < ng; ig++) GIndp(ig) = ig;  

   for (jg = 0;  jg < ng; jg++) 
      for (ig = 0;  ig < ng; ig++)  {
         i = ig + jg*ng;
         MGIndp(i) = i;  
      }
   
   if ( nOp > 1 )  { 
      SxMatrix3<Double> aMat = geom->aMat; 
      SxMatrix3<Double> bMat = (aMat.inverse());
      const SymOpList &sOp = geom->symOp;
      SxList<SxMatrix3<Double> > symOp;   symOp.resize(nOp); 
//    (aMat^bMat).print(); 
      
      for (iOp = 0; iOp < nOp; iOp++)  {
         //symOp(iOp) =  ( aMat ^  SxMatrix3<Double>(sOp(iOp))  ^ bMat).inverse(); 
         symOp(iOp) =  ( aMat ^  sOp(iOp) ^ bMat );
      }
      
      L = 0; 
      for (ig = 0; ig < ng;  ig++) { 
         if ( GIndx(ig) == -1 ) {
            Gp = gBasis->getG(ig); 
            GIndx(ig) = GIndp(L++) = ig;  
            for (iOp = 0; iOp < nOp; iOp++) { 
               G = symOp(iOp) ^ Gp; 
               ik = lookup( G, gBasis); 
               if (ik > ig && GIndx(ik) == -1)  {  
                  GIndx(ik) = ig; 
                  nStars(ig)++; nStars(ik)++; 
               } 
            }   
         }
      }   
      nGIndp = L;    GIndp.resize(nGIndp, true);    

//    cout << "nGIndp = " << nGIndp <<  endl; 
//    cout << "nGIndx = " << GIndx  <<  endl; 
//    cout << "GIndp  = " << GIndp; 

      SxMatrix<Int> idx(nOp,ng);
      for (ig = 0; ig < ng;  ig++) { 
         Gp = gBasis->getG(ig); 
         for (iOp = 0; iOp < nOp; iOp++)
            idx(iOp,ig) = lookup ( symOp(iOp) ^ Gp, gBasis);
      }
      
      L = 0;
      for (jg = 0; jg < ng;  jg++)  {
         for (ig = 0; ig < ng;  ig++) { 
            currentIdx = ig + jg*ng; // TODO: to vector class
            if (MGIndx(currentIdx) == -1)  {
               MGIndx(currentIdx) = MGIndp(L++) = currentIdx;
               for (iOp = 0; iOp < nOp; iOp++)  { // exclude 1
                  i  = idx(iOp,ig);
                  j  = idx(iOp,jg);
                  ij = i + j*ng;  // TODO: goes into the vector class
                  if (ij > currentIdx && MGIndx(i,j) == -1)  {  
                     MGIndx(i,j) = currentIdx; 
                     MnStars(i,j)++; MnStars(ig,jg)++;
                  }
               }
            } 
         } 
      }

      int MnGIndp = L;    MGIndp.resize(MnGIndp, true);    
//    sxprintf ("MnStars(%d) = %d\n", 2, MnStars(2));
//    cout << "1st row:\n" << MGIndx.colRef(0) << endl;
//    cout << "2nd row:\n" << MGIndx.colRef(1) << endl;
   } 
}*/


//----------------------------------------------------------------------------
bool SxGkBasis::equal(const SxVector3<TPrecG> &v1 , 
                      const SxVector3<TPrecG> &v2) const
{
   const Real8 DLT = 1e-7;
   return (  fabs(v1(0) - v2(0)) < DLT 
         &&  fabs(v1(1) - v2(1)) < DLT   
         &&  fabs(v1(2) - v2(2)) < DLT  ); 
}        

//----------------------------------------------------------------------------
bool SxGkBasis::equal(const PrecG &s1, const PrecG &s2) const
{
   const Real8 DLT = 1e-7;
   return ( fabs(s1-s2) < DLT );

}

 //----------------------------------------------------------------------------
int SxGkBasis::lookup (const SxVector3<TPrecG> &gVec, const SxGBasis *gk) const
{
   const Real8 DLT = 1e-7;
   SxDiracVec<TPrecG> gk2;      // |G+k| table
   PrecG            gVec2;
   int              ng, il, ih, i;         
   bool             not_found; 
   // --- Initializations ...
   gVec2 = (gVec^gVec);
   // --- retrieve some information about basis
   ng  = gk->ng;                // get the No. of G+k vectors;
   gk2 = gk->g2;                // get list of ALL |G+k|^2

   // --- Search for the shell (binary search tree algorithm)...
   if ( !(fabs (gVec2 - gk2(ng-1)) < DLT) && gVec2 > gk2(ng-1) ) {
      return -1 ;  //Vector is outside all shells.
   }
   il = 0; ih = ng-1;   
   do {
      i = (il+ih)/2;
      if (gVec2 < gk2(i)) ih = i-1;
      else             il = i+1;   
      not_found = !equal(gVec2, gk2(i));
   }
   while 	( il <= ih && not_found );
   
   if (not_found)  {
      return -1; // shell not found, otherwise:
   }
   // --- Search through the shell for the correct vector!   
   il = i; ih = i+1;
   
   // --- Go down until the begining of the shell    
   while ( il>0 && equal(gVec2,gk2(il)) && !equal(gVec,gk->getG(il)) ) 
      il--;
   if ( equal(gVec,gk->getG(il)) )  {
      return il ;    
   }
   
   // --- Go up until the end of the shell    
   while ( ih<ng-1 && equal(gVec2,gk2(ih)) && !equal(gVec,gk->getG(ih)) ) 
      ih++;
   if ( equal(gVec,gk->getG(ih)) ) {
      return ih ;    
   }

   // --- Not found
   return -1;
}


   
//----------------------------------------------------------------------------
int SxGkBasis::lookup (int ik, int igk, int iq, int igq)
{
   SX_CHECK (ik  >= 0 && ik < nk, ik, nk);
   SX_CHECK (iq  >= 0 && iq < nk, iq, nk);
   SX_CHECK  (igk >= 0, igk);
   SX_CHECK  (igq >= 0, igq);
   return GinGk(GplusGinG( GkInG(igk,ik), GkInG(igq,iq) ), iq );
}



//------------------------------------------------------------------------------
// Symetrization in reciprocal space. 
//------------------------------------------------------------------------------ 
void SxGkBasis::symmetrize_g (SxDiracVec<TPrecCoeffG> &vector) 
{ 
   int i, ig, ng = (int)vector.getSize(); 
   SxDiracVec<TPrecCoeffG> v(ng, 0.); 
   
   for (ig=0; ig<ng; ig++) { 
      i = GIndx(ig); 
      v(i) += vector(ig); 
   } 
   for (ig=0; ig<ng; ig++) { 
      i = GIndx(ig); 
      vector(ig) = v(i)/(PrecCoeffG)nStars(i); 
   } 
} 


void SxGkBasis::symmetrize_g (SxMatrix<TPrecCoeffG> &matrix) 
{
   SX_CHECK  (matrix.getSize() > 0, matrix.getSize());
   SX_CHECK (matrix.nRows() == matrix.nCols(), matrix.nRows(), matrix.nCols());

   int j, jg, i, idx, ig, n = (int)matrix.nRows(); 
   SxMatrix<TPrecCoeffG> m(n,n); m.set (0.);

   for (jg=1; jg < n; jg++) { 
      for (ig=1; ig < n; ig++) { 
         idx = MGIndx(ig,jg); 
         j   = idx / gBasis->ng; i = idx - j*gBasis->ng;
//         sxprintf ("%d,%d -> %d,%d, ", ig, jg, i,j);
//            cout << "|G|="    << gBasis->getG(i).absSqr() 
//                 << ", |G'|="  << gBasis->getG(j).absSqr()
//                 << ", |G*|="  << gBasis->getG(ig).absSqr()
//                 << ", |G'*|=" << gBasis->getG(jg).absSqr() << endl;
//
         m(i,j) += matrix(ig,jg); 
         SX_CHECK (i != 0);
         SX_CHECK (j != 0);
      }
   } 
   for (jg=1; jg < n; jg++) { 
      for (ig=1; ig < n; ig++) { 
         idx = MGIndx(ig,jg); 
         j   = idx / gBasis->ng; i = idx - j*gBasis->ng;
         matrix(ig,jg) = m(i,j)/(PrecCoeffG)MnStars(i,j); 
      }
   } 
} 



void SxGkBasis::write (const SxBinIO &io) const
{
   SX_MPI_LEVEL("waves-k");
   int ik, ig;
   SxVector<Int> nGk(nk);
   int nAllGk = 0, ng;
   
   // --- get nGk lengths
   for (ik = 0; ik < nk; ik++)
   {
      if (SxLoopMPI::myWork(ik))
      {
         nGk(ik) = ng = (gBasisList)(ik)->ng;
         nAllGk += ng;
      }
      else
      {
         nGk(ik) = 0;
      }

   }
   nAllGk = SxLoopMPI::sum(nAllGk);
   SxLoopMPI::sum (nGk);

   // --- support for parallel NetCDF4 IO -- testing is needed

   // Idea for parallelization:
   // -- copy data to contiguous buffers on each rank
   // -- use index arrays to handle the field-global and rank-local offsets
   // -- write with the NetCDF4 library using these index arrays

   // --- containers for k-points: gVec, n123 & source and dest iterators
   int nkRank = 0;
   int nGkRank = 0;
   for (ik = 0; ik < nk; ik++)
   {
      if (SxLoopMPI::myWork(ik))
      {
         nkRank += 1;
         nGkRank += nGk(ik);
      }
   }
   SxMatrix<TPrecG> gVectors(nGkRank, 3);
   SxMatrix<TPrecG> kVectors(nkRank, 3);
   SxVector<Int> n123(nGkRank);
   // --- offsets for writing at the correct location in the file
   SxVector<Int> gkOffset(nk); // global (row) offset
   SxVector<Int> gkRankOffset(nk); // local offset for the gk arrays
   SxVector<Int> kRankOffset(nk); // local offset for the k arrays

   // --- put kpoints into containers accross the MPI ranks
   int igTot = 0;
   int igRank = 0;
   int ikRank = 0;
   for (ik = 0; ik < nk; ik++)
   {
      ng = nGk(ik);
      gkOffset(ik) = igTot;
      //
      if (SxLoopMPI::myWork(ik))
      {
         kRankOffset(ik) = ikRank;
         gkRankOffset(ik) = igRank;
      }
      else
      {
         // dangerous, but useful for debug purposes
         kRankOffset(ik) = -1;
         gkRankOffset(ik) = -1;
      }
      //
      for (ig = 0; ig < ng; ig++, igTot++)
      {
         if (SxLoopMPI::myWork(ik))
         {
            n123(igRank) = gBasisList(ik)->n123(0)(ig);
            for (int dim = 0; dim < 3; dim++)
               gVectors(igRank, dim) = gBasisList(ik)->gVec(ig, dim);
            igRank++;
         }
      }
      if (SxLoopMPI::myWork(ik))
      {
         for (int dim = 0; dim < 3; dim++)
         {
            kVectors (ikRank, dim) = kVec(ik)(dim);
         }
         ikRank++;
      }
   }

   try  {
      // --- create dimensions
      io.addDimension ("nk",     nk);
      io.addDimension ("xyz",    3);
      io.addDimension ("nAllGk", nAllGk);
      
      // --- write data
      if ((io.ncMode == SxBinIO::WRITE_HEADER) || (SxLoopMPI::me() == 0))
         io.write ("eCut", gkCut);  // scalar, not parallel

      if ((io.ncMode == SxBinIO::WRITE_HEADER) || (SxLoopMPI::me() == 0))
         io.write ("meshDim", gBasisList(0)->fft3d(0).mesh, "xyz"); // 3-integer vector, not parallel

      if (!io.contains("nGk") || (io.ncMode == SxBinIO::WRITE_DATA))
      {
//         io.write ("nGk", nGk, "nk");
         for (ik = 0; ik < nk; ik++)
         {
            if ((io.ncMode == SxBinIO::WRITE_HEADER) || (SxLoopMPI::me() == 0))
            {
               io.write ("nGk", nGk, "nk");
               if (io.ncMode == SxBinIO::WRITE_HEADER)
                  break;
            }
         }
      }

//      io.write ("gkVec", gVectors, "nAllGk", "xyz");
      for (ik = 0; ik < nk; ik++)
      {
         if (SxLoopMPI::myWork(ik))
         {
            io.write ("gkVec", gVectors, "nAllGk", "xyz", gkOffset(ik), nGk(ik), gkRankOffset(ik));
            if (io.ncMode == SxBinIO::WRITE_HEADER)
               break;
         }
      }

//      io.write ("fftIdx", n123, "nAllGk" /*, rowOffset, colOffset */);
      for (ik = 0; ik < nk; ik++)
      {
         if (SxLoopMPI::myWork(ik))
         {
            io.write ("fftIdx", n123, "nAllGk", gkOffset(ik), nGk(ik), gkRankOffset(ik));
            if (io.ncMode == SxBinIO::WRITE_HEADER)
               break;
         }
      }

//      io.write ("kWeights", weights, "nk");
      for (ik = 0; ik < nk; ik++)
      {
         if (SxLoopMPI::myWork(ik))
         {
            io.write ("kWeights", weights, "nk", /* trivial: kOffset(ik) is simply */ ik, 1, /* kRankOffset(ik) */ ik);
            if (io.ncMode == SxBinIO::WRITE_HEADER)
               break;
         }
      }

      // redundant: can be generated from gkVec
//      io.write ("kVec", kVectors, "nk", "xyz");
      for (ik = 0; ik < nk; ik++)
      {
         if (SxLoopMPI::myWork(ik))
         {
            io.write ("kVec", kVectors, "nk", "xyz", /* trivial: kOffset(ik) is simply */ ik, 1, kRankOffset(ik));
            if (io.ncMode == SxBinIO::WRITE_HEADER)
               break;
         }
      }

      if ( ((io.ncMode == SxBinIO::WRITE_HEADER) || (SxLoopMPI::me() == 0))
            && (foldingMP(0) > 0 && foldingMP(1) > 0 && foldingMP(2) > 0) )
      {
         io.write ("foldingMP", foldingMP, "xyz"); // 3-integer vector, not parallel
      }

   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}



void SxGkBasis::read (const SxBinIO &io, bool initFFT, bool saveMemory)
{
   cout << "Reading |G+k> basis ..." << endl;

   SxKPoints::read (io);

   try  {
      // clean everything we don't have (and probably don't need)
      gBasis = NULL;

      GkInG = SxMatrix<Int> ();
      GplusGinG = SxMatrix<Int> ();
      GinGk = SxMatrix<Int> ();
      GminusGinG = SxMatrix<Int> ();

      nGkMax = -1;
      nGkMin = -1;
      nGMax = -1;
      GIndx = SxDiracVec< Int > ();
      GIndp = SxDiracVec< Int > ();
      nStars = SxDiracVec< Int > ();
      MGIndx = SxMatrix< Int > ();
      MnStars = SxMatrix< Int > ();
      MGIndp = SxDiracVec< Int > ();
      nGIndp = -1;

      // read cutoff energy
      io.read ("eCut", &gkCut);
      
      // --- read G bases data
      SxVector<Int> nGk(nk);
      io.read ("nGk", &nGk, nk);
      SxVector3<Int> mesh;
      io.read ("meshDim", &mesh);
      
      // --- set up G bases
      gBasisList.resize (nk);
      SxPtr<SxGBasis> g;
      int start = 0;

      SxDiracMat<TPrecG> gVecMatrixRef;

      for (int ik = 0; ik < nk; ik++)  {
         SxGBasis::MemoryUsage memMode = saveMemory ? SxGBasis::SaveMemory
                                                    : SxGBasis::SaveTime;

         if (!SxLoopMPI::myWork(ik))  {
            memMode = SxGBasis::SaveMemory;
         }
         gBasisList(ik) = g = SxPtr<SxGBasis>::create (memMode);

         if (SxLoopMPI::myWork(ik))  {
            g->read (io, nGk(ik), start,
                     initFFT ? mesh : SxVector3<Int> (0,0,0));
            if (!initFFT)  {
               // set up the mesh size
               g->fft3d.resize(1);
               g->fft3d(0).mesh = mesh;
               g->fft3d(0).meshSize = mesh.product ();
            }
         }

         // advance start index
         start += nGk(ik);

         cout << "k = " << kVec(ik) << " done." << endl;
      }

#ifndef NDEBUG
      int nAllGk = io.getDimension ("nAllGk");
      SX_CHECK (start == nAllGk, start, nAllGk);
#endif
      
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}


//------------------------------------------------------------------------------
void SxGkBasis::print () const
{
   SxString x, y, z, w, g, l;
   SxVector3<Double> k;
   // --- Monkhorst-Pack folding
   if (foldingMP(0) > 1 || foldingMP(1) > 1 || foldingMP(2) > 1) {
      SxCell recCell = getTau ().cell.getReciprocalCell ();
      cout << "| Using Monkhorst-Pack mesh\n";
      cout << "|    folding:                         ";
      cout << foldingMP(0) << " x "
           << foldingMP(1) << " x "
           << foldingMP(2) << "\n";
      cout << "| Primitive k-points:                 in units of b1,b2,b3\n";
      sxprintf ("|  -ik-     -x-      -y-       -z-    |  -weight-\n");
      for (int ik=0; ik < kVecPrimitive.getSize(); ik++)  {
         k = recCell.carToRel (kVecPrimitive(ik));  // to rel. coordinates
         sxprintf ("| %4d: % 8.6f % 8.6f % 8.6f | %8.6f\n", 
                 ik+1, k(0), k(1), k(2), weightsPrimitive(ik));
      }
      cout << SX_SEPARATOR;

      if (kVecMP.getSize () <= 1000)  {
         cout << "| Unsymmetrized k-points:             in cartesian coordinates\n";
         sxprintf ("|  -ik-     -x-      -y-       -z-    |  -weight-\n");
         for (int ik=0; ik < kVecMP.getSize(); ik++)  {
            k = kVecMP(ik);  // convert to cart. coords
            sxprintf ("| %4d: % 8.6f % 8.6f % 8.6f | %8.6f\n", 
                       ik+1, k(0), k(1), k(2), weightsMP(ik));
         }
         cout << SX_SEPARATOR;
      }
   } 
   cout << "| Symmetrized k-points:               in cartesian coordinates\n";
   // all k-Points were provided by user
   sxprintf ("|  -ik-     -x-      -y-       -z-    |  -weight-    -nG-    -label-\n");
   
   
   for (int ik=0; ik < nk; ik++)  {
      k = kVec(ik);
      if (kLabels.getSize() == nk) 
         l = kLabels(ik);
      else
         l = "";
      sxprintf ("| %4d: % 8.6f % 8.6f % 8.6f | %8.6f    %-6d   %s\n", 
               ik+1, k(0), k(1), k(2), weights(ik),
               gBasisList(ik)->ng, l.ascii());
   }
   cout << SX_SEPARATOR;
   cout << endl;
   SX_MPI_LEVEL("waves-k");
   double ngSum = 0.;
   for (int ik = 0; ik < nk; ++ik)  {
      if (SxLoopMPI::myWork(ik))  {
         ngSum += weights(ik) * gBasisList(ik)->ng;
      }
   }
   ngSum = SxLoopMPI::sum (ngSum);
   double cutVol = ngSum * getTau ().cell.getReciprocalCell ().volume;
   double eCutEff = pow(cutVol * 3. / FOUR_PI, 2./3.);
   sxprintf ("Ecut(eff) = %.12f\n",eCutEff);
}


void SxGkBasis::setNComp(int nComp)
{
   int ik;
   for (ik = 0; ik < nk; ik++)
      gBasisList(ik)->setNComp(nComp);
}

double SxGkBasis::getMinGk () const
{
   double result;
   result = gBasisList(0)->g2(0);
   for (int ik = 1; ik < nk; ik++)   {
      result = result > gBasisList(ik)->g2(0) ? 
         gBasisList(ik)->g2(0) : result; 
   }

   return sqrt(result);
}

double SxGkBasis::getMaxGk () const
{
   double result;
   int dim = (int)gBasisList(0)->g2.getSize ();
   result = gBasisList(0)->g2(dim - 1);
   for (int ik = 1; ik < nk; ik++)   {
      dim = (int)gBasisList(ik)->g2.getSize ();
      result = result < gBasisList(ik)->g2(dim-1) ? 
         gBasisList(ik)->g2(dim-1) : result; 
   }

   return sqrt(result);
}

SxArray<ssize_t> SxGkBasis::getKSortIdx () const
{
   //int nk = getNk();
   SxArray<ssize_t> result (nk);
   SxVector<Double> kNorm(nk);
   for(int ik = 0; ik < nk; ik++) kNorm(ik) = (*this)(ik).getK().norm();
   result = kNorm.getSortIdx();

   return result;
}

void SxGkBasis::registerMemoryObservers ()
{
   TRACK_MEMORY (GkInG);
   TRACK_MEMORY (GplusGinG);
   TRACK_MEMORY (GinGk);
   TRACK_MEMORY (GminusGinG);
   TRACK_MEMORY (GIndx);
   TRACK_MEMORY (GIndp);
   TRACK_MEMORY (nStars);
   TRACK_MEMORY (MGIndx);
   TRACK_MEMORY (MnStars);
   TRACK_MEMORY (MGIndp);
}
