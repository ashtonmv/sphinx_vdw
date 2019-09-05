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

#include <SxGauss.h>
#include <SxBasis.h>
#include <SxProjector.h>

SxGauss::SxGauss ()
   : gkSetPtr (),
     gPtr     (NULL),
     realSpaceMesh (SxGauss::None)
{
   // empty
}

SxGauss::SxGauss (SxPtr<SxGkBasis> gkPtr_,
                  const SxCell &cell_, const SxSymbolTable *table)
   : gkSetPtr (gkPtr_),
     gPtr     (NULL),
     cell     (cell_),
     realSpaceMesh (SxGauss::Radial)
{
   SX_CHECK (table);


   SxSymbolTable         *init = NULL, *gauss = NULL;
   int                    iGauss;
   Coord                  pos;
   double                 width;
   SxString               label;
   bool                   relPos;
   
   meshDim = SxVector3<Int> ((*gkSetPtr)(0).fft3d(0).mesh);

   // --- parse input file
   try  {
      init = table->getGroup("initialization");

      for (gauss  = init->getGroup("gauss"), iGauss = 0;
           gauss != NULL;
           gauss  = gauss->nextSibling ("gauss"), iGauss++)
      {
         pos    =  Coord (gauss->get("coords")->toList());
         width  =  gauss->get("width")->toReal();
         label  = (gauss->contains("label"))
                ?  gauss->get("label")->toString()
                :  "";
         relPos = (gauss->contains("relative"))
                ?  gauss->get("relative")->toAttribute()
                :  false;

         // --- initialization.gauss.relative
         if (relPos)  pos = cell ^ pos;

         positions << pos;
         widths << width;
         labels << label;
      }
      nGauss = iGauss;
   }  catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }

   cout << "nGauss = " << nGauss << endl;

   // --- initialise |R> basis
   R.set (meshDim, cell);

}

void SxGauss::setGBasis (SxGBasis &gBasis)
{
   gPtr = &gBasis;
}

void SxGauss::setRealSpaceMesh (SxGauss::RealSpaceMesh realSpaceMesh_)
{
   realSpaceMesh = realSpaceMesh_;
}

void SxGauss::print () const
{
   if (nGauss == 0)  {
      sxprintf ("|   No initial guesses at the Wannier functions used.\n");
      return;
   }

   sxprintf ("|   -i-   -centers  [in Cartesian coords]-    -width-    "
             "-annotation-\n|\n");
   
   Coord pos;
   SxString x, y, z, w, l;

   for (int iGauss = 0; iGauss < nGauss; iGauss++)  {
      pos = positions(iGauss);
      x   = SxString (pos(0), "% 8.6f");
      y   = SxString (pos(1), "% 8.6f");
      z   = SxString (pos(2), "% 8.6f");
      w   = SxString (widths(iGauss), "%8.5f");
      l   = labels(iGauss);
      sxprintf ("| %4d:    %9s %9s %9s    %9s     %s\n",
                iGauss+1,  x.ascii(), y.ascii(), z.ascii(), w.ascii(),
                l.ascii());
   }
}

void SxGauss::computeGaussiansRadialMesh ()
{
   int                 i, nPoints = 600;
   double              rMin, rMax, ldr, width;
   Coord               pos;
   SxDiracVec<Double>  mesh(nPoints), gauss(nPoints);

   SxArray<SxDiracVec<Double> >  radMeshes(nGauss);
   SxArray<double>               logDr(nGauss);

   gaussiansRadialMesh.resize (nGauss);

   for (int iGauss = 0; iGauss < nGauss; iGauss++)  {
      width  = widths(iGauss);
      pos    = positions(iGauss);

      // --- create radial basis
      rMin = 0.0001 / width;   // the more localised, the closer to 0
      rMax = 1.07238 * width;  // Gaussian has value exp(-4.6) = 0.01 there,
                               // and I think, the rest can be neglected for
                               // most applications.
      rMax *= 2;
      ldr  = log(rMax/rMin) / double(nPoints - 1);
      for (i = 0; i < nPoints; i++)  mesh(i) = rMin * exp(ldr * i);
      SX_CHECK (fabs(mesh(nPoints-1)-rMax) < 1.e-12, mesh(nPoints-1), rMax);
      radMeshes(iGauss).resize (nPoints);
      radMeshes(iGauss).set (mesh);
      logDr(iGauss) = ldr;

      // --- create Gaussian centered aroud (0,0,0)
      gauss  = exp ( -((2. * mesh / width).sqr()) );
      gauss *= 64. / (sqrt (PI * PI * PI) * width * width * width);
      // TODO: check the norm in G space, is it necessary at all?

      // TODO: ugly!!!
      gauss.handle->auxData.is = iGauss;  // iGauss plays the role of iSpecies
      gauss.handle->auxData.l  = 0;       // angular momentum not needed
      gauss.handle->auxData.m  = 0;       // magnetic quantum number not needed

      gaussiansRadialMesh(iGauss).resize (nPoints);
      gaussiansRadialMesh(iGauss).set (gauss);
   }

   // --- initialise radial basis |r> and cell periodic |G+0> basis
   rPtr = SxPtr<SxRadBasis>::create (radMeshes, logDr);

   // --- set |r> basis to the Gaussians
   for (int iGauss = 0; iGauss < nGauss; iGauss++)
      gaussiansRadialMesh(iGauss).setBasis (&*rPtr);

   // --- compute structure factors required for shifting the Gaussians
   computeStructureFactorsG ();
}

void SxGauss::computeGaussiansFFTMesh ()
{
   int             iMesh, nMesh = (*gkSetPtr)(0).fft3d(0).meshSize;
   int             m1, m2, m3;
   Coord           s, r, pos, m, a;
   double          t, width;
   SxArray<Coord>  mesh = computeMeshPoints (meshDim, RelVec (1,1,1), cell);
   cout << "nMeshPoints = " << nMesh << endl;

   gaussiansFFTMesh.resize (nGauss);

   cout << "| Creating Gauss bells ... ";  cout.flush ();

   for (int iGauss = 0; iGauss < nGauss; iGauss++)  {
      if (iGauss+1 < nGauss)  cout << iGauss+1 << ", ";
      else                    cout << nGauss   << " and ";
      cout.flush();
      gaussiansFFTMesh(iGauss).resize (nMesh);
      pos   = positions(iGauss);
      width = widths(iGauss);

      // --- consider - at least approximately - periodic boundaries
      gaussiansFFTMesh(iGauss).set (0.);
      for (m1 = -2; m1 <= 2; m1++)  {
         for (m2 = -2; m2 <= 2; m2++)  {
            for (m3 = -2; m3 <= 2; m3++)  {
               m = SxVector3<Double> ((double)m1, (double)m2, (double)m3);
               a = cell ^ m;

               for (iMesh = 0; iMesh < nMesh; iMesh++)  {
                  r = mesh(iMesh) - pos - a;
                  s = 2. * r / width;
                  t = - s.normSqr ();
                  gaussiansFFTMesh(iGauss)(iMesh) += exp( t );
                  /*
                  if (iGauss==1 && iMesh==1)  {
                  cout << "mesh(iMesh) = " << mesh(iMesh) << endl;
                  cout << "pos = " << pos << endl;
                  cout << "a = " << a << endl;
                  cout << "r = " << r << endl;
                  cout << "s = " << s << endl;
                  cout << "t = " << t << endl;
                  cout << "gaussiansFFTMesh(iGauss)(iMesh) = "
                       << gaussiansFFTMesh(iGauss)(iMesh) << endl;
                  }
                  */
               }  // :iMesh
            }  // :m3
         }  // :m2
      }  // :m1
   }  // :iGauss

   cout << "done\n";

   // --- set |R> basis to the Gaussians
   for (int iGauss = 0; iGauss < nGauss; iGauss++)
      gaussiansFFTMesh(iGauss).setBasis(&R);
}

void SxGauss::computeGaussiansOld ()
{
   cout << SX_SEPARATOR;
   cout << "| Old initialisation" << endl;
   cout << "|       Sometimes S/PHI/nX is not at all willing to do the" << endl;
   cout << "|       things you want it to do. Then there is no way other"<<endl;
   cout << "|       than a step by step comparison with a working code!"<< endl;
   cout << SX_SEPARATOR;
   cout << endl;

   gaussians.resize (nGauss);

   SxVector3<Double> a1 = SxVector3<Double> (cell.row(0));
   SxVector3<Double> a2 = SxVector3<Double> (cell.row(1));
   SxVector3<Double> a3 = SxVector3<Double> (cell.row(2));

   cout << SX_SEPARATOR;
   cout << "| unit cell vectors" << endl;
   cout << SX_SEPARATOR;
   sxprintf ("|   a1 = (%5.2g, %5.2g, %5.2g)\n", a1(0), a1(1), a1(2));
   sxprintf ("|   a2 = (%5.2g, %5.2g, %5.2g)\n", a2(0), a2(1), a2(2));
   sxprintf ("|   a3 = (%5.2g, %5.2g, %5.2g)\n", a3(0), a3(1), a3(2));

   const SxCell::CellType cellType (cell);
   const double aLat = cellType.getA ();
   cellType.print ();

   const double sigma = sqrt (3.) * aLat / 32.;
                                    // sigma: spread of Gaussian bell-curves

   // --- centers of Gaussian bell-curves
   cout << endl;
   cout << SX_SEPARATOR;
   const double s = aLat * 0.125;  // TODO: generalisation needed (fcc only)!
   const SxVector3<Double> transl = 0.25 * (a1 + a2 + a3);
   SxArray<Coord>          tau(nGauss);  // centers of Gaussian bell-curves
   sxprintf ("| preparing \"bonding\" initialisation ...\n");
   tau(0) = transl - SxVector3<Double> ( s, s, s);  // bonding
   tau(1) = transl - SxVector3<Double> ( s,-s,-s);  // bonding
   tau(2) = transl - SxVector3<Double> (-s, s,-s);  // bonding
   tau(3) = transl - SxVector3<Double> (-s,-s, s);  // bonding
   if (nGauss == 5)  {
      tau(4) = SxVector3<Double> (0,0,0);  // Ga
   }

   // --- write on screen
   cout << "| initialise Wannier functions "
        << "by projecting \"" << "..." << "\" onto Gaussians" << endl;
   cout << SX_SEPARATOR;
   sxprintf ("| centers of Gaussians:\n");
   for (int iGauss = 0; iGauss < nGauss; iGauss++)  {
      sxprintf ("|   tau( %d ) = ", iGauss);
      sxprintf ("(%5.2g,%5.2g,%5.2g)\n",
                tau(iGauss)(0), tau(iGauss)(1), tau(iGauss)(0));
   }

   // --- create Gaussians
   int                 ir, nr = (*gkSetPtr)(0).fft3d(0).meshSize;
   int                 m1, m2, m3;
   double              enumerator = 1. / (2. * sigma * sigma);
   double              sumCells, exponent;
   SxVector3<Double>   translation;  // translation into neighbour cell
   SxVector3<Double>   rSpot;
   SxVector3<Int>      reptn = SxVector3<Int> (1,1,1);
   SxArray<Coord>      bigMesh = computeMeshPoints (meshDim, reptn, cell);

   for (int iGauss = 0; iGauss < nGauss; iGauss++)  {
      gaussians(iGauss).resize (nr);

      for (ir = 0; ir < nr; ir++)  {
         // consider -- at least approximately --
         // periodic boundary conditions
         sumCells = 0.;
         for (m1 = -2; m1 <= 2; m1++)  {
            for (m2 = -2; m2 <= 2; m2++)  {
               for (m3 = -2; m3 <= 2; m3++)  {
                  translation  = (double)m1*a1 + (double)m2*a2 + (double)m3*a3;
                  rSpot        = bigMesh(ir) - (tau(iGauss) + translation);
                  exponent     = (rSpot * rSpot).sum() * enumerator;
                  sumCells    += exp (-exponent);  // normalisation comes later
               }
            }
         }

         gaussians(iGauss)(ir) = sumCells;
      }  // :ir
   }  // :iGauss
   cout << SX_SEPARATOR;
   cout << endl;

   // --- normalise Gaussians
   double               normQuad;  // norm square
   SxDiracVec<TPrecPhi> f;         // f considered as C^n vector, not as fctn.

   for (int iGauss = 0; iGauss < nGauss; iGauss++)  {
      f = SxDiracVec<TPrecPhi> (gaussians(iGauss));
      normQuad = (f * f).sum();  // C^n norm
      gaussians(iGauss) /= sqrt (normQuad);
      cout << "norm = " << sqrt (normQuad) << endl;
   }

   // --- write Gaussians to output file(s)
   SxString filename = "gOld";
   try  {
      SxBinIO io;
      for (int iGauss = 0; iGauss < nGauss; iGauss++)  {
         io.open      (filename+"-"+(iGauss+1)+".sxb", SxBinIO::BINARY_WRITE_ONLY);
         io.writeMesh (gaussians(iGauss), cell, meshDim);
         io.setMode   (SxBinIO::WRITE_DATA);
         io.writeMesh (gaussians(iGauss), cell, meshDim);
         io.close ();
      }
   }  catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }
}

void SxGauss::computeGaussiansGSpace ()
{
   SX_CHECK (gPtr);
   SxGBasis &G = *gPtr;

   gaussiansGSpace.resize (nGauss); 

   if (realSpaceMesh == Radial)  {
      if (gaussiansRadialMesh.getSize() == 0)  computeGaussiansRadialMesh ();
      SX_CHECK (gaussiansRadialMesh.getSize() == nGauss,
                gaussiansRadialMesh.getSize(), nGauss);
      SX_CHECK      (&*rPtr);
      //SxRadBasis  &r = *rPtr;

      // --- enable projections
      //G.registerRBasis (R);
      //R.registerGBasis (G);

      // --- create periodic Gaussians in G space
      //     -
      //     Please note, that the following projections generate periodic
      //     functions out of the non-periodic Gauss functions given on the
      //     radial grid, since we perform actually a Fouries SERIES expansion
      //     (the G's are discrete) rather than a Fourier transformation.
      for (int iGauss = 0; iGauss < nGauss; iGauss++)
         gaussiansGSpace(iGauss) = ( G | gaussiansRadialMesh(iGauss) );

      // --- shift Gaussians to the specified positions
      for (int iGauss = 0; iGauss < nGauss; iGauss++)
         gaussiansGSpace(iGauss) *= structureFactorsG(iGauss);

      return;
   }

   if (realSpaceMesh == FFT)  {
      if (gaussiansFFTMesh.getSize() == 0)  computeGaussiansFFTMesh ();
      SX_CHECK (gaussiansFFTMesh.getSize() == nGauss,
                gaussiansFFTMesh.getSize(), nGauss);
      SX_CHECK      (&R);

      // --- enable projections
      //G.registerRBasis (R);
      //R.registerGBasis (G);

      // --- bring Gaussians to G space
      for (int iGauss = 0; iGauss < nGauss; iGauss++)
         gaussiansGSpace(iGauss) = ( G | gaussiansFFTMesh(iGauss) );

      return;
   }

   // --- should not happen
   SX_EXIT;
}

void SxGauss::computeStructureFactorsG ()
{
   SX_CHECK (gPtr);
   SxGBasis                &G = *gPtr;
   int                      ig, ng = G.ng;
   double                   phase;
   SxDiracVec<TPrecCoeffG>::Iterator strFacIt;
   Coord                    pos;

   structureFactorsG.resize(nGauss);

   for (int iGauss = 0; iGauss < nGauss; iGauss++)  {
      structureFactorsG(iGauss).resize(ng);
      pos = positions(iGauss);

      strFacIt = structureFactorsG(iGauss).begin();
      for (ig = 0; ig < ng; ig++, strFacIt++)  {
         phase     = -(G.getG(ig) ^ pos);
         *strFacIt = PrecCoeffG (cos(phase), sin(phase));
      }
   }
}

void SxGauss::computeGaussiansGkSpace ()
{
   SX_CHECK (gkSetPtr);
   SX_CHECK (&*rPtr);
   //SxRadBasis &r = *rPtr;
   SxGkBasis  &gkSet = *gkSetPtr;
   SxGBasis   *gkPtr;
   int         ik, nk = gkSet.nk;
   int         iSpin, nSpin = 1;  // spin polarisation not implemented
   int         iGauss;
   SxDiracVec<TPrecCoeffG> strFacs;

   SX_CHECK (realSpaceMesh == Radial, realSpaceMesh, Radial);

   if (gaussiansRadialMesh.getSize() == 0)  computeGaussiansRadialMesh ();
   SX_CHECK (gaussiansRadialMesh.getSize() == nGauss,
             gaussiansRadialMesh.getSize(), nGauss);

   // --- Gaussians projected on the |G+k> basises stored like Bloch functions
   gaussiansGkSpace = SxPW (nGauss, nSpin, gkSetPtr);

   // --- create periodic Gaussians in G space
   //     -
   //     Please note, that the following projections generate periodic
   //     functions out of the non-periodic Gauss functions given on the
   //     radial grid, since we perform actually a Fouries SERIES expansion
   //     (the G's are discrete) rather than a Fourier transformation.
   for (ik = 0; ik < nk; ik++)  {
      gkPtr = &gkSet(ik);

      for (iSpin = 0; iSpin < nSpin; iSpin++)  {
         for (iGauss = 0; iGauss < nGauss; iGauss++)  {

            // --- < G+k | projection
            gaussiansGkSpace(iGauss, iSpin, ik)
               <<= ( *gkPtr | gaussiansRadialMesh(iGauss) );

            // --- shift to correct position by structure factors
            strFacs = computeStructureFactorsGk(ik, iGauss);
            gaussiansGkSpace(iGauss, iSpin, ik) *= strFacs;
         }  // :iGauss
      }  // :iSpin
   }  // :ik
}

SxDiracVec<TPrecCoeffG>
SxGauss::computeStructureFactorsGk (int ik, int iGauss)
{
   SX_CHECK (gkSetPtr);
   SxGBasis                &gk = (*gkSetPtr)(ik);
   int                      ig, ng = gk.ng;
   double                   phase;
   SxDiracVec<TPrecCoeffG>  strFacs(ng);
   SxDiracVec<TPrecCoeffG>::Iterator strFacIt;
   Coord                    pos = positions(iGauss);

   strFacIt = strFacs.begin();
   for (ig = 0; ig < ng; ig++, strFacIt++)  {
      phase     = -(gk.getG(ig) ^ pos);
      *strFacIt = PrecCoeffG (cos(phase), sin(phase));
   }

   return strFacs;
}
   
void SxGauss::write ()
{
   // --- radial mesh case
   if (realSpaceMesh == Radial)  {
      SX_CHECK (gaussiansRadialMesh.getSize() > 0,
                gaussiansRadialMesh.getSize());

      // --- Gaussians given on a radial mesh can't be plotted directly,
      //     but have to be brought to the G space first (then they are
      //     periodic and at the right position) and then transformed to
      //     their phinaxable FFT mesh representation.
      if (gaussiansGSpace.getSize() == 0)  computeGaussiansGSpace ();
      SxArray<SxDiracVec<TPrecPhi> >  gaussiansRSpace(nGauss);
      for (int iGauss = 0; iGauss < nGauss; iGauss++)
         gaussiansRSpace(iGauss) = ( R | gaussiansGSpace(iGauss) );

      // --- write to output file(s)
      SxString filename = "gauss";
      for (int iGauss = 0; iGauss < nGauss; iGauss++)  {
         try  {
            SxBinIO io (filename+"-"+(iGauss+1)+".sxb", SxBinIO::BINARY_WRITE_ONLY);
            io.writeMesh (gaussiansRSpace(iGauss), cell, meshDim);
            io.setMode (SxBinIO::WRITE_DATA);
            io.writeMesh (gaussiansRSpace(iGauss), cell, meshDim);
            io.close ();
         }  catch (SxException e)  {
            e.print ();
            SX_QUIT;
         }
      }

      return;
   }

   // --- FFT mesh case
   if (realSpaceMesh == FFT)  {
      SX_CHECK (gaussiansFFTMesh.getSize() > 0, gaussiansFFTMesh.getSize());
      SxString filename = "g";
      for (int iGauss = 0; iGauss < nGauss; iGauss++)  {
         try  {
            SxBinIO io (filename+"-"+(iGauss+1)+".sxb", SxBinIO::BINARY_WRITE_ONLY);
            io.writeMesh (gaussiansFFTMesh(iGauss), cell, meshDim);
            io.setMode (SxBinIO::WRITE_DATA);
            io.writeMesh (gaussiansFFTMesh(iGauss), cell, meshDim);
            io.close ();
         }  catch (SxException e)  {
            e.print ();
            SX_QUIT;
         }
      }
      return;
   }

   // --- should not happen
   SX_EXIT;
}

void SxGauss::project (const SxPW &waves, int iBottom, int iTop)
{
   SX_EXIT;  // to be accomplished
   SX_CHECK (gkSetPtr);
   int nSpin   = waves.getNSpin();
   int nStates = iTop - (iBottom - 1);
   SX_CHECK  (nStates > 0, nStates);
//   SxPW projections;

   // --- could happen in non-debug cases
   if (nStates != nGauss)  {
      sxprintf ("ERROR:  The number of Gaussians provided must be identical\n"
                "        to the number of bands you want to consider, i.e.\n"
                "\n"
                "        iTop - (iBottom - 1) .\n");
      SX_QUIT;
   }

   projections = SxPW (nStates, nSpin, gkSetPtr);

   // ...
}

void SxGauss::projectOld (const SxPW &waves, int iBottom, int iTop)
{
   SX_CHECK       (gkSetPtr);
   SX_CHECK   (realSpaceMesh == FFT, realSpaceMesh);
   int          ik, nk = waves.getNk();
   int          iSpin, nSpin = waves.getNSpin();
   int          iGauss;
   int          i, nStates = iTop - (iBottom - 1);
   SX_CHECK  (nStates > 0 && nStates <= waves.getNStates(),
              nStates, waves.getNStates());
   SxGkBasis   &gkSet = *gkSetPtr;
   SX_CHECK  (gaussiansFFTMesh.getSize() == nGauss,
              gaussiansFFTMesh.getSize(), nGauss);
   //SxDirTensor3<Complex16> D1;

   PrecCoeffG              spr;
   SxDiracVec<TPrecCoeffG> gaussGk, u0i;

   // --- could happen in non-debug cases
   if (nStates != nGauss)  {
      sxprintf ("ERROR:  The number of Gaussians provided must be identical\n"
                "        to the number of bands you want to consider, i.e.\n"
                "\n"
                "        iTop - (iBottom - 1) .\n");
      SX_QUIT;
   }

   // --- resize projection trafo
   trafoProjection.resize(nk);
   for (ik = 0; ik < nk; ik++)  {
      trafoProjection(ik).resize (nSpin);
      for (iSpin = 0; iSpin < nSpin; iSpin++)  {
         trafoProjection(ik)(iSpin).reformat (nStates,nStates);
      }
   }

   // --- calculate projections
   projections = SxPW (nStates, nSpin, gkSetPtr);

   for (ik = 0; ik < nk; ik++)  {

      for (iSpin = 0; iSpin < nSpin; iSpin++)  {
         for (iGauss = 0; iGauss < nGauss; iGauss++)  {
            projections(iGauss, iSpin, ik).set (0.);
            gaussGk = ( gkSet(ik) | gaussiansFFTMesh(iGauss) );

            for (i = iBottom; i <= iTop; i++)  {
               u0i = waves(i, iSpin, ik);
               spr = (u0i ^ gaussGk).chop();
               projections(iGauss, iSpin, ik) += u0i * spr;

               trafoProjection(ik)(iSpin)(i,iGauss) = spr;
               //D1.set (spr, ik, i, iGauss);
            }  // :i
         }  // :iGauss
      }  // :iSpin
   }  // :ik

   //D1.writeDirTensor3 ("D1.sxb");
}

void SxGauss::orthonormalizeProjections (SxOrthoMethod method)
{
   if (method == Loewdin)  {
      int ik, nk = projections.getNk();
      int iSpin, nSpin = projections.getNSpin();
      int nStates = projections.getNStates();
      if (nSpin > 1)  SX_QUIT;  // spin polarisation not yet implemented

      // --- resize array for Loewdin trafo
      trafoLoewdin.resize(nk);
      for (ik = 0; ik < nk; ik++)  {
         trafoLoewdin(ik).resize (nSpin);
         for (iSpin = 0; iSpin < nSpin; iSpin++)  {
            trafoLoewdin(ik)(iSpin).reformat (nStates,nStates);
         }
      }

      // --- orthonormalisation
      projections.orthonormalize (method, &trafoLoewdin);
   }  else  {
      projections.orthonormalize (method);
   }
}

SxArray<SxArray<SxDiracMat<TPrecCoeffG> > >
SxGauss::getProjectionTrafo () const
{
   SX_CHECK (trafoProjection.getSize() > 0, trafoProjection.getSize());

   return trafoProjection;
}

SxArray<SxArray<SxDiracMat<TPrecCoeffG> > >
SxGauss::getLoewdinTrafo () const
{
   SX_CHECK (trafoLoewdin.getSize() > 0, trafoLoewdin.getSize());

   return trafoLoewdin;
}

SxArray<Coord>
SxGauss::computeMeshPoints (const RelVec &dim, const RelVec &repeat,
                            const SxCell &cell_)
{
   SX_CHECK (cell_.volume > 0., cell_.volume);
   SX_CHECK (dim(0) > 0, dim(0));
   SX_CHECK (dim(1) > 0, dim(1));
   SX_CHECK (dim(2) > 0, dim(2));
   SX_CHECK (repeat(0) > 0, repeat(0));
   SX_CHECK (repeat(1) > 0, repeat(1));
   SX_CHECK (repeat(2) > 0, repeat(2));

   SxMesh3D        dimOut = dim * repeat;
   Coord           r, rRel, trans;
   SxArray<Coord>  meshOut( dimOut.product() );

   double delta1 = 1./(double)dim(0);
   double delta2 = 1./(double)dim(1);
   double delta3 = 1./(double)dim(2);

   int x, y, z, i, j, k, idx, xOff, yOff, zOff;

   // --- access points of the mesh
   for (z = 0; z < dim(2); z++)  {
      for (y = 0; y < dim(1); y++)  {
         for (x = 0; x < dim(0); x++)  {

            rRel = Coord ((double)x * delta1,
                          (double)y * delta2,
                          (double)z * delta3);

            r    = cell_ ^ rRel;

            // --- translation to neighboured unit cells
            for (k = 0; k < repeat(2); k++)  {
               zOff = z + k * dim(2);
               for (j = 0; j < repeat(1); j++)  {
                  yOff = y + j * dim(1);
                  for (i = 0; i < repeat(0); i++)  {
                     xOff = x + i * dim(0);

                     rRel  = Coord ((double)i, (double)j, (double)k);
                     trans = cell_ ^ rRel;
                     idx   = (int)dimOut.getMeshIdx (xOff, yOff, zOff,
                                                     SxMesh3D::Positive);
                     meshOut(idx) = r + trans;

                  }
               }
            }

         }
      }
   }

   return meshOut;
}

