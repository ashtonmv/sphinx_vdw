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

#include <SxTimer.h>
#include <SxEAM.h>
#include <fstream>
SxEAM::SxEAM (const SxAtomicStructure &str, 
                          const SxSymbolTable *table) {
   sxprintf ("This is SxEAM::SxEAM\n");
   speciesData = SxSpeciesData (&*table);
   nAtoms = str.getNAtoms ();
   species.resize (nAtoms);
   nSpecies = str.getNSpecies ();
   cout << nSpecies << endl;
   

   //--- up to 2 different species implemented yet
   if (nSpecies > 2) {
      cout << "Only up to 2 different species implemented in SxEAM" << endl;
      SX_QUIT;
   } else {
      if (nSpecies == 1) {
         V.resize (1);  F.resize (1); rho.resize (1);
         dV.resize (1);  dF.resize (1); drho.resize (1);
      } 
      if (nSpecies == 2) {
         V.resize (3);  F.resize (3); rho.resize (3);
         dV.resize (3);  dF.resize (3); drho.resize (3);
      } 
   }
   
   int counter = 0;
   for (int is = 0; is < str.getNSpecies (); is++) {
      SxString elem = speciesData.chemName(is);
      for (int ia = 0; ia < str.getNAtoms(is); ia++)  {
         species(counter) = elem;
         counter++;
      }
   }

   
   setupSplines ((table -> getGroup("eamPot") -> getGroup("params")));

   if (table -> getGroup("eamPot") -> contains("noLinkCell"))
      noLinkCell = table -> getGroup("eamPot") -> get ("noLinkCell") -> toBool ();
   else
      noLinkCell = false;   
  
   neigh.resize (nAtoms);
   dist.resize (nAtoms);
   pw.resize (nAtoms);
   cell.resize (nAtoms);

   if (!noLinkCell) {
      length = 6;
      setLinkCellMethodParams (str);
   }
   update (str);

   sxprintf("TOTAL ENERGY (Hartree): %.11f\n", getTotalEnergy ());

  /* 
   SxAtomicStructure test; test.copy (tau); 

   SxMatrix3<Double> mat; mat.set(0.);
   
   for (double aLat = 3.45;  aLat < 8.0; aLat += 0.01) {
      
      double a0 = tau.cell(0)(0);
      double a = 2.*A2B*aLat;
      double scale = a/a0;
      mat (0, 0) = mat(1, 1) = mat(2, 2) = a;
      test.cell.set(mat);
      test.set (scale*test.coordRef ());
      update (test);

      cout << "MURN: " << aLat << " " << getTotalEnergy () << endl;


   }
*/
  /*  
   SxAtomicStructure test; test.copy (tau);

   for (int i = 0; i < nAtoms; i++) {
      for (int j = 0; j < 3; j++) {
         test.ref(i)(j) += 0.05*rand()/(double) 0x7fffffff;
      }
   }
   

     
   SxAtomicStructure f1 = getForces (test, table);
   SxAtomicStructure f2 = getNumericalForces (test, 0.01);
   SxAtomicStructure f3 = getNumericalForces (test, 0.0001);

   SxAtomicStructure f = f1 - f2;
  cout << "LARS3" << endl; fflush(stdout); 

   cout << f1 << endl;
   cout << f2 << endl;
   cout << f3 << endl; 
   cout << f << endl; 
   QUIT;
   */
  
   /*
   for (int i = 0; i < 100; i++) {
      cout <<  "UPDATE" << i << endl;
      SxAtomicStructure f = getForces (test, table);
      
   }
   
   for (int i = 0; i < neigh.getSize (); i++) 
      cout << neigh(i).getSize () << endl;

      */
}

SxEAM::~SxEAM ()
{
   // empty
}


void SxEAM::setLinkCellMethodParams (const SxAtomicStructure &tau)
{
   int i, j;
   double min =  1e20;
   double max = -1e20;
   for ( i = 0; i < nAtoms; i++ )  {
      for (j = 0; j < 3; j++) {
         if ( tau(i)(j) < min ) min = tau(i)(j);
         if ( tau(i)(j) > max ) max = tau(i)(j);
      }
   }
   cout << "min " << min << endl;
   cout << "max " << max << endl;
   meshOrigin.set(min-(cutoff+2*length));
   double meshLength = max + 2*cutoff + 4*length - meshOrigin(0);
   n = int(meshLength/length);
   length = meshLength/n;


   cout << "n: " << n << endl;
   cout << "length: " << length << endl;
   cout << "meshOrigin: " << meshOrigin << endl;

   int nHalf = (int) n/2;
   SxVector3<Double> RmeshHalf;
   RmeshHalf.set(length*(nHalf-1));
   RmeshHalf += meshOrigin;

   SxArray<SxVector3<Int> > a1, a2;
   a1.resize(64);
   a2.resize(64);
   for (i = 0; i < 64; i++) {
		a1(i)(0) =  (int) i/32;
		a1(i)(1) =  (int) (i - a1(i)(0)*32)/16;
		a1(i)(2) =  (int) (i - a1(i)(0)*32 - a1(i)(1)*16)/8;
		a2(i)(0) =  (int) (i - a1(i)(0)*32 - a1(i)(1)*16 - a1(i)(2)*8)/4;
		a2(i)(1) =  (int) (i - a1(i)(0)*32 - a1(i)(1)*16 - a1(i)(2)*8 
                        - a2(i)(0)*4)/2;
		a2(i)(2) =  (int) (i - a1(i)(0)*32 - a1(i)(1)*16 - a1(i)(2)*8 
                        - a2(i)(0)*4 - a2(i)(1)*2);
   }

   SxList<SxVector3<Int> > neighborsReference;
   neighborsReference.resize(0);
   SxVector3<Double> Rmesh, v;
   int x, y, z;
   double norm;
   double minNorm;
   for (z = 0; z < n; z++) {
      for (y = 0; y < n; y++) {
         for (x = 0; x < n; x++) {
            minNorm = 1e20;
            Rmesh(0) = length*x + meshOrigin(0);
            Rmesh(1) = length*y + meshOrigin(1);
            Rmesh(2) = length*z + meshOrigin(2);
            for (i = 0; i < 64; i++) {
               v = Rmesh + length*a1(i) - (RmeshHalf + length*a2(i));
               if ((norm = v.norm()) < minNorm) minNorm = norm;
            }
            if ((minNorm < cutoff) && ((nHalf-1)!=x || (nHalf-1)!=y 
                     || (nHalf-1)!=z))
               neighborsReference.append(SxVector3<Int> (x,y,z));
         }
      }
   }
   cout << "size of neighborsReference " << 
      neighborsReference.getSize() << endl;

   SxList<SxVector3<Int> > origMeshTmp;
   SxList<SxArray<SxVector3<Int> > > neighborsTmp;
   SxArray<SxVector3<Int> > shiftedReference = neighborsReference;
   SxVector3<Double> midVec;
   origMeshTmp.resize(0);
   neighborsTmp.resize(0);
   for (z = 0; z < n; z++) {
      for (y = 0; y < n; y++) {
         for (x = 0; x < n; x++) {
            Rmesh = length * (SxVector3<Int> (x,y,z)) + meshOrigin;
            for (i = 0; i < 64; i=i+8) {
					v = Rmesh + length*a1(i);
               v = tau.cell.inverse()^v;
					if ((v(0) >= 0. && v(0) < 1.) && (v(1) >= 0. && v(1) < 1.)
                     && v(2) >= 0. && v(2) < 1.) {
                  origMeshTmp.append(SxVector3<Int> (x,y,z));
                  for (j = 0; j < shiftedReference.getSize(); j++)
                     shiftedReference(j) = neighborsReference(j) + 
                        (SxVector3<Double> (x-nHalf+1,y-nHalf+1,z-nHalf+1));
                  neighborsTmp.append(shiftedReference);
                  break;
               }
            }
         }
      }
   }
   origMesh = origMeshTmp;
   neighbors = neighborsTmp;
}


bool SxEAM::isRegistered (const SxSymbolTable *cmd) const
{
   SX_CHECK (cmd);
   SxString str = cmd->getName ();
   return ( str=="forces" || str=="HesseMatrix" || str=="TTensor" );
}


int SxEAM::getInteractionType (int i, int j) 
{
   if (speciesData.chemName.getSize () == 1) {
      return 0;
   } else {
      if (i == j) {
         if  (species(i) == speciesData.chemName(0))  return 0;
         if  (species(i) == speciesData.chemName(1))  return 1;
         SX_EXIT;
      } else {
         if ( (species(i) == speciesData.chemName(0)) &&
               (species(j) == speciesData.chemName(0)) ) return 0;
         if ( (species(i) == speciesData.chemName(1)) &&
               (species(j) == speciesData.chemName(1)) ) return 1;
         return 2;
      }
   }
}
           



double SxEAM::getDist (int i, int j, const SxVector3<Double> &R) 
{
   SxVector3<Double> v = tau.ref(i) - tau.ref(j) - R;

   if (fabs(v(0)) > cutoff) return fabs (v(0));
   if (fabs(v(1)) > cutoff) return fabs (v(1));
   if (fabs(v(2)) > cutoff) return fabs (v(2));
	
   return v.norm ();

}

void SxEAM::updateNeighsLinkCellMethod () 
{
   int cX, cY, cZ, i, j;
  // int sCut = (int) (cutoff/tau.cell(0, 0)) + 1;
 
   int sCut = 1;
   for (i = 0; i < 3; i++) {
      int a = (int) (cutoff/tau.cell(i, i)) + 1;
      if (a > sCut) sCut = a;
   }
	



   SxVector3<Double> v, v1, vMesh, R, vCheck, vReduced;

   SxArray<SxList<int> > neighList;
   SxArray<SxList<int> > pwList;
   SxArray<SxList<double> > distList;
   SxArray<SxList<SxVector3<Double> > > cellList;

   neighList.resize(nAtoms);
   distList.resize(nAtoms);
   pwList.resize(nAtoms);
   cellList.resize(nAtoms);

   SxArray<SxArray<SxArray<SxList<int> > > > meshAtomNr, meshAtomNrPeriodic;
   SxArray<SxArray<SxArray<SxList<SxVector3<Double> > > > > 
                              meshAtomCoord, meshAtomCoordPeriodic;
   SxArray<SxArray<SxArray<SxList<SxVector3<Double> > > > > 
                              meshAtomCell, meshAtomCellPeriodic;
   meshAtomNr.resize(n);
   meshAtomNrPeriodic.resize(n);
   meshAtomCoord.resize(n);
   meshAtomCoordPeriodic.resize(n);
   meshAtomCell.resize(n);
   meshAtomCellPeriodic.resize(n);
   for (i = 0; i < n; i++) {
      meshAtomNr(i).resize(n);
      meshAtomNrPeriodic(i).resize(n);
      meshAtomCoord(i).resize(n);
      meshAtomCoordPeriodic(i).resize(n);
      meshAtomCell(i).resize(n);
      meshAtomCellPeriodic(i).resize(n);
		for (j = 0; j < n; j++) {
         meshAtomNr(i)(j).resize(n);
         meshAtomNrPeriodic(i)(j).resize(n);
         meshAtomCoord(i)(j).resize(n);
         meshAtomCoordPeriodic(i)(j).resize(n);
         meshAtomCell(i)(j).resize(n);
         meshAtomCellPeriodic(i)(j).resize(n);
      }
   }
   
   SxCell inv = tau.cell.inverse();

  
    for (cZ = -sCut; cZ <= sCut; cZ++) {
      for (cY = -sCut; cY <= sCut; cY++) {
         for (cX = -sCut; cX <= sCut; cX++) { 
            R = cX*tau.cell(0) + cY*tau.cell(1) + cZ*tau.cell(2);
            for (i = 0; i < nAtoms; i++) {
               v = R + tau(i);
               vMesh = (v - meshOrigin)/length;
               vReduced = inv^v;

              // cout << vReduced << " " << n << endl;

               if (0 <= vMesh(0) && vMesh(0) < n 
                     && 0 <= vMesh(1) && vMesh(1) < n
                     && 0 <= vMesh(2) && vMesh(2) < n) {
                  if ((vReduced(0) >= -0.00001 && vReduced(0) < 0.99999999) 
                        && (vReduced(1) >= -0.00001 && vReduced(1) < 0.99999999)
                        && (vReduced(2) >= -0.00001 && vReduced(2) < 0.99999999)) {

         // <--------
         meshAtomNr((int)vMesh(0))((int)vMesh(1))((int)vMesh(2)).append(i);
         meshAtomCoord((int)vMesh(0))((int)vMesh(1))((int)vMesh(2)).append(v);
         meshAtomCell((int)vMesh(0))((int)vMesh(1))((int)vMesh(2)).append(R);
         // <--------
                  } else {
   // <--------------
   meshAtomNrPeriodic((int)vMesh(0))((int)vMesh(1))((int)vMesh(2)).append(i);
   meshAtomCoordPeriodic((int)vMesh(0))((int)vMesh(1))((int)vMesh(2)).append(v);
   meshAtomCellPeriodic((int)vMesh(0))((int)vMesh(1))((int)vMesh(2)).append(R);
   // <--------------
                  }
               }
            }
         }
      }
   }
   
   for (i = 0; i < nAtoms; i++) {
      distList(i). resize (0);
      pwList(i). resize (0);
      neighList(i).resize (0);
      cellList(i). resize (0);
   }
   
   int a1, a2, k, l;
   double distance;
   SxArray<int> currentAtomNrs, neighborsAtomNrs;
   SxArray<int> currentAtomNrsPeriodic, neighborsAtomNrsPeriodic;
   SxArray<SxVector3<Double> > currentAtomCoords, neighborsAtomCoords;
   SxArray<SxVector3<Double> > currentAtomCoordsPeriodic, neighborsAtomCoordsPeriodic;
   SxArray<SxVector3<Double> > currentAtomCells, neighborsAtomCells;
   SxArray<SxVector3<Double> > currentAtomCellsPeriodic, neighborsAtomCellsPeriodic;


   for (i = 0; i < origMesh.getSize(); i++) {
      int a = origMesh(i)(0);
      int b = origMesh(i)(1);
      int c = origMesh(i)(2);

      currentAtomNrs = meshAtomNr(a)(b)(c);
      currentAtomCoords = meshAtomCoord(a)(b)(c);
      currentAtomCells =  meshAtomCell(a)(b)(c);
      currentAtomNrsPeriodic =  meshAtomNrPeriodic(a)(b)(c);
      currentAtomCoordsPeriodic =   meshAtomCoordPeriodic(a)(b)(c);
      currentAtomCellsPeriodic =   meshAtomCellPeriodic(a)(b)(c);
		
      for (j = 0; j < currentAtomNrs.getSize(); j++) {
         a1 = currentAtomNrs(j);
         v1 = currentAtomCoords(j);
         for (k = 0; k < currentAtomNrs.getSize(); k++) {
            if ( k != j ) {
               a2 = currentAtomNrs(k);
               v = currentAtomCoords(k);
               R = currentAtomCells(k);
               vCheck = v1-v;
               if (fabs(vCheck(0)) < cutoff && fabs(vCheck(1)) < cutoff
                     && fabs(vCheck(2)) < cutoff) {
                  if ((distance = sqrt(vCheck.absSqr().sum())) < cutoff) {
                     neighList(a1).append (a2);
                     distList(a1).append (distance);
                     cellList(a1).append (vCheck);
                     pwList(a1).append (getInteractionType (a1, a2));
                  }
               }
            }
         }
         for (k = 0; k < currentAtomNrsPeriodic.getSize(); k++) {
            a2 = currentAtomNrsPeriodic(k);
            v = currentAtomCoordsPeriodic(k);
            R = currentAtomCellsPeriodic(k);
            vCheck = v1-v;


            if (fabs(vCheck(0)) < cutoff && fabs(vCheck(1)) < cutoff
                  && fabs(vCheck(2)) < cutoff) {
               if ((distance = sqrt(vCheck.absSqr().sum())) < cutoff) {
                  neighList(a1).append (a2);
                  distList(a1).append (distance);
                  cellList(a1).append (vCheck);
                  pwList(a1).append (getInteractionType (a1, a2));
               }
            }
         }


         for (l = 0; l < neighbors(i).getSize(); l++) {
            int a = neighbors(i)(l)(0);
            int b = neighbors(i)(l)(1);
            int c = neighbors(i)(l)(2);

            // <--------
            neighborsAtomNrs = meshAtomNr(a)(b)(c);
            neighborsAtomCoords = meshAtomCoord(a)(b)(c);
            neighborsAtomCells = meshAtomCell(a)(b)(c);
            neighborsAtomNrsPeriodic = meshAtomNrPeriodic(a)(b)(c);
            neighborsAtomCoordsPeriodic = meshAtomCoordPeriodic(a)(b)(c);
            neighborsAtomCellsPeriodic = meshAtomCellPeriodic(a)(b)(c);

            // <--------
            for (k = 0; k < neighborsAtomNrs.getSize(); k++) {
               a2 = neighborsAtomNrs(k);
               v = neighborsAtomCoords(k);
               R = neighborsAtomCells(k);
               vCheck = v1-v;
               if (fabs(vCheck(0)) < cutoff && fabs(vCheck(1)) < cutoff
                     && fabs(vCheck(2)) < cutoff) {
                  if ((distance = sqrt(vCheck.absSqr().sum())) < cutoff) {
                     neighList(a1).append (a2);
                     distList(a1).append (distance);
                     cellList(a1).append (vCheck);
                     pwList(a1).append (getInteractionType (a1, a2));
                  }
               }
            }
            for (k = 0; k < neighborsAtomNrsPeriodic.getSize(); k++) {
               a2 = neighborsAtomNrsPeriodic(k);
               v = neighborsAtomCoordsPeriodic(k);
               R = neighborsAtomCellsPeriodic(k);
               vCheck = v1-v;
               if (fabs(vCheck(0)) < cutoff && fabs(vCheck(1)) < cutoff
                     && fabs(vCheck(2)) < cutoff) {
                  if ((distance = sqrt(vCheck.absSqr().sum())) < cutoff) {
                     neighList(a1).append (a2);
                     distList(a1).append (distance);
                     cellList(a1).append (vCheck);
                     pwList(a1).append (getInteractionType (a1, a2));
                  }
               }
            }
         }
      }
   }
   for (i = 0; i < nAtoms; i++) {
		neigh(i) = neighList(i);
		dist(i) = distList(i);
      pw(i) = pwList(i);
		cell(i) = cellList(i);
   }
}

void SxEAM::updateNeighs () 
{
   double distance;
   SxVector3<Double> R;
	
   SxArray<SxList<int> > neighList;
   SxArray<SxList<double> > distList;
   SxArray<SxList<SxVector3<Double> > > cellList;
   neighList.resize(nAtoms);
   distList.resize(nAtoms);
   cellList.resize(nAtoms);
   
   int i, j;

   int sCut = (int) (cutoff/tau.cell(0, 0)) + 1;
	
   for (i = 0; i < nAtoms; i++) {
      distList(i). resize (0);
      neighList(i).resize (0);
      cellList(i). resize (0);
   }
   for (i = 0; i < nAtoms; i++) {
      for (int cX = -sCut; cX <= sCut; cX++) {
         for (int cY = -sCut; cY <= sCut; cY++) {
            for (int cZ = -sCut; cZ <= sCut; cZ++) {
		
                  R = cX*tau.cell(0) + cY*tau.cell(1) 
                    + cZ*tau.cell(2); 
                  for (j = 0; j < nAtoms; j++) {
                     if ( ((i != j)) 
                           && ((distance = getDist(i, j, R)) < cutoff )) {
                        neighList(i).append (j);
                        distList(i).append (distance);
                        cellList(i).append ( tau.ref(i) - tau.ref(j) - R );
                                              
                     } 
                  }
               } 
            }
         }
      }
   for (i = 0; i < nAtoms; i++) {
		neigh(i) = neighList(i);
		dist(i) = distList(i);
		cell(i) = cellList(i);
   }   
}

void SxEAM::update (const SxAtomicStructure &newTau)
{
   double R;
   
   tau = SxAtomicStructure (newTau, SxAtomicStructure::Copy);

   if (noLinkCell) updateNeighs(); else updateNeighsLinkCellMethod();

   //--- updating charge densities on atom positions
   rhoSum.resize (0);
   rhoSum.resize (nAtoms);
   
   int iType;
   for (int i = 0; i < nAtoms; i++) {
      rhoSum(i) = 0.;
		for (int j = 0; j < neigh(i).getSize (); j++) {
			R = dist(i)(j);
         iType = getInteractionType (neigh(i)(j), neigh(i)(j));
         rhoSum(i) += getRho(R, iType);
		}
	}
}

SxList<double> SxEAM::readData(const SxString &fileName)
{
   SxList<double> list;
   ifstream filestr;
   double d1, d2;
   filestr.open(fileName.ascii());
   filestr >> d1;
   filestr >> d2;
   list.append(d1);
   list.append(d2);
   
   while (filestr) {
      filestr >> d1;
      if (filestr) {
         filestr >> d2;
         list.append(d1);
         list.append(d2);
      }
   }
   filestr.close();
   return list;
}

void SxEAM::setupSplines (const SxSymbolTable *cmd) 
{
   cout << "Setting up splines ..." << endl;
   SxList<double> readList;
   SxArray<double> Vx, Vy, y2;
   double dx, dy, xMax, xMin, x;
   int noI = (int)V.getSize ();
   int nSP;

   
   //--- number of points for tabulating potentials
   //nPoints = 1000000;
   nPoints = 100000;

   //---read in V(r) in eV(Angstroem)
   if (cmd->contains("V"))
      readList = cmd -> get ("V") -> toList ();
   else 
      readList = readData(cmd->get("VFile")->toString());

   nSP = (int)readList.getSize () / 2 / noI;
   
  
   for (int nI = 0; nI < noI; nI++) {
      Vx.resize (readList.getSize () / 2 / noI);
      Vy.resize (readList.getSize () / 2 / noI);
      y2.resize (readList.getSize () / 2 / noI);

      for (int i = 0; i < Vx.getSize (); i++) {
         Vx(i) = readList(2*(nSP*nI + i)    );
         Vy(i) = readList(2*(nSP*nI + i) + 1);
      }

      dy = Vy(1) - Vy(0);
      dx = Vx(1) - Vx(0);

      spline (Vx, Vy, (int)Vx.getSize (), dy/dx, 0., &y2);

      V(nI).resize (2);
      V(nI)(0).resize (nPoints+1);
      V(nI)(1).resize (nPoints+1);
   
      xMax = Vx(Vx.getSize () - 1);
      cutoff = xMax * A2B;
   //cutoff = xMax;
      xMin = Vx(0);
      dx = (xMax -  xMin)/(double)nPoints;
   
      for (int i = 0; i <= nPoints; i++) {
         x = xMin + dx*(double)i;
         V(nI)(0)(i) = x;
         V(nI)(1)(i) = splint (Vx, Vy, y2, (int)Vx.getSize (), x);

      }
   }
   
   //---read in Rho(r) in eV(Angstroem)
   if (cmd->contains("RHO"))
      readList = cmd -> get ("RHO") -> toList ();
   else 
      readList = readData(cmd->get("RHOFile")->toString());

   nSP = (int)readList.getSize () / 2 / noI;

   for (int nI = 0; nI < noI; nI++) {
      Vx.resize (readList.getSize () / 2 / noI);
      Vy.resize (readList.getSize () / 2 / noI);
      y2.resize (readList.getSize () / 2 / noI);

      for (int i = 0; i < Vx.getSize (); i++) {
         Vx(i) = readList(2*(nSP*nI + i)    );
         Vy(i) = readList(2*(nSP*nI + i) + 1);
      }

      dy = Vy(1) - Vy(0);
      dx = Vx(1) - Vx(0);

      spline (Vx, Vy, (int)Vx.getSize (), dy/dx, 0., &y2);
   
      rho(nI).resize (2);
      rho(nI)(0).resize (nPoints+1);
      rho(nI)(1).resize (nPoints+1);
   
      xMax = Vx(Vx.getSize () - 1);
      if ( xMax*A2B > cutoff ) cutoff = xMax*A2B;
      xMin = Vx(0);
      dx = (xMax -  xMin)/(double)nPoints;
   
      for (int i = 0; i <= nPoints; i++) {
         x = xMin + dx*(double)i;
         rho(nI)(0)(i) = x;
         rho(nI)(1)(i) = splint (Vx, Vy, y2, (int)Vx.getSize (), x);

     // cout << "RHO: " << rho(0)(i) << " " << rho(1)(i) << endl;

      }
   }

   
   //---read in F(r) in eV(Angstroem)
   if (cmd->contains("F"))
      readList = cmd -> get ("F") -> toList ();
   else 
      readList = readData(cmd->get("FFile")->toString());
   
   nSP = (int)readList.getSize () / 2 / noI;

   for (int nI = 0; nI < noI; nI++) {
      Vx.resize (readList.getSize () / 2 / noI);
      Vy.resize (readList.getSize () / 2 / noI);
      y2.resize (readList.getSize () / 2 / noI);

      for (int i = 0; i < Vx.getSize (); i++) {
         Vx(i) = readList(2*(nSP*nI + i)    );
         Vy(i) = readList(2*(nSP*nI + i) + 1);
      }

      dy = Vy(1) - Vy(0);
      dx = Vx(1) - Vx(0);

      spline (Vx, Vy, (int)Vx.getSize (), dy/dx, 0., &y2);


      F(nI).resize (2);
      F(nI)(0).resize (nPoints+1);
      F(nI)(1).resize (nPoints+1);
   
      xMax = Vx(Vx.getSize () - 1);
      xMin = Vx(0);
      dx = (xMax -  xMin)/(double)nPoints;
   
      for (int i = 0; i <= nPoints; i++) {
         x = xMin + dx*(double)i;
         F(nI)(0)(i) = x;
         F(nI)(1)(i) = splint (Vx, Vy, y2, (int)Vx.getSize (), x);

         //  cout << "FF: " << F(0)(i) << " " << F(1)(i) << endl;

      }
   } 

   //--- tabulating derivatives
   for (int nI = 0; nI < noI; nI++) {
      dV(nI).resize (2);
      dV(nI)(0).resize (V(nI)(0).getSize ());
      dV(nI)(1).resize (V(nI)(1).getSize ()); 
     
      drho(nI).resize (2);
      drho(nI)(0).resize (rho(nI)(0).getSize ()); 
      drho(nI)(1).resize (rho(nI)(1).getSize ()); 
     
      dF(nI).resize (2);
      dF(nI)(0).resize (F(nI)(0).getSize ()); 
      dF(nI)(1).resize (F(nI)(1).getSize ()); 
   

      for (int i = 1; i < (dV(nI)(0).getSize () - 1); i++) {
         dV(nI)(0)(i) = V(nI)(0)(i);
         dy = V(nI)(1)(i + 1) - V(nI)(1)(i - 1);
         dx = V(nI)(0)(i + 1) - V(nI)(0)(i - 1);
         dV(nI)(1)(i) = dy/dx;
   } 
      dV(nI)(0)(0) = V(nI)(0)(0);
      dV(nI)(0)(dV(nI).getSize() - 1) = V(nI)(0)(dV(nI).getSize () - 1);
      dV(nI)(1)(0) = dV(nI)(1)(1);
      dV(nI)(1)(dV(nI).getSize() - 1) = dV(nI)(1)(dV(nI).getSize () - 2);
      
      for (int i = 1; i < (dF(nI)(0).getSize () - 1); i++) {
         dF(nI)(0)(i) = F(nI)(0)(i);
         dy = F(nI)(1)(i + 1) - F(nI)(1)(i - 1);
         dx = F(nI)(0)(i + 1) - F(nI)(0)(i - 1);
         dF(nI)(1)(i) = dy/dx;
      } 
      dF(nI)(0)(0) = F(nI)(0)(0);
      dF(nI)(0)(dF(nI).getSize() - 1) = F(nI)(0)(dF(nI).getSize () - 1);
      dF(nI)(1)(0) = dF(nI)(1)(1);
      dF(nI)(1)(dF(nI).getSize() - 1) = dF(nI)(1)(dF(nI).getSize () - 2);
      
      for (int i = 1; i < (drho(nI)(0).getSize () - 1); i++) {
         drho(nI)(0)(i) = rho(nI)(0)(i);
         dy = rho(nI)(1)(i + 1) - rho(nI)(1)(i - 1);
         dx = rho(nI)(0)(i + 1) - rho(nI)(0)(i - 1);
         drho(nI)(1)(i) = dy/dx;
      } 
      drho(nI)(0)(0) = rho(nI)(0)(0);
      drho(nI)(0)(drho(nI).getSize() - 1) = rho(nI)(0)(drho(nI).getSize () - 1);
      drho(nI)(1)(0) = drho(nI)(1)(1);
      drho(nI)(1)(drho(nI).getSize() - 1) = drho(nI)(1)(drho(nI).getSize () - 2);
   }
   
 
  // for (int i = 0; i < F(0)(0).getSize (); i++) 
     // if ((i < 100) || (i > (drho(0)(0).getSize () - 100)))
   //      cout << F(1)(0)(i) << " " << F(1)(1)(i) << " " << dF(1)(1)(i) << endl;
  // QUIT;
 

}


double SxEAM::getEnergy () const
{
	return totalEnergy;
}

double SxEAM::getPotentialEnergy () const
{
	return totalEnergy;
}


SxSpeciesData SxEAM::getSpeciesData () const
{
   return speciesData;
}



SxAtomicStructure SxEAM::getNumericalForces (const SxAtomicStructure &tauIn, const double &dx) {
   SxAtomicStructure Forces (tauIn, SxAtomicStructure::Copy);
	SxAtomicStructure undevCoord (tauIn, SxAtomicStructure::Copy); 
	SxAtomicStructure devCoord (tauIn, SxAtomicStructure::Copy); 
	
   double undevEnergy;


   update(tauIn);
	
	undevEnergy = getTotalEnergy ();


	for (int i = 0; i < nAtoms; i++) {
		for (int j = 0; j < 3; j++) {
			Forces.ref(i)(j) = 0.;
			devCoord.ref(i)(j) += dx;
			update(devCoord);
			Forces.ref(i)(j) 
            -= (getTotalEnergy () - undevEnergy)/dx/2.;
			devCoord.ref(i)(j) -= 2.*dx;
			update(devCoord);
			Forces.ref(i)(j) += (getTotalEnergy() - undevEnergy)/dx/2.;
			devCoord.ref(i)(j) += dx;

		}
	}
	return Forces;
}

SxArray<SxAtomicStructure> SxEAM::getNumericalHesseMatrix 
(const SxAtomicStructure &tauIn, const double &dx) {
 
   SxArray<SxAtomicStructure> Hesse;
   Hesse.resize (3*nAtoms);
   SxAtomicStructure A, B;
	SxAtomicStructure devCoord (tauIn, SxAtomicStructure::Copy); 
   SxSymbolTable *table = NULL; // CF 2016-07-07: where should table come from ???

	for (int i = 0; i < nAtoms; i++) {
		for (int j = 0; j < 3; j++) {
			devCoord.ref(i)(j) += dx;
			A = getForces (devCoord, table);
			devCoord.ref(i)(j) -= 2.*dx;
			B = getForces (devCoord, table);
         Hesse(3*i+j).copy ((1./dx/2.)*(A-B));
         //Hesse(3*i+j) =  (1./dx/2.)*(A-B);
         devCoord.ref(i)(j) += dx;
		}
	}
	return Hesse;
}

SxArray<SxArray<SxAtomicStructure> > SxEAM::getNumericalTTensor (const SxAtomicStructure &tauIn, const double &dx) 
{
   double threshold = 1e-9;

   SxArray<SxArray<SxAtomicStructure> > TTensor;
   TTensor.resize(3*nAtoms);
   SxArray<SxAtomicStructure> A, B, Delta;
   SxAtomicStructure devCoord (tauIn, SxAtomicStructure::Copy);
   int i,j,k,l,m;

   sxprintf("nAtoms=%d",nAtoms); fflush(stdout);
   for (i = 0; i < nAtoms; i++) {
        sxprintf("  %d",i+1); fflush(stdout);
   	for (j = 0; j < 3; j++) {
   		devCoord.ref(i)(j) += dx;
   		A = getNumericalHesseMatrix (devCoord,dx);
   		devCoord.ref(i)(j) -= 2.*dx;
   		B = getNumericalHesseMatrix (devCoord,dx);
         //Delta = A-B;

         Delta.resize(3*nAtoms);
         for (k=0;k<Delta.getSize();k++) {
            Delta(k) = A(k) - B(k);
   			Delta(k) = (1./dx/2.) * Delta(k);
            for (l=0; l<nAtoms; l++) {
               for (m=0; m<3; m++) {
                  if (Delta(k).ref(l)(m) < threshold) Delta(k).ref(l)(m)=0.;
               }
            }
         }
        // cout << "Delta " << Delta(24)( << endl;
         TTensor(3*i+j) = Delta;
         devCoord.ref(i)(j) += dx;
   	}
   }
   cout << endl;
   return TTensor;
}

void SxEAM::exportForces (const SxAtomicStructure &forces, const SxString &fileName)
{
   int i;
   ofstream filestr;
   filestr.open(fileName.ascii(), ifstream::trunc);
   for (i=0; i<nAtoms; i++) {
      filestr << forces(i)(0) << " " << forces(i)(1) << " " << forces(i)(2) << endl;
   }
   filestr.close();
}

void SxEAM::exportHesse (const SxArray<SxAtomicStructure> &hesse, const SxString &fileName)
{
   int i, j;
   ofstream filestr;
   filestr.open(fileName.ascii(), ifstream::trunc);
   for (i=0; i<3*nAtoms; i++) {
      for (j=0; j<nAtoms; j++) {
         filestr << -hesse(i)(j)(0) << " " << -hesse(i)(j)(1) << " " << -hesse(i)(j)(2) << " ";
      }
      filestr << endl;
   }
   filestr.close();
}

void SxEAM::exportTTensor (const SxArray<SxArray<SxAtomicStructure> > &TTensor, const SxString &fileName)
{
   int i, j, k;
   ofstream filestr;
   filestr.open(fileName.ascii(), ifstream::trunc);
	for (k=0;k<3*nAtoms;k++) {
      for (i=0; i<3*nAtoms; i++) {
         for (j=0; j<nAtoms; j++) {
            filestr << TTensor(k)(i)(j)(0) << " " << TTensor(k)(i)(j)(1) << " " << TTensor(k)(i)(j)(2) << " ";
         }
      }
      filestr << endl;
   }
   filestr.close();
}

SxAtomicStructure SxEAM::getForces (const SxAtomicStructure &tauIn,
                                                const SxSymbolTable *dummy)
{  // SVN_HEAD; // this line MUST NOT be removed
   SxAtomicStructure forces (tauIn, SxAtomicStructure::Copy);
   update (tauIn);
  
   

   
   double R = 0.;
   SxVector3<Double> e;
   int n = 0;
   int iType1;
   int iType2;
   int iType3;
   

   for (int i = 0; i < nAtoms; i++) {
      forces.ref(i).set(0.);
      for (int j = 0; j < neigh(i).getSize (); j++) {
         n = neigh(i)(j);
         R = dist(i)(j);
         e = (1./R)*cell(i)(j);

         iType1 = getInteractionType(i, n);
         iType2 = getInteractionType(n, n);
         iType3 = getInteractionType(i, i);

       
         forces.ref(i) 
            -= (getdV(R, iType1) 
             + getdF(rhoSum(n), iType2) *getdRho(R, iType3) 
             + getdF(rhoSum(i), iType3) *getdRho(R, iType2))*e;
      }
   }
   totalEnergy = getTotalEnergy ();
   
//   cout << tauIn << endl;
   return ((1./HA2EV)*forces);
}


double SxEAM::getV (const double &r, int isP)
{
   double rB = r/A2B;
   double idxD = ((rB - V(isP)(0)(0))/(V(isP)(0)(1)-V(isP)(0)(0)));
   int idx = (int) idxD;
   double interpol = idxD - (double)idx; 
   return ((1. - interpol)*V(isP)(1)(idx) + interpol*V(isP)(1)(idx + 1));
}

double SxEAM::getdV (const double &r, int isP)
{
   double rB = r/A2B;
   double idxD = ((rB - dV(isP)(0)(0))/(dV(isP)(0)(1)-dV(isP)(0)(0)));
   int idx = (int) idxD;
   double interpol = idxD - (double)idx;
   //double interpol = 0.;
   return (((1. - interpol)*dV(isP)(1)(idx) + interpol*dV(isP)(1)(idx + 1))/A2B);
}

double SxEAM::getRho (const double &r, int isP)
{
   double rB = r/A2B;
   //double last = rho(isP)(0)(rho(isP)(0).getSize()-1);
   double idxD = ((rB - rho(isP)(0)(0))/(rho(isP)(0)(1)  -rho(isP)(0)(0)));
   int idx = (int) idxD;
   double interpol = idxD - (double)idx; 
   return ((1. - interpol)*rho(isP)(1)(idx) + interpol*rho(isP)(1)(idx + 1));
}

double SxEAM::getdRho (const double &r, int isP)
{
   double rB = r/A2B;
   double idxD = ((rB - drho(isP)(0)(0))/(drho(isP)(0)(1) - drho(isP)(0)(0)));
   int idx = (int) idxD;
   double interpol = idxD - (double)idx; 
   //double interpol = 0.; 
   return (((1. - interpol)*drho(isP)(1)(idx) + interpol*drho(isP)(1)(idx + 1))/A2B);
}

double SxEAM::getF (const double &rh, int isP)
{
   double rB = rh;
   double idxD = ((rB - F(isP)(0)(0))/(F(isP)(0)(1)-F(isP)(0)(0)));
   int idx = (int) idxD;
   if (idx >= F(isP)(0).getSize ()) return (F(isP)(1)(F(isP)(0).getSize () - 1));
   double interpol = idxD - (double)idx; 
   return ((1. - interpol)*F(isP)(1)(idx) + interpol*F(isP)(1)(idx + 1));
}

double SxEAM::getdF (const double &rh, int isP)
{
   double rB = rh;
   double idxD = ((rB - dF(isP)(0)(0))/(dF(isP)(0)(1) - dF(isP)(0)(0)));
   int idx = (int) idxD;
   if (idx >= dF(isP)(0).getSize ()) return 0.;
   double interpol = idxD - (double)idx;
   //double interpol = 0.;
   return ((1. - interpol)*dF(isP)(1)(idx) + interpol*dF(isP)(1)(idx + 1));
}


double SxEAM::getTotalEnergy () 
{
   double energy, R;
   int i, j;
	energy = 0.;
   int iType;

	for (i = 0; i < nAtoms; i++) {
		for (j = 0; j < neigh(i).getSize (); j++) {
         iType = getInteractionType(i, neigh(i)(j));
			R = dist(i)(j);
         energy += 0.5*getV (R, iType);
		}
      iType = getInteractionType(i, i);
      energy += getF (rhoSum(i), iType);
	}
	return (energy/HA2EV);
}		

void SxEAM::execute (const SxSymbolTable *cmd, bool calc)
{
   cout << "Executing EAM Potential" << endl << endl;
   SxString s = cmd->getName();
   double disp;
   update (tau);
   SxAtomicStructure str=tau;

   exportForces(tau,"cartesian_coords");
   cout << "cartesian_coords in bohr exported." << endl << endl;

   if (s=="forces") {
      if (cmd->contains("numerical")) {
         disp = (cmd->contains("disp")) ? cmd->get("disp")->toReal() : 0.01;
         cout << "Calculating numerical forces with displacement " 
              << disp << " Angstrom..." << endl;
         exportForces(getNumericalForces(tau, disp),"Forces");
         cout << "Forces in hartree/bohr exported." << endl << endl;
      } else {
         cout << "Calculating analytical forces..." << endl;
         exportForces(getForces(tau),"Forces");
         cout << "Forces in hartree/bohr exported." << endl << endl;
      }
   }
   
   if (s=="HesseMatrix") {
      disp = (cmd->contains("disp")) ? cmd->get("disp")->toReal() : 0.01;
      
      cout << "Calculating numerical Hesse matrix with displacement " 
           << disp << " Angstrom (forces are analytical)..." << endl;
      exportHesse(getNumericalHesseMatrix(tau, disp),"HesseMatrix");
      cout << "HesseMatrix in hartree/(u*bohr^2) exported." << endl << endl;
   }
   
   if (s=="TTensor") {
      disp = (cmd->contains("disp")) ? cmd->get("disp")->toReal() : 0.01;
      
      cout << "Calculating numerical TTensor with displacement " 
           << disp << " Angstrom (forces are analytical)..." << endl;
      SxArray<SxArray<SxAtomicStructure> > TTensor;
      TTensor = getNumericalTTensor(tau,disp);
      exportTTensor(TTensor,"TTensor");
      cout << "TTensor in hartree/(u*bohr^3) exported." << endl << endl;
   }

   return;
}


void SxEAM::spline 
(const SxArray<double> &x, const SxArray<double> &y, int n, 
 double yp1, double ypn, SxArray<double> *y2) 
{
   int i, k;
   double p, qn, sig, un;

   SxArray<double> u (x.getSize ());

   if (yp1 > 0.99e30)
      (*y2)(0) = u(0) = 0.;
   else {
      (*y2)(0) = -0.5;
      u(0) = (3.0/(x(1) - x(0)))
           * ((y(1) - y(0))/(x(1) - x(0)) - yp1);
   }

   for (i = 1; i <= (n-2); i++) {
      sig = (x(i) - x(i-1))/(x(i+1) - (x(i-1)));
      p   = sig* (*y2)(i-1) + 2.;
      (*y2)(i) = (sig - 1.)/p;
      u(i)  = (y(i+1) - y(i))/(x(i+1) - x(i)) - (y(i) - y(i-1))/(x(i) - x(i-1));
      u(i)  = (6.0*u(i)/(x(i+1) - x(i-1)) -sig*u(i-1))/p;
   }

   if (ypn > 0.99e30)
      qn = un = 0.0;
   else {
      qn = 0.5;
      un = (3.0/(x(n-1) - x(n-2)))*(ypn - (y(n-1) - y(n-2))/(x(n-1) - x(n-2)));
   }
   (*y2)(n-1) = (un - qn*u(n-2))/(qn* (*y2)(n-2) + 1.);

   for (k = (n-2); k>= 0; k--) 
      (*y2)(k) = (*y2)(k)*(*y2)(k + 1) + u(k);

}

double SxEAM::splint 
(const SxArray<double> &xa, const SxArray<double> &ya, 
 const SxArray<double> &y2a, int n, double x)
{
   int klo, khi, k;
   double h,b,a;

   klo = 0;
   khi = n-1;

   while ( (khi-klo) > 1) {
      k = (khi + klo) >> 1;
      if (xa(k) > x) khi = k;
      else klo = k;
   }

   h = (xa(khi) - xa(klo));

   if (h == 0.) {
      cout << "Bad xa input to routine splint" << endl;
      SX_EXIT;
   }

   a = (xa(khi) - x)/h;
   b = (x - xa(klo))/h;

   double y = a*ya(klo) + b*ya(khi)+((a*a*a - a)*y2a(klo) 
            + (b*b*b - b)*y2a(khi))*(h*h)/6.0;
   return y;
}


