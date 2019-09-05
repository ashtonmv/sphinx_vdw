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
#include <SxTBAtomAtom.h>
#include <SxDFTBFilesReader.h>

SxTBAtomAtom::SxTBAtomAtom ()
{
   // empty   
}

SxTBAtomAtom::~SxTBAtomAtom ()
{
   // empty   
}

void SxTBAtomAtom::init (int iS, int jS, const SxTBSpecies &inSpeciesData,
                                         const bool        noERepul)
{
   // --- read from file
   const int BUFLEN = 1024;
   char      buffer[BUFLEN];
   FILE     *fp;

   withoutERepulsive = noERepul;

   SxVector<Double> epsL, hubL; // L: for orbital quantum number (l)
   int              iOrb, l, iPar, i;

   //iSName:    slater-koster 1st species name (skElementName)
   //jSName:     slater-koster 2nd species name (skElementName)
   //path:       slater-koster files path
   //fileFormat: slater-koster files format (1=orig. dftb code, 2=sphinx)
   SxString iSName     = inSpeciesData.skElementName(iS);
   SxString jSName     = inSpeciesData.skElementName(jS);
   SxString path       = inSpeciesData.skFilesPath; 
   int      fileFormat = inSpeciesData.skFilesFormat;

   lMax1 = inSpeciesData.lMax(iS);
   lMax2 = inSpeciesData.lMax(jS);

   lOrb1 = getlOrb (lMax1);
   lOrb2 = getlOrb (lMax2);

   isSameSpecies = false;
   if (iS == jS)  {
      epsL.resize (lMax1+1, 0.);
      hubL.resize (lMax1+1, 0.);
      nValElect = 0;
      isSameSpecies = true;
   }

   // --- early file formats made for sphingx
   if (fileFormat == 2) { 

      if (lMax1 > 1 || lMax2 > 1) {
         cout << "\n invalid value(s) of lMax" << endl;
         cout << " only up to lMax = p is supported for this sk-files format" << endl;
         cout << " Try www.DFTB.org SK-files formats; skFilesFormat = 1"      << endl;
         SX_EXIT;
      }

      fileName = path + "/" + iSName + jSName + ".dat"; 
      isTransposed = false;

      if ( !(fp = fopen (fileName.ascii(), "r")) )  {
         sxprintf ("| Can't open file: %s \n", fileName.ascii());
         fileName = path + "/" + jSName + iSName + ".dat"; 
         sxprintf ("|  Thus, I'll try: %s \n", fileName.ascii());
         isTransposed = true;
         if ( !(fp = fopen (fileName.ascii(), "r")) )  {
            sxprintf ("\nCan't open file %s\n", fileName.ascii());
            SX_EXIT;
         }
      }
      cout << "| Reading file : " << fileName.ascii() << endl;

      fgets (buffer, BUFLEN, fp);  // to skip line

      int lMax1SK, lMax2SK;
      if (!isTransposed)
         fscanf (fp, "%le %d %d %d", &delta, &nPoints, &lMax1SK, &lMax2SK);
      else
         fscanf (fp, "%le %d %d %d", &delta, &nPoints, &lMax2SK, &lMax1SK);
      fgets (buffer, BUFLEN, fp);  // to skip rest of line

      if (lMax1 != lMax1SK || lMax2 != lMax2SK) {
         cout << "\n  contraditions in lMax (max. angular momenta) values:\n"
              << "     an lMax value given in species data is different from\n"
              << "     the one given in the slater-koster file\n";
         SX_QUIT;
      }

      // --- if the file is for two atoms of same type: read the number of 
      //     valence elctorns, and the atomic eigenvalues  
      double numb;
      if (isSameSpecies)  {
         for (i = 0; i < 2; i++)  fgets (buffer, BUFLEN, fp); //skip 2 lines
         fscanf (fp, "%d", &nValElect); // read number of valence electrons
         if (inSpeciesData.valenceCharge(iS) != nValElect) {
            cout << "\n contradiction in the valence charge provided" << endl;
            cout << " in the species data and in the sk-files\n" << endl;
            SX_EXIT;
         }
         for (l = 0; l <= lMax1; l++)  {
            fscanf (fp, "%le", &numb);
            epsL(l) = numb;
         }
         for (l = 0; l <= lMax1; l++)  {
            fscanf (fp, "%le", &numb);
            hubL(l) = numb;
         }
         fgets (buffer, BUFLEN, fp);  // to skip rest of line
      }

      for (i = 0; i < 3; i++)  fgets (buffer, BUFLEN, fp); //skip 3 lines

      // --- In these files there are always 5 values, independent of lMax
      //     (unused values are set to zero in the files).
      //     The order is: 
      //     Hss0 Hsp0 Hps0 Hpp1 Hpp0 Sss0 Ssp0 Sps0 Spp1 Spp0 
      //     (where 0=sigma, 1=pi, 2=delta, H=Hamiltonian, S=Overlap)
      int nPar = 5;  
      slkoHam.reformat (nPoints, nPar);  
      slkoOvl.reformat (nPoints, nPar);

      double ds;
      int    ip, il;
      SxVector<Double> distance(nPoints);

      for (ip = 0; ip < nPoints; ip++)  {
         fscanf (fp, "%le", &ds);
         distance(ip) = ds;
         for (iPar = 0; iPar < nPar; iPar++)  {
            fscanf (fp, "%le", &numb);
            slkoHam(ip, iPar) = numb;
         }
         for (iPar = 0 ; iPar < nPar; iPar++)  {
            fscanf (fp, "%le", &numb);
            slkoOvl(ip, iPar) = numb;
         }
         fgets (buffer, BUFLEN, fp);  // to skip rest of line
      }

      if (!withoutERepulsive)  {
         // --- read the repulsive energy coefficients part 
         char type [20];
         fscanf (fp, "%d", &nLines);
         fscanf (fp, "%le", &dCut);
         fscanf (fp, "%s", type);
         fgets  (buffer, BUFLEN, fp); 

         repType = type;

         if (repType == "dataFit")  {

            // --- new format: repulsive energies on 
            //     a grid of distances
            repData.resize (nLines);
            for (il = 0; il < nLines; il++)  {
               repData(il).resize (2);
               fscanf (fp, "%le", &numb); // distance
               repData(il)(0) = numb;
               fscanf (fp, "%le", &numb); // eRep
               repData(il)(1) = numb;
            }

         }  else if (repType == "functionFit")  {

            // --- old format: function coefficients
            interval.reformat (nLines, 2); //each interval made of 2 numbers
            repulsiveCoef1.reformat (nLines-1, 4); // 4 coefficients 
            repulsiveCoef2.resize (6);             // 6 coefficients
            func.resize (3);
            // --- complete reading and storing
            for (ip = 0; ip < 3; ip++)  {
               fscanf (fp, "%le", &numb);
               func (ip) = numb;
            }
            fgets (buffer, BUFLEN, fp);  
            for (il = 0; il < (nLines-1) ; il++)  {
               for (i = 0; i < 2; i++)  {
                  fscanf (fp, "%le", &numb);
                  interval(il, i) = numb;
               }
               for (i = 0; i < 4; i++)  {
                  fscanf (fp, "%le", &numb);
                  repulsiveCoef1(il, i) = numb;
               }
               fgets (buffer, BUFLEN, fp);  
            }
            // --- read the last line 
            for (i = 0; i < 2; i++)  {
               fscanf (fp, "%le", &numb);
               interval(nLines-1, i) = numb;
            }
            for (i = 0; i < 6; i++)  {
               fscanf (fp, "%le", &numb);
               repulsiveCoef2(i) = numb;
            }

         }  else  {
            cout << "\n can't read repulsive energy data" << endl;
            cout << " could be old slater-koster files format! " << endl;
            SX_QUIT;
         }
      }
      fclose(fp);

      // determine cutoff for the repulsive energy 
      cutoff = distance(nPoints-3);//the use of spline needs the next two values  

   //--- format used for the DFTB code (Frauenheim's group)
   } else if (fileFormat == 1) {

      fileName1 = path + "/" + iSName + jSName + "_dftb.dat"; 
      fileName2 = path + "/" + jSName + iSName + "_dftb.dat"; 
      isTransposed = false;

      // --- open the first file
      if ( !(fp = fopen (fileName1.ascii(), "r")) )  {
         sxprintf ("| Can't open file: %s \n", fileName1.ascii());
         SX_EXIT;
      }
      cout << "| Reading file : " << fileName1.ascii() << endl;

      // --- Determine whether the file has the new "extended format".
      //     Not supported for now (no f orbitals implemented here)
      char ident [5];
      SxString identifier;
      fscanf (fp, "%s", ident);
      identifier = ident;
      if (identifier == "@") {
         sxprintf ("\n the current version does not support the format");
         sxprintf ("\n of the slater-koster file %s", fileName1.ascii());
         sxprintf ("\n this new extended format might be supported later\n");
         SX_EXIT;
      } else {
         rewind (fp);
      }

      fscanf (fp, "%le %d", &delta, &nPoints);
      fgets (buffer, BUFLEN, fp);  // to skip rest of line
      
      // --- if the file is for two atoms of same type: read the number of 
      // --- valence elctorns, and the atomic eigenvalues  
      double numb;
      if (isSameSpecies)  {

         // --- These files always have values (could be 0) for d p s orbitals
         //     hence, we always read the same number of values
         for (l = 2; l >= 0; l--)  { 
            fscanf (fp, "%le", &numb);
            if (l <= lMax1) epsL(l) = numb; // store only relevant data
         }
         fscanf (fp, "%le", &numb);  // reads a useless number 
         for (l = 2; l >= 0; l--)  { // hubbardU parameters
            fscanf (fp, "%le", &numb);
            if (l <= lMax1) hubL(l) = numb;
         }
         nValElect = 0;
         for (l = 2; l >= 0; l--)  { // occupation numbers
            fscanf (fp, "%le", &numb);
            if (l <= lMax1) nValElect += (int) numb;
         }
         fgets (buffer, BUFLEN, fp);  // to skip rest of line

         if (inSpeciesData.valenceCharge(iS) != nValElect) {
            cout << "\n contradiction in the valence charge provided" << endl;
            cout << " in the species data and in the sk-files\n" << endl;
            cout << inSpeciesData.valenceCharge(iS) << " Vs. " << nValElect << endl;
            SX_EXIT;
         }
      }
      fgets (buffer, BUFLEN, fp); // skip next line. Not used input 

      // --- In these files there are always 10 values, independent of lMax,
      //     in the following order: 
      //     Hdd0 Hdd1 Hdd2 Hpd0 Hpd1 Hpp0 Hpp1 Hsd0 Hsp0 Hss0 
      //     Sdd0 Sdd1 Sdd2 Spd0 Spd1 Spp0 Spp1 Ssd0 Ssp0 Sss0
      //     (where 0=sigma, 1=pi, 2=delta, H=Hamiltonian, S=Overlap)
      int nParComp = 10;
      SxMatrix<Double> slkoHamComp, slkoOvlComp;
      SxMatrix<Double> slkoHamCompFile2, slkoOvlCompFile2;
      slkoHamComp.reformat (nPoints, nParComp);  
      slkoOvlComp.reformat (nPoints, nParComp);
      slkoHamCompFile2.reformat (nPoints, nParComp);  
      slkoOvlCompFile2.reformat (nPoints, nParComp);

      // the complete set of spd interactions contains 14 elements
      // saved in the following order:
      // (0)ss_sigma  (1)sp_sigma (2)ps_sigma (3)pp_pi  (4)pp_sigma 
      // (5)sd_sigma  (6)ds_sigma (7)pd_pi    (8)dp_pi  (9)pd_sigma (10)dp_sigma 
      // (11)dd_delta (12)dd_pi   (13)dd_sigma
      int nPar = 14;
      slkoHam.reformat (nPoints, nPar);  
      slkoOvl.reformat (nPoints, nPar);

      int    ip, il;
      SxVector<Double>  distance(nPoints);
      SxDFTBFilesReader reader (fp);

      distance.set(0);
      double tempDist = 0.0;

      for (ip = 0; ip < nPoints; ip++)  {
         // I am not sure if the first distance is zero or delta!
         // most probably it is delta, but I used 0 in some old pot files!
         tempDist     += delta;
         distance(ip) += tempDist;
         for (iPar = 0; iPar < nParComp; iPar++)  {
            numb = reader.getNextReal ();
            slkoHamComp(ip, iPar) = numb;
         }
         for (iPar = 0 ; iPar < nParComp; iPar++)  {
            numb = reader.getNextReal ();
            slkoOvlComp(ip, iPar) = numb;
         }
         fgets (buffer, BUFLEN, fp);  // to skip rest of line
      }

      // --- read the repulsive energy coefficients part 
      if (!withoutERepulsive)  {
         
         char type [20];
         fscanf (fp, "%s", type);
         fscanf (fp, "%d", &nLines);
         fscanf (fp, "%le", &dCut);
         fgets  (buffer, BUFLEN, fp); 

         repType = type;

         if (repType == "Spline")  {

            // --- old format: function coefficients
            interval.reformat (nLines, 2); //each interval made of 2 numbers
            repulsiveCoef1.reformat (nLines-1, 4); // 4 coefficients 
            repulsiveCoef2.resize (6);             // 6 coefficients
            func.resize (3);
            // --- complete reading and storing
            for (ip = 0; ip < 3; ip++)  {
               fscanf (fp, "%le", &numb);
               func (ip) = numb;
            }
            fgets (buffer, BUFLEN, fp);  
            for (il = 0; il < (nLines-1) ; il++)  {
               for (i = 0; i < 2; i++)  {
                  fscanf (fp, "%le", &numb);
                  interval(il, i) = numb;
               }
               for (i = 0; i < 4; i++)  {
                  fscanf (fp, "%le", &numb);
                  repulsiveCoef1(il, i) = numb;
               }
               fgets (buffer, BUFLEN, fp);  
            }
            // --- read the last line 
            for (i = 0; i < 2; i++)  {
               fscanf (fp, "%le", &numb);
               interval(nLines-1, i) = numb;
            }
            for (i = 0; i < 6; i++)  {
               fscanf (fp, "%le", &numb);
               repulsiveCoef2(i) = numb;
            }

         }  else  {
            cout << "\n unknown repulsive energy type = " << repType << endl;
            cout << " can't read repulsive energy data"              << endl;
            cout << " could be an old slater-koster files format! "     << endl;
            SX_QUIT;
         }
      }
      fclose(fp);

      // determine cutoff for the repulsive energy 
      cutoff = distance(nPoints-3);//the use of spline needs the next two values  

      if (!isSameSpecies) {

         // --- open the second file: we need some matrix elements
         if ( !(fp = fopen (fileName2.ascii(), "r")) )  {
            sxprintf ("| Can't open file: %s \n", fileName2.ascii());
            SX_EXIT;
         }
         cout << "| Reading file : " << fileName2.ascii() << endl;

         // skip 2 lines
         for (i = 0; i < 2; i++) fgets (buffer, BUFLEN, fp);

         for (ip = 0; ip < nPoints; ip++)  {
            for (iPar = 0; iPar < nParComp; iPar++)  {
               numb = reader.getNextReal ();
               slkoHamCompFile2(ip, iPar) = numb;
            }
            for (iPar = 0 ; iPar < nParComp; iPar++)  {
               numb = reader.getNextReal ();
               slkoOvlCompFile2(ip, iPar) = numb;
            }
            fgets (buffer, BUFLEN, fp);  // to skip rest of line
         }
         fclose(fp);

      } else { // same species
         slkoHamCompFile2 = slkoHamComp;
         slkoOvlCompFile2 = slkoOvlComp;
      }

      // --- an ugly conversion for matrix elements, 
      //     due to different ordering in sphinx and dftb format
      for (ip = 0; ip < nPoints; ip++)  {
         // Hamiltonian matrix elements
         // --- ss
         slkoHam(ip,0)  = slkoHamComp(ip, 9);        // ss_sigma
         // --- sp
         slkoHam(ip,1)  = slkoHamComp(ip, 8);        // sp_sigma
         slkoHam(ip,2)  = slkoHamCompFile2(ip, 8);   // ps_sigma
         // --- pp
         slkoHam(ip,3)  = slkoHamComp(ip, 6);        // pp_pi
         slkoHam(ip,4)  = slkoHamComp(ip, 5);        // pp_sigma
         // --- sd
         slkoHam(ip,5)  = slkoHamComp(ip, 7);        // sd_sigma
         slkoHam(ip,6)  = slkoHamCompFile2(ip, 7);   // ds_sigma
         // --- pd
         slkoHam(ip,7)  = slkoHamComp(ip, 4);        // pd_pi
         slkoHam(ip,8)  = slkoHamCompFile2(ip, 4);   // dp_pi
         slkoHam(ip,9)  = slkoHamComp(ip, 3);        // pd_sigma
         slkoHam(ip,10) = slkoHamCompFile2(ip, 3);   // dp_sigma
         // --- dd
         slkoHam(ip,11) = slkoHamComp(ip, 2);        // dd_delta
         slkoHam(ip,12) = slkoHamComp(ip, 1);        // dd_pi
         slkoHam(ip,13) = slkoHamComp(ip, 0);        // dd_sigma

         // Overlap matrix elements
         // --- ss
         slkoOvl(ip,0)  = slkoOvlComp(ip, 9);        // ss_sigma
         // --- sp                                              
         slkoOvl(ip,1)  = slkoOvlComp(ip, 8);        // sp_sigma
         slkoOvl(ip,2)  = slkoOvlCompFile2(ip, 8);   // ps_sigma
         // --- pp                                              
         slkoOvl(ip,3)  = slkoOvlComp(ip, 6);        // pp_pi
         slkoOvl(ip,4)  = slkoOvlComp(ip, 5);        // pp_sigma
         // --- sd
         slkoOvl(ip,5)  = slkoOvlComp(ip, 7);        // sd_sigma
         slkoOvl(ip,6)  = slkoOvlCompFile2(ip, 7);   // ds_sigma
         // --- pd
         slkoOvl(ip,7)  = slkoOvlComp(ip, 4);        // pd_pi
         slkoOvl(ip,8)  = slkoOvlCompFile2(ip, 4);   // dp_pi
         slkoOvl(ip,9)  = slkoOvlComp(ip, 3);        // pd_sigma
         slkoOvl(ip,10) = slkoOvlCompFile2(ip, 3);   // dp_sigma
         // --- dd
         slkoOvl(ip,11) = slkoOvlComp(ip, 2);        // dd_delta
         slkoOvl(ip,12) = slkoOvlComp(ip, 1);        // dd_pi
         slkoOvl(ip,13) = slkoOvlComp(ip, 0);        // dd_sigma

      }

   } else {

      cout << "\n unknown slater-koster file format " << endl;
      SX_QUIT;

   }

   ham.reformat (lOrb1, lOrb2);
   ovl.reformat (lOrb1, lOrb2);

   int counter = 0;
   if (isSameSpecies)  {
      hubbardU.resize (lOrb1, 0.);
      epsAtom.resize  (lOrb1, 0.);
      for (l = 0; l <= lMax1; l++)  {
         for (iOrb = 0; iOrb < (2*l+1); iOrb++)  {
            hubbardU(counter) = hubL(l);  // extending over m_l quantum number 
            epsAtom (counter) = epsL(l);  // i.e. over all orbitals
            counter += 1;
         }
      }
   }

   // --- default values
   dCos.set(0.); dSin.set(0.);
   dist = 0.;    dist2 = 0.;
   ham.set(0.);  ovl.set(0.);

   /* // TEST for future SX file format
   SxParser parser;
   SxParser::Table slaterKoster = parser.read ("sk.sx");

   SxSymbolTable *data, *elements = NULL;  
   int l1, l2, coupling;
   try  {
      for ( data  = slaterKoster->getGroup ("data");
            data != NULL;
            data  = data->nextSibling ("data"))
      {
         l1       = data->get("l1")->toInt();
         l2       = data->get("l2")->toInt();
         coupling = data->get("coupling")->toInt();
         for ( elements  = data->getGroup ("elements");
               elements != NULL;
               elements  = elements->nextSibling ("elements"))
         {
            //r << (elements->get("r")->toReal());
            //H << (elements->get("H")->toReal());
            //S << (elements->get("S")->toReal());
         }
      }
   }  catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   */

}

double SxTBAtomAtom::getERepulsive (double distIn) const
{
   int    i, j;
   double repEn = 0; 
   double deltaX, fac, deltaAB, deltaBC, fA, fB, fC;
   bool   found;
   double distance = distIn;

   if (withoutERepulsive)  return 0;

   if (repType == "functionFit" || "Spline" )  {
      if (distance < 0.01 || distance >= dCut)
         return 0.;
      else if (distance < interval(0, 0) )
         return ( exp(-func(0)*distance + func(1)) + func(2) );
      else  { 
         found = false;
         for (i = 0; i < nLines; i++)  {
            if (distance >= interval(i,0) && distance < interval(i,1))  {
               found  = true;
               deltaX = distance - interval(i,0);
               fac    = 1;
               repEn  = 0.;
               if (i < nLines-1)  {
                  for (j = 0; j < 4; j++)  {
                     repEn += repulsiveCoef1(i,j)* fac; 
                     fac   *= deltaX;
                  }
               }  else  {
                  for (j = 0; j < 6; j++)  {
                     repEn += repulsiveCoef2(j)* fac; 
                     fac   *= deltaX;
                  }
               }
               break;  
            }
         }
         if (!found) {
            cout << " in getERepulsive : " ;
            cout << " I can't find your interval !!" << endl;
            SX_EXIT;
         }
         return repEn;
      }

   }  else if (repType == "dataFit")  {

      if (distance < repData(0)(0))
         // TODO: this will cause a discontinuity!, but in a not interesting region
         // should be modefied soon
         return 10; // somme large value !
      else if (distance > repData(nLines-3)(0)) // last few repData should = 0 !
         return 0;
      else  {
         found = false;
         for (i = 0; i < (nLines-2); i++)  {
            if (distance >= repData(i)(0) && distance < repData(i+1)(0))  {
               found  = true;
               deltaX = distance - repData(i)(0);
               deltaAB = repData(i+1)(0) - repData(i)(0);
               deltaBC = repData(i+2)(0) - repData(i+1)(0);
               fA = repData(i)(1);
               fB = repData(i+1)(1);
               fC = repData(i+2)(1);
               // Taylor series expansion
               repEn = fA + (fB - fA)*deltaX/deltaAB 
                  + (fC + fA - 2.*fB)*(deltaX/deltaAB)*(deltaX/deltaBC-1.)/2.0; 
               break;  
            }
         }
         if (!found) {
            cout << " in getERepulsive : " ;
            cout << " I can't find your interval !!" << endl;
            SX_EXIT;
         }
         return repEn;
      }

   }  else  { 
      cout << "\n unknown repulsive energy type = " << repType << endl;
      SX_QUIT;
      return 0;

   }

}

double SxTBAtomAtom::getRepulsiveEGrad (double distIn) const 
{
   int    i, j;
   double repGrad = 0;
   double deltaX, fac;
   bool found;
   double distance = distIn;

   if (withoutERepulsive)  return 0;

   if (repType == "functionFit" || "Spline")  {
      if (distance < 0.01 || distance >= dCut)
         return 0.;
      else if (distance < interval(0, 0) )
         return ( -func(0)*exp(-func(0)*distance + func(1)));
      else  { 
         found = false;
         for (i = 0; i < nLines; i++)  {
            if (distance >= interval(i,0) && distance < interval(i,1))  {
               found   = true;
               deltaX  = distance - interval(i,0);
               fac     = 1;
               repGrad = 0.;
               if (i < nLines-1)  {
                  for (j = 1; j < 4; j++)  {
                     repGrad += j*repulsiveCoef1(i,j)* fac; 
                     fac     *= deltaX;
                  }
               }  else  {
                  for (j = 1; j < 6; j++)  {
                     repGrad += j*repulsiveCoef2(j)* fac; 
                     fac     *= deltaX;
                  }
               }
               break;  
            }
         }
         if (!found) {
            cout << " in getRepulsiveEGrad : " ;
            cout << " I can't find your interval !!" << endl;
            SX_EXIT;
         }
         return repGrad;
      }
   }  else if (repType == "dataFit")  {
      // numerical derivative in this case 
      repGrad = ( getERepulsive (distance+0.01) 
                - getERepulsive (distance-0.01) )/0.02; 
      return repGrad;
   }  else  { 
      cout << "\n unknown repulsive energy type = " << repType << endl;
      SX_QUIT;
      return 0;
   }

}
   
int SxTBAtomAtom::getNOrb1 () const
{
   return lOrb1; 
} 
   
int SxTBAtomAtom::getNOrb2 () const
{
   return lOrb2; 
} 

int SxTBAtomAtom::getNElect () const
{
   SX_CHECK (isSameSpecies, isSameSpecies);
   return nValElect; 
} 

int SxTBAtomAtom::getlOrb (int lMax) const
{
   int lOrb, l, n;
   lOrb = 0;
   for (l = 0; l <= lMax; l++)  {
      n = 2*l + 1;
      lOrb += n;
   }
   return lOrb; 
} 

void SxTBAtomAtom::setDistVec (const SxVector3<Double> &dVecIn)
{
   int    i;

   dist = dVecIn.norm();
   dCos = dVecIn/dist;
   for (i = 0; i < 3; i++)
      dSin(i) = sqrt(1 - dCos(i)*dCos(i)); 
   
   index    = (int)((dist + 1e-8)/delta) - 1;//files start from dist=0.02 not 0.
   fraction = (dist - (index + 1)*delta)/delta;

   // --- kind of validation of the distance
   if (dist < 0.4)  {   // in Bohr. This value depends on SK files !
      cout << "\n two atoms are very close !\n";
      SX_QUIT;
   }
}
   
void SxTBAtomAtom::setElements ()
{
   nCol = (int)slkoHam.nCols ();
   vCoeff.resize (nCol);
   sCoeff.resize (nCol);

   if (index + 2 >= nPoints)  {
      //cout << "Distance between 2 atoms > largest value in the corresponding "
      //     << "SK file,\nthis will set the corresponding elements of "  
      //     << "hamiltonian and overlap to 0.\n\n";
      vCoeff.set(0.);
      sCoeff.set(0.);
      ham.set(0.);
      ovl.set(0.);
      return;
   }

   vCoeff = spline (slkoHam, index, fraction);
   sCoeff = spline (slkoOvl, index, fraction);

   ham = getSubMatrix (vCoeff);
   ovl = getSubMatrix (sCoeff);

}
 
void SxTBAtomAtom::setDiagElements ()
{
   SX_CHECK (isSameSpecies);
   
   int il, jl;
   for (il = 0; il < lOrb1; il++)  {
      for (jl = 0; jl < lOrb2; jl++)  {
         if (il == jl)  {
            ham (il, jl) = epsAtom (il);
            ovl (il, jl) = 1.;
         }  else  {
            ham (il, jl) = 0.;
            ovl (il, jl) = 0.;
         }
      }
   }

}

SxVector<Double> SxTBAtomAtom::spline (const SxMatrix<Double> &slko, 
                                       int ind, double frac) const
{
   //SX_CHECK  (frac >= -1e-12 && frac <= 1., frac);
   int    iCol;
   double f0, f1, f2;
   SxVector<Double> sk(nCol);

   for (iCol = 0; iCol < nCol; iCol++)  {
      f0       = slko (ind, iCol);
      f1       = slko (ind + 1, iCol);
      f2       = slko (ind + 2, iCol);
      sk(iCol) = f0 + (f1-f0)*frac + (f2+f0-2.0*f1)*frac*(frac-1.)/2.; 
   }
   return sk;
}

SxMatrix<Double> SxTBAtomAtom::getSubMatrix 
                              (const SxVector<Double> &coeff) const
{
   SxMatrix<Double> inter (lOrb1, lOrb2);

   SxArray<SxMatrix<Double> >  matElement(3);
   // matElement(type(sigma=0, pi=1, delta=2..))(1st orbital, 2nd orbital)
   // e.g. matElement (0)(0,1) = Vsp-sigma, matElement (1)(2,2) = Vdd-pi

   for (int i = 0; i < matElement.getSize(); i++)  {
      matElement(i).reformat (lMax1+1, lMax2+1); 
      matElement(i).set(0);
   }

   // --- store
   int l1, l2, i;
   for (l1 = 0; l1 <= lMax1; l1++)  {
      for (l2 = 0; l2 <= lMax2; l2++)  {
         if      (l1 == 0 && l2 == 0) matElement(0)(0,0) = coeff(0);
         else if (l1 == 0 && l2 == 1) matElement(0)(0,1) = coeff(1);
         else if (l1 == 1 && l2 == 0) matElement(0)(1,0) = coeff(2);
         else if (l1 == 1 && l2 == 1)  {
            matElement(0)(1,1) = coeff(4);
            matElement(1)(1,1) = coeff(3);
         }  
         else if (l1 == 0 && l2 == 2) matElement(0)(0,2) = coeff(5); // s d interactions  
         else if (l1 == 2 && l2 == 0) matElement(0)(2,0) = coeff(6); // d s interactions
         else if (l1 == 1 && l2 == 2) {  // p d interactions
            matElement(0)(1,2) = coeff(9);
            matElement(1)(1,2) = coeff(7);
         }  
         else if (l1 == 2 && l2 == 1)  { // d p interactions
            matElement(0)(2,1) = coeff(10);
            matElement(1)(2,1) = coeff(8);
         }  
         else if (l1 == 2 && l2 == 2)  { // d d interactions
            matElement(0)(2,2) = coeff(13);
            matElement(1)(2,2) = coeff(12);
            matElement(2)(2,2) = coeff(11);
         }  
         else  {
            cout << "\n invalid value(s) of lMax\n";
            SX_EXIT;
         }
      }	
   }
   // now we can transpose them if necessary
   if (isTransposed)
      for (i = 0; i < matElement.getSize(); i++)
         matElement(i) = matElement(i).transpose();

   double M = dCos(0);
   double N = dCos(1);
   double L = dCos(2);

   // --------------------------------------------------
   // --- see Ref.: Phys. Rev. B 69, 233101 (2004)
   // --------------------------------------------------
   for (l1 = 0; l1 <= lMax1; l1++)  {
      for (l2 = 0; l2 <= lMax2; l2++)  {
         
         if (l1 == 0 && l2 == 0)  {
            inter(0, 0) = matElement(0)(0,0);
         }  else if (l1 == 0 && l2 == 1)  {
            inter(0, 1) = matElement(0)(0,1) * dCos(0);
            inter(0, 2) = matElement(0)(0,1) * dCos(1);
            inter(0, 3) = matElement(0)(0,1) * dCos(2);
         }  else if (l1 == 1 && l2 == 0)  {
            inter(1, 0) = -matElement(0)(1,0) * dCos(0);
            inter(2, 0) = -matElement(0)(1,0) * dCos(1);
            inter(3, 0) = -matElement(0)(1,0) * dCos(2);
         }  else if (l1 == 1 && l2 == 1)  {
            inter(1, 1) = matElement(0)(1,1) * dCos(0)*dCos(0) 
                        + matElement(1)(1,1) * dSin(0)*dSin(0);
            inter(1, 2) = matElement(0)(1,1) * dCos(0)*dCos(1)
                        + matElement(1)(1,1) *-dCos(0)*dCos(1);
            inter(1, 3) = matElement(0)(1,1) * dCos(0)*dCos(2)
                        + matElement(1)(1,1) *-dCos(0)*dCos(2);
            
            inter(2, 1) = matElement(0)(1,1) * dCos(0)*dCos(1)
                        + matElement(1)(1,1) *-dCos(1)*dCos(0);
            inter(2, 2) = matElement(0)(1,1) * dCos(1)*dCos(1)
                        + matElement(1)(1,1) * dSin(1)*dSin(1);
            inter(2, 3) = matElement(0)(1,1) * dCos(1)*dCos(2)
                        + matElement(1)(1,1) *-dCos(1)*dCos(2);
            
            inter(3, 1) = matElement(0)(1,1) * dCos(0)*dCos(2)
                        + matElement(1)(1,1) *-dCos(2)*dCos(0);
            inter(3, 2) = matElement(0)(1,1) * dCos(1)*dCos(2) 
                        + matElement(1)(1,1) *-dCos(2)*dCos(1);
            inter(3, 3) = matElement(0)(1,1) * dCos(2)*dCos(2) 
                        + matElement(1)(1,1) * dSin(2)*dSin(2);

         }  else if (l1 == 0 && l2 == 2)  {
            inter(0, 4) = sqrt (3.)*L*M*matElement(0)(0,2);
            inter(0, 5) = sqrt (3.)*M*N*matElement(0)(0,2);
            inter(0, 6) = ((-1 + 3*pow (N,2))*matElement(0)(0,2))/2.;
            inter(0, 7) = sqrt (3.)*L*N*matElement(0)(0,2);
            inter(0, 8) = (sqrt (3.)*(L - M)*(L + M)*matElement(0)(0,2))/2.;

         }  else if (l1 == 2 && l2 == 0)  {
            inter(4, 0) = sqrt (3.)*L*M*matElement(0)(2,0);
            inter(5, 0) = sqrt (3.)*M*N*matElement(0)(2,0);
            inter(6, 0) = ((-1 + 3*pow (N,2))*matElement(0)(2,0))/2.;
            inter(7, 0) = sqrt (3.)*L*N*matElement(0)(2,0);
            inter(8, 0) = (sqrt (3.)*(L - M)*(L + M)*matElement(0)(2,0))/2.;

         }  else if (l1 == 1 && l2 == 2)  {

            inter(1, 4) = sqrt (3.)*L*pow (M,2)*matElement(0)(1,2) - (L*(pow (L,2) 
                     + pow (M,2)*(-1 + 2*pow (N,2)))*matElement(1)(1,2))/(-1 + pow (N,2));
            inter(1, 5) = (N*(sqrt (3.)*pow (M,2)*(-1 + pow (N,2))*matElement(0)(1,2) 
                     + (-pow(L,2) + pow (M,2)*(1 - 2*pow (N,2)))*matElement(1)(1,2)))/(-1 + pow (N,2));
            inter(1, 6) = (M*((-1 + 3*pow (N,2))*matElement(0)(1,2) 
                     - 2*sqrt (3.)*pow (N,2)*matElement(1)(1,2)))/2.;
            inter(1, 7) = L*M*N*(sqrt (3.)*matElement(0)(1,2) - 2*matElement(1)(1,2));
            inter(1, 8) = (M*(pow (M,2)*(-(sqrt (3.)*(-1 + pow (N,2))*matElement(0)(1,2)) 
                        + 2*pow (N,2)*matElement(1)(1,2)) + pow (L,2)*(sqrt (3.)*(-1 + pow (N,2))*matElement(0)(1,2) 
                        - 2*(-2 + pow (N,2))*matElement(1)(1,2))))/(2.*(-1 + pow (N,2)));

            inter(2, 4) = L*M*N*(sqrt (3.)*matElement(0)(1,2) - 2*matElement(1)(1,2));
            inter(2, 5) = M*(pow (N,2)*(sqrt (3.)*matElement(0)(1,2) 
                     - 2*matElement(1)(1,2)) + matElement(1)(1,2));
            inter(2, 6) = (N*((-1 + 3*pow (N,2))*matElement(0)(1,2) 
                     - 2*sqrt (3.)*(-1 + pow (N,2))*matElement(1)(1,2)))/2.;
            inter(2, 7) = L*(pow (N,2)*(sqrt (3.)*matElement(0)(1,2) 
                     - 2*matElement(1)(1,2)) + matElement(1)(1,2));
            inter(2, 8) = ((L - M)*(L + M)*N*(sqrt (3.)*matElement(0)(1,2) - 2*matElement(1)(1,2)))/2.;

            inter(3, 4) = sqrt (3.)*pow (L,2)*M*matElement(0)(1,2) 
               - (M*(pow (M,2) + pow (L,2)*(-1 + 2*pow (N,2)))*matElement(1)(1,2))/(-1 + pow (N,2));
            inter(3, 5) = L*M*N*(sqrt (3.)*matElement(0)(1,2) - 2*matElement(1)(1,2));
            inter(3, 6) = (L*((-1 + 3*pow (N,2))*matElement(0)(1,2) 
                     - 2*sqrt (3.)*pow (N,2)*matElement(1)(1,2)))/2.;
            inter(3, 7) = (N*(sqrt (3.)*pow (L,2)*(-1 + pow (N,2))*matElement(0)(1,2) 
                     + (-pow(M,2) + pow (L,2)*(1 - 2*pow (N,2)))*matElement(1)(1,2)))/(-1 + pow (N,2));
            inter(3, 8) = (L*(pow (L,2)*(sqrt (3.)*(-1 + pow (N,2))*matElement(0)(1,2) - 2*pow (N,2)*matElement(1)(1,2)) 
                     + pow (M,2)*(-(sqrt (3.)*(-1 + pow (N,2))*matElement(0)(1,2)) 
                        + 2*(-2 + pow (N,2))*matElement(1)(1,2))))/(2.*(-1 + pow (N,2)));

         }  else if (l1 == 2 && l2 == 1)  {

            inter(4, 1) = -(sqrt (3.)*L*pow (M,2)*matElement(0)(2,1)) 
               + (L*(pow (L,2) + pow (M,2)*(-1 + 2*pow (N,2)))*matElement(1)(2,1))/(-1 + pow (N,2));
            inter(4, 2) = L*M*N*(-(sqrt (3.)*matElement(0)(2,1)) + 2*matElement(1)(2,1));
            inter(4, 3) = -(sqrt (3.)*pow (L,2)*M*matElement(0)(2,1)) 
               + (M*(pow (M,2) + pow (L,2)*(-1 + 2*pow (N,2)))*matElement(1)(2,1))/(-1 + pow (N,2));

            inter(5, 1) = (N*(-(sqrt (3.)*pow (M,2)*(-1 + pow (N,2))*matElement(0)(2,1)) 
                     + (pow (L,2) + pow (M,2)*(-1 + 2*pow (N,2)))*matElement(1)(2,1)))/(-1 + pow (N,2));
            inter(5, 2) = -(M*(pow (N,2)*(sqrt (3.)*matElement(0)(2,1) 
                        - 2*matElement(1)(2,1)) + matElement(1)(2,1)));
            inter(5, 3) = L*M*N*(-(sqrt (3.)*matElement(0)(2,1)) + 2*matElement(1)(2,1));

            inter(6, 1) = (M*((1 - 3*pow (N,2))*matElement(0)(2,1) 
                     + 2*sqrt (3.)*pow (N,2)*matElement(1)(2,1)))/2.;
            inter(6, 2) = (N*((1 - 3*pow (N,2))*matElement(0)(2,1) 
                     + 2*sqrt (3.)*(-1 + pow (N,2))*matElement(1)(2,1)))/2.;
            inter(6, 3) = (L*((1 - 3*pow (N,2))*matElement(0)(2,1) 
                     + 2*sqrt (3.)*pow (N,2)*matElement(1)(2,1)))/2.;

            inter(7, 1) = L*M*N*(-(sqrt (3.)*matElement(0)(2,1)) + 2*matElement(1)(2,1));
            inter(7, 2) = -(L*(pow (N,2)*(sqrt (3.)*matElement(0)(2,1) 
                        - 2*matElement(1)(2,1)) + matElement(1)(2,1)));
            inter(7, 3) = (N*(-(sqrt (3.)*pow (L,2)*(-1 + pow (N,2))*matElement(0)(2,1)) 
                     + (pow (M,2) + pow (L,2)*(-1 + 2*pow (N,2)))*matElement(1)(2,1)))/(-1 + pow (N,2));

            inter(8, 1) = (M*(pow (M,2)*(sqrt (3.)*(-1 + pow (N,2))*matElement(0)(2,1) 
                        - 2*pow (N,2)*matElement(1)(2,1)) + pow (L,2)*(-(sqrt (3.)*(-1 + pow (N,2))*matElement(0)(2,1)) 
                        + 2*(-2 + pow (N,2))*matElement(1)(2,1))))/(2.*(-1 + pow (N,2)));
            inter(8, 2) = -((L - M)*(L + M)*N*(sqrt (3.)*matElement(0)(2,1) - 2*matElement(1)(2,1)))/2.;
            inter(8, 3) = (L*(pow (L,2)*(-(sqrt (3.)*(-1 + pow (N,2))*matElement(0)(2,1)) 
                        + 2*pow (N,2)*matElement(1)(2,1)) + pow (M,2)*(sqrt (3.)*(-1 + pow (N,2))*matElement(0)(2,1) 
                        - 2*(-2 + pow (N,2))*matElement(1)(2,1))))/(2.*(-1 + pow (N,2)));

         }  else if (l1 == 2 && l2 == 2)  {

            inter(4, 4) = ((-1 + pow (N,2))*(3*pow (L,2)*pow (M,2)*(-1 + pow (N,2))*matElement(0)(2,2) 
                     - (pow (L,4) + pow (M,4) + 2*pow (L,2)*pow (M,2)*(-1 + 2*pow (N,2)))*matElement(1)(2,2)) 
                  + (pow (M,2) + pow (L,2)*pow (N,2))*(pow (L,2) + pow (M,2)*pow (N,2))*matElement(2)(2,2))/pow(-1 + pow (N,2),2);
            inter(4, 5) = (L*N*(pow (L,2)*(-matElement(1)(2,2) + matElement(2)(2,2)) 
                     + pow (M,2)*(3*(-1 + pow (N,2))*matElement(0)(2,2) + (3 - 4*pow (N,2))*matElement(1)(2,2) 
                        + pow (N,2)*matElement(2)(2,2))))/(-1 + pow (N,2));
            inter(4, 6) = (sqrt (3.)*L*M*((-1 + 3*pow (N,2))*matElement(0)(2,2) 
                     + matElement(2)(2,2) + pow (N,2)*(-4*matElement(1)(2,2) + matElement(2)(2,2))))/2.;
            inter(4, 7) = (M*N*(pow (M,2)*(-matElement(1)(2,2) + matElement(2)(2,2)) 
                     + pow (L,2)*(3*(-1 + pow (N,2))*matElement(0)(2,2) 
                        + (3 - 4*pow (N,2))*matElement(1)(2,2) + pow (N,2)*matElement(2)(2,2))))/(-1 + pow (N,2));
            inter(4, 8) = (L*(L - M)*M*(L + M)*(3*matElement(0)(2,2) - 4*matElement(1)(2,2) + matElement(2)(2,2)))/2.;

            inter(5, 4) = (L*N*(pow (L,2)*(-matElement(1)(2,2) + matElement(2)(2,2)) 
                     + pow (M,2)*(3*(-1 + pow (N,2))*matElement(0)(2,2) + (3 - 4*pow (N,2))*matElement(1)(2,2) 
                        + pow (N,2)*matElement(2)(2,2))))/(-1 + pow (N,2));
            inter(5, 5) = (pow (L,2)*(-(pow (N,2)*matElement(1)(2,2)) + (-1 + pow (N,2))*matElement(2)(2,2)) 
                  + pow (M,2)*(-matElement(1)(2,2) + pow (N,2)*(-1 + pow (N,2))*(3*matElement(0)(2,2) - 4*matElement(1)(2,2) 
                        + matElement(2)(2,2))))/(-1 + pow (N,2));
            inter(5, 6) = (sqrt (3.)*M*N*((-1 + 3*pow (N,2))*matElement(0)(2,2) 
                     + (2 - 4*pow (N,2))*matElement(1)(2,2) + (-1 + pow (N,2))*matElement(2)(2,2)))/2.;
            inter(5, 7) = L*M*(matElement(1)(2,2) - matElement(2)(2,2) + pow (N,2)*(3*matElement(0)(2,2) 
                     - 4*matElement(1)(2,2) + matElement(2)(2,2)));
            inter(5, 8) = (M*N*(pow (L,2)*(3*(-1 + pow (N,2))*matElement(0)(2,2) 
                        + (6 - 4*pow (N,2))*matElement(1)(2,2) + (-3 + pow (N,2))*matElement(2)(2,2)) 
                     - pow (M,2)*(3*(-1 + pow (N,2))*matElement(0)(2,2) + (2 - 4*pow (N,2))*matElement(1)(2,2) 
                        + (1 + pow (N,2))*matElement(2)(2,2))))/(2.*(-1 + pow (N,2)));


            inter(6, 4) = (sqrt (3.)*L*M*((-1 + 3*pow (N,2))*matElement(0)(2,2) 
                     + matElement(2)(2,2) + pow (N,2)*(-4*matElement(1)(2,2) + matElement(2)(2,2))))/2.;
            inter(6, 5) = (sqrt (3.)*M*N*((-1 + 3*pow (N,2))*matElement(0)(2,2) 
                     + (2 - 4*pow (N,2))*matElement(1)(2,2) + (-1 + pow (N,2))*matElement(2)(2,2)))/2.;
            inter(6, 6) = (pow (1 - 3*pow (N,2),2)*matElement(0)(2,2) + 3*(-1 + pow (N,2))*(-4*pow (N,2)*matElement(1)(2,2) 
                     + (-1 + pow (N,2))*matElement(2)(2,2)))/4.;
            inter(6, 7) = (sqrt (3.)*L*N*((-1 + 3*pow (N,2))*matElement(0)(2,2) 
                     + (2 - 4*pow (N,2))*matElement(1)(2,2) + (-1 + pow (N,2))*matElement(2)(2,2)))/2.;
            inter(6, 8) = (sqrt (3.)*(L - M)*(L + M)*((-1 + 3*pow (N,2))*matElement(0)(2,2) 
                     + matElement(2)(2,2) + pow (N,2)*(-4*matElement(1)(2,2) + matElement(2)(2,2))))/4.;

            inter(7, 4) = (M*N*(pow (M,2)*(-matElement(1)(2,2) + matElement(2)(2,2)) 
                     + pow (L,2)*(3*(-1 + pow (N,2))*matElement(0)(2,2) + (3 - 4*pow (N,2))*matElement(1)(2,2) 
                        + pow (N,2)*matElement(2)(2,2))))/(-1 + pow (N,2));
            inter(7, 5) = L*M*(matElement(1)(2,2) - matElement(2)(2,2) + pow (N,2)*(3*matElement(0)(2,2) 
                     - 4*matElement(1)(2,2) + matElement(2)(2,2)));
            inter(7, 6) = (sqrt (3.)*L*N*((-1 + 3*pow (N,2))*matElement(0)(2,2) 
                     + (2 - 4*pow (N,2))*matElement(1)(2,2) + (-1 + pow (N,2))*matElement(2)(2,2)))/2.;
            inter(7, 7) = (pow (M,2)*(-(pow (N,2)*matElement(1)(2,2)) + (-1 + pow (N,2))*matElement(2)(2,2)) 
                  + pow (L,2)*(-matElement(1)(2,2) + pow (N,2)*(-1 + pow (N,2))*(3*matElement(0)(2,2) 
                        - 4*matElement(1)(2,2) + matElement(2)(2,2))))/(-1 + pow (N,2));
            inter(7, 8) = (L*N*(pow (M,2)*(-3*(-1 + pow (N,2))*matElement(0)(2,2) 
                        + (-6 + 4*pow (N,2))*matElement(1)(2,2) - (-3 + pow (N,2))*matElement(2)(2,2)) 
                     + pow (L,2)*(3*(-1 + pow (N,2))*matElement(0)(2,2) + (2 - 4*pow (N,2))*matElement(1)(2,2) 
                        + (1 + pow (N,2))*matElement(2)(2,2))))/(2.*(-1 + pow (N,2)));

            inter(8, 4) = (L*(L - M)*M*(L + M)*(3*matElement(0)(2,2) 
                     - 4*matElement(1)(2,2) + matElement(2)(2,2)))/2.;
            inter(8, 5) = (M*N*(pow (L,2)*(3*(-1 + pow (N,2))*matElement(0)(2,2) 
                        + (6 - 4*pow (N,2))*matElement(1)(2,2) + (-3 + pow (N,2))*matElement(2)(2,2)) 
                     - pow (M,2)*(3*(-1 + pow (N,2))*matElement(0)(2,2) + (2 - 4*pow (N,2))*matElement(1)(2,2) 
                        + (1 + pow (N,2))*matElement(2)(2,2))))/(2.*(-1 + pow (N,2)));
            inter(8, 6) = (sqrt (3.)*(L - M)*(L + M)*((-1 + 3*pow (N,2))*matElement(0)(2,2) 
                     + matElement(2)(2,2) + pow (N,2)*(-4*matElement(1)(2,2) + matElement(2)(2,2))))/4.;
            inter(8, 7) = (L*N*(pow (M,2)*(-3*(-1 + pow (N,2))*matElement(0)(2,2) 
                        + (-6 + 4*pow (N,2))*matElement(1)(2,2) - (-3 + pow (N,2))*matElement(2)(2,2)) 
                     + pow (L,2)*(3*(-1 + pow (N,2))*matElement(0)(2,2) + (2 - 4*pow (N,2))*matElement(1)(2,2) 
                        + (1 + pow (N,2))*matElement(2)(2,2))))/(2.*(-1 + pow (N,2)));
            inter(8, 8) = (3*pow (pow (L,2) - pow (M,2),2)*matElement(0)(2,2) - (4*(4*pow (L,2)*pow (M,2) 
                        + pow (pow (L,2) - pow (M,2),2)*pow (N,2))*matElement(1)(2,2))/(-1 + pow (N,2)) 
                  + ((pow (L,4)*pow (1 + pow (N,2),2) + pow (M,4)*pow (1 + pow (N,2),2) 
                        - 2*pow (L,2)*pow (M,2)*(1 - 6*pow (N,2) 
                           + pow (N,4)))*matElement(2)(2,2))/pow(-1 + pow (N,2),2))/4.;

         }  else  {
            cout << endl;
            cout << "\n invalid value(s) of lMax\n";
            SX_EXIT;
         }
         
      }
   }

   return inter;
}  

double SxTBAtomAtom::getHam (int iOrb, int jOrb) const
{
   SX_CHECK (iOrb <= lOrb1 && iOrb >= 0, iOrb, lOrb1);
   SX_CHECK (jOrb <= lOrb2 && jOrb >= 0, jOrb, lOrb2);
   return ham(iOrb, jOrb);
}    

double SxTBAtomAtom::getOvl (int iOrb, int jOrb) const
{
   SX_CHECK (iOrb <= lOrb1 && iOrb >= 0, iOrb, lOrb1);
   SX_CHECK (jOrb <= lOrb2 && jOrb >= 0, jOrb, lOrb2);
   return ovl(iOrb, jOrb);
}    

SxMatrix<Double> SxTBAtomAtom::getAtomAtomHam () const
{
   return ham;
}    

SxMatrix<Double> SxTBAtomAtom::getAtomAtomOvl () const
{
   return ovl;
}

double SxTBAtomAtom::getNeighborsCutoff () const
{
   return cutoff;
}

double SxTBAtomAtom::getHubbardU (int iOrb) const
{
   SX_CHECK  (isSameSpecies, isSameSpecies);
   SX_CHECK (iOrb <= lOrb1 && iOrb >= 0, iOrb, lOrb1);

   return hubbardU(iOrb);
}

