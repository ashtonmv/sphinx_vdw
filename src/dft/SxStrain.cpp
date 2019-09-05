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


/* Implementation of continuum theory strain field calculation as done in:

ref1:
Craig Pryor: Eight-band calculations of strained InAs/GaAs quantum dots
compared with one-, four- and six-band approximations,
Phys. Rev. \textbf{B}, Vol. 57, 12 (1998)
*/

#include <SxStrain.h>
#include <SxProjector.h>
#include <SxGBasis.h>
#include <SxFSAction.h>

SxStrain::SxStrain (const SxPWSet &waves, const SxRBasis* rPtr)
   : SxHamiltonian (),
     rBasisPtr (rPtr),
     G1(waves(0,0,0).getBasis<SxGBasis> (),1)
// SxStrain
{
   SX_CHECK (rBasisPtr);
   cout << "\n Strain field calculation\n";
   gVecs.resize(1);
   step = 0;
   Diag.resize(1);
   SX_CHECK (G1.structPtr);
   R1.set (G1.fft3d(0).mesh, G1.structPtr->cell);
   R1.setNComp(1);
   G1.registerRBasis(R1);
   R1.registerGBasis(G1);
}

SxStrain::~SxStrain ()
{

}

void SxStrain::printEnergies () const
{
}

double SxStrain::interpolate (double y1, double y2, double pIn, double bow)
{
   return ( y2 - pIn * (y2 - y1) - y2 * bow * pIn * (1 - pIn));
}

const SxSymbolTable *
SxStrain::getHamiltonianGroup (const SxSymbolTable *table)
{
   SX_CHECK (table);
   SxSymbolTable *hGroup = NULL;
   try  { hGroup = table->getGroup("StrainField"); }
   catch (const SxException&) { /*EXIT;*/ }
   return hGroup;
}


PrecEnergy SxStrain::getETrial ()
{
   return 0.;
}

/* Seems to be unused, CF 2017-06-02
double SxStrain::mixing (double x1, double y1, double x2, double y2, double x)
{
   double y2I = interpolate(y2, y1, p, 0.);
   double m;
   m = (y1 + (x + x2) * (y2I - y1) / (x2 - x1));
   return m;
}
*/

SxMeshR SxStrain::parameters(int iParam)
{
   SxMeshR param;
   int iMat;
   param = matParam(0, iParam) * materials(0);
   for (iMat = 1; iMat < nMat; iMat++)
      param += matParam(iMat, iParam) * materials(iMat);
   return param;
}

void SxStrain::read (const SxSymbolTable *table)
{
   generalized = false;
   generalizedLat = false;
   nonzeroC.resize(6);
   for (int i = 0; i < 6; i++)  {
      nonzeroC(i).resize(6);
      for (int j = 0; j < 6; j++)
         nonzeroC(i)(j) = false;
   }


   // --- read the material parameter files ---

   SxList<SxString> matNames; // names of involved materials
   SxList<SxMatrix<Double> > cValues; // list of elastic constant tensors
   SxList<SxVector<Double> > latValues; // list of lattice constant tensors

   SxList<double> paramVals;
   SxMap<SxString,int> paramMap; // each parameter gets an integer index

   nMat = 0;
   containsVacuum = false;
   readExternalStrains = false;
   nParam = 0;
	SxSymbolTable *paramSet = NULL, *param = NULL, *materialMap = NULL;

   // --- loop over all materials parameter sets
   SxMap<SxString, int>::Iterator it;
   vacMat = -1;
   for (paramSet  = table->getGroup("parameterSet");
        paramSet != NULL;
        paramSet  = paramSet->nextSibling ("parameterSet"))
   {
      // --- first material defines parameter map name->integer
      if (nMat == 0)  {
         for (param  = paramSet->getGroup("parameter");
               param != NULL;
               param  = param->nextSibling ("parameter"))
         {
            SxString name = param->get("name")->toString ();
            paramMap(name) = nParam;
            paramNames << name;
            nParam++;
         }
      }

      SxString material = paramSet->get("material")->toString ();
      matNames << material;
      SxVector<Double> tmpLat;
      tmpLat.resize(3);
      if (paramSet->contains("latticeConstants"))  {
         generalizedLat = true;
         tmpLat = SxMatrix<Double>(paramSet->get("latticeConstants")->toList());
         cout << "Lattice constants of material " + material + ":" << tmpLat << endl;
         latValues << tmpLat;
      }
      if (paramSet->contains("elasticTensor"))  {
         generalized = true;
         SxVector<Double> tmp;
         tmp.resize(36);
         tmp = SxMatrix<Double>(paramSet->get("elasticTensor")->toList());
         SxMatrix<Double> tmpMatrix(6,6);
         for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++) {
               tmpMatrix(i,j) = tmp(i*6 + j);
               if (tmpMatrix(i,j) != 0) {nonzeroC(i)(j) = true;}
               // --- issue warning if Cji != Cij
               if ((i > j) && (tmpMatrix(i,j) != tmpMatrix(j,i)))  {
                  cout << "WARNING: " << material << "'s C" << i+1 << j+1
                       << " != C" << j+1 << i+1 << "!" << endl;
               }
            }
         cout << "elastic tensor of material " + material + ":" << endl;
         cout << tmpMatrix << endl;
         cValues << tmpMatrix;
      }

      // --- browse parameter set for all parameter names
      for (it = paramMap.begin(); it != paramMap.end(); it++)  {
         if (paramSet->contains("vacuum"))  {
            vacMat = nMat;
            containsVacuum = true;
         }
         bool found = false; // parameter found in material parameter group?
         param = paramSet->getGroup("parameter");
         while ((!found) && (param != NULL)) {
            SxString name = param->get("name")->toString();
            if (it.getKey() == name)
            {
               double val = param->get("value")->toReal();
               paramVals << val;
               found = true;
            }
            param = param->nextSibling ("parameter");
         }
         if (!found)  {
            cout << "Error: Parameter " << it.getKey() <<
               " was not found in the " << material << " parameter file."
               << endl << "EXITING." << endl;
            SX_EXIT;
         }
      }
      nMat++;
   }

   matParam = SxMatrix<Double> (nMat, nParam, paramVals);
   int iMat, iParam, idx;

   cout << "####### Input materials and parameters ########" << endl;
   cout << "Mat.\t  |  ";
   for (iMat = 0; iMat < nMat; iMat++)
      cout << matNames(iMat) << "\t| ";
   cout << endl;
   for (iParam = 0; iParam < nParam; iParam++)  {
      cout << paramNames(iParam) << " \t| ";
      for (iMat = 0; iMat < nMat; iMat++)  {
         cout << matParam(iMat, iParam) << " \t| ";
      }
      cout << endl;
   }

   // --- read the material map file ---
   cout << "Load map... ";
   SxString inputASCII = "";
   SxString inputBinary = "";
   mesh3d = rBasisPtr->getMesh();
   SxMesh3D mesh = mesh3d;
   cout << mesh3d << endl;

   int xMax = mesh3d(0);
   int yMax = mesh3d(1);
   int zMax = mesh3d(2);
   int x, y, z;
   float value;
   rSize = xMax * yMax * zMax;
   one.resize(rSize);
   zero.resize(rSize);
   one.set (1.);
   zero.set (0.);

   materialMap  = table->getGroup("materialMap");

   short dataType = -1;
   if (materialMap->contains("asciiFile"))  {
      inputASCII = materialMap->get("asciiFile")->toString();
      dataType = 1;
   } else if (materialMap->contains("binaryFile"))  {
      inputBinary = materialMap->get("binaryFile")->toString();
      dataType = 0;
   }
   cout << "success." << endl;
   cout << "check for composition consistency..." << endl;
   materials.resize(nMat);
   cout << nMat << " materials found." << endl;
   for (iMat = 0; iMat < nMat; iMat++)
      materials(iMat).resize(rSize);
   int ok = 0; // nr of mesh points that have a total composition of 1
   FILE *file = 0;
   SxMatrix3<Double> cell = rBasisPtr->cell;
   if (dataType == 0)  {
      cout << "read binary data" << endl;
      SxBinIO io (inputBinary, SxBinIO::BINARY_READ_ONLY);
      materials = io.readMesh (&cell, &mesh);
      io.close();
   } else if (dataType == 1)  {
      cout << "read ASCII data" << endl;
      ifstream ifile(inputASCII.ascii());
      if (!ifile)   {
         cout << endl << "file "
              << inputASCII.ascii() << " does not exist. EXITING." << endl;
         SX_EXIT;
      }

      try { file = fopen(inputASCII.ascii(), "r+");} catch (const SxException&) { };
      SxVector3<Int> pos;
      for (x = 0; x < xMax; x++)  {
         pos(0) = x;
         for (y = 0; y < yMax; y++)  {
            pos(1) = y;
            for (z = 0; z < zMax; z++)  {
               pos(2) = z;
               int meshIdx = (int)mesh.getMeshIdx(pos, SxMesh3D::Positive);
               double checksum = 0.;
               for (iMat = 0; iMat < nMat; iMat++)  {
                  fscanf(file,"%f\t", &value);
                  materials(iMat)(meshIdx) = value;
                  checksum += value;
               }
               if (fabs(checksum - 1.) < 1.e-5) ok++;
               else if (fabs(checksum - 1.) > 1.e-5)
                  cout << "checksum error at (" << x << ", " << y << ", "
                  << z << "), diff = " << fabs(1. - checksum) << endl;
            }
         }
      }
   }
   // --- check if input map fits mesh size
   cout << "rSize: " << rSize << ", materials(0)'s size: " << materials(0).getSize() << endl;
   if (rSize != materials(0).getSize())  {
      cout << "Mesh size does not fit input map size! EXITING." << endl;
      SX_EXIT;
   }
   cout << "Read strain group" << endl;
   SxSymbolTable *strainGroup = table->getGroup("StrainField");
   /* TODO: COMMENT OUT LATER - THIS FEATURE IS SCARCELY NEEDED AND
      NOT MADE FOR THE RELEASE VERSION OF SPHINX. */
   if (strainGroup->containsGroup("extStrains"))  {
      SxSymbolTable *extStrain = NULL;
      cout << "Read external strains to compute Vpol." << endl;
      readExternalStrains = true;
      extStrain = strainGroup->getGroup("extStrains");
      SxArray<SxString> strainInput; // 3 diagonal, 3 off-diagonal strains
      SxArray<SxMeshR> strainR;
      dataType = -1; // 0 - binary, 1 - ascii
      strainInput.resize(6);
      eIJ.resize(6);
      strainInput(0) = extStrain->get("eXX")->toString();
      strainInput(1) = extStrain->get("eYY")->toString();
      strainInput(2) = extStrain->get("eZZ")->toString();
      strainInput(3) = extStrain->get("eXY")->toString();
      strainInput(4) = extStrain->get("eXZ")->toString();
      strainInput(5) = extStrain->get("eYZ")->toString();
      for (idx = 0; idx < 6; idx++)  {
         if (strainInput(idx).contains(".sxb")) dataType = 0;
         else if (strainInput(idx).contains(".dat")) dataType = 1;
         else {
            cout << strainInput(idx)
               << " has unknown data format. Please provide .dat or.sxb file."
               << " Exiting." << endl;
            SX_EXIT;
         };
         if (dataType == 0)  {
         // --- read binary strain file
            SxBinIO io (strainInput(idx), SxBinIO::BINARY_READ_ONLY);
            strainR.resize (1);
            strainR = io.readMesh (&cell, &mesh);
            eIJ(idx)= strainR(0);
            cout << "open " << strainInput(idx) << endl;
            io.close();

         } else if (dataType == 1)  {
         // --- read ascii strain file
            ifstream ifile(strainInput(idx).ascii());
            if (!ifile)   {
               cout << endl << "file "
                    << strainInput(idx).ascii() << " does not exist. EXITING." << endl;
               SX_EXIT;
            }

            file = fopen(strainInput(idx).ascii(), "r+");
            SxVector3<Int> pos;
            eIJ(idx).resize(rSize);
            for (x = 0; x < xMax; x++)  {
               pos(0) = x;
               for (y = 0; y < yMax; y++)  {
                  pos(1) = y;
                  for (z = 0; z < zMax; z++)  {
                     pos(2) = z;
                     int meshIdx = (int)mesh.getMeshIdx(pos, SxMesh3D::Positive);
                     fscanf(file,"%f\n", &value);
                     eIJ(idx)(meshIdx) = value;
                  }
               }
            } // end x,y,z loops
         }

         if (rSize != eIJ(idx).getSize())  {
            cout << "Strain field size does not match mesh size! EXITING."
                 << endl;
            SX_EXIT;
         }

      } // end loop over strains
//   cout << "test strain: " << endl;
//   for (int idx = 0; idx < eIJ(0).getSize(); idx++)
//      cout << eIJ(0)(idx) << endl;
//   SX_EXIT;
   }
   /* TODO: COMMENT OUT UNTIL HERE. */

   if (generalized)  {
      cout << "Elastic tensor: " << endl;
      cout << "+------------------------------------------+" << endl << " ";
      for (int i = 0; i < 6; i++)  {
         for (int j = 0; j < 6; j++)  {
            if (nonzeroC(i)(j))  {
               if (i < j)  { cout << "C" << i+1 << j+1 << "\t";}
               else { cout << "C" << j+1 << i+1 << "\t";}
            } else {
               cout << " 0 " << "\t";
            }
         }
         cout << endl << " ";
      }
      cout << "+------------------------------------------+" << endl;



      cout << "open Polarisation vector file ";
      SxString polFile;
         polFile = strainGroup->get("polFile")->toString();
      cout << polFile << "... ";
      ifstream ifile(polFile.ascii());
      if (!ifile)   {
         cout << endl << "file "
              << polFile.ascii() << " does not exist. EXITING." << endl;
         SX_EXIT;
      }
      try { file = fopen(polFile.ascii(), "r+");} catch (const SxException&) { };
      char ch;
      SxString polString;
      inFile = "";

      while (!feof(file))  {
         ch = (char)fgetc(file);
         if ((ch != '\n') && (ch != '\t') && (ch != ' '))
            inFile = inFile + ch;
      }
      inFile = inFile.subString(0, inFile.getSize() - 2);
      cout << "success." << endl;
   }

   SxVector3<Double> reference;
   secondOrderPiezo = false;
   AParams = false;
   reference = SxVector3<Double>(strainGroup->get("reference")->toList());
   moreIO = true;
   if (strainGroup->contains("lessIO"))
      moreIO = false;
   crystalStruct = 1;
   rotated = false;
   if (strainGroup->contains("rotated"))  {
      rotated = true;
      crystalStruct = 3;
      }
   if (strainGroup->contains("externalCharge"))  {
      dataType = -1;
      SxString extFile = strainGroup->get("externalCharge")->toString();
      if (extFile.contains(".sxb")) dataType = 0;
      else if (extFile.contains(".dat")) dataType = 1;
      else {
         cout << extFile
            << " has unknown data format. Please provide .dat or.sxb file."
            << " Exiting." << endl;
         SX_EXIT;
      };
      if (dataType == 0)  {
         SxBinIO io (extFile, SxBinIO::BINARY_READ_ONLY);
         extChg.resize (1);
         extChg = (io.readMesh (&cell, &mesh))(0);
         io.close();
      } else if (dataType == 1)  {
         // --- read ascii polarisation file
         ifstream ifile(extFile.ascii());
         if (!ifile)   {
            cout << endl << "file "
                 << extFile.ascii() << " does not exist. EXITING." << endl;
            SX_EXIT;
         }

         file = fopen(extFile.ascii(), "r+");
         SxVector3<Int> pos;
         extChg.resize(rSize);
         for (x = 0; x < xMax; x++)  {
            pos(0) = x;
            for (y = 0; y < yMax; y++)  {
               pos(1) = y;
               for (z = 0; z < zMax; z++)  {
                  pos(2) = z;
                  int meshIdx = (int)mesh.getMeshIdx(pos, SxMesh3D::Positive);
                  fscanf(file,"%f\n", &value);
                  extChg(meshIdx) = value;
               }
            }
         }
      }
      if (rSize != extChg.getSize())  {
         cout << "External charge file size does not match mesh size! EXITING."
              << endl;
         SX_EXIT;
      }
      cout << "total external charge: " << extChg.sum() << endl;
   }
   wgt(0) = 1.; wgt(1) = 1.; wgt(2) = 1.;
   if (strainGroup->contains("weight"))  {
      wgt = SxVector3<Double>(strainGroup->get("weight")->toList());
      for (int i = 0; i < 3; i ++)
         wgt(i) = 1./wgt(i);
      }
   cout << "weighting derivatives with: " << wgt << endl;
   SxArray<SxMeshR> ref;
   ref.resize(3);
   oldC.resize(6);
   sigma0X.resize(rSize);
   sigma0Y.resize(rSize);
   sigma0Z.resize(rSize);
   for (int i = 0; i < 3; i++)  {
      ref(i).resize(rSize);
      for (idx = 0; idx < rSize; idx++)
         ref(i)(idx) = reference(i);
   }
   P.resize(3);
   P2.resize(3);
   P2a.resize(3);
   for (int i = 0; i < 3; i++)  {
      P(i).resize(rSize);
      P2(i).resize(rSize);
      P2a(i).resize(rSize);
   }

   c15given = false;
   for (iParam = 0; iParam < nParam; iParam++)  {
      // --- determine lattice mismatch
      if (!generalizedLat)  {
         if (paramNames(iParam) == "aLat") {
           sigma0X = (ref(0) - parameters(iParam))/parameters(iParam);
           sigma0Y = (ref(1) - parameters(iParam))/parameters(iParam);
           sigma0Z = (ref(2) - parameters(iParam))/parameters(iParam);
         }
         if (paramNames(iParam) == "cLat")  {
           sigma0Z = (ref(2) - parameters(iParam))/parameters(iParam);
         }
      }

      // --- lattice constants
      if (generalizedLat)  {
         latConst.resize(3);
         for (int i = 0; i < 3; i++)  {
            latConst(i).resize(materials(0).getSize());
            latConst(i) = latValues(0)(i) * materials(0);
            for (iMat = 1; iMat < nMat; iMat++)  {
               latConst(i) = latConst(i) + latValues(iMat)(i) * materials(iMat);
            }
         }
         sigma0X = (ref(0) - latConst(0)) / (latConst(0));
         sigma0Y = (ref(1) - latConst(1)) / (latConst(1));
         sigma0Z = (ref(2) - latConst(2)) / (latConst(2));
      }
      // --- elastic constants
      if (generalized)  {
         C.resize(6);
         for (int i = 0; i < 6; i++)  {
            C(i).resize(6);
            for (int j = 0; j < 6; j++)  {
               if (nonzeroC(i)(j))  {
                  C(i)(j).resize(materials(0).getSize());
                  C(i)(j) = cValues(0)(i,j) * materials(0);
                  for (iMat = 1; iMat < nMat; iMat++)  {
                     C(i)(j) = C(i)(j) + cValues(iMat)(i,j) * materials(iMat);
                  }
               } else {
                  C(i)(j).resize(1); C(i)(j)(0) = 0.0;
               }
            }
         }

      }
      if (paramNames(iParam) == "c11") oldC(0) = parameters(iParam);
      if (paramNames(iParam) == "c12") oldC(1) = parameters(iParam);
      if (paramNames(iParam) == "c13") {
         crystalStruct = 2;
         oldC(3) = parameters(iParam);
      }
      if (paramNames(iParam) == "c15")  {
         c15given = true;
         oldC(5) = parameters(iParam);
      }
      if (paramNames(iParam) == "c33") oldC(4) = parameters(iParam);
      if (paramNames(iParam) == "c44") oldC(2) = parameters(iParam);
      if (paramNames(iParam) == "epsilon") epsilon = parameters(iParam);
      if (paramNames(iParam) == "e14") e14 = parameters(iParam);
      if (paramNames(iParam) == "e15") e15 = parameters(iParam);
      if (paramNames(iParam) == "e31") e31 = parameters(iParam);
      if (paramNames(iParam) == "e33") e33 = parameters(iParam);
      if (paramNames(iParam) == "B114") {
         B114 = parameters(iParam);
         secondOrderPiezo = true;
      }

      // --- calculate elastic constants in 111-rotated cell
      if (paramNames(iParam) == "B124") B124 = parameters(iParam);
      if (paramNames(iParam) == "B156") B156 = parameters(iParam);
      if (paramNames(iParam) == "Psp") pSpont = parameters(iParam);
      if (paramNames(iParam) == "A1") {
         A1 = parameters(iParam);
         secondOrderPiezo = true;
         AParams = true;
      }	
      if (paramNames(iParam) == "A2") A2 = parameters(iParam);
   }
   if (rotated)  {
      cout << "calculate rotated parameters" << endl;
      c11S = (1./2.) * (oldC(0) + oldC(1)) + oldC(2);
      c33S = (1./3.) * (oldC(0) + 2. * oldC(1)) + (4./3.) * oldC(2);
      c44S = (1./3.) * (oldC(0) - oldC(1)) + (1./3.) * oldC(2);
      c66S = (1./6.) * (oldC(0) - oldC(1)) + (2./3.) * oldC(2);
      c12S = (1./6.) * (oldC(0) + 5. * oldC(1)) - (1./3.) * oldC(2);
      c13S = (1./3.) * (oldC(0) + 2. * oldC(1)) - (2./3.) * oldC(2);
      if (c15given) {c15S = oldC(5);}
      else
      c15S = (sqrt(2.)/6.) * (oldC(1) - oldC(0)) + (sqrt(2.)/3.) * oldC(2);
      cout << "done. " << endl;
      cout << "values of rotated parameters:" << endl;
      cout << "C11': " << c11S.minval() << " to " << c11S.maxval() << endl;
      cout << "C33': " << c33S.minval() << " to " << c33S.maxval() << endl;
      cout << "C44': " << c44S.minval() << " to " << c44S.maxval() << endl;
      cout << "C66': " << c66S.minval() << " to " << c66S.maxval() << endl;
      cout << "C12': " << c12S.minval() << " to " << c12S.maxval() << endl;
      cout << "C13': " << c13S.minval() << " to " << c13S.maxval() << endl;
      cout << "C15': " << c15S.minval() << " to " << c15S.maxval() << endl;
      cout << "-----------------------------" << endl;
      crystalStruct = 3;
   }
   if (generalized) { crystalStruct = 0;}

   // required for preconditioner
   cP.resize(5);
   for (int i = 0; i < 5; i++)
      cP(i) = oldC(i).sum() / rSize;

   sigma03.resize(3);

   sigma03(0) = sigma0X;
   sigma03(1) = sigma0Y;
   sigma03(2) = sigma0Z;

   dV = rBasisPtr->dOmega;
   cout << "dV: " << dV << endl;

   if (crystalStruct == 0)
      cout << "Using generalised elastic tensor and polarisation" << endl;
   else if (crystalStruct == 1)
      cout << "Calculating in Zincblend material" << endl;
   else if (crystalStruct == 2)
      cout << "Calculating in Wurtzite material" << endl;
   else if (crystalStruct == 3)
      cout << "Calculating in (111) Zincblend material" << endl;

   // --- used for debugging purposes, to be removed later
/*   if (moreIO)  {   // uncomment if required for debugging purposes
      SxMeshR help = sigma0X;
      SxMeshR help2 = sigma0Z;
      SxBinIO ioXX ("helpNew1.sxb", SxBinIO::BINARY_WRITE_ONLY);
      ioXX.writeMesh (help, cell, mesh);
      ioXX.setMode (SxBinIO::WRITE_DATA);
      ioXX.writeMesh (help, cell, mesh);
      ioXX.close();

      SxBinIO ioXX2 ("helpNew2.sxb", SxBinIO::BINARY_WRITE_ONLY);
      ioXX2.writeMesh (help2, cell, mesh);
      ioXX2.setMode (SxBinIO::WRITE_DATA);
      ioXX2.writeMesh (help2, cell, mesh);
      ioXX2.close();
   }*/
//   cout << "successfully finished read. DEBUGGING. EXITING HERE." << endl; SX_EXIT;
}

SxArray<PsiR> SxStrain::gIntegrate (const PsiG &psiG)
{
   const SxRBasis &rBasis = *rBasisPtr;
   int i;
   SxArray<PsiR> gradR(3);

   for (i=0; i < 3; i++)  {
      gradR(i) = (rBasis | (gG2(i)) * I * psiG ) * wgt(i);
   }
   return gradR;
}


SxArray<PsiR> SxStrain::firstDerRComp (const PsiR &psiR) const
{
   SX_CHECK (psiR.getBasisPtr () == &R1);

   SxArray<SxDiracVec<Double> > fdPar;
   fdPar.resize(3);
   SxArray<PsiR> out;
   out.resize(3);

   PsiG psiG = G1 | psiR;
   for (int i=0; i < 3; i++)  {
      //fdPar(i) = (R1 | (psiG * I * gVecs(i).getComponent(0))).real ()*wgt(i);
      fdPar(i) = (R1 | ( psiG * gVecs(i).getComponent(0)) ).imag () * (-wgt(i));
      out(i) = fdPar(i).getComponent(0);
      out(i).setBasis(R1);
   }
   return out;
}

// translating Voigt-Notation Cab to stiffness tensor Cijkl
bool SxStrain::nonZeroElement(int i, int j, int k, int l) const
{
   int a,b;
   if (i == j) {a = i;} // C1*, C2*, C3*
   else {
      if (i + j == 3) {a = 3;} // C4*
      if (i + j == 2) {a = 4;} // C5*
      if (i + j == 1) {a = 5;} // C6*
   }
   if (k == l) {b = k;} // C*1, C*2, C*3
   else {
      if (k + l == 3) {b = 3;} // C*4
      if (k + l == 2) {b = 4;} // C*5
      if (k + l == 1) {b = 5;} // C*6
   }
   return nonzeroC(a)(b);
}

SxMeshR SxStrain::elConst(int i, int j, int k, int l) const
{
   int a,b;
   if (i == j) {a = i;} // C1*, C2*, C3*
   else {
      if (i + j == 3) {a = 3;} // C4*
      if (i + j == 2) {a = 4;} // C5*
      if (i + j == 1) {a = 5;} // C6*
   }
   if (k == l) {b = k;} // C*1, C*2, C*3
   else {
      if (k + l == 3) {b = 3;} // C*4
      if (k + l == 2) {b = 4;} // C*5
      if (k + l == 1) {b = 5;} // C*6
   }
   return (C(a)(b));

}

double delta(int i, int j)
{
   if (i == j) return 1.0;
   else return 0.0;
}

SxMeshR SxStrain::strain(int i, int j)
{
   SxMeshR rtn = eXX;
   if ((i == 0) && (j == 0)) {rtn = eXX;}
   if ((i == 1) && (j == 1)) {rtn = eYY;}
   if ((i == 2) && (j == 2)) {rtn = eZZ;}

   if (((i == 0) && (j == 1)) || ((i == 1) && (j == 0))) {rtn = eXY;}
   if (((i == 0) && (j == 2)) || ((i == 2) && (j == 0))) {rtn = eXZ;}
   if (((i == 1) && (j == 2)) || ((i == 2) && (j == 1))) {rtn = eYZ;}
   return rtn;
}

double SxStrain::getFreeEnergy()
{
   double rtn = 0.;
   int i, j, k, l;
   i = 0;
   SxMeshR sum = zero;
   if (generalized)  { // NEW CODE
      for (j = 0; j < 3; j ++)
         for (k = 0; k < 3; k ++)
            for (l = 0; l < 3; l ++)  {
               if (nonZeroElement(i,j,k,l))  {
                  sum = sum + (elConst(i,j,k,l) * strain(i,j) * strain(k,l));
               }
            }
   } else {
      cout << "total free energy not available in deprecate elasticity format." << endl;
      rtn = -1;
   }
   rtn = 0.5 * sum.sum() * dV * 0.00092489629; // converts GPa * bohrradius^3 to eV
   fMax = sum.maxval(); fMin = sum.minval();
   return rtn;
}

SxArray<PsiR> SxStrain::uInR (const PsiG u) const
{
   SxArray<PsiR> uR(3);
   const SxRBasis &rBasis = *rBasisPtr;
   SxDiracVec<Double> uRallComponents = rBasis | u;
   for (int i = 0; i < 3; i++)  {
      // 10^7 is to make calculation numerically possible.
      // order of magnitude of u's and other contributions differ at about
      // 10^7 -> without this factor minimisation doesn't work for numerical
      // reasons
      uR(i) = 1.e7 * uRallComponents.getComponent(i);
      uR(i).setBasis(&R1);
   }
   return uR;
}

PsiR SxStrain::gradient (const PsiG &u) const
{
   SX_CHECK (rBasisPtr);

   int i,j,k,l;
   int i1, i2;

//   SxArray<PsiR> uR;
   uStore = u;
   SxArray<PsiR> uR = uInR(u);
   PsiR grad(*rBasisPtr); grad.set (0.);
   PsiR dF(*rBasisPtr); dF.set (0.);
   SxArray<PsiR> fd(3);
   if (!generalized)
      for (i = 0; i < 3; i++)  {
         fd(i) = firstDerRComp(uR(i))(i);  //dUi/dRi
      }

   if (!generalized)
      const_cast<SxStrain*>(this)->calculateSigma(1.e9*u, uR);

   for (i = 0; i < 3; i++)  {
      i1 = (i + 1) % 3;
      i2 = (i + 2) % 3;
      if (generalized)  { // NEW CODE
         dF = 0.* uR(0);
         for (j = 0; j < 3; j ++)
            for (k = 0; k < 3; k ++)
               for (l = 0; l < 3; l ++)  {
                  if (nonZeroElement(i,j,k,l))  {
                     dF = dF +
                     firstDerRComp (elConst(i,j,k,l)
                     * (firstDerRComp(uR(k))(l) + delta(k,l)*sigma03(l)))(j);
                  }
               }
      } else
      if (crystalStruct == 1)  { // Zincblende

         dF = firstDerRComp(
                  oldC(0) * ( fd(i) + sigma03(i) )
                  + oldC(1) * (fd(i1) + fd(i2) + sigma03(i1) + sigma03(i2))
                  )(i)

         + 1.0 * firstDerRComp(oldC(2)
               * (firstDerRComp(uR(i1))(i)
                  + firstDerRComp(uR(i))(i1)))(i1)
         + 1.0 * firstDerRComp(oldC(2)
               * (firstDerRComp(uR(i2))(i)
                  + firstDerRComp(uR(i))(i2)))(i2)
            ;
      }


      /*
       Cxxxx = Cyyyy                                  = C11 = oldC(0)
       Cxxyy = Cyyxx                                  = C12 = oldC(1)
       Cxyxy = Cyxyx = Cxzxz = Czxzx = Cyzyz = Czyzy  = C44 = oldC(2)
       Cxxzz = Czzxx = Cyyzz = Czzyy                  = C13 = oldC(3)
       Czzzz                                          = C33 = oldC(4)
       C                                              = C15 = oldC(5)
       */
      else if (crystalStruct == 2)  { // Wurtzite

         if (i == 0)
            dF = firstDerRComp(
                  oldC(0) * ( fd(0) + sigma03(0) )
                  + oldC(1) * (fd(1) + sigma03(1))
                  + oldC(3) * (fd(2) + sigma03(2))
                  )(i)
         + 1. * firstDerRComp(oldC(2)
               * (firstDerRComp(uR(i2))(i)
                  + firstDerRComp(uR(i))(i2)))(i2)
            + firstDerRComp(0.5*(oldC(0) - oldC(1))
              * (firstDerRComp(uR(0))(1)
                + firstDerRComp(uR(1))(0)))(1)
         ;
         if (i == 1)
            dF = firstDerRComp(
                  oldC(0) * ( fd(1) + sigma03(1) )
                  + oldC(1) * (fd(0) + sigma03(0))
                  + oldC(3) * (fd(2) + sigma03(2))
                  )(i)
         + 1. * firstDerRComp(oldC(2)
               * (firstDerRComp(uR(i1))(i)
                  + firstDerRComp(uR(i))(i1)))(i1)
         + firstDerRComp(0.5*(oldC(0) - oldC(1))
              * (firstDerRComp(uR(0))(1)
                + firstDerRComp(uR(1))(0)))(0)
               ;

         if (i == 2)
            dF = firstDerRComp(
                  oldC(4) * ( fd(2) + sigma03(2) )
                  + oldC(3) * (fd(0) + sigma03(0))
                  + oldC(3) * (fd(1) + sigma03(1))
                  )(i)
         + 1. * firstDerRComp(oldC(2)
               * (firstDerRComp(uR(i1))(i)
                  + firstDerRComp(uR(i))(i1)))(i1)
         + 1. * firstDerRComp(oldC(2)
               * (firstDerRComp(uR(i2))(i)
                  + firstDerRComp(uR(i))(i2)))(i2)
               ;

      }
      else if (crystalStruct == 3)  { // Zincblende, 111
         if (i == 0)
            dF = firstDerRComp(
                  c11S * ( fd(0) + sigma03(0) )
                  + c12S * (fd(1) + sigma03(1))
                  + c13S * (fd(2) + sigma03(2))
                  )(0)
            + firstDerRComp(c66S
               * (firstDerRComp(uR(1))(0)
                  + firstDerRComp(uR(0))(1)))(1)
            + firstDerRComp(c44S
               * (firstDerRComp(uR(2))(0)
                  + firstDerRComp(uR(0))(2)))(2)
            + firstDerRComp(c15S * (firstDerRComp(uR(2))(0) + firstDerRComp(uR(0))(2)))(0)
            - firstDerRComp(c15S * (firstDerRComp(uR(2))(1) + firstDerRComp(uR(1))(2)))(1)
            + firstDerRComp(c15S * (fd(0) - fd(1) + sigma03(0) - sigma03(1)))(2);

         if (i == 1)
            dF = firstDerRComp(
                  c11S * ( fd(1) + sigma03(1) )
                  + c13S * (fd(2) + sigma03(2))
                  + c12S * (fd(0) + sigma03(0))
                  )(1)
            + firstDerRComp(c44S
               * (firstDerRComp(uR(2))(1)
                  + firstDerRComp(uR(1))(2)))(2)
            + firstDerRComp(c66S
               * (firstDerRComp(uR(0))(1)
                  + firstDerRComp(uR(1))(0)))(0)
            - firstDerRComp(c15S * (firstDerRComp(uR(2))(1) + firstDerRComp(uR(1))(2)))(0)
            - firstDerRComp(c15S * (firstDerRComp(uR(2))(0) + firstDerRComp(uR(0))(2)))(1)
            - firstDerRComp(c15S * (firstDerRComp(uR(0))(1) + firstDerRComp(uR(1))(0)))(2);

         if (i == 2)
            dF = firstDerRComp(
                  c33S * ( fd(2) + sigma03(2) )
                  + c13S * (fd(0) + fd(1) + sigma03(0) + sigma03(1))
                  )(2)
            + firstDerRComp(c44S
               * (firstDerRComp(uR(0))(2)
                  + firstDerRComp(uR(2))(0)))(0)
            + firstDerRComp(c44S
               * (firstDerRComp(uR(1))(2)
                  + firstDerRComp(uR(2))(1)))(1)
	    + firstDerRComp(c15S * (fd(0) + sigma03(0) - fd(1) - sigma03(1)))(0)
            - firstDerRComp(c15S * (firstDerRComp(uR(0))(1) +  firstDerRComp(uR(1))(0)))(1);

      }
      if ((step == 1) && (containsVacuum)) {
         if (i == 0)
         cout << "remove unphysical displacements in vacuum..." << endl;
         dF = -uR(i) * materials(vacMat);
      }
      if (containsVacuum)
         dF -= uR(i) * materials(vacMat);

      grad.setComponent(i,-1. * dF);
   }
   return grad;
}

SxString SxStrain::replaceUnknown(SxString in)
{
   SxString rtn = in;
   // tokenise
   in = in.substitute("+", "|");
   in = in.substitute("-", "|");
   in = in.substitute("*", "|");
   in = in.substitute("/", "|");
   in = in.substitute("(", "|");
   in = in.substitute(")", "|");
   in = in.substitute("{", "|");
   in = in.substitute("}", "|");
   in = in.substitute(",", "|");
   SxList<SxString> split = in.tokenize('|');

   int idx;
   for (idx = 0; idx <  split.getSize(); idx++)  {
      if (  // keyword, double, function or known material parameter
          (split(idx) == "eXX") ||
          (split(idx) == "eYY") ||
          (split(idx) == "eZZ") ||
          (split(idx) == "eXY") ||
          (split(idx) == "eXZ") ||
          (split(idx) == "eYZ") ||
          (split(idx) == "Tr") ||
          (split(idx) == "Pi") ||
          (split(idx) == "Sin") ||
          (split(idx) == "Cos") ||
          (split(idx) == "Tan") ||
          (split(idx) == "Exp") ||
          (split(idx) == "Root") ||
          (split(idx) == "Sqrt") ||
          (split(idx).isDouble()) ||
          (paramNames.contains(split(idx)))
         ) {/* empty */}
      else { // check wether unknown elements are defined above
         if (inFile.contains(split(idx) + "=") < 1)  { // if not: exit
            cout << "Element " << split(idx) << " was not found anywhere.\n";
            SX_EXIT;
         } else  { // if yes, replace element by its expression
            SxString replace = "(" + inFile.right(split(idx) + "=").left(";") + ")";
            if (replace.contains("=") > 0)  { // if expression lacks semicolon, exit
               cout << "Error: semicolon missing after " + split(idx) + "=\n";
               SX_EXIT;
            }
            rtn = rtn.substitute(split(idx), replaceUnknown(replace));
         }
      }
   }
   return rtn;
}

int SxStrain::firstOutsideBracket(SxString expr, char op)  {
   int pos = -1;
   int bLevel = 0;
   int i;

   for (i = 0; i < expr.getSize(); i++)  {
      if (expr(i) == '(') bLevel++;
      if (expr(i) == '{') bLevel++;
      if (expr(i) == ')') bLevel--;
      if (expr(i) == '}') bLevel--;
      if ((expr(i) == op) && (bLevel == 0) && (pos == -1))
         pos = i;
   }

   return pos;
}

double SxStrain::parseFunction (SxString in)  {
   double rtn = 1.;
   double tmp;

   SxString left, right, newIn;
   // check for '+'
   int pos = -1;
   int i;
   pos = firstOutsideBracket(in, '+');
   if (pos > -1)  {
      left = in.subString(0, pos -1);
      right = in.subString(pos +1);
      if (right.getSize() == 0) {
         cout << "Syntax error: '+' without right side argument." << endl;
         SX_EXIT;
      }
      if (left.getSize() == 0) {
         rtn = parseFunction(right);
      } else {
         rtn = parseFunction (left)
             + parseFunction (right);
      }
   }
   // check for '-'
   if (pos == -1)  { // no '+' outside brackets found
      pos = firstOutsideBracket(in, '-');
      if (pos > -1)  {
         left = in.subString(0, pos -1);
         right = in.subString(pos +1);
         if (right.getSize() == 0) {
            cout << "Syntax error: '-' without right side argument." << endl;
            SX_EXIT;
         }

         if (left.getSize() == 0) {
            rtn = - parseFunction(right);
         } else {
            rtn = parseFunction (left)
                - parseFunction (right);
         }
      }
   }
   // check for '*'
   if (pos == -1)  { // no '+', '-' outside brackets found
      pos = firstOutsideBracket(in, '*');
      if (pos > -1)  {
         left = in.subString(0, pos -1);
         right = in.subString(pos +1);
         if ((left.getSize() < 1) || (right.getSize() < 1))  {
            cout << "Syntax error: '*' with only one argument." << endl;
            SX_EXIT;
         }
         rtn = parseFunction (left)
             * parseFunction (right);
      }
   }
   // check for '/'
   if (pos == -1)  { // no '+', '-', '*' outside brackets found
      pos = firstOutsideBracket(in, '/');
      if (pos > -1)  {
         left = in.subString(0, pos -1);
         right = in.subString(pos +1);
         if ((left.getSize() < 1) || (right.getSize() < 1))  {
            cout << "Syntax error: '/' with only one argument." << endl;
            SX_EXIT;
         }
         tmp = parseFunction (right);
         if (tmp != 0)  {rtn = parseFunction (left) / tmp;}
         else {cout << "division by zero detected." << endl; SX_EXIT;}
      }
   }
   // check for '()'
   if (pos == -1)  { // no mathematical operator found
      if (in(0) == '(')  {
         newIn = "";
         for (i = 1; i < in.getSize()-1; i++)   {// remove 1st and last char --> '(' and ')'
            newIn = newIn + in(i);
         }
         rtn = parseFunction (newIn);
         pos = 1;
      }
   }

   if (pos == -1)  {
      if (in.isDouble()) {rtn = in.toDouble(); pos = 1;}
   }
   if (pos == -1)  { // no mathematical operator or bracket found, elem is no double either
      SxString inShort = in;
      inShort.resize(5, true);
      // --- Sqrt function
      if (inShort.contains("Sqrt{") > 0)  {
         newIn = in.right("Sqrt{");
         newIn.resize(newIn.getSize() - 1, true);
         tmp = parseFunction(newIn);
         if (tmp >= 0) {rtn = sqrt(tmp);}
         else {cout << "Sqrt of negative number detected." << endl; SX_EXIT;}
      }
      // --- Pi value
      if (inShort.contains("Pi") > 0) {
            rtn = PI;
         }
      // --- Sin function
      if (inShort.contains("Sin{") > 0)  {
         newIn = in.right("Sin{");
         newIn.resize(newIn.getSize() - 1, true);
         rtn = sin(parseFunction(newIn));
      }
      // --- Cos function
      if (inShort.contains("Cos{") > 0)  {
         newIn = in.right("Cos{");
         newIn.resize(newIn.getSize() - 1, true);
         rtn = cos(parseFunction(newIn));
      }
      // --- Tan function
      if (inShort.contains("Tan{") > 0)  {
         newIn = in.right("Tan{");
         newIn.resize(newIn.getSize() - 1, true);
         rtn = tan(parseFunction(newIn));
      }
      // --- Exp function
      if (inShort.contains("Exp{") > 0)  {
         newIn = in.right("Exp{");
         newIn.resize(newIn.getSize() - 1, true);
         rtn = exp(parseFunction(newIn));
      }
      // --- Root function
      if (inShort.contains("Root{") > 0)  {
         double value,e;
         newIn = in.right("Root{");
         newIn.resize(newIn.getSize() - 1, true);
         if (newIn.contains(",") < 1)  {
            cout << "root exponent missing. Usage: Root{a,b}" << endl;
            SX_EXIT;
         } else {
            value = parseFunction(newIn.left(","));
            e = 1./parseFunction(newIn.right(","));
         }
         if ((value >= 0) && (e != 0)) {
            rtn = pow(value,e);
         } else {
            cout << "Negative number or zero exponent detected." << endl;
            SX_EXIT;
         }
      }

   }

   return rtn;
}

PsiR SxStrain::parseElement(SxString in)  {
   PsiR rtn;
   rtn = zero;
   SxString left, right, newIn;
   // check for '+'
   int pos = -1;
   int i, iParam;
   pos = firstOutsideBracket(in, '+');
   if (pos > -1)  {
      left = in.subString(0, pos -1);
      right = in.subString(pos +1);
      if (right.getSize() == 0) {
         cout << "Syntax error: '+' without right side argument." << endl;
         SX_EXIT;
      }
      if (left.getSize() == 0) {
         rtn = parseElement(right);
      } else {
         rtn = parseElement (left)
             + parseElement (right);
      }
   }
   // check for '-'
   if (pos == -1)  { // no '+' outside brackets found
      pos = firstOutsideBracket(in, '-');
      if (pos > -1)  {
         left = in.subString(0, pos -1);
         right = in.subString(pos +1);
         if (right.getSize() == 0) {
            cout << "Syntax error: '-' without right side argument." << endl;
            SX_EXIT;
         }

         if (left.getSize() == 0) {
            rtn = - parseElement(right);
         } else {
            rtn = parseElement (left)
                - parseElement (right);
         }
      }
   }
   // check for '*'
   if (pos == -1)  { // no '+', '-' outside brackets found
      pos = firstOutsideBracket(in, '*');
      if (pos > -1)  {
         left = in.subString(0, pos -1);
         right = in.subString(pos +1);
         if ((left.getSize() < 1) || (right.getSize() < 1))  {
            cout << "Syntax error: '*' with only one argument." << endl;
            SX_EXIT;
         }
         rtn = parseElement (left)
             * parseElement (right);
      }
   }
   // check for '/'
   if (pos == -1)  { // no '+', '-', '*' outside brackets found
      pos = firstOutsideBracket(in, '/');
      if (pos > -1)  {
         left = in.subString(0, pos -1);
         right = in.subString(pos +1);
         if ((left.getSize() < 1) || (right.getSize() < 1))  {
            cout << "Syntax error: '/' with only one argument." << endl;
            SX_EXIT;
         }

         rtn = parseElement (left)
                / parseElement (right);
      }
   }
   // check for '()'
   if (pos == -1)  { // no mathematical operator found
      if (in(0) == '(')  {
         newIn = "";
         for (i = 1; i < in.getSize()-1; i++)   {// remove 1st and last char --> '(' and ')'
            newIn = newIn + in(i);
         }
         rtn = parseElement (newIn);
         pos = 1;
      }
   }
   if (pos == -1)  { // no mathematical operator & no brackets found
      if (in.isDouble()) {rtn = in.toDouble() * one; pos = 1;}
      if (in == "eXX") {rtn = eXX; pos = 1;}
      if (in == "eYY") {rtn = eYY; pos = 1;}
      if (in == "eZZ") {rtn = eZZ; pos = 1;}
      if (in == "eXY") {rtn = eXY; pos = 1;}
      if (in == "eXZ") {rtn = eXZ; pos = 1;}
      if (in == "eYZ") {rtn = eYZ; pos = 1;}
      if (in == "Tr") {rtn = eXX + eYY + eZZ; pos = 1;}
      for (iParam = 0; iParam < paramNames.getSize(); iParam++)  {
         if (in == paramNames(iParam)) {rtn = parameters(iParam); pos = 1;}
      }
   }
   if (pos == -1)  {
      rtn = parseFunction(in) * one;
   }
   return rtn;
}

SxArray<PsiR> SxStrain::parsePolarisation()  {
//   cout << "Start parser..." << endl;
   if (inFile.contains("Px=") < 1) {cout << "Error: Px missing.\n"; SX_EXIT;}
   if (inFile.contains("Py=") < 1) {cout << "Error: Py missing.\n"; SX_EXIT;}
   if (inFile.contains("Pz=") < 1) {cout << "Error: Pz missing.\n"; SX_EXIT;}

   SxString px = inFile.right("Px=").left(";");
   if (px.contains("=") > 0)  {
      cout << "Error: semicolon missing after Px. EXITING." << endl;
      SX_EXIT;
   }
   SxString py = inFile.right("Py=").left(";");
   if (py.contains("=") > 0)  {
      cout << "Error: semicolon missing after Py. EXITING." << endl;
      SX_EXIT;
   }
   if (inFile.right("Pz=").contains(";") < 1)  {
      cout << "Error: semicolon missing after Pz. EXITING." << endl;
      SX_EXIT;
   }
   SxString pz = inFile.right("Pz=").left(";");
   px = replaceUnknown(px);
   py = replaceUnknown(py);
   pz = replaceUnknown(pz);
   cout << "Px = " + px << endl;
   cout << "Py = " + py << endl;
   cout << "Pz = " + pz << endl;
   SxArray<PsiR> rtn;
   rtn.resize(3);
   rtn(0) = parseElement(px) / epsilon;
   rtn(1) = parseElement(py) / epsilon;
   rtn(2) = parseElement(pz) / epsilon;
   return rtn;
}

void SxStrain::calculateSigma (const SxArray<PsiR> &uR)
{
   // --- calculating piezoelectric field
   int i;
   const SxRBasis &rBasis = *rBasisPtr;
   const SxGBasis &gBasis = *dynamic_cast<const SxGBasis *>(uStore.getBasisPtr());

   PsiR peComps = (rBasis | uStore);
   cout << "employ generalised polarisation..." << endl;
   if (readExternalStrains) cout << "read strains from external files instead..." << endl;
   eXX = firstDerRComp(uR(0))(0) + sigma0X;
   eYY = firstDerRComp(uR(1))(1) + sigma0Y;
   eZZ = firstDerRComp(uR(2))(2) + sigma0Z;

   eXY = .5 * (firstDerRComp(uR(1))(0)
         + firstDerRComp(uR(0))(1));
   eXZ = .5 * (firstDerRComp(uR(2))(0)
         + firstDerRComp(uR(0))(2));
   eYZ = .5 * (firstDerRComp(uR(2))(1)
         + firstDerRComp(uR(1))(2));
   fTot = getFreeEnergy();
   double au = 0.017478; // C/m^2 -> e / bohrradius^2
   P = parsePolarisation();
   cout << "max. of polarisation: " << P(0)(0) << ", " << P(1)(0) << ", " << P(2)(0) << endl;

   for (i = 0; i < 3; i++)  {
      peComps.setComponent(i, au*P(i));
   }

   double NM = 18.897261;
   PsiR chgR;
//   cout << "extChg's size: " << extChg.getSize() << endl;
   PsiG vG;
   if (extChg.getSize() > 0)  {
      chgR = 0.0 * peComps.getComponent(0) + extChg;
//      cout << "chgR's size: " << chgR.getSize() << endl;
      vG = (gBasis | (chgR / epsilon)).getComponent(0)
           * wgt(0) * wgt(1) * wgt(2) / (NM*scaledG2);
//      vG = (gBasis | (chgR / epsilon)).getComponent(0) / gBasis.g2;
      vG(0) = 0;
//      cout << "vG's size: " << vG.getSize() << endl;
//      cout << "gSize: " << uStore.getSize() / 3 << endl;
   }

   if (generalized)  {
      vP = -FOUR_PI * HA2EV * (
        gIntegrate(gBasis|peComps)(0).getComponent(0)
      + gIntegrate(gBasis|peComps)(1).getComponent(1)
      + gIntegrate(gBasis|peComps)(2).getComponent(2));
      if (extChg.getSize() > 0) {
         PsiR chgPot = (rBasis | vG).getComponent(0);
         vP = vP + chgPot;
      }
   } else { vP = -(FOUR_PI*HA2EV)*(gIntegrate(gBasis|peComps)(0).getComponent(0)
      + gIntegrate(gBasis|peComps)(1).getComponent(1)
      + gIntegrate(gBasis|peComps)(2).getComponent(2)
      );
   }
   Tr = eXX + eYY + eZZ;
}

void SxStrain::calculateSigma (const PsiG &u, const SxArray<PsiR> &uR)
{
   SxMeshR sigmaXXP, sigmaYYP, sigmaZZP, C1, C2, C3, C4, C5;

   double au = 0.017478; // C/m^2 -> e / bohrradius^2
   eXX = firstDerRComp(uR(0))(0) + sigma0X;
   eYY = firstDerRComp(uR(1))(1) + sigma0Y;
   eZZ = firstDerRComp(uR(2))(2) + sigma0Z;

   eXY = .5 * (firstDerRComp(uR(1))(0)
         + firstDerRComp(uR(0))(1));
   eXZ = .5 * (firstDerRComp(uR(2))(0)
         + firstDerRComp(uR(0))(2));
   eYZ = .5 * (firstDerRComp(uR(2))(1)
         + firstDerRComp(uR(1))(2));

   /* TODO: COMMENT OUT AFTER USAGE - NOT FOR RELEASE VERSION  */
   if (readExternalStrains)  {
      cout << "use external strains for polarization..." << endl;
      eXX = eIJ(0);
      eYY = eIJ(1);
      eZZ = eIJ(2);
      eXY = eIJ(3);
      eXZ = eIJ(4);
      eYZ = eIJ(5);
   }
   /* TODO: END COMMENT */
   sigmaXXP = eXX;
   sigmaYYP = eYY;
   sigmaZZP = eZZ;


   sigmaXY = eXY;
   sigmaXZ = eXZ;
   sigmaYZ = eYZ;

/*   sigmaXXRot = (1./6.)*sigmaXXP + (1./2.)*sigmaYYP + (1./3.)*sigmaZZP
    - (1./sqrt(3.))*sigmaXY + (sqrt(2.)/3.)*sigmaXZ - sqrt(2./3.)*sigmaYZ;
   sigmaYYRot = (1./6.)*sigmaXXP + (1./2.)*sigmaYYP + (1./3.)*sigmaZZP
    + (1./sqrt(3.))*sigmaXY + (sqrt(2.)/3.)*sigmaXZ + sqrt(2./3.)*sigmaYZ;
   sigmaZZRot = (2./3.)*sigmaXXP + (1./3.)*sigmaZZP - sqrt(2.)*(2./3.)*sigmaXZ;
   sigmaXYRot = (1./6.)*sigmaXXP - (1./2.)*sigmaYYP + (1./3.)*sigmaZZP
    + sqrt(2./9.)*sigmaXZ;
   sigmaXZRot = (1./3.)*(sigmaZZP - sigmaXXP) - (1./(3.*sqrt(2.)))*sigmaXZ
    + sqrt(1./3.)*sigmaXY - sqrt(1./6.)*sigmaYZ;
   sigmaYZRot = (1./3.)*(sigmaZZP - sigmaXXP) - (1./(3.*sqrt(2.)))*sigmaXZ
    + sqrt(1./3.)*sigmaXY - sqrt(1./6.)*sigmaYZ;*/
   Tr = sigmaXXP + sigmaYYP + sigmaZZP;

   // --- calculating piezoelectric field
   int i;
   const SxRBasis &rBasis = *rBasisPtr;
   const SxGBasis &gBasis = *dynamic_cast<const SxGBasis *>(u.getBasisPtr());

   PsiR peComps = (rBasis | u);
   if (crystalStruct == 0)  {
      cout << "employ generalised polarisation..." << endl;
      eXX = firstDerRComp(uR(0))(0) + sigma0X;
      eYY = firstDerRComp(uR(1))(1) + sigma0Y;
      eZZ = firstDerRComp(uR(2))(2) + sigma0Z;

      eXY = .5 * (firstDerRComp(uR(1))(0)
            + firstDerRComp(uR(0))(1));
      eXZ = .5 * (firstDerRComp(uR(2))(0)
            + firstDerRComp(uR(0))(2));
      eYZ = .5 * (firstDerRComp(uR(2))(1)
            + firstDerRComp(uR(1))(2));
   /* TODO: COMMENT OUT AFTER USAGE - NOT FOR RELEASE VERSION  */
   if (readExternalStrains)  {
      cout << "use external strains for polarization..." << endl;
      eXX = eIJ(0);
      eYY = eIJ(1);
      eZZ = eIJ(2);
      eXY = eIJ(3);
      eXZ = eIJ(4);
      eYZ = eIJ(5);
   }
   /* TODO: END COMMENT */

      P = parsePolarisation();
   } else if (crystalStruct == 1)  { // zincblend materials
      P (0) = 2. * e14 * sigmaYZ / (epsilon);
      P (1) = 2. * e14 * sigmaXZ / (epsilon);
      P (2) = 2. * e14 * sigmaXY / (epsilon);
      if (secondOrderPiezo)  {
      P2(0) = 2. * B114 * sigmaXXP * sigmaYZ / epsilon
                  + 2. * B124 * sigmaYZ * (sigmaYYP + sigmaZZP) / epsilon
                  + 4. * B156 * sigmaXZ * sigmaXY / epsilon;
      P2(1) = 2. * B114 * sigmaYYP * sigmaXZ / epsilon
                  + 2. * B124 * sigmaXZ * (sigmaXXP + sigmaZZP) / epsilon
                  + 4. * B156 * sigmaXY * sigmaYZ / epsilon;
      P2(2) = 2. * B114 * sigmaZZP * sigmaXY /epsilon
                  + 2. * B124 * sigmaXY * (sigmaXXP + sigmaYYP) /epsilon
                  + 4. * B156 * sigmaXZ * sigmaYZ / epsilon;

      P2a(0) = 2.* (A1 * (sigmaXXP + sigmaYYP + sigmaZZP) * sigmaYZ
                    + A2 * sigmaYZ *
                    ( sigmaXXP - 0.5 * (sigmaYYP + sigmaZZP) )
                    ) / epsilon
                  + (4./sqrt(6.)) * B156
                    * (sigmaXZ * sigmaXY + sigmaYZ * sigmaXY
                       - 2.*sigmaYZ*sigmaXZ) / epsilon;
      P2a(1) = 2.* (A1 * (sigmaXXP + sigmaYYP + sigmaZZP) * sigmaXZ
                    + A2 * sigmaXZ *
                    ( sigmaYYP - 0.5 * (sigmaXXP + sigmaZZP) )
                    ) / epsilon
                  + (4./sqrt(2.)) * B156
                    * (sigmaYZ * sigmaXY - sigmaXZ * sigmaXY) / epsilon;
      P2a(2) = 2.* (A1 * (sigmaXXP + sigmaYYP + sigmaZZP) * sigmaXY
                    + A2 * sigmaXY
                      * ( sigmaZZP - 0.5 * (sigmaXXP + sigmaYYP) )
                    ) / epsilon
                  + (4./sqrt(3.)) * B156 *
                    (sigmaXZ * sigmaXY + sigmaYZ * sigmaXY
                     + sigmaYZ * sigmaXZ) / epsilon;
      }

   } else if (crystalStruct == 2)  { // wurtzite materials
      P (0) =  e15 * sigmaXZ * 2.0 / (epsilon);
      P (1) =  e15 * sigmaYZ * 2.0 / (epsilon);
      P (2) =  (e31 * (sigmaXXP + sigmaYYP)
            + e33 * (sigmaZZP) + pSpont) / (epsilon);
   } else if (crystalStruct == 3)  { // zincblende (111)
      Tr = sigmaXXP + sigmaYYP + sigmaZZP;
      C1 = (1./3.) * (sigmaZZP - sigmaXXP) - (1./(3.*sqrt(2.))) * sigmaXZ;
      C2 = (1./sqrt(6.)) * sigmaYZ - (1./sqrt(3.)) * sigmaXY;
      C3 = 2. * sigmaXY - sqrt(2.) * sigmaYZ;
      C4 = 0.25 * (sigmaYYP - sigmaXXP) + sqrt(0.5) * sigmaXZ;
      C5 = sqrt(3./2.) * sigmaYZ + sqrt(3./4.) * sigmaXY;
      P (0) = ((-sqrt(1./3.)) * e14
                     * (-sqrt(2.) * (sigmaYYP - sigmaXXP) + 2. * sigmaXZ)
                     / epsilon);
      P (1) = ((-sqrt(1./3.)) * e14 * 2.
                     * (sigmaYZ - sqrt(2.) * sigmaXY) / epsilon);
      P (2) = ((2./sqrt(3.)) * e14 * (sigmaZZP
                      - 0.5 * (sigmaXXP + sigmaYYP)) / epsilon);
      // --- simplified version to simulate wurtzite
//      P (0) = ((-sqrt(1./3.)) * e14 * (2. * sigmaXZ) / epsilon);
//      P (1) = ((-sqrt(1./3.)) * e14 * 2. * (sigmaYZ) / epsilon);
//      P (2) = ((2./sqrt(3.)) * e14 * (sigmaZZP - 0.5 * (sigmaXXP + sigmaYYP)) / epsilon);
      if (secondOrderPiezo)  {
//      cout << "employ 2nd-order polarisation..." << endl;
      P2 (0) = 2. * (-sqrt(1./3.) * (B114 + B124)
                  * (2. * sigmaXZ - sqrt(2.) * (sigmaYYP - sigmaXXP))
                  * (sigmaXXP + sigmaYYP + sigmaZZP)) / epsilon
                  ;
      P2 (1) = 2. * (-sqrt(1./3.) * (B114 + B124)
                  * (2. * sigmaYZ - 2. * sqrt(2.) * sigmaXY)
                  * (sigmaXXP + sigmaYYP + sigmaZZP)) / epsilon
                  ;
      P2 (2) = 2. * (-sqrt(1./3.) * (B114 + B124)
                  * (sigmaXXP + sigmaYYP - 2. * sigmaZZP)
                  * (sigmaXXP + sigmaYYP + sigmaZZP)) / epsilon
                  ;
      P2a (0) = (
                  (-sqrt(1./3.)) * A1 * Tr
                  * (-sqrt(2.) * (sigmaYYP - sigmaXXP) + 2. * sigmaXZ)
                  // --- A2 term
                  + 2. * A2 * (
                  (sqrt(1./6.) * (sigmaYYP - sigmaXXP) + sqrt(4./3.) * sigmaXZ)
                  * (0.5 * (sigmaZZP - sigmaYYP) + sqrt(1./8.) * sigmaXZ)
                  + sqrt(1./6.) * (sqrt(1./2.) * sigmaYZ + 0.5 * sigmaXY) * C3)
                  // --- B156 term
                  + 2. * B156 * sqrt(2./3.) * (C1
                  * (0.5 * (sigmaXXP - sigmaYYP) + sqrt(0.5) * sigmaXZ)
                  + C2 * C2)
                  ) / epsilon
                  ;
      P2a (1) = (
                  (-sqrt(2./3.)) * A1 * Tr * 2.
                   * (sigmaYZ - sqrt(2.) * sigmaXY)
                  // --- A2 term
                  + 2. * A2 * (sqrt(1./6.) * C4 * C3
                  + sqrt(1./18.) * C5
                  * (2. * (sigmaZZP - sigmaXXP) - sqrt(2.) * sigmaXZ) )
                  // --- B156 term
                  + 2. * B156 * sqrt(2.) * C2
                  * ((1./6.) * sigmaXXP - 0.5 * sigmaYYP
                    + (1./3.) * sigmaZZP + (sqrt(2.)/3.) * sigmaXZ)
                  ) / epsilon
                  ;
      P2a (2) = (
                  (2./sqrt(3.)) * A1 * Tr
                  * (sigmaZZP - 0.5 * (sigmaXXP + sigmaYYP))
                  // --- A2 term
                  + 2. * A2 * (sqrt(1./3.) * C4
                  * (sigmaYYP - sigmaXXP - sqrt(2.) * sigmaXZ)
                  + C5 * C3 )
                  // --- B156 term
                  + 2. * B156 * (1./sqrt(3.)) * (
                  C1 * (sigmaZZP - sigmaYYP + sqrt(0.5) * sigmaXZ) - C2 * C2)
                  ) / epsilon
                  ;
      }
   }
   // no second-order piezoelectricity for wurtzite included
   if (crystalStruct == 2) secondOrderPiezo = false;
   if (!secondOrderPiezo)  {
      for (i = 0; i < 3; i++)  {
         peComps.setComponent(i, -au*P(i));
      }
   }
   else  {
      for (i = 0; i < 3; i++)  {
         if (AParams)  {
            peComps.setComponent(i, -au*P(i) - au*P2a(i));
         }  else  {
            peComps.setComponent(i, -au*P(i) - au*P2a(i));
         }
      }
   }

   vP = (FOUR_PI*HA2EV)*((gIntegrate(gBasis|peComps)(0).getComponent(0)
      + gIntegrate(gBasis|peComps)(1).getComponent(1)
      + gIntegrate(gBasis|peComps)(2).getComponent(2)
      ));
/*   if (crystalStruct == 3)  { // for test purposes
      sigmaXZ = sigmaXZRot;
      sigmaYZ = sigmaYZRot;
      sigmaXY = sigmaXYRot;
      sigmaXX = sigmaXXP;
      sigmaYY = sigmaYYP;
      sigmaZZ = sigmaZZP;
   }*/
   Tr = eXX + eYY + eZZ;
}

void SxStrain::initGvec (const SxGBasis &gBasis)
{
   // setup of basis set for u's (3 component psi)
   // and sigma's/P's (1 component psi)

   PsiG temporary;
   gVecs.resize(3);
   gG2.resize(3);
//      PsiG scaledG2;  // common gBasis.g2 does not work once k-scaling is active
   scaledG2.resize(G1.ng);
   for (int idx = 0; idx < G1.ng; idx++)  {
      scaledG2(idx) = sqr((gBasis).gVec.colRef(0)(idx) * wgt(0)) +
                   sqr((gBasis).gVec.colRef(1)(idx) * wgt(1)) +
                   sqr((gBasis).gVec.colRef(2)(idx) * wgt(2));
   }
   scaledG2(0) = 0.0; // to avoid division by zero, meaningless later on as temporary(0)=0
   temporary.resize(G1.ng);
   int gSize = (int)gBasis.getNElements ();
   for (int i = 0; i < 3; i++)  {
      gVecs(i).resize(gSize);
      gG2(i).resize(gSize);
      gVecs(i).setBasis(gBasis);
      gG2(i).setBasis(gBasis);
      for (int iComp = 0; iComp < 3; iComp++)  {
         gVecs(i).setComponent(iComp, gBasis.gVec.colRef(i));
         temporary(0) = 0.;
         for (int idx = 1; idx < G1.ng; idx++)
            temporary(idx) = -G1.gVec(idx,i) / scaledG2(idx);
         gG2(i).setComponent(iComp, temporary);
      }
   }
   /*
   cout << "set up gG2s..." << endl;
   cout << "gSize: " << gSize << endl;
   cout << "scaledG2's length: " << scaledG2.getSize() << endl;
   gG2s.resize(gSize/3);
   gG2s.setBasis(psiG.getBasisPtr());
   gG2s = -(gBasis).gVec / scaledG2;
   for (int idx = 0; idx < gG2s.getSize(); idx++)
      cout << idx << ": " << gG2s(idx) << endl;
   gG2s(0) = 0;
   cout << "success." << endl;
   */
   if (Diag.getSize() == 1)  { // set up diagonal elements for preconditioner

      cout << "set up preconditioner elements... " << endl;
      SxArray<SxArray<PsiR> > sd;
      sd.resize(3);
      for (int i = 0; i < 3; i++)  {
         sd(i).resize(3);
         for (int j = 0; j < 3; j++)
            sd(i)(j) = gVecs(i) * gVecs(j);
      }
      SxArray<PsiG> PDiag(3);

      Diag.resize(gBasis.getNElements ());
      Diag.setBasis (gBasis);
      for (int i = 0; i < 3; i++)
         PDiag(i) = cP(0) * (sd(i)(i)).getComponent(i);
      if (crystalStruct == 2)
         PDiag(2) = cP(4) * (sd(2)(2)).getComponent(2);
      for (int i = 0; i < 3; i++)
         Diag.setComponent(i, PDiag(i));
      cout << "done." << endl;
   }
}

PsiG SxStrain::operator* (const PsiG &psiG) const
{
   const SxGBasis &gBasis = psiG.getBasis<SxGBasis> ();
   // if not yet done: initialise gVecs for first and second derivatives
   // NEEDS A BETTER SOLUTION (?)
   if (gVecs.getSize() == 1)
      const_cast<SxStrain*>(this)->initGvec (gBasis);
   step++;
   PsiG out = (gBasis|gradient(psiG));

   return out;

}


PrecEnergy SxStrain::getEnergy (const SxPsiSet &psiSet,
                                const SxFermi &)
{
   SX_CHECK (dynamic_cast<const SxPWSet *> (&psiSet));
   if (generalized) calculateSigma(uInR(uStore));

   return energy;
}

SxDiracVec<TPrecCoeffG::TReal>
SxStrain::preconditioner (const PsiG &psi,
                                 Preconditioner type) const
{
   // TODO: ugly!!!
   if (type == 0) {/*get rid of compiler warnings*/}


   int nCol = (int)psi.nCols ();
   SX_CHECK (nCol > 0, nCol);
   double kin = 0.;
   SxLaplacian L;
   // average kinetic energy over all states (usually 1)
   for (int i = 0; i < nCol; ++i)
      kin += (psi.colRef (i) | L | psi.colRef(i) );
   kin /= nCol;
//   cout << "kin: " << kin << endl;

   SxDiracVec<TPrecCoeffG::TReal> x, x2, x3, x4, n, K;
   x.resize(psi.getSize());
   x = Diag / ((Diag ^ psi.absSqr()).chop());
   x2 = x.sqr();
   x3 = x.cub();
   x4 = x2.sqr();
   n  = 27. + 18.*x + 12.*x2 + 8.*x3;
   K  = n / (n + 16.*x4);
   K /= kin;

   return K;

}

void SxStrain::writeMeshASCII (const SxString &name, const PsiR &data) const
{

   SxString fileName = name + ".dat";
   ofstream fileN ((fileName).ascii (), fstream::out | fstream::trunc);
   for (int x = 0; x < mesh3d(0); x++)  {
      for (int y = 0; y < mesh3d(1); y++)  {
         for (int z = 0; z < mesh3d(2); z++)  {
            int idx = (int)mesh3d.getMeshIdx(x, y, z, SxMesh3D::Positive);
            fileN << data(idx).re << endl;
         }
      }
   }
   cout << name << " written to " << fileName << endl;
}

void SxStrain::writeMeshSxb (const SxString &name, const SxMeshR &data) const
{
   SxString file = name + ".sxb";
   cout << "writing " << file << endl;;
   SX_CHECK (rBasisPtr);
   const SxCell &cell = rBasisPtr->cell;
   SxBinIO io (file, SxBinIO::BINARY_WRITE_ONLY);
   io.writeMesh (data, cell, mesh3d);
   io.setMode (SxBinIO::WRITE_DATA);
   io.writeMesh (data, cell, mesh3d);
   io.close();
}

void SxStrain::writePotentials () const
{
   if (step > 0) {

      cout << "enter writeRho: " << endl;

      // --- write ascii data
      fstream::openmode mode = fstream::out | fstream::trunc;
      if (eXX.getSize() > 0)  {
         int idx, x, y, z;
         SxString fileName;
         SxString file;
         fileName="data111.dat";
         PsiR sigma;
         cout << endl;
         double val;
         // --- line scan along 111 through center
         {
            ofstream fileN ((fileName).ascii (), mode);
            double eXXVal, eYYVal, eZZVal, eXYVal, eXZVal, eYZVal, trVal;
            int max = mesh3d(0);
            if (mesh3d(1) < max) max = mesh3d(1);
            if (mesh3d(2) < max) max = mesh3d(2);
            for (z = 0; z < max; z++)  {
               idx = (int)mesh3d.getMeshIdx(z, z, z, SxMesh3D::Positive);
               val = vP(idx);
               eXXVal = eXX(idx);
               eYYVal = eYY(idx);
               eZZVal = eZZ(idx);
               eXYVal = eXY(idx);
               eXZVal = eXZ(idx);
               eYZVal = eYZ(idx);
               trVal = eXXVal + eYYVal + eZZVal;

               fileN << SxString(z) << "\t" << SxString(val) << "\t"
                     << SxString (eXXVal) << "\t" << SxString(eYYVal) << "\t"
                     << SxString(eZZVal) << "\t"
                     << SxString (eXYVal) << "\t" << SxString(eXZVal) << "\t"
                     << SxString(eYZVal) << "\t"
                     << SxString (trVal) << endl;
            }
            cout << "strain and Vp along 111 written to " << fileName << endl;
         }
         // --- line scan along 001 through center
         fileName = "data001.dat";
         {
            double eXXVal, eYYVal, eZZVal, eXYVal, eXZVal, eYZVal, trVal;
            ofstream fileN ((fileName).ascii (), mode);
            x = mesh3d(0)/2;
            y = mesh3d(1)/2;
            for (z = 0; z < mesh3d(2); z++)  {
               idx = (int)mesh3d.getMeshIdx(x, y, z, SxMesh3D::Positive);
               val = vP(idx);
               eXXVal = eXX(idx);
               eYYVal = eYY(idx);
               eZZVal = eZZ(idx);
               eXYVal = eXY(idx);
               eXZVal = eXZ(idx);
               eYZVal = eYZ(idx);
               trVal = eXXVal + eYYVal + eZZVal;

               fileN << SxString(z) << "\t" << SxString(val) << "\t"
                     << SxString (eXXVal) << "\t" << SxString(eYYVal) << "\t"
                     << SxString(eZZVal) << "\t"
                     << SxString (eXYVal) << "\t" << SxString(eXZVal) << "\t"
                     << SxString(eYZVal) << "\t"
                     << SxString (trVal) << endl;
            }
            cout << "strain and Vp along 001 written to " << fileName << endl;
         }
         // --- line scan along 100 through center
         fileName = "data100.dat";
         {
            double eXXVal, eYYVal, eZZVal, eXYVal, eXZVal, eYZVal, trVal;
            ofstream fileN ((fileName).ascii (), mode);
            val = -100.;
            y = mesh3d(1)/2;
            z = mesh3d(2)/2;
            for (x = 0; x < mesh3d(0); x++)  {
               idx = (int)mesh3d.getMeshIdx(x, y, z, SxMesh3D::Positive);
               val = vP(idx);
               eXXVal = eXX(idx);
               eYYVal = eYY(idx);
               eZZVal = eZZ(idx);
               eXYVal = eXY(idx);
               eXZVal = eXZ(idx);
               eYZVal = eYZ(idx);
               trVal = eXXVal + eYYVal + eZZVal;

               fileN << SxString(x) << "\t" << SxString(val) << "\t"
                     << SxString (eXXVal) << "\t" << SxString(eYYVal) << "\t"
                     << SxString(eZZVal) << "\t"
                     << SxString (eXYVal) << "\t" << SxString(eXZVal) << "\t"
                     << SxString(eYZVal) << "\t"
                     << SxString (trVal) << endl;
            }
            cout << "strain and Vp along 100 written to " << fileName << endl;
         }
         // --- line scan along 110 through center
         fileName = "data110.dat";
         {
            double eXXVal, eYYVal, eZZVal, eXYVal, eXZVal, eYZVal, trVal;
            ofstream fileN ((fileName).ascii (), mode);
            val = -100.;
            z = mesh3d(2)/2;
            for (x = 0; x < mesh3d(0); x++)  {
               idx = (int)mesh3d.getMeshIdx(x, x, z, SxMesh3D::Positive);
               val = vP(idx);
               eXXVal = eXX(idx);
               eYYVal = eYY(idx);
               eZZVal = eZZ(idx);
               eXYVal = eXY(idx);
               eXZVal = eXZ(idx);
               eYZVal = eYZ(idx);
               trVal = eXXVal + eYYVal + eZZVal;

               fileN << SxString(x) << "\t" << SxString(val) << "\t"
                     << SxString (eXXVal) << "\t" << SxString(eYYVal) << "\t"
                     << SxString(eZZVal) << "\t"
                     << SxString (eXYVal) << "\t" << SxString(eXZVal) << "\t"
                     << SxString(eYZVal) << "\t"
                     << SxString (trVal) << endl;
            }
            cout << "strain and Vp along 110 written to " << fileName << endl;
         }

         if (moreIO)  {
            writeMeshASCII ("strainXX", eXX);
            writeMeshASCII ("strainYY", eYY);
            writeMeshASCII ("strainZZ", eZZ);
            writeMeshASCII ("strainXY", eXY);
            writeMeshASCII ("strainXZ", eXZ);
            writeMeshASCII ("strainYZ", eYZ);
            writeMeshASCII ("vPol", vP);
            writeMeshASCII ("Px", P(0));
            writeMeshASCII ("Py", P(1));
            writeMeshASCII ("Pz", P(2));
         }
      }
      // --- write .sxb strain fields

      cout << "write .sxb files..." << endl;
      if ((eXX.getSize() > 0) && (moreIO)) {
         writeMeshSxb ("epsilonXX", eXX);
         writeMeshSxb ("epsilonYY", eYY);
         writeMeshSxb ("epsilonZZ", eZZ);
         writeMeshSxb ("epsilonXY", eXY);
         writeMeshSxb ("epsilonXZ", eXZ);
         writeMeshSxb ("epsilonYZ", eYZ);

         writeMeshSxb ("vPol", vP);
         cout << "maximum of polarisation potential: " << vP.maxval() << endl;
         cout << "minimum of polarisation potential: " << vP.minval() << endl;
      }
      cout << "Total free energy: " << fTot << " eV. " << "Max: "
           << fMax << " eV. " << "Min: " << fMin << " eV." << endl;
   }
}


