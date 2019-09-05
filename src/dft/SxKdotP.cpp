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


/* Implementation of n band k dot p perturbation theory.
*/

#include <SxKdotP.h>
#include <SxProjector.h>
#include <SxGBasis.h>
#include <SxFSAction.h>
#include <SxFileIO.h>

SxKdotP::SxKdotP (const SxPWSet &waves,
                        const RhoR &rhoRIn)
   : SxHamiltonian (),
     SxRho (rhoRIn),
     wavesPtr (&waves)
// SxKdotP
{
   SX_CHECK (rBasisPtr);
   sdC.resize(3);
   fdC.resize(3);
   
   gBasisPtr = &(rBasisPtr->getGBasis ());
   
   gVecs.resize(1);
   D.resize(1);
   
}

SxArray<PsiR> SxKdotP::firstDerivative (const PsiG &psiG)
{
   const SxRBasis &R = *rBasisPtr;
   SxArray<PsiR> gradR(3);
   int i;

   for (i=0; i < 3; i++)  {
      gradR(i) = I * wgt(i) * (R | ( (-gVecs(i)) * psiG ) );
   }
   return gradR;
}

SxArray<SxArray<PsiR> > SxKdotP::secondDerivative (const PsiG &psiG)
{
   const SxRBasis &R = *rBasisPtr;
   SxArray<SxArray<PsiR> > grad2R(3);

   for (int i = 0; i < 3; i++)  {
      grad2R(i).resize(3);
      grad2R(i)(i) = wgt(i) * wgt(i) *
         ( R | ((  gVecs(i) * gVecs(i)) * psiG) );
   }
      
   grad2R(0)(1) = wgt(0) * wgt(1) *
      ( R | ((  gVecs(0) * (gVecs(1)) * psiG) ) );
   grad2R(0)(2) = wgt(0) * wgt(2) *
      ( R | ((  gVecs(0) * (gVecs(2)) * psiG) ) );
   grad2R(1)(2) = wgt(1) * wgt(2) *
      ( R | ((  gVecs(1) * (gVecs(2)) * psiG) ) );
   
   grad2R(1)(0) = grad2R(0)(1);
   grad2R(2)(0) = grad2R(0)(2);
   grad2R(2)(1) = grad2R(1)(2);
   return grad2R;
}

PrecEnergy SxKdotP::getETrial()
{
   return eTrial;
}

SxMeshR SxKdotP::returnMesh(SxString inString)
{
   if (inString == "epsilonR")
   return epsilonR;
   else return zero;
}

// TODO: This is an ugly routine. Does 2 different things when called, distinguished by
// the counter variable. It is required for optical spectra that rely on the Hamiltonian
// but not used for the electronic structure calculation itself at all.
// Problem is that all communication with the Hamiltonian from anywhere else
// (e.g. a tool) is done via the HamPtr, so whatsoever is required needs a virtual
// function in SxHamiltonian.h. Maybe it is reasonable to allow for some arbitrary
// communication channel in SxHamiltonian that can be used to transfer various
// data between a Hamiltonian and other tools & add-ons?
void SxKdotP::formDerivatives()
{
   if (counter == 0)  {
      accurateInterfaces = false;
      cout << "Switch to simplified interfaces for optical spectra..." << endl;
      counter ++; 
   } else
   {
      cout << "form k-derivative of Hamiltonian" << endl;
      derivatives = true;
      // --- replace with templates to ensure no parameters are substituted twice
      if (inFile.contains("k2") > 0) inFile = inFile.substitute("k2", "###K2###");
      if (inFile.contains("kx2") > 0) inFile = inFile.substitute("kx2", "###KX2###");
      if (inFile.contains("ky2") > 0) inFile = inFile.substitute("ky2", "###KY2###");
      if (inFile.contains("kz2") > 0) inFile = inFile.substitute("kz2", "###KZ2###");
      if (inFile.contains("kxy") > 0) inFile = inFile.substitute("kxy", "###KXY###");
      if (inFile.contains("kxz") > 0) inFile = inFile.substitute("kxz", "###KXZ###");
      if (inFile.contains("kyz") > 0) inFile = inFile.substitute("kyz", "###KYZ###");
      if (inFile.contains("kx") > 0) inFile = inFile.substitute("kx", "###KX###");
      if (inFile.contains("ky") > 0) inFile = inFile.substitute("ky", "###KY###");
      if (inFile.contains("kz") > 0) inFile = inFile.substitute("kz", "###KZ###");
      if (inFile.contains("k") > 0) inFile = inFile.substitute("k", "###K###");
      // --- replace with k-derivatives, multiplied by light polarisation a
      if (inFile.contains("###K2###")) inFile = inFile.substitute("###K2###", "(2*kx*ax+2*ky*ay+2*kz*az)");
      if (inFile.contains("###KX2###")) inFile = inFile.substitute("###KX2###", "2*kx*ax");
      if (inFile.contains("###KY2###")) inFile = inFile.substitute("###KY2###", "2*ky*ay");
      if (inFile.contains("###KZ2###")) inFile = inFile.substitute("###KZ2###", "2*kz*az");
      if (inFile.contains("###KXY###")) inFile = inFile.substitute("###KXY###", "(ky*ax+kx*ay)");
      if (inFile.contains("###KXZ###")) inFile = inFile.substitute("###KXZ###", "(kz*ax+kx*az)");
      if (inFile.contains("###KYZ###")) inFile = inFile.substitute("###KYZ###", "(kz*ay+ky*az)");
      if (inFile.contains("###K###")) inFile = inFile.substitute("###K###", "(ax+ay+az)");
      if (inFile.contains("###KX###")) inFile = inFile.substitute("###KX###", "(ax)");
      if (inFile.contains("###KY###")) inFile = inFile.substitute("###KY###", "(ay)");
      if (inFile.contains("###KZ###")) inFile = inFile.substitute("###KZ###", "(az)");
      cout << "build updated tree: " << endl;
      SxString hamString = ((inFile.right("Hamiltonian")).right("=[")).left("];");
      SxList<SxString> cols = hamString.tokenize(']');
      int i, j;
      for (i = 0; i < nComp; i++)  {
         cols(i) = cols(i).right("[");
         SxList<SxString> rows = cols(i).tokenize(',');
         int nRows = (int)rows.getSize();
         // --- split columns in single elements
         for (j = 0; j < nRows; j++)  {
            expression(i)(j).resize(0);
            // --- evaluate single element
            buildTree(i, j, rows(j), inFile);
         }
      }
      cout << "success." << endl;
   }
}

void SxKdotP::update (const SxFermi &fermi)
{
   fermiPtr = &fermi;
   // --- external potential
   eExt = 0.;
}

void SxKdotP::compute (const SxFermi &fermi,
                       bool tauChanged,
                       bool rhoChanged)
{
   if (tauChanged && rhoChanged) {/* just to get rid of compiler-warning*/};
   update(fermi);
}

PsiG SxKdotP::operator* (const PsiG &psiG)
{
   // --- test: show <psiI|H|psiJ>
   SX_CHECK (psiG.getBasisPtr());

   const SxGBasis &G = *dynamic_cast<const SxGBasis *>(psiG.getBasisPtr());

   // initialize gVecs for first and second derivatives
   if (gVecs.getSize() == 1)  {
      gVecs.resize(3);
      for (int i = 0; i < 3; i++)  {
         gVecs(i).resize(psiG.getSize());
         gVecs(i).setBasis(psiG.getBasisPtr());
         for (int iComp = 0; iComp < nComp; iComp++)
            gVecs(i).setComponent(iComp, G.gVec.colRef(i));
      }
      // setup operator -iG/G^2 to compute potential from external charges

      PsiG scaledG2;

//      gG2.resize(3);
      scaledG2.resize(psiG.getSize());
      for (int iComp = 0; iComp < nComp; iComp++)
         for (int idx = 0; idx < psiG.getSize()/nComp; idx++)  {
            scaledG2(idx + iComp * psiG.getSize() / nComp)
                       = sqr(G.gVec.colRef(0)(idx) * wgt(0)) +
                         sqr(G.gVec.colRef(1)(idx) * wgt(1)) +
                         sqr(G.gVec.colRef(2)(idx) * wgt(2));
         }
      scaledG2(0) = 0.0; // to avoid division by zero, meaningless later on as temporary(0)=0
/*      temporary.resize(psiG.getSize() / nComp);
      for (int i = 0; i < 3; i++)  {
//         gG2(i).resize(psiG.getSize());
//         gG2(i).setBasis(psiG.getBasisPtr());
         for (int iComp = 0; iComp < nComp; iComp++)  {
            for (int idx = 0; idx < psiG.getSize() / nComp; idx++){
               temporary(idx)
                  = -(G).gVec.colRef(i)(idx) / scaledG2(idx);
            }
            temporary(0) = 0.;
//            gG2(i).setComponent(iComp, temporary);
         }
      }*/
      if (nCharges > 0)  {
         cout << "compute potential from external charges" << endl;
         vChg.resize(rSize*nComp);
         double NM = 18.897261;
         const SxRBasis &R = *rBasisPtr;
         SxMatrix3<Double> cell = R.cell;
         SxMesh3D mesh = R.getMesh();
         double d = 1;
         for (int i = 0; i < 3; i++)
            d = d * (cell(i) / (mesh(i) * NM))(i);
         PsiR chgR = (R | psiG);//.getComponent(0);
         PsiG vG;
         if (totalCharge.getSize() > 0)  {
            for (int iComp = 0; iComp < nComp; iComp++)  {
            chgR.setComponent(iComp,totalCharge/epsilonR); // ugly...
            }
//            cout << "gsf" << endl;
            vG = (G | (chgR))
                 * wgt(0) * wgt(1) * wgt(2) / (d*NM);//*scaledG2);
            for (int idx = 0; idx < vG.getSize(); idx++)
                if (scaledG2(idx).re != 0) {vG(idx) = vG(idx)/scaledG2(idx);}
                else {vG(idx) = 0;} 
//            cout << "vG's size: " << vG.getSize() << endl;
//            cout << "scaledG2's size: " << scaledG2.getSize() << endl;
            vG(0) = 0;
            vChg = (R | vG).getComponent(0);
            chargePotential.resize(vChg.getSize());
            for (int idx = 0; idx < vChg.getSize(); idx++)
               chargePotential(idx) = vChg(idx)/HA2EV;
            cout << "Charge Potential: " << chargePotential.minval() * HA2EV
                 << " eV to " << chargePotential.maxval() * HA2EV << " eV."
                 << endl << "Potential difference: " 
                 << (chargePotential.maxval() - chargePotential.minval() ) * HA2EV
                 << " eV." << endl;
            cout << "write charge potential " << endl;
            
            ofstream vChgFile;
            fstream::openmode mode = fstream::out | fstream::trunc;
            vChgFile.open (("vChg.dat"), mode);
            double valChg;
            for (int x = 0; x < mesh(0); x++)  {
               for (int y = 0; y < mesh(1); y++)  {
                  for (int z = 0; z < mesh(2); z++)  {
                     ssize_t idx = mesh.getMeshIdx(x,y,z, SxMesh3D::Positive);
                     valChg = chargePotential(idx) * HA2EV;
                     vChgFile << SxString(valChg) << endl;
                  }
               }
            }
            vChgFile.close ();
            cout << "done. " << endl;
         }
      }
   }
   return (G | hamNxN (psiG));
}

void SxKdotP::writeMeshASCII (const SxString &name, const PsiR &data) const
{
   
   const SxRBasis &R = *rBasisPtr;
   SxMesh3D mesh = R.getMesh();
   SxString fileName = name + ".dat";
   ofstream fileN ((fileName).ascii (), fstream::out | fstream::trunc);
   for (int x = 0; x < mesh(0); x++)  {
      for (int y = 0; y < mesh(1); y++)  {
         for (int z = 0; z < mesh(2); z++)  {
            int idx = (int)mesh.getMeshIdx(x, y, z, SxMesh3D::Positive);
            fileN << data(idx).re << endl;
         }
      }
   }
   cout << name << " written to " << fileName << endl;
}

SxMatrix<Complex16> SxKdotP::hMatrixBS(SxVector3<Double> kVal, SxVector3<Int> r)
{
   SxComplex16 e;
   const SxRBasis &R = *rBasisPtr;
   SxMesh3D mesh = R.getMesh();
   int coordIdx = (int)mesh.getMeshIdx(r, SxMesh3D::Positive);

   SxMatrix<Complex16> Ham(nComp, nComp);
   SxMatrix<Complex16 > hamiltonian(nComp, nComp);
   int iComp, jComp;
   for (iComp = 0; iComp < nComp; iComp++)  {
      for (jComp = 0; jComp < nComp; jComp++)  {
         e = evaluateTreeBS(iComp, jComp, 0, coordIdx, kVal);
         hamiltonian(iComp, jComp) = e;
      }
   }

   return hamiltonian;
}

// --- provide bandstructure between given k's
void SxKdotP::showBandstructure (SxString outFile, SxVector3<Int> r)
{

   int iStep, iBand;
   SxMatrix<Complex16> hMat;
   SxBinIO::deleteFile(outFile);
/*   SxString outFile2d = "bs2d.dat";
   SxBinIO::deleteFile(outFile2d);
   cout << "produce 3d band properties at k = 0..." << endl;
   const SxRBasis &R = *rBasisPtr;
   SxMesh3D mesh = R.getMesh();
   SxVector3<Int> pos;
   SxVector3<Double> null;
   null(0) = 0; null(1) = 0; null(2) = 0;
   pos(2) = r(2);
   for ( x = 0; x < mesh(0); x++)  {
      cout << "x: " << x << endl;
      pos(0) = x;
      for ( y = 0; y < mesh(1); y++)  {
         pos(1) = y;
         idx = mesh.getMeshIdx(pos, SxMesh3D::Positive);
         hMat = hMatrixBS(null, pos);
         SxString outLine = "";
         for (iBand = 0; iBand < nComp; iBand++)
            outLine = outLine + "   "
             + SxString((hMat.eigenvalues())(iBand).re * HA2EV);
         outLine = outLine + "\n";
         outLine.appendToFile(outFile2d);
      }
   }
   cout << "done." << endl;
*/         
   try  {
      SxFileIO outio;
      outio.open (outFile, "a");
      for (iStep = 0; iStep < bsPts.getSize(); iStep++) {
         hMat = hMatrixBS(bsPts(iStep), r);
      
         SxString outLine = "" + SxString(bsPts(iStep)(0)) + "  " + 
            SxString(bsPts(iStep)(1)) + "  " + SxString(bsPts(iStep)(2));
         for (iBand = 0; iBand < nComp; iBand++)
            outLine = outLine + "   "
               + SxString((hMat.eigenvalues())(iBand).re * HA2EV);
         outLine = outLine + "\n";
         outio.write (outLine);
      }
      outio.close ();
   } catch (const SxException &e) {
      e.print ();
      SX_EXIT;
   }
   cout << "+------- Band structure calculation finished. -------" << endl;
   if (onlyBS) {
      cout << "Only band structure calculation requested.\n"
           << "Exiting now." << endl; 
      exit(0); 
   }
}

SxArray<PsiR> SxKdotP::rFirstDerivative(const PsiR &partial)  {
   PsiG psiG = (*wavesPtr)(0,0,0);
   int i, iComp;
   const SxRBasis &R = *rBasisPtr;
   PsiR psiRFull = (R | psiG);
   const SxGBasis &G = *dynamic_cast<const SxGBasis *>(psiG.getBasisPtr());
   PsiR pTemp;
   SxArray<PsiR> fdParTmp, fdPar;
   fdPar.resize(3);
   pTemp = 0.* psiRFull;
   for (iComp = 0; iComp < nComp; iComp++)
      pTemp.setComponent(iComp,partial);
   fdParTmp = firstDerivative(G|pTemp);
   fdPar(0) = fdParTmp(0).getComponent(0);
   fdPar(1) = fdParTmp(1).getComponent(0);
   fdPar(2) = fdParTmp(2).getComponent(0);
   for (i = 0; i < 3; i++)  {
      fdPar(i) = -I*fdParTmp(i).getComponent(0);
   } 
   return fdPar;
}

// --- the Hamiltonian itself
PsiR SxKdotP::hamNxN (const PsiG &psiG)
{
   SX_CLOCK(Timer::NxNham);
   SX_CHECK (psiG.getBasisPtr());
   const SxRBasis &R = *rBasisPtr;
   PsiR psiRFull = (R | psiG);
   PsiR Hcol = 0. * psiRFull;
   PsiR H;
   H.resize(Hcol.getSize());
   int iComp, jComp, i, j;
   PsiR elem;
   bool precond = false;
   SxArray<PsiR> fd; // first derivatives dPsiG/d(x,y,z)
   SxArray<SxArray<PsiR> > sd; // second derivatives d/d(x,y,z)(dPsiG/d(x,y,z)

   // initialize preconditioner and potential from input charges
   if (D.getSize() == 1)  {
      SX_CLOCK(Timer::init);
      gSize = (int)psiG.getComponent(0).getSize();
      gZero.resize(gSize);
      gOne.resize(gSize);
      int idx;
      for (idx = 0; idx < gSize; idx++)  {
         gZero(idx) = 0.;
         gOne(idx) = 1.;
      }
      // --- resize fdC and sdC, now both can be used by hamElem(i,j,prec)
      for (i = 0; i < 3; i++)  {
         fdC(i).resize(gSize);
         for (j = 0; j < 3; j++)
            sdC(i)(j).resize(gSize);
      }

      SxArray<SxArray<PsiG> > sdP; // second derivatives in preconditioner
      precond = true;
      D.resize(psiG.getSize());
      D = 0. * psiG;
      sdP.resize(3);
      for (i = 0; i < 3; i++)  {
         sdP(i).resize(3);
         for (j = 0; j < 3; j++)
            sdP(i)(j) = gVecs(i) * gVecs(j);
      }
      for (iComp = 0; iComp < nComp; iComp++)  {
         for (i = 0; i < 3; i++)  {
            for (j = 0; j < 3; j++)
            if (i == j) {sdC(i)(j) <<= sdP(i)(j).getComponent(iComp);}
            else {sdC(i)(j) <<= 0.*sdP(i)(j).getComponent(iComp);}

         }
         elem = evaluateTree(iComp,iComp, 0, precond).abs();
         gSize = (int)D.getSize()/nComp;
         for (i = 0; i < gSize; i++) {
         D(i + iComp*gSize) = elem(i);
         }

      }
      precond = false;

      // --- resize fdC and sdC, now both can be used by hamElem(i,j,prec)
      for (i = 0; i < 3; i++)  {
         fdC(i).resize(rSize);
         for (j = 0; j < 3; j++)  {
            sdC(i)(j).resize(rSize);
         }
      }

   }   // timings: 0.0%

   {
   SX_CLOCK(Timer::fstDer);
   fd = firstDerivative(psiG);
   }   // timing: 7.8%
   {
   SX_CLOCK(Timer::sndDer);
   sd = secondDerivative(psiG);
   }   // timing: 15.2 %
   iMult = 0;
   for (iComp = 0; iComp < nComp; iComp++)  {
      psiR = psiRFull.getComponent(iComp);
      for (i = 0; i < 3; i++)  {
         SX_CLOCK(Timer::derivatives);
         fdC(i) = -I * fd(i).getComponent(iComp);
         for (j = 0; j < 3; j++) {
            sdC(i)(j) = sd(i)(j).getComponent(iComp);
         }
      }

      for (jComp = 0; jComp < nComp; jComp++)  {
         {

         SX_CLOCK(Timer::evalTree);
         elem = evaluateTree(jComp,iComp, 0, precond);
         }
         Hcol.setComponent(jComp, elem);
      }
      if (iComp == 0)
         H <<= Hcol;
      else H += Hcol;
   }
   firstStep = false;
   return H;
}

const SxSymbolTable *
SxKdotP::getHamiltonianGroup (const SxSymbolTable *table)
{
   SX_CHECK (table);
   SxSymbolTable *hGroup = NULL;
   try  { hGroup = table->getGroup("kpHamiltonian"); }
   catch (const SxException&) { /*EXIT;*/ }
   return hGroup;
}

int SxKdotP::getNEmptyStates () {
   return nEmptyStates;
}

void SxKdotP::validateNStates (int nStatesIn)
{
   int nStatesMax = wavesPtr->getGkBasis().nGkMin;

   if (nStatesIn > nStatesMax)  {
      sxprintf ("Error: input parameter(s) \"nEmptyStates\" or "
              "\"nEmptyStatesChi\" in (PW)Hamil-\n"
              "tonian group chosen too large. Please reduce it/them or "
              "try with a\n"
              "larger cutoff \"eCut\". The number of states you want to "
              "compute must\n"
              "       not exceed %d.\n", nStatesMax);
      SX_QUIT;
   }
}

SxMeshR SxKdotP::parameters(int iParam)
{
   SxMeshR p;
   int iMat;
   if (!speedOpt)  {
      p = matParam(0, iParam) * materials(0);
      if (bowParam(0, iParam) != 0)
            p -= bowParam(0, iParam)
               * (1. - materials(0)) * materials(0);
      for (iMat = 1; iMat < nMat; iMat++)  {
         p = p + matParam(iMat, iParam) * materials(iMat);
         if (bowParam(iMat, iParam) != 0)
            p -= bowParam(iMat, iParam)
               * (1. - materials(iMat)) * materials(iMat);
      }
      return p;
   }
   else
      return pMem(iParam);
}

int SxKdotP::firstOutsideBracket(SxString expr, char op)
{
   int pos = -1;
   int bLevel = 0;
   int i;

   for (i = 0; i < expr.getSize(); i++)  {
      if (expr(i) == '(') bLevel++;
      if (expr(i) == ')') bLevel--;
      if ((expr(i) == op) && (bLevel == 0) && (pos == -1))
         pos = i;
   }

   return pos;
}

SxString SxKdotP::containsUnknown(SxList<SxString> inList)
{
   SxString unknownElem = "###";
   int iElem;
   for (iElem = 0; iElem < inList.getSize(); iElem++)  {
      SxString elem = inList(iElem);
      if (
            ( // --- is one of the known operators
                (elem == "k2") ||
                (elem == "kx2") ||
                (elem == "ky2") ||
                (elem == "kz2") ||

                (elem == "kxy") ||
                (elem == "kxz") ||
                (elem == "kyz") ||

                (elem == "k") ||
                (elem == "kx") ||
                (elem == "ky") ||
                (elem == "kz") ||

                (elem == "e") ||
                (elem == "eXX") ||
                (elem == "eYY") ||
                (elem == "eZZ") ||

                (elem == "eXY") ||
                (elem == "eXZ") ||
                (elem == "eYZ") ||

                (elem == "Vp") ||
                (elem == "Vchg") ||
                (elem == "Vext")

         )  || // --- or known parameter
            (paramNames.contains(elem))
            || // --- or complex number
            (elem == "i")
            ||
            (
             (elem.substitute("i","").substitute("I","")
              ).isDouble()
             )
            || // --- or function Sqrt{}
            (
             (elem.trim().substitute("_{", "", 1)
               .substitute("}", "", 1)
              ).isDouble()
            )
            || // --- or function Sin{}
            (
             (elem.trim().substitute("!{", "", 1)
               .substitute("}", "", 1)
              ).isDouble()
            )
            || // --- or function Cos{}
            (
             (elem.trim().substitute("?{", "", 1)
               .substitute("}", "", 1)
              ).isDouble()
            )
            || // --- or function Tan{}
            (
             (elem.trim().substitute(".{", "", 1)
               .substitute("}", "", 1)
              ).isDouble()
            )
      )  {
      } else  {
         unknownElem = elem;
      }
   }   
  return unknownElem;
}

SxString SxKdotP::replaceUnknown(SxString elem, SxString ham)
{
   elem = elem.substitute(" ", "");
   SxString elemOut = elem;
   /* Blanks are required to make sure that unknown
    strings are not replaced as parts of other strings.
    Blanks have to be removed afterwards again.*/
   elemOut = elemOut.substitute("+", " + ");
   elemOut = elemOut.substitute("-", " - ");
   elemOut = elemOut.substitute("*", " * ");
   elemOut = elemOut.substitute("/", " / ");
   elemOut = elemOut.substitute("(", " ( ");
   elemOut = elemOut.substitute(")", " ) ");

   // --- split element in single contributions
   elem = elem.substitute("+", "|");
   elem = elem.substitute("-", "|");
   elem = elem.substitute("*", "|");
   elem = elem.substitute("/", "|");
   elem = elem.substitute("(", "|");
   elem = elem.substitute(")", "|");

   // --- set up list of single contributions
   SxList<SxString> split = elem.tokenize('|');
   
   SxString unknownElem = containsUnknown(split);

   if (unknownElem != "###")  { // --- unknown element detected
      SxString replacement;
      SxString search = unknownElem + "=";
      int srch = (int)ham.find(search);
      if ( (srch > -1 ) && (srch < ham.getSize()) )  {
         // --- check wether all terms A = ... contain semicolon
         replacement = "(" + ham.right(unknownElem + "=").left(";") + ")";
         if (replacement.contains("=") > 0) {
             cout << "Error: '" << unknownElem << "=' lacks semicolon."
                  << endl;
             SX_EXIT;
         }
         elemOut = elemOut.substitute(unknownElem, replacement);
      }
      else  {
      cout << "Error: element " << unknownElem
         << " is not found anywhere. EXITING."<< endl;
      cout << paramNames << endl;
      SX_EXIT;
      }
      elemOut = replaceUnknown(elemOut.trim(), ham);

   }
   return elemOut;
}

int SxKdotP::resolveExpr(SxString expr, int col, int row)
{
   int pos = -1;
   // --- check for +
   SxString left, right;
   int l = 0, r  = 0, rtn = 0;
   pos = firstOutsideBracket(expr, '+');
   if (pos > -1)  {
      left = expr.subString(0, pos - 1);
      right = expr.subString(pos + 1);
      expression(col)(row) << "+";
      rtn = (int)expression(col)(row).getSize() - 1;
      leftPtr(col)(row).resize(rtn + 1);
      rightPtr(col)(row).resize(rtn + 1);
      l = resolveExpr(left, col, row);
      r = resolveExpr(right, col, row);
      leftPtr(col)(row)(rtn) = l;
      rightPtr(col)(row)(rtn) = r;
   } 
   // --- check for -
   if (pos == -1) {
      pos = firstOutsideBracket(expr, '-');
      if (pos > -1)  {
         left = expr.subString(0, pos - 1);
         if (left.trim() == "") left = "0";
         right = expr.subString(pos + 1);
         expression(col)(row) << "-";
         rtn = (int)expression(col)(row).getSize() - 1;
         leftPtr(col)(row).resize(rtn + 1);
         rightPtr(col)(row).resize(rtn + 1);
         l = resolveExpr(left, col, row);
         r = resolveExpr(right, col, row);
         leftPtr(col)(row)(rtn) = l;
         rightPtr(col)(row)(rtn) = r;
      } 
   }
   // --- check for *
   if (pos == -1) {
      pos = firstOutsideBracket(expr, '*');
      if (pos > -1)  {
         left = expr.subString(0, pos - 1);
         right = expr.subString(pos + 1);
         expression(col)(row) << "*";
         rtn = (int)expression(col)(row).getSize() - 1;
         leftPtr(col)(row).resize(rtn + 1);
         rightPtr(col)(row).resize(rtn + 1);
         l = resolveExpr(left, col, row);
         r = resolveExpr(right, col, row);
         leftPtr(col)(row)(rtn) = l;
         rightPtr(col)(row)(rtn) = r;
      }
   }
   // --- check for /
   if (pos == -1) {
      pos = firstOutsideBracket(expr, '/');
      if (pos > -1)  {
         left = expr.subString(0, pos - 1);
         right = expr.subString(pos + 1);
         expression(col)(row) << "/";
         rtn = (int)expression(col)(row).getSize() - 1;
         leftPtr(col)(row).resize(rtn + 1);
         rightPtr(col)(row).resize(rtn + 1);
         l = resolveExpr(left, col, row);
         r = resolveExpr(right, col, row);
         leftPtr(col)(row)(rtn) = l;
         rightPtr(col)(row)(rtn) = r;
      } 
   }
   // if none of the above is found: search for brackets
   if (pos == -1)  {
      int lBrck = -1, rBrck = -1; // left bracket, right bracket's position
      int length = (int)expr.getSize();
      int i;
      for (i = 0; i < length; i++)  {
         if ((lBrck == -1) && (expr(i) == '('))
            lBrck = i;
      }
      for (i = length-1; i > 0; i--)  {
         if ((rBrck == -1) && (expr(i) == ')'))
            rBrck = i;
      }
      if ((lBrck > -1) && (rBrck > -1))  {
         expression(col)(row) << "()";
         rtn = (int)expression(col)(row).getSize() - 1;
         leftPtr(col)(row).resize(rtn + 1);
         rightPtr(col)(row).resize(rtn + 1);

         leftPtr(col)(row)(rtn) = rtn + 1;
         rightPtr(col)(row)(rtn) = -1;
         resolveExpr(expr.subString(lBrck + 1, rBrck - 1), col, row);
         pos = 0;
      } 
   }
   // --- if still no element identified: expr is a leaf
   if (pos == -1) {
      expression(col)(row) << expr.trim();
      rtn = (int)expression(col)(row).getSize() - 1;
      leftPtr(col)(row).resize(rtn + 1);
      rightPtr(col)(row).resize(rtn + 1);
      leftPtr(col)(row)(rtn) = -1;
      rightPtr(col)(row)(rtn) = -1;
   }
   return rtn;
}

SxComplex16 SxKdotP::evaluateTreeBS(int col, int row, int pos, int coordIdx,
      SxVector3<Double> kVec)
{
   SxComplex16 rtn = 0;
   SxComplex16 kx = kVec(0);
   SxComplex16 ky = kVec(1);
   SxComplex16 kz = kVec(2);

   SxString expr = expression(col)(row)(pos).substitute("!", "");
   int pIdx;
   // --- evaluate +, -, *, /
   if (expr == "+")
      rtn = evaluateTreeBS(col, row, leftPtr(col)(row)(pos), coordIdx, kVec)
         + evaluateTreeBS(col, row, rightPtr(col)(row)(pos), coordIdx, kVec);
   else if (expr == "-")
      rtn = evaluateTreeBS(col, row, leftPtr(col)(row)(pos), coordIdx, kVec)
         - evaluateTreeBS(col, row, rightPtr(col)(row)(pos), coordIdx, kVec);
   else if (expr == "*")
      rtn = evaluateTreeBS(col, row, leftPtr(col)(row)(pos), coordIdx, kVec)
         * evaluateTreeBS(col, row, rightPtr(col)(row)(pos), coordIdx, kVec);
   else if (expr == "/")
      rtn = evaluateTreeBS(col, row, leftPtr(col)(row)(pos), coordIdx, kVec)
         / evaluateTreeBS(col, row, rightPtr(col)(row)(pos), coordIdx, kVec);
   // --- evaluate ()
   else if (expr == "()")
      rtn = evaluateTreeBS(col, row, leftPtr(col)(row)(pos), coordIdx, kVec);

   // --- evaluate leafs

   // --- k^2 - operators, also mixed op's
   else if (expr == "k2") rtn = kx * kx + ky * ky + kz * kz;
   else if (expr == "kx2") rtn = kx * kx;
   else if (expr == "ky2") rtn = ky * ky;
   else if (expr == "kz2") rtn = kz * kz;
   else if (expr == "kxy") rtn = kx * ky;
   else if (expr == "kxz") rtn = kx * kz;
   else if (expr == "kyz") rtn = ky * kz;

   // --- linear k - operators

   else if (expr == "k") rtn = kx + ky + kz;
   else if (expr == "kx") rtn = kx;
   else if (expr == "ky") rtn = ky;
   else if (expr == "kz") rtn = kz;

   // --- strains

   else if (expr == "e") rtn = eIJ(0)(coordIdx)
      + eIJ(1)(coordIdx) + eIJ(2)(coordIdx);
   else if (expr == "eXX") rtn = eIJ(0)(coordIdx);
   else if (expr == "eYY") rtn = eIJ(1)(coordIdx);
   else if (expr == "eZZ") rtn = eIJ(2)(coordIdx);
   else if (expr == "eXY") rtn = eIJ(3)(coordIdx);
   else if (expr == "eXZ") rtn = eIJ(4)(coordIdx);
   else if (expr == "eYZ") rtn = eIJ(5)(coordIdx);

   // --- external or polarization potential
   else if (expr == "Vp") rtn = vP(coordIdx);
   else if (expr == "Vext") rtn = vExt(coordIdx);
   else if (expr == "Vchg") rtn = chargePotential(coordIdx);

   // --- function
   else if (expr.contains("_{") > 0)  {
      if ((expr.right("_{").left("}")).isDouble())  {
         rtn = sqrt((expr.right("_{").left("}")).toDouble());
      }
   }
   else if (expr.contains("!{") > 0)  {
      if ((expr.right("!{").left("}")).isDouble())  {
         rtn = sin((expr.right("!{").left("}")).toDouble());
      }
   }
   else if (expr.contains("?{") > 0)  {
      if ((expr.right("?{").left("}")).isDouble())  {
         rtn = cos((expr.right("?{").left("}")).toDouble());
      }
   }
   else if (expr.contains(".{") > 0)  {
      if ((expr.right(".{").left("}")).isDouble())  {
         rtn = tan((expr.right(".{").left("}")).toDouble());
      }
   }
   // --- complex factor
   else if ((expr.contains("i") > 0) && (expr.left("i").isDouble()))  {
      rtn = I * expr.left("i").toDouble();
   }
   else if (expr == "i")  {
      rtn = I;
   }
   else if (expr.isDouble()) rtn = expr.toDouble();

   // --- known parameter
   else  {
      for (pIdx = 0; pIdx < paramNames.getSize(); pIdx++)
         if (expr == paramNames(pIdx)) rtn = parameters(pIdx)(coordIdx);
   }

   return rtn;
}

int SxKdotP::resolveOpKey(SxString opString)
{
   if (opString == "k") return 0;
   else if (opString == "kx") return 1;
   else if (opString == "ky") return 2;
   else if (opString == "kz") return 3;
   else if (opString == "k2") return 4;
   else if (opString == "kx2") return 5;
   else if (opString == "ky2") return 6;
   else if (opString == "kz2") return 7;
   else if (opString == "kxy") return 8;
   else if (opString == "kxz") return 9;
   else if (opString == "kyz") return 10;
   else return -1;
}

PsiR SxKdotP::returnAccurate(int i)
{
   PsiR rtn = zero;
   SX_CLOCK(Timer::kk2s);
   if ((opKey(i)) == 0)  { // k = kx + ky + kz
      rtn = par(i) * (fdC(0) + fdC(1) + fdC(2))
          + 0.5 * kIPar(i) * psiR;
   }
   if ((opKey(i)) == 1)  { // kx
      rtn = (par(i) * fdC(0) + 0.5 * kIPar(i) * psiR);
   }
   if ((opKey(i)) == 2)  { // ky
      rtn = (par(i) * fdC(1) + 0.5 * kIPar(i) * psiR);
   }
   if ((opKey(i)) == 3)  { // kz
      rtn = (par(i) * fdC(2) + 0.5 * kIPar(i) * psiR);
   }

   if ((opKey(i)) == 4)  { // k2 = kx2 + ky2 + kz2
      rtn = par(i) * (sdC(0)(0) + sdC(1)(1) + sdC(2)(2))
          + (kIPar(i) * fdC(0) + kJPar(i) * fdC(1) + kKPar(i) * fdC(2));
   }
   if ((opKey(i) == 5))  { // kx2
      rtn = par(i) * sdC(0)(0) + kIPar(i) * fdC(0);
   }
   if ((opKey(i) == 6))  { // ky2
      rtn = par(i) * sdC(1)(1) + kIPar(i) * fdC(1);
   }
   if ((opKey(i) == 7))  { // kz2
      rtn = par(i) * sdC(2)(2) + kIPar(i) * fdC(2);
   }
   if ((opKey(i) == 8))  { // kxy = kx*ky
      rtn = par(i) * sdC(0)(1) + 0.5 * (kIPar(i) * fdC(1) + kJPar(i) * fdC(0));
   }

   if ((opKey(i) == 9))  { // kxz = kx * kz
      rtn = par(i) * sdC(0)(2) + 0.5 * (kIPar(i) * fdC(2) + kJPar(i) * fdC(0));
   }

   if ((opKey(i) == 10))  { // kyz = ky * kz
      rtn = par(i) * sdC(1)(2) + 0.5 * (kIPar(i) * fdC(2) + kJPar(i) * fdC(1));
   }
   return rtn;
}

void SxKdotP::determineDerivatives(int i)
{
   if (opKey(i) == 0)  {// k
      kIPar(i) = rFirstDerivative(par(i))(0)
            + rFirstDerivative(par(i))(1)
            + rFirstDerivative(par(i))(2);
      kJPar(i) = zero;
      kIJPar(i) = zero;
      }
   if (opKey(i) == 1)  {// kx
      kIPar(i) = rFirstDerivative(par(i))(0);
      kJPar(i) = zero;
      kIJPar(i) = zero;
      }
   if (opKey(i) == 2)  {// ky
      kIPar(i) = rFirstDerivative(par(i))(1);
      kJPar(i) = zero;
      kIJPar(i) = zero;
      }
   if (opKey(i) == 3)  {// kz
      kIPar(i) = rFirstDerivative(par(i))(2);
      kJPar(i) = zero;
      kIJPar(i) = zero;
      }
   if (opKey(i) == 4)  {// k2
      kIPar(i) = rFirstDerivative(par(i))(0);
      kJPar(i) = rFirstDerivative(par(i))(1);
      kKPar(i) = rFirstDerivative(par(i))(2);
      kIJPar(i) = rFirstDerivative(rFirstDerivative(par(i))(0))(0)
                + rFirstDerivative(rFirstDerivative(par(i))(1))(1)
                + rFirstDerivative(rFirstDerivative(par(i))(2))(2);
/*      if (help.getSize() == 1)  {
      cout << "initialize help..." << endl;
      help.resize(rSize);
      help2.resize(rSize);
      help     = kIJPar(i);
      help2    = par(i);
      }*/
      }
   if (opKey(i) == 5)  {// kx2
      kIPar(i) = rFirstDerivative(par(i))(0);
      kJPar(i) = kIPar(i);
      kIJPar(i) = rFirstDerivative(rFirstDerivative(par(i))(0))(0);
      }
   if (opKey(i) == 6)  {// ky2
      kIPar(i) = rFirstDerivative(par(i))(1);
      kJPar(i) = kIPar(i);
      kIJPar(i) = rFirstDerivative(rFirstDerivative(par(i))(1))(1);
      }
   if (opKey(i) == 7)  {// kz2
      kIPar(i) = rFirstDerivative(par(i))(2);
      kJPar(i) = kIPar(i);
      kIJPar(i) = rFirstDerivative(rFirstDerivative(par(i))(2))(2);
      }
   if (opKey(i) == 8)  {// kxy
      kIPar(i) = rFirstDerivative(par(i))(0);
      kJPar(i) = rFirstDerivative(par(i))(1);
      kIJPar(i) = rFirstDerivative(rFirstDerivative(par(i))(0))(1);
      }
   if (opKey(i) == 9)  {// kxz
      kIPar(i) = rFirstDerivative(par(i))(0);
      kJPar(i) = rFirstDerivative(par(i))(2);
      kIJPar(i) = rFirstDerivative(rFirstDerivative(par(i))(0))(2);
      }
   if (opKey(i) == 10)  {// kyz
      kIPar(i) = rFirstDerivative(par(i))(1);
      kJPar(i) = rFirstDerivative(par(i))(2);
      kIJPar(i) = rFirstDerivative(rFirstDerivative(par(i))(1))(2);
      }
}

void SxKdotP::printElement(int pos, int col, int row)
{
   SxString expr = expression(col)(row)(pos);
   if (expr == "()")  {
      cout << "(";
      printElement(leftPtr(col)(row)(pos), col, row);
      cout << ")";
   }
   else if (expr == "*")  {
      printElement(leftPtr(col)(row)(pos), col, row);
      cout << " * ";
      printElement(rightPtr(col)(row)(pos), col, row);
   }
   else if (expr == "/")  {
      printElement(leftPtr(col)(row)(pos), col, row);
      cout << " / ";
      printElement(rightPtr(col)(row)(pos), col, row);
   }
   else if (expr == "+")  {
      printElement(leftPtr(col)(row)(pos), col, row);
      cout << " + ";
      printElement(rightPtr(col)(row)(pos), col, row);
   }
   else if (expr == "-")  {
      printElement(leftPtr(col)(row)(pos), col, row);
      cout << " - ";
      printElement(rightPtr(col)(row)(pos), col, row);
   } else {
     cout << expr;
   }

}

PsiR SxKdotP::evaluateTree(int col, int row, int pos, bool prec)
{
   /* This is the main part of the input file parser.
      Though it always did an excellent job for my purposes, I am aware that
      this parser can certainly be made more efficient and also I am quite
      sure there are certain redundancies that could be removed. However, as
      it was tested for a wide range of systems and Hamiltonians, I hesitate
      to do any changes to it myself.
      Whoever wants to make the parser more efficient, I will at least
      provide a couple of hints & ideas:

   - The Hamiltonian itself is nowhere calculated. Instead, H|Psi> is
     what is commonly calculated. Means that the parser somehow has to
     distinguish between operator-like multiplications (parameter * kx2
          or 0.5 * parameter * (kx + ky) or so ) and potential-like
          contributions such as Ecb or Vp or Vext. The latter ones only need
     multiplication with |Psi>, whereas the others require operators to
     be applied.

   - Please be aware that the product rule commonly applies to products
     involving a real-space parameter and operators (this is the
     'accurateInterfaces'-tag in the code). As an example, the effective
     mass kinetic contribution reads:
       (d/dr) (1/me) (d/dr)|Psi> rather than (1/me) * (d^2/dr^2) |Psi>
     This made it so far impossible to use sums of operators directly
     in a bracket, i.e. in the current version you need to use
       (parameter * kx2 + parameter * ky2) rather than
             parameter * (kx2 + ky2). Bit cumbersome, but it does the job.
     I fully agree it can certainly be implemented better.

   - So far, complex operations such as Sqrt{} are hardly available. In
     particular, the Sqrt{} function can so far be applied only to real
     numbers. However, it might be useful to have this one also for some
     parameters. Likewise, one might consider implementing Sqr{} as well.
     I have never found any use for Sin{}, Cos{} and Tan{} and such, but
     if there is need for it, I suppose it can be done quite fast.

      If you feel you want to donate a better parser, I hope this helps. For
      further questions, do not hesitate to ask: o.marquardt@yahoo.ie.     
   */
   PsiR rtn = zero;
   if (!prec)  { // --- Hamiltonian call not used for preconditioner
      int pIdx;
      // --- in case no operator is present in a summand, multiply with psiR
      bool noOperator = false;
      SxString expr = expression(col)(row)(pos);
      if (expr.contains("!"))  {
         noOperator = true;
         expr = expr.left("!");
      }
      // --- evaluate +, -, *, /
      if (expr == "+") {
         rtn = evaluateTree(col, row, leftPtr(col)(row)(pos), prec)
            + evaluateTree(col, row, rightPtr(col)(row)(pos), prec);
         }
      else if (expr == "-") {
         rtn = evaluateTree(col, row, leftPtr(col)(row)(pos), prec)
            - evaluateTree(col, row, rightPtr(col)(row)(pos), prec);}
      else if (expr == "*")  {
         /* in case the accurateInterface option is set, the simplified
            interface treatment d/dr (P(r) * Psi(r)) = P(r) * dPsi(r)/dr
            is replaced by the correct one:
            d/dr (P(r) * Psi(r)) = dP(r)/dr * Psi(r) + dPsi(r)/dr * P(r)
         */
         if ( firstStep && accurateInterfaces &&
         ( ( (expression(col)(row)(leftPtr(col)(row)(pos)) == "()")
         && (!isRealNumber(rightPtr(col)(row)(pos), col, row ) )
         && (containsOperator(leftPtr(col)(row)(pos), col, row) ) )
         || ( (expression(col)(row)(rightPtr(col)(row)(pos)) == "()")
         && (!isRealNumber(leftPtr(col)(row)(pos), col, row ) )
         && (containsOperator(rightPtr(col)(row)(pos), col, row) ) ) ) ) {
            cout << "Error in element: " << col << ", " << row << ": ";
            printElement(pos,col,row);
            cout << endl;
            cout << "Multiplication with bracket that contains operator "
                 << "detected.\nThis construct is not allowed in "
                 << "accurateInterfaces mode. EXITING." << endl;
            SX_EXIT;
         }
          
         if (  accurateInterfaces
               && ((expression(col)(row)(leftPtr(col)(row)(pos)).contains("k"))
               || (expression(col)(row)(rightPtr(col)(row)(pos)).contains("k")))
            )  {
            if (expression(col)(row)(leftPtr(col)(row)(pos)).contains("k"))  {
               if (opMultStr(iMult) == "~")  { // if not yet defined
                  opMultStr(iMult) = "%"; // defined now
                  // --- resolve operator key
                  opKey(iMult)
                    = resolveOpKey(
                      expression(col)(row)(leftPtr(col)(row)(pos)));
                  // --- identify parameter P(r)
                  par(iMult)
                    = evaluateTree(col, row, rightPtr(col)(row)(pos), prec);
                  determineDerivatives(iMult);
               }
               rtn = returnAccurate(iMult);
               iMult++;
            } else 
            if (expression(col)(row)(rightPtr(col)(row)(pos)).contains("k"))  {
               if (opMultStr(iMult) == "~")  { // if not yet defined
                  opMultStr(iMult) = "%"; // defined now
                  // --- resolve operator key
                  opKey(iMult)
                    = resolveOpKey(
                    expression(col)(row)(rightPtr(col)(row)(pos)));
                  // --- identify parameter P(r)
                  par(iMult)
                    = evaluateTree(col, row, leftPtr(col)(row)(pos), prec);
                  determineDerivatives(iMult);
               }
               rtn = returnAccurate(iMult);   // timing: 16.0
               iMult++;
            }
          
         } else {
            rtn = evaluateTree(col, row, leftPtr(col)(row)(pos), prec)
               * evaluateTree(col, row, rightPtr(col)(row)(pos), prec);
            }
         }
      else if (expr == "/") {

         rtn = evaluateTree(col, row, leftPtr(col)(row)(pos), prec)
            / evaluateTree(col, row, rightPtr(col)(row)(pos), prec);
         }
      // --- evaluate ()
      else if (expr == "()") {
         rtn = evaluateTree(col, row, leftPtr(col)(row)(pos), prec);
      }
      // --- evaluate leafs

      // --- k^2 - operators, also mixed op's
/*      else if (expr == "k2")  rtn = sdC(0)(0) + sdC(1)(1) + sdC(2)(2);
      else if (expr == "kx2") rtn = sdC(0)(0);
      else if (expr == "ky2") rtn = sdC(1)(1);
      else if (expr == "kz2") rtn = sdC(2)(2);
      else if (expr == "kxy") rtn = sdC(0)(1);
      else if (expr == "kxz") rtn = sdC(0)(2);
      else if (expr == "kyz") rtn = sdC(1)(2);
      // --- linear k - operators

      else if (expr == "k") rtn = fdC(0) + fdC(1) + fdC(2);
      else if (expr == "kx") rtn = fdC(0);
      else if (expr == "ky") rtn = fdC(1);
      else if (expr == "kz") rtn = fdC(2);*/

      // --- strains

      else if (expr == "e")   {SX_CLOCK(Timer::extField); rtn = eIJ(0) + eIJ(1) + eIJ(2);}
      else if (expr == "eXX") {SX_CLOCK(Timer::extField); rtn = eIJ(0);}
      else if (expr == "eYY") {SX_CLOCK(Timer::extField); rtn = eIJ(1);}
      else if (expr == "eZZ") {SX_CLOCK(Timer::extField); rtn = eIJ(2);}
      else if (expr == "eXY") {SX_CLOCK(Timer::extField); rtn = eIJ(3);}
      else if (expr == "eXZ") {SX_CLOCK(Timer::extField); rtn = eIJ(4);}
      else if (expr == "eYZ") {SX_CLOCK(Timer::extField); rtn = eIJ(5);}

      // --- external, charge or polarization potential
      else if (expr == "Vp") {SX_CLOCK(Timer::extField); rtn = vP;}
      else if (expr == "Vext") {SX_CLOCK(Timer::extField); rtn = vExt;}
      else if (expr == "Vchg") {SX_CLOCK(Timer::extField); rtn = chargePotential;}

      // --- function
      else if (expr.contains("_{") > 0)  {
         if ((expr.right("_{").left("}")).isDouble())  {
            rtn = sqrt((expr.right("_{").left("}")).toDouble()) * one;
         }
      }
      else if (expr.contains("!{") > 0)  {
         if ((expr.right("!{").left("}")).isDouble())  {
            rtn = sin((expr.right("!{").left("}")).toDouble()) * one;
         }
      }
      else if (expr.contains("?{") > 0)  {
         if ((expr.right("?{").left("}")).isDouble())  {
            rtn = cos((expr.right("?{").left("}")).toDouble()) * one;
         }
      }
      else if (expr.contains(".{") > 0)  {
         if ((expr.right(".{").left("}")).isDouble())  {
            rtn = tan((expr.right(".{").left("}")).toDouble()) * one;
         }
      }

      // --- complex factor
      else if ((expr.contains("i") > 0) && (expr.left("i").isDouble()))  {
         SX_CLOCK(Timer::complexNr);
         rtn = I * expr.left("i").toDouble() * one;
      }   // timing: 0.5
      else if (expr == "i")  {
         SX_CLOCK(Timer::complexNr);
         rtn = I * one;
      }   // timing: 0.5
      else if (expr.isDouble()) {SX_CLOCK(Timer::doubleNr); rtn = expr.toDouble() * one;}   // timing: 3.7

      // --- known parameter
      else  {
         SX_CLOCK(Timer::known);
         for (pIdx = 0; pIdx < paramNames.getSize(); pIdx++)
            if (expr == paramNames(pIdx)) rtn = parameters(pIdx);
         }   //timing: 9.7

      // --- multiply constants with psiR
      if (noOperator)  {
         SX_CLOCK(Timer::constant);
         rtn = rtn * psiR;
         if (derivatives) rtn = zero;
      }   // timing: 1.0
   } // --- end (!prec)-loop  timing: 55.5%
   else  { // --- if used for preconditioner
      if (col != row)
         rtn = gZero; // zero preconditioner element in reciprocal space
      else  {
         rtn = gZero;
         int pIdx;
         // --- in case no operator present in a summand, multiply with psiR
         bool noOperator = false;
         SxString expr = expression(col)(row)(pos);
         if (expr.contains("!"))  {
            noOperator = true;
            expr = expr.left("!");
         }
         // --- evaluate +, -, *, /
         if (expr == "+")
            rtn = evaluateTree(col, row, leftPtr(col)(row)(pos), prec)
               + evaluateTree(col, row, rightPtr(col)(row)(pos), prec);
         else if (expr == "-")
            rtn = evaluateTree(col, row, leftPtr(col)(row)(pos), prec)
               - evaluateTree(col, row, rightPtr(col)(row)(pos), prec);
         else if (expr == "*")
            rtn = evaluateTree(col, row, leftPtr(col)(row)(pos), prec)
               * evaluateTree(col, row, rightPtr(col)(row)(pos), prec);
         else if (expr == "/")
            rtn = evaluateTree(col, row, leftPtr(col)(row)(pos), prec)
               / evaluateTree(col, row, rightPtr(col)(row)(pos), prec);
         // --- evaluate ()
         else if (expr == "()") {
              rtn = evaluateTree(col, row, leftPtr(col)(row)(pos), prec); 
            }
         // --- evaluate leafs

         // --- k^2 - operators, also mixed op's
         else if (expr == "k2") rtn = sdC(0)(0) + sdC(1)(1) + sdC(2)(2) ;
         else if (expr == "kx2") rtn = sdC(0)(0);
         else if (expr == "ky2") rtn = sdC(1)(1);
         else if (expr == "kz2") rtn = sdC(2)(2);
         else if (expr == "kxy") rtn = sdC(0)(1);
         else if (expr == "kxz") {rtn = sdC(0)(2);

         }
         else if (expr == "kyz") rtn = sdC(1)(2);

         // --- linear k - operators

         else if (expr == "k") rtn = fdC(0) + fdC(1) + fdC(2);
         else if (expr == "kx") rtn = fdC(0);
         else if (expr == "ky") rtn = fdC(1);
         else if (expr == "kz") rtn = fdC(2);

         // --- strains - no contribution to preconditioner

         else if (expr == "e") rtn = gZero;         
         else if (expr == "eXX") rtn = gZero;
         else if (expr == "eYY") rtn = gZero;
         else if (expr == "eZZ") rtn = gZero;
         else if (expr == "eXY") rtn = gZero;
         else if (expr == "eXZ") rtn = gZero;
         else if (expr == "eYZ") rtn = gZero;

         // --- additional potentials - no contribution to preconditioner
         else if (expr == "Vp") rtn = gZero;
         else if (expr == "Vext") rtn = gZero;
         else if (expr == "Vchg") rtn = gZero;

         // --- complex factor
         else if (expr.contains("_{") > 0)  {
            if ((expr.right("_{").left("}")).isDouble())  {
               rtn = sqrt((expr.right("_{").left("}")).toDouble()) * gOne;
            }
         }
         else if (expr.contains("!{") > 0)  {
            if ((expr.right("!{").left("}")).isDouble())  {
               rtn = sin((expr.right("!{").left("}")).toDouble()) * gOne;
            }
         }
         else if (expr.contains("?{") > 0)  {
            if ((expr.right("?{").left("}")).isDouble())  {
               rtn = cos((expr.right("?{").left("}")).toDouble()) * gOne;
            }
         }
         else if (expr.contains(".{") > 0)  {
            if ((expr.right(".{").left("}")).isDouble())  {
               rtn = tan((expr.right(".{").left("}")).toDouble()) * gOne;
            }
         }
         else if ((expr.contains("i") > 0) && (expr.left("i").isDouble()))  {
            rtn = I * expr.left("i").toDouble() * gOne;
         }
         else if (expr == "i")  
            rtn = I * gOne;
         else if (expr.isDouble()) rtn = expr.toDouble() * gOne;

         // --- known parameter
         else 
            for (pIdx = 0; pIdx < paramNames.getSize(); pIdx++) {
               if (expr == paramNames(pIdx)) 
                  rtn = matParam(precMaterial,pIdx) * gOne;
               }

         // --- in preconditioner: if no operator is present, element = 0
         if (noOperator)
            rtn = gZero;
      }

   rtn.resize(gSize);
   }   // timing: 0.0%
//   for (int idx=0; idx < rtn.getSize(); idx++)  {
//      cout << idx << ": " << rtn(idx) << endl;
//   }
   return rtn;
}

bool SxKdotP::isRealNumber(int pos, int col, int row)
{
bool rtn = false;
   SxString expr = expression(col)(row)(pos);
   if (expr == "()")  {
      rtn = isRealNumber(leftPtr(col)(row)(pos), col, row);
   }
   else if (expr == "*")  {
      rtn = isRealNumber(leftPtr(col)(row)(pos), col, row) 
      && isRealNumber(rightPtr(col)(row)(pos), col, row);
   }
   else if (expr == "/")  {
      rtn = isRealNumber(leftPtr(col)(row)(pos), col, row) 
      && isRealNumber(rightPtr(col)(row)(pos), col, row);
   }
   else if (expr == "+")  {
      rtn = isRealNumber(leftPtr(col)(row)(pos), col, row) 
      && isRealNumber(rightPtr(col)(row)(pos), col, row);
   }
   else if (expr == "-")  {
      rtn = isRealNumber(leftPtr(col)(row)(pos), col, row) 
      && isRealNumber(rightPtr(col)(row)(pos), col, row);
   } else if (
     ((expr.substitute("i","")).substitute("_{",""))
      .substitute("}","").isDouble())  {
      rtn = true;
   }

return rtn;
}

bool SxKdotP::containsOperator(int pos, int col, int row)
{
   bool rtn = false;
   SxString expr = expression(col)(row)(pos);
//   cout << "expr: " << expr << endl;
   if (expr == "()")  {
      rtn = containsOperator(leftPtr(col)(row)(pos), col, row);
   }
   if (expr == "*")  {
      rtn = containsOperator(leftPtr(col)(row)(pos), col, row) 
      || containsOperator(rightPtr(col)(row)(pos), col, row);
   }
   if (expr == "/")  {
      rtn = containsOperator(leftPtr(col)(row)(pos), col, row) 
      || containsOperator(rightPtr(col)(row)(pos), col, row);
   }
   if (expr == "+")  {
      rtn = containsOperator(leftPtr(col)(row)(pos), col, row) 
      || containsOperator(rightPtr(col)(row)(pos), col, row);
   }
   if (expr == "-")  {
      rtn = containsOperator(leftPtr(col)(row)(pos), col, row) 
      || containsOperator(rightPtr(col)(row)(pos), col, row);
   }
   if ((expr == "k")
   || (expr == "kx")
   || (expr == "ky")
   || (expr == "kz")
   || (expr == "k2")
   || (expr == "kx2")
   || (expr == "ky2")
   || (expr == "kz2")
   || (expr == "kxy")
   || (expr == "kxz")
   || (expr == "kyz"))  {
      rtn = true;
   };
   return rtn;
}

int SxKdotP::whatIsElement(int pos, int col, int row)
{
   int out = -1; // 0: operator, 1: potential, 2: not clear

   /* bracket:
      If element in bracket is operator, mark bracket with "?"
      and remove "?" from element in bracket to prevent double
      operator countings.
      If element in bracket is a potential, mark bracket with "!"
      and remove "!" from element in bracket to prevent double
      multiplication with psi.
    */
   if (expression(col)(row)(pos) == "()")  {
      out = whatIsElement(leftPtr(col)(row)(pos), col, row);
      if ((out == 0))  {
         expression(col)(row)(pos) += "?";
         expression(col)(row)(leftPtr(col)(row)(pos)) =
         expression(col)(row)(leftPtr(col)(row)(pos))
            .substitute("?", "", 1);
/*    // DEBUG 11.5.2015
         expression(col)(row)(rightPtr(col)(row)(pos)) =
         expression(col)(row)(rightPtr(col)(row)(pos))
            .substitute("?", "", 1);*/
         out = 0;
      }
      else  {
         if (expression(col)(row)(leftPtr(col)(row)(pos)).contains("!"))  {
            expression(col)(row)(pos) += "!";
            expression(col)(row)(leftPtr(col)(row)(pos)) =
               expression(col)(row)(leftPtr(col)(row)(pos))
               .substitute ("!", "", 1);
            out = 1;
         }
         else  {
            out = 0;
         }
      }
   }
   
   /* multiplication, division:
      If one of the elements is an operator, mark expression with "?"
      and remove "?" and "!" from left and right branches.
      If left and right branches are potentials, mark expression with "!"
      and remove "!" from left and right branches.
      */
   else if ((expression(col)(row)(pos) == "*")
         || (expression(col)(row)(pos) == "/"))  {
      // --- multiplication yields operator
      if ((whatIsElement(leftPtr(col)(row)(pos), col, row) == 0)
         || (whatIsElement(rightPtr(col)(row)(pos), col, row) == 0))  {
         out = 0;
         expression(col)(row)(pos) += "?";
         expression(col)(row)(pos) = expression(col)(row)(pos)
            .substitute("!", "");
         expression(col)(row)(leftPtr(col)(row)(pos)) = 
         (expression(col)(row)(leftPtr(col)(row)(pos))
            .substitute ("?", "", 1)).substitute("!", "", 1);
         expression(col)(row)(rightPtr(col)(row)(pos)) = 
         (expression(col)(row)(rightPtr(col)(row)(pos))
            .substitute ("?", "", 1)).substitute("!", "", 1);
      }
      // --- multiplication yields potential
      else if ((whatIsElement(leftPtr(col)(row)(pos), col, row) == 1)
         || (whatIsElement(rightPtr(col)(row)(pos), col, row) == 1))  {
         out = 2;
         if((whatIsElement(leftPtr(col)(row)(pos), col, row) == 1)
               && (whatIsElement(rightPtr(col)(row)(pos), col, row) == 1))  {
            out = 1;
            expression(col)(row)(pos) += "!";
         }
         expression(col)(row)(leftPtr(col)(row)(pos)) = 
         expression(col)(row)(leftPtr(col)(row)(pos))
          .substitute ("!", "", 1);
         expression(col)(row)(rightPtr(col)(row)(pos)) = 
         expression(col)(row)(rightPtr(col)(row)(pos))
          .substitute ("!", "", 1);
      }
      else out = 2;
   }
   /*
      addition, subtraction:
      If left and right branches are operators, mark expression with "?"
      and remove "?" from left and right branches.
      If left and right branches are potentials, mark expression with "!"
      and remove "!" from left and right branches.
      If left and right branches are not of the same kind, do nothing
      and leave left and right branch marks unaffected.
      */
   else if ((expression(col)(row)(pos) == "+")
         || (expression(col)(row)(pos) == "-"))  {
      if ((whatIsElement(leftPtr(col)(row)(pos), col, row) == 0)
         && (whatIsElement(rightPtr(col)(row)(pos), col, row) == 0))  {
         out = 0;
         expression(col)(row)(pos) += "?";
         expression(col)(row)(leftPtr(col)(row)(pos)) = 
         expression(col)(row)(leftPtr(col)(row)(pos))
         .substitute ("?", "", 1);
         expression(col)(row)(rightPtr(col)(row)(pos)) = 
         expression(col)(row)(rightPtr(col)(row)(pos))
         .substitute ("?", "", 1);
      }
      else if ((whatIsElement(leftPtr(col)(row)(pos), col, row) == 1)
         && (whatIsElement(rightPtr(col)(row)(pos), col, row) == 1)) {
         if (
               (expression(col)(row)(leftPtr(col)(row)(pos)).contains("!"))
               && 
               (expression(col)(row)(rightPtr(col)(row)(pos)).contains("!"))
            ) {
            out = 1;
            expression(col)(row)(pos) += "!";
            expression(col)(row)(leftPtr(col)(row)(pos)) = 
               expression(col)(row)(leftPtr(col)(row)(pos))
               .substitute ("!", "", 1);
            expression(col)(row)(rightPtr(col)(row)(pos)) = 
               expression(col)(row)(rightPtr(col)(row)(pos))
               .substitute ("!", "", 1);
      }}
      else {
         out = 2;
         expression(col)(row)(pos) = expression(col)(row)(pos)
            .substitute("!", "", 1)
            .substitute("?", "", 1);
      }

   }
   else if (
         (expression(col)(row)(pos) == "k2") ||
         (expression(col)(row)(pos) == "kx2") ||
         (expression(col)(row)(pos) == "ky2") ||
         (expression(col)(row)(pos) == "kz2") ||
         (expression(col)(row)(pos) == "kxy") ||
         (expression(col)(row)(pos) == "kxz") ||
         (expression(col)(row)(pos) == "kyz") ||
         (expression(col)(row)(pos) == "k") ||
         (expression(col)(row)(pos) == "kx") ||
         (expression(col)(row)(pos) == "ky") ||
         (expression(col)(row)(pos) == "kz") 
            )  {
      out = 0;
      expression(col)(row)(pos) += "?";

   }
      else (out = 1);
   if ((out == 1) && (!expression(col)(row)(pos).contains("!"))
         && (!expression(col)(row)(pos).contains("?"))
         && (pos > 0))  {
      expression(col)(row)(pos) += "!";

   }

   if (expression(col)(row)(pos).contains("?")) out = 0;
   if (expression(col)(row)(pos).contains("!")) out = 1;
   return out;

}

void SxKdotP::removeOpMarks(int pos, int col, int row)
{

   expression(col)(row)(pos) = expression(col)(row)(pos)
      .substitute("?", "");
   if (leftPtr(col)(row)(pos) > 0)
      removeOpMarks(leftPtr(col)(row)(pos), col, row);
   if (rightPtr(col)(row)(pos) > 0)
      removeOpMarks(rightPtr(col)(row)(pos), col, row);
}

void SxKdotP::correctOperators(int col, int row)
{
   whatIsElement(0, col, row);
   removeOpMarks(0, col, row);
}

void SxKdotP::countOperatorMultiplications(SxString in)
{
   int length = (int)in.getSize();
   for (int i = 0; i < length-1; i++)
      if ((in(i) == '*') && (in(i+1) == 'k'))
   nOpMult++;
}

void SxKdotP::buildTree (int col, int row, SxString elem, SxString ham)
{
   // --- substitute all unknown variables in Hamiltonian element, if possible
   SxString elemFull = replaceUnknown(elem, ham);
   // --- only in accurateInterface mode: check for operator in brackets
   if (accurateInterfaces)  {
      countOperatorMultiplications(elemFull.substitute(" ", ""));
   }
   // --- routine that actually sets up the tree
   resolveExpr(elemFull, col, row);
   /* Problem: operators kx, ky, kz automatically apply to Psi.
    Constants need multiplication with Psi if not somehow multiplied
    with an operator.*/
   correctOperators(col,row);
   // --- print out tree, UNCOMMENT FOR DEBUGGING PURPOSES
/*   cout << "List of expressions, tree " << col << ", " << row << ": " << endl;
   for (int i = 0; i < expression(col)(row).getSize(); i++)  {
      cout << i << ":\t" << expression(col)(row)(i);
      if (leftPtr(col)(row)(i) > -1)
         cout << "\t" << leftPtr(col)(row)(i);
      if (rightPtr(col)(row)(i) > -1)
         cout << "\t" << rightPtr(col)(row)(i);
     cout <<endl;
   }*/
//   SX_EXIT;
}

void SxKdotP::read (const SxSymbolTable *table)
{

   accurateInterfaces = true;
   derivatives = false;
   firstStep = true;
   counter = 0;
   // --- read the material parameter files ---
   nMat = 0;
   nParam = 0;
   speedOpt = false;
   moreIO   = true;
   precMaterial = 0;
   SxSymbolTable *paramSet = NULL, *param = NULL, *materialMap = NULL,
                 *hamiltonian = NULL, *strain = NULL, *bandstructure = NULL,
                 *path = NULL, *extCharge = NULL;

   // --- loop over all materials parameter sets
   SxMap<SxString, int>::Iterator it;
   iEpsR = -1;
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
            if ((name == "k2") 
               || (name == "kx2")
               || (name == "ky2")
               || (name == "kz2")
               || (name == "kxy")
               || (name == "kxz")
               || (name == "kyz")
               || (name == "k")
               || (name == "kx")
               || (name == "ky")
               || (name == "kz")
               || (name == "kz")
               || (name == "e")
               || (name == "eXX")
               || (name == "eYY")
               || (name == "eZZ")
               || (name == "eXY")
               || (name == "eXZ")
               || (name == "eYZ")
               || (name == "Vp")
               || (name == "Vchg")
               || (name == "Vext"))  {
               cout << "parameter name is reserved keyword " << name
                    << "! EXITING." << endl;
               SX_EXIT;
            }
       
            paramMap(name) = nParam;
            paramNames << name;
            if (name == "epsilon") iEpsR = nParam;
            nParam++;
         }
      }

      SxString material = paramSet->get("material")->toString ();
      matNames << material;
      if (paramSet->contains("useForPreconditioner"))
         precMaterial = nMat;

      // --- browse parameter set for all parameter names
      for (it = paramMap.begin(); it != paramMap.end(); it++)  {
         bool found = false; // parameter found in material parameter group?
         param = paramSet->getGroup("parameter");
         while ((!found) && (param != NULL)) {
            SxString name = param->get("name")->toString();
            if (it.getKey() == name)
            {
               double val = param->get("value")->toReal();
               paramVals << val;
               found = true;
               double bowing = 0.;
               if (param->contains("bowing"))  {
                  cout << "bowing employed" << endl;
                  bowing = param->get("bowing")->toReal();
               }
               bowingVals << bowing;

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
   
   // --- read external charge file names and prefactors //TODO:NEW FEATURE
   nCharges = 0;
   short dataType = -1;
   SxMeshR chargeTmp;

   if (table->containsGroup("charge"))  {
      for (extCharge  = table->getGroup("charge");
            extCharge != NULL;
            extCharge  = extCharge->nextSibling ("charge"))

      {
         SxString filename = extCharge->get("file")->toString ();
         double chgPrefactor = extCharge->get("prefactor")->toReal ();
         chargeFiles << filename;
         chargePrefactors << chgPrefactor;
         nCharges++;
      }
   }

   matParam = SxMatrix<Double> (nMat, nParam, paramVals);
   bowParam = SxMatrix<Double> (nMat, nParam, bowingVals);

   cout << "####### Input materials and parameters ########" << endl;
   cout << "Mat.\t  |  ";
   for (int iMat = 0; iMat < nMat; iMat++)
      cout << matNames(iMat) << "\t| ";
   cout << endl;
   for (int iParam = 0; iParam < nParam; iParam++)  {
      cout << paramNames(iParam) << " \t| ";
      for (int iMat = 0; iMat < nMat; iMat++)  {
         cout << matParam(iMat, iParam) << " \t| ";
      }
      cout << "\n";
   }

   // --- read the material map file ---
   const SxRBasis &R = *rBasisPtr;
   SxString inputASCII = "";
   SxString inputBinary = "";
   SxMesh3D mesh = R.getMesh();
   
   int xMax = mesh(0);
   int yMax = mesh(1);
   int zMax = mesh(2);
   int x, y, z, i, j;
   const int meshSize = xMax * yMax * zMax;
   float value;
   rSize = xMax * yMax * zMax;
   help.resize(1);

   for (i = 0; i < 3; i++)  {
      sdC(i).resize(3);
      fdC(i).resize(rSize);
      for (j = 0; j < 3; j++)
         sdC(i)(j).resize(rSize);
   }
   
   materialMap  = table->getGroup("materialMap");
   
   if (materialMap->contains("asciiFile"))  {
      inputASCII = materialMap->get("asciiFile")->toString();
      dataType = 1;
   } else if (materialMap->contains("binaryFile"))  {
      inputBinary = materialMap->get("binaryFile")->toString();
      dataType = 0;
   }

   hamiltonian  = table->getGroup("kpHamiltonian");
   if (hamiltonian->contains("simplifiedInterfaces"))
      accurateInterfaces = false;
   // --- optimize for speed?
   if (hamiltonian->contains("speedOpt"))
      speedOpt = true;

   if (hamiltonian->contains("lessIO"))
      moreIO = false;

   materials.resize(nMat);
   for (int iMat = 0; iMat < nMat; iMat++)
      materials(iMat).resize(rSize);
   psiR.resize(rSize);
   zero.resize(rSize);
   one.resize(rSize);
   zero.set (0.);
   one.set (1.);
   int ok = 0; // nr of mesh points that have a total composition of 1
   FILE *file =0;
   SxMatrix3<Double> cell = R.cell;
   if (dataType == 0)  {
      SxBinIO io (inputBinary, SxBinIO::BINARY_READ_ONLY);
      materials = io.readMesh (&cell, &mesh);
      io.close();
   } else if (dataType == 1)  {
      int line = 0;
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
               line++;
               double checksum = 0.;
               for (int iMat = 0; iMat < nMat; iMat++)  {
                  fscanf(file,"%f\t", &value);
                  materials(iMat)(meshIdx) = value;
                  checksum += value;
               }
               if (fabs(checksum - 1.) < 1.e-5) ok++;
               else if (fabs(checksum - 1.) > 1.e-5)
                  cout << "checksum error at (" << x << ", " << y << ", "
                  << z << "), diff = " << fabs(1. - checksum) << ", line: "
                  << line << endl;
            }
         }
      }
   }
   // --- check if size of map fits mesh dimensions
   if (materials(0).getSize() != meshSize)  {
      cout << "Mesh size does not fit map size! EXITING." << endl;
      SX_EXIT;
   }
   if (speedOpt)  {
      cout << "optimized for speed -> higher memory consumption" << endl;
      pMem.resize(nParam);
      speedOpt = false; // allows to use SxMeshR SxKdotP::parameters(iParam)
      for (int iParam = 0; iParam < nParam; iParam++)  {
         pMem(iParam) = parameters(iParam);
      }
      speedOpt = true;
   }
   else  {
      cout << "optimized for memory -> higher computational time" << endl;
   }

   int iComp;

   // --- read energy where to search for states
   eTrial = hamiltonian->get("eTrial")->toReal();
   int nBands = hamiltonian->get("nBands")->toInt();
   cout << "eTrial = " << eTrial << " Hartree, nBands = " << nBands << endl;

   cout << "open Hamiltonian file ";
   SxString hamFile;
      hamFile = hamiltonian->get("hamFile")->toString();
   cout << hamFile << "... ";
   ifstream hfile(hamFile.ascii());
   if (!hfile)   {
      cout << endl << "file "
           << hamFile.ascii() << " does not exist. EXITING." << endl;
      SX_EXIT;
   }
   try { file = fopen(hamFile.ascii(), "r+");} catch (const SxException&) { };
   cout << "success." << endl;
   cout << "Build Hamiltonian tree...";
   char c;
   SxString hamString;
   inFile = "";

   while (!feof(file))  {
      c = (char)fgetc(file);
      if ((c != '\n') && (c != '\t') && (c != ' '))
         inFile = inFile + c;
   }
   inFile = inFile.subString(0, inFile.getSize() - 2);
   cout << "success." << endl;

   SxString comment;

   while (inFile.contains("/*") > 0)  {
      if (inFile.right("/*").contains("*/") < 1)  {// no end of comment
            cout << "Comment not closed in Hamiltonian! Exiting." << endl;
            SX_EXIT;
         }
      comment = inFile.right("/*").left("*/");
      inFile = inFile.substitute("/*" + comment + "*/", "", 1);
   }

   inFile = inFile.substitute("Sqrt", "_");
   inFile = inFile.substitute("Sin", "!");
   inFile = inFile.substitute("Cos", "?");
   inFile = inFile.substitute("Tan", ".");
   // --- allow treatment of A - B - C terms
   inFile = inFile.substitute("-", "+(-1.)*"); 
   inFile = inFile.substitute("=+", "="); 
   inFile = inFile.substitute(",+", ","); 
   inFile = inFile.substitute("[+", "[");
   // --- check wether parameter names are overwritten
   for (int pIdx = 0; pIdx < paramNames.getSize(); pIdx++)  {
      if (inFile.contains(paramNames(pIdx) + "=") > 0) {
         cout << "Error: Parameter " << paramNames(pIdx)
              << " must not be overwritten.\nPlease rename the parameter "
              << "in either the parameter file or the Hamiltonian." << endl;
         SX_QUIT;
      }
   }
   // ------------------------------------------------

   hamString = ((inFile.right("Hamiltonian")).right("=[")).left("];");
   if (nBands != int(sqrt(1. * (hamString.contains(",") + 1))))  {
      cout << "Number of bands is not consistent with provided Hamiltonian. "
           << endl << "nBands = " << nBands << endl
           << "Elements in Hamiltonian: " << hamString.contains(",") + 1
           << ", should be " << nBands * nBands << ". EXITING." << endl;
      SX_EXIT;
   }
   
   // --- substitute keywords by symbols
   SxList<SxString> cols = hamString.tokenize(']');
   int nCols = (int)cols.getSize();
   // --- split Hamiltonian in columns
   nComp = nCols;
   cout << nComp << "-band Hamiltonian used here..." << endl;
   // --- set up matrix with Hamiltonian elements
   expression.resize(nComp);
   leftPtr.resize(nComp);
   rightPtr.resize(nComp);
   for (iComp = 0; iComp < nComp; iComp++)  {
      expression(iComp).resize(nComp);
      leftPtr(iComp).resize(nComp);
      rightPtr(iComp).resize(nComp);
   }
   nOpMult = 0;
   cout << "cols: " << nCols << endl;
   for (i = 0; i < nCols; i++)  {
      cols(i) = cols(i).right("[");
      SxList<SxString> rows = cols(i).tokenize(',');
      int nRows = (int)rows.getSize();
      if ( i == 0) cout << "rows: " << nRows << endl;
      if (nRows != nCols) {
         cout << "Hamiltonian is no square matrix. EXITING." << endl;
         cout << "nRows: " << nRows << endl;
         cout << "nCols: " << nCols << endl;
         SX_EXIT;
      }
      // --- split columns in single elements
      for (j = 0; j < nRows; j++)  {
         // --- evaluate single element
         buildTree(i, j, rows(j), inFile);
      }
   }
   cout << "tree set up." << endl;
   if (accurateInterfaces)  {
      cout << "Resize arrays for operator multiplications" << endl;
      opMultStr.resize(nOpMult); // string that represents factor 
      opKey.resize(nOpMult); // int that represents operator:
      //0: k, 1: kx, 2: ky, 3: kz, 4: k2, 5: kx2, 6: ky2, 7: kz2
      par.resize(nOpMult);
      kIPar.resize(nOpMult);
      kJPar.resize(nOpMult);
      kIJPar.resize(nOpMult);
      kKPar.resize(nOpMult);
      for (i = 0; i < nOpMult; i++)  {
         opMultStr(i) = "~"; // not initialized
      }
   }
   // --- input for bandstructure plot
   bs = false;
   onlyBS = false;
   outputPar = 0;
   outPar = ""; // empty
   if (hamiltonian->contains("outputParameter"))  {
      outPar = hamiltonian->get("outputParameter")->toString();
      for (int pIdx = 0; pIdx < paramNames.getSize(); pIdx++)
         if (outPar == paramNames(pIdx)) {
            outputPar = pIdx;
         }
   }
   outMesh = parameters(outputPar);
   SxBinIO ioRho (outPar+".sxb",SxBinIO::BINARY_WRITE_ONLY);
   ioRho.writeMesh (parameters(outputPar), cell, mesh);
   ioRho.setMode (SxBinIO::WRITE_DATA);
   ioRho.writeMesh (parameters(outputPar), cell, mesh);
   ioRho.close();

   // --- read polarization potential.
   SxArray<SxMeshR> vPolR, vExtR, strainR, chargeR;
   SxString polFile;
   vP = zero;
   if (hamiltonian->contains("polarization"))  {
      cout << "read polarization potential..." << endl;

      dataType = -1;
      polFile = hamiltonian->get("polarization")->toString();
      if (polFile.contains(".sxb")) dataType = 0;
      else if (polFile.contains(".dat")) dataType = 1;
      else {
         cout << polFile
            << " has unknown data format. Please provide .dat or.sxb file."
            << " Exiting." << endl;
         SX_EXIT;
      };

      if (dataType == 0)  {
         SxBinIO io (polFile, SxBinIO::BINARY_READ_ONLY);
         vPolR.resize (1);
         vPolR = io.readMesh (&cell, &mesh);
         vP = (vPolR(0)*(1./HA2EV));
         io.close();

      } else if (dataType == 1)  {
         // --- read ascii polarization file
         ifstream pfile(polFile.ascii());
         if (!pfile)   {
            cout << endl << "file "
                 << polFile.ascii() << " does not exist. EXITING." << endl;
            SX_EXIT;
         }

         file = fopen(polFile.ascii(), "r+");
         vP.resize(rSize);
         for (x = 0; x < xMax; x++)  {
            for (y = 0; y < yMax; y++)  {
               for (z = 0; z < zMax; z++)  {
                  ssize_t meshIdx = mesh.getMeshIdx(x,y,z, SxMesh3D::Positive);
                  fscanf(file,"%f\n", &value);
                  vP(meshIdx) = value / HA2EV;
               }
            }
         }
      }
      // --- check if mesh size reflects provided data
      if (meshSize != vP.getSize())  {
         cout <<
            "Polarisation potential size does not match mesh size! EXITING."
            << endl;
         SX_EXIT;
      }
   }

   // --- read external potential.
   SxString extFile;
   vExt = zero;
   if (hamiltonian->contains("extPotential"))  {
      dataType = -1;
      extFile = hamiltonian->get("extPotential")->toString();
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
         vExtR.resize (1);
         vExtR = io.readMesh (&cell, &mesh);
         vExt = (vExtR(0)*(1./HA2EV));
         io.close();
      } else if (dataType == 1)  {
         // --- read ascii polarization file
         ifstream efile(extFile.ascii());
         if (!efile)   {
            cout << endl << "file "
                 << extFile.ascii() << " does not exist. EXITING." << endl;
            SX_EXIT;
         }

         file = fopen(extFile.ascii(), "r+");
         SxVector3<Int> pos;
         vExt.resize(rSize);
         for (x = 0; x < xMax; x++)  {
            pos(0) = x;
            for (y = 0; y < yMax; y++)  {
               pos(1) = y;
               for (z = 0; z < zMax; z++)  {
                  pos(2) = z;
                  int meshIdx = (int)mesh.getMeshIdx(pos, SxMesh3D::Positive);
                  fscanf(file,"%f\n", &value);
                  vExt(meshIdx) = value / HA2EV;
               }
            }
         }
      }
      if (meshSize != vExt.getSize())  {
         cout << "External potential size does not match mesh size! EXITING."
              << endl;
         SX_EXIT;
      }
   }

   if (iEpsR > -1)  {
      epsilonR = parameters(iEpsR);
   }
   else  {
      epsilonR.resize(1);
   }
   // --- read strain fields.
   if (hamiltonian->containsGroup("strain"))  {
      cout << "read strain fields..." << endl;
      strain = hamiltonian->getGroup("strain");
      SxArray<SxString> strainInput; // 3 diagonal, 3 off-diagonal strains
      dataType = -1; // 0 - binary, 1 - ascii
      strainInput.resize(6);
      eIJ.resize(6);
      strainInput(0) = strain->get("eXX")->toString();
      strainInput(1) = strain->get("eYY")->toString();
      strainInput(2) = strain->get("eZZ")->toString();
      strainInput(3) = strain->get("eXY")->toString();
      strainInput(4) = strain->get("eXZ")->toString();
      strainInput(5) = strain->get("eYZ")->toString();
      for (int idx = 0; idx < 6; idx++)  {
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
            ifstream sfile(strainInput(idx).ascii());
            if (!sfile)   {
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
      
         if (meshSize != eIJ(idx).getSize())  {
            cout << "Strain field size does not match mesh size! EXITING."
                 << endl;
            SX_EXIT;
         }
  
      } // end loop over strains 
   }
   // --- read external charges
   for (int idx = 0; idx < nCharges; idx++)  {
      if (chargeFiles(idx).contains(".sxb")) dataType = 0;
         else if (chargeFiles(idx).contains(".dat")) dataType = 1;
         else {
            cout << chargeFiles(idx)
               << " has unknown data format. Please provide .dat or.sxb file."
               << " Exiting." << endl;
            SX_EXIT;
         };
         if (dataType == 0)  {
         // --- read binary external charge file
            SxBinIO io (chargeFiles(idx), SxBinIO::BINARY_READ_ONLY);
            chargeR.resize(1);
            chargeR = io.readMesh (&cell, &mesh);
            chargeTmp = (chargeR(0));
            cout << "open " << chargeFiles(idx) << endl;
            io.close();

         } else if (dataType == 1)  {
         // --- read ascii external charge file
            ifstream cfile(chargeFiles(idx).ascii());
            if (!cfile)   {
               cout << endl << "file "
                    << chargeFiles(idx).ascii() << " does not exist. EXITING." << endl;
               SX_EXIT;
            }

            file = fopen(chargeFiles(idx).ascii(), "r+");
            SxVector3<Int> pos;
            chargeTmp.resize(rSize);
            for (x = 0; x < xMax; x++)  {
               pos(0) = x;
               for (y = 0; y < yMax; y++)  {
                  pos(1) = y;
                  for (z = 0; z < zMax; z++)  {
                     pos(2) = z;
                     int meshIdx = (int)mesh.getMeshIdx(pos, SxMesh3D::Positive);
                     fscanf(file,"%f\n", &value);
                     chargeTmp(meshIdx) = value;
                  }
               }
            } // end x,y,z loops
         }
      
         if (meshSize != chargeTmp.getSize())  {
            cout << "external charge map size does not match mesh size! EXITING."
                 << endl;
            SX_EXIT;
         }
         if (epsilonR.getSize() == 1)  {
            cout << "Error: no permittivity parameters given. "
                 << "Parameter `epsilon` required in all material files, "
                 << "alternatively, remove charges. EXITING." << endl;
            SX_EXIT;
         }
         if (idx == 0)  {totalCharge = chargeTmp * chargePrefactors(idx);}
         else {totalCharge = totalCharge + chargeTmp * chargePrefactors(idx);}
      }
   
   // --- TODO:REMOVE AFTER DEBUGGING
   cout << "found the following charge density files: " << endl;
   for (i = 0; i < nCharges; i++)  {cout << chargeFiles(i) << " * " << chargePrefactors(i) << endl;}
   cout << "total charge in system: " << totalCharge.sum() << endl;
//   SX_EXIT;
   // --- TODO:'til 'ere

   // --- show bandstructure
   wgt(0) = 1.;
   wgt(1) = 1.;
   wgt(2) = 1.;
   if (hamiltonian->contains("weight"))  {
      wgt = SxVector3<Double>(hamiltonian->get("weight")->toList());
      for (i = 0; i < 3; i ++)
         wgt(i) = 1./wgt(i);
      }
   cout << "weighting k-vectors with: " << wgt << endl;
   if (hamiltonian->containsGroup("bandstructure"))  {
      int iStep = 0;
      SxString bsOutFile;
      SxVector3<Int> rCoord;
      bandstructure = hamiltonian->getGroup("bandstructure");
      cout << "+------ Band structure calculation requested. -------" << endl;
      if (bandstructure->contains("onlyBS")) onlyBS = true;
      bs = true;
      bsOutFile = bandstructure->get("outFile")->toString();
      rCoord = SxVector3<Int>(bandstructure->get("rCoord")->toIntList());
      for (i = 0; i < 3; i++)
         if (rCoord(i) >= mesh(i))  {
            cout << "| Real space coordinate is not inside mesh." << endl;
            SX_EXIT;
         }
      cout << "| Calculate band structure at real space grid point " << rCoord
           << endl;
      cout << "| Material composition at this point: " << endl << "| ";
      double composition;
      int meshIdx = (int)mesh.getMeshIdx(rCoord, SxMesh3D::Positive);
      for (int iMat = 0; iMat < nMat; iMat++)  {
         composition = materials(iMat)(meshIdx);
         cout << matNames(iMat) << ": " << composition << "\t";
      }
      cout << endl;

      SxVector3<Double> step;
      cout << "| Set up path for band structure" << endl;
      for (path  = bandstructure->getGroup("path");
               path != NULL;
               path  = path->nextSibling ("path"))
         {
            bsStart   = SxVector3<Double>(path->get("from")->toList());
            bsEnd     = SxVector3<Double>(path->get("to")->toList());
            stepsBs   = path->get("steps")->toInt();
            step = 1. / stepsBs * (bsEnd - bsStart);
            // --- set up full path for bandstructure
            for (iStep = 0; iStep < stepsBs; iStep++)  {
               k = bsStart + iStep * step;
               bsPts << k;
            }
         }
      k = bsStart + iStep * step;
      bsPts << k;
      cout << "| Calculate band structure... " << endl;
      showBandstructure(bsOutFile, rCoord);
   }
}

void SxKdotP::printEnergies () const
{
   cout << "printEnergies()" << endl;
   cout.precision(3);
   const SxPWSet &waves = getWavesRef();
   PsiG comp, psiI, psiI2;
   PsiR compR, psiIR, compR2;
   int i, iComp, nStates;
   const SxRBasis &R = *rBasisPtr;
   nStates = waves.getNStates(); 
   cout << "composition of the wave functions of each state:" << endl;
   SxArray<SxDiracVec<Double> > rho;
   SxArray<PsiR> rComp;
   SxMatrix3<Double> cell = R.cell;
   SxMesh3D mesh = R.getMesh();
   rho.resize(nStates);

   if (moreIO)  
   for (i=0; i < nStates; i++)  {
      rho(i) = zero;
      cout << " " << i << " :\t";

      psiI = waves(i,0,0);
      psiIR = (R | psiI);
      rComp.resize(nComp);

      for (iComp = 0; iComp < nComp; iComp++)  {
           
         comp = psiI.getComponent(iComp);
         compR = psiIR.getComponent(iComp);
         rComp(iComp) = compR;
         
         cout << dot(comp, comp).re << "\t";
         // --- calculate charge density
         rho(i) = rho(i) + compR.absSqr();
      }

      SxBinIO ioRho ("rho-" + SxString(i) + ".sxb",
            SxBinIO::BINARY_WRITE_ONLY);
      ioRho.writeMesh (rho(i), cell, mesh);
      ioRho.setMode (SxBinIO::WRITE_DATA);
      ioRho.writeMesh (rho(i), cell, mesh);
      ioRho.close();

      // --- open output files
      ofstream file, wfile, file1d, outParFile, sum1d, help1d, vChgFile, vChg1d;
      fstream::openmode mode = fstream::out | fstream::trunc;
      file.open (("rho_" + SxString(i) + ".dat").ascii (), mode);
      vChgFile.open (("vChg.dat"), mode);
      vChg1d.open (("vChg1d.dat"), mode);
      wfile.open(("psi_" + SxString(i) + ".dat").ascii (), mode);
      file1d.open (("rho1d_" + SxString(i) + ".dat").ascii (), mode);
      help1d.open ("help1d.dat", mode);
      sum1d.open (("sum1d_" + SxString(i) + ".dat").ascii (), mode);
      if (outPar != "") outParFile.open ((outPar + ".dat").ascii (), mode);

      int x,y,z,idx;
      double val, val2, val0, valChg;
      val0=0;z=0;
      SxVector3<Int> pos;
      SxString waveStr;
      cout << "size of chargePotential: " << chargePotential.getSize() << endl;
      if ((chargePotential.getSize() > 1)) { cout << "write charge potential" << endl;}
      else {cout << "no charge potential" << endl;}
      for ( x = 0; x < mesh(0); x++)  {
         pos(0) = x;
         for ( y = 0; y < mesh(1); y++)  {
            pos(1) = y;
            for ( z = 0; z < mesh(2); z++)  {
               pos(2) = z;
               idx = (int)mesh.getMeshIdx(pos, SxMesh3D::Positive);
               val = rho(i)(idx);
               if (z == 0) val0 = val;
               file << SxString(val) << endl;
               // --- DEBUGGING ONLY, REMOVE LATER
               if ((chargePotential.getSize() > 1)) {
                  valChg = chargePotential(idx) * HA2EV;
                  vChgFile << SxString(valChg) << endl;
                  if ((y == mesh(0) - x) && (z == 1))  vChg1d << SxString(valChg) << endl;
               }
               // --- 'til 'ere
               waveStr = "";
               for (iComp = 0; iComp < nComp; iComp++)
                  waveStr = waveStr
                          + SxString((rComp(iComp)(idx)).re)
                          + "   "
                          + SxString((rComp(iComp)(idx)).im)
                          + "   ";
               wfile << waveStr << endl;
//               if ((x == mesh(0)/2) && (y == mesh(1)/2))  {
                  file1d << SxString(z) << "\t" << SxString(val) << endl;
                  val2 = outMesh(idx);
                  if (outPar != "") {
                     outParFile << SxString(val2) << endl;
                  }

                  if (help.getSize() > 1)
                     help1d << SxString(help(idx)) << "\t" 
                            << SxString(help2(idx)) << endl;
//               }
            }
         }
      }
      file1d << SxString(z) << "\t" << SxString(val0) << endl;
      for ( z = 0; z < mesh(2); z++)  {
         pos(2) = z;
         val = 0.;
         for ( x = 0; x < mesh(0); x++)  {
            pos(0) = x;
            for ( y = 0; y < mesh(1); y++)  {
               pos(1) = y;
               idx = (int)mesh.getMeshIdx(pos, SxMesh3D::Positive);
               val = val + rho(i)(idx);
            }
         }
         sum1d << SxString(z) << "\t" << SxString(val) << endl;
      }
      cout << endl;
      // --- close files
      file.close ();
      vChgFile.close ();
      vChg1d.close ();
      wfile.close ();
      file1d.close ();
      help1d.close ();
      sum1d.close ();
      if (outPar != "") outParFile.close ();
   }
   cout << endl;
}

SxRho &SxKdotP::getRho ()
{
   return *this;
}

void SxKdotP::computeRho (const SxPsiSet &wavesIn, const SxFermi &fermi)
{  
   SX_CHECK (dynamic_cast<const SxPW *> (&wavesIn));
   const SxPW &wavesRef = *dynamic_cast<const SxPW *> (&wavesIn);
   SxRho::computeRho (fermi.focc, wavesRef);
}

PrecEnergy SxKdotP::getEnergy (const SxPsiSet &psiSet,
                                  const SxFermi &fermi)
{
   SX_CHECK (dynamic_cast<const SxPWSet *>(&psiSet));   
   wavesPtr = static_cast<const SxPWSet *> (&psiSet);
   update (fermi);
   
   return eKin; 
}


void SxKdotP::set (const SxPsiSet &psiSet, const SxFermi &fermi)
{
   SX_CHECK (dynamic_cast<const SxPWSet *> (&psiSet));
   wavesPtr = dynamic_cast<const SxPWSet *>(&psiSet);
   compute (fermi, true, true);
   fermiPtr = &fermi;
}

void SxKdotP::readRho (const SxBinIO &io)
{
   SxRho::readRho (io);
}

void SxKdotP::normalizeRho ()
{
   SxRho::normalizeRho ();
}

void SxKdotP::writeRho (const SxString &filename) const
{
   SxRho::writeRho (filename);
}

SxDiracVec<TPrecCoeffG::TReal> 
SxKdotP::preconditioner (const PsiG &psi,
                                 Preconditioner ) const
{
   SxDiracVec<TPrecCoeffG::TReal> x, x2, x3, x4, n, K;
   x.resize(psi.getSize());
   double kin = ((D ^ psi.absSqr()).chop());
   x = D / kin;
   x2 = x.sqr();
   x3 = x.cub();
   x4 = x2.sqr();
   n  = 27. + 18.*x + 12.*x2 + 8.*x3;
   K  = n / (n + 16.*x4);

   return K; 
}
