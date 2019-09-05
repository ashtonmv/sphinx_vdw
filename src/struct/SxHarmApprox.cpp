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
#include <SxHarmApprox.h>
#include <SxArray.h>
#include <SxComplex.h>


SxHarmApprox::SxHarmApprox ()
{
   // potential = NULL;
   nRLat = nTau = 0;
   forceNoSym = false;
}


SxHarmApprox::SxHarmApprox (bool forceNoSym_)
{
   // potential = NULL;
   nRLat = nTau = 0;
   forceNoSym = forceNoSym_;
}


SxHarmApprox::~SxHarmApprox ()
{
   // empty
}


bool SxHarmApprox::equal(const Coord &a, const Coord &b)
{
   if ((a-b).norm() < 1e-8)
      return true;
   else 
      return false;
}


int SxHarmApprox::getBoxIndex(const Coord &a, 
      const SxArray<Coord > &box)
{
   int i;
   for (i = 0; i < box.getSize(); i++)  {
      if (equal(box(i),a))
         return i;
   } 
   return -1;
}


Coord SxHarmApprox::getRLat(const Coord &a, const SxCell &cell)
{
   SxVector3<Double> rLat;
   rLat = cell.inverse()^a;
   for (int i = 0; i < 3; i++)
      rLat(i) = floor(rLat(i) + 1e-6);
   rLat = cell^rLat;
   return rLat;
}


void SxHarmApprox::checkStr(const SxAtomicStructure &str_, 
                            const SxCell &smallCell_)
{
   for (int i = 0; i < 3; i++)  {
      // check if str is invariant against translations 
      // by smallCell basis vectors
      if (str_ + smallCell_(i) != str_)  {
         sxprintf("\nTest failed!\n"
               "Please check your structure and small cell!\n");
         SX_EXIT;
      }
      sxprintf("Moving along %d. lattice vector is fine.\n", i+1);
   }        
   cout << endl;
}


void SxHarmApprox::setIndexesAndMoved()
{
   SX_CHECK(equal(str(0,0),SxVector3<Double>(0.,0.,0.)));
   
   int js, ja;
   Coord rLat, tau;
   int index, size;
   // resize rLatIndex, tauIndex and moved so that each atom of each 
   // species of the structure can get its own index and moved value
   rLatIndex.resize(str.getNSpecies());
   tauIndex.resize(str.getNSpecies());
   moved.resize(str.getNSpecies());
   for (js = 0; js < str.getNSpecies(); js++)  {
    	rLatIndex(js).resize(str.getNAtoms(js));  
    	tauIndex(js).resize(str.getNAtoms(js));  
    	moved(js).resize(str.getNAtoms(js));  
   } 
   // set first atom indexes and moved value, also first value of boxes
   rLatIndex(0)(0) = tauIndex(0)(0) = 0;
   moved(0)(0) = 1;
   rLatBox.resize(1);
   tauBox.resize(1);
   rLatBox(0) = tauBox(0) = str.getAtom(0,0);
   
   for (js = 0; js < str.getNSpecies(); js++)  {
      for (ja = 0; ja < str.getNAtoms(js); ja++)  {
         rLat = getRLat(str(js,ja), smallCell);
         tau = str(js,ja) - rLat;         
         // assigning moved
         if (equal(rLat, str(0,0)))  {
            moved(js)(ja) = 1;
         }  else  {
            moved(js)(ja) = 0;
         }
         // assigning rLatIndex
         index = getBoxIndex(rLat, rLatBox);
         if (index == -1)  {
            size = (int)rLatBox.getSize();
            rLatBox.resize(size+1,true);
            rLatBox(size) = rLat;
            rLatIndex(js)(ja) = size;
         }  else  {
            rLatIndex(js)(ja) = index;
         }
         //assigning tauIndex
         index = getBoxIndex(tau, tauBox);
         if (index == -1)  {
            size = (int)tauBox.getSize();
            tauBox.resize(size+1,true);
            tauBox(size) = tau;
            tauIndex(js)(ja) = size;
         }  else  {
            tauIndex(js)(ja) = index;
         }        
      }
   }
   nRLat = (int)rLatBox.getSize();
   nTau  = (int)tauBox.getSize();
   
  SX_CHECK(nTau == str.getNSpecies());
}


void SxHarmApprox::centerStr()
{
  int js, ja;
  int i, j, k;
  Coord rLat, tau, newRLat, smallestRLat;
  for (js = 0; js < str.getNSpecies(); js++)  {
     for (ja = 0; ja < str.getNAtoms(js); ja++)  {
        rLat = rLatBox(rLatIndex(js)(ja));
        tau = tauBox(tauIndex(js)(ja));
        smallestRLat = rLat;
        for (i = -1; i < 2; i++)  {
           for (j = -1; j < 2; j++)  {
              for (k = -1; k < 2; k++)  {
                 newRLat = rLat + double(i)*str.cell(0) 
                    + double(j)*str.cell(1) + double(k)*str.cell(2);
                 // 1e-5 added so that equal rLat won't be changed
                 if ((newRLat.norm()-smallestRLat.norm()+1e-5) < 0)   
                    smallestRLat = newRLat;
              }
           }
        }
        str.ref(js,ja) = smallestRLat + tau;
        rLatBox(rLatIndex(js)(ja)) = smallestRLat;
     }
  }
  // saving centered structure for visualization
  saveStr(SxString("centeredStr.sx"),str);
}


void SxHarmApprox::set(/*SxPotential *potential_,*/
      const SxSymbolTable *symbolTable, const SxCell &smallCell_, 
      double deviation_)
{
   SxSpeciesData species_(symbolTable);
   SxAtomicStructure str_(symbolTable);
   // potential = potential_; TODO
   deviation = deviation_;
   species = species_;
   // check if structure is well repeated
   checkStr(str_, smallCell_);
   smallCell = smallCell_;
   str = SxAtomicStructure (str_, SxAtomicStructure::Copy);
   // for using all symmetries
   str -= str.getAtom(0,0);
   setIndexesAndMoved();
   centerStr();
   
   dynMatR.resize(nRLat);
   int r = 3*nTau;
   for (int i = 0; i < dynMatR.getSize(); i++)  
      dynMatR(i).reformat(r,r);
   
}


void SxHarmApprox::printStr()
{
   int js, ja;
   SxString mov;
   Coord a, b, c;
   for (js = 0; js < str.getNSpecies(); js++)  {
      for (ja = 0; ja < str.getNAtoms(js); ja++)  {
         if (moved(js)(ja) == 1) {mov = "true";} else {mov = "false";}
         a = str(js,ja);
         b = rLatBox(rLatIndex(js)(ja));
         c = tauBox(tauIndex(js)(ja));
         sxprintf("js = %d, ja = %d,  coord: (%f,%f,%f\nrLatIndex = %d,  "
               "tauIndex = %d,  moved = %s\nrLat from rLatBox: (%f,%f,%f)\n"
               "tau from tauBox: (%f,%f,%f)\n\n", js, ja, a(0), a(1), a(2),
               rLatIndex(js)(ja), tauIndex(js)(ja), mov.ascii(), 
               b(0), b(1), b(2), c(0), c(1), c(2));  
      }
   }
}


void SxHarmApprox::printDynMatR(bool detail)
{
   SxMatrix<Double> sum(3*nTau, 3*nTau);
   sum.set(0.);
   for (int i = 0; i < nRLat; i++)  {
      if (detail)  {
         cout << "rLat " << i << ":  " << rLatBox(i) << endl;
         cout << "sub dynamical matrix in R space:" << endl;
         cout << endl;
         dynMatR(i).print();
      }
      sum += dynMatR(i);
   }
   cout << "sum matrix: " << endl;
   sum.print();
}


void SxHarmApprox::saveDynMatR (const SxString &outFile)
{
   int r, c, nRLatIndex;
   FILE *fp;
   fp = fopen(outFile.ascii(), "w"); 
   for (nRLatIndex = 0; nRLatIndex < nRLat; nRLatIndex++)  {
      for (r = 0; r < dynMatR(nRLatIndex).nRows(); r++)  {
         for (c = 0; c < dynMatR(nRLatIndex).nCols(); c++)  {
            fprintf(fp, "%.9f\t ", dynMatR(nRLatIndex)(r,c));
         }
         if (r == 0)  {
            fprintf(fp, "\t %f  %f  %f\n", rLatBox(nRLatIndex)(0),
                 rLatBox(nRLatIndex)(1), rLatBox(nRLatIndex)(2)); 
         } else fprintf(fp,"\n");
      }
      fprintf(fp,"\n");
   }
   
   fclose(fp);
}


void SxHarmApprox::loadDynMatR (const SxString &dynMatRFile)
{  
   // some variables needed for reading the dynMatR file
   int r, c, i = 0, listSize;
   const int BUFLEN = 10240; 
   char buffer[BUFLEN];
   SxString line;
   SxList<SxString> list;
   FILE *fp = NULL;
   
   if (!(fp = fopen(dynMatRFile.ascii(), "r")))  {
      cout << "Can't open file " << dynMatRFile.ascii() << ".\n";
      SX_EXIT;
   }
   
   line = fgets(buffer, BUFLEN, fp);
      
   while (!feof(fp))  {
      list = line.tokenize (' ');
      listSize = (int)list.getSize(); 
      
      dynMatR.resize(i+1,true);
      rLatBox.resize(i+1,true);
      dynMatR(i).reformat(listSize-4,listSize-4);
      dynMatR(i).set(0.);
      // read out dynMatR
      for (c = 0; c < listSize-4; c++)  {
         dynMatR(i)(0,c) = list(c).toDouble();
      }
      // read out rLat
      for (c = listSize-3; c < listSize; c++)  {
         rLatBox(i)(c-listSize+3) = list(c).toDouble();
      }   
      // read out rest rows of dynMatR
      for (r = 1; r < listSize-4; r++)  {
         line = fgets(buffer, BUFLEN, fp);
         list = line.tokenize (' ');
         for (c = 0; c < listSize-4; c++)  {
            dynMatR(i)(r,c) = list(c).toDouble();
         }            
      }
      i++;
      // forwarding one line
      line = fgets(buffer, BUFLEN, fp);
      // for next loop
      line = fgets(buffer, BUFLEN, fp);
   }
   nRLat = (int)rLatBox.getSize();
   nTau  = int(floor((double)dynMatR(0).nCols()/3. + 1e-06));
   fclose (fp);
}


void SxHarmApprox::saveStr(const SxString &file, const SxAtomicStructure &str_)
{
   FILE *fp;
   fp = fopen(file.ascii(),"w");
   if (fp)  {
      str_.fprint(fp);
      fclose(fp);
   } else {
      cout << "Warning: could not write to '" << file << "'." << endl;
   }
}


SxAtomicStructure SxHarmApprox::getForcesDFT(const SxAtomicStructure &str_,
      const SxSymbolTable *table)
{
   SxPotential *pot = new SxHamSolver(str_,table);
   try {
      SxSymbolTable *initialGuess = table->getGroup("initialGuess");
      pot->execute(initialGuess);
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }   
   SxSymbolTable *cmd = table->getGroup("main")->begin();

   SX_EXIT;
// pot->applyVDWCorrection = false;
   return pot->getSymForces(str_, cmd);
}


void SxHarmApprox::computeDynMatR1d(int alpha, const SxSymbolTable *symbolTable)
{
   SX_CHECK(alpha == 0 || alpha == 1 || alpha == 2);
   
   int js, ja, is, ia;
   int nSpecies = str.getNSpecies();
   int beta;  // {x,y,z}
   double m1, m2;
   SxAtomicStructure forces;

   double sum = 0.;
   
   // TODO is from 0 !!!!!!!!!!!!!!!!!
   for (is = 0; is < nSpecies; is++)  {
    	for (ia = 0; ia < str.getNAtoms(is); ia++)  {
         if (moved(is)(ia) == 1)  {
            m1 = species.ionicMass(is); // atomic units
            
            // moved cause of symmorph symmetries
            str -= str.getAtom(is,ia);

            sxprintf("is: %d, ia: %d, coords: (%f,%f,%f)\n\n",
                  is, ia, str(is,ia)(0), str(is,ia)(1), str(is,ia)(2));
            
            str.ref(is,ia)(alpha) += deviation;
            str.updateSymmetries ();
            for (int g = 0; g < str.cell.symGroupPtr->getNSymmorphic(); g++)
               str.cell.symGroupPtr->getSymmorphic(g).print();
            forces = getForcesDFT(str,symbolTable); // [H/B]
            cout << forces << endl;
            for (ja = 0; ja < forces.getNAtoms(0); ja++)  {
               sum += forces(0,ja)(0);
            }
            cout << "sum: " << sum << endl;
            str.ref(is,ia)(alpha) -= deviation;
            str += str.getAtom(is,ia);
            
            for (js = 0; js < nSpecies; js++)  {
               for (ja = 0; ja < str.getNAtoms(js); ja++)  {
                  m2 = species.ionicMass(js); // atomic units
                  for (beta = 0; beta < 3; beta++)  
     //  Physics !              
     dynMatR(rLatIndex(js)(ja))  // [H/B/B/au]
        (3*tauIndex(is)(ia) + alpha, 3*tauIndex(js)(ja) + beta) 
        = forces(js,ja)(beta)/deviation/sqrt(m1*m2); 
     //             
               }
            }
            saveDynMatR(SxString("dynMatR.tmp"));
         }
      }  
   }
}


SxArray<SxArray<SymMat> > SxHarmApprox::getDynMatRSymOp()
{
   int i, j;
   SxAtomicStructure s(str, SxAtomicStructure::Copy);
   SxCell cell(s.cell);
   // enlarging cell, cause isolated str must be checked
   s.cell = 2.*cell;
   SxArray<SymMat> sOp = s.getSymmetries();
   // symOp that convert x to y or x to z
   SxArray<SxArray<SymMat> > xyzSOp(2);
   Coord test(1.,0.,0.);
   for (j = 0; j < sOp.getSize(); j++)  {
      for (i = 0; i < 2; i++)  {
         if (equal(sOp(j).row(i+1),test))  {  
            xyzSOp(i).resize(xyzSOp(i).getSize()+1,true);
            xyzSOp(i)(xyzSOp(i).getSize()-1) = sOp(j);
         }
      }
   }
   cout << "Dynamical matrix symmetries:" << endl;
   for (j = 0; j < 2; j++)  {
      sxprintf("x -> %c:\n", 121+j);
      if (xyzSOp(j).getSize() != 0)  {
         for (i = 0; i < xyzSOp(j).getSize(); i++)  {
               xyzSOp(j)(i).print();
         }
         cout << endl;
      }  else  {
         sxprintf("No Symmetries!\n %c deviation must be also calculated!\n\n",
               121+j);
      }
   }
   return xyzSOp;
}


Coord SxHarmApprox::getSubDynMatRrow(int n, int r, int c, int beta)
{
   Coord row;
   for (int i = 0; i < 3; i++)  {
      row(i) = dynMatR(n)(3*r + beta, 3*c + i);
   }
   return row;
}


void SxHarmApprox::setSubDynMatRrow(int n, int r, int c, int beta, 
      const Coord row)
{
   for (int i = 0; i < 3; i++)  {
      dynMatR(n)(3*r + beta, 3*c + i) = row(i);
   }
}


void SxHarmApprox::fillDynMatR(int beta, const SymMat &S)
{
   // TODO introduce filling y -> z, y -> x, ... 
   int i, j, r, c;
   // subDrow -> row from sub dynamical matrix in R space
   Coord symRLat, subDrow;
   for (j = 0; j < nRLat; j++)  {
      symRLat = S^rLatBox(j);
      for (i = 0; i < nRLat; i++)  {
         if (equal(symRLat,rLatBox(i)))  {
            for (r = 0; r < nTau; r++)  {
               for (c = 0; c < nTau; c++)  {
                  subDrow = getSubDynMatRrow(j,r,c,0);
                  subDrow = S^subDrow;
                  setSubDynMatRrow(i,r,c,beta,subDrow);
               }
            }
         }
      } 
   }
}


void SxHarmApprox::computeDynMatR(const SxSymbolTable *symbolTable) 
{
   int beta;
   SxArray<SxArray<SymMat> > xyzSOp = getDynMatRSymOp();
   cout << "\nDeviation in x direction\n";
   computeDynMatR1d(0, symbolTable);
   for (beta = 1; beta < 3; beta++)  {
      if (xyzSOp(beta - 1).getSize() != 0 && !forceNoSym)  {
         // taking the first one, cause it doesn't matter which one
         fillDynMatR(beta, xyzSOp(beta - 1)(0));
      }  else  {
         sxprintf("\nDeviation in %c direction\n", 120 + beta);
         computeDynMatR1d(beta, symbolTable);
      }
   }
}


// TODO sym for test purposes, also only for lattices with inv sym
void SxHarmApprox::symmetrizeDynMatR()
{
   int index1, index2; // i, beta, gamma, nDiff = 0;
   SxAtomicStructure s(str, SxAtomicStructure::Copy);
   SxCell cell(s.cell);
   Coord coord;
   //double diff = 0.;
   SxMatrix3<Double> mat1, mat2;
   // enlarging cell, cause isolated str must be checked
   s.cell = 2.*cell;
   SxArray<SymMat> sOp = s.getSymmetries();   
 //  for (i = 0; i < sOp.getSize(); i++)  {
      for (index1 = 0; index1 < nRLat; index1++)  {
         coord = -1.*rLatBox(index1);
         index2 = getBoxIndex(coord, rLatBox);
         if (index1 != index2)  {
            dynMatR(index1) = 0.5*(dynMatR(index1) + dynMatR(index2));
            dynMatR(index2) = dynMatR(index1);
         }
         
/*         coord = sOp(i)^rLatBox(index1);
         index2 = getBoxIndex(coord, rLatBox);
         if (index1 != index2)  {
            for (beta = 0; beta < 3; beta++)  {
               for (gamma = 0; gamma < 3; gamma++)  {
                  mat1(beta,gamma) = dynMatR(index1)(beta,gamma);
                  mat2(beta,gamma) = dynMatR(index2)(beta,gamma);
               }
            }
            mat2 = sOp(i)^(mat2^sOp(i));
            mat1 = 0.5*(mat1 + mat2);
            for (beta = 0; beta < 3; beta++)  {
               for (gamma = 0; gamma < 3; gamma++)  {
                  dynMatR(index1)(beta,gamma) = 
                     dynMatR(index2)(beta,gamma) = mat1(beta,gamma);
               }
            }           

            
            for (beta = 0; beta < 3; beta++)  {
               diff += (mat1 - mat2).col(beta).sum();
            }
            if (diff > 5e-04)  {
               cout << " mat1:" << mat1 << endl;
               cout << " mat2:" << mat2 << endl;
               cout << " mat-:" << mat1-mat2 << endl;
               cout << " diff:" << diff << endl << endl;
               nDiff++;
            }
            diff = 0.;
         }*/
      }
  // }

   cout << "Dynamical matrix was symmetrized !" << endl;
}


void SxHarmApprox::extendDynMatR()
{
   // TODO only valid for bravais latices, extend !!
	int nR, index, i;
   Coord rLat;
   SxArray<SymMat> sOp = str.getSymmetries(); 
   
   for (i = 0; i < sOp.getSize(); i++)  {
      for (nR = 0; nR < nRLat; nR++)  {
         rLat = sOp(i)^rLatBox(nR);
         index = getBoxIndex(rLat,rLatBox);
         if (index == -1)  {
            dynMatR.resize(dynMatR.getSize() + 1, true);
            rLatBox.resize(rLatBox.getSize() + 1, true);
            dynMatR(nR) = 0.5*dynMatR(nR);
            dynMatR(dynMatR.getSize() - 1) = dynMatR(nR);
            rLatBox(rLatBox.getSize() - 1) = rLat;         
         } 
      }
   }
   nRLat = (int)rLatBox.getSize();   
   saveDynMatR(SxString("extendedDynMatR.dat"));
   cout << "Dynamical matrix was extended !" << endl;
}
   

double SxHarmApprox::getRCut()
{
   int i, j, k, n = 1;
   Coord rLat;
   Coord outsideRLat(1e50,1e50,1e50);
   bool found = false;
   while(!found)  {
      for (i = -n; i < n+1; i++)  {
         for (j = -n; j < n+1; j++)  {
            for (k = -n; k < n+1; k++)  {
               rLat = double(i)*smallCell(0) + double(j)*smallCell(1) 
                  + double(k)*smallCell(2);
               if (getBoxIndex(rLat, rLatBox) == -1 && 
                     rLat.norm() < outsideRLat.norm())  {
                  outsideRLat = rLat;
                  found = true;
               }
            }
         }
      }
      n++;
   }
   // smallest rLat outside the structure minus some small number
   // not to include itself
   double r = outsideRLat.norm() - 1e-05;
  // cout << "rCut: " << r << endl;
  // cout << "outsideRLat: " << outsideRLat << endl;

   // saving sphere structure for visualization
   int is, ia;
   SxAtomicStructure s;
   s.cell = str.cell;
   for (is = 0; is < str.getNSpecies(); is++)  {
      s.newSpecies();
      for (ia = 0; ia < str.getNAtoms(is); ia++)  {
         rLat = rLatBox(rLatIndex(is)(ia));
         if (rLat.norm() < r)  {
            s.addAtom(str(is,ia));
         }
      }
   }
   s.endCreation();

   str = SxAtomicStructure(s, SxAtomicStructure::Copy);
   
   saveStr(SxString("sphereStr.sx"),s);
   return r;
}


void SxHarmApprox::cut()
{
   int nR, index;
   SxArray<SxMatrix<Double> > dynMatRnew;
   SxArray<Coord> rLatBoxNew;
   for (nR = 0; nR < nRLat; nR++)  {
      index = getBoxIndex(-1.*rLatBox(nR),rLatBox);
      if (index != -1)  {
         dynMatRnew.resize(dynMatRnew.getSize() + 1, true);
         rLatBoxNew.resize(rLatBoxNew.getSize() + 1, true);
         dynMatRnew(dynMatRnew.getSize() - 1) = dynMatR(nR);
         rLatBoxNew(rLatBoxNew.getSize() - 1) = rLatBox(nR);
      }  else  {
         dynMatRnew(0) += dynMatR(nR);
         cout << rLatBox(nR) << "   was taken out." << endl;
      }
   }
   dynMatR.resize(dynMatRnew.getSize());
   dynMatR = dynMatRnew;
   rLatBox.resize(rLatBoxNew.getSize());
   rLatBox = rLatBoxNew;
   nRLat = (int)rLatBox.getSize();

   cout << "Dynamical matrix was cut !" << endl;
}


void SxHarmApprox::cutSphereDynMatR()
{
   int nR;
   //rCut = getRCut();
   rCut = 13.5;
   SxArray<SxMatrix<Double> > dynMatRnew;
   SxArray<Coord> rLatBoxNew;
   for (nR = 0; nR < nRLat; nR++)  {
      if (rLatBox(nR).norm() < rCut)  {
         dynMatRnew.resize(dynMatRnew.getSize() + 1, true);
         rLatBoxNew.resize(rLatBoxNew.getSize() + 1, true);
         dynMatRnew(dynMatRnew.getSize() - 1) = dynMatR(nR);
         rLatBoxNew(rLatBoxNew.getSize() - 1) = rLatBox(nR);
      }  else  {
         dynMatRnew(0) += dynMatR(nR);
         cout << rLatBox(nR) << "   was taken out." << endl;
      }
   }
   dynMatR.resize(dynMatRnew.getSize());
   dynMatR = dynMatRnew;
   rLatBox.resize(rLatBoxNew.getSize());
   rLatBox = rLatBoxNew;
   nRLat = (int)rLatBox.getSize();

   cout << "Dynamical matrix was sphered !" << endl;
}

bool SxHarmApprox::setLoadFillSave(const SxSymbolTable *symbolTable, 
      const SxCell &smallCell_, const SxString &inputDynMatR,
      const SxString &outputDynMatR)
{
   // TODO user can specify from which row the dynMatR is filled
   set(symbolTable, smallCell_, 0);
   loadDynMatR(inputDynMatR.ascii());
   SxArray<SxArray<SymMat> > xyzSOp = getDynMatRSymOp();
   fillDynMatR(1,xyzSOp(0)(0));
   fillDynMatR(2,xyzSOp(1)(0));
   saveDynMatR(outputDynMatR);
   return true;
}


SxArray<SxVector<Complex16> > SxHarmApprox::computeEigFreq(
      const SxArray<Coord> &kPoints)
{
   SX_CHECK(nRLat > 0 && nTau > 0);
   SX_CHECK(rLatBox.getSize() == nRLat);
   SX_CHECK(dynMatR.getSize() == rLatBox.getSize());

   int i, j, nk, nR, counter = 0;
   SxMatrix<Complex16> dynMatK_(3*nTau, 3*nTau);
   SxMatrix<Complex16>::Eigensystem eig;
   SxArray<SxVector<Complex16> > eigFreq(kPoints.getSize());
   Coord k;
   for (j = 0; j < eigFreq.getSize(); j++)  
      eigFreq(j).resize(3*nTau);   
   dynMatK_.set(0.);
   
   for (nk = 0; nk < kPoints.getSize(); nk++)  {
      k = kPoints(nk);
     // cout << "k: " << k << endl;
      for (nR = 0; nR < nRLat; nR++)  {
         dynMatK_ -= dynMatR(nR)*SxComplex16
            (cos(k^rLatBox(nR)), sin(k^rLatBox(nR)));
      }
      eig = dynMatK_.eigensystem();
     // eig.print();
      for (i = 0; i < 3*nTau; i++)  {
         if (fabs(eig.vals(i).im) > 1e-04)  {
            counter++;
         }
         eigFreq(nk)(i) = sqrt(eig.vals(i).re); // [sqrt(H/B/B/au)]
         eigFreq(nk)(i) *= 96.826; // conversion to [THz]
      }
     // cout << endl;
      dynMatK_.set(0.);
   }
   if (counter > 0)  {
      sxprintf("\nImaginary part of eigenvalue was %d times (%f%%) bigger than 1e-04 !\n"
         "It will be cut, please check your dynamical matrix !\n", counter, 
         100.*double(counter)/(3.*double(kPoints.getSize())));
   }

   return eigFreq;  // [THz]
}
