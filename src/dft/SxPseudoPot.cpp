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
#include <SxPseudoPot.h>
#include <SxConstants.h>
#include <SxRadBasis.h>

SxPseudoPot::SxPseudoPot ()
{
   nlcc = false;
}


SxPseudoPot::SxPseudoPot (const SxSymbolTable *group)
   : SxSpeciesData ()
{
   int l, iSpecies;
   SxSymbolTable *species;
   SxList<int>    intList;
   SxList<double> realList;
   SxString potType;
   try  {
      if (group->getName () != "pseudoPot")
         group = group->getGroup("pseudoPot");
      SxSpeciesData::readSpecies (group);

      // --- pseudoPot
      // --- pseudoPot.species
      nSpecies = group->getGroup("species")->getNItems ("species");
      lMax.resize           (nSpecies);
      lLoc.resize           (nSpecies);
      realSpace.resize      (nSpecies);
      rGauss.resize         (nSpecies);
      pseudoFocc.resize     (nSpecies);
      foccAtomicRho.resize  (nSpecies);
      pseudoPotFiles.resize (nSpecies);
      for ( species  = group->getGroup  ("species"), iSpecies=0;
            species != NULL;
            species  = species->nextSibling ("species"), iSpecies++ )
      {
         // pseudoPot.species.potential
         pseudoPotFiles(iSpecies) = species->get("potential")->toString();
         // pseudoPot.species.lMax
         lMax(iSpecies) = species->get("lMax")->toInt();
         // pseudoPot.species.lLoc
         lLoc(iSpecies) = species->get("lLoc")->toInt();
         // pseudoPot.species.realSpaceProjectors
         realSpace(iSpecies) = species->contains("realSpaceProjectors")
                             ? species->get("realSpaceProjectors")
                               ->toAttribute ()
                             : false;
         // pseudoPot.species.rGauss
         rGauss(iSpecies) = species->get("rGauss")->toReal();
         // pseudoPot.species.lcaoOrbitals
         intList = species->get("lcaoOrbitals")->toIntList();
         pseudoFocc(iSpecies).resize (lMax(iSpecies)+1);
         for (l=0; l <= lMax(iSpecies); l++)
            pseudoFocc(iSpecies)(l) = intList.contains(l);
         if (species->contains("potType")) potType = species->get("potType")->toString();
         else potType = "FHI";
         // pseudoPot.species.atomicRhoOcc
         if (species->contains("atomicRhoOcc"))  {
            SxSymbol *sym = species->get("atomicRhoOcc");
            int ind, iSpin, nSpin = sym->getRank();
            if (nSpin == 0)  nSpin = 1; // for hydrogen/helium
            foccAtomicRho(iSpecies).resize(lMax(iSpecies)+1);
            for (l=0; l <= lMax(iSpecies); l++)
               foccAtomicRho(iSpecies)(l).resize (nSpin);
            realList  = species->get("atomicRhoOcc")->toList();
            if (realList.getSize() != nSpin*(lMax(iSpecies)+1))  {
               cout << "In species ";
               if (elementName(iSpecies).getSize () > 0)
                  cout << elementName(iSpecies);
               else if (chemName(iSpecies).getSize () > 0)
                  cout << chemName(iSpecies);
               else cout << (iSpecies+1);
               sxprintf (": %d <--> %d\n", (int)realList.getSize(),
                         nSpin*(lMax(iSpecies)+1));
               cout << "Corrupt atomicRhoOcc!" << endl;
               SX_QUIT;
            }
            for (iSpin = 0; iSpin < nSpin; iSpin++)  {
               for (l=0; l <= lMax(iSpecies); l++)  {
                  ind = l + iSpin*(lMax(iSpecies)+1);
                  foccAtomicRho(iSpecies)(l)(iSpin) = realList(ind);
               }
            }
         }
      }

   }  catch (SxException e) {
      e.print();
      SX_EXIT;
   }

   SxArray<SxString> psiFiles(nSpecies); // missing!
   if (potType == "FHI") initPseudoPotentials ();
   else if (potType == "Siesta") initPseudoPotSiesta (psiFiles);
   else {
      cout << "Unknown potential type, S/PHI/nX quits here!" << endl;
      SX_QUIT;
   }
}   

SxPseudoPot::~SxPseudoPot ()
{
   // empty
}


/* In old SPHInX versions, the pseudopotentials were read as floats.
   This changes energies by 1e-5 to 1e-4 Hartree/atom.
   To allow for a direct comparison with old calculations, define the
   SX_PP_READFLOAT preparser macro.
*/
// #define SX_PP_READFLOAT
void SxPseudoPot::initPseudoPotentials ()
{
   // TODO: separate pseudopotential reader...
   const int BUFLEN = 1024;
   char buffer[BUFLEN];

   FILE *fp;
   char const *filename = NULL;
#ifdef SX_PP_READFLOAT
#warning "float truncation enabled - this is deprecated!!"
   cout << SX_SEPARATOR;
   cout << "| WARNING: float-truncation in PP activated!" << endl;
   cout << SX_SEPARATOR;
   float zv, radf, psi, pot, rhoCore, lDr;
#else
   double zv, radf, psi, pot, rhoCore, lDr;
#endif
   int i, j, n, l, lMaxIdx, lLocIdx, iSpecies;
   int nVarsRead, iPoint;
   Real8 norm;

   SxList<Real8> radList, psiList, potList;  // :r
   nlcc = false;

   rad.resize (nSpecies);
   pseudoPsi.resize (nSpecies);
   pseudoPot.resize (nSpecies);
   radCore.resize   (nSpecies);
   logDr.resize   (nSpecies);
   for (iSpecies=0; iSpecies < nSpecies; iSpecies++)  {
      radCore(iSpecies).resize(2);
      filename = pseudoPotFiles (iSpecies).ascii();
      if ( !(fp = fopen (filename, "r")) )  {
         sxprintf ("Can't open file %s\n", filename);
         SX_EXIT;
      }

      // --- read valence charge
#ifdef SX_PP_READFLOAT
      nVarsRead = fscanf (fp, "%e %d", &zv, &lMaxIdx);
#else
      nVarsRead = fscanf (fp, "%le %d", &zv, &lMaxIdx);
#endif
      if (nVarsRead != 2)  {
         cout << "Failed to read valence charge and max l.\n'"
              << filename << "' is not a valid pseudotential file.\n";
         SX_QUIT;
      }
      if (lMaxIdx > lMax(iSpecies)+1)  {
         cout << "WARNING: omitting l > " << lMax(iSpecies) 
              << " from pseudopotential file '" << filename
              << "'." << endl;
      } else if (lMaxIdx <= lMax(iSpecies))  {
         cout << "Requested l=0.." << lMax(iSpecies) << " from '"
              << filename << "', but it contains only " << lMaxIdx
              << " l-channels.\nCheck the pseudopotential metadata!\n";
         SX_QUIT;
      }
      fgets (buffer, BUFLEN, fp);  // skip rest of line

      // --- skip following 10 lines
      for (i=0; i<10; i++)  fgets (buffer, BUFLEN, fp);

      pseudoPsi(iSpecies).resize (lMax(iSpecies) + 1);
      pseudoPot(iSpecies).resize (lMax(iSpecies) + 1);
      // --- read orbitals (2nd column=rad, 3th column=psi,
      //                    4th column=potential)
      for (l=0; l < lMaxIdx; l++)  {
         // --- number of entries of next orbital
#ifdef SX_PP_READFLOAT
         nVarsRead = fscanf (fp, "%d %f", &n, &lDr);
#else
         nVarsRead = fscanf (fp, "%d %lf", &n, &lDr);
#endif
         if (nVarsRead != 2 || n <= 0)  {
            cout << "Failed to read l=" << l << ".\nPseudotential file '"
                 << filename << "' is corrupted.\n";
            SX_QUIT;
         }
         lDr = log(lDr);
         if (l <= lMax(iSpecies))
            logDr(iSpecies) = lDr; // radial mesh is checked below
         fgets (buffer, BUFLEN, fp);  // skip rest of line
         for (i=0; i<n; i++)  {
#ifdef SX_PP_READFLOAT
            nVarsRead = fscanf (fp, "%d %e %e %e", &iPoint, &radf, &psi, &pot);
#else
            nVarsRead = fscanf (fp, "%d %le %le %le", &iPoint, &radf, &psi, &pot);
#endif
            if (nVarsRead != 4 || iPoint != (i+1))  {
               cout << "Failed to read radial point " << (i+1)
                    << " for l=" << l << ".\nPseudotential file '"
                    << filename << "' is corrupted.\n";
               SX_QUIT;
            }
            radList.append (radf);
            psiList.append (psi/radf);
            potList.append (pot);
         }
         norm = sqrt(( SxVector<TReal8> (radList).cub()
                     * SxVector<TReal8> (psiList).sqr()).integrate(lDr)
                    );
////         for (j=0; j<n; j++)  psiList(j) *= exp (-rad/rc)^2;
         if (fabs (norm) > 1e-10)
            for (j=0; j<n; j++)  psiList(j) /= norm;

         if (l <= lMax(iSpecies))  {
            if (rad(iSpecies).getSize () == 0)
               rad(iSpecies) = radList;
            else if (rad(iSpecies).getSize () != radList.getSize ()
                     || (rad(iSpecies) - SxDiracVec<Double>(radList))
                        .normSqr () > 1e-12 * (double)radList.getSize ())  {
               cout << "In pseudopotential file " << filename << ":" << endl;
               cout << "Varying radial meshes are not supported within a "
                       "species." << endl;
               cout << "If necessary, interpolate the offending channel l="
                    << l << endl;
               SX_QUIT;
            }

            pseudoPsi(iSpecies)(l) = psiList;
            // Set up auxData
            SxDiracAux &auxData = pseudoPsi(iSpecies)(l).handle->auxData;
            auxData.is = iSpecies;
            auxData.n = 0;
            auxData.l = l;
            pseudoPot(iSpecies)(l) = potList;
         }

         radList.removeAll ();
         psiList.removeAll ();
         potList.removeAll ();
      }

      // --- NLCC
      clearerr (fp);
      int ir = 0;
      while ( !feof (fp) )  {
         // TODO: old/new pseudopotentials (2/4 columns)
#ifdef SX_PP_READFLOAT
         nVarsRead = fscanf (fp, "%e %e", &radf, &rhoCore);
#else
         nVarsRead = fscanf (fp, "%le %le", &radf, &rhoCore);
#endif
         if (nVarsRead != 2 && nVarsRead != EOF)  {
            cout << "Failed to read point " << (radCore(iSpecies)(0).getSize ()+1)
                 << " for pseudo core density.\nPseudotential file '"
                 << filename << "' is corrupted.\n";
            SX_QUIT;
         }
         fgets (buffer, BUFLEN, fp);  // skip rest of line
         if ( !feof (fp) )  {
            radCore(iSpecies)(0).append (radf);
            radCore(iSpecies)(1).append (rhoCore);
            if (fabs(radf-rad(iSpecies)(ir)) > 1e-12 * fabs(radf))  {
               cout << "In PP-file " << filename << ":" << endl;
               cout << "Pseudo-core density must have same mesh as psi's.\n";
               SX_QUIT;
            }
            ir++;
         }
      }
      if (radCore(iSpecies)(0).getSize() > 0)
         nlcc = true;

      fclose (fp);
   }

   // --- substract local part from pseudopotential
   for (iSpecies=0; iSpecies < nSpecies; iSpecies++)  {
      lLocIdx = lLoc(iSpecies);
      for (l=0; l<=lMax(iSpecies); l++)  {
         if (l != lLocIdx )  {
            pseudoPot(iSpecies)(l) -=  pseudoPot(iSpecies)(lLocIdx);
         }
      } // l
   } // iSpecies
}

void SxPseudoPot::initPseudoPotSiesta (SxArray<SxString> /*psiFiles*/)
{
   // TODO: separate pseudopotential reader...
   const int BUFLEN = 1024;
   char buffer[BUFLEN];

   FILE *fp;
   char const *filename = NULL;
#ifdef SX_PP_READFLOAT
#warning "float truncation enabled - this is deprecated!!"
   cout << SX_SEPARATOR;
   cout << "| WARNING: float-truncation in PP activated!" << endl;
   cout << SX_SEPARATOR;
#endif

   rad.resize (nSpecies);
   pseudoPsi.resize (nSpecies);
   pseudoPot.resize (nSpecies);
   radCore.resize   (nSpecies);
   logDr.resize   (nSpecies);
   for (int iSpecies=0; iSpecies < nSpecies; iSpecies++)  {
      radCore(iSpecies).resize(2);
      filename = pseudoPotFiles (iSpecies).ascii();
      if ( !(fp = fopen (filename, "r")) )  {
         sxprintf ("Can't open file %s\n", filename);
         SX_EXIT;
      }

      // --- read header
      SxString element, XC, rel, pseudoCore;
      int nVarsRead = fscanf (fp, "%s %s %s %s", 
            element.elements, 
            XC.elements, 
            rel.elements, 
            pseudoCore.elements);
      if (nVarsRead != 4)  {
         cout << "Failed to read chemname, XC, relativistic treatment, pseudo core.\n'"
              << filename << "' is not a valid pseudotential file.\n";
         SX_QUIT;
      }

      if (!(element = elementName(iSpecies))) {
         cout << "Wrong element! Is " << element << "should be " << elementName(iSpecies) << endl;
         SX_QUIT;
      }
      fgets (buffer, BUFLEN, fp);  // skip rest of line
      fgets (buffer, BUFLEN, fp);  // skip rest of line
      int nDown, nUp, nGrid;
      double r0, zVal;
#ifdef SX_PP_READFLOAT
      nVarsRead = fscanf (fp, "%i %i %i %e %e %f", &nDown, &nUp, &nGrid, &r0, &logDr(iSpecies), &zVal);
#else
      nVarsRead = fscanf (fp, "%i %i %i %le %le %lf", &nDown, &nUp, &nGrid, &r0, &logDr(iSpecies), &zVal);
#endif
      if (nUp > 0) {
         cout << "Spin in pseudoPotential not yet implemented!" << endl;
         SX_EXIT
      }
      rad(iSpecies).resize(nGrid);
      pseudoPot(iSpecies).resize(nDown);
      pseudoPsi(iSpecies).resize(nDown);
      int lines = nGrid / 4;
      int lastEntrys = nGrid % 4;

      fgets (buffer, BUFLEN, fp);  // skip rest of line
      // read radial Basis
      int counter = 0;
      SxDiracVec<TReal8> &r = rad(iSpecies);
      for (int iLine = 0; iLine < lines; iLine++,counter+=4)  {
#ifdef SX_PP_READFLOAT
            nVarsRead = fscanf (fp, "%e %e %e %e", &r(counter), &r(counter+1), &r(counter+2), &r(counter+3));
#else
            nVarsRead = fscanf (fp, "%le %le %le %le", &r(counter), &r(counter+1), &r(counter+2), &r(counter+3));
#endif
      }
      switch (lastEntrys)  {
         case 1:
#ifdef SX_PP_READFLOAT
            nVarsRead = fscanf (fp, "%e", &r(counter));
#else
            nVarsRead = fscanf (fp, "%le", &r(counter));
#endif
            break;
         case 2:
            #ifdef SX_PP_READFLOAT
            nVarsRead = fscanf (fp, "%e %e", &r(counter), &r(counter+1));
#else
            nVarsRead = fscanf (fp, "%le %le", &r(counter), &r(counter+1));
#endif
            break;
         case 3:
            #ifdef SX_PP_READFLOAT
            nVarsRead = fscanf (fp, "%e %e %e", &r(counter), &r(counter+1), &r(counter+2));
#else
            nVarsRead = fscanf (fp, "%le %le %le", &r(counter), &r(counter+1), &r(counter+2));
#endif
            break;
         default:
            break;
      }

      // read PseudoPotential
      for (int iDown = 0; iDown < nDown; iDown++)  {
         counter = 0;
         int l;
         fgets (buffer, BUFLEN, fp);  // skip rest of line
         nVarsRead = fscanf (fp, "%i", &l);
         pseudoPot(iSpecies)(l).resize(nGrid);
         SxDiracVec<TReal8> &pp = pseudoPot(iSpecies)(l);
         for (int iLine = 0; iLine < lines; iLine++,counter+=4)  {
#ifdef SX_PP_READFLOAT
            nVarsRead = fscanf (fp, "%e %e %e %e", &pp(counter), &pp(counter+1), &pp(counter+2), &pp(counter+3));
#else
            nVarsRead = fscanf (fp, "%le %le %le %le", &pp(counter), &pp(counter+1), &pp(counter+2), &pp(counter+3));
#endif
         }
         switch (lastEntrys)  {
            case 1:
#ifdef SX_PP_READFLOAT
               nVarsRead = fscanf (fp, "%e", &pp(counter));
#else
               nVarsRead = fscanf (fp, "%le", &pp(counter));
#endif
               break;
            case 2:
#ifdef SX_PP_READFLOAT
               nVarsRead = fscanf (fp, "%e %e", &pp(counter), &pp(counter+1));
#else
               nVarsRead = fscanf (fp, "%le %le", &pp(counter), &pp(counter+1));
#endif
               break;
            case 3:
#ifdef SX_PP_READFLOAT
               nVarsRead = fscanf (fp, "%e %e %e", &pp(counter), &pp(counter+1), &pp(counter+2));
#else
               nVarsRead = fscanf (fp, "%le %le %le", &pp(counter), &pp(counter+1), &pp(counter+2));
#endif
               break;
            default:
               break;
         }
      }

      // read core Charge
      /*
      radCore(iSpecies)(0) = 1.0 * rad(iSpecies);
      radCore(iSpecies)(1).resize(nGrid);
      counter = 0;
      nlcc = true;
      SxDiracVec<TReal8> &rc = radCore(iSpecies);
      for (int iLine = 0; iLine < lines; iLine++,counter+=4)  {
#ifdef SX_PP_READFLOAT
            nVarsRead = fscanf (fp, "%e %e %e %e", &rc(1)(counter), &rc(1)(counter+1), &rc(1)(counter+2), &rc(1)(counter+3));
#else
            nVarsRead = fscanf (fp, "%le %le %le %le", &rc(1)(counter), &rc(1)(counter+1), &rc(1)(counter+2), &rc(1)(counter+3));
#endif
      }
      switch (lastEntrys)  {
         case 1:
#ifdef SX_PP_READFLOAT
            nVarsRead = fscanf (fp, "%e", &rc(counter));
#else
            nVarsRead = fscanf (fp, "%le", &rc(counter));
#endif
            break;
         case 2:
            #ifdef SX_PP_READFLOAT
            nVarsRead = fscanf (fp, "%e %e", &rc(counter), &rc(counter+1));
#else
            nVarsRead = fscanf (fp, "%le %le", &rc(counter), &rc(counter+1));
#endif
            break;
         case 3:
            #ifdef SX_PP_READFLOAT
            nVarsRead = fscanf (fp, "%e %e %e", &rc(counter), &rc(counter+1), &rc(counter+2));
#else
            nVarsRead = fscanf (fp, "%le %le %le", &rc(counter), &rc(counter+1), &rc(counter+2));
#endif
            break;
         default:
            break;
      }
      */

      fclose(fp);

      // read pseudoPsi missing

   }

   // --- substract local part from pseudopotential
   for (int iSpecies=0; iSpecies < nSpecies; iSpecies++)  {
      int lLocIdx = lLoc(iSpecies);
      for (int l = 0; l <= lMax(iSpecies); l++)  {
         if (l != lLocIdx )  {
            pseudoPot(iSpecies)(l) -=  pseudoPot(iSpecies)(lLocIdx);
         }
      } // l
   } // iSpecies
}

void SxPseudoPot::print () const
{
   int iSpecies;
   Real8 lDr;
   for (iSpecies=0; iSpecies < nSpecies; iSpecies++)  {
      sxprintf ("|    %-15s  ", elementName(iSpecies).ascii());
      cout << "rGauss=" << rGauss(iSpecies)
           << ", zv="   << valenceCharge(iSpecies)
           << ", loc="  << lLoc(iSpecies);
      if (radCore(iSpecies)(1).getSize() > 0)  {
         lDr = logDr(iSpecies);
         cout << ", nCoreEl=" <<
            (
             (SxVector<TReal8> (radCore(iSpecies)(1))
            * SxVector<TReal8> (radCore(iSpecies)(0)).cub()).integrate(lDr)
           );
      }
      if (realSpace(iSpecies)) cout << ", rsp";
      cout << endl;
      cout << "|             Potential:   " << "file:/"
                                            << pseudoPotFiles(iSpecies)
                                            << endl;
   }
}

void SxPseudoPot::reConfinePsi (double sigma) 
{
   int iSpecies, il;
   double norm;

   SxDiracVec<TReal8> psi;
   SxDiracVec<TReal8> radFunc;

   for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)  {
      for (il = 0; il < pseudoPsi(iSpecies).getSize(); il++)  {
         
         psi = pseudoPsi(iSpecies)(il); // reference
         radFunc = rad(iSpecies);

         // mutiply psi by Gaussian : Exp(-r^2/(2*sigma^2))
         psi *= exp ( (- radFunc.sqr()/(2.0*sigma*sigma)) );

         // --- re-normalize
         norm = sqrt((radFunc.cub() * psi.sqr()).integrate(logDr(iSpecies)));
         if (fabs (norm) > 1e-10) psi /= norm;

      }
   }
   
}

int SxPseudoPot::getLMax () const
{
   int lmax = -1, ll;
   for (int iSpecies = 0; iSpecies < lMax.getSize (); ++iSpecies)
      if (lmax < (ll = lMax(iSpecies))) lmax = ll;
   SX_CHECK (lmax >= 0);
   return lmax;
}

SxArray<SxArray<SxDiracVec<TReal8> > > SxPseudoPot::getPseudoPsi ()
{
   SxArray<SxArray<SxDiracVec<Double> > > result (nSpecies);
   for (int is = 0; is < nSpecies; is++)   {
      result(is).resize(lMax(is)+1);
      for (int l = 0; l <= lMax(is); l++)   {
         result(is)(l) = pseudoPsi(is)(l).getCopy ();
      }
   }
   return result;
}

SxDiracVec<TReal8> SxPseudoPot::getPseudoPsi (int is, int l)
{
   return pseudoPsi(is)(l).getCopy ();
}

SxDiracVec<TPrecCoeffG>
SxPseudoPot::getAtomRhoG(const SxGBasis &G, int iSpecies, int iSpin) const
{
   SX_CHECK(G.structPtr);
   double volume = G.structPtr->cell.volume;
   SxRadBasis r (rad, logDr);  // |r>
   const SxDiracVec<TReal8> &psRad = rad(iSpecies);
   SxDiracVec<TReal8> radFun (psRad.getSize ());
   radFun.set(0.);
   radFun.setBasis (&r);
   for (int l=0; l <= lMax(iSpecies); l++)  {
      double focc = 0.;
      if (iSpin >= 0) {
         focc =  foccAtomicRho(iSpecies)(l)(iSpin);
      } else {
         SX_LOOP(jSpin) focc += foccAtomicRho(iSpecies)(l)(jSpin);
      }
      radFun += focc * pseudoPsi(iSpecies)(l).sqr();
   }
   //cout << "AO density should be projected Dirac-like" << endl;
   PsiG psiAtom  = r.toPWBasis  (psRad,radFun, G, 0, logDr(iSpecies));
   psiAtom *= SxYlm::getYlmNormFactor (0,0) * SQRT_1_4PI // = Y(0,0)
           * FOUR_PI / sqrt(volume); // normalization
   return psiAtom;
}
