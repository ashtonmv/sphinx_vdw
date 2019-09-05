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

#include <SxCLI.h>
#include <SxConfig.h>
#include <SxHarmonicPotential.h>
#include <SxParser.h>
#include <SxRho.h>
#include <SxStructOpt.h>
#include <SxGrid.h>
#include <SxFileParser.h>


inline double mySqrt(double x)
{
   return (x > 0.) ? sqrt(x) : -sqrt(-x);
}

int main (int argc, char **argv)
{
   SxCLI cli(argc, argv);
   SxString inFile =
      cli.option ("-i|--input","file","S/PHI/nX input file")
      .toString ("input.sx");

   SxFFT::plannerCLI (cli);

   SxString forceFile = cli.option ("--param", "file ",
                                    "force file for parametrization")
           .toString ("");
   bool param = forceFile.getSize () > 0;
   bool useSVD = cli.option ("--svd",
                                 "use singular-value decomposition for fitting")
                 .toBool ();

   SxString outFile = cli.option ("-o|--output","file","output file")
                      .toString ("");

   double incompleteShift = cli.option ("--incompleteShift",
                                        "value", "shift modes not covered by the fitting to a large value")
                            .toDouble (0.);

   cli.option ("--check", "file ", "force file. "
               "Code will print out errors in computing these forces");
   if (!param)  {
      forceFile = cli.last ().toString ("");
   } else {
      cli.last ().toString ("");
   }
   bool vaspInput = cli.option ("--vasp","parametrization from VASP data: lines of <a1> <a2> <a3> <fx> <fy> <fz>").toBool ();
   bool vaspRel = false;
   if (cli.last ().hasValue ())  {
      SxString val = cli.last ().getValue ();
      vaspRel =  (val == "relative");
   }

   SxString printHesse
      = cli.option ("--printHesse", "file", "print Hesse to file").toString ("");

   int phonons = cli.newGroup ("phonons");
   Coord kVec ( cli.option ("--kVec", "vector 3",
                            "wave vector for phonons (relative)")
                .toList3 ());

   cli.setGroup (cli.generalGroup);
   cli.option ("--betaSoft", "double", "beta for G summation");
   if (cli.last ().exists ())
      SxHarmonicPotential::betaSoft = cli.last ().toDouble ();

   SxList<int> atomList
      = cli.option ("--atoms", "index list",
                    "phonons from input file: project vibrations on these atoms"
                    ).required (false).toIdxList ();

   SxArray<int> frozenAtoms
      = cli.option ("--frozen", "index list",
                    "keep these atoms frozen"
                    ).required (false).toIdxList ();

   bool checkHesse = cli.option ("-C|--HesseCheck",
                                 "check the Gamma Hesse matrix numerically")
                     .toBool ();
   int paramQ = cli.newGroup ("parametrize Q tensors from force/field");
   Coord field (cli.option ("--field", "vector", "field in eV/bohr")
                 .toList3 ());
   field /= HA2EV;

   cli.finalize ();
   SxFFT::quickFFTPlanner ();

   // --- read input file
   SxParser parser;
   SxParser::Table table = parser.read (inFile);

   SxAtomicStructure structure (&*table);

   SxAtomicStructure refStr;
   if (table->getGroup ("harmonicPotential")->containsGroup ("structure"))  {
      refStr = SxAtomicStructure (table->getGroup ("harmonicPotential")
                                     ->getGroup ("structure"));
      refStr.replaceInfo (structure.atomInfo);
   } else {
     refStr = structure;
   }
   SxHarmonicPotential pot(&*table, refStr);

   if (table->containsGroup ("vExt"))
   {
      try {
         SxString file = table->getGroup ("vExt")->get ("file")->toString ();
         SxBinIO io(file, SxBinIO::BINARY_READ_ONLY);
         SxMesh3D mesh;
         io.read ("dim", &mesh);
         SxRBasis R(mesh, structure.cell);
         SxMeshR V = SxRho (io, &R)(0);
         V = -V; // sign convention
         pot.setExtraPot (V, false);
      } catch (SxException e) {
         e.print ();
         SX_EXIT;
      }
   }

   if (checkHesse && forceFile.getSize () == 0)  {
      SxMatrix<Complex16> H = pot.getHesse (Coord(0,0,0), structure);
      double d = 1e-4;
      double error = 0.;
      SX_LOOP(iDof)  {
         SxAtomicStructure disp(structure,SxAtomicStructure::Copy);
         disp.coordRef ()(iDof) += d;
         SxAtomicStructure df = pot.getForces (disp);
         disp.coordRef ()(iDof) -= 2. * d;
         df -= pot.getForces (disp);
         df /= -2. * d;
         cout << "iDof = " << iDof << endl;
         SxVector<Double> Hcol = H.colRef (iDof).real ();
         cout << "H analytic: " << Hcol << endl;
         cout << "H numeric:  " << df.coordRef () << endl;
         cout << "delta:      " << (Hcol - df.coordRef ()) << endl;
         double norm2 =  (Hcol - df.coordRef ()).normSqr ();
         cout << "error:      " << norm2 << endl;
         error+= norm2;
      }
      cout << SX_SEPARATOR;
      cout << "Total error: " << error << endl;
      return 0;
   }

   if (cli.groupAvailable (paramQ))  {
      ssize_t nGenQ = pot.matSymQ.getSize ();
      SxArray<SxMatrix<Double> >  xx(nGenQ), fx(nGenQ);
      int line = 0;
      SxList<int> whichAtoms;
      while (!feof(stdin))  {
         cout << "Enter atom id & force (Hartree/bohr)" << endl;
         int ia;
         SxVector<Double> f(3);
         line++;
         int nRead = fscanf (stdin, "%d %lf %lf %lf ", &ia, &f(0), &f(1), &f(2));
         if (nRead != 4)  {
            if (nRead == 0 && feof(stdin)) break;
            cout << "Failed to read " << line+1 << ". force."  << endl;
            cout << nRead << endl;
            SX_EXIT;
         }
         ia--;
         if (ia >= refStr.getSize () || ia < 0)
            cout << "Illegal atom " << (ia+1) << endl;
         int g = pot.genQ (ia);
         cout << "Atom " << ia + 1;
         if (g < 0)  {
            cout << " has no charge tensor" << endl;
            continue;
         } else {
            cout << " has charge tensor " << (g+1) << endl;
         }
         whichAtoms << ia;
         SxMatrix<Double> dfdqp = pot.paramQ (ia, field);
         if (xx(g).getSize () == 0)  {
            ssize_t np = pot.matSymQ(g).nCols ();
            xx(g).reformat (np, np); xx(g).set (0.);
            fx(g).resize (np); fx(g).set (0.);
         }
         xx(g) += dfdqp ^ dfdqp.transpose ();
         fx(g) -= dfdqp ^ f;
      }
      SX_LOOP (g)  {
         if (xx(g).getSize () == 0) continue;
         SxMatrix<Double>::Eigensystem eig = xx(g).eigensystem ();
         ssize_t np = pot.matSymQ(g).nCols ();
         SxVector<Double> res(np), space(9), ft = eig.vecs.transpose () ^ fx(g);
         res.set (0.); space.set (1.);
         SX_LOOP (i)  {
            if (fabs(eig.vals(i).re) > 1e-10)  {
               res.plus_assign_ax (ft(i)/eig.vals(i).re, eig.vecs.colRef(i));
               space += (pot.matSymQ(g) ^ eig.vecs.colRef(i)).absSqr ();
            }
            space -= pot.matSymQ(g).colRef (i).absSqr ();
         }
         res = pot.matSymQ(g) ^ res;

         cout << "// charge tensor " << (g+1) << " ("
              << pot.matSymQ(g).nCols () << " free parameters)" << endl;
         cout << "chargeTensor = [";
         for (int a = 0; a < 3; ++a)  {
            sxprintf ("[% 6.4f, % 6.4f, % 6.4f]%2s // completeness:",
                      res(3 * a + 0), res(3 * a + 1), res(3*a + 2),
                      (a == 2) ? "];" : ", ");
            for (int b = 0; b < 3; ++b)
               cout << ' ' << space(3 * a + b);
            cout << endl;
         }
         {
            SxMatrix<Double> Q(3,3);
            SX_LOOP2(a,b) Q(a,b) = res(3 * a + b);
            eig = Q.eigensystem ();
            SX_LOOP(i)
               cout << "// Q=" << eig.vals(i).re
                    << " for " << eig.vecs.colRef(i)
                    << endl;
         }
         {
            SxMatrix3<Double> Q;
            SX_LOOP2(a,b) Q(a,b) = res(3 * a + b);

            for (int i = 0; i < whichAtoms.getSize (); ++i)  {
               int ia = whichAtoms(i);
               if (pot.genQ(ia) != g) continue;
               SxMatrix3<Double> Qrot;
               if (pot.genQSymId(ia) >= 0)  {
                  SymMat S = pot.getSym(pot.genQSymId(ia));
                  Qrot = S ^ Q ^ S.transpose ();
               } else {
                  Qrot = Q;
               }
               cout << "// " << (ia+1) << " "
                    << (-field ^ Qrot) << endl;
            }
         }
      }
      return 0;
   }

   if (cli.groupAvailable (phonons))  {
      cout.precision (10);
      //if (kVec.normSqr () < 1e-20)  {
      //   cout << "No phonons at Gamma point!" << endl;
      //   SX_QUIT;
      //}
      SxCell recCell = structure.cell.getReciprocalCell ();
      kVec = recCell.relToCar (kVec);
      kVec = recCell.getMapped(kVec, SxCell::WignerSeitz);
      cout << "kVec = " << kVec << endl;
      SxMatrix<Complex16> H = pot.getHesse (kVec, structure);
      //cout << H << endl;

      /*
      {
         FILE *fp = fopen ("diagH.dat", "w");
         for (int ia = 0; ia < structure.getNAtoms (); ++ia)  {
            int is = structure.getISpecies (ia);
            fprintf (fp, "%s%d", pot.chemName(is).ascii (),
                         ia - structure.atomInfo->offset(is) + 1);
            for (int i = 0; i < 3; ++i)
               for (int j = 0; j < 3; ++j)
                  fprintf (fp, " %.8f", H(3*ia + i, 3*ia + j).re);
            fprintf (fp, "\n");
         }
         fclose (fp);
      }
      */
      if (printHesse.getSize () > 0)  {
         FILE *fp = sxfopen (printHesse, "w");
         if (H.imag ().normSqr () < 1e-12 * H.normSqr ())  {
            SX_LOOP(i)  {
               SX_LOOP(j) fprintf (fp, "%.8f\t", H(i,j).re);
               fprintf (fp, "\n");
            }
         } else {
            SX_LOOP(i)  {
               SX_LOOP(j) fprintf (fp, "(%.8f %.8f)\t", H(i,j).re, H(i,j).im);
               fprintf (fp, "\n");
            }
         }
         fclose (fp);
      }

      // --- symmetric form
      for (int i = 0; i < H.nCols (); ++i)  {
         int is = structure.getISpecies (i / 3);
         H.colRef (i) /= sqrt(pot.ionicMass (is));
      }
      H = H.transpose ();
      for (int i = 0; i < H.nCols (); ++i)  {
         int is = structure.getISpecies (i / 3);
         H.colRef (i) /= sqrt(pot.ionicMass (is));
      }
      H = H.transpose ();

      SxMatrix<Complex16>::Eigensystem eig = H.eigensystem ();
      for (int i = 0; i < eig.vals.getSize (); ++i)  {
         cout << (mySqrt(eig.vals(i).re) * AU2CM);
         if (fabs(eig.vals(i).im)  > 1e-12)  {
            cout << "+ " << eig.vals(i).im << "i";
         }
         cout << ": {(" << endl;
         for (int ja = 0; ja < structure.getNAtoms (); ++ja)  {
            SxComplex16 phase = exp( I * (kVec ^ structure(ja)));
            for (int j = 0; j < 3; ++j)  {
               SxComplex16 c = (eig.vecs(3 * ja + j, i) * phase);
               cout << eig.vecs(3 * ja +j, i) << " " << c;
               if (c.absSqr () > 1e-12) cout << " " << (c/c.abs ());
               cout << endl;
            }
         }
         cout << ")}" << endl;
      }
      return 0;
   }

   if (table->containsGroup ("phonon") && !cli.groupAvailable (phonons))  {
      SxCell recCell = structure.cell.getReciprocalCell ();
      SxKPoints kPoints;
      try {
         const SxSymbolTable *phononGroup = table->getGroup ("phonon");
         kPoints.read (recCell, phononGroup, "k", "kPoint", "kPoints");
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }

      FILE *fp = fopen ("phonon.dat", "w");
      if (!fp)  {
         cout << "Failed to open 'phonon.dat' for writing." << endl;
         SX_QUIT;
      }
      fprintf (fp, "# phonon modes in cm^-1\n");

      FILE *fp2 = NULL;
      if (atomList.getSize () > 0)  {
         fp2 = fopen ("phononWeight.dat", "w");
         fprintf (fp2, "# phonon weigts for atoms");
         SxList<int>::Iterator it;
         for (it = atomList.begin (); it != atomList.end (); it++)  {
            int ia = *it;
            if (ia >= 0 && ia < structure.getNAtoms ())  {
               fprintf (fp2, " %d", (ia+1));
            } else {
              cout << "Atom index out of range: " << (ia + 1) << endl;
              SX_QUIT;
            }
         }
         fprintf (fp2, "\n");
      }

      SxArray<bool> isFrozen;
      ssize_t nFree = 0;
      if (frozenAtoms.getSize () > 0)  {
         if (!pot.noCoulomb)  {
            cout << "Frozen atoms for Coulomb Hessian is not implemented"
                 << endl;
            SX_EXIT;
         }
         // --- set up bool array of "frozen" status
         isFrozen.resize (structure.getNAtoms ());
         isFrozen.set (false);
         SX_LOOP(iFrozen) isFrozen(frozenAtoms(iFrozen)) = true;
         SX_LOOP(ia) if(!isFrozen(ia)) nFree++;

         SX_LOOP(iTl)  {
            if (isFrozen(iTl))  {
               // remove neighbors of frozen atoms
               pot.neighborIdx(iTl).resize (0);
               pot.hesse(iTl).resize (0);
            } else {
               // set Hessian to frozen atoms to zero
               for (int iN = 0; iN < pot.neighborIdx(iTl).getSize (); ++iN)
                  if (isFrozen(pot.neighborIdx(iTl)(iN)))
                     pot.hesse(iTl)(iN).set (0.);
            }
         }
      }

      for (int ik = 0; ik < kPoints.getNk (); ++ik)  {
         kVec = kPoints.getK (ik);
         cout << "kVec = " << kVec
              << " = " << recCell.getMapped(kVec, SxCell::WignerSeitz)
              << endl;
         kVec = recCell.getMapped(kVec, SxCell::WignerSeitz);
         /*
         if (recCell.getMapped (kVec, SxCell::Origin).normSqr () < 1e-20)  {
            cout << "No phonons at Gamma point!" << endl;
            SX_QUIT;
         }
         */
         fprintf (fp, "%d", ik + 1);
         if (fp2) fprintf (fp2, "%d", ik + 1);
         SxMatrix<Complex16> H = pot.getHesse (kVec, structure);
         cout << "H hermit = " << (H - H.adjoint ()).normSqr () << endl;
         //SxMatrix<Complex16> Hel = pot.getHesseEl (kVec, structure);
         //cout << H << endl;
         if (frozenAtoms.getSize () > 0)  {
            SxMatrix<Complex16> Hfree(3 * nFree, 3 * nFree);
            SxVector<Double> mI(nFree);

            // --- get selected elements from Hessian
            const SxArray<double> &M = pot.ionicMass;
            ssize_t iTlF = 0;
            SX_LOOP2(is,ia)  {
               ssize_t iTl = structure.getIAtom(is, ia);
               if (isFrozen(iTl)) continue;
               mI(iTlF) = 1./sqrt(pot.ionicMass(is));

               ssize_t jTlF = 0;
               SX_LOOP2(js,ja)  {
                  ssize_t jTl = structure.getIAtom(js, ja);
                  if (isFrozen(jTl)) continue;

                  for (ssize_t i = 0; i < 3; i++)
                     for (ssize_t j = 0; j < 3; j++)
                        Hfree(iTlF*3 + i, jTlF*3 + j) = H(iTl*3 + i, jTl*3 + j)
                                                      / sqrt(M(is)*M(js));

                  jTlF++;
               }
               iTlF++;
            }
            H = Hfree;
         } else {
            // --- symmetric form
            for (int i = 0; i < H.nCols (); ++i)  {
               int is = structure.getISpecies (i / 3);
               H.colRef (i) /= sqrt(pot.ionicMass (is));
            }
            H = H.transpose ();
            for (int i = 0; i < H.nCols (); ++i)  {
               int is = structure.getISpecies (i / 3);
               H.colRef (i) /= sqrt(pot.ionicMass (is));
            }
            H = H.transpose ();
         }
         cout << "H hermit = " << (H - H.adjoint ()).normSqr () << endl;

         SxMatrix<Complex16>::Eigensystem eig = H.eigensystem ();
         {
            SxMatrix<Complex16> U = eig.vecs.adjoint () ^ eig.vecs;
            for (int i = 0; i < U.nCols (); ++i) U(i,i) -= 1.;
            cout << "U: " << U.normSqr () << endl;
         }
         for (int i = 0; i < eig.vals.getSize (); ++i)  {
            double omega = mySqrt(eig.vals(i).re) * AU2CM;
            fprintf (fp, "\t%9.4f", omega);

            /*
            cout << omega;
            if (fabs(eig.vals(i).im)  > 1e-12)  {
               cout << "+ " << eig.vals(i).im << "i";
            }
            cout << ": " << eig.vals(i).re * sqr(AU2CM)
                 << " " << dot(eig.vecs.colRef(i), Hel ^ eig.vecs.colRef(i)).re
                 * sqr(AU2CM)
                 << endl;
            */
            /*
            cout << ": {(" << endl;
            for (int ja = 0; ja < structure.getNAtoms (); ++ja)  {
               SxComplex16 phase = exp( I * (kVec ^ structure(ja)));
               for (int j = 0; j < 3; ++j)  {
                  SxComplex16 c = (eig.vecs.colRef(i)(3 * ja + j) * phase);
                 if (c.absSqr () > 1e-12)  {
                    //cout << eig.vecs(3 * ja +j, i) << " " << c;
                    //cout << " " << (c/c.abs ());
                    cout << c.absSqr ();
                    cout << endl;
                 } else {
                    cout << "0" << endl;
                 }
               }
            }
            cout << ")}" << endl;
            */
            if (fp2)  {
               double proj = 0.;
               SxList<int>::Iterator it;
               for (it = atomList.begin (); it != atomList.end (); it++)  {
                  int offset = 3 * *it;
                  for (int j = 0; j < 3; ++j)
                     proj += eig.vecs(offset + j, i).absSqr ();
               }
               fprintf (fp2, "\t%11.8f", proj);
            }
         }
         fprintf (fp, "\n");
         if (fp2) fprintf (fp2, "\n");
      }
      fclose (fp);
      if (fp2) fclose (fp2);
      printTiming ();
      return 0;
   }

   if (table->containsGroup ("main"))  {
      SxStructOpt structOpt (&pot, structure);
      bool hasForces = table->getGroup ("structure")
                            ->getGroup ("species")
                            ->getGroup ("atom")
                            ->contains ("force");
      if (hasForces)  {
         // --- optimize with external forces
         SxAtomicStructure resForces = structure.getNewStr ();
         // --- set external forces for pot to 0.
         resForces.set (Coord(0,0,0));
         pot.setExtraForce (resForces);
         // read external forces
         resForces = SxAtomicStructure (&*table, "force");
         // substract non-zero internal forces for this structure
         resForces -= pot.getForces (structure);
         // set this as external force for harmonic potential
         pot.setExtraForce (resForces);
      }

      SxSymbolTable *main = table->getGroup("main");
      SxSymbolTable *cmd  = NULL;
      for (cmd  = main->begin();
           cmd != NULL;
           cmd  = cmd->nextSibling())
      {
         if (SxStructOpt::isRegistered (cmd))
            structOpt.execute  (cmd);
      }
      structure.print (pot);
   } else if (table->getGroup ("harmonicPotential")
              ->containsGroup ("structure"))
   {  const SxSymbolTable *strGroup
         = table->getGroup ("harmonicPotential")->getGroup ("structure");

      bool hasForces = strGroup->getGroup ("species")
                               ->getGroup ("atom")
                               ->contains ("force");
      if (hasForces)  {
         SxAtomicStructure str(strGroup);
         str.replaceInfo (structure.atomInfo);
         // --- optimize with external forces
         SxAtomicStructure resForces = str.getNewStr ();
         // --- set external forces for pot to 0.
         resForces.set (Coord(0,0,0));
         pot.setExtraForce (resForces);
         // read external forces
         resForces = SxAtomicStructure (strGroup, "force");
         // substract non-zero internal forces for this structure
         resForces -= pot.getForces (str);
         // set this as external force for harmonic potential
         pot.setExtraForce (resForces);
      }
   }

   cout << "|d|^2 = " << (structure - refStr).absSqr ().sum () << endl;

   SxAtomicStructure forces = pot.getForces (structure, NULL);
   forces.print (pot);

   //pot.getForces (structure + Coord(1,2,3), NULL).print (pot);

   if (forceFile.getSize () > 0)  {
      SxList<SxAtomicStructure> configs, forceList, refStructs;
      SxList<int> nConfigs;

      refStructs << structure;
      nConfigs << 0;

      if (vaspInput)  {
         SxFileParser fp(forceFile);
         fp.skipWhite ();
         while (!feof(fp.fp))  {
            SxAtomicStructure config = structure.getNewStr (),
                              force  = structure.getNewStr ();
            SX_LOOP(ia)  {
               Coord x;
               SX_LOOP (iDir) fp >> x(iDir);
               if (vaspRel) x = structure.cell.relToCar (x);
               else         x *= A2B;
               // remove wrap-arounds
               x -= structure(ia);
               structure.cell.map (&x, SxCell::Origin);
               x += structure(ia);
               config.ref(ia) = x;
               SX_LOOP (iDir)
                  force.ref(ia)(iDir) = fp.getDouble ();
            }
            //config.coordRef () *= A2B;
            force.coordRef () /= HA2EV * A2B;
            fp.skipWhite ();
            nConfigs.last ()++;
            configs   << config;
            forceList << force;
         }
         cout << "Read " << configs.getSize () << " configurations" << endl;

      } else {
         // --- SPHInX format
         table = parser.read (forceFile, "std/forces.std");

         // --- loop over configurations
         for (const SxSymbolTable *item = table->begin (); item;
              item = item->nextSibling ())
         {
            // --- ignore everything which is not a structure
            if (item->getName () != "structure") continue;

            // --- read structure from file
            structure = SxAtomicStructure (item, "coords");

            // --- new configuration or new reference structure ?
            bool hasForces = false;
            try {
               hasForces = item->getGroup ("species")->getGroup ("atom")
                           ->contains ("force");
            } catch (SxException e)  {
               e.print ();
               SX_EXIT;
            }

            if (!hasForces)  {
               // --- this is a new reference structure
               refStructs << structure;
               nConfigs << 0;
               continue;
            }

            SxAtomicStructure &refStruct = refStructs.last ();
            if (refStruct.getNSpecies () != structure.getNSpecies ()
                || (refStruct.atomInfo->nAtoms - structure.atomInfo->nAtoms)
                   .absSqr ().sum () > 0
                || (refStruct.cell - structure.cell).absSqr ().sum () > 1e-10)
            {
               cout << "Structure inconsistency detected!" << endl;
               cout << "Configuration is incompatible to reference structure"
                    << endl;
               cout << "Please check content of " << forceFile << endl;
               SX_QUIT;
            }

            // --- read forces from file
            forces    = SxAtomicStructure (item, "force");
            structure.replaceInfo (refStruct.atomInfo);
            forces.replaceInfo    (refStruct.atomInfo);

            // --- add to list of configurations
            configs << structure;
            forceList << forces;
            nConfigs.last ()++;
         }
      }

      SxMatrix<Double> fxFull, xxFull;
      if (checkHesse)  {
         if (refStructs.getSize () > 1) {
            cout << "Not possible to check Hesse for multiple cells" << endl;
            SX_QUIT;
         }
         ssize_t nDof = refStructs(0).getNAtoms () * 3;
         fxFull.reformat (nDof, nDof);
         xxFull.reformat (nDof, nDof);
         fxFull.set (0.);
         xxFull.set (0.);
      }



      for (int iRef = 0; iRef < refStructs.getSize (); ++iRef)
         refStructs(iRef).updateSymGroup ();

      int nParam = pot.getNParam ();
      SxVector<Double> fx(nParam);
      SxMatrix<Double> xx(nParam, nParam);
      fx.set (0.);
      xx.set (0.);

      // --- set up parametrization problem from configurations
      SxAtomicStructure refForce = -pot.getElForces (refStructs(0));
      for (int iConfig = 0, iRef = 0;
           iConfig < configs.getSize ();
           ++iConfig, nConfigs(iRef)--)
      {
         // --- new reference structure?
         if (nConfigs(iRef) == 0)  {
            iRef++;
            table = parser.read (inFile);
            pot = SxHarmonicPotential (&*table, refStructs(iRef));
            refForce = -pot.getElForces (refStructs(iRef));
         }

         // --- get reference, configuration, forces
         SxAtomicStructure &refStruct = refStructs(iRef);
         structure = configs(iConfig);
         forces = forceList(iConfig);
         SxAtomicStructure disp = structure - refStruct;

         cout << "Configuration " << (iConfig + 1) << endl;
         for (int ia = 0; ia < disp.getNAtoms (); ++ia)  {
            if (disp.getAtom(ia).normSqr () > 1e-10)  {
               cout << "Displaced atom " << (ia+1) << " by "
                    << disp.getAtom (ia) << endl;
            }
         }

         // --- check if this is the reference configuration
         if (disp.coordRef ().normSqr () < 1e-10 * refStruct.getNAtoms ()
             && nConfigs(iRef) > 1)
         {
            cout << "Found forces for reference structure no. "
                 << (iRef + 1) << endl;
            cout << "These forces will be subtracted from "
                 << "following configurations" << endl;
            if (param)  {
               // this is the non-zero force of the reference structure
               refForce = forces - pot.getElForces (refStruct);
            } else {
               refForce.copy (forces);
               pot.setExtraForce (SxAtomicStructure ());
               pot.setExtraForce (forces - pot.getForces (structure));
            }
            // proceed to next configuration
            continue;
         }

         if (param)  {
            // subtract forces of reference structure
            forces -= refForce;
            // subtract electrostatic forces
            forces -= pot.getElForces (structure);
            cout << "harmonic forces" << endl;
            forces.print (pot);
         } else if (checkHesse)  {
            // subtract forces of reference structure
            forces -= refForce;
            SX_LOOP2(i,j)  {
               fxFull(j,i) -= forces.coords(j) * disp.coords(i);
               xxFull(j,i) += disp.coords(j) * disp.coords(i);
            }
         } else {
            SxAtomicStructure df = forces - pot.getForces (structure);
            df.print (pot);
            cout << "force norm=" << sqrt((forces-refForce).absSqr ().sum ()) << endl;
            cout << "df norm=" << sqrt(df.absSqr ().sum ()) << endl;
            continue;
         }


         // --- exploit symmetries of reference structure
         const SxSymGroup &syms = *refStruct.cell.symGroupPtr;
         SxGrid grid(refStruct, 10);
         int nSyms = syms.getSize ();
         if (nSyms > 100)  {
            cout << "Do not use translational symmetries..." << endl;
            nSyms = syms.getNPrimitive ();
         }
         nSyms = 1; cout << "No refStruct symmetries exploited" << endl;
         cout << "Applying " << nSyms << " symmetries..." << endl;
         for (int iSym = 0; iSym < nSyms; ++iSym)  {
            // --- get reordering
            SxConstPtr<SxAtomInfo> rotInfo;
            rotInfo = refStruct.match (grid, syms(iSym) ^ refStruct);

            // rotated displacement pattern
            SxAtomicStructure rotS(refStruct, SxAtomicStructure::Copy);
            rotS += (syms.getRot(iSym) ^ disp).replaceInfo (rotInfo);
            // rotated forces
            SxAtomicStructure rotForce (refStruct.getNewStr ());
            rotForce <<= (syms.getRot(iSym) ^ forces).replaceInfo (rotInfo);
            // --- update parametrization matrices
            SxMatrix<Double> dfdg = pot.paramGrad (rotS);
            fx += rotForce.coordRef () ^ dfdg;
            xx += dfdg.adjoint () ^ dfdg;
         }

      }
      if (checkHesse)  {
         SxMatrix<Double> HesseFull = xxFull.solve (fxFull.transpose ());
         SxMatrix<Double> model = pot.getHesse (Coord(0,0,0), refStructs(0));
         cout << model.normSqr () << endl;

         for (int ia = 0; ia < structure.getNAtoms (); ++ia)  {
            for (int ja = ia; ja < structure.getNAtoms (); ++ja)  {
               SxMatrix3<Double> fullHij, modelHij;
               SX_LOOP2(i,j) fullHij(i,j) = HesseFull(3 * ia + i, 3 * ja + j);
               SX_LOOP2(i,j) modelHij(i,j) = model(3 * ia + i, 3 * ja + j);
               SxMatrix3<Double> dH = fullHij - modelHij;
               if (dH.absSqr ().sum () > 1e-7)  {
                  cout << "ia=" << (ia+1) << ' ' << "ja=" << (ja+1)
                       << ": " << endl;
                  cout << "dH=" << dH << endl;
                  cout << "fullH=" << fullHij << endl;
                  cout << SxMatrix<Double>(fullHij).eigenvalues () << endl;
                  cout << "model=" << modelHij << endl;
               }
            }
         }
      }
      if (!param) return 0;

      // --- symmetry constraints
      SxMatrix<Double> freeParam = pot.getParamSpace ();
      xx = freeParam.adjoint () ^ xx ^ freeParam;
      fx = freeParam.adjoint () ^ fx;

      int nFreeParam = (int)freeParam.nCols ();
      SxVector<Complex16> res(nFreeParam);
      SxVector<Double> space(nParam);
      res.set (0.);
      space.set (1.);
      if (useSVD)  {
         res = xx.solve (fx);
      } else {
         // --- solve the parametrization problem
         SxMatrix<Complex16>::Eigensystem eig
            = SxMatrix<Complex16>(xx).eigensystem ();
         cout << eig.vals.real () << endl;

         for (int i = 0; i < nFreeParam; ++i)  {
            if (fabs(eig.vals(i).re) > 1e-10)  {
               SxVector<Complex16> u = eig.vecs.colRef(i);
               space += (freeParam ^ u).absSqr ();
               res.plus_assign_ax (dot(u,SxVector<Complex16>(fx))/eig.vals(i).re, u);
            } else {
               res.plus_assign_ax (incompleteShift, eig.vecs.colRef(i));
            }
            space -= freeParam.colRef(i).absSqr ();
         }
      }
      res = freeParam ^ res;
      int nParamBulk = nParam;
      if (pot.nGenSets.getSize () > 0)
         nParamBulk -= 9 * pot.nGenSets.sum ();
      ofstream paramFile;
      if (outFile.getSize () > 0)
         paramFile.open (outFile.ascii ());
      ostream &paramOut = (outFile.getSize () > 0) ? paramFile : cout;

      paramOut << "bulkHesse = [" << endl;
      for (int i = 0; i < nParamBulk; i += 9)  {
         if (i > 0) paramOut << "," << endl;
         paramOut << "// Hesse matrix " << (i/9 + 1) << endl;
         paramOut << "[";
         for (int b = 0; b < 3; ++b)  {
            paramOut << "["
                 << res(i + b).re << ", "
                 << res(i+3+b).re << ", "
                 << res(i+6+b).re << "]"
                 << (b < 2 ? "," : "]")
                 << " // completeness: "
                 << space(i + b) << ", "
                 << space(i+3+b) << ", "
                 << space(i+6+b)
                 << endl;
         }
      }
      paramOut << "];" << endl << endl;
      for (int iSet = 0, i=nParamBulk; iSet < pot.nGenSets.getSize (); ++iSet) {
         paramOut << "defectHesse = [" << endl;
         paramOut << "// defect neighbor type " << (iSet+1);
         paramOut << " with " << pot.nGenSets(iSet) << " independent interactions" << endl;
         for (int j = 0; j < pot.nGenSets(iSet); j++, i+=9)  {
            if (j > 0) paramOut << "," << endl;
            paramOut << "// @ " << pot.getGenPos (i/9) << endl;
            paramOut << "[";
            for (int b = 0; b < 3; ++b)  {
               paramOut << res(i + b).re << ", "
                    << res(i+3+b).re << ", "
                    << res(i+6+b).re
                    << (b < 2 ? "," : "]")
                    << " // completeness: "
                    << space(i + b) << ", "
                    << space(i+3+b) << ", "
                    << space(i+6+b)
                    << endl;
            }
         }
         paramOut << "];" << endl;
      }
      if (outFile.getSize () > 0) paramFile.close ();

      if (refStructs.getSize () == 1 && (space - 1.).normSqr () < 1e-12)  {
         pot.setHesse (res);

         int nDispl = 0;
         double Esum = 0., E2sum = 0.;
         for (int iConfig = 0; iConfig < configs.getSize (); ++iConfig) {
            const SxAtomicStructure &config = configs(iConfig);
            if ((config - refStructs(0)).coords.normSqr () > 1e-12) {
               pot.getForces (config);
               double E = pot.getEnergy ();
               Esum += E;
               E2sum += E * E;
               nDispl++;
            }
         }
         if (nDispl > 0)  {
            sxprintf ("Average energy above reference = %.12f Hartree\n",
                      Esum / nDispl);
            sxprintf ("<(E-Eavg)^2> = %.12f Hartree^2\n",
                      (E2sum - Esum * Esum / nDispl) / nDispl);
         }

         if (printHesse.getSize () > 0)  {
            SxMatrix<Double> H = pot.getHesse (Coord(0.,0.,0.), refStructs(0));
            FILE *fp = sxfopen (printHesse, "w");
            SX_LOOP(i)  {
               SX_LOOP(j) fprintf (fp, "%.8f\t", H(i,j));
               fprintf (fp, "\n");
            }
            fclose (fp);
         }
      }

      return 0;
   }

}
