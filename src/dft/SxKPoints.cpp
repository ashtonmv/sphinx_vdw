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

#include <SxKPoints.h>
#include <SxConfig.h>
#include <SxGrid.h>

SxKPoints::SxKPoints (const SxCell &cell,
                      const SxList<SxVector3<TPrecG> > &kp_,
                      const SxList<double> &weights_,
                      const SxVector3<Int> &folding_,
                      bool useInvSymmetry_)
   : useInvSymmetry(useInvSymmetry_)
{
   SX_CHECK (cell.volume > 0.);
   init (cell.getReciprocalCell (), kp_, weights_, folding_);
}

SxKPoints::SxKPoints (const SxCell        &cell,
                      const SxSymbolTable *table,
                            bool           useInvSymmetry_)
   : useInvSymmetry(useInvSymmetry_)
{
   SX_CHECK (table);
   SX_CHECK (cell.volume > 0.);
   
   // --- basis
   const SxSymbolTable *basis = NULL;
   SxString pointSpec = "k";
   try {
      if (table->containsGroup ("basis")) {
         basis = table->getGroup ("basis");
      } else if (table->getName () == "basis")  {
         basis = table;
      } else if (table->getName () == "qMesh")  {
         basis = table; pointSpec = "q";
      }  else if (table->getName() == "qPath")  {
         basis = table->parent; pointSpec = "q";
//#     ifdef WITH_TB
      } else if (table->containsGroup ("tbBasis")) {
         basis = table->getGroup ("tbBasis");
      } else if (table->getName () == "tbBasis")  {
         basis = table;
//#     endif
      } else if (table->containsGroup ("phonon")) {
         basis = table->getGroup ("phonon");
      }
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   if (!basis)  {
      cout << "Basis group is missing. " << pointSpec << "-points "
           << "cannot be initialized.\n";
      SX_QUIT;
   }
   SxString kPointId  = pointSpec + "Point";
   SxString kPointsId = (pointSpec == "k") ?  "kPoints" : "qPath";

   read (cell.getReciprocalCell (), basis, pointSpec, kPointId, kPointsId);
}

SxKPoints::SxKPoints (const SxKPoints &in)
   : nk(in.nk), foldingMP (-1,-1,-1)
{
   // silently ignore uninitialized "in"
   if (nk <= 0) return;

   // resize
   kVec.resize (nk);
   weights.resize (nk);
   bool labels = (in.kLabels.getSize () == nk);
   if (labels) kLabels.resize (nk);
   
   // Monkhorst-Pack initialization stuff is NOT copied

   // copy vectors, labels
   for (int ik = 0; ik < nk; ++ik)  {
      kVec(ik) = in.kVec(ik);
      if (labels && in.kLabels(ik).getSize () > 0)
         kLabels(ik) = in.kLabels(ik);
   }
   kVecMP = kVec;

   // copy weights
   weights <<= in.weights;
}

SxKPoints::~SxKPoints ()
{
   // empty
}

bool SxKPoints::isHighSymPoint (const SxCell &recCell, const Coord &k)
{
   double kNorm = k.norm ();
   // Gamma point: kNorm^3 is negligible compared to volume
   if (kNorm * kNorm * kNorm < 1e-18 * recCell.volume) return true;
   // check that k is inside Wigner-Seitz cell (within 1e-4)
   double kMapNorm = recCell.getMapped (k, SxCell::WignerSeitz).norm ();
   if (kNorm > 1.0001 * kMapNorm)   {
      cout << "k-Vector is outside Wigner-Seitz cell!" << endl;
      return false;
   }
   // check that k * 1.0001 is outside Wigner-Seitz cell
   kMapNorm = recCell.getMapped (k * 1.0001, SxCell::WignerSeitz).norm ();
   if (kMapNorm > kNorm) {
      cout << "Slightly outwards shiftet k-Vector is still inside Wigner-Seitz cell!" << endl;
      return false;
   }
   SxArray<SymMat> latSym = recCell.getLatticeSymmetries ();
   SymMat unity (1.,0.,0.,0.,1.,0.,0.,0.,1.), S1, S2;
   for (int iSym = 0; iSym < latSym.getSize (); ++iSym)  {
      // (S - 1)
      S1 = latSym(iSym) - unity;
      S2 = latSym(iSym) + unity;
      // unity (trivial symmetry)
      if (S1.absSqr ().sum () < 1e-6) continue;
      // inversion (trivial symmetry for WS boundary)
      if (S2.absSqr ().sum () < 1e-6) continue;
      // Symmetry maps k-vector on its own
      if (recCell.getMapped (S1 ^ k, SxCell::Origin).norm () < 1e-6 * kNorm) return true;
      // Symmetry maps k-vector to -k, which is trivial symmetry for WS boundary
      if (recCell.getMapped (S2 ^ k, SxCell::Origin).norm () < 1e-6 * kNorm) return true;
   }
   cout << "No Symmetry found that maps k-Vector on its own!" << endl;
   return false;
}

void SxKPoints::read (const SxCell &recCell,
                      const SxSymbolTable *table,
                      const SxString &pointSpec,
                      const SxString &kPointId, 
                      const SxString &kPointsId)
{
   SX_CHECK (table);
   SX_CHECK (recCell.volume > 0.);
   
   kLabels.resize (0);
   
   SxList<SxVector3<TPrecG> > kpIn;        // :ik,:{xyz}
   SxList<double>             weightsIn;   // :ik
   SxVector3<Int>             foldingIn;   // {xyz}

   // --- basis
   try {
      // --- basis.folding
      if (table->contains ("folding")) {
         foldingIn = SxVector3<Int> (table->get("folding")->toIntList());
      } else {
         foldingIn = 1; // 1x1x1
      }
      // --- basis.kPoint
      if (table->containsGroup (kPointId))  {
         SxSymbolTable *kp = NULL;
         int ik;
         SxVector3<Double> k;
         for (kp  = table->getGroup (kPointId), ik=0;
              kp != NULL;
              kp  = kp->nextSibling (kPointId), ik++)
         {
            // --- basis.kPoint.coords
            k = SxVector3<Double> (kp->get("coords")->toList());
            // --- table.relative
            bool relCoords = false;
            if (kp->contains("relative")) 
               relCoords = kp->get("relative")->toBool();

            if (relCoords)  {
               k = recCell.relToCar (k);
            }  else  {
               if (foldingIn.product() > 1)  {
                  cout << SX_SEPARATOR;
                  sxprintf ("Warning: Your %s-points are provided "
                        "in cartesian coordinates.\n", pointSpec.ascii());
                  sxprintf ("         Primitive %s-points for the Monkhorst-Pack folding are "
                        "usually\n", pointSpec.ascii());
                  sxprintf ("         given in relative coordinates, i.e. in "
                        "terms of (b1,b2,b3).\n");
                  cout << SX_SEPARATOR;
               }
            }
            kpIn << k;
            // --- basis.kPoint.weight
            if (kp->contains ("weight"))
               weightsIn << (kp->get ("weight")->toReal());
            // --- basis.kPoint.label
            if (kp->contains ("label"))
               kLabels.append (kp->get ("label")->toString());
            else
               kLabels.append ("");
               
         }
         int nWeights = int(weightsIn.getSize ());
         if (nWeights != kpIn.getSize () && nWeights != 0)  
         {
            cout << "Inconsistent " << pointSpec << "-weights!\n"
                 << "If weights are specified, they must be specified for all"
                    " points.\n";
            SX_QUIT;
         }
      }
      // --- basis.kPoints
      SxSymbolTable *kp = NULL;

      if (table->containsGroup (kPointsId))  {
         kp = table->getGroup (kPointsId);

         SxVector3<Double> k0, k1, dK, kInt;
         SxString label0, label1;
         bool relCoords = false;
         int ik, nPoints;
         
         // --- basis.kPoints.relative
         if (kp->contains ("relative"))
            relCoords = kp->get("relative")->toAttribute();
         
         // --- basis.kPoints.from
         SxSymbolTable *from = kp->getGroup ("from");
         // --- basis.kPoints.from.coords
         k0 = SxVector3<Double> (from->get ("coords")->toList());
         // --- basis.kPoints.from.label
         if (from->contains ("label"))
            label0 = from->get("label")->toString ();
         
         if (from->contains("relative") )  {
            // local relative overrides
            if (from->get("relative")->toAttribute ())
               k0 = recCell.relToCar (k0);
         } else if (relCoords)  {
            k0 = recCell.relToCar (k0);
            cout << SX_SEPARATOR;
            sxprintf ("Warning: Your %s-points are provided "
                    "in relative coordinates, i.e.\n"
                    "         in terms of (b1,b2,b3).\n", pointSpec.ascii());
            sxprintf ("%s-points for band structure calculations are "
                    "usually\n",pointSpec.ascii ());
            sxprintf ("         given in cartesian coordinates.\n");
            cout << SX_SEPARATOR;
         }

         if (!isHighSymPoint (recCell, k0))  {
            cout << "WARNING: " << k0;
            if (label0.getSize ()) cout << " (" << label0 << ")";
            cout << " is not a high symmetry point on Brillouin zone boundary"
                 << endl;
            label0 += "?!";
         }
   
         // --- add the first k-point to the list      
         kpIn      << k0;
         kLabels   << label0;
            
         // --- basis.kPoints.to
         SxSymbolTable *to = NULL;
         for (to  = kp->getGroup("to");
              to != NULL;
              to  = to->nextSibling ("to"))
         { 
            // --- basis.kPoints.to.coords
            k1 = SxVector3<Double> (to->get ("coords")->toList());
            // --- basis.kPoints.to.label
            if (to->contains ("label"))
               label1 = to->get("label")->toString ();

            // --- handle relative coordinates
            if (to->contains("relative") )  {
               // local relative overrides
               if (to->get("relative")->toAttribute ())
                  k1 = recCell.relToCar (k1);
            } else if (relCoords)  {
               k1 = recCell.relToCar (k1);
            }

            if (! isHighSymPoint (recCell, k1))  {
               cout << "WARNING: " << k1;
               if (label1.getSize ()) cout << " (" << label1 << ")";
               cout << " is not on Wigner-Seitz Brillouin zone boundary!\n";
               label1 += "?!";
            }

            // --- basis.kPoints.to.nPoints
            if (to->contains ("nPoints"))
               nPoints = to->get ("nPoints")->toInt();
            else
               nPoints = int(ceil ((k1 - k0).norm ()
                               / to->get ("dK")->toReal ()));
            
            if (nPoints > 0)  {
               // --- interpolate
               dK = (k1 - k0) / (double)nPoints;
               
               for (ik=1; ik < nPoints; ik++)  {
                  kInt = k0 + ik * dK;
                  kpIn      << kInt;
                  kLabels   << "";
              //  cout << kInt << endl;
               }
               kpIn      << k1;
               kLabels   << label1;
            } else {
               if (label1.getSize () > 0)  {
                  kLabels.last () += "=" + label1;
               }
            }

            // --- reset data for the next k-point
            k0     = k1;
            label0 = label1;
            label1 = "";
         }
         
         if (foldingIn.product () > 1)  {
            
            cout << SX_SEPARATOR;
            cout << "| Path folding for 2D and 1D projected band structures" 
                 << endl;
            cout << "| WARNING: no Monkhorst-Pack folding" << endl;
            cout << "| for normal bandstructures, set folding=[1,1,1];\n";
            cout << SX_SEPARATOR;
            int nFoldExtra = foldingIn.product () - 1;
            
            
            Coord scale;
            for (int dim = 0; dim < 3; ++dim)  {
               // pathes are folded from 0 ... (n-1)/2(n-1)
               // rather than 0 ... (n-1)/n 
               foldingIn(dim) -= 1; 
               scale(dim) = (foldingIn(dim) == 0)
                            ? 0. 
                            : 0.5 / double(foldingIn(dim));
            }
            
            // set up shift vectors
            SxVector3<Int> i;
            SxArray<Coord> shift(nFoldExtra);
            int idx = 0;
            for (i(0) = 0; i(0) <= foldingIn(0); i(0)++)
               for (i(1) = 0; i(1) <= foldingIn(1); i(1)++)
                  for (i(2) = 0; i(2) <= foldingIn(2); i(2)++)  {
                     if (i.sum () == 0) continue;
                     shift(idx++) = recCell.relToCar (scale * i);
                  }
            SX_CHECK (idx == nFoldExtra);
            
            // append shifted k-vectors to list
            int nkIn = int(kpIn.getSize ()), iShift;
            SxList<SxVector3<TPrecG> >::Iterator kIt = kpIn.begin ();
            for (ik = 0; ik < nkIn; ++ik, ++kIt)
               for (iShift = 0; iShift < nFoldExtra; ++iShift)
                  kpIn << (*kIt + shift(iShift));

            foldingIn.set(-1); // request symmetry reduction in init
         }
         
         weightsIn.resize (0);
      }
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   if (weightsIn.getSize () == 0)  {
      // equal weights
      PrecG weight = 1. / double(kpIn.getSize());
      for (int ik = 0; ik < kpIn.getSize (); ++ik)
         weightsIn << weight;
   }

   init (recCell, kpIn, weightsIn, foldingIn);
}

void SxKPoints::init (const SxCell &recCell,
                      const SxList<SxVector3<TPrecG> > &kp_,
                      const SxList<double> &weights_,
                      const SxVector3<Int> &folding_)
{
   kVecPrimitive    = kp_;
   weightsPrimitive = weights_;
   foldingMP        = folding_;
   
   int ik;

   // --- allocate memory
   nk = int(kp_.getSize ());
   kVec.resize (nk);
   weights.resize (nk);

   // --- copy kPoints array.
   for (ik=0; ik < nk; ik++)  {
      kVec(ik)    = SxVector3<TPrecG> (kp_(ik)(0), kp_(ik)(1), kp_(ik)(2));
      weights(ik) = weights_(ik);
   }

   // --- fold k-points
   if (foldingMP(0) == -1 && foldingMP(1) == -1 && foldingMP(2) == -1)  {
      cout << "Symmetry reduction ..." << endl;
      weightsMP = weights;
      kVecMP = kVec;
      symmetryReduce (recCell);
   } else if (foldingMP(0) != 1 || foldingMP(1) != 1 || foldingMP(2) != 1)  {
      generateMonkPack (recCell);
   } else  {
      kVecMP = kVec;
   }

   // --- Check norm of k-point weights
   if ( fabs(weights.sum() - 1.) > 1e-8)  {
      sxprintf ("Error: Sum of kpoint weights is incorrect.\n");
      sxprintf ("       sum = %15.12f (should be 1)\n", weights.sum());
      SX_EXIT;
   }
   SxAutoLevel ().append ("waves-k",getNk ());
}

void operator <<= (SxArray<Coord> &pos, const SxAtomicStructure &str)
{
   pos.resize (str.getSize ());
   for (int i = 0; i < str.getSize (); ++i) pos(i) = str(i);
}

void SxKPoints::generateMonkPack (const SxCell &recCell)
{
   SxVector3<Int>    i;
   SxVector3<TPrecG> k, kGridVec;
   SxArray<SxVector3<TPrecG> > kOffset(nk);
   int nFold = foldingMP.product();
   SX_CHECK (nFold > 1, nFold);

   kVecMP.resize (nFold * nk);
   weightsMP.resize (nFold * nk);

   // --- create Monkhorst-Pack mesh
   SxVector3<TReal8> folding = SxVector3<TReal8> (foldingMP);

   // Set up k vector offset
   for (int ik = 0; ik < nk; ik++)
      kOffset(ik) = recCell.relToCar (recCell.carToRel (kVec(ik)) / folding);

   // FORTRAN order to be compatible with fhi98start
   int ik, idx = 0;
   for (i(0)=0; i(0) < foldingMP(0); i(0)++)  {
      for (i(1)=0; i(1) < foldingMP(1); i(1)++)  {
         for (i(2)=0; i(2) < foldingMP(2); i(2)++)  {
            kGridVec = recCell.relToCar (i / folding);
            for (ik=0; ik < nk; ik++, idx++)  {
               kVecMP.ref(idx) = kGridVec + kOffset(ik);
               weightsMP(idx)  = weights(ik) / nFold;
            }
         }
      }
   }
   if (recCell.symGroupPtr)  {
      symmetryReduce (recCell);
   } else  { 
      kVec <<= kVecMP;
      weights = weightsMP;
      nk *= nFold;
   }
}

void SxKPoints::symmetryReduce (const SxCell &recCell)
{
   SX_CHECK (recCell.symGroupPtr);
   const SxArray<SymMat> symOp = recCell.symGroupPtr->getSymmorphic ();
   int iOp, nOp = int(symOp.getSize());
      
   // --- find generating k-points
   int ik, jk, nkMP = kVecMP.getSize ();
   SxVector3<Double> k, kS;
   
   SxArray<int>    foundIdx (nkMP);
   SxArray<bool>   found (nkMP); found.set(false);
   SxVector<TPrecWeights> newWeights (nkMP, 0.);

   if (useInvSymmetry)  {
      sxprintf ("| Time reversal symmetry used for reducing the k-point "
                "set ...\n");
   }  else  {
      sxprintf ("| Time reversal symmetry NOT used for reducing the k-point "
                "set ...\n");
   }

   int nFound=0;
   kVecMP.cell = recCell;
   kVecMP.epsEqual = 1e-6;
   SxGrid grid (kVecMP, foldingMP);
   for (ik=0; ik < kVecMP.getSize(); ik++)  {
      
      if ( found (ik) ) continue; // this k-point is done already 
      
      foundIdx(nFound)   = ik;
      newWeights(nFound) = weightsMP(ik);
      k = kVecMP(ik);         // k  in cart. coords

      for (iOp=0; iOp < nOp; iOp++)  {

         kS = symOp(iOp) ^ k;          // k' = S k   in cart. coords

         jk = kVecMP.find (kS, grid);
         if (jk >= 0 && jk != ik)  {
            if ( !found(jk) ) newWeights(nFound) += weightsMP(jk);
            found(jk) = true;
         }

         if (useInvSymmetry)  {
            jk = kVecMP.find(-kS, grid);
            if (jk >= 0 && jk != ik)  {
               if ( !found(jk) ) newWeights(nFound) += weightsMP(jk);
               found(jk) = true;
            }
         }
      }
      nFound++;
   }
   nk = nFound;

   // --- use symmetrized Monkhorst-Pack mesh from now on
   kVec.resize (nk); newWeights.resize (nk, true);
   for (ik=0; ik < nk; ik++)  kVec(ik) = kVecMP(foundIdx(ik));
   weights = newWeights;

  
/* // --- write unsymmetrized k-points
   FILE *fp = NULL;
   fp = fopen ("kp.sx", "w");
   if (fp)  {
      for (ik=0; ik < weightsMP.getSize(); ik++)  {
         fprintf (fp, "kPoint { coords = [%15.12f, %15.12f, %15.12f]; "
                               "weight = %15.12f; }\n",
                  kVecMP(ik)(0), kVecMP(ik)(1), kVecMP(ik)(2), weightsMP(ik));
      }
      fclose (fp);
   }
   

   // --- write symmetrized k-points
   fp = fopen ("kp-sym.sx", "w");
   if (fp)  {
      for (ik=0; ik < nk; ik++)  {
         fprintf (fp, "kPoint { coords = [%15.12f, %15.12f, %15.12f]; "
                               "weight = %15.12f; }\n",
                  kVec(ik)(0), kVec(ik)(1), kVec(ik)(2), weights(ik));
      }
      fclose (fp);
   }  */
}


SxVector3<TPrecG> SxKPoints::getK (int ik) const
{
   return SxVector3<TPrecG> (kVec(ik)(0), kVec(ik)(1), kVec(ik)(2));
}

void SxKPoints::read (const SxBinIO &io)
{
   SX_CHECK (io.mode == SxBinIO::BINARY_READ_ONLY);
   try {
      // read nk
      nk = io.getDimension ("nk");

      // --- read kVec
      SxMatrix<TPrecG> kVectors (nk, 3);
      io.read ("kVec", &kVectors, nk, 3);

      kVec.resize (nk);
      for (int ik = 0; ik < nk; ik++)
      {
//         kVec(ik)(0) = kVectors(ik,0);
//         kVec(ik)(1) = kVectors(ik,1);
//         kVec(ik)(2) = kVectors(ik,2);
         kVec(ik) = kVectors.row (ik).toVector3 ();
      }

      // --- read weights
      weights.resize (nk);
      io.read ("kWeights", &weights, nk);

      // clean everything we don't have (and probably don't need)
      kLabels.resize (nk);
      kVecMP.resize (0);
      weightsMP.resize (0);
      kVecPrimitive.resize (0);
      weightsPrimitive.resize (0);
      if (io.contains ("foldingMP"))  {
         io.read ("foldingMP", &foldingMP);
         cout << "foldingMP = " << foldingMP << endl;
      } else {
         foldingMP.set (-1); // no idea
      }

   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}

