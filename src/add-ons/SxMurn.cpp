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

#include <SxMurn.h>
#include <SxString.h>
#include <SxCLI.h>
#include <SxVector.h>
#include <SxList.h>
#include <SxBinIO.h>
#include <SxFileIO.h>

#ifndef SX_STANDALONE

SxMurn::SxMurn (const SxList<double> &volumesIn, 
                const SxList<double> &energiesIn)
{
   SX_CHECK(volumesIn.getSize () == energiesIn.getSize (),
            volumesIn.getSize (), energiesIn.getSize ());
         
   volumes = SxVector<Double>(volumesIn);
   energies = SxVector<Double>(energiesIn);
   init ();
}



SxMurn::SxMurn (const SxVector<Double> &volumesIn, 
                const SxVector<Double> &energiesIn)
 : volumes(volumesIn),
   energies(energiesIn)
{
   init ();
}

void SxMurn::init ()
{
   if (volumes.getSize () != 0)  {
      minVolume = volumes(0);
      minEnergy = energies(0);
      for (int i = 0; i < volumes.getSize (); i++)  {
         if (energies(i) < minEnergy)  {
            minEnergy = energies(i);
            minVolume = volumes(i);
         }
      }
   } else {
      minVolume = 0.;
      minEnergy = 0.;
   }
    
   bulkModulusDerivative = 4.;
}
    
double SxMurn::getEnergy (double vol)
{
   return minEnergy + bulkModulus * innerFunction (vol);
}

double SxMurn::getBulkModulus ()  {
   // 1 Hartree / Bohr^3 = 29421.01 GPa
   return bulkModulus * 29421.01;
}

double SxMurn::getMinVolume ()  {
   return minVolume;
}

double SxMurn::getMinEnergy ()  {
   return minEnergy;
}

double SxMurn::innerFunction (double vol)
{
   SX_CHECK (vol > 0., vol);
   SX_CHECK (minVolume > 0., minVolume);

   double lnRatio = log (minVolume / vol);
   // General formula not defined for B0' = 0 or 1. 
   if (bulkModulusDerivative < 1e-5)
      return (minVolume - vol * ( lnRatio + 1. ) );
   if (fabs(1. - bulkModulusDerivative) < 1e-5)
      return (vol + minVolume * lnRatio );
    
   double result = (vol - minVolume);
   result += (exp(lnRatio * bulkModulusDerivative) - 1)
             * vol / bulkModulusDerivative;
   result /= bulkModulusDerivative - 1.;
   return result;
}

double SxMurn::getDeviation ()
{
   SX_CHECK(volumes.getSize () == energies.getSize (),
            volumes.getSize (), energies.getSize ());
   int nPoints = int(volumes.getSize ());
   SX_CHECK (nPoints != 0);
   
   // vector of 1.'s
   SxVector<Double> one (nPoints); one.set (1.);

   // vector: innerFunction (vol)
   SxVector<Double> fit (nPoints);
   for (int i = 0; i < nPoints; i++)
      fit(i) = innerFunction (volumes(i));

   // calculate analytically optimized bulkModulus
   bulkModulus = (((energies.sum () / nPoints) * one - energies) * fit).sum ();
   bulkModulus /= (((fit.sum () / nPoints) * one - fit) * fit). sum ();
   bulkModulus = fabs(bulkModulus);
   // calculate analytically optimized minimum energy 
   minEnergy = (energies - bulkModulus * fit).sum () / nPoints;
      
   return (energies - minEnergy * one - bulkModulus * fit).absSqr ().sum (); 
}

void SxMurn::computeFit ()
{
   double tol = 1e-9;
   double oldBmdDeviation = getDeviation (); 
   double newDeviation;
   double oldMinVolume;

   double bmdStep = 0.1;
   double volStep = minVolume * 0.01;
   bool trySmaller = true;

   int nSteps = 0;
   
   do  {
      // change bulkModulus
      bulkModulusDerivative += bmdStep;

      // --- optimization of V0
      double oldVolDeviation = getDeviation ();
      volStep = (volStep > minVolume * 0.001) ? minVolume * 0.01 : 10*volStep;
      oldMinVolume = minVolume;
      bool trySmallerStep = true;
      do  {
         nSteps++;
         if (nSteps >= 10000)  {
            cout << "Fitting failed" << endl;
            SX_QUIT;
         }
         // change volume
         minVolume += volStep;
         while (minVolume <= 0.)  {
            volStep *= 0.5; // must be 1/2!
            minVolume-=volStep;
         }
         newDeviation = getDeviation ();
         // DEBUG print
         /*
         sxprintf("%u %f %f %f %f %g %g %g\n", nSteps, minEnergy, minVolume,
                bulkModulus, bulkModulusDerivative, volStep, bmdStep, 
                newDeviation);
         */
         if (newDeviation < oldVolDeviation)  {
            oldVolDeviation = newDeviation;
            volStep *= 1.3;
            trySmallerStep = true;
            continue;
         }
         minVolume -= volStep;
         if (oldVolDeviation == 0.) break; 
         if (fabs(oldVolDeviation - newDeviation)/oldVolDeviation < tol) break;
         if (trySmallerStep)
            volStep /= 4.;
         else
            volStep =- volStep;
         trySmallerStep = !trySmallerStep;
      } while ((fabs(volStep/minVolume) > tol));

      newDeviation = oldVolDeviation;
      if (newDeviation == 0.) break;
      if (newDeviation < oldBmdDeviation)  {
         oldBmdDeviation = newDeviation;
         //if (fabs(bmdStep) < 0.1) 
         bmdStep *=1.3;
         trySmaller = true;
         continue;
      }
      // go back
      bulkModulusDerivative -= bmdStep;
      minVolume = oldMinVolume; 
      //if (fabs(oldBmdDeviation - newDeviation)/oldBmdDeviation < tol) break;
      if (trySmaller)  {
         bmdStep /= 4.;
         trySmaller = false;
      } else {
         bmdStep = -bmdStep;
         trySmaller = true;
      }
   } while (fabs(bmdStep) > tol);
}

int getOrderOfMagnitude (double val)  {
   if (val == 0.) return 0;
   return int (floor(log10 (val)));
}
         
void SxMurn::writeToFile(SxString fileName, int nPoints, bool pressurePrint)
{
   SX_CHECK (volumes.getSize () != 0);
   SX_CHECK (nPoints > 1, nPoints);
   SX_CHECK(fileName.getSize () != 0);

   double from = minVolume, to = minVolume;

   for (int i = 0; i < volumes.getSize (); i++)  {
      if (from > volumes(i)) from = volumes(i);
      if (to < volumes(i)) to = volumes(i);
   }
   double step = (to - from);
   SX_CHECK (step > 0., step);
   from -= step * 0.05;
   to   += step * 0.05;
   step *= 1.1 / double(nPoints);

   SxVector<Double> vol (nPoints);
   SxVector<Double> energy (nPoints);
   SxVector<Double> pressure (nPoints);
   double pressPrefactor = - bulkModulus / bulkModulusDerivative
                                       // 1 Hartree / Bohr^3 = 29421.01 GPa
                                         * 29421.01;

   for (int iVol = 0; iVol < nPoints; iVol++)  {
      vol(iVol) = from + double(iVol) * step;
      energy(iVol) = getEnergy(vol(iVol));
      if (pressurePrint)  {
         pressure(iVol) = pressPrefactor *
          ( 1. - exp (-log (vol(iVol) / minVolume) * bulkModulusDerivative));
      }
   } 
   
   int xDecimals = -getOrderOfMagnitude (step) + 2;
   if (xDecimals < 0) xDecimals = 0;
   int yDecimals = -getOrderOfMagnitude ((energy.maxval () - energy.minval ())
                                         / double (nPoints) ) + 3;
   
   try  {
      SxBinIO io(fileName, SxBinIO::ASCII_WRITE_ONLY);
      fprintf(io.fp, "# Murnaghan fit by sxmurn \n");
      fprintf(io.fp, "# bulk modulus B0 = %f GPa\n", getBulkModulus ());
      fprintf(io.fp, "# bulk modulus derivative B0' = %f\n", 
              bulkModulusDerivative);
      fprintf(io.fp, "# minimum energy E0 = %f\n", getMinEnergy ());
      fprintf(io.fp, "# optimal volume V0 = %f\n", getMinVolume ());
      if (!pressurePrint)  {
         io.writeXYPlot(vol, energy, xDecimals, yDecimals);
      } else {
         for (int i = 0; i < nPoints; ++i)  {
            fprintf(io.fp, "%20.12f %20.12f %20.12f\n",
                    vol(i), energy(i), pressure(i) );
         }
      }
      io.close ();
   } catch (SxException e) {
      e.print ();
   }
}


#else /* SX_STANDALONE */


int main (int argc, char **argv)
{
   SxCLI cli (argc,argv);
   cli.preUsageMessage = 
      "Performs a fit to the Murnaghan state equation. This reads\n\n"
      "|               B0 V                    V0    | V0| B0'\n"
      "| E(V) = E0 + --------- * ( B0' * (1 - ---) + |---|    - 1 )\n"
      "|             B0'(B0'-1)                V     | V |\n\n"
      "E0, B0, B0', and V0 are the parameters being determined.";
   cli.authors = "C. Freysoldt, M. Friak";
  
   SxString inFile = cli.option
      ("-i|--input", "file",
       "input file, contains lines of <volume> <energy>.\n"
       "The volume is expected in cubic-bohr, the energy in Hartree."
      ).toString("");
   int nPoints = cli.option
      ("-n|--points","number",
       "number of points to sample"
      ).toInt(100,2);
   SxString outFile = cli.option 
      ("-o|--output", "file",
       "output file (xmgrace-readable)"
      ).toString("murn.dat");
   bool startVals = cli.option ("--start","enter start values").toBool ();

   bool printPressure = cli.option ("-p|--pressurePrint",
                                    "print first derivative (pressure)")
                        .toBool ();
      
   
   cli.version ("1.1");
   cli.finalize ();
   
   initSPHInXMath ();

   SxList<double> vList, eList;
   if (inFile.getSize () != 0)  {
      try {
         SxString data = SxFileIO::readBinary (inFile, -1);
         data.stripComments();
         SxList<SxString> lines = data.tokenize('\n');
         SxList<SxString>::Iterator lineIt;
         for (lineIt = lines.begin(); lineIt != lines.end(); ++lineIt)  {
            double V, E;
            int nRead = sscanf((*lineIt).ascii (), "%lf%lf", &V, &E);
            if (nRead == 2)  {
               if (V <= 0.)  {
                  cout << "Line '" << *lineIt 
                       << "': Volume must be larger than 0!" << endl;
                  SX_QUIT;
               }
               vList << V;
               eList << E;
            } else {
               cout << "Failed to parse line '" << (*lineIt) << "'. Skipped.\n";
            }
         }
      } catch (SxException e) {
         e.print ();
         SX_QUIT;
      }
   } else {
      cout << SxString("No input file given. If you have the input data in a "
                       "file you should use the -i option. Now, the input "
                       "data is read from the keyboard (stdin).")
              .wrap () << endl;
      cout << "Enter volume (cubic-bohr) and energy (Hartree) pairs!" << endl;
      cout << "Any non-number will end input." << endl;
      while (!feof(stdin))  {
         double vol, energy;
         int nRead = fscanf(stdin,"%lf%lf", &vol, &energy);
         if (nRead == 2)  {
            if ((vol <= 0.))  {
               cout << "Line " << vol << " " << energy;
               cout << ": Volume must be larger than 0!" << endl;
               SX_QUIT;
            }
            vList.append (vol);
            eList.append (energy);
         } else {
            break;
         }
      }
   }

   if (vList.getSize () < 2)  {
      cout << "Too little data given" << endl;
      SX_QUIT;
   }
   
   SxMurn murn (vList, eList);
   if (startVals)  {
      cout << "Enter V0: ";
      cin >> murn.minVolume;
      cout << "Enter B0': ";
      cin >> murn.bulkModulusDerivative;
      cout << "\nDeviation: " << murn.getDeviation () << endl;
   }

   murn.computeFit ();
   murn.writeToFile(outFile, nPoints, printPressure);
   if (startVals) cout << "Final deviation: " << murn.getDeviation () << endl;
   sxprintf ("# Murnaghan fit by sxmurn \n");
   sxprintf ("# bulk modulus B0 = %f GPa\n", murn.getBulkModulus ());
   sxprintf ("# bulk modulus derivative B0' = %f\n", 
             murn.bulkModulusDerivative);
   sxprintf ("# minimum energy E0 = %f\n", murn.getMinEnergy ());
   sxprintf ("# optimal volume V0 = %f\n", murn.getMinVolume ());
   
   return 0;
}


#endif /* SX_STANDALONE */
