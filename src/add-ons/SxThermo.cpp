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

#include <SxThermo.h>
#include <SxParser.h>
#include <iostream>
#include <fstream>
#include <SxCLI.h>
#include <SxBinIO.h>

#ifndef SX_STANDALONE

SxVector<Double> SxThermo::logVec (const SxVector<Double> &y) const
{
   int nElt = int(y.getSize ()),
       nIg = 0;//no. ignored freq.
   SxVector<Double> res (nElt, 0.); //initialize with 0.
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
   for (ssize_t i = 0; i < y.getSize (); ++i)
      if (y(i) > 1e-4)
         res(i) = log(y(i));
      else
         nIg++;

   // rescale weighing
   if (nIg < nElt)  { res /= (1. - nIg * 1. / nElt); }
   else {
      cout << "All freq. too small!" << endl;
      SX_EXIT;
   }
   if (nIg * 10000. > nElt)
   { cout << " Fract.ign.freq: " << (nIg * 1. / nElt); }
   return res;
}

SxMatrix<Double> SxThermo::readFreq (const SxString &input,
                                           double &vol) const
{  cout << "readFreq " << input << endl;
   SX_CHECK (input != "");
   SxMatrix<Double> freq;//:(nMode, nQ)
   try  {
     // ---read sxb file
      SxBinIO io (input, SxBinIO::BINARY_READ_ONLY);

      int nQ = io.getDimension ("nQ"),
          nMode = io.getDimension ("nMode");
      cout << "nMode" << nMode << "nQ" << nQ << endl;

      // make freq proper size
      freq.reformat (nMode, nQ);
      io.read ("frequencies", &freq, nMode, nQ);
      io.read ("volume", &vol);
   }  catch (SxException e)  {
      e.print();
      SX_QUIT;
   }
   cout << "input read" << endl;
   return freq;
}

SxArray<SxMatrix<Double> > SxThermo::getXtraPhon (
      const SxArray<SxMatrix<Double> > &allFreq,
      const SxArray<double> &vol,
      const bool lin) const
{
   cout << "getXtraPhon" << endl;

   // ---extrapolate phonons to other vols, linear or exponential
   int nTotVol = int(vol.getSize ()),
       iVol, nVol = int(allFreq.getSize ()),
       nFreq = int(allFreq(0).getSize ()),
       nRow = int(allFreq(0).nRows ());
   for (iVol = 1; iVol < nVol; iVol++)  {
      if (allFreq(iVol).getSize () != nFreq)  {
         cout << "ERROR in SxThermo::getXtraPhon\n"
              << "no. frequencies differ for the volumes!\n"
              << "vol no. " << iVol << " has " << allFreq(iVol).getSize ()
              << " frequencies, while volume 0 had " << nFreq << endl
              << "you MUST use identical qPoint grids!";
         SX_EXIT;
      }
   }
   // set average input volume (linear or logarithmic)
   double avVol = (lin ? SxVector<Double> (vol) (SxIdx(0, nVol - 1))
                       : logVec (SxVector<Double> (vol) (SxIdx(0, nVol - 1))))
                  .sum () / nVol;

   // set average frequencies (linear or logarithmic)
   SxMatrix<Double> avFreq (nRow, nFreq / nRow);
   avFreq.set (0.);
   if (lin)  {
      for (iVol = 0; iVol < nVol; iVol++)  { avFreq += allFreq(iVol); }
   } else  { //logarithmic
      for (iVol = 0; iVol < nVol; iVol++)  {
         avFreq += logVec (allFreq(iVol).abs ());
      }
   }
   avFreq /= nVol;

   // set average to 1st param
   SxArray<SxMatrix<Double> > freqParam (2);
   freqParam(0).copy (avFreq);

   // set average sqr vol or log sqr vol
   double vol2 = (lin ? SxVector<Double> (vol) (SxIdx(0, nVol - 1))
                      : logVec (SxVector<Double> (vol) (SxIdx(0, nVol - 1))))
                  .normSqr ();

   // set 1st moment of freq (linear or logarithmic)
   avFreq.set (0.);
   if (lin)  {
      for (iVol = 0; iVol < nVol; iVol++)  {
         avFreq += allFreq(iVol) * (vol(iVol) - avVol);
      }
   }  else  {//logarithmic
      for (iVol = 0; iVol < nVol; iVol++)  {
         avFreq += logVec (allFreq(iVol).abs ()) * (log (vol(iVol)) - avVol);
      }
   }

   // set 1st moment over vol sqr to 2nd param
   avFreq /= vol2;
   freqParam(1).copy (avFreq);

   // --extrapolate freq to the new volumes
   SxArray<SxMatrix<Double> > freqNew (allFreq);
   freqNew.resize (nTotVol, true);
   for (iVol = nVol; iVol < nTotVol; iVol++)
      { freqNew(iVol) = freqParam(0) + freqParam(1) * (vol (iVol) - avVol); }
   if (! lin)  { //logarithmic
      for (iVol = nVol; iVol < nTotVol; iVol++)
         { freqNew(iVol) = exp (freqNew(iVol)); } //*sign (allFreq(iVol)?!
   }
   cout << "done" << endl;
   return freqNew;
}

void SxThermo::printDOS (const SxArray<double> &vol,
                         const SxArray<SxMatrix<Double> > &freq) const
{
   cout << "printDOS" << endl;

   // check array sizes
   SX_CHECK (vol.getSize () == freq.getSize (),
             vol.getSize (), freq.getSize ());

   //--- calculate dosses.
   double dE = double(freq(0).getSize ()) / 400.; //energy interval (meV)
   int iVol, nVol = int(vol.getSize ()),
       nInt = int (40 * dE); //on average 10 freq per interval
   SxArray<SxArray<double> > dos (nVol);
   SxVector<Double>::Iterator freqIt;
   for (iVol = 0; iVol < nVol; iVol++)  {

      // make some checks on the input
      SX_CHECK (freq(iVol).getSize () > 0);
      SX_CHECK (vol(iVol) > 0.);

      dos(iVol).resize (nInt);
      dos(iVol).set (0);
      for (freqIt = freq(iVol).begin (); freqIt != freq(iVol).end (); freqIt++)
      {  if (*freqIt < 0.)  { dos(iVol)(0) += dE; }
         else if (*freqIt >= 40 - 1./dE)  { dos(iVol)(nInt - 1) += dE; }
         else { dos(iVol)(int (*freqIt * dE)) += dE; }
      } //dos ranges from 0-40 meV
   }
   cout << "dos determined" << endl;

   //--- write to output
   ofstream dosOut;
   dosOut.open ("phononDos.out");
   dosOut << "#phononDOS (1/meV/atom) for volumes (A^3)";
   for (iVol = 0; iVol < nVol; iVol++)
   {  dosOut << (vol(iVol) / cube (A2B)) << " "; }
   for (int iE = 0; iE < nInt; iE++)  {
      dosOut << endl << (iE + .5) / dE;
      for (iVol = 0; iVol < nVol; iVol++)
      { dosOut << "\t" << (3 * dos(iVol)(iE) / (double)freq(iVol).getSize ()); }
   }
   dosOut << endl;
   dosOut.close ();
   cout << "dos printed" << endl;
}

SxVector<Double> SxThermo::getFreeEn (const SxMatrix<Double> &freq,
                                      const int maxT,
                                      const double dT,
                                      const bool ignore) const
{
   cout << "getFreeEn" << endl;

   int nMode = int(freq.nRows ()),
       nQ  = (int)freq.nCols (), //!
       nT = int(maxT / dT) + 1; //no. of temperature steps
   cout << "nQ" << nQ << "nMode" << nMode << endl;
   double K2meV = KB * JL2HA * HA2EV * 1000., //scaling
          qWeight, kT;
   SxVector<Double> freeEn (nT),
                    boltzmannT;
   SX_CHECK (nQ > 0 && nMode > 0);

   // --zero point vibration, remove minusses
   qWeight = 3.e-3 / HA2EV / nQ / nMode; //scaling, includes meV->H
   freeEn(0) = freq.abs ().sum () * qWeight / 2.;
   // Warning: soft systems are too stable, but cannot rescale here

   // --temperature dependent term
   cout << "nT=" << nT << endl;
   for (int iT = 1; iT < nT; iT++)  {
      (cout << ".").flush ();
      kT = K2meV * iT * dT;

      // Boltzmann factor; minus signs of soft modes removed
      if (ignore)  { boltzmannT = exp (- freq / kT); }
      else { boltzmannT = exp (- freq.abs () / kT); }

      // start with zero point energy and add T dependent term
      freeEn(iT) = freeEn(0) + kT * (logVec (1. - boltzmannT)).sum () * qWeight;
   }
   cout << endl;
   cout << "done" << endl;
   return freeEn;  // H per atom
}

SxArray<SxVector<Double> > SxThermo::readXtraF (const SxString &input) const {

   cout << "read xtraF" << endl;

   // check non-empty input
   SX_CHECK (! (input == ""));
   SxList<SxVector<Double> > res; //:nVol,nT
   try  {
      // ---input to parser
      SxParser parser;
      SxParser::Table parTab = parser.read (input);

      // read free energy data
      SxSymbolTable *volPtr = NULL;
      for (volPtr = parTab->getGroup ("vol"); volPtr != NULL;
           volPtr = volPtr->nextSibling ("vol"))
      { res << SxVector<Double> (volPtr->get ("xtraF")->toList ()); }
   }  catch (SxException e)  {
      e.print();
      SX_QUIT;
   }

   return res;// in H per atom
}

void SxThermo::printF (const SxArray<SxVector<Double> > &allFreeEn,//:nVol,nT
                       const SxArray<double> &vol,
                       const double dT,
                             ofstream &output) const
{
   cout << "printF" << endl;
   int iVol, nVol = (int)vol.getSize ();

   // check amount of free energies
   SX_CHECK (nVol == allFreeEn.getSize ());

   //-- C_V is the second derivative * -T using iterators
   int iT, nT = int(allFreeEn(0).getSize ());
   double F0, F1;
   SxArray<SxVector<Double> > Cv (nVol);//:nT-2
   SxVector<Double>::Iterator freeIt,
                              CvIt;
   for (iVol = 0; iVol < nVol; iVol++)  {
      Cv(iVol).resize (nT - 2);
      freeIt = allFreeEn(iVol).begin ();
      F0 = *freeIt++;
      F1 = *freeIt++;
      CvIt = Cv(iVol).begin ();
      for (iT = 2; iT < nT; iT++)  {
         *CvIt++ = iT * (2. * F1 - *freeIt - F0);
         F0 = F1;
         F1 = *freeIt++;
      }
      Cv(iVol) /= dT * KB * JL2HA; //H
   }

   //-- write output
   output << "#No murn data or too few volume: only harmonic results\n"
          << "#(Vib.) free en. (meV/atom) & C_V (k_B/atom) vs. T (K)\n"
             "#for volumes (A^3)";
   for (iVol = 0; iVol < nVol; iVol++)
   {  output << (vol(iVol) / cube (A2B)) << "\t"; }
   for (iT = 0; iT < nT; iT++)  {
      output << endl << (iT * dT);
      for (iVol = 0; iVol < nVol; iVol++)  {
         output << "\t" << (allFreeEn(iVol)(iT) * HA2EV * 1000.) << "\t"
                << (iT > 0 && iT < nT - 1 ? Cv(iVol)(iT - 1) : 0);
      }
   }
   output << endl;
}

SxMurn SxThermo::readMurnDat (const SxString &input,
                              const double &BM,
                              const double &BMD) const
{
   cout << "read murnDat" << endl;

   //check that input not empty
   SX_CHECK (! (input == ""));

   // ---input to parser
   SxParser murnPars;
   SxParser::Table murnParTab = murnPars.read (input);

   // read murnaghan data
   SxSymbolTable *murnDat = murnParTab->getGroup ("murndata");
   SxArray<double> vols = SxArray<double> (murnDat->get ("volumes")->toList ());
   SxArray<double> ens = SxArray<double> (murnDat->get ("energies")->toList ());

   // check if we have enough data
   int nVol = int(vols.getSize ());
   SX_CHECK (nVol == ens.getSize ());
   if (nVol < 3)  {
      cout << "Error in SxThermo::readMurnDat! " << endl
           << "No. provided volumes is only " << nVol << endl
           << "At least 3 are required!" << endl;
      SX_EXIT;
   }  else if (nVol < 5)  {
      cout << "Warning from SxThermo::readMurnDat! "  << endl
           << "No. provided volumes is only " << nVol << endl
           << "Fits may be bad!" << endl;
   }
   // set SxMurn
   SxMurn inMurn (vols, ens);

   // compute fit
   inMurn.computeFit ();

   // --write some output
   cout << "inMurn.volume " << inMurn.minVolume << "inMurn.bulkMod ";

   //reset BM/BMD when provided
   if (BM > 0.)  {
      cout << "reset to " << BM;
      inMurn.bulkModulus = BM * cube (BOHRRADIUS) * JL2HA * 1.e+9;
   } else {
      cout << inMurn.getBulkModulus (); //GPa
   }
   if (BMD > 0.)  {
      cout << " BM' reset to " << BMD << endl;
      inMurn.bulkModulusDerivative = BMD;
   } else {
      cout << "inMurn.BM' " << inMurn.bulkModulusDerivative << endl;
   }
   return inMurn;
}

SxArray<Coord> SxThermo::paramFreeEn (
      const SxArray<SxVector<Double> > &allFreeEn,//:nVol,nT
      const SxArray<double> &vol,
      const bool linF) const //set only linear fit
{
   cout << "parametrise free energy" << endl;

   int iVol, nVol = int(vol.getSize ()),
       nT = int(allFreeEn(0).getSize ());
   SX_CHECK (nT > 0);
   SX_CHECK (nVol == allFreeEn.getSize ());

   // at least three volumes needed
   if (nVol < 3)  {
      cout << "please provide more volumes!";
      SX_EXIT;
   }

   //---analytically solve 2nd order fit for free en at all T
   // calculate moments of volume
   double avVol = SxVector<Double> (vol).sum () / nVol,
          avVol2 = (SxVector<Double> (vol) - avVol).normSqr () / nVol,
          avVol3 = (SxVector<Double> (vol) - avVol).cub ().sum () / nVol,
          avVol4 = (SxVector<Double> (vol) - avVol).sqr ().sqr ().sum () / nVol;
   if (linF)  {avVol3 = 0; avVol4 = 0;} //only linear fit

   // calculate average free energy
   SxVector<Double> avF (nT, 0.);
   for (iVol = 0; iVol < nVol; iVol++)  { avF += allFreeEn(iVol); }

   // scale
   avF /= nVol;

   // calc. 1st and 2nd moments of F
   SxVector<Double> avF1 (nT, 0.),
                    avF2 (nT, 0.);
   for (iVol = 0; iVol < nVol; iVol++)  {
      avF1 += (vol(iVol) - avVol) * (allFreeEn(iVol) - avF);
      if (! linF)  { avF2 += sqr (vol(iVol) - avVol) * (allFreeEn(iVol) - avF); }
   }
   // scale
   avF1 /= nVol;
   avF2 /= nVol;

   // -put everything together
   double den = avVol4 - sqr (avVol2) - sqr (avVol3) / avVol2;
   SxVector<Double> par0 = avF + (avF1 * avVol3 - avF2 * avVol2) / den,
                    par1 = (sqr (avVol3) * avF1 - avVol3 * avVol2 * avF2)
                         / den / sqr (avVol2) + avF1 / avVol2,
                    par2 = (avVol2 * avF2 - avVol3 * avF1) / den / avVol2,
                    dev (nT, 0.);

   // check linear fit
   if (linF)  {
      if (par2.normSqr () > 1.e-6)  {
         cout << "Error in SxThermo::paramFreeEn\n"
              << "linear fit has quadratic contribution!" << par2.norm ()
              << nVol << endl;
         SX_EXIT;
      }
   }
   //reset par: avVol -> 0.
   par0 += par2 * sqr (avVol) - par1 * avVol;
   par1 -= 2. * par2 * avVol;

   // calc. average sqrt deviation from fit
   for (iVol = 0; iVol < nVol; iVol++)  {
      dev += (allFreeEn(iVol) - par0 - par1 * vol(iVol) - par2 * sqr (vol(iVol))
             ).sqr ();
   }
   // print to standard out
   cout << "sqrt average deviation from fit:"
        << (sqrt (dev) / sqrt (double(nVol))) << endl;

   // set output: 0th, 1st & 2nd order parameters (mev, B^3)
   SxArray<Coord> param (nT);
   SxVector<Double>::Iterator par0It = par0.begin (),
                              par1It = par1.begin (),
                              par2It = par2.begin ();
   for (int iT = 0; iT < nT; iT++, par0It++, par1It++, par2It++)
   { param (iT) = Coord (*par0It, *par1It, *par2It); }
   return param;
}

SxArray<Coord> SxThermo::paramHeat (const SxArray<Coord> &paramF,//:nT
                                    const double dT) const
{
   cout << "ParamHeat" << endl;

   // volume parametrization of the heat capacity C_V (k_B/atom)
   int nT = int(paramF.getSize ());
   double KT = dT * KB * JL2HA; //K->H
   SxArray<Coord> Cv (nT - 2);

   // take second derivative * -T
   for (int iT = 1; iT < nT - 1; iT++)
   { Cv(iT-1) = iT * (2. * paramF(iT) - paramF(iT + 1) - paramF(iT - 1)) / KT; }
   return Cv;
}

SxArray<SxMurn> SxThermo::fitMurn (SxMurn &inMurn,
                     const SxArray<Coord> &paramF, //:nT
                     const SxArray<double> &vol,
                     const double &press) const
{
   cout << "fitMurn" << endl;

   //--- make all murn fits
   int iVol, nVol = int(vol.getSize ()),
       nT = int(paramF.getSize ());

   // check that there are enough volumes
   SX_CHECK (nVol > 5);
   SX_CHECK (nT > 0);
   double V;
   SxArray<double> energs (nVol);
   SxArray<SxMurn> allMurn (nT);
   for (int iT = 0; iT < nT; iT++)  {
      for (iVol = 0; iVol < nVol; iVol++)  {
         V = vol(iVol);

         // add pressure and vibrational contributions to total energy
         energs(iVol) = inMurn.getEnergy (vol(iVol)) + press * vol(iVol)
                      + (paramF(iT) ^ Coord (1., V, sqr (V)));
      }
      allMurn(iT) = SxMurn (vol, energs);
      allMurn(iT).computeFit ();
   }
   return allMurn;
}

double SxThermo::getEqVol (SxMurn &inMurn,
                     const Coord &paramF,
                     const double &press) const
{
   cout << "getEqVol ";

   /* analytical solution (assuming F almost linear) would be:
    * dE0/dV = BM0/BMD (1-(V0/V)^BMD), dE0/dV+dF/dV+press = 0
   double BM0 = inMurn.bulkModulus,
          BMD = inMurn.bulkModulusDerivative,
          num = 1. + BMD * (paramF(1) + press) / BM0,
          anV = V / (exp (log (num) / BMD) + 2. * paramF(2) * V * num/ BM0);
    */
   // find minVol numerically
   int iStep = 0;
   double E0, E = inMurn.getMinEnergy (), //current energy
          V = inMurn.getMinVolume (), //initial volume
          dV = .2, //volume step
          slope; //slope of free energy
    while (fabs (dV) > 5.e-7 && iStep < 1.e+5)  {// correct aLat within .01 A
      E0 = E;
      if (V + dV < 0.) { break; }
      V += dV;
      E = inMurn.getEnergy (V);
      slope = (E - E0) / dV + paramF(1) + 2. * paramF(2) * V;
      if ((slope + press) * dV > 0.) { dV *= -.5; }
      iStep++;
   }
   // restrict max. no. steps
   if (iStep == 1e5)  {
      cout << "Warning: minimum volume not found! Current volume:"
           << V << " step:" << dV << endl;
   }
   return V;
}

SxArray<SxMurn> SxThermo::paramMurn (SxMurn &inMurn,
                       const SxArray<Coord> &paramF,//:nT
                       const double &press) const
{
   cout << "paramMurn" << endl;

   int iT, nT = int(paramF.getSize ());
   SX_CHECK (nT > 0);
   double V, V0 = inMurn.getMinVolume (),
          BM, BM0 = inMurn.bulkModulus,
          BMD = inMurn.bulkModulusDerivative;
   SxArray<SxMurn> allMurn (nT);

   //--- assign minV, minE, bulkModulus & Derivative
   for (iT = 0; iT < nT; iT++)  {

      // set the equilibrium vol
      V = getEqVol (inMurn, paramF(iT), press);
      allMurn(iT).minVolume = V;

      // set the equilibrium Energy
      allMurn(iT).minEnergy = inMurn.getEnergy (V) + press * V
                            + (paramF(iT) ^ Coord (1., V, sqr (V)));

      // set the new BM (V d^2F/dV^2) from vol. change  E_tot + vib. contribution
      BM = BM0 * exp (BMD * log (V0 / V)) + 2. * paramF(iT)(2) * V;
      allMurn(iT).bulkModulus = BM;

      // set new BM' (dB/dp, this is taken constant in a murn fit)
      allMurn(iT).bulkModulusDerivative = BMD - 2. * paramF(iT)(2) * V / BM;
   }
   return allMurn;
}

SxArray<double> SxThermo::getCp (SxArray<SxMurn> &allMurn, //:nT
                            const SxArray<Coord> &Cv, //:nT-2
                            const double dT) const //:nT
{
   cout << "getCp" << endl;

   // check the no. temperatures
   int iT, nT = int(Cv.getSize ()) + 2;
   SX_CHECK (allMurn.getSize () == nT, allMurn.getSize (), nT);

   //--- get Cp from Cv @ V + lattice expansion
   double V, V0 = allMurn(0).minVolume,
          E, E0 = allMurn(0).minEnergy,
          KT = dT * KB * JL2HA; //K->H
   SxArray<double> Cp (nT - 2);
   for (iT = 1; iT < nT - 1; iT++)  {
       V = allMurn(iT).minVolume;
       E = allMurn(iT).minEnergy;
       Cp(iT - 1) = (Cv(iT - 1) ^ Coord (1., V, sqr (V))) - (iT - .5) / KT *
          (E + E0 - allMurn(iT - 1).getEnergy (V) - allMurn(iT).getEnergy (V0));
       E0 = E;
       V0 = V;
   }
   return Cp;
}

void SxThermo::printHeat (const SxArray<Coord> &Cv, //:nT-2
                          const SxArray<Coord> &CvXF, //:nT-2 or 0
                          const SxArray<double> &vol,
                          const SxArray<SxVector<Double> > &Cp, //:nPres,nT-2
                          const double &unitPress, //H per B^3
                          const double dT) const
{
   cout << "printHeat" << endl;

   // check no. temperatures
   bool xf = (CvXF.getSize () > 0);
   int iT, nT = int(Cv.getSize () + 2);

   // check no. temperatures
   SX_CHECK (nT == Cp(0).getSize () + 2,
             nT, Cp(0).getSize ());
   if (xf)  {SX_CHECK (CvXF.getSize () == nT - 2);}

   // check no. volumes
   int iVol, nVol = int(vol.getSize ());
   SX_CHECK (nVol > 0);

   // check no. pressures
   int iPres, nPres = int(Cp.getSize ());
   SX_CHECK (nPres > 0);

   // open output
   ofstream output;
   output.open ("C_p.out");

   //--- print C_p
   output << "#T (K), C_p (k_B/atom) for pressures " ;
   for (iPres = 0; iPres < nPres; iPres++)
   { output << (iPres * unitPress / cube (BOHRRADIUS) / JL2HA) << " "; }
   output << "Pa\n" ;
   for (iT = 1; iT < nT - 1; iT++)  {
      output << (iT * dT);
      for (iPres = 0; iPres < nPres; iPres++)  {
         output << "\t" << Cp(iPres)(iT - 1);
      }
      output << endl;
   }
   cout << "C_p printed" << endl;

   // close old, open new output
   output.close ();
   output.open ("C_V.out");

   //--- print C_V
   if (xf)  {output << "#T (K), C_V (k_B/atom), xtraF part for volumes ";}
   else  {output << "#T (K), C_V (k_B/atom) for volumes ";}
   for (iVol = 0; iVol < nVol; iVol++)
   {  output << (vol(iVol) / cube (A2B)) << " "; }
   output << "\\AA^3\n";
   double V;
   for (iT = 1; iT < nT - 1; iT++)  {
      output << (iT * dT);
      for (iVol = 0; iVol < nVol; iVol++)  {
         V = vol(iVol);
         output << "\t" << (Cv(iT - 1) ^ Coord (1., V, sqr (V)));
         if (xf)  {output << "\t" << (CvXF(iT - 1) ^ Coord (1., V, sqr (V)));}
      }
      output << endl;
   }
   output.close ();
   cout << "C_V printed\n";
}

void SxThermo::printThermo (SxArray<SxMurn> &allMurn,
                            SxArray<Coord> &param,
                            SxArray<Coord> &paramXF,
                                   const double dT,
                                   ofstream &output) const
{
   cout << "printThermo\n";
   bool xf = (paramXF.getSize () > 0);
   int nT = int(allMurn.getSize ());
   SX_CHECK (param.getSize () == nT);
   if (xf)  {SX_CHECK (paramXF.getSize () == nT);}

   output <<"#T(K) quasi-harm.freeEn.(meV/atom), vib. part,";
   if (xf)  {output << ", xtra part";}
   output << "; atom vol.(\\AA^3) \alpha (1e-5/K) BM (GPa) BM'\n";
   double expand, minV, minE;
   Coord Vc;
   for (int iT = 0; iT < nT; iT++)  {
      minV = allMurn(iT).minVolume;
      Vc = Coord (1., minV, sqr (minV));
      minE = allMurn(iT).minEnergy;
      output << (iT * dT) << "\t" << setprecision(16) << setw(8) << (minE * HA2MEV) << " ";

      // print also contributions to F
      if (xf)  {
          output << ((param(iT) - paramXF(iT)) ^ Vc * HA2MEV)
                << " " << (paramXF(iT) ^ Vc * HA2MEV) << "\t";
      }
      output << (param(iT) ^ Vc * HA2MEV)
             << "\t"   << (minV / cube (A2B)) << "\t";
      //linear expansion coefficient
      if (iT == 0)  { output << 0. << "\t"; }
      else  {
         expand = cbrt (minV / allMurn(iT - 1).minVolume) - 1.;
         output << (expand * 1.e+5 / dT) << "\t" ;
      }
      // BM and BM'
      output << allMurn(iT).getBulkModulus () << "\t"
             << allMurn(iT).bulkModulusDerivative << endl;
    }
}

#else /* SX_STANDALONE */

int main (int argc, char **argv)  {

   // command line interface: collect identical options
   SxCLI cli (argc, argv);
   SxArray<SxString> phononInput = cli.option("-p|--phonon","file",
                                   "(multiple option): phonon files").toList ();
   SxString murnInput = cli.option ("-m|--murn", "file",
                                    "murnaghan input file").toString (""),
            xtraFinput =  cli.option ("-xF|--xtraF", "file",
                             "input file with electronic/magnetic contribution")
                             .toString ("");
   bool fit = cli.option ("-f|--fitMurn", "make murn fit (T) instead of"
                          " numerical determination").toBool (),
        linF = cli.option ("--linF", "make linear fit F").toBool (),
        ignore = cli.option ("-i|--ignore", "ignore imag. freq. instead of "
                             "taking the absolute values").toBool (),
        lin = cli.option ("-l|--linear", "set linear instead of exponential"
                          "interpolation of the phonons.").toBool ();
   int maxT = cli.option ("-T|--maxT", "maximal temperature (K)").toInt (500);
   double dT = cli.option ("-dT|--dT", "temperature step (K)").toDouble (100);
   int nT = int(maxT / dT) + 1;
   double unitPress = cli.option ("-dP|--presStep", "pressure steps (bar)")
                                 .toDouble (.1),
          BM = cli.option ("-BM|--bulkModulus", "overwrite BM0 (GPa)")
                          .toDouble (0.),
          BMD = cli.option ("-bmd|--modulusDerivative", "overwrite BMD0")
                           .toDouble (0.);
   SxArray<double> extraVol = cli.option ("-xV|--extraVol", "comma separated"
               " list of volumes (B^3/atom) for which the phonons are"
               " inter/extra-polated. WARNING equal BZ-grids must be used"
               " for the phonons!").toDoubleList ();
   cli.finalize ();

   // check if input provided
   if (phononInput.getSize () == 0)  {
      cout << "Error in SxThermo: no input provided!\n";
      SX_EXIT;
   }

   // open output file
   ofstream output;
   output.precision(16);
   output.open ("thermo.out");
   output << "#Input files are " << phononInput << " " << murnInput << " "
          << xtraFinput << "\n#Maximal temperature " << maxT << "K, step " << dT
          << "\n#Pressure step " << unitPress << "bar; fitMurn?" << fit
          << "; ignore freq?" << ignore << endl;

   SxThermo T;

   // --read frequencies for different volumes
   int iVol, nVol = int(phononInput.getSize ());
   SxArray<double> vol (nVol);
   SxArray<SxMatrix<Double> > allFreq (nVol);//:(nMode, nQ)
   for (iVol = 0; iVol < nVol; iVol++)
   { allFreq(iVol) = T.readFreq (phononInput(iVol), vol(iVol)); }

   output << "#vol are (\\AA^3/atom) " ;
   for (iVol = 0; iVol < nVol; iVol++)
   { output << (vol(iVol) / cube (A2B)) << " "; }

   if (extraVol.getSize () > 0 && nVol > 1)  {

      //-- enlarge vol
      vol.resize (nVol + extraVol.getSize (), true);
      for (iVol = 0; iVol < extraVol.getSize (); iVol++)
      {  vol(iVol + nVol) = extraVol(iVol); }

      // reset nVol
      nVol = int(vol.getSize ());

      // extrapolate frequencies to other volumes (linear/exponential)
      // Warning: exponential extrapolation only to positive frequencies
      allFreq = T.getXtraPhon (allFreq, vol, lin); //auto resize
      output << "#XtraPhon for vols (A^3/atom): " ;
      for (iVol = 0; iVol < extraVol.getSize (); iVol++)
      { output << (extraVol(iVol) / cube (A2B)) << " "; }
   }  else  { output << "#No XtraPhon"; }
   output << endl;
   cout << "@vol are (B^3/atom)" << vol << endl;

   // print the phonon DOS for all volumes
   T.printDOS (vol, allFreq);

   // get the vibrational free energy for all vol (H/atom)
   SxArray<SxVector<Double> > allFreeEn (nVol);//:nT
   for (iVol = 0; iVol < nVol; iVol++)
   { allFreeEn(iVol) = T.getFreeEn (allFreq(iVol), maxT, dT, ignore); }

   // read extra free energy contribution when provided
   SxArray<SxVector<Double> > xtraF (0);
   if (! (xtraFinput == ""))  {

      //assume all vectors have same size
      xtraF = T.readXtraF (xtraFinput);

      // check if proper dimensions
      if (xtraF.getSize () != nVol || xtraF(0).getSize () != nT)  {
         cout << "Warning dimensions of xtraF not correct!" << endl
              << "nVol is " << xtraF.getSize () << " and should be " << nVol
              << "\n nT is " << xtraF(0).getSize () << " should be " << nT
              << "\n xtraF will be ignored\n";
      }  else  {

         //add xtraF to allFreeEn
         for (iVol = 0; iVol < nVol; iVol++) { allFreeEn(iVol) += xtraF(iVol); }
      }
   }

   // no murnaghan fit possible
   if (nVol < 3 || murnInput == "")  {

      // write free energies (& Cv) of input volumes to output
      T.printF (allFreeEn, vol, dT, output);

      //close output & end program
      output.close ();
      return 0;
   }

   // read murnaghan data (total energies & vols)
   SxMurn inMurn = T.readMurnDat (murnInput, BM, BMD);

   // parametrise the free energy with volume, contains the extraPhon
   SxArray<Coord> param = T.paramFreeEn (allFreeEn, vol, linF), paramXF (0);
   if (xtraF.getSize () > 0)  {paramXF = T.paramFreeEn (xtraF, vol, linF);}

   // parametrise the heat capacity with vol (per atom)
   SxArray<Coord> Cv = T.paramHeat (param, dT), CvXF (0);
   if (paramXF.getSize () > 0)  {CvXF = T.paramHeat (paramXF, dT);}

   if (fit)  {
      //--- reset the volumes for the murnaghan fits to the double of the range minVol -> maxVol)
      double minV = SxVector<Double> (vol).minval (),
             maxV = SxVector<Double> (vol).maxval ();
      vol.resize (10);
      for (iVol = 0; iVol < 10; iVol++)
      { vol(iVol) = (3. * minV - maxV) / 2. + (maxV - minV) * iVol * 2. / 9.; }
   }
   unitPress *=  1.e5 * cube (BOHRRADIUS) * JL2HA; //scale to H/B^3
   double press;
   SxArray<SxMurn> allMurn; //:nT
   SxArray<SxVector<Double> > Cp (5); //:nT-2
   for (int iPres = 4; iPres >= 0; iPres--)  {

      //-set and print pressure
      press = iPres * unitPress;
      cout << "pressure is" << press << endl;
      if (fit)  {
         // do murnaghan fits at all temperatures
         allMurn = T.fitMurn (inMurn, param, vol, press);
      } else {
         // parametrise the murn data
         allMurn = T.paramMurn (inMurn, param, press);
      }

      // get the heat capacity at constant pressure
      Cp(iPres) = T.getCp (allMurn, Cv, dT);
   }

   // reset the heat capacity volumes from the lattice expansion
   vol.resize (5);
   for (iVol = 0; iVol < 5; iVol++)
   { vol(iVol) = allMurn((iVol * (nT - 1)) / 4).minVolume; }

   // print the heat capacities
   T.printHeat (Cv, CvXF, vol, Cp, unitPress, dT);

   // write (other) thermodynamic results to output
   T.printThermo (allMurn, param, paramXF, dT, output);

   // close output & end program
   output.close ();
   return 0;
}
#endif /* SX_STANDALONE */
