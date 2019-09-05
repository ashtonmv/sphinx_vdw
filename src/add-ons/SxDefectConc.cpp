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


#include <SxDefectConc.h>

#include <SxTimer.h>

enum DefectTimers {
   GetCharge, ReactionEvolve
};

SX_REGISTER_TIMERS(DefectTimers)
{
   regTimer (GetCharge, "Get charge");
   regTimer (ReactionEvolve, "reaction");
}

#ifndef SX_STANDALONE
int SxDefectConc::findSpecies (const SxString &name) const
{
   int is = int(elementNames.findPos (name));
   if (is < 0) throw SxException
      (("Could not find element '" + name + "' in species list").ascii (),
       __FILE__, __LINE__);
   return is;
}

int SxDefectConc::findDefect (const SxString &name) const
{
   SxList<Defect>::ConstIterator it;
   int id = 0;
   for (it = begin (); it != end (); ++it, ++id)  {
      if ((*it).name == name) return id;
   }
   throw SxException
      (("Could not find defect'" + name + "' in defect list").ascii (),
       __FILE__, __LINE__);
   // return -1;
}

SxDefectConc::SxDefectConc (const SxSymbolTable *table)
{
   SX_CHECK (table);
   // SxStack<double> muList;
   SxList<double> muList; 
   double eMassByhBarSqr = 1.;
   try {
      // --- find kBoltzmann 
      if (table->contains ("kBoltzmann"))  {
         // directly ...
         kBoltzmann = table->get ("kBoltzmann")->toReal ();
      } else if (table->contains ("energyUnit"))  {
         // ... or via energy unit
         SxString unit = table->get ("energyUnit")->toString ();
         if (unit == "eV")           kBoltzmann = 8.61734216109122e-5;
         else if (unit == "meV")     kBoltzmann = 8.61734216109122e-2;
         else if (unit == "kJ/mol" ||
                  unit == "kJ")      kBoltzmann = 8.314472e-3;
         else if (unit == "kcal/mol" ||
                  unit == "kcal")    kBoltzmann = 1.985877519824e-3;
         else if (unit == "Hartree") kBoltzmann = 3.16681515923377e-6;
         else {
            cout << "Cannot interprete energyUnit=\"" << unit << "\".\n";
            SX_QUIT;
         }
      } else {
         throw SxException ("Need kBoltzmann or energyUnit");
      }
      
      temperature = table->contains("temperature")
                  ? table->get("temperature")->toReal ()
                  : 300.;
      siteConc    = table->contains("siteConcentration")
                  ? table->get ("siteConcentration")->toReal ()
                  : 1.;
      
      // --- read species and chemical potentials
      const SxSymbolTable *chemPot = table->getGroup ("chemicalPotentials");
      SxStack<int> mobileList;
      int is = 0;
      for (const SxSymbolTable *species = chemPot->getGroup ("species");
           species ; species = species->nextSibling ("species"), ++is)  {
         elementNames << species->get("element")->toString ();
         muList << species->get ("value")->toReal ();
         if (   species->contains("mobile") 
             && species->get("mobile")->toAttribute ())
            mobileList << is;
      }
      mu = muList;
      mobileSpecies = mobileList;

      const SxSymbolTable *defect = table->getGroup ("defect");
      for ( ; defect ; defect = defect->nextSibling ("defect"))
         defects << Defect(defect, *this);

      const SxSymbolTable *band = table->containsGroup ("band")
                                ? table->getGroup ("band")
                                : NULL;
      if (table->contains("eMassByhBarSqr"))
         eMassByhBarSqr = table->get("eMassByhBarSqr")->toReal ();
      else if (band)  {
         SxString eUnit = table->get("energyUnit")->toString ();
         if (eUnit == "eV")           eMassByhBarSqr = 0.036749326;
         else if (eUnit == "meV")     eMassByhBarSqr = 3.6749326e-5;
         else if (eUnit == "kJ/mol" ||
                  eUnit == "kJ")      eMassByhBarSqr = 3.8087989e-4;
         else if (eUnit == "kcal/mol" ||
                  eUnit == "kcal")    eMassByhBarSqr = 0.0015946679;
         else if (eUnit == "Hartree") eMassByhBarSqr = 1.;
         else {
            cout << "Unknown energy unit '" << eUnit << "'" << endl;
            cout << "Specify eMassByhBarSqr in units of length^{-2} energy^-1."
                 << endl;
            SX_EXIT;
         }
         SxString lUnit = table->get("lengthUnit")->toString ();
         if (lUnit == "m")       eMassByhBarSqr *= 3.5710649e+20;
         else if (lUnit == "cm") eMassByhBarSqr *= 3.5710649e+16;
         else {
            cout << "Unknown length unit '" << lUnit << "'" << endl;
            cout << "Specify eMassByhBarSqr in units of length^{-2} energy^-1."
                 << endl;
            SX_EXIT;
         }
      }
      for ( ; band ; band=band->nextSibling("band"))
         bands << Band(band, eMassByhBarSqr);

      dielecConstant = table->get("dielecConstant")->toReal ();
      if (table->contains("eps0"))  {
         dielecConstant *= table->get ("eps0")->toReal ();
      } else {
         dielecConstant /= FOUR_PI; // in atomic units
         SxString eUnit = table->get("energyUnit")->toString ();
         if (eUnit == "eV")           dielecConstant /= 27.211383;
         else if (eUnit == "meV")     dielecConstant /= 27211.383;
         else if (eUnit == "kJ/mol" ||
                  eUnit == "kJ")      dielecConstant /= 2625.4996;
         else if (eUnit == "kcal/mol" ||
                  eUnit == "kcal")    dielecConstant /= 627.08981;
         else if (eUnit == "Hartree") /* nothing */;
         else {
            cout << "Unknown energy unit '" << eUnit << "'" << endl;
            cout << "Specify eps0 in units of e^2 length^-1 energy^-1."
                 << endl;
            SX_EXIT;
         }
         SxString lUnit = table->get("lengthUnit")->toString ();
         if (lUnit == "m")       dielecConstant /= 5.2917721e-11;
         else if (lUnit == "cm") dielecConstant /= 5.2917721e-09;
         else {
            cout << "Unknown length unit '" << lUnit << "'" << endl;
            cout << "Specify eps0 in units of e^2 length^-1 energy^-1."
                 << endl;
            SX_EXIT;
         }
         cout << "eps=" << dielecConstant << endl;
         
      }
      if (table->containsGroup("reaction")) reactions.read (table,*this);
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
}

void SxDefectConc::print () const
{
   SxList<SxString>::ConstIterator specIt = elementNames.begin ();
   cout << "Species: " << (*specIt++);
   for (; specIt != elementNames.end (); ++specIt)
      cout << ", " << *specIt;
   cout << endl << "Defects: " << endl;
   SxList<Defect>::ConstIterator it = defects.begin ();
   for (; it != defects.end (); ++it)
      (*it).print (elementNames);
   cout << "Bands: " << endl;
   SxList<Band>::ConstIterator bIt = bands.begin ();
   for (; bIt != bands.end (); ++bIt)
      (*bIt).print (getkT ());

}

void SxDefectConc::Defect::print (const SxList<SxString> &elements) const
{
   cout << '"' << name << "\" (";
   if (multiplicity != 1) cout << multiplicity << "x; ";
   cout << "q=" << charge(0);
   for (int ics = 1; ics < charge.getSize (); ++ics)
      cout << ',' << charge(ics);
   cout << ';';
   for (int is = 0; is < sumFormula.getSize (); ++is)  {
      if (sumFormula(is) != 0)
         sxprintf (" %+d %s", (int)round(sumFormula(is)), elements(is).ascii ());
   }
   cout << ')' << endl;      
}

void SxDefectConc::printConcentrations (double eFermi) const
{
   double kT = getkT();
   cout << "Fermi level: " << eFermi << " (kT=" << kT << ")" << endl;
   cout << "mu = " << mu << endl;
   cout << endl;
   
   SxList<SxDefectConc::Defect>::ConstIterator it;
   int is, id;
   SxVector<Double> cS(getNSpecies ()), c;
   cS.set (0.);
   for (it = defects.begin (), id = 0; it != defects.end (); ++it, ++id)  {
      c = (*it).getConcentration (kT, mu, eFermi) * siteConc;
      cout << (*it).name << ": " << c << endl;
      cS += (*it).sumFormula * c.sum ();
   }
   cout << endl;
   cout << "Free carriers:" << getFreeChargeDensity (eFermi) << endl;
   cout << "Charge density:" << getChargeDensity (eFermi) << endl;

   // --- print total species concentrations
   for (is = 0; is < getNSpecies (); ++is)
      cout << "c(" << elementNames(is) << ")= " << cS(is) << endl;
}

SxDefectConc::Defect::Defect ()
   : multiplicity(0)
{
   // empty
}

SxDefectConc::Defect::Defect (const SxSymbolTable *table,
                              const SxDefectConc &parent)
{
   SX_CHECK (table);
   const SxString qLabel       = "charge",
                  eLabel       = "energy",
                  stateLabel   = "state",
                  multLabel    = "configurations",
                  speciesLabel = "species",
                  elemLabel    = "element",
                  nameLabel    = "name",
                  nAtomLabel   = "nAtoms";

   int nSpecies = parent.getNSpecies ();
   sumFormula.resize (nSpecies);
   sumFormula.set (0);
   try {
      int nChargeStates = table->getGroup (stateLabel)->getNItems (stateLabel),
          ics = 0;
      E0.resize (nChargeStates);
      charge.resize (nChargeStates);
      
      for (const SxSymbolTable *csTable = table->getGroup (stateLabel);
           csTable; csTable = csTable->nextSibling (stateLabel), ++ics)  { 
         E0(ics)      = csTable->get (eLabel)->toReal ();
         charge(ics)  = csTable->contains(qLabel) 
                      ? csTable->get(qLabel)->toReal ()
                      : 0.;
      }
      SX_CHECK (ics == nChargeStates, ics, nChargeStates);
      multiplicity = table->contains(multLabel)
                   ? table->get(multLabel)->toInt ()
                   : 1;
      name         = table->contains (nameLabel)
                   ? table->get(nameLabel)->toString ()
                   : "Defect " + SxString (parent.getNDefects () + 1);

      // --- read defect composition
      int is = -1;
      const SxSymbolTable *species = table->containsGroup (speciesLabel)
                                   ? table->getGroup (speciesLabel)
                                   : NULL;
      for (; species ; species = species->nextSibling (speciesLabel))  {
         if (species->contains (elemLabel))
            is = parent.findSpecies (species->get(elemLabel)->toString ());
         else
            ++is;
         if (is >= nSpecies)
            throw SxException (("Too many species for '" + name + "'").ascii (),
                  __FILE__, __LINE__);
         sumFormula(is) = species->get(nAtomLabel)->toInt ();
      }
      if (table->contains("diffusionBarrier"))  {
         diffusionBarrier   = table->get("diffusionBarrier")->toReal ();
         diffusionPrefactor = table->get("diffusionPrefactor")->toReal ();
         diffusionCharge    = table->get("diffusionCharge")->toReal ();
      } else {
         diffusionBarrier = 0.;
         diffusionPrefactor = -1.;
      }
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
}

SxDefectConc::Band::Band (const SxSymbolTable *table, double eMassByhBarSqr)
{
   SX_CHECK(table);
   SX_CHECK (eMassByhBarSqr > 0., eMassByhBarSqr);
   try {
      E0         = table->get("energy")->toReal ();
      mass       = table->get("mass")->toReal () * eMassByhBarSqr;
      name       = table->contains("name")
                 ? table->get("name")->toString ()
                 : SxString((mass > 0.) ? "electron" : "hole");
      degeneracy = table->contains("degeneracy")
                 ? table->get("degeneracy")->toInt ()
                 : 1;
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}

double SxDefectConc::Band::getChargeDensity (double eFermi, double kT) const
{
   double x = kT * fabs(mass) / TWO_PI;
   if ( (mass > 0.) && (E0-eFermi > kT))  {
         // electrons
         return  x * sqrt(x) * degeneracy * 2. // spin factor 2
                   * polyLog3_2 (-exp(-(E0 - eFermi) / kT));
   } else if (mass < 0. && (eFermi - E0 > kT))  {
         // holes
         return -x * sqrt(x) * degeneracy * 2. // spin factor 2
                   * polyLog3_2 (-exp((E0 - eFermi) / kT));
   }
   // --- integrate by hand
   // Fermi energy relative to band edge
   double eF = (mass > 0.) ? eFermi - E0 : E0 - eFermi;
   // starting energy where occupation = 1-1e-14
   double e = eF + kT * log(1e-14);
   if (e < 0.) e = 0.;
   // analytic integral with occupation=1
   double sum = pow(2. * fabs(mass) * e, 1.5) / (3. * PI_2);
   // numerically integrate partial occupation part
   double de = 1e-2 * kT;
   double eMax = eF - kT * log(1e-14);
   double beta = 1. / kT,
          fac  =  fabs(mass) * sqrt(2. * fabs(mass)) * de / (3. * PI_2);
   bool odd = true;
   sum -= fac * sqrt(e) / (1. + exp(beta * (e - eF)));
   for (; e < eMax || odd; e+=de, odd = !odd)
      sum += (odd ? 2. : 4.) * fac * sqrt(e) / (1. + exp(beta * (e - eF)));
   sum += fac * sqrt(e) / (1. + exp(beta * (e - eF)));
   SX_CHECK_NUM(sum);
   SX_CHECK (sum < 1e50 && sum > 0., sum);
   return degeneracy * (mass > 0. ? -sum : sum);
}

void SxDefectConc::Band::print (double kT) const
{
   cout << name << ": effective DOS (@E=" << E0 << "): ";
   double x = kT * fabs(mass) / TWO_PI;
   cout << 2. * degeneracy * x * sqrt(x) << endl;
}

SxVector<Double>
SxDefectConc::Defect::getConcentration (double kT,
                                        const SxVector<Double> mu,
                                        double fermiEnergy) const
{
   return double(multiplicity) * exp (-getEnergy(mu, fermiEnergy) / kT);
}

double SxDefectConc::getFreeChargeDensity (double eFermi) const
{
   SxList<Band>::ConstIterator it = bands.begin ();
   double conc = 0;
   for (; it != bands.end (); it++)
      conc += (*it).getChargeDensity (eFermi, getkT ());
   return conc;
}

double SxDefectConc::getChargeDensity (double eFermi) const
{
   double conc = 0., kT = getkT ();
   SxList<Defect>::ConstIterator it;
   for (it = defects.begin (); it != defects.end (); ++it)  {
      conc += ((*it).charge * (*it).getConcentration (kT, mu, eFermi)).sum ();
   }
   //cout << "eFermi = " << eFermi << "=>" << (conc * siteConc) << endl;
   return conc * siteConc + getFreeChargeDensity (eFermi);
}

double SxDefectConc::getChargeDensity (double eFermi,
                                       const SxVector<Double> &concDefect) const
{
   SX_CHECK (concDefect.getSize () == defects.getSize (),
             concDefect.getSize (),  defects.getSize ());
   double conc = 0., kT = getkT ();
   SxList<Defect>::ConstIterator it;
   int id;
   for (it = defects.begin (), id = 0; it != defects.end (); ++it, ++id)  {
      conc += (*it).getChargeDensity (kT, eFermi, concDefect(id));
   }
   return conc + getFreeChargeDensity (eFermi);
}

SxDefectConc::Gap SxDefectConc::getGapRange () const
{
   SX_CHECK (bands.getSize () >= 2, bands.getSize ());
   SxList<Band>::ConstIterator bIt = bands.begin ();
   Gap res((*bIt).E0, (*bIt).E0);
   double E;
   for (++bIt; bIt != bands.end (); ++bIt)  {
      if ((E = (*bIt).E0) < res.from) res.from = E;
      if ((E = (*bIt).E0) > res.to)   res.to = E;
   }
   return res;
}

double SxDefectConc::getEFermi () const
{
   Gap range = getGapRange ();
   range.from -= 10. * getkT ();
   range.to   += 10. * getkT ();

   double fMin = getChargeDensity (range.from),
          fMax = getChargeDensity (range.to);
   if (fMin * fMax > 0.)  {
      cout << "Can't find range in Fermi level search" << endl;
      SX_EXIT;
   }
   for (int i = 0; i < 200; i++)  {
      double eF = 0.5 * (range.from + range.to);
      double n = getChargeDensity (eF);
      if (fabs(n) < 1e-10) return eF;
      if (n < 0.)  {
         range.to = eF;
         fMax = n;
      } else {
         range.from = eF;
         fMax = n;
      }
      if (fabs(range.to - range.from) < 1e-10) return eF;
   }
   SX_EXIT;
}

SxVector<Double> 
SxDefectConc::getConcSpecies (const SxVector<Double> &defectConc,
                              const SxArray<bool> &active) const
{
   // --- get species concentrations
   SxVector<Double> cS(getNSpecies ()); 
   cS.set (0.);
   int nDefects = getNDefects (), id;
   SxList<Defect>::ConstIterator dIt = defects.begin ();
   for (id = 0; id < nDefects; ++id, ++dIt)
      if (active(id))
         cS += (*dIt).sumFormula * defectConc(id);
   return cS;
}

SxVector<Double> 
SxDefectConc::equilibrateDefects (double eFermi,
                                  const SxArray<bool> &active,
                                  const SxVector<Double> &initialConc,
                                  bool saveMu)
{
   SxVector<Double> cS = getConcSpecies (initialConc, active);
   SxVector<Double> chemPot = chemicalEquilibrium (eFermi, active, cS);
   if (saveMu) mu = chemPot;
   // --- return final defect distribution
   return getDefConc (eFermi, chemPot, active, initialConc);
}
   
SxVector<Double> 
SxDefectConc::getDefConc(double eFermi,
                         const SxVector<Double> &chemPot,
                         const SxArray<bool> &active,
                         const SxVector<Double> &initialConc) const
{
   int id, nDefects = int(defects.getSize ());
   SxList<Defect>::ConstIterator dIt;

   SxVector<Double> res(nDefects);
   double kT = getkT ();
   for (id = 0, dIt = defects.begin (); id < nDefects; ++id, ++dIt)  {
      if (active(id))  {
         res(id) = (*dIt).getConcentration (kT, chemPot, eFermi).sum ()
                 * siteConc;
      } else {
         res(id) = initialConc(id);
      }
   }
   return res;
}

SxVector<Double> 
SxDefectConc::chemicalEquilibrium (double eFermi,
                                   const SxArray<bool> &active,
                                   const SxVector<Double> &cSIn,
                                   const SxArray<int> ignoreSpecies) const
{
   int nDefects = getNDefects (), nSpecies = getNSpecies (), id, is;
   SxList<Defect>::ConstIterator dIt;
   double kT = getkT ();
   SxVector<Double> chemPot(nSpecies), cSrel(nSpecies), cS;
   chemPot.set (0.);
   cS.copy (cSIn);

   SxMatrix<Double> one(nSpecies,nSpecies);
   one.set (0.); for (is = 0; is < nSpecies; ++is) 
      one(is,is) = 1;

   SxMatrix<Double> deriv(nSpecies, nSpecies);
   SxVector<Double> dMu;

   // --- check which species are used at all
   SxArray<bool> speciesUsed(nSpecies);
   speciesUsed.set (false);
   for (id = 0, dIt = defects.begin (); id < nDefects; ++id, ++dIt)  {
      if (active(id))  {
         for (is = 0; is < nSpecies; ++is)  {
            if (fabs((*dIt).sumFormula(is)) > 1e-8)  {
               speciesUsed(is) = true;
            }
         }
      }
   }

   SxVector<Double> cDef0(nDefects);
   for (id = 0, dIt = defects.begin (); id < nDefects; ++id, ++dIt)  {
      if (active(id))  {
         cDef0(id) = (*dIt).getConcentration (kT, chemPot, eFermi).sum ()
                     * siteConc;
      }
   }
   
   for (int iIgnore = 0; iIgnore < ignoreSpecies.getSize (); ++iIgnore)  {
      is = ignoreSpecies(iIgnore);
      speciesUsed(is) = false;
      chemPot(is) = mu(is);
   }
   
   chemPot<<= mu;
   
   // --- set unused species concentrations to arbitrary value
   for (is = 0; is < nSpecies; ++is)
      if (!speciesUsed(is)) cS(is) = 1.;

   //cout << "eFermi = " << eFermi << endl;
   //cout << "species conc=" << cS << endl;
      
   // --- loop to find chemical potentials
   SxArray<SxVector<Double> > cSX(nDefects);
   int it;
   for (it = 0; it < 180; ++it)  {
      if (it == 90) {
         for (is = 0; is < nSpecies; ++is)
            if (speciesUsed(is)) chemPot(is) = 0.; 
      }
      //cout << endl << "it=" << (it+1) << endl;
      //cout << "mu=" << chemPot << endl;
      
      // --- determine defect and species concentrations
      cSrel.set (0.);
      for (id = 0, dIt = defects.begin (); id < nDefects; ++id, ++dIt)  {
         if (active(id))  {
            double conc = cDef0(id) * exp(dot((*dIt).sumFormula, chemPot) / kT);
            //cout << (*dIt).name << ": " << conc << endl;
            cSX(id) = (*dIt).sumFormula * conc / cS;
            cSrel += cSX(id);
         }
      }
      for (is = 0; is < nSpecies; ++is)
         if (!speciesUsed(is)) cSrel(is) = 1.;
      
      
      //cout << "dcs=" << cSrel << endl;
      //cout << "|dcs-1|=" << (cSrel-1.).norm () << endl;

      if ((cSrel-1.).normSqr () < 1e-16 * nSpecies) break;

      deriv.set (0.);
      double minVal = cSrel.minval (),
             maxVal = cSrel.maxval ();

      if (minVal < 0.)  {
         // --- wrong sign: move chemical potential into right direction
         dMu.resize (nSpecies);
         dMu.set (0.);
         for (is = 0; is < nSpecies; ++is)
            if (cSrel(is) < 0.)
               dMu(is)= (cS(is) > 0) ? 2.*kT : -2.*kT;
      } else if ((maxVal > 2 || minVal < 0.5))  {
         //cout << "|log(cRel)| min" << endl;
         // --- signs are ok, but ratio is large: minimize |log(cRel)|
         for (id = 0, dIt = defects.begin (); id < nDefects; ++id, ++dIt)  {
            if (active(id))  {
               for (is = 0; is < nSpecies; ++is)
                  deriv.colRef(is) += ((*dIt).sumFormula(is)/kT) * cSX(id) 
                                    / cSrel;
            }
         }
         // --- set unused species diagonal
         for (is = 0; is < nSpecies; ++is)  {
            if (!speciesUsed(is))  {
               for (int js = 0; js < nSpecies; ++js)
                  deriv(is,js) = deriv(js,is) = 0.;
               deriv(is,is) = 1.;
            }
         }
         //cout << "deriv" << deriv << endl;
         
         SxVector<Double> evalsSqr = deriv.eigenvalues ().absSqr ();
         //cout << "deriv.evals" << deriv.eigenvalues () << endl;
         double smallVal = 1e-4 * max(evalsSqr.maxval (), 1e-2/kT);
         if (evalsSqr.minval () > smallVal)  {
            dMu = -deriv.inverse () ^ log(cSrel);
         } else {
            double mag = smallVal;
            dMu = - ((deriv.transpose () ^ deriv) + mag * one).inverse () ^ deriv.transpose () ^ log(cSrel);
         }
         
         if (dMu.normSqr () > sqr(30. * kT * nSpecies))  {
            dMu *= 30. * kT * nSpecies / dMu.norm ();
         }
         //cout << "trial" << (log(cSrel) + (deriv ^ dMu)) << endl;
         if (dMu.normSqr () < 1e-8 * log(cSrel).normSqr ())  {
            cout << "dMu=" << dMu << endl;
            // dMu vanishes, but gradient does not: 
            // move mu into right direction
            cout << "failed to find finite dMu: fixed step" << endl;
            for (is = 0; is < nSpecies; ++is)  {
               if (fabs(cSrel(is) - 1.) > 1e-3)  {
                  dMu(is)= (cS(is) > 0.) ? 2.*kT : -2.*kT;
               }
            }
         }
      } else {
         //cout << "|1-cRel| min" << endl;
         // --- minimize |cRel - 1|
         for (id = 0, dIt = defects.begin (); id < nDefects; ++id, ++dIt)  {
            if (active(id))  {
               for (is = 0; is < nSpecies; ++is)
                  deriv.colRef(is) += ((*dIt).sumFormula(is)/kT) * cSX(id);
            }
         }
         // --- set unused species diagonal
         for (is = 0; is < nSpecies; ++is)  {
            if (!speciesUsed(is))  {
               for (int js = 0; js < nSpecies; ++js)
                  deriv(is,js) = deriv(js,is) = 0.;
               deriv(is,is) = 1.;
            }
         }
         //cout << "deriv" << deriv << endl;

         SxVector<Double> evalsSqr = deriv.eigenvalues ().absSqr ();
         //cout << "deriv.evals" << deriv.eigenvalues () << endl;
         double smallVal = 1e-4 * max(evalsSqr.maxval (), 1./kT);
         if (evalsSqr.minval () > smallVal)  {
            dMu = deriv.inverse () ^ (1. - cSrel);
         } else {
            double mag = 1e-2 * smallVal;
            dMu = ((deriv.transpose () ^ deriv) + mag * one).inverse () 
                  ^ (deriv.transpose () ^ (1. - cSrel));
         }
         
         //cout << "trial" << (cSrel + (deriv ^ dMu)) << endl;
         //cout << "now |dcs|=" << (1. - cSrel).norm ()
         //     << ";next |dcs|=" << (cSrel - 1. + (deriv ^ dMu)).norm () << endl;
         if (dMu.normSqr () < 1e-8 * (1. - cSrel).normSqr ()
             || (cSrel + (deriv ^ dMu) - 1.).norm () > 0.95 * (1. - cSrel).norm ())  {
            cout << "dMu=" << dMu << endl;
            // dMu vanishes, but gradient does not: 
            cout << "failed to find finite dMu: line search" << endl;
            double curvP, curvM, concP, concM, dFac;
            double minCurv = sqrt(evalsSqr.maxval ()) * 1e-4;
            SxVector<Double> searchDir(nSpecies);
            searchDir.set (0.);
            
            /*
            SxMatrix<Double>::Eigensystem eig = deriv.eigensystem ();
            double ovlp = -1.;
            for (is = 0; is < nSpecies; ++is)  {
               if (evalsSqr(is) < minCurv 
                   && (fabs(dot(searchDir, eig.vecs.colRef(is))) > ovlp))  {
                  searchDir = eig.vecs.colRef(is);
                  ovlp = fabs(dot(searchDir, eig.vecs.colRef(is)));
                  break;
               }
            }
            */
            searchDir = 1. - cSrel;
            //searchDir.normalize ();
            dMu = (5. * kT) * searchDir;
            int iMu;
            for (iMu = 1; iMu < 100; ++iMu)  {
               // --- line search
               curvP = 0.; curvM = 0.;
               dMu *= 1.2;
               for (id = 0, dIt = defects.begin (); id < nDefects;
                    ++id, ++dIt)  {
                  if (active(id))  {
                     concP = (*dIt).getConcentration (kT, chemPot + dMu, eFermi)
                             .sum () * siteConc;
                     concM = (*dIt).getConcentration (kT, chemPot - dMu, eFermi)
                             .sum () * siteConc;
                     //cout << (*dIt).name << ": " << conc << endl;
                     dFac = dot ((*dIt).sumFormula, searchDir)
                          * dot((*dIt).sumFormula / cS, searchDir);
                     curvP += dFac * concP;
                     curvM += dFac * concM;
                  }
               }
               //cout << "iMu=" << iMu
               //     << "\t " << curvP << '\t' << curvM << endl;
               if (fabs(curvP) > minCurv)  {
                  break;
               } else if (fabs(curvM) > minCurv) {
                  dMu = -dMu;
                  break;
               }

            }
            SX_CHECK (iMu < 100);
         }
      }
      //cout << "dMu=" << dMu << endl;

      for (is = 0; is < nSpecies; ++is)
         if (!speciesUsed(is))  dMu(is) = 0.;

      // --- update chemical potential
      chemPot += dMu;
   }
   if (it >= 180)  {
      const_cast<SxDefectConc*>(this)->mu = chemPot;
      printConcentrations (eFermi);
      cout << "failed to determine chemical equilibrium" << endl;
      cout << "eFermi = " << eFermi << endl;
      cout << "target concentrations:" << cS << endl;
      cout << "actual rel. concentrations:" << cSrel << endl;
      cout << "last mu: " << chemPot << endl;
      cout << "last dMu: " << dMu << endl;
      if ((cSrel-1.).normSqr () > 1e-4)  {
         SX_EXIT; 
      }
   }
   //cout << "it(chem)=" << it << endl;
   return chemPot;
   
}

void SxDefectConc::Reaction::computeConstants (const SxDefectConc &defconc,
                                               double eFermi, double ekt)
{
   SX_CHECK (defconc.getNDefects () == coefficients.getSize (),
             defconc.getNDefects (), coefficients.getSize ());
   int iDef;
   SxList<Defect>::ConstIterator def = defconc.begin ();
   rateConstantForward = rateConstantBackward 
                       = prefactor * exp(-(barrier + charge * eFermi)/ekt)
                       * defconc.siteConc;
   for (iDef = 0; iDef < defconc.getNDefects (); ++iDef, ++def)  {
      if (fabs(coefficients(iDef)) > 1e-6)  {
         double tot = 0.;
         const SxVector<Double> &E0 = (*def).E0,
                                &q = (*def).charge;
         for (int i = 0; i < E0.getSize (); ++i)
            tot += exp( -(E0(i) + q(i) * eFermi)/ekt);
         tot *= (*def).multiplicity * defconc.siteConc;
         double rf = pow(tot, -fabs(coefficients(iDef)));
         if (coefficients(iDef) > 0.)
            rateConstantBackward *= rf;
         else
            rateConstantForward *= rf;
      }
   }
}

SxVector<Double> SxDefectConc::Reaction::getRate (const SxVector<Double> &conc)
{
   SX_CHECK (conc.getSize () == coefficients.getSize (),
             conc.getSize (), coefficients.getSize ());
   double rateForward = rateConstantForward,
          rateBackward = rateConstantBackward;
   for (int iDef = 0; iDef < conc.getSize (); ++iDef)  {
      if (coefficients(iDef) > 1e-6)
         rateBackward *= pow(conc(iDef), coefficients(iDef));
      else if (coefficients(iDef) < -1e-6)
         rateForward *= pow(conc(iDef), fabs(coefficients(iDef)));
   }
   return (rateForward - rateBackward) * coefficients;
}



void SxDefectConc::ReactionSystem::read (const SxSymbolTable *table,
                                         const SxDefectConc &defconc)
{
   SxList<Reaction> reacList;
   int nDef = defconc.getNDefects ();
   SxVector<Double> speciesSum(defconc.getNSpecies ());
   try  {
      const SxSymbolTable *reacTable = table->getGroup("reaction");
      for (; reacTable; reacTable = reacTable->nextSibling ("reaction"))  {
         Reaction reac;
         reac.coefficients = SxVector<Double> (nDef);
         reac.coefficients.set (0.);
         reac.charge = reacTable->get("transitionCharge")->toReal ();
         reac.barrier = reacTable->get("barrierEnergy")->toReal ();
         reac.prefactor = reacTable->get("prefactor")->toReal ()
                        / defconc.siteConc; // see below

         speciesSum.set (0.);
         cout << "Reaction:\nEducts:";
         for (const SxSymbolTable *educt = reacTable->getGroup ("educt");
              educt; educt = educt->nextSibling ("educt")) {
            SxString name = educt->get("name")->toString ();
            cout << '"'<< name << "\" ";
            int id = defconc.findDefect (name);
            reac.coefficients (id) -= 1.;
            speciesSum -= defconc.defects (id).sumFormula;
            reac.prefactor *= defconc.siteConc; // dimensionless prefactor
         }
         cout << endl << "Products:";
         for (const SxSymbolTable *product = reacTable->getGroup ("product");
              product; product = product->nextSibling ("product")) {
            SxString name = product->get("name")->toString ();
            cout << '"'<< name << "\" ";
            int id = defconc.findDefect (name);
            reac.coefficients (id) += 1.;
            speciesSum += defconc.defects (id).sumFormula;
         }
         cout << endl;
         if (speciesSum.norm () > 1e-6)  {
            cout << "Reaction is not balanced!" << endl;
            cout << "Species balance:" << speciesSum << endl;
            SX_QUIT;
         }
         reacList << reac;
      }
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   reactions = reacList;
}

bool SxDefectConc::ReactionSystem::setup (const SxVector<Double> &concIn,
                                          const SxDefectConc &defconc,
                                          double eFermi,
                                          double ekt)
{
   conc.copy (concIn);
   for (int id = 0; id < conc.getSize (); ++id) 
      if (conc(id) < 0.) return false;
   for (int i = 0; i < reactions.getSize (); ++i)  {
      reactions(i).computeConstants (defconc, eFermi, ekt);
   }
   return true;
}


void SxDefectConc::ReactionSystem::evolve (double time)
{
   SX_CLOCK (ReactionEvolve);
   SX_CHECK (time > 0., time);

   /*
   currentTime = 0.;
   cout << currentTime;
   for (int id = 0; id < conc.getSize (); ++id)
      cout << '\t' << conc(id);
   cout << endl;
   */
#ifndef NDEBUG
   for (int id = 0; id < conc.getSize (); ++id) 
      SX_CHECK (conc(id) >= 0., id, conc(id));
#endif

   evolve (time, 1);
}

void SxDefectConc::ReactionSystem::evolve (double dt, int nSteps)
{
   SX_CHECK (nSteps > 0, nSteps);
   SX_CHECK (dt > 0., dt);
   int id, nDef = int(conc.getSize ());
   SxVector<Double> rate, newRate, newConc;
   newConc.copy (conc);
   for (int i = 0; i < nSteps; ++i)  {
      // --- get reaction rates
      if (newRate.getSize () == 0)  {
         rate.resize (nDef);
         rate.set (0.);
         for (int ir = 0; ir < reactions.getSize (); ++ir)
            rate += reactions(ir).getRate (conc);
      } else {
         rate = newRate;
         newRate = SxVector<Double> ();
      }

      // compute new concentrations
      newConc.plus_assign_ax (dt, rate);
      
      // --- check if any concentration drops below zero ...
      for (id = 0; id < nDef; ++id)
         if (newConc(id) < 0.) break;
      if (id < nDef)  {
         // ... finer time step
         evolve (0.1 * dt, 10);
         newConc.copy (conc);
         continue;
      } 

      // --- compute rate constant at new conc
      newRate.resize (nDef);
      newRate.set (0.);
      for (int ir = 0; ir < reactions.getSize (); ++ir)
         newRate += reactions(ir).getRate (newConc);
      
      // --- check if any relevant rate changes by more than 5 %
      for (id = 0; id < nDef; ++id)  {
         bool bigRateChange = fabs(newRate(id) - rate(id)) 
                              > 0.025 * (fabs(newRate(id)) + fabs(rate(id)));
         bool relevantChange = fabs(rate(id) * dt) > 1e-6 * conc(id);
         if (bigRateChange && relevantChange)  {
            //cout << "# " << id << '\t' << rate(id) << '\t' << newRate(id)
            //     << '\t' << (rate(id) * dt / conc(id)) << endl;
            break;
         }
      }
      if (id < nDef)  {
         // ... finer time step
         evolve (0.2 * dt, 5);
         newConc.copy (conc);
         newRate.resize (0);
         continue;
      }
      /*
      currentTime += dt;
      cout << currentTime;
      for (id = 0; id < nDef; ++id)
         cout << '\t' << newConc(id);
      cout << endl;
      */
      // write concentrations to conc
      conc <<= newConc;
   }

}

double SxDefectConc::equilibrate (const SxSymbolTable *cmd)
{
   //SX_CHECK (cmd->getName () == "equilibrate");
   int is, nSpecies = getNSpecies ();
   SxVector<Double> concS(nSpecies);
   concS.set (0.);
   SxString element;
   SxStack<int> ignoreList;
   double eFermi;
   try {
      const SxSymbolTable *speciesGroup;
      for (speciesGroup = cmd->getGroup ("species");
           speciesGroup; 
           speciesGroup = speciesGroup->nextSibling ("species"))
      {
         element = speciesGroup->get("element") ->toString ();
         is = findSpecies (element);
         if (is < 0)  {
            cout << "Unknown element " << element 
                 << " in equilibrate." << endl;
            SX_QUIT;
         }
         if (speciesGroup->contains("concentration"))  {
            concS(is) = speciesGroup->get("concentration")->toReal ();
         } else {
            mu(is) = speciesGroup->get("mu")->toReal ();
         }
      }
      for (is = 0; is < nSpecies; ++is)  {
         if (concS(is) == 0.)  {
            ignoreList << is;
         }
      }
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   int nDefects = int(defects.getSize ());
   SxArray<int> ignoreSpecies(ignoreList);
   SxArray<bool> active (nDefects);
   active.set (true);
   if (cmd->contains("eFermi"))  {
      eFermi = cmd->get("eFermi")->toReal ();
      mu = chemicalEquilibrium (eFermi, active, concS,
                                                ignoreSpecies);
   } else {
      Gap gap = getGapRange ();
      double eMin = gap.from - 10. * getkT (),
             eMax = gap.to + 10. * getkT (), E;
      try {
         if (cmd->contains ("fermiRange"))  {
            SxList<double> range = cmd->get("fermiRange")->toList ();
            if (range.getSize () != 2)  {
               cout << "Invalid fermiRange: " << range << endl;
               SX_EXIT;
            }
            eMin = range(0);
            eMax = range(1);
         }
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
      double rho;
      do {
         E = 0.5 * (eMin + eMax);
         mu = chemicalEquilibrium (E, active, concS,
                                                   ignoreSpecies);
         rho = getChargeDensity (E);
         if (rho > 0.)  {
            eMin = E;
         } else  {
            eMax = E;
         }
      } while ( (eMax - eMin) > 1e-8);
      eFermi = E;
   }
   cout << SX_SEPARATOR;
   printConcentrations (eFermi);
   return eFermi;
}

#else

/** \brief Solve Poisson equation in 1D

    @param density density
    @param eps    absolute dielectric constant,
                  i.e. \f$\varepsilon_r \varepsilon_0\f$
    @param dV     voltage drop (units: [energy]/[e])
    @param dl     step length (= l/n)
  
  */
SxVector<Double> solvePoisson (const SxVector<Double> &density,
                               const SxVector<Double> &eps,
                               double dV,
                               double dl)
{
   double prefactor = dl * dl;

   int n = int(density.getSize ());
   SX_CHECK (eps.getSize () == n, eps.getSize (), n);
   int i;
   SxVector<Double> res(n);
   SxVector<Double>::Iterator nIt, epsIt, resIt;
   nIt   = density.begin ();
   epsIt = eps.begin ();
   resIt = res.begin ();
   // integrate differential equation
   // d/dz 
   double field = 0., pot = 0.;
   for (i = 0; i < n; ++i)  {
      pot   += field;
      *resIt++ = pot;
      field -= prefactor * *nIt++ / *epsIt++;
   }
   // meet boundary conditions
   field = (dV-pot) / (n - 1);
   resIt = res.begin ();
   for (i = 0; i < n; ++i)
      *resIt++ += field * i;

   return res;
}

enum AdjustMode { KeepField, KeepInitial };
SxVector<Double> adjustPotential (const SxDefectConc &defects,
                                  const SxMatrix<Double> &cDefect,
                                  SxVector<Double> *Vptr,
                                  enum AdjustMode mode)
{
   SX_CHECK (Vptr);
   double step = 1.;
   int n = int(Vptr->getSize ()), id;
   SX_CHECK (cDefect.nCols () == n, cDefect.nCols (), n);
   SX_CHECK (cDefect.nRows () == defects.getNDefects (),
             cDefect.nRows (), defects.getNDefects ());
      
   SxVector<Double> density(n);
   SxList<SxDefectConc::Defect>::ConstIterator it;

   double qTot = 0., qTotOld, dqBydV;
   int iLoop;
   for (iLoop = 0; iLoop < 1000; ++iLoop)  {
      //cout << "V(0)=" << (*Vptr)(0);
      //cout << "V(n-1)=" << (*Vptr)(n-1);
      // --- compute new density with present deltaV
      for (int i = 0; i < n; ++i)  {
         density(i) = defects.getFreeChargeDensity ((*Vptr)(i));
      }
      for (it = defects.begin (), id = 0;
           it != defects.end ();
           ++it, ++id)  {
         for (int i = 0; i < n; ++i)  {
            density(i) += (*it).getChargeDensity (defects.getkT (), (*Vptr)(i),
                                                  cDefect(id, i));
         }
      }

      qTotOld = qTot;
      qTot = density.sum ();
      //cout << " => qTot = " << qTot << endl;
      if (iLoop == 0)  {
         step = -defects.getkT () * qTot 
              / sqrt(density.normSqr () * double(density.getSize ()));
      } else {
         // --- adjust potential to enforce charge neutrality
         double maxStep = 2. *fabs(step);
         //if (maxStep > 0.1) maxStep = 0.1 * pow(0.99, double(iLoop % 10));
         dqBydV = (qTot-qTotOld)/step;
         //cout << "dq/dV=" << dqBydV << endl;
         step = -qTot / dqBydV;
         if (fabs(step) > maxStep)
            step *= maxStep / fabs(step);
      }
      if (mode == KeepField)  {
         *Vptr += step;
      } else if (mode == KeepInitial) {
         double f = step/n;
         for (int i = 0; i < n; ++i)
            (*Vptr)(i) += i * f;
      } else {
         SX_EXIT;
      }
      if (fabs(step) < 1e-8) break;
   }
   SX_CHECK (iLoop < 1000);
   return density;
}

class ProfileCharge
{
   public:
      enum Status { OK, Overflow, Underflow } status;
      double charge;
      ProfileCharge (Status stat = OK) : status(stat), charge(0.) {}
      ProfileCharge (double c) : status(OK), charge(c) {}
      bool positive ()  {
         return status == Underflow || (status == OK && charge > 0.); 
      }
      void print () const 
      {
         if (status == OK)
            cout << charge << endl;
         else
            cout << ((status == Underflow) ? "Underflow" : "Overflow") << endl;
      }
};

/** \brief   
  */
ProfileCharge getCharge (const SxDefectConc &defects,
                         const SxMatrix<Double> &cDefect,
                         SxVector<Double> *Vptr,
                         const SxVector<Double> &eps,
                         double initialField,
                         double dl,
                         int startAt)
{
   SX_CLOCK(GetCharge);
   SX_CHECK (Vptr);
   SxVector<Double> &V = *Vptr;
   double prefactor = dl * dl;
   SxDefectConc::Gap gap = defects.getGapRange ();
   //gap.from -= 10. * defects.getkT ();
   //gap.to   += 10. * defects.getkT ();

   int n = int(V.getSize ());
   SX_CHECK (eps.getSize () == n, eps.getSize (), n);
   int i,id;
   // integrate differential equation
   // d/dz 
   double field = initialField, pot = V(startAt) - initialField, density;
   double charge = 0.;
   SxList<SxDefectConc::Defect>::ConstIterator it;
   for (i = startAt; i < n; ++i)  {
      pot += field;
      V(i) = pot;
      if (pot < gap.from) return ProfileCharge::Underflow;
      if (pot > gap.to)   return ProfileCharge::Overflow;
      density = defects.getFreeChargeDensity (V(i));
      for (it = defects.begin (), id = 0; it != defects.end (); ++it, ++id)  {
         density += (*it).getChargeDensity (defects.getkT (), V(i),
                                            cDefect(id, i));
      }
      // cout << i << '\t' << pot << '\t' << field << '\t' << density << endl;
      charge += density;
      field -= prefactor * density / eps(i);
   }
   return ProfileCharge (charge);
}

void findPinnedPotential(SxVector<Double> *VPtr,
                            const SxVector<Double> &eps,
                            double V0,
                            double dl,
                            SxDefectConc &defects,
                            const SxMatrix<Double> &cDefect,
                            int startAt = 0,
                            double addCharge = 0.)
{
   SX_CHECK(VPtr);
   SxVector<Double> &V = *VPtr;
   SX_CHECK (cDefect.nCols () == V.getSize (),
             cDefect.nCols (), V.getSize ());
   SX_CHECK (cDefect.nRows () == defects.getNDefects (),
             cDefect.nRows (), defects.getNDefects ());
   int iLast = int(V.getSize ()) - 1;

   ProfileCharge res;
   
   // --- get 1st point: initialField = 0
   double startField = 0.;
   if (startAt > 0) startField = V0 - V(startAt-1);
   V(startAt) = V0;
   res = getCharge (defects, cDefect, VPtr, eps, startField, dl, startAt);
   if (startAt > 0) res.charge += addCharge;
   bool posField = res.positive ();
   //cout << startField << ": "; res.print ();
   
   // --- get 2nd point
   double iniField;
   if (fabs(startField) > 1e-6)  {
      iniField = startField + fabs(startField) * (posField ? 1e-1 : -1e-1);
   } else {
      iniField = (posField ? 1. : -1.) / double(V.getSize ());
      if (startAt > 0) iniField *= 1e-5;
   }

   int it;
   for (it = 1000; it; --it)  {
      V(startAt) = V0;
      res = getCharge (defects, cDefect, VPtr, eps, iniField, dl, startAt);
      if (startAt > 0) res.charge += addCharge;
      //cout << iniField << ": "; res.print ();
      if (posField)  {
         if (!res.positive ()) break;
      } else {
         if (res.positive ()) break;
      }

      if (fabs(startField) < 1e-6)
         iniField *= 2.;
      else
         iniField += (iniField - startField);

   } 
   if (it == 0)  {
      cout << "Failed to converge potential" << endl;
      SX_EXIT;
   }

   // --- set start interval
   double lower, upper;
   if (posField)  {
      lower = startField; upper = iniField;
   } else {
      lower = iniField; upper = startField;
   }


   double lowerV, upperV;
   bool lowerOK, upperOK;
   V(startAt) = V0;
   res = getCharge (defects, cDefect, VPtr, eps, lower, dl, startAt);
   if (startAt > 0) res.charge += addCharge;
   lowerOK = res.status == ProfileCharge::OK;
   lowerV = V(iLast);
   SX_CHECK (res.positive ());

   V(startAt) = V0;
   res = getCharge (defects, cDefect, VPtr, eps, upper, dl, startAt);
   if (startAt > 0) res.charge += addCharge;
   upperOK = res.status == ProfileCharge::OK;
   upperV = V(iLast);
   SX_CHECK (!res.positive ());

   // --- root finding in interval 
   do {
      V(startAt) = V0;
      iniField = 0.5 * (upper + lower);
      res = getCharge (defects, cDefect, VPtr, eps, iniField, dl, startAt);
      if (startAt > 0) res.charge += addCharge;
      /*
      cout << iniField << ": ";
      if (res.status == ProfileCharge::OK)
         cout << res.charge << endl;
      else  {
         cout << ((res.status == ProfileCharge::Underflow)
                  ? "Underflow"
                  : "Overflow");
         cout << endl;
      }
      */

      if (res.positive ())  {
         lower = iniField;
         lowerOK = (res.status == ProfileCharge::OK);
         lowerV = V(iLast);
      } else {
         upper = iniField;
         upperOK = (res.status == ProfileCharge::OK);
         upperV = V(iLast);
      }
      //cout << '(' << lower << "," << upper << ')' << endl;
   } while ((upper - lower) > 1e-10 * fabs(iniField));

   // --- check that V is converged until end, or restart where it diverges
   if (!lowerOK || !upperOK || fabs(lowerV - upperV) > 1e-8)  {
      SxVector<Double> lowV(V.getSize ());
      lowV(startAt) = V0;
      getCharge (defects, cDefect, &lowV, eps, lower, dl, startAt);
      V(startAt) = V0;
      getCharge (defects, cDefect, VPtr, eps, upper, dl, startAt);

      // --- find point where upper / lower potentials start to diverge
      int i, id;
      SxList<SxDefectConc::Defect>::ConstIterator defIt;
      SxDefectConc::Gap gap = defects.getGapRange ();
      for (i = startAt + 1; i < V.getSize (); ++i)  {
         if (fabs(lowV(i) - V(i)) > 1e-8) break;
         if (lowV(i) < gap.from) break;
         if (V(i) > gap.to) break;
         V(i) = 0.5 * (lowV(i) + V(i));
         // --- update total charge until here
         addCharge += defects.getFreeChargeDensity (V(i));
         for (defIt = defects.begin (), id = 0;
              defIt != defects.end (); ++defIt,++id)  {
            addCharge += (*defIt).getChargeDensity (defects.getkT (), V(i),
                                                    cDefect(id, i));
         }
      }
      if (i + 1 == V.getSize ())  {
         if (i == startAt + 1)  {
            cout << V << endl;
            cout << "Failed to converge potential!" << endl;
            SX_EXIT;
         }
         i--;
      }
      lowV.resize (0);

      /*
      cout << "Restart at " << i << ':' << V(i) 
           << '/' << (V(i) - V(i-1))
           << endl;
      */

      // compute potential from here on
      findPinnedPotential (VPtr, eps, V(i), dl, defects, cDefect, i, addCharge);
   }
}

void computePinnedPotential(SxVector<Double> *VPtr,
                            double V0,
                            double dl,
                            SxDefectConc &defects,
                            const SxMatrix<Double> &cDefect)
{
   SX_CHECK(VPtr);
   SxVector<Double> &V = *VPtr;
   SX_CHECK (cDefect.nCols () == V.getSize (),
             cDefect.nCols (), V.getSize ());
   SX_CHECK (cDefect.nRows () == defects.getNDefects (),
             cDefect.nRows (), defects.getNDefects ());

   SxVector<Double> eps(V.getSize ());
   eps = defects.dielecConstant;

   ProfileCharge res;
   findPinnedPotential (VPtr, eps, V0, dl, defects, cDefect);
}

void distributeH (const SxVector<Double> &V,
                  SxDefectConc &defects,
                  SxMatrix<Double> &cDefect)
{
   SxArray<bool> active(defects.getNDefects ());
   active.set (true);
   double kT = defects.getkT ();
   SxVector<Double> cS, chemPot;
   double cHOld = -1., cHtarget = 0., dMu = 0.;
   for (int icl = 0; icl < 100; ++icl)  {
      double cH = 0.;
      for (int ir = 0; ir < V.getSize (); ++ir)  {
         cS = defects.getConcSpecies (cDefect.colRef (ir), active);
         if (icl == 0)  {
            cHtarget += cS(defects.mobileSpecies(0));
         }
         chemPot = defects.chemicalEquilibrium (V(ir), active, cS, defects.mobileSpecies);

         cDefect.colRef (ir) <<= defects.getDefConc (V(ir), chemPot, active, cDefect.colRef(ir));
         cS = defects.getConcSpecies (cDefect.colRef (ir), active);
         cH += cS(defects.mobileSpecies(0));
      }
      cout << "nH/nHtarget - 1 = " << (cH/cHtarget - 1.) << endl;
      if (fabs(cH/cHtarget - 1.) < 1e-7) break;
      if (icl == 0)  {
         dMu = (cH > cHtarget) ? - kT : kT;
      } else {
         double dcdmu = (cH - cHOld) / dMu;
         double maxdMu = 2. * fabs(dMu);
         cout << "dc/dMu = " << dcdmu << endl;
         dMu = (cHtarget - cH) / dcdmu;
         if (fabs(dMu) > maxdMu) dMu *= maxdMu/ fabs(dMu);
      }
      cout << "dMu = " << dMu << endl;
      defects.mu(defects.mobileSpecies(0)) += dMu;
      cHOld = cH;
   }
   cout << "mu=" << defects.mu << endl;
}

#include <SxRhoMixer.h>
#include <SxRho.h>
void computePotential(SxVector<Double>       *VPtr,
                      SxVector<Double> *densityPtr,
                      double dV,
                      double dl,
                      SxDefectConc &defects,
                      const SxMatrix<Double> &cDefect,
                      double maxRes,
                      enum AdjustMode mode)
{
   SX_CHECK(VPtr);
   SX_CHECK(densityPtr);
   SxVector<Double> &V = *VPtr, &density = *densityPtr;
   SX_CHECK (V.getSize () == density.getSize (),
             V.getSize (), density.getSize ());
   SX_CHECK (cDefect.nCols () == V.getSize (),
             cDefect.nCols (), V.getSize ());
   SX_CHECK (cDefect.nRows () == defects.getNDefects (),
             cDefect.nRows (), defects.getNDefects ());

   SxArray<bool> active(defects.getNDefects ());
   active.set (true);


   SxVector<Double> eps(V.getSize ());
   eps = defects.dielecConstant;

   SxVector<Double> newV, R, Rold;
   double residue;
   SxRhoMixer mixer(SxRhoMixer::Pulay, 0.8, 40);
   //mixer.preconditioner.scaling = 0.3;
   mixer.setNormModus (SxRhoMixer::RenormOff);
   VALIDATE_VECTOR(V);
   density = adjustPotential (defects, cDefect, &V, mode);
   int iLoop;
   for (iLoop = 0; iLoop < 100; ++iLoop)  {

      //mixer.addRhoIn (V);
      mixer.addRhoIn (density);
      //distributeH (V, defects, const_cast<SxMatrix<Double>&>(cDefect));
      if (mode == KeepField)  {
         newV = solvePoisson (density, eps, dV, dl) + V(0);
      } else {
         newV = solvePoisson (density, eps, V(V.getSize ()-1)-V(0), dl) + dV;
      }
      density = adjustPotential (defects, cDefect, &newV, mode);
      //cout << density << endl;
      //cout << newV << endl;
      //mixer.addRhoOut(newV);
      mixer.addRhoOut(density);
      //V = toVector (mixer.getMixedRho ()(0));
      density = toVector (mixer.getMixedRho ().getRef<SxRho> ()(0));
      //if (mode == KeepInitial) V += dV - V(0);
      //residue = mixer.getNormR () / sqrt (V.getSize ());
      //cout << "R(rho)=" << mixer.getNormR ()/ sqrt(density.getSize ()) << endl;
      residue = sqrt((V - newV).normSqr () /double(V.getSize ()));
      V = newV;
      
      //cout << "R(V) = " << residue << endl;
      //density = adjustPotential (defects, cDefect, &V, mode);
      
      if (residue < maxRes && iLoop > 1) break;
   }
   SX_CHECK (iLoop < 100);
}

void output (const SxDefectConc &defects,
             double l,
             const SxVector<Double> &V,
             const SxMatrix<Double> cDefect,
             const SxString &file,
             const SxString &header = SxString ())
{
   SxBinIO io;
   SxVector<Double> x(V.getSize ());
   double dl = l/double(V.getSize ());
   for (int i = 0; i < V.getSize (); ++i)  {
      x(i) = i * dl;
   }
   try {
      io.open (file, SxBinIO::ASCII_WRITE_ONLY);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   if (header.getSize () > 0)
      fprintf (io.fp, "%s", header.ascii ());
   SxList<SxDefectConc::Defect>::ConstIterator it = defects.begin ();
   for (int id = 0; id < defects.getNDefects (); ++id, ++it)  {
      fprintf(io.fp, "# %s\n", (*it).name.ascii ());
      io.writeXYPlot (x, cDefect.row (id));
      fprintf(io.fp, "&\n");
   }
   fprintf(io.fp, "# potential\n");
   io.writeXYPlot (x, V);
   fprintf(io.fp, "&\n");
   io.close ();
}

void diffusion (SxDefectConc &defects,
                SxMatrix<Double> &cDefect,
                SxVector<Double> *Vptr,
                double l,
                double time,
                enum AdjustMode mode,
                double dV,
                double dumpTime,
                const SxString &dumpFileBase)
{
   SX_CHECK (Vptr);
   SxVector<Double> &V = *Vptr;
   int n = int(Vptr->getSize ()), id, nDefects = int(defects.getNDefects ());
   SX_CHECK (cDefect.nCols () == n, cDefect.nCols (), n);
   SX_CHECK (cDefect.nRows () == defects.getNDefects (),
             cDefect.nRows (), defects.getNDefects ());
      
   SxList<SxDefectConc::Defect>::ConstIterator it;
   double kT = defects.getkT ();

   SxVector<Double> density(n);
   density.set (0.);

   // TODO: get from input file
   SxArray<bool> active(defects.getNDefects ());
   active.set (true);
   
   // --- find mobile defects
   SxArray<int> mobileDefects;
   SxVector<Double> diffusionConstants;
   SxVector<Double> driftConstants;

   {
      SxList<int> mobileDefectList;
      SxList<double> diffusionConstantList;
      SxList<double> driftList;
      for (id = 0, it = defects.begin (); id < nDefects; ++it, ++id)  {
         if ((*it).diffusionPrefactor > 0.)  {
            mobileDefectList << id;
            diffusionConstantList << (*it).getDiffusionConstant (kT);
            driftList << ((*it).diffusionCharge / kT);
            cout << (*it).name << ": D=" << (*it).getDiffusionConstant (kT)
                 << endl;
         }
      }
      SX_CHECK (mobileDefectList.getSize () > 0);
      mobileDefects = mobileDefectList;
      diffusionConstants = diffusionConstantList;
      driftConstants = driftList;
   }

   // --- set time step
   const double maxEffectiveD = 0.1;
   double dl = l/n;
   double dt = maxEffectiveD / diffusionConstants.maxval () * sqr(dl);

   diffusionConstants *= dt / sqr(dl);

   // --- time evolution
   int iTime, nTime = (int) ceil(time / dt);
   int dumpAt = (dumpTime > 0.) ? (int)ceil(dumpTime / dt) : 0;
   int iDump = 0;
   double last, current, next, D, fieldDrift, dVnext, dVlast;
   SxVector<Double>::Iterator cIt, nextIt;
   cout << dumpTime << endl;
   cout << "dt = "    << dt << endl;
   cout << "nTime = " << nTime << endl;
   cout << "nTimeDump = " << dumpAt << endl;
   for (iTime = 0; iTime < nTime; ++iTime)  {
      //cout << "iTime=" << iTime << endl;
      // --- compute potential and chemical equilibrium
      if (mode == KeepInitial)
         computePinnedPotential(Vptr, dV, dl, defects, cDefect);
      else
         computePotential (Vptr, &density, dV, dl,
                           defects, cDefect, 1e-5*kT, mode);

      // --- write output
      if (dumpAt && iTime % dumpAt == 0)  {
         SxString filename = dumpFileBase + iDump + ".dat";
         cout << "Save state at t=" << (iTime * dt) << " to "
              << filename << endl;
         output (defects, l, *Vptr, cDefect, filename,
                 SxString("# time = ") + (iTime * dt) + "\n");
         iDump++;
      }
      /*
      // --- chemical equilibrium
      for (int ir = 0; ir < V.getSize (); ++ir)  {
         cDefect.colRef (ir) <<=
            defects.equilibrateDefects (V(ir), active, cDefect.colRef(ir),
                                        true);
      }
      */

      // --- reactions
      if (defects.reactions.getSize () > 0)  {
         for (int ir = 0; ir < V.getSize (); ++ir)  {
            bool ok = defects.reactions.setup (cDefect.colRef(ir), defects,
                                               V(ir), kT);
            // skip evolution if any conc is negative
            if (ok)  {
               defects.reactions.evolve (dt);
               cDefect.colRef (ir) <<= defects.reactions.conc;
            } else {
               cout << "Skipping ir=" << ir << endl;
            }
         }
      }
      
      // --- diffusion step
      for (int idm = 0; idm < mobileDefects.getSize (); ++idm)  {
         id = mobileDefects(idm);
         D = diffusionConstants(idm);
         fieldDrift = 0.5 * D * driftConstants(idm);
         
         // lower boundary
         last = cDefect(id, 0);
         dVnext = V(1) - V(0);
         cDefect(id,0) += D * (cDefect(id, 1) - cDefect(id, 0))
                       + fieldDrift * dVnext * (cDefect(id,1) + cDefect(id,0));
         
         // start at il=1
         cIt = cDefect.begin (); cIt += id + nDefects;
         current = *cIt;
         nextIt = cIt;
         for (int il = 1; il < n - 1; ++il)  {
            nextIt += nDefects;
            next = *nextIt;
            dVlast = dVnext;
            dVnext = V(il+1) - V(il);
            *cIt += D * (last + next - 2. * current)
                 + fieldDrift * (  dVnext * (next + current)
                                 - dVlast * (last + current));
            //cout << "drift=" << (fieldDrift/current * (  dVnext * (next + current)
            //                     - dVlast * (last + current))) << endl;
                                  
            last = current;
            current = next;
            cIt = nextIt;
         }
         // upper boundary
         dVlast = dVnext;
         *cIt += D * (last - current)
               - fieldDrift * (dVlast * (last + current));

         // --- out-diffusion
         //cDefect(id,0) = 1e10;

         // not yet implemented
         
      }
      
   }
}

void chargeVsPot (const SxSymbolTable *cmd,
                  SxDefectConc &defects)
{
   int is, nSpecies = defects.getNSpecies ();
   SxVector<Double> concS(nSpecies);
   concS.set (0.);
   SxString element;
   double from = 0., to = 0.;
   try {
      const SxSymbolTable *speciesGroup;
      for (speciesGroup = cmd->getGroup ("species");
           speciesGroup; 
           speciesGroup = speciesGroup->nextSibling ("species"))
      {
         element = speciesGroup->get("element") ->toString ();
         is = defects.findSpecies (element);
         if (is < 0)  {
            cout << "Unknown element " << element 
                 << " in equilibrate." << endl;
            SX_QUIT;
         }
         if (speciesGroup->contains("concentration"))  {
            concS(is) = speciesGroup->get("concentration")->toReal ();
         } else {
            defects.mu(is) = speciesGroup->get("mu")->toReal ();
         }
         SxList<double> range = cmd->get("range")->toList ();
         if (range.getSize () != 2)  {
            cout << "Invalid range" << endl;
            SX_QUIT;
         }
         from = range(0);
         to   = range(1);
      }
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   SxStack<int> ignoreList;
   for (is = 0; is < nSpecies; ++is)  {
      if (concS(is) == 0.)  {
         ignoreList << is;
      }
   }

   int nDefects = defects.getNDefects ();
   SxArray<int> ignoreSpecies(ignoreList);
   SxArray<bool> active (nDefects);
   active.set (true);
   double kT = defects.getkT ();
   FILE *fp = fopen ("nVsPot.dat", "w");
   if (!fp)  {
      cout << "Unable to open nVsPot.dat" << endl;
      SX_QUIT;
   }
   for (double eFermi = from; eFermi < to; eFermi += 0.1 * kT)  {
      defects.mu = defects.chemicalEquilibrium (eFermi, active, concS, 
                                                ignoreSpecies);
      fprintf (fp, "%.8f\t%15.8g\n", eFermi,
               defects.getChargeDensity (eFermi));
      // defects.printConcentrations (eFermi);
   }
   fclose (fp);
}

class SxXYPlot {
   public:
      SxXYPlot (const SxString &fileName);

      SxVector<Double> xVals, yVals;

      double getY (double x) const;
};

SxXYPlot::SxXYPlot (const SxString &fileName)
{
   FILE *fp = fopen (fileName.ascii (), "r");
   if (!fp)  {
      cout << "Could not open " << fileName << " for reading." << endl;
      SX_QUIT;
   }
   SxStack<double> xIn, yIn;
   int l = 0;
   while (!feof(fp))  {
      l++;
      char c[2];
      if (fscanf (fp, " %[#]", c) == 1)  {
         cout << "Suppressing line " << l << endl;
      } else {
         double xVal, yVal;
         if (fscanf (fp, " %lf %lf", &xVal, &yVal) == 2)  {
            xIn << xVal;
            yIn << yVal;
         } else {
            break;
         }
      }
      // --- read until end of line
      while (!feof(fp))  {
         if (fgetc (fp) == '\n') break;
      }
   }
   cout << "Read " << l << " lines." << endl;
   fclose (fp);
   xVals = xIn;
   yVals = yIn;
   SxArray<ssize_t> idx = xVals.getSortIdx ();
   xVals.sortByIdx (idx);
   yVals.sortByIdx (idx);
}

double SxXYPlot::getY (double x) const
{
   int n = int(xVals.getSize ());
   SX_CHECK (n > 1, n);
   int l = 0, u = n - 1;
   SX_CHECK (xVals(0) <= x, xVals(0), x);
   SX_CHECK (xVals(u) >= x, xVals(u), x);
   while (u - l > 1)  {
      int i = (l + u) / 2;
      if (fabs(xVals(i) - x) < 1e-10 * x) return yVals(i);
      if (xVals(i) > x)  {
         u = i;
      } else {
         l = i;
      }
   }
   double t = (x - xVals(l)) / (xVals(u) - xVals(l));
   return t * yVals(u) + (1. - t) * yVals(l);
}

enum PlotType { Charge, Defect, Species };
void pinningCurve (const SxSymbolTable *cmd,
                   SxDefectConc &defects)
{
   int is, nSpecies = defects.getNSpecies ();
   SxVector<Double> concS(nSpecies);
   concS.set (0.);
   SxString element;
   double pinFermi = 0., startFermi = 0.;
   bool startValue = false;
   SxStack<PlotType> plotType;
   SxStack<int> idx;
   SxPtr<SxXYPlot> qVPtr;
   try {
      const SxSymbolTable *speciesGroup;
      for (speciesGroup = cmd->getGroup ("species");
           speciesGroup; 
           speciesGroup = speciesGroup->nextSibling ("species"))
      {
         element = speciesGroup->get("element") ->toString ();
         is = defects.findSpecies (element);
         if (is < 0)  {
            cout << "Unknown element " << element 
                 << " in equilibrate." << endl;
            SX_QUIT;
         }
         if (speciesGroup->contains("concentration"))  {
            concS(is) = speciesGroup->get("concentration")->toReal ();
         } else {
            defects.mu(is) = speciesGroup->get("mu")->toReal ();
         }
      }
      if (cmd->contains("alignFile"))  {
         qVPtr = SxPtr<SxXYPlot>::create (cmd->get("alignFile")->toString ());
      } 
      pinFermi = cmd->get("pinFermi")->toReal ();
      if ((startValue = cmd->contains ("startFermi")))  {
         startFermi = cmd->get("startFermi")->toReal ();
      }
      // --- read in plot data
      const SxSymbolTable *plotGroup = cmd->getGroup ("plot");
      for (plotGroup = cmd->getGroup ("plot");
           plotGroup; plotGroup = plotGroup->nextSibling ("plot"))
      {
         if (plotGroup->contains ("charge"))  {
            plotType << Charge;
            idx << -1;
         } else if (plotGroup->contains("element"))  {
            plotType << Species;
            idx << defects.findSpecies (plotGroup->get("element")->toString ());
         } else if (plotGroup->contains("name"))  {
            plotType << Defect;
            idx << defects.findDefect (plotGroup->get("name")->toString ());
         } else {
            cout << "Illegal plot type" << endl;
            SX_EXIT;
         }
      }
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   
   SxArray<PlotType> toPlot(plotType);
   SxArray<int> plotIdx (idx);

   SxStack<int> ignoreList;
   for (is = 0; is < nSpecies; ++is)  {
      if (concS(is) == 0.)  {
         ignoreList << is;
      }
   }

   int nDefects = defects.getNDefects ();
   SxArray<int> ignoreSpecies(ignoreList);
   SxArray<bool> active (nDefects);
   active.set (true);
   double kT = defects.getkT ();
   FILE *fp = fopen ("pinFermi.dat", "w");
   if (!fp)  {
      cout << "Unable to open pinFermi.dat" << endl;
      SX_QUIT;
   }
   double eFermi, dV;
   double prefactor = 1. / defects.dielecConstant; 
   double field;
   double d, x, debyeLength = 0.;
   double n0;
   if (startValue) {
      eFermi = startFermi;
      dV = 0.01 * (pinFermi > startFermi ? kT : -kT);
      // equilibrate
      defects.mu = defects.chemicalEquilibrium (eFermi, active, concS, 
                                                ignoreSpecies);
      cout << SX_SEPARATOR;
      // get charge density
      double rho = defects.getChargeDensity(eFermi);
      // set up initial values
      d = sqrt (-dV / (prefactor * rho));
      x = 0.;
      field = 0.;
   } else {
      eFermi = defects.equilibrate (cmd);
      cout << SX_SEPARATOR;
      dV = 0.01 * (pinFermi > eFermi ? kT : -kT);
      // starting value for d: 0.1 * Debye length
      double scale = 0.01;
      do {
         cout << "dn/dV for delta V=" << (scale * kT) << endl;
         defects.mu = defects.chemicalEquilibrium (eFermi+scale*kT,
                                                   active, concS, 
                                                   ignoreSpecies);
         n0  = defects.getChargeDensity (eFermi+scale*kT);
         defects.mu = defects.chemicalEquilibrium (eFermi-scale*kT,
                                                   active, concS, 
                                                   ignoreSpecies);
         n0 -= defects.getChargeDensity (eFermi-scale*kT);
         n0 /= - 2. * scale;

         if (scale < 1.) 
            scale *=2.; 
         else 
            scale += 1.;
      } while (n0 == 0.);
      cout << "n0=" << n0 << endl;
      debyeLength = sqrt (defects.dielecConstant * kT/ n0);
      d = debyeLength;
      cout << "Debye length = " << debyeLength << endl;
      // adjust Fermi energy
      defects.mu = defects.chemicalEquilibrium (eFermi,
                                                active, concS, 
                                                ignoreSpecies);
      double dEFermi = defects.getChargeDensity (eFermi) / n0 * kT;
      eFermi += dEFermi;
      cout << "eFermi  = " << eFermi << endl;
      cout << "dEFermi = " << dEFermi << endl;
      x = -100. * debyeLength;
   }

   // --- sample pinning curve
   SxVector<Double> cS(defects.getNSpecies ()),
                    c(defects.getNDefects ());
   double eFermi0 = eFermi;
   double Q = 0., Qt = 0.;
   for (; (eFermi - pinFermi) / dV <= 0.; )  {
      defects.mu = defects.chemicalEquilibrium (eFermi, active, concS, 
                                                ignoreSpecies);
      // --- compute concentrations
      SxList<SxDefectConc::Defect>::ConstIterator it;
      int id;
      cS.set (0.);
      for (it = defects.begin (), id = 0; it != defects.end (); ++it, ++id)  {
         c(id) = (*it).getConcentration (kT, defects.mu, eFermi).sum ()
               * defects.siteConc;
         cS += (*it).sumFormula * c(id);
      }
      double rho = defects.getChargeDensity(eFermi);

      // --- plot output
      fprintf(fp, "%.8g\t%.8f", x, eFermi);
      for (int io = 0; io < toPlot.getSize (); ++io)  {
         if (toPlot(io) == Charge)  {
            fprintf (fp, "\t%12.8g", rho);
         } else if (toPlot(io) == Species)  {
            fprintf (fp, "\t%12.8g", cS(plotIdx(io)));
         } else if (toPlot(io) == Defect)  {
            fprintf (fp, "\t%12.8g", c(plotIdx(io)));
         }
      }
      fprintf (fp, "\n");

      if (qVPtr)  {
         double QtOld = Qt;
         Qt = Q + qVPtr->getY (eFermi);
         cout << x << '\t' << Q << '\t' << Qt << endl;
         if (Qt * QtOld <= 0. && QtOld != 0.)  {
            double t = Qt / (Qt - QtOld);
            fprintf (fp, "# match: d=%8.4g V=%f Q=%8.4g\n",
                     x - d*t, eFermi - t*dV, Q - t * d * rho);
            fclose (fp);
            return;
         }
      }

      // --- integrate
      if (x < 0.)  {
         SX_CHECK (!startValue);
         x += d;
         double V = dV * exp (x / debyeLength);
         eFermi = eFermi0+V;
         field = V / debyeLength;
      } else {
         field -= prefactor * d * rho;
         d = dV / field;
         x+= d;
         eFermi += dV;
      }
      Q += d * rho;
   
   }
   fclose (fp);
}


#include <SxCLI.h>
#include <SxParser.h>

int main (int argc, char **argv)
{
   initSPHInXMath ();
   SxCLI cli(argc, argv);
   SxString inputFile 
      = cli.option ("-i|--input", "input file", "defect description file")
        .toString ("defects.sx");
   /*
   double l = cli.option ("-l|--length", "number", "length in cm")
              .toDouble (1e-6);
   double dV = cli.option ("-V|--deltaV", "voltage", "potential drop (eV)")
              .toDouble (1);
   int    n = cli.option ("-n|--points","integer","number of points")
              .toInt (1000,10);
   double time = cli.option ("-t|--time", "number", "time in s")
              .toDouble (1);
   */
   cli.finalize ();

   SxParser parser;
   SxParser::Table table = parser.read (inputFile);
   SxDefectConc defects (&*table);

   defects.print ();
   cout << SX_SEPARATOR;

   double eFermi = defects.getEFermi (), kT = defects.getkT ();

   defects.printConcentrations (eFermi);

   int nDefects = defects.getNDefects ();
   //int is, nSpecies = defects.getNSpecies ();

   /*
   SxMatrix<Double> cDefect(nDefects, n);
   double c;
   SxList<SxDefectConc::Defect>::ConstIterator it;
   for (it = defects.begin (); it != defects.end (); ++it, ++id)  {
      c = (*it).getConcentration (kT, defects.mu, eFermi).sum ()
          * defects.siteConc;
      for (int i = 0; i < n; ++i) cDefect(id, i) = c;
   }

   SxVector<Double> V(n);
   V.set (eFermi - dV / 2.);

   diffusion (defects, cDefect, &V, l, time);

   for (id = 0; id < defects.getNDefects (); ++id)  {
      cout << cDefect.row (id) << endl;
   }
   cout << V << endl;
   */

   /*
   SxVector<Double> V(n), density(n);
   V.set (eFermi - dV / 2.);
   density.set (0.);
   
   computePotential (&V, &density, dV, l/n,
                     defects, cDefect, 1e-3*kT);
   
   SxArray<bool> active(defects.getNDefects ());
   active.set (true);
   for (int ir = 0; ir < V.getSize (); ++ir)  {
      cDefect.colRef (ir) <<= defects.equilibrateDefects (V(ir), active, 
                              cDefect.colRef(ir), true);
   }
   computePotential (&V, &density, dV, l/n,
                     defects, cDefect, 1e-3*kT);
   */

   /*
   cout << V << endl;
   cout << density << endl;
   SxVector<Double> free(n);
   for (int i = 0; i < n; ++i)
      free(i) = defects.bands(0).getChargeDensity (V(i), defects.getkT ());
      //free(i) = defects.getFreeChargeDensity (V(i));
   cout << free << endl;
   */


   /* extract potentials
      grep '^{' sxdefectconc.log | sed -e 's/[,{}]/ /g' | awk '{for (i = 1; i <= NF; i++) print $i; print "&"; }' | xmgrace -
   */

   /*
   SxDefectConc::Band band;

   band.E0 = 0.;
   band.mass = -0.23 *1.31234228e15;
   band.degeneracy = 1;

   for (double e = -1; e < 1.; e+= 0.1*kT)
      cout << e << '\t' << band.getChargeDensity (e, kT) << endl;
   */

   /*
   SxArray<bool> active(defects.getNDefects ());
   active.set (true);
   SxVector<Double> conc(defects.getNDefects ());

   int id;
   for (id = 0, it = defects.begin (); it != defects.end (); ++it, ++id)  {
      conc(id) = (*it).getConcentration (kT, defects.mu, eFermi).sum ()
               * defects.siteConc;
   }

   conc = defects.equilibrateDefects (eFermi, active, conc);

   for (id = 0, it = defects.begin (); it != defects.end (); ++it, ++id)  {
      cout << (*it).name << ": " << conc(id) << endl;
   }
   */

   /*
   SxVector<Double> conc(defects.getNDefects ());

   int id;
   SxList<SxDefectConc::Defect>::ConstIterator it;
   for (id = 0, it = defects.begin (); it != defects.end (); ++it, ++id)  {
      conc(id) = (*it).getConcentration (kT, defects.mu, eFermi).sum ()
               * defects.siteConc;
   }

   cout << "#<<<" << endl;
   defects.reactions.setup (conc, defects, eFermi, kT);
   defects.reactions.evolve (10);
   defects.reactions.setup (conc, defects, eFermi, 1.1 * kT);
   defects.reactions.evolve (10);
   conc <<= defects.reactions.conc;
   defects.reactions.setup (conc, defects, eFermi, kT);
   defects.reactions.evolve (10);
   */




   const SxSymbolTable *main, *cmd;
   try {
      main = table->getGroup("main");
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   
   SxMatrix<Double> cDefect;
   double dV = 0.;
   enum AdjustMode mode = KeepField;
   int n;
   double length = 0.;
   SxVector<Double> V, density;

   for (cmd = main->begin (); cmd; cmd = cmd->nextSibling ())  {
      if (cmd->getName () == "equilibrate")  {
         eFermi = defects.equilibrate (cmd);
      }
      if (cmd->getName () == "chargeVsPot")  {
         chargeVsPot (cmd, defects);
      }
      if (cmd->getName () == "pinningCurve")  {
         pinningCurve (cmd, defects);
      }
      if (cmd->getName () == "profile")  {
         // --- parse input
         try  {
            n = cmd->get ("nPoints")->toInt ();
            length = cmd->get ("length")->toReal ();
            if (cmd->contains("pinFermi"))  {
               mode = KeepInitial;
               dV = cmd->get("pinFermi")->toReal ();
            } else if (cmd->contains("voltage")) {
               mode = KeepField;
               dV = cmd->get("voltage")->toReal ();
            }

         } catch (SxException e)  {
            e.print ();
            SX_EXIT;
         }

         // --- setup initial concentrations
         cDefect.reformat (nDefects, n);
         double c;
         SxList<SxDefectConc::Defect>::ConstIterator it;
         int id = 0;
         if (cmd->containsGroup ("part"))  {
            try  {
               SxSymbolTable *part = cmd->getGroup ("part");
               for (; part; part=part->nextSibling("part")) {
                  int from, to;
                  SxArray<int> range = part->get("range")->toIntList ();
                  if (range.getSize () != 2)  {
                     cout << "Illegal range size" << endl;
                     SX_QUIT;
                  }
                  from = range(0)-1;
                  to = range(1)-1;
                  if (from < 0)  {
                     cout << "Illegal range: from < 1!" << endl;
                     SX_QUIT;
                  }
                  if (to > n)  {
                     cout << "Illegal range: to > nPoints" << endl;
                     SX_QUIT;
                  }
                  if (from >= to)  {
                     cout << "Illegal range: from >= to" << endl;
                     SX_QUIT;
                  }
                  cout << SX_SEPARATOR;
                  cout << "Range " << from << "..." << to << endl;
                  eFermi = defects.equilibrate (part);
                  for (id = 0, it = defects.begin (); it != defects.end ();
                       ++it, ++id)  {
                     c = (*it).getConcentration (kT, defects.mu, eFermi).sum ()
                         * defects.siteConc;
                     for (int i = from; i <= to; ++i) cDefect(id, i) = c;
                  }
               }
            } catch (SxException e)  {
               e.print ();
               SX_EXIT;
            }
         } else {
            for (it = defects.begin (); it != defects.end (); ++it, ++id)  {
               c = (*it).getConcentration (kT, defects.mu, eFermi).sum ()
                   * defects.siteConc;
               for (int i = 0; i < n; ++i) cDefect(id, i) = c;
            }
         }
         // gradient for defect 1
         //for (int i = 0; i < n; ++i) cDefect(0, i) *= 0.9 + i*0.2/n;

         // --- setup potential
         if (V.getSize () != n)  {
            V.resize (n);
            density.resize (n);
            if (mode == KeepField)  {
               V.set (eFermi - dV / 2.);
            } else {
               V.set (eFermi);
            }
         }
         density.set (0.);
         double dl = length / n;

         //bool equilibrate = (defects.mobileSpecies.getSize () > 0);

         /*
         SxVector<Double> oldV;
         for (int ipl = 0; ipl < 100; ++ipl)  {
            computePotential (&V, &density, dV, dl,
                              defects, cDefect, 1e-3*kT, mode);
            cout << V << endl;

            distributeH (V, defects, cDefect);

            if (ipl > 0)
              cout << "delta V=" << ((oldV - V).normSqr ()/n) << endl;
            if (ipl > 0 && (oldV - V).normSqr () < sqr(1e-2 * kT) * n)
               break;
            oldV.copy (V);
         }
         */
         if (mode == KeepInitial)
            computePinnedPotential(&V, dV, dl, defects, cDefect);
         else
            computePotential (&V, &density, dV, dl,
                              defects, cDefect, 1e-5*kT, mode);

         // --- output
         output (defects, length, V, cDefect, "cDefect.dat");
         /*
         for (id = 0; id < defects.getNDefects (); ++id)  {
            cout << cDefect.row (id) << endl;
         }
         cout << V << endl;
         */
      }
      
      if (cmd->getName () == "diffusion")  {
         double time;
         double dumpTime = -1.;
         SxString fileBase = "diffusion-";
         try {
            time = cmd->get("time")->toReal ();
            if (cmd->contains("dumpTime"))  {
               dumpTime = cmd->get("dumpTime")->toReal ();
            }
            if (cmd->contains ("file"))  {
               fileBase = cmd->get("file")->toString ();
            }
         } catch (SxException e)  {
            e.print ();
            SX_EXIT;
         }
         if (cDefect.getSize () == 0)  {
            cout << "Need profile {} group to initialize diffusion {}." << endl;
            SX_QUIT;
         }
         diffusion (defects, cDefect, &V, length, time, mode, dV,
                    dumpTime, fileBase);
         // --- output
         output (defects, length, V, cDefect,"cDefect.dat");
      }
         
   }
   //printTiming ();
}
#endif

