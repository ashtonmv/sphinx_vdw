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
#include <SxElemDB.h>
#include <SxList.h>
#include <SxParser.h>

SxElemDB::SxElemDB ()
{
   SxParser parser;
   try {
      SxParser::Table table = parser.read ("species/elements.sx");
      read(&*table);
   } catch (SxException e) {
      cout << "WARNING: failed to load element database" << endl;
   }
}


SxElemDB::SxElemDB (const SxSymbolTable *table)
{
   read (table);
}

SxElemDB::SxElemDB (const SxConstPtr<SxSymbolTable> &table)
{
   read (&*table);
}

SxElemDB::~SxElemDB ()
{
   // empty
}

void SxElemDB::read (const SxSymbolTable *table)
{
   SX_CHECK (table);
   SX_CHECK (name.getSize() == 0);
   try  {
      SxSymbolTable *group = table->getGroup ("chemElements", true);
      SxSymbolTable *elem  = group->getGroup ("elem", true);

      // --- count number of elements
      int i=0, nElem = 0;
      while (elem)  {
         elem = elem->nextSibling ("elem");
         ++nElem;
      }
      elem  = group->getGroup ("elem", true);

      // --- allocate memory
      chemSymbol.resize(nElem);
      name.resize(nElem);
      atomicNumber.resize(nElem);
      weight.resize(nElem);
      covalentRadius.resize(nElem);
      atomicRadius.resize(nElem);
      gyromagneticRatio.resize(nElem);
      red.resize(nElem);
      green.resize(nElem);
      blue.resize(nElem);
      polarizability.resize(nElem);
      C6.resize(nElem);
      vdwRadius.resize(nElem);

      SxList<double> color;
//    while (elem && i < nElem)  {
      while (elem)  {
         color = elem->get("CPKcolor", true)->toList();
         chemSymbol(i)     = (elem->get ("symbol", true)->toString()
                             .toUpper().ascii());
         name(i)           = (elem->get ("name", true)->toString().ascii());
         atomicNumber(i)   = (elem->get ("number", true)->toInt());
         weight(i)         = (elem->get ("weight", true)->toReal());
         covalentRadius(i) = (elem->get ("covalentRadius", true)->toReal());
         atomicRadius(i)   = (elem->get ("atomicRadius", true)->toReal());
         if (elem->contains ("gyromagneticRatio"))  {
            gyromagneticRatio(i) = 
               (elem->get ("gyromagneticRatio", true)->toReal());
         }  else  {
            gyromagneticRatio(i) = 0.;
         }
         red(i)            = (color(0)/255.);
         green(i)          = (color(1)/255.);
         blue(i)           = (color(2)/255.);
         
         polarizability(i) = (elem->get ("polarizability", true)->toReal());
         C6(i) = (elem->get ("C6", true)->toReal());
         vdwRadius(i) = (elem->get ("vdwRadius", true)->toReal());

         elem = elem->nextSibling ("elem");
         ++i;
      }
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}




SxString SxElemDB::getChemSymbol (double atomicMass)  const
{
   SX_CHECK (weight.getSize() == chemSymbol.getSize(),
             weight.getSize(),   chemSymbol.getSize());

   ssize_t i, nElem = weight.getSize();
   for (i=0; i < nElem; ++i)
      if ( fabs(weight(i) - atomicMass) < 0.1 )  return chemSymbol(i);
   return "??";  // nothing found
}


   

int SxElemDB::getAtomicNumber (int idx) const
{
   if (idx < atomicNumber.getSize())  return atomicNumber(idx);
   else                               return 0; // not found
}


int SxElemDB::getAtomicNumber (const SxString &chemName) const
{
   SX_CHECK (chemSymbol.getSize() == atomicNumber.getSize(),
             chemSymbol.getSize(),   atomicNumber.getSize());

   SxString inName = chemName.toUpper();

   ssize_t i, nElem = chemSymbol.getSize();
   for (i=0; i < nElem; ++i)
      if (chemSymbol(i).toUpper() == inName)  return atomicNumber(i);
   return 0;  // not found
}

int SxElemDB::getRow (const SxString &chemName) const
{
   int Z = getAtomicNumber (chemName);
   if (Z == 0) return 0;
   if (Z <= 2) return 1;
   if (Z <= 10) return 2;
   if (Z <= 18) return 3;
   if (Z <= 36) return 4;
   if (Z <= 54) return 5;
   if (Z <= 86) return 6;
   if (Z <= 118) return 7;
   return -1;
}

double SxElemDB::getAtomicWeight (const SxString &chemName) const 
{
	return weight(getIdx(getAtomicNumber(chemName)));
}

double SxElemDB::getAtomicWeight (int idx) const
{
   if (idx < weight.getSize())  return weight(idx);
   else                         return 0;  // not found
}

int SxElemDB::getIdx (int atomNumber) const
{
   int i, nElem = int(atomicNumber.getSize());
   for (i=0; i < nElem; ++i)
      if (atomicNumber(i) == atomNumber)  return i;

   sxprintf ("SxElemDB::getIdx: Can't find element %d\n", atomNumber);
   return -1;  // not found
}

int SxElemDB::getIdx (const SxString &chemName_) const
{
   SX_CHECK (chemSymbol.getSize() == atomicNumber.getSize(),
             chemSymbol.getSize(),   atomicNumber.getSize());

   SxString inName = chemName_.toUpper();

   int i, nElem = int(chemSymbol.getSize());
   for (i=0; i < nElem; ++i)
      if (chemSymbol(i).toUpper() == inName)  return i;
   return -1;  // not found
}

double SxElemDB::getRed (int i) const
{
   return red(i);
}


double SxElemDB::getGreen (int i) const
{
   return green(i);
}


double SxElemDB::getBlue (int i) const
{
   return blue(i);
}

SxVector3<Double> SxElemDB::getRGB (int i) const
{
   return SxVector3<Double> (red(i), green(i), blue(i));
}

double SxElemDB::getCovalentRadius (const SxString &chemName_) const
{
   SX_CHECK (chemSymbol.getSize() == atomicNumber.getSize(),
             chemSymbol.getSize(),   atomicNumber.getSize());

   SxString inName = chemName_.toUpper();

   ssize_t i, nElem = chemSymbol.getSize();
   for (i=0; i < nElem; ++i)
      if (chemSymbol(i).toUpper() == inName)  return covalentRadius(i);
   return 0;  // not found
}

double SxElemDB::getCovalentRadius (int i) const
{
   return covalentRadius(i);
}

double SxElemDB::getAtomicRadius (int i) const
{
   return atomicRadius(i);
}

SxString SxElemDB::getChemSymbol (int i) const
{
   return chemSymbol(i);
}


SxString SxElemDB::getName (int i) const
{
   return name(i);
}

double SxElemDB::getGyrRatio (const SxString &chemName_) const
{
   SX_CHECK (chemSymbol.getSize() == gyromagneticRatio.getSize(),
             chemSymbol.getSize(),   gyromagneticRatio.getSize());

   SxString inName = chemName_.toUpper();

   ssize_t i, nElem = chemSymbol.getSize();
   for (i=0; i < nElem; ++i)
      if (chemSymbol(i).toUpper() == inName)  return gyromagneticRatio(i);
   return 0;  // not found
}

void SxElemDB::set (int idx_, int atomNumber_, const SxString &chemSymbol_,
                    const SxString &name_, double weight_, double rCov_,
                    double rAtom_, double gyrAt_, double red_, double green_, 
                    double blue_)
{
   chemSymbol(idx_)     = chemSymbol_;
   name(idx_)           = name_;
   atomicNumber(idx_)   = atomNumber_;
   weight(idx_)         = weight_;
   covalentRadius(idx_) = rCov_;
   atomicRadius(idx_)   = rAtom_;
   gyromagneticRatio(idx_) = gyrAt_;
   red(idx_)            = red_;
   green(idx_)          = green_;
   blue(idx_)           = blue_;
}


void SxElemDB::append (int atomNumber_, const SxString &chemSymbol_,
                       const SxString &name_, double weight_, double rCov_,
                       double rAtom_, double gyrAt_, double red_, double green_, 
                       double blue_)
{
   ssize_t i, nElem;
   i = nElem = chemSymbol.getSize();
   ++nElem;

   chemSymbol.resize (nElem, true);     chemSymbol(i)     = chemSymbol_;
   name.resize (nElem, true);           name(i)           = name_;
   atomicNumber.resize (nElem, true);   atomicNumber(i)   = atomNumber_;
   weight.resize (nElem, true);         weight(i)         = weight_;
   covalentRadius.resize (nElem, true); covalentRadius(i) = rCov_;
   atomicRadius.resize (nElem, true);   atomicRadius(i)   = rAtom_;
   gyromagneticRatio.resize (nElem, true);   gyromagneticRatio(i) = gyrAt_;
   red.resize (nElem, true);            red(i)            = red_;
   green.resize (nElem, true);          green(i)          = green_;
   blue.resize (nElem, true);           blue(i)           = blue_;
}


void SxElemDB::remove (int idx)
{
   SxList<SxString> strList;
   SxList<int>      intList;
   SxList<double>   dblList;

   ssize_t i, nElem = chemSymbol.getSize();
   for (i=0; i < nElem; ++i)  if (i != idx)  strList << chemSymbol(i);
   chemSymbol = strList;
   strList.removeAll ();

   for (i=0; i < nElem; ++i)  if (i != idx)  strList << name(i);
   name = strList;
   strList.removeAll ();

   for (i=0; i < nElem; ++i)  if (i != idx)  intList << atomicNumber(i);
   atomicNumber = intList;
   intList.removeAll ();
   
   for (i=0; i < nElem; ++i)  if (i != idx)  dblList << weight(i);
   weight = dblList;
   dblList.removeAll ();
   
   for (i=0; i < nElem; ++i)  if (i != idx)  dblList << covalentRadius(i);
   covalentRadius = dblList;
   dblList.removeAll ();
   
   for (i=0; i < nElem; ++i)  if (i != idx)  dblList << atomicRadius(i);
   atomicRadius = dblList;
   dblList.removeAll ();

   for (i=0; i < nElem; ++i)  if (i != idx)  dblList << gyromagneticRatio(i);
   gyromagneticRatio = dblList;
   dblList.removeAll ();
   
   for (i=0; i < nElem; ++i)  if (i != idx)  dblList << red(i);
   red = dblList;
   dblList.removeAll ();
   
   for (i=0; i < nElem; ++i)  if (i != idx)  dblList << green(i);
   green = dblList;
   dblList.removeAll ();
   
   for (i=0; i < nElem; ++i)  if (i != idx)  dblList << blue(i);
   blue = dblList;
   dblList.removeAll ();
}




void SxElemDB::print () const
{
   sxprintf ("\nElements:\n");
   for (int i=0; i < chemSymbol.getSize(); i++)  {
      sxprintf ("Symbol:         %s (%d)\n", 
              chemSymbol(i).ascii(), atomicNumber(i));
      sxprintf ("Name:           %s\n", name(i).ascii());
      sxprintf ("weight:         %g\n", weight(i));
      sxprintf ("cov.   radius:  %g\n", covalentRadius(i));
      sxprintf ("atomic radius:  %g\n", atomicRadius(i));
      sxprintf ("gyromagnetic Ratio:  %g\n", gyromagneticRatio(i));
      sxprintf ("color:          (%g,%g,%g)\n", red(i), green(i), blue(i));
      sxprintf ("\n");
   }
   
}

const SxElemDB& SxElemDB::getElemDB ()
{
   static SxPtr<SxElemDB> dbPtr;
   if (!dbPtr) dbPtr = SxPtr<SxElemDB>::create ();
   return *dbPtr;
}
