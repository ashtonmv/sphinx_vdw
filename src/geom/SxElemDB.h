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


#ifndef SX_ELEMENTS_H
#define SX_ELEMENTS_H

#include <SxSymbolTable.h>
#include <SxGeom.h>
#include <SxPtr.h>
#include <SxString.h>
#include <SxVector3.h>
#include <SxArray.h>

/** This class contains the most relevant data taken from the 
    periodic system of elements (PSE).
    The data themselve are read in from a file called elements.sx
    which is being shipped with the SPHInX package.
    This file must reside in the search path.

    Usually this class should be instantiated once and it is accessable
    globally. 

    \par Initialization
    The default constructor reads the database file
    share/sphinx/species/elements.sx which is distributed with the
    SPHInX package.
    Alternative files can be selected by using the constructor from a
    symbol table.
\code
   #include <SxElemDB.h>
   #include <SxParser.h>
   ...
   SxParser parser;
   SxParser::Table elemTable = parser.read (myElementsFile);
   SxElemDB elemDB (elemTable);

\endcode
    
    \par Examples
\code
 int covRad = elemDB->getCovalentRadius (8);
\endcode
    or
\code
 int z = elemDB->getAtomicNumber("O");
 int covRad = elemDB->getCovalentRadius (z);
\endcode
    @author Sixten Boeck
  */
class SX_EXPORT_GEOM SxElemDB
{
   public:

      SxElemDB ();
      SxElemDB (const SxSymbolTable *);
      SxElemDB (const SxConstPtr<SxSymbolTable> &);
      ~SxElemDB ();

   protected:
      void read (const SxSymbolTable *);

   public:
      /** Returns the atomic number from the given chemical symbol.
          This routine is case-insensitive. */
      int getAtomicNumber      (const SxString &chemName) const;
      int getAtomicNumber      (int idx) const;

      /// Return the row for the given chemical symbol
      int getRow (const SxString &chemName) const;
      /** Return the index of the entry in the database from an atomic
          number. If the elements.sx file is sorted corectly these
          numbers are equal. */
      int getIdx               (int atomNumber) const;
      int getIdx               (const SxString &chemName_) const;
      /** Returns the red component of the CPK color of the i-th element.
          The color values range between 0 and 1 */
      double getRed            (int i) const;
      /** Returns the green component of the CPK color of the i-th element 
          The color values range between 0 and 1 */
      double getGreen          (int i) const;
      /** Returns the blue component of the CPK color of the i-th element
          The color values range between 0 and 1 */
      double getBlue           (int i) const;
		
      /** Returns the CPK color of all 3 color components at once between
          0...1 */
      SxVector3<Double> getRGB (int i) const;
      /** Returns the covalent radius of the i-th element */
      double getCovalentRadius (const SxString &chemName_) const;
      double getCovalentRadius (int i) const;
      /** Returns the atomic radius of the i-th element */
      double getAtomicRadius (int i) const;
      /** Returns the chemical Symbol for a given input mass */
      SxString getChemSymbol (double atomicMass) const;
      /** Returns the chemical symbol of the i-th element */
      SxString getChemSymbol   (int i) const;
      /** Returns the (engl.) name  of the i-th element */
      SxString getName         (int i) const;

      /** Returns the polarizability of the i-th element */
      double getPolarizability (const SxString &chemName_) const;
      double getPolarizability (int i) const;

      /** Returns the dispersion coefficient (C6) of the i-th element */
      double getC6 (const SxString &chemName_) const;
      double getC6 (int i) const;

      /** Returns the vdW radius of the i-th element */
      double getVdwRadius (const SxString &chemName_) const;
      double getVdwRadius (int i) const;

		 /** Returns the atomic weight for a given chemical symbol
          (does not work for isotopes) 
          \author Lars Ismer */
		double getAtomicWeight (const SxString &chemName) const; 
		double getAtomicWeight (int idx) const; 
      double getGyrRatio (const SxString &chemName_) const;



      int getSize () const { return int(name.getSize()); }

      void set (int idx_, int atomNumber_, const SxString &chemSymbol_,
                const SxString &name_, double weight_, double rCov_,
                double rAtom_, double gyrAt_, double red_, double green_, 
                double blue_);

      void append (int atomNumber_, const SxString &chemSymbol_,
                const SxString &name_, double weight_, double rCov_,
                double rAtom_, double gyrAt_, double red_, double green_,
                double blue_);

      void remove (int idx);
      

      void print () const;

      static const SxElemDB &getElemDB ();

   protected:
      /** The standardized chemical symbol. It is one or two characters
          long, for example "H" - Hydrogen, "Ga" - Gallium, etc. */
      SxArray<SxString> chemSymbol;
      /** The (english) name of the elements, e.g. "Hydrogen" */
      SxArray<SxString> name;
      /** The atomic number, 1 - Hydrogen, 103 - Lawrencium */
      SxArray<int>     atomicNumber;
      /** The atomic weight based upon Carbon */
      SxArray<double>  weight;
      /** The covalent radii in Bohr */
      SxArray<double>  covalentRadius;   // [Bohr]
      /** The atomic radii in Bohr */
      SxArray<double>  atomicRadius;   // [Bohr]
      /** The gyromagnetic ratio in MHz/T **/
      SxArray<double> gyromagneticRatio; // [MHz/T]
      /** The atomic radii in Bohr. Elements with a coordination number
          12 the metallic radii are taken */
      /** The red color channel for the CPK model, 0 ... 1 */
      SxArray<double>  red;
      /** The green color channel for the CPK model, 0 ... 1 */
      SxArray<double>  green;
      /** The blue color channel for the CPK model, 0 ... 1 */
      SxArray<double>  blue;
      /** The polarizability in Bohr^3 */
      SxArray<double>  polarizability;
      /** The dispersion coefficient C6 in Hartree * Bohr^6 */
      SxArray<double> C6;
      /** The vdW radius in Bohr */
      SxArray<double> vdwRadius;
};

#endif /* SX_ELEMENTS_H */
