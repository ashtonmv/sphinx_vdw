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

#ifndef _SX_SPECIES_DATA_H_
#define _SX_SPECIES_DATA_H_

#include <SxArray.h>
#include <SxSymbolTable.h>
#include <SxString.h>
#include <SxParser.h>
#include <SxMap.h>
#include <SxBinIO.h>
#include <SxGeom.h>

/** \brief Data about species involved in the system

    \b SxSpeciesData = SPHInX Species Data

    This class provides data which is independent of the method being 
    used. All method-specific stuff goes into a derived class.

    \todo Throw out all structure optimization stuff from here.

    \author  Sixten Boeck, boeck@mpie.de */
class SX_EXPORT_GEOM SxSpeciesData
{
   public:
      SxSpeciesData ();
      SxSpeciesData (const SxSymbolTable *table);

      /** Init from the list of unique chemNames. */
      SxSpeciesData (const SxList<SxString> &chemNames_);

      virtual ~SxSpeciesData ();

      void readSpecies (const SxSymbolTable *table);

      /** \brief Guess included species files
        
        \param parser a parser where the table has been read from
        \param table  top-level or species level of symbol tree
        
        This routine tries to guess the included species files.
        This is done by locating the origin of the "valenceCharge"
        entry. The return value is an array of filenames (including
        partial pathes that are needed in addition to the parser search
        path).

        \example
        \code
SxParser parser;
SxParser::Table table = parser.read (inFile);
SxArray<SxString> speciesFiles;
speciesFiles = SxSpeciesData::guessIncludeFiles (parser, table);
SxAtomicStructure structure (&*table);
SxSpeciesData sData(&*table);
structure.fprint (stdout, sData.elementNames, speciesFiles);
        \endcode
        
        */
      static SxArray<SxString> guessIncludeFiles (const SxParser& parser,
                                                  SxParser::Table table);

      /** \brief Read chemical symbols from input file
          
          The chemical symbols are the minimum species information needed
          to print out structures. Usually, they are part of the potential
          group. However, in minimal input structures, they may be placed
          directly into the structure group.

          This routine finds the chemical symbols in both cases, preferring
          the potential group over the structure group.

        */
      static SxArray<SxString> getElements (const SxSymbolTable *);

      /// returns number of species
      int getNSpecies () const;

      /** \brief Find species number by chemName, elementName, or number

          This tries to identify a species by name or number. If
          the species is not found, -1 is returned.
          \example
          For GaAs, the following strings would identify arsenic 
          (second species):
          -# "As"
          -# "Arsenic"
          -# "2"
          all of this would return 1 (as we start from 0)

        */
      int find (const SxString &name) const;
      
      /** @brief Standard PSE name of the chemical element
       
          The international identifies for chemical elements are usually
          1 or 2 characters, like "H" for Hydrogen, "He" for Helium etc. */
      SxArray<SxString>   chemName;
      /** @brief User provided name for the species */
      SxArray<SxString>   elementName;
      /** @brief Nuclear (all-electron) charge */
      SxArray<double>     nuclearCharge;
      /** @brief Valence charge */
      SxArray<double>     valenceCharge;
      /** @brief Damping parameter for damped Newton relaxation

          \sa SxStructOpt::dampedNewton */
      SxArray<double>     dampingMass;
      /** @brief Reciprocal mass for damped Newton relaxation

          \sa SxStructOpt::dampedNewton */
      SxArray<double>     reciprocalMass;

      /** @brief Ionic mass */
      SxArray<double>     ionicMass;

      /** \brief Find appropriate group for SxSpeciesData
        \note This looks for the registered group names and structure {}
        */
      static const SxSymbolTable* 
      findTable (const SxSymbolTable *tableIn);

      /// Produce best possible name for a species
      SxString prettyName (int is) const;
      
   public:
      /// Auxiliary class for printing in input-file format
      class SX_EXPORT_GEOM PrintInfo {
         public:
            /// Input group name
            SxString groupName;
            /// Species include files
            SxArray<SxString> includeFiles;
            /// Setup routine
            void setup (const SxParser &parser, const SxSymbolTable *table);

            /// Empty constructor
            PrintInfo () { /* empty */ }
            /// Constructor
            PrintInfo (const SxParser &parser, const SxSymbolTable *table)
            {
               setup (parser, table); 
            }
            
            /// Check if printing info is available
            bool isValid () const { return groupName.getSize () > 0; }
            /// Print to file in input format
            void fprint (FILE *output) const;
      };
      /** \brief Read chemNames from binary file.
        */
      void read (const SxBinIO &);
};

#endif /* _SX_SPECIES_DATA_H_ */
