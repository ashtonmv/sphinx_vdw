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

#ifndef _SX_FORST_H_
#define _SX_FORST_H_

#include <SxClassic.h>
#include <SxPotential.h>
#include <SxAtomicStructure.h>
#include <SxSymbolTable.h>
#include <SxForstData.h>

/** \brief ...

    \b SxForst = S/PHI/nX ...

    ....

    \author Timo Jacob, jacob@fhi-berlin.mpg.de */
class SX_EXPORT_CLASSIC SxForst : public SxPotential
{
   public:

      /** \brief ...

          ... 
       */
      SxForst ();

      /** \brief ...

          ... 
       */
      SxForst (const SxAtomicStructure &, const SxSymbolTable *);

      /** \brief Destructor
       */
      virtual ~SxForst ();

      /** \brief ...

          ...
       */
      virtual bool isRegistered (const SxSymbolTable *) const;

      /** \brief ...

          ...
          */
      virtual void execute (const SxSymbolTable *, bool calc=true);

      /** \brief ...

          ...
          */
      virtual SxAtomicStructure getForces (const SxAtomicStructure &,
                                           const SxSymbolTable *cmd=NULL);

      
      /** \brief ...

          ...
          */
      virtual PrecEnergy getEnergy () const;

      virtual SxSpeciesData getSpeciesData () const;

   protected:

      SxForstData speciesData;


};

#endif /* _SX_FORST_H_ */
