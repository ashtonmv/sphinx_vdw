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

#ifndef _SX_PARTIAL_RHO_H_
#define _SX_PARTIAL_RHO_H_
#include <SxRho.h>
#include <SxFermi.h>
#include <SxPW.h>
#include <SxPtr.h>
#include <SxExt.h>

/** \brief Partial densities

  \b SxClass = S/PHI/nX ...

  Addition calculates partial charge density resolved
  with respect to one-particle energies
  or with respect to k and i
  @ingroup Add-Ons
  @author Alexey Dick, dick@fhi-berlin.mpg.de
  @author Christoph Freysoldt freyso@fhi-berlin.mpg.de

*/
class SX_EXPORT_EXT SxPartialRho : private SxRho
{
   public:
      
      /// Pointer to waves
      SxPtr<SxPW> waves;
      /// Pointer to energies/occupations
      SxPtr<SxFermi> fermi;
      /** \brief Whether to use the occupation number
          
           If set to false, full occupation is assumed for
           all states.
        */
      bool useFocc;

      /** \brief Whether to symmetrize the density */
      bool symmetrize;

      /** \brief Constructor */
      SxPartialRho (const SxPtr<SxPW>    &wavesIn, 
                    const SxPtr<SxFermi> &fermiIn);

      /** \brief Destructor */
      ~SxPartialRho () { /* empty */}

      /// Write partial density to file
      void writeRho (const SxString &fileName) const
      {
         SxRho::writeRho (fileName);
      }

      /// Get computed rho
      const RhoR& getRho () const
      {
         return rhoR;
      }
      
      /// Compute partial density in energy interval
      void compute(double eMin, double eMax);

      /**
        The electron density is calculated from all states in the states' list 
        for all k-points in the k-point list. 
        The k-points are weighted according to the k-point weight 
        (from the |G+k> basis). The states are occupied according to the 
        Fermi distribution if useFocc is true.
        Otherwise, full occupation (1 for spin-polarized waves, or 2) is
        assumed.
        In this way, the unoccupied states can be visualized, too.
       */
      void compute(const SxList<int> &kPoints, const SxList<int> &states);
};

#endif /* _SX_PARTIAL_RHO_H_ */
