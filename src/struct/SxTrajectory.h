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

#ifndef _SX_TRAJECTORY_H_
#define _SX_TRAJECTORY_H_

#include <SxMatrix.h>
#include <SxVector.h>
#include <SxBinIO.h>
#include <SxStruct.h>


/** \brief sxtrajectory 

    \b SxHessianOps = class for i/o management of trajectories
       input-data: forces, structures, velocities, energies etc. 
                   from md and struct-opt
       output-files: binary trajectory file 
                     several ascii trajectory files       
    \ingroup   group_structure
    \author    Lars Ismer, ismer@fhi-berlin.mpg.de */

class SX_EXPORT_STRUCT SxTrajectory
{
   public:

      //---------------------------------------------------------------------
      /**\brief Constructors and Destructors */
      //---------------------------------------------------------------------
      SxTrajectory ();
      ~SxTrajectory ();
};

#endif 
