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

#ifndef _SX_FOCK_TIMER_H_
#define _SX_FOCK_TIMER_H_
#include <SxTimer.h>

namespace Timer {
   enum FockTimer {
      FockCompute,
      FockApply
   };
}

SX_REGISTER_TIMERS (Timer::FockTimer)
{
   using namespace Timer;
   regTimer (FockCompute, "Fock computation");
   regTimer (FockApply,   "Fock * psi");
}

#endif /* _SX_FOCK_TIMER_H_ */
