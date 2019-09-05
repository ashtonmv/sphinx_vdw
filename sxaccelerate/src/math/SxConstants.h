// ---------------------------------------------------------------------------
//
//           The general purpose cross platform C/C++ framework
//
//                       S x A c c e l e r a t e
//
//           Home:       https://www.sxlib.de
//           License:    Apache 2
//           Authors:    see src/AUTHORS
//
// ---------------------------------------------------------------------------
#ifndef _SX_CONSTANTS_H_
#define _SX_CONSTANTS_H_

#include <SxComplex.h>

/** \brief complex \b i */
#define I             SxComplex16 (0.,1.)

/** \brief \f$\pi\f$ */
#define PI            3.14159265358979323846264338327950

/** \brief \f$2 \pi\f$ */
#define TWO_PI        6.28318530717958647692528676655900

/** \brief \f$4 \pi\f$ */
#define FOUR_PI      12.56637061435917295385057353311800

/** \brief \f$2 \pi^2\f$ */
#define TWO_PI_2     19.73920880217871723766898199975226

/** \brief \f$4 \pi^2\f$ */
#define FOUR_PI_2    39.47841760435743447533796399950452

/** \brief \f$\pi^2 \f$ */
#define PI_2          9.86960440108935861883449099987613

/** \brief \f$\sqrt{\pi}\f$ */
#define SQRT_PI       1.77245385090551602729816748334114

/** \brief \f$\sqrt{\frac{1}{4\pi}}\f$ */
#define SQRT_1_4PI    .282094791773878143474039725780386

/** \brief \f$\sqrt{2}\f$ */
#define SQRT2         1.41421356237309504880168872420970

/** \brief \f$\sqrt{3}\f$ */
#define SQRT3         1.73205080756887729352744634150587

/** \brief \f$\sqrt{5}\f$ */
#define SQRT5         2.23606797749978969640917366873127

/** \brief \f$\sqrt{12}\f$ */
#define SQRT12        3.46410161513775458705489268301174

/** \brief \f$\sqrt{15}\f$ */
#define SQRT15        3.87298334620741688517926539978240

/** \brief \f$\sqrt{30}\f$ */
#define SQRT30        5.47722557505166113456969782800802

/** \brief \f$\cos{30}\f$  */
#define COS30         0.86602540378443864676372317075294

/** \brief Hartree ->eV/meV */
#define HA2EV            27.211383
#define HA2MEV            27211.383

/** \brief eV/meV->Hartree */
#define EV2HA            1.0/27.211383
#define MEV2HA           1.0/27211.383

/** \brief Just a large value */
#define SX_HUGE       1e50

/** \brief prefactor \frac{180}{\pi} used when converting from RAD to DEG */
#define RAD2DEG       180./PI

/** \brief prefactor \frac{\pi}{180} used when converting from DEG to RAD */
#define DEG2RAD       PI/180.

/** \brief Boltzmann constant in J/K */
#define KB            1.3806504e-23

/** \brief Planck constant in SI units*/
#define HPLANCK       6.62606896e-34

/** velocity of light in vacuum in m/s */
#define CVEL          299792458

/** S/PHI/nX time units (atom units) to seconds*/
// sqrt(mass_H bohrradius^2/hartree)
#define AU2S          1.0327499e-15

/** \brief Joule -> Hartree */
#define JL2HA         2.2937128e+17

/** \brief Hartree -> kcal/mol */

#define HA2KCALPM     627.08975

/** \brief 1/cm -> kcal/mol */

#define CM2KCALPM     0.0028572312

/**\ brief Joule -> kcal/mol */

#define JL2KCALPM     1.4383638e+20

/**\brief  1/(2pi * atomic time unit) to cm^-1*/
 
#define AU2CM         5140.4871

/**\brief Bohrradius in m*/
#define BOHRRADIUS    5.2917721e-11

/**\brief angstrom in bohr*/
#define A2B           1.8897261

/**\brief atomic mass unit -> Kg (=AU2S^2/JL2HA/BOHRRADIUS^2) */ 
#define U2KG          1.66053886e-27

/**\brief TWO_PI over atomic time unit -> meV (=1000*HPLANCK*JL2HA*HA2EV/TWO_PI/AU2S) */
#define AU2MEV        637.33912

#define SPIN_UP       0
#define SPIN_DOWN     1

#endif // _SX_CONSTANTS_H_
