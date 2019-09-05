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

#ifndef _SX_MURN_H_
#define _SX_MURN_H_

#include <SxString.h>
#include <SxVector.h>
#include <SxList.h>
#include <SxExt.h>

/** \brief Murnaghan fit

    \b SxMurn = S/PHI/nX Murnaghan Fit
    
    The Murnaghan state equation resolved to energy reads
    \f[E = E_0 + \frac{B_0V}{B_0'(B_0' - 1)}\left[ 
    B_0'\left(1 - \frac{V_0}{V}\right) + \left(\frac{V_0}{V}\right)^{B_0'} - 1
    \right] \f]
 
    The add-on performs a fit of the parameters \f$B_0\f$, \f$B_0'\f$, 
    \f$V_0\f$, and \f$E_0\f$ to a number of points (volume/energy).
    
    \par Fitting procedure
    In this add-on, the bulk modulus \f$B_0\f$ and the minimum energy \f$E_0\f$
    are optimized analytically, while the derivative \f$B_0'\f$ and the minimum
    volume are obtained by two (very primitive) single line optimizations.
    \n
    The quadratic deviation is defined as
    \n
    \f[\Delta = \sum_n (E(V_n) - E_n)^2 \f]
    \n
    where \f$V_n\f$ / \f$E_n\f$ are the given points, and E(V) is the Murnaghan 
    function.
    \n
    Putting all the complicated stuff into a function f, the thing looks easier
    \f[E(V) = E_0 + B_0\cdot f_{V_0,B_0'} (V) \f]
    \n
    If \f$E_0\f$ is the optimal parameter, the partial derivative of 
    \f$\Delta\f$ with respect to \f$E_0\f$ must vanish. Thus
    \n
    \f[0 = \frac{\partial \Delta}{\partial E_0} 
    =  2 \sum_n \left(E(V_n) - E_n\right) \cdot \sum_n 
    \underbrace{\frac{\partial E(V_n)}{\partial E_0}}_1
    = 2[\underbrace{\sum_n (E_0}_{n\cdot E_0} + B_0 \cdot f(V_n) - E_n)] \cdot n 
    \f]
    \f[ \Rightarrow E_0 = \frac {\sum_n E_n - B_0 f(V_n)}{n} \f]
    \n
    If \f$B_0\f$ is the optimal parameter, the partial derivative of 
    \f$\Delta\f$ with respect to \f$B_0\f$ must vanish. Thus
    \n
    \f[0 = \frac{\partial \Delta}{\partial B_0} 
      = \frac{\partial (\sum_n (E(V_n) - E_n)^2)}{\partial B_0}  
      = \frac{\partial \left(\sum_n (E_0 + B_0\cdot f(V_n) - E_n)^2\right)}{\partial B_0}
      = 2 \sum_n f(V_n) \cdot (E_0 + B_0\cdot f(V_n) - E_n)
    \f]
    Dividing by 2 and inserting the expression for \f$E_0\f$ we arrive at
    \f[
    0 =  \sum_n f(V_n)\cdot \left(\frac {\sum_n E_n - B_0 f(V_n)}{n} 
         + B_0\cdot f(V_n) - E_n\right)
      =  \sum_n f(V_n)\cdot \left[\frac{\sum E_n}{n} - E_n - B_0 \cdot\left(
         \frac{\sum_n f(V_n)}{n} - f(V_n)\right)\right]
      =  \sum_n f(V_n)\cdot \left[\frac{\sum E_n}{n} - E_n\right] 
         - B_0 \sum_n f(V_n) \cdot\left[\frac{\sum_n f(V_n)}{n} - f(V_n)\right]
    \f]
    \f[ \Rightarrow B_0 = \frac
    {\sum_n f(V_n)\cdot \left(\frac{\sum E_n}{n} - E_n\right)}
    {\sum_n f(V_n) \cdot\left(\frac{\sum_n f(V_n)}{n} - f(V_n)\right)}
    \f]
    With this optimized \f$B_0\f$, the corresponding \f$E_0\f$ is calculated as
    derived above.
    f depends implicitly on both \f$B_0'\f$ and \f$V_0\f$. The line optimization
    works with the analytically optimized values for \f$E_0\f$ and \f$B_0\f$.
 
   \ingroup  group_addons
   \author   Christoph Freysoldt, mailto:freyso@fhi-berlin.mpg.de */
class SX_EXPORT_EXT SxMurn  
{
   public:
      /// Volume at minimum \f$V_0\f$
      double minVolume;
      /// Bulk modulus \f$B_0\f$
      double bulkModulus;
      /// Bulk modulus derivative wrt pressure \f$B_0'\f$
      double bulkModulusDerivative;
      /// Minimum energy \f$E_0\f$
      double minEnergy;

      /// Given volumes
      SxVector<Double> volumes;
      /// Given energies
      SxVector<Double> energies;

      /// Constructor
      SxMurn ()
         : minVolume(0.),
           bulkModulus(0.), 
           bulkModulusDerivative(0.),
           minEnergy (0.)
           { /* empty */ } 
      /// Constructor
      SxMurn (const SxList<double> &, const SxList<double> &);
      /// Constructor
      SxMurn (const SxVector<Double> &, const SxVector<Double> &);
      /// Destructor
      ~SxMurn () {/* empty */}

      /// Fit E0, V0, B0, B0' to data
      void computeFit ();

      /// Print result (xmgrace-format)
      /** 
        @param FileName output file
        @param nPoints number of points to sample
      */
      void writeToFile (SxString FileName, int nPoints, bool pressurePrint);

      /// Get energy for volume according to Murn function
      double getEnergy (double);

      /// Get \f$\sum (E_{V_0,B_0'}(V_n) - E_n)^2 \f$
      /** \f$E_0\f$ and \f$B_0\f$ are optimized analytically as described in the
          introduction.
        */
      double getDeviation ();

      /// Returns bulk modulus in GPa
      double getBulkModulus ();
      /// Returns volume at minimum
      double getMinVolume ();
      /// Return minimum energy
      double getMinEnergy ();
      
   private:
      /// This is f from the introduction: the complicated part of the equation
      double innerFunction (double);
   protected:
      /// Internal initialization routine
      void init ();
};



#endif /* _SX_MURN_H_ */
