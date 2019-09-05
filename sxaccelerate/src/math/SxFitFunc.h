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

#ifndef _SX_FIT_FUNC_H_
#define _SX_FIT_FUNC_H_

#include <SxMath.h>
#include <SxVector.h>

/** \brief Fitted function class


    \author Christoph Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_MATH SxFitFunc
{
   protected:
      /// The function parameters
      SxVector<Double> param;
      /// The number of variables expected
      int nVar;

      /// Constructor
      SxFitFunc (int nVarIn = -1) : nVar(nVarIn) { }

   public:
      // --- interface to the actual implementation
      /** \brief The evaluation function
          This function returns the value of the fitted function
          for the variable vector var.
        */
      virtual double evaluate (const SxVector<Double> &var) const = 0;
      /** \brief Derivative with respect to variables
          @param var input variables
          @param valPtr if non-zero, the value of the function must
                        be computed, too, and written to valPtr.
          This function returns the derivative of the fitted function
          and possibly the function value.
        */
      virtual SxVector<Double> getDeriv (const SxVector<Double> &var,
                                         double *val = NULL) const = 0;
      /** \brief Derivative with respect to parameters
          @param var input variables
          @param valPtr if non-zero, the value of the function must
                        be computed, too, and written to valPtr.
          This function returns the change in the function value
          with respect to the function parameters, and possibly the
          function value.
        */
      virtual SxVector<Double> getParamDeriv (const SxVector<Double> &var,
                                              double *val = NULL) const = 0;

      /** This function is called when the parameters change */
      virtual void paramChange () { }

      /// Stop criterium for fitting
      virtual bool fitOK (const SxVector<Double> &dParam) const = 0;

      /// Initial trial step
      virtual double paramFitStep () { return 0.1; }

   public:
      // --- Service functions

      /// Get number of parameters
      ssize_t getNParam () const { return param.getSize (); }
      /// Get number of variables
      ssize_t getNVar () const { return nVar; }
      /// Set parameters
      void setParam (const SxVector<Double> &newParam)
      {
         SX_CHECK(newParam.getSize () == getNParam () || getNParam () == 0,
                  newParam.getSize (), getNParam ());
         param.resize (newParam.getSize ());
         param <<= newParam;
         paramChange ();
      }

      double operator() (const SxVector<Double> &var) const
      {
         return evaluate (var);
      }

      /// Fit from a set of values via conjugate-gradient
      void fitFromValsCG (const SxArray<SxVector<Double> > &xList,
                          const SxArray<double> &yList);
   protected:
      /// Auxiliary function for fitFromVals
      SxVector<Double> getParamGrad (const SxArray<SxVector<Double> > &xList,
                                     const SxArray<double> &yList);

};

class SxFitPolynomial : public SxFitFunc
{
   protected:
      int order;
   public:
      /** \brief Constructor
          @param number of variables
          @param polynomial order
        */
      SxFitPolynomial (int nVar, int order);
      // --- interface to the actual implementation
      virtual double evaluate (const SxVector<Double> &var) const;
      virtual SxVector<Double> getDeriv (const SxVector<Double> &var,
                                         double *val = NULL) const;
      virtual SxVector<Double> getParamDeriv (const SxVector<Double> &var,
                                              double *val = NULL) const;
      //virtual void paramChange () { }
      virtual bool fitOK (const SxVector<Double> &dParam) const;

      double relAcc;

      inline SxVector<Double> &paramRef ()
      {
         return param;
      }

      inline int getOrder () const { return order; }

      /** \brief Apply a scale transformation
          Rewrite coefficients for a new polynomial \f$\tilde P\f$
          \f[
          \tilde P(\alpha_i x_i) = P(x_i)
          \f\
          \@param scale The factors by which each variable will be scaled.
        */
      void rescale (const SxVector<Double> &scale);
   protected:
      enum GradType { None, VarGrad, ParamGrad} ;
      /// The internal routine
      SxVector<Double> evaluate (const SxVector<Double> &var,
                                 double *val,
                                 GradType gradientType) const;
   public:
      /// Evaluation of Horner scheme
      double evaluateHorner0 (const SxVector<Double> &var) const;
      /// Evaluation of Horner scheme
      SxVector<Double> evaluate3 (const SxVector<Double> &var, int levels) const;
};

#endif /* _SX_FIT_FUNC_H_ */
