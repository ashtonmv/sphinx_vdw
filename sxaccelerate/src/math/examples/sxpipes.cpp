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

#include <SxUtil.h>
#include <SxList.h>
#include <SxString.h>
#include <SxOperator.h>


class Shift : public SxOperatorBase<double>
{
   public:
      double shift;

      Shift (double i) : shift(i)  { }

      virtual ~Shift () { }

   protected:
      virtual double operator* (const double &in) const  {
         return shift + in;
      }
      // in place transformation for operator*=
      virtual void applyInPlace (double &data) const
      {
         data += shift;
      }

      SXOPERATOR_GETCOPY(Shift,double);
};



class Scale : public SxOperatorBase<double>
{
   public:
      double scale;

      Scale (double s) : scale (s) { }

      virtual ~Scale () { }

   protected:
      virtual double operator* (const double &in) const  {
         return scale * in;
      }
      virtual void applyInPlace (double &data) const { SX_EXIT; }
      
      SXOPERATOR_GETCOPY(Scale,double);
};




/** \examples sxpipes.cpp

    This example demonstrates how transformations can be constructed that
    support pipelining.
    Two simple operators are defined: Shift and Scale.

    \author Sixten Boeck, boeck@sfhingx.de */
int main ()
{
   // --- initialize SFHIngX
   cout << "\tResult\tDescription" << endl;
   cout << SEPARATOR;

   // --- define new operators
   Shift T(2);   // translation
   Scale S(5);   // scaling

   // --- apply operators directly
   cout << "\t" << ( T | 1. )  << "\t 2+1     = 3"  << endl;
   cout << "\t" << ( S | 6. )  << "\t 5*6     = 30" << endl;

   double x = 7;
   // in place transformation (implemented in applyInPlace)
   x *= T;
   cout << "\t" << x           << "\t 7+2     = 9"  << endl;

   // --- order of operations: pipelining
   SxOperator<double> O1;
   O1 << T << S;   // apply: 1st: translation, 2nd: scale
   cout << "\t" << ( O1 | 1. ) << "\t (2+1)*5 = 15" << endl;
   cout << "\t" << ( O1 | 6. ) << "\t (2+6)*5 = 40" << endl;

   // --- cleaning up an old pipeline
   O1 = S;
   cout << "\t" << ( O1 | 6. ) << "\t 5*6     = 30" << endl;

   SxOperator<double> O2;
   O2 << S << T;   // apply: 1st: scale, 2nd: translation
   cout << "\t" << ( O2 | 1. ) << "\t (5*1)+2 = 7"  << endl;
   cout << "\t" << ( O2 | 6. ) << "\t (5*6)+2 = 32" << endl;
   
   // --- reverse ordering using |=
   O2  = T;
   O2 |= S; // apply: 1st: scale, 2nd: translation
   // equivalent to T | ( S | ... )
   cout << "\t" << ( O2 | 1. ) << "\t (5*1)+2 = 7"  << endl;
   cout << "\t" << ( O2 | 6. ) << "\t (5*6)+2 = 32" << endl;

   // --- pipes using operator references

   // use SxPtr
   SxPtr<Scale> pO3 = SxPtr<Scale>::create (3.);
   
   Scale & O3 = *pO3; // only for clarity

   SxOperator<double> O4;
   // put operator reference(O3) into pipe
   O4 << T << pO3;

   // note that  we added a reference to O3 to O4
   // also useful for "heavy" operators

   // side remark: O4 has an autopointer pointing to O3
   int &nRef = *pO3.refCounter; // number of additional references
   cout << "number of refs to O3: " << (nRef+1) << endl;

   // the content of pO3 (the autopointer) is therefore irrelevant for O4
   pO3 =  SxPtr<Scale> ();
   cout << "number of refs to O3: " << (nRef+1) << endl;
   // what matters is the content of O3 (managed from within O4 now)

   // O4 will execute O3
   cout << "\t" << ( O4 | 1. ) << "\t 3*(2+1) = 9" << endl;
   // IMPORTANT: change O3 => change O4
   O3.scale = 4.;
   cout << "\t" << ( O4 | 1. ) << "\t 4*(2+1) = 12" << endl;

   // clean up O4
   O4 = T;
   // note: O3 and nRef are invalid now

   return 0;
}
