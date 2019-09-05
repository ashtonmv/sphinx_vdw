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

#ifndef _SX_N_H_
#define _SX_N_H_

#include <SxGQExpr.h>

/** \brief Graph Query Node Symbol Class

    \b N = SPHInX Graph Query Node Symbol Class

     N class represents a node symbol
     within a graph query. It's
     constructor takes the property
     name as parameter.

\code
N("prop1") == "val1";

or

N("prop2") != "val2";

or

N("prop3").any();

\endcode
 */
namespace sx {
class N
{
   public:
      N ();
      N (SxString prop_);

      SxString getPropName () const;
      SxString getCName () const;

      N &c (SxString name = "");
      bool isCaptured () const;
      N &any ();
      operator SxPtr<SxGQExprBase>() const;

      template<class T>
      N &eq(const T &propVal);

      template<class T>
      N &neq(const T &propVal);

   protected:
      SxString propName;
      SxVariant val;
      SxString captureName;
      SxGQExprBase::OpType op;
      bool captureEnabled;
};

#include <SxN.hpp>

}
#endif /*_SX_N_H_*/
