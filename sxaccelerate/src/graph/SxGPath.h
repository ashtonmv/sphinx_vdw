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

#ifndef _SX_G_PATH_H_
#define _SX_G_PATH_H_

#include<SxGraph.h>
#include<SxGProps.h>

/** \brief Graph Path class

    \b SxGPath = SPHInX Graph Path Class

    This class represents an individual path/selection
    and it allows to fetch the nodes based on idx.

 */
class SxGPath
{
   public:
      typedef SxPtr<SxUniqueList<ssize_t> > Selection;
      typedef SxPtr<SxList<Selection> >  SelSet;

      SxGPath ();
      SxGPath (const SxPtr<SxGraph<SxGProps> > &gPtr_,
               const Selection &sel_);

     ~SxGPath ();

      ssize_t getSize () const;

      SxGProps &operator() (ssize_t idx);

      template<class Fn>
      void foreach (Fn fn);

      template<class Fn>
      void foreach (Fn fn) const;

   protected:
      SxPtr<SxGraph<SxGProps> > gPtr;
      Selection sel;
};

#include <SxGPath.hpp>

#endif /* _SX_G_PATH_H_ */
