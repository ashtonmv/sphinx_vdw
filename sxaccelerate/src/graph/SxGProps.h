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

#ifndef _SX_G_PROPS_H_
#define _SX_G_PROPS_H_

#include <SxVariant.h>
#include <SxString.h>
#include <SxMap.h>
#include <SxHashFunction.h>
#include <SxPair.h>

/** \brief Graph Node for property graphs

    \b SxGProps = SPHInX Properties storage Class

     SxGProps class represents a graph
     node for a property graph. A node
     in a property graph can store multiple
     key-value pairs called properties.
 */
class SxGProps
{
   public:
      SxGProps ();
      SxGProps (ssize_t id_);

      ssize_t getId () const;
      SxVariant &getProperty (const SxString &key);
      const SxVariant &getProperty (const SxString &key) const;

      const SxMap<SxString,SxVariant> &getProperties () const;

      bool operator== (const SxGProps &n1) const;

      template<class T>
      void setProperty (const SxString &key, const T &v)
      {
         props(key) = SxVariant(v);
      }

      bool hasProperty (const SxString &key) const;

      operator ssize_t() const;

      friend std::ostream &operator<< (std::ostream &os, const SxGProps &in)
      {
         os << in.getId ();
         return os;
      }

   protected:
      ssize_t id;
      SxMap<SxString,SxVariant> props;
};

size_t sxHash (const SxGProps &in);
#endif /*_SX_G_PROPS_H_*/
