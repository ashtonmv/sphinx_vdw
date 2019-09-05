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

#include <SxGProps.h>

size_t sxHash (const SxGProps &in) {
   return SxHashFunction::hash (in.getId ());
}


SxGProps::SxGProps () : id(-1) { }
SxGProps::SxGProps (ssize_t id_) : id(id_) {}

ssize_t SxGProps::getId () const
{
   return id;
}

SxVariant &SxGProps::getProperty (const SxString &key)
{
   if(props.containsKey (key))  return props(key);
   else                         SX_CHECK (false, key);
}

const SxVariant &SxGProps::getProperty (const SxString &key) const
{
   if(props.containsKey (key))  return props(key);
   else                         SX_CHECK (false, key);
}

const SxMap<SxString,SxVariant> &
SxGProps::getProperties () const
{
   return props;
}

bool SxGProps::operator== (const SxGProps &n1) const
{
   return id == n1.id;
}

bool SxGProps::hasProperty (const SxString &key) const
{
   return this->props.containsKey (key);
}

SxGProps::operator ssize_t() const {
   return id;
} 
