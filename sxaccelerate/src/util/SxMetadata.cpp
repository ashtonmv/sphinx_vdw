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

#include <SxMetadata.h>

void SxMetadata::removeAll ()
{
   for (int i = 0; i < metaMaps.getSize (); ++i)
      delete metaMaps(i);
   metaMaps.resize (0);
}

SxMetadata& SxMetadata::operator= (const SxMetadata &in)
{
   removeAll ();
   metaMaps.resize (in.metaMaps.getSize ());
   for (int i = 0; i < in.metaMaps.getSize (); ++i)
      metaMaps(i) = in.metaMaps(i)->getCopy ();
   return *this;
}

SxMetaBase* SxMetadata::findMap (size_t typeSize, const type_info &typeId) const
{
   for (int i = 0; i < metaMaps.getSize (); ++i)
   {
      // shortcut failure to save on dynamic casts
      if (metaMaps(i)->typeSize != typeSize) continue;

      // this implementation seems much faster than dynamic_cast
      if (metaMaps(i)->typeMatch (typeId)) return metaMaps(i);

   }
   // not found? => new one
   ssize_t n = metaMaps.getSize ();
#ifndef NDEBUG
   if (metaMaps.getSize () > 10)  {
      // so many different types used as keys?! Unify the key types!
      // probably switching systematically to SxString/const char*
      // keys is best solution if you have so many different key types.
      // You are getting slow anyway.
      cout << "ERROR: too many key types in metadata" << endl;
      SX_EXIT;
   }
#endif
   metaMaps.resize (n + 1,true);
   return NULL;
}
