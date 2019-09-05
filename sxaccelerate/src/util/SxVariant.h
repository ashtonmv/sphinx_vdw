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

#ifndef _SX_VARIANT_H_
#define _SX_VARIANT_H_

#include <SxUtil.h>
#include <SxPtr.h>
#include <SxRegex.h>

namespace SxVariantType {

   enum DataType { Undefined, Any, Int, Double, Bool, String, List, Group };

   template<class T, size_t N=0> class ScalarType { };

   template<> class ScalarType<int> {
      public:
         enum { Type = Int };
         typedef int64_t AtomType;
         inline static int64_t minVal () { return INT64_MIN ;}
         inline static int64_t maxVal () { return INT64_MAX; }
   };
   template<> class ScalarType<long> {
      public:
         enum { Type = Int };
         typedef int64_t AtomType;
         inline static int64_t minVal () { return INT64_MIN ;}
         inline static int64_t maxVal () { return INT64_MAX; }
   };
   template<> class ScalarType<long long> {
      public:
         enum { Type = Int };
         typedef int64_t AtomType;
         inline static int64_t minVal () { return INT64_MIN ;}
         inline static int64_t maxVal () { return INT64_MAX; }
   };
   template<> class ScalarType<bool> {
      public:
         enum { Type = Bool };
         typedef uint8_t AtomType;
         inline static uint8_t minVal () { return 0; }
         inline static uint8_t maxVal () { return 1; }
   };

   template<> class ScalarType<float> {
      public:
         enum { Type = Double };
         typedef double AtomType;
         inline static float minVal () { return -HUGE_VAL; }
         inline static float maxVal () { return  HUGE_VAL; }
   };
   template<> class ScalarType<double> {
      public:
         enum { Type = Double };
         typedef double AtomType;
         inline static double minVal () { return -HUGE_VAL; }
         inline static double maxVal () { return  HUGE_VAL; }
   };

   template<> class ScalarType<const char *> {
      public:
         enum { Type = String };
         typedef char AtomType;
         inline static char minVal () { return CHAR_MIN; }
         inline static char maxVal () { return CHAR_MAX; }
   };
   template<> class ScalarType<char[]> {
      public:
         enum { Type = String };
         typedef char AtomType;
         inline static char minVal () { return CHAR_MIN; }
         inline static char maxVal () { return CHAR_MAX; }
   };
   template<size_t N> class ScalarType<char[N]> {
      public:
         enum { Type = String };
         typedef char AtomType;
         inline static char minVal () { return CHAR_MIN; }
         inline static char maxVal () { return CHAR_MAX; }
   };
   template<> class ScalarType<SxString> {
      public:
         enum { Type = String };
         typedef char AtomType;
         inline static char minVal () { return CHAR_MIN; }
         inline static char maxVal () { return CHAR_MAX; }
   };
}

/** \brief Variant of simplest data types

    This class provides a variant of simple data types. Using this variant
    lists or arrays with elements of varying types can be created at runtime.

    \author Sixten Boeck */
class SX_EXPORT_UTIL SxVariant
{
   public:

      enum CastType { ExactCast, ImplicitCast, ExplicitCastRequired,
                      IncompatibleTypes };

      SxVariant ();

      SxVariant (const SxVariant &in);
      void operator= (const SxVariant &in);

      template<class T> SxVariant (const T &in, const SxString &tag="");

      template<class A, class B, class T> SxVariant (const T &in,
                                                     const A &min,
                                                     const B &max,
                                                     const SxString &tag="");

      ~SxVariant ();

      bool  operator== (const SxVariant &) const;
      bool  operator!= (const SxVariant &) const;
      bool  operator< (const SxVariant &) const;

      // --- Set value (with type/range/regex check)
      template<class T> void operator= (const T in);
                        void operator= (const SxString &in);

      template<class T> void set (const T in);
                        void set (const SxString &in);
                        void set (const SxVariant &, bool allowTypecast=false);
                        void set (const SxList<SxVariant> &in);

      // --- List
      void append (const SxVariant &in);

      int     getRank () const;
      ssize_t getListSize () const;

      // --- List iterators
      SxList<SxVariant>::Iterator begin ();
      SxList<SxVariant>::ConstIterator begin () const;
      SxList<SxVariant>::Iterator end ();
      SxList<SxVariant>::ConstIterator end () const;

      // --- Limits
      template<class A, class B>
      SxVariant &setLimits (const A &min, const B &max);

      template<class A, class B>
      SxVariant &getLimits (A *min, B *max);

      template<class A> SxVariant &setMin (const A &min);
      template<class A> SxVariant &setMax (const A &max);

      void resetLimits ();

      SxVariant &setRegex (const SxString &, const SxString & = "");
      SxString getRegex () const;

      // --- Match
      bool match (const SxVariant &in, bool withLimits=false) const;
      bool matchLimits (const SxVariant &in) const;
      CastType castDiff (const SxVariant &in) const;

      // --- Values 
      int64_t getInt () const;
      double getDouble () const;
      bool   getBool () const;
      const SxString &getString () const;

      // --- Values with cast
      int64_t toInt () const;
      double toDouble () const;
      bool   toBool ()  const;
      SxString toString () const;
      SxString toByteArray () const;

      // --- explicit type initialization
      void setType (int type);

      bool isInitialized () const;
      int getType () const;
      SxString getTypeName () const;
      int getDataType () const;
      const SxString &getTag () const;
      void setTag (const SxString &);

      SxString printToString () const;
      SxString getTypeArgs () const;
      SxString getDescription (bool withTypeArgs=false) const;


      static SxString getTypeStr (int type);
      static int      getTypeId  (const SxString &str);

   protected:

      int type;
      int dataType;
      bool initialized;
      void *data;

      void *minPtr, *maxPtr;
      ssize_t minElem;
      ssize_t maxElem;
      SxPtr<SxRegex> regex;
      SxString regexPattern, regexFlags;

      SxString tag;

      template<class T> void setTypeT (const T in);
                        void setTypeT (const SxString &in);


      bool hasStringLimits () const;
      bool isValidString (const SxString &) const;
      void validate (const SxString &) const;

      void undefTest (const SxVariant &a, const SxVariant &b) const;
      void undefTest (const SxVariant &a) const;

      void reset ();

      int64_t minInt () const;
      int64_t maxInt () const;
      double minDouble () const;
      double maxDouble () const;
};

#ifdef WIN32
   SX_EXPORT_UTIL std::wostream &operator<< (std::wostream &, const SxVariant &);
#endif
SX_EXPORT_UTIL std::ostream &operator<< (std::ostream &, const SxVariant &);
SX_EXPORT_UTIL SxVariant operator+ (const SxVariant &, const SxVariant &);

#include <SxVariant.hpp>

#endif /* _SX_VARIANT_H_ */
