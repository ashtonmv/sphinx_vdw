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

#include <SxVariant.h>

void SxVariant::setType (int type_)
{
   if (type != type_)  {
      reset ();
      type = type_;
      switch (type)  {
         case SxVariantType::Undefined:
         case SxVariantType::Any:
            dataType = SxVariantType::Undefined;
            break;
         case SxVariantType::Int:
            dataType = SxVariantType::Int;
            data   = (void *)new int64_t;
            minPtr = (void *)new int64_t;
            maxPtr = (void *)new int64_t;
            *(int64_t *)data   = 0;
            *(int64_t *)minPtr = 0;
            *(int64_t *)maxPtr = 0;
            break;
         case SxVariantType::Double:
            dataType = SxVariantType::Double;
            data   = (void *)new double;
            minPtr = (void *)new double;
            maxPtr = (void *)new double;
            *(double *)data   = 0.0;
            *(double *)minPtr = 0.0;
            *(double *)maxPtr = 0.0;
            break;
         case SxVariantType::Bool:
            dataType = SxVariantType::Bool;
            data   = (void *)new bool;
            minPtr = (void *)new uint8_t;
            maxPtr = (void *)new uint8_t;
            *(bool *)data      = false;
            *(uint8_t *)minPtr = 0.0;
            *(uint8_t *)maxPtr = 0.0;
            break;
         case SxVariantType::String:
            dataType = SxVariantType::String;
            data   = (void *)new SxString;
            minPtr = (void *)new char;
            maxPtr = (void *)new char;
            *(char *)minPtr = 0;
            *(char *)maxPtr = 0;
            break;
         case SxVariantType::List:
            dataType = SxVariantType::List;
            data   = (void *)new SxList<SxVariant>;
            minPtr = (void *)new int;
            maxPtr = (void *)new int;
            *(int *)minPtr = SxVariantType::Undefined;
            *(int *)maxPtr = SxVariantType::Undefined;
            initialized = true;
            break;
         case SxVariantType::Group:
            break;
         default: SX_EXIT;
      }
      resetLimits ();
   }
}

void SxVariant::reset ()
{
   if (data)  {
      switch (dataType)  {
         case SxVariantType::Undefined:
         case SxVariantType::Any:
            SX_EXIT;
            break;
         case SxVariantType::Int:
            delete ((int64_t *)data);
            delete ((int64_t *)minPtr);
            delete ((int64_t *)maxPtr);
            break;
         case SxVariantType::Double:
            delete ((double *)data);
            delete ((double *)minPtr);
            delete ((double *)maxPtr);
            break;
         case SxVariantType::Bool:
            delete ((bool *)data);
            delete ((uint8_t *)minPtr);
            delete ((uint8_t *)maxPtr);
            break;
         case SxVariantType::String:
            delete ((SxString *)data);
            delete ((char *)minPtr);
            delete ((char *)maxPtr);
            break;
         case SxVariantType::List:
            delete ((SxList<SxVariant> *)data);
            delete ((int *)minPtr);
            delete ((int *)maxPtr);
            break;
         default: SX_EXIT;
      }
      data = NULL;
      minPtr = NULL;
      maxPtr = NULL;
      initialized = false;
      type = SxVariantType::Undefined;
      dataType = SxVariantType::Undefined;
   }
}

void SxVariant::resetLimits ()
{
   if (data)  {
      switch (dataType)  {
         case SxVariantType::Undefined:
         case SxVariantType::Any:
            break;
         case SxVariantType::Int:
            *(int64_t *)minPtr = SxVariantType::ScalarType<int>::minVal ();
            *(int64_t *)maxPtr = SxVariantType::ScalarType<int>::maxVal ();
            break;
         case SxVariantType::Double:
            *(double *)minPtr = SxVariantType::ScalarType<double>::minVal ();
            *(double *)maxPtr = SxVariantType::ScalarType<double>::maxVal ();
            break;
         case SxVariantType::Bool:
            *(uint8_t *)minPtr = SxVariantType::ScalarType<bool>::minVal ();
            *(uint8_t *)maxPtr = SxVariantType::ScalarType<bool>::maxVal ();
            break;
         case SxVariantType::String:
            *(char *)minPtr = SxVariantType::ScalarType<SxString>::minVal ();
            *(char *)maxPtr = SxVariantType::ScalarType<SxString>::maxVal ();
            break;
         case SxVariantType::List:
            *(int *)minPtr = SxVariantType::Undefined;
            *(int *)maxPtr = SxVariantType::Undefined;
            break;
         default: SX_EXIT;
      }
      minElem = 0;
      maxElem = 0;
   }
}

void SxVariant::operator= (const SxVariant &in)
{
   if (&in == this)  return;

   setType (in.type);

   switch (dataType) {
      case SxVariantType::Undefined:
      case SxVariantType::Any:
         break;
      case SxVariantType::Int:
         *(int64_t *)data   = *(int64_t *)in.data;
         *(int64_t *)minPtr = *(int64_t *)in.minPtr;
         *(int64_t *)maxPtr = *(int64_t *)in.maxPtr;
         break;
      case SxVariantType::Double:
         *(double *)data   = *(double *)in.data;
         *(double *)minPtr = *(double *)in.minPtr;
         *(double *)maxPtr = *(double *)in.maxPtr;
         break;
      case SxVariantType::Bool:
         *(bool *)data      = *(bool *)in.data;
         *(uint8_t *)minPtr = *(uint8_t *)in.minPtr;
         *(uint8_t *)maxPtr = *(uint8_t *)in.maxPtr;
         break;
      case SxVariantType::String:
         *(SxString *)data  = *(SxString *)in.data;
         *(char *)minPtr = *(char *)in.minPtr;
         *(char *)maxPtr = *(char *)in.maxPtr;
         break;
      case SxVariantType::List:
         *(SxList<SxVariant> *)data = *(SxList<SxVariant> *)in.data;
         *(int *)minPtr = *(int *)in.minPtr;
         *(int *)maxPtr = *(int *)in.maxPtr;
         break;
      default: SX_EXIT;
   }

   initialized = in.initialized;
   minElem = in.minElem;
   maxElem = in.maxElem;
   regex   = in.regex;
   regexPattern = in.regexPattern;
   regexFlags   = in.regexFlags;

   tag = in.tag;
}

bool SxVariant::isValidString (const SxString &in) const
{
   if (minElem && in.getSize () < minElem)
      return false;
   else if (maxElem && in.getSize () > maxElem)
      return false;
   if (regex.getPtr() != NULL)  {
      try { 
         return (regex->match(in).getSize() == 1);
      }  catch (SxException e) {
         return false;
      }
   }

   return true;
}

void SxVariant::validate (const SxString &in) const
{
   if (minElem && in.getSize () < minElem)
      SX_THROW("String is below minimal length " + SxString(minElem));
   else if (maxElem && in.getSize () > maxElem)
      SX_THROW("String exceeds maximum length " + SxString(maxElem));
}

bool SxVariant::hasStringLimits () const
{
   if (regex.getPtr() != NULL || minElem || maxElem) return true;
   else                                              return false;
}

bool SxVariant::matchLimits (const SxVariant &in) const
{
   using namespace SxVariantType;

   int typeB = in.getType ();

   if (type == Int && typeB == Int)  {
      int64_t val = in.getInt ();
      return !(val < minInt() || val > maxInt());
   } else if (type == Int && typeB == Double)  {
      int64_t val = static_cast<int64_t>(in.getDouble ());
      return !(val < minInt() || val > maxInt());
   } else if (type == Double && typeB == Double)  {
      double val = in.getDouble ();
      return !(val <  minDouble() || val > maxDouble());
   } else if (type == Double && typeB == Int)  {
      double val = static_cast<double>(in.getInt ());
      return !(val <  minDouble() || val > maxDouble());
   } else if (type == String && typeB == String)  {
      if (hasStringLimits ())
         return isValidString (in.toString ());
      else
         return true;
   }

   return false;
}

bool SxVariant::match (const SxVariant &in, bool withLimits_) const
{
   using namespace SxVariantType;

   if (type == SxVariantType::Any)  {
      return true;
   } else if (type == in.getType ())  {
      if (withLimits_)  return matchLimits(in);
      else              return true;

   } else if (   (type == Int && in.getType () == Double)
              || (type == Double && in.getType () == Int))
   {
      if (withLimits_)  return matchLimits(in);
      else              return true;
   }

   return false;
}

SxVariant::CastType SxVariant::castDiff (const SxVariant &in) const
{
   using namespace SxVariantType;
   if (type == in.type)       return  ExactCast;
   if (type == Any)           return  ExactCast;
   if (in.type == Undefined)  return  IncompatibleTypes;

   if (type == Int)  {
      switch (in.type)  {
         case Double : return ImplicitCast;
         case Bool   : return ImplicitCast;
         case String : return ExplicitCastRequired;
         case List   : return ExplicitCastRequired;
      }
   } else if (type == Double)  {
      switch (in.type)  {
         case Int    : return ExplicitCastRequired;
         case Bool   : return ImplicitCast;
         case String : return ExplicitCastRequired;
         case List   : return ExplicitCastRequired;
      }
   } else if (type == Bool)  {
      switch (in.type)  {
         case Int    : return ImplicitCast;
         case Double : return ImplicitCast;
         case String : return ExplicitCastRequired;
         case List   : return ExplicitCastRequired;
      }
   } else if (type == String)  {
      switch (in.type)  {
         case Int    : return  ExplicitCastRequired;
         case Double : return  ExplicitCastRequired;
         case List   : return  ExplicitCastRequired;
      }
   }

   SX_EXIT;
}


SxString SxVariant::getTypeName () const
{
   switch (type) {
      case SxVariantType::Undefined:
         return "undefined";
      case SxVariantType::Any:
         return "any";
      case SxVariantType::Int:
         return "int";
      case SxVariantType::Double:
         return "float";
      case SxVariantType::Bool:
         return "bool";
      case SxVariantType::String:
         return "string";
      case SxVariantType::List:
         return "list";
      case SxVariantType::Group:
         return "group";
      default: SX_EXIT;
   }
   return "unknown";
}

SxString SxVariant::getTypeArgs () const
{
   SxString res;

   SxString typeStr;
   SxString rangeStr;
   SxString minStr;
   SxString maxStr;

   typeStr = getTypeName ();
   switch (type)  {
      case SxVariantType::Int:
         if (*(int64_t *)minPtr != SxVariantType::ScalarType<int>::minVal())
            minStr = SxString(*(int64_t *)minPtr);
         if (*(int64_t *)maxPtr != SxVariantType::ScalarType<int>::maxVal()) 
            maxStr = SxString(*(int64_t *)maxPtr);
         break;
      case SxVariantType::Double:
         if (*(double *)minPtr != SxVariantType::ScalarType<double>::minVal())
            minStr = SxString(*(double *)minPtr);
         if (*(double *)maxPtr != SxVariantType::ScalarType<double>::maxVal()) 
            maxStr = SxString(*(double *)maxPtr);
         break;
      case SxVariantType::String:
         if (minElem > 0)
            minStr = SxString(minElem);
         if (maxElem > 0)
            maxStr = SxString(maxElem);
         break;
      case SxVariantType::Bool:
         minStr = SxString ("false");
         maxStr = SxString ("true");
         break;
      case SxVariantType::List:
      case SxVariantType::Group:
      case SxVariantType::Undefined:
      case SxVariantType::Any:
         break;
      default: SX_EXIT;
   }

   if (minStr != "" || maxStr != "")  {
      rangeStr = "(";
      if (minStr != "")
         rangeStr += "min " + minStr;
      if (maxStr != "")  {
         if (minStr != "")
            rangeStr += ", ";
         rangeStr += "max " + maxStr;
      }
      rangeStr += ")";
   }
   if (regex.getPtr() != NULL)
      rangeStr += "[/" + regexPattern + "/" + regexFlags + "]";

   return typeStr + rangeStr;
}

SxString SxVariant::getDescription (bool withTypeArgs) const
{
   SxString res;

   if (tag != "")  {
      if (withTypeArgs)  res = tag + ": ";
      else               res = "<" + tag + ">";
   }

   SxString argsStr;
   if (withTypeArgs)  argsStr = getTypeArgs ();

   return res + argsStr;
}

SxString SxVariant::printToString () const
{
   SxList<SxString> lines;
   SxList<SxVariant>::ConstIterator it;

   if (type != SxVariantType::List && !initialized)
      return "uninitialized";

   switch (type)  {
      case SxVariantType::Undefined:
         return "undefined";
      case SxVariantType::Any:
         return "any";
      case SxVariantType::Int:
         return SxString(*(int64_t *)data);
      case SxVariantType::Double:
         return SxString(*(double *)data);
      case SxVariantType::Bool:
         return ( (*(bool *)data)==true ?
                  SxString("true") : SxString("false") );
      case SxVariantType::String:
         return "\'" + (*(SxString *)data) + "\'";
      case SxVariantType::List:
         it = ((SxList<SxVariant> *)data)->begin ();
         while (it.isValid ())  {
            lines.append (it->printToString ());
            ++it;
         }
         return "[" + SxString::join (lines, ", ") + "]";
      default: SX_EXIT;
   }   

   return "unknown";
}


void SxVariant::undefTest (const SxVariant &a, const SxVariant &b) const
{
   undefTest (a);
   undefTest (b);
}

void SxVariant::undefTest (const SxVariant &a) const
{
   if (a.getType () == SxVariantType::Undefined)
      SX_THROW ("Undefined type '" + a.getTag () + "'");

   if (a.getType () == SxVariantType::Any)
      SX_THROW ("Undefined any type '" + a.getTag () + "'");

   if (!a.isInitialized ())
      SX_THROW ("Uninitialized value '" + a.getTag () + "'");
}

// --------------------------------------------------------------------------

SxVariant operator+ (const SxVariant &a, const SxVariant &b)
{
   using namespace SxVariantType;

   int typeA = a.getType();
   int typeB = b.getType();
   if (typeA == Int)  {
      if (typeB == Int)     return a.toInt()         + b.toInt();
      if (typeB == Double)  return (double)a.toInt() + b.toDouble();
   } else if (typeA == Double)  {
      if (typeB == Int)     return a.toDouble() + (double)b.toInt();
      if (typeB == Double)  return a.toDouble() + b.toDouble();
   } else if (typeA == String)  {
      if (typeB == String)  return a.toString() + b.toString();
   }
   return SxVariant();
}

// --------------------------------------------------------------------------

SxString SxVariant::getTypeStr (int type) {
   switch (type) {
      case SxVariantType::Undefined:
         return "undefined";
      case SxVariantType::Any:
         return "any";
      case SxVariantType::Int:
         return "int";
      case SxVariantType::Double:
         return "float";
      case SxVariantType::Bool:
         return "bool";
      case SxVariantType::String:
         return "string";
      case SxVariantType::List:
         return "list";
      case SxVariantType::Group:
         return "group";
      default: SX_EXIT;
   }
   return "unknown";
}

int SxVariant::getTypeId (const SxString &str) {
   if (str == "group") {
      return SxVariantType::Group;
   } else if (str == "list") {
      return SxVariantType::List;
   } else if (str == "int") {
      return SxVariantType::Int;
   } else if (str == "float") {
      return SxVariantType::Double;
   } else if (str == "string") {
      return SxVariantType::String;
   } else if (str == "bool") {
      return SxVariantType::Bool;
   } else if (str == "group") {
      return SxVariantType::Group;
   } else {
      return SxVariantType::Undefined;
   }
}
