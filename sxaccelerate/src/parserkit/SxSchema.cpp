#include <SxSchema.h>

SxSchema::SxSchema () { }

SxSchema::~SxSchema () { }

SxSchema::SxSchema (const SxPtr<SxGraph<SxGProps> > &schemaG_)
{
   SX_CHECK (schemaG_.getPtr ());

   schemaG  = schemaG_;
   simplifySchemaAst ();
}

ssize_t SxSchema::getInt (const SxGraph<SxGProps>::Iterator &sIt,
                          const SxString &key,
                          const ssize_t &default_)
{
   auto it = get (sIt, "__sx_Key", key);
   if (it.isValid ())
      return it->getProperty ("__sx_Value").getInt ();
   else
      return default_;
}

double SxSchema::getDouble (const SxGraph<SxGProps>::Iterator &sIt,
                            const SxString &key,
                            const double &default_)
{
   auto it = get (sIt, "__sx_Key", key);
   if (it.isValid ())
      return it->getProperty ("__sx_Value").getDouble ();
   else
      return default_;
}

SxString SxSchema::getString (const SxGraph<SxGProps>::Iterator &sIt,
                              const SxString &key,
                              const SxString &default_)
{
   auto it = get (sIt, "__sx_Key", key);
   if (it.isValid ())
      return it->getProperty ("__sx_Value").getString ();
   else
      return default_;
}

void SxSchema::validationError (const SxString &msg,
                                const SxString &sTag,
                                const SxString &dTag)
{
   SX_THROW (dTag+": "+msg+" declared at "+sTag);
}


bool SxSchema::validateArray (const SxGraph<SxGProps>::Iterator &sIt,
                              const SxGraph<SxGProps>::Iterator &dIt)
{
   auto itemTypeIt = get(sIt, "__sx_Key", "itemsType");
   SX_CHECK (itemTypeIt.isValid ());

   Type itype = (Type)itemTypeIt->getProperty ("type").getInt ();
   ssize_t nItemsFound = 0;

   if (itype == Type::Group) {
      for (ssize_t l = 0; l < dIt.getSizeOut (); ++l) {
         auto elemIt = dIt.out (l);
         if (elemIt->getProperty ("__sx_Value").getType () != itype) {
            validationError ("type mismatch "
                             + elemIt->getProperty ("__sx_Value").getTypeName ()
                             + " != " + SxVariant::getTypeStr(itype),
                             itemTypeIt->getProperty("__sx_Value").getTag (),
                             elemIt->getProperty ("__sx_Value").getTag ());
         }
         if (validateObj (itemTypeIt, elemIt))
            nItemsFound++;
      }

   } else if (itype == Type::List) {
      auto subArrayIt = get(sIt, "__sx_Key", "subarray");
      for (ssize_t l = 0; l < dIt.getSizeOut (); ++l) {
         auto elemIt = dIt.out (l);
         if (elemIt->getProperty ("__sx_Value").getType () != itype) {
            validationError ("type mismatch "
                             + elemIt->getProperty ("__sx_Value").getTypeName ()
                             + " != " + SxVariant::getTypeStr(itype),
                             itemTypeIt->getProperty ("__sx_Value").getTag (),
                             elemIt->getProperty ("__sx_Value").getTag ());
         }
         if (validateArray (subArrayIt, elemIt))
            nItemsFound++;
      }

   } else if (itype == Type::Int) {
      for (auto elemIt = dIt->getProperty ("__sx_Value").begin ();
           elemIt.isValid (); ++elemIt) {

         if (elemIt->getType () != itype) {
            validationError ("type mismatch "
                             + elemIt->getTypeName ()
                             + " != " + SxVariant::getTypeStr(itype),
                             itemTypeIt->getProperty ("__sx_Value").getTag (),
                             elemIt->getTag ());
         }

         ssize_t val = elemIt->getInt ();

         if (!itemTypeIt->getProperty ("__sx_Value").matchLimits (val)) {
            ssize_t minVal = 0;
            ssize_t maxVal = 0;
            itemTypeIt->getProperty ("__sx_Value").getLimits (&minVal, &maxVal);
            validationError (SxString("value of type 'int' is out of range (")
                             + minVal + ", "
                             + maxVal + ")",
                             itemTypeIt->getProperty ("__sx_Value").getTag (),
                             elemIt->getTag ());
         }
         nItemsFound++;
      }

   } else if (itype == Type::Double) {
      for (auto elemIt = dIt->getProperty ("__sx_Value").begin ();
           elemIt.isValid (); ++elemIt) {

         if (elemIt->getType () != itype) {
            validationError ("type mismatch "
                             + elemIt->getTypeName ()
                             + " != " + SxVariant::getTypeStr(itype),
                             itemTypeIt->getProperty ("__sx_Value").getTag (),
                             elemIt->getTag ());
         }
         double val = elemIt->getDouble ();

         if (!itemTypeIt->getProperty ("__sx_Value").matchLimits (val)) {
            double minVal = 0;
            double maxVal = 0;
            itemTypeIt->getProperty ("__sx_Value").getLimits (&minVal, &maxVal);
            validationError (SxString ("value of type 'float' is out of range (")
                             + minVal + ", "
                             + maxVal + ")",
                             itemTypeIt->getProperty ("__sx_Value").getTag (),
                             elemIt->getTag ());
         }
         nItemsFound++;
      }

   } else if (itype == Type::Bool) {
      for (auto elemIt = dIt->getProperty ("__sx_Value").begin ();
           elemIt.isValid (); ++elemIt) {

         if (elemIt->getType () != itype) {
            validationError ("type mismatch "
                             + elemIt->getTypeName ()
                             + " != " + SxVariant::getTypeStr(itype),
                             itemTypeIt->getProperty ("__sx_Value").getTag (),
                             elemIt->getTag ());
         }

         nItemsFound++;
      }

   } else if (itype == Type::String) {

      for (auto elemIt = dIt->getProperty ("__sx_Value").begin ();
         elemIt.isValid (); ++elemIt) {

         if (elemIt->getType () != itype) {
            validationError ("type mismatch "
                             + elemIt->getTypeName ()
                             + " != " + SxVariant::getTypeStr(itype),
                             itemTypeIt->getProperty ("__sx_Value").getTag (),
                             elemIt->getTag ());
         }
         SxString val = elemIt->getString ();

         if (!itemTypeIt->getProperty ("__sx_Value").matchLimits (val)) {
            ssize_t minLen = 0;
            ssize_t maxLen = 0;
            itemTypeIt->getProperty ("__sx_Value").getLimits (&minLen, &maxLen);
            validationError (SxString("string length is out of range (")
                             + minLen + ", "
                             + maxLen + ")",
                             itemTypeIt->getProperty ("__sx_Value").getTag (),
                             elemIt->getTag ());
         }
         nItemsFound++;
      }

   }

   SxString opt = sIt->hasProperty ("optional")?
                  sIt->getProperty ("optional").getString () :
                  "";

   if (nItemsFound == 0 && (opt == "" || opt == "false")) {
      validationError ("array is empty",
                       sIt->getProperty ("__sx_Value").getTag (),
                       dIt->getProperty ("__sx_Value").getTag ());
   }

   ssize_t n        = sIt->hasProperty ("nItems")?
   sIt->getProperty ("nItems").getInt () :
   max<ssize_t>();
   ssize_t maxItems = sIt->hasProperty ("maxItems")?
   sIt->getProperty ("maxItems").getInt () :
   max<ssize_t>();
   ssize_t minItems = sIt->hasProperty ("minItems")?
   sIt->getProperty ("minItems").getInt () :
   min<ssize_t>();
   if (n != max<ssize_t>()) {
      if (nItemsFound > n) {
         validationError ("too many array items",
                          sIt->getProperty ("__sx_Value").getTag (),
                          dIt->getProperty ("__sx_Value").getTag ());
      }
      if (nItemsFound < n) {
         validationError ("too few array items",
                          sIt->getProperty ("__sx_Value").getTag (),
                          dIt->getProperty ("__sx_Value").getTag ());
      }
   }
   if (nItemsFound > maxItems) {
      validationError ("too many array items",
                       sIt->getProperty ("__sx_Value").getTag (),
                       dIt->getProperty ("__sx_Value").getTag ());
   }
   if (nItemsFound < minItems) {
      validationError ("too few array items",
                       sIt->getProperty ("__sx_Value").getTag (),
                       dIt->getProperty ("__sx_Value").getTag ());
   }
   return true;
}



bool SxSchema::validateObj (const SxGraph<SxGProps>::Iterator &sIt,
                            const SxGraph<SxGProps>::Iterator &dIt)
{
   SxMap<ssize_t, SxList<ssize_t> > foundGroups;
   SxMap<SxString, ssize_t> nameToIdx;
   for (ssize_t i = 0; i < sIt.getSizeOut (); ++i) {
      auto chIt = sIt.out (i);
      Type type = (Type)chIt->getProperty ("type").getInt ();

      if (type == Type::Group) {
         const SxVariant &key = chIt->getProperty ("__sx_Key");
         nameToIdx (key.getString ()) = chIt.getIdx ();
         foundGroups(chIt.getIdx ()); // insert empty matched list

         for (ssize_t k = 0; k < dIt.getSizeOut (); ++k) {
            auto dchIt = dIt.out (k);
            if (dchIt->getProperty ("__sx_Key") == key) {
               foundGroups(chIt.getIdx ()).append (dchIt.getIdx ());
            }
         }
      } else if (type == Type::List) {
         const SxVariant &key = chIt->getProperty ("__sx_Key");
         bool keyFound = false;
         for (ssize_t k = 0; k < dIt.getSizeOut (); ++k) {
            auto dchIt = dIt.out (k);
            if (dchIt->getProperty ("__sx_Key") == key) {
               keyFound = true;
               validateArray (chIt, dchIt);
            }
         }
         SxString opt = chIt->hasProperty ("optional")?
            chIt->getProperty ("optional").getString () :
            "";
         if (!keyFound && opt != "true") {
            validationError ("key '"
                             + key.getString ()
                             + "' of type 'array' not found in '"
                             + dIt->getProperty ("__sx_Key").getString ()
                             + "'", chIt->getProperty ("__sx_Value").getTag (),
                             dIt->getProperty ("__sx_Value").getTag ());
         }

      } else if (type == Type::Int) {
         SxString k = chIt->getProperty ("__sx_Key").getString ();

         if (!dIt->hasProperty (k)) {
            validationError ("key '"
                             + k
                             + "' of type 'int' not found in '"
                             + dIt->getProperty ("__sx_Key").getString ()
                             + "'", chIt->getProperty ("__sx_Value").getTag (),
                             dIt->getProperty ("__sx_Value").getTag ());
         }

         if (dIt->getProperty (k).getType () != Type::Int) {
            validationError ("type mismatch "
                             + dIt->getProperty (k).getTypeName ()
                             + " != int",
                             chIt->getProperty ("__sx_Value").getTag (),
                             dIt->getProperty (k).getTag ());
         }

         nameToIdx (k) = chIt.getIdx ();

         ssize_t val = dIt->getProperty (k).getInt ();

         if (!chIt->getProperty ("__sx_Value").matchLimits (val)) {
            ssize_t minVal = 0;
            ssize_t maxVal = 0;
            chIt->getProperty ("__sx_Value").getLimits (&minVal, &maxVal);
            validationError (SxString("value of type 'int' is out of range (")
                             + minVal + ", "
                             + maxVal + ")",
                             chIt->getProperty ("__sx_Value").getTag (),
                             dIt->getProperty (k).getTag ());
         }

      } else if (type == Type::Double) {
         SxString k = chIt->getProperty ("__sx_Key").getString ();

         if (!dIt->hasProperty (k)) {
            validationError ("key '"
                             + k
                             + "' of type 'float' not found in '"
                             + dIt->getProperty ("__sx_Key").getString ()
                             + "'",
                             chIt->getProperty ("__sx_Value").getTag (),
                             dIt->getProperty ("__sx_Value").getTag ());
         }

         if (dIt->getProperty (k).getType () != Type::Double) {
            validationError ("type mismatch "
                             + dIt->getProperty (k).getTypeName ()
                             + " != float",
                             chIt->getProperty ("__sx_Value").getTag (),
                             dIt->getProperty (k).getTag ());
         }

         nameToIdx (k) = chIt.getIdx ();

         double val = dIt->getProperty (k).getDouble ();

         if (!chIt->getProperty ("__sx_Value").matchLimits (val)) {
            double minVal = 0;
            double maxVal = 0;
            chIt->getProperty ("__sx_Value").getLimits (&minVal, &maxVal);
            validationError (SxString("value of type 'float' is out of range (")
                             + minVal + ", "
                             + maxVal + ")",
                             chIt->getProperty ("__sx_Value").getTag (),
                             dIt->getProperty (k).getTag ());
         }

      } else if (type == Type::Bool) {
         SxString k = chIt->getProperty ("__sx_Key").getString ();

         if (!dIt->hasProperty (k)) {
            validationError ("key '"
                             + k
                             + "' of type 'bool' not found in '"
                             + dIt->getProperty ("__sx_Key").getString ()
                             + "'",
                             chIt->getProperty ("__sx_Value").getTag (),
                             dIt->getProperty ("__sx_Value").getTag ());
         }

         if (dIt->getProperty (k).getType () != Type::Bool) {
            validationError ("type mismatch "
                             + dIt->getProperty (k).getTypeName ()
                             + " != bool",
                             chIt->getProperty ("__sx_Value").getTag (),
                             dIt->getProperty (k).getTag ());
         }

         nameToIdx (k) = chIt.getIdx ();

      } else if (type == Type::String) {
         SxString k = chIt->getProperty ("__sx_Key").getString ();

         if (!dIt->hasProperty (k)) {
            validationError ("key '"
                             + k
                             + "' of type 'string' not found in '"
                             + dIt->getProperty ("__sx_Key").getString ()
                             + "'",
                             chIt->getProperty ("__sx_Value").getTag (),
                             dIt->getProperty ("__sx_Value").getTag ());
         }

         if (dIt->getProperty (k).getType () != Type::String) {
            validationError ("type mismatch "
                             + dIt->getProperty (k).getTypeName ()
                             + " != string",
                             chIt->getProperty ("__sx_Value").getTag (),
                             dIt->getProperty (k).getTag ());
         }

         nameToIdx (k) = chIt.getIdx ();

         SxString val = dIt->getProperty (k).getString ();

         if (!chIt->getProperty ("__sx_Value").matchLimits (val)) {
            ssize_t minLen = 0;
            ssize_t maxLen = 0;
            chIt->getProperty ("__sx_Value").getLimits (&minLen, &maxLen);
            validationError (SxString("string length is out of range (")
                             + minLen + ", "
                             + maxLen + ")",
                             chIt->getProperty ("__sx_Value").getTag (),
                             dIt->getProperty (k).getTag ());
         }

      } else {
         validationError ("type '" + SxVariant::getTypeStr(type)
                          + "' not recognized",
                          chIt->getProperty ("__sx_Value").getTag (),
                          chIt->getProperty ("__sx_Value").getTag ());
      }

   }

   auto groupIt  = foundGroups.getKeys ().begin ();
   auto matchIt = foundGroups.getValues ().begin ();
   for (; groupIt.isValid (); (++groupIt,++matchIt)) {
      auto chIt = schemaG->begin (*groupIt);
      const SxVariant &key = chIt->getProperty ("__sx_Key");
      ssize_t nItemsFound = 0;
      for (auto mIt = matchIt->begin (); mIt != matchIt->end (); ++mIt) {
         auto dchIt = dataG->begin (*mIt);
         if (   dchIt->getProperty ("__sx_Key") == key
             && dchIt->getProperty ("__sx_Value").getType () == Type::Group) {
            validateObj (chIt, dchIt);
            nItemsFound++;
         }
      }

      // check needed groups
      SxString needStr = chIt->hasProperty ("needs")?
                         chIt->getProperty ("needs").getString()
                       : "";
      if (needStr != "" && nItemsFound > 0) {
         SxList<SxString> needList = needStr.tokenize (',');

         bool isPresent = true;
         for (auto nIt = needList.begin (); nIt != needList.end (); ++nIt) {
            if (!nameToIdx.containsKey (*nIt)) {
               isPresent = false;
            }
         }

         if (!isPresent) {
            validationError ("Presence of group '"
                             + key.getString ()
                             + "' requires also '" + needStr
                             + "'",
                             chIt->getProperty ("__sx_Value").getTag (),
                             dIt->getProperty ("__sx_Value").getTag ());
         }
      }


      // check not groups
      SxString notStr = chIt->hasProperty ("not")?
                        chIt->getProperty ("not").getString ()
                      : "";

      if (notStr != "" && nItemsFound > 0) {
         SxList<SxString> notList = notStr.tokenize (',');

         bool isPresent = false;
         for (auto notIt = notList.begin (); notIt != notList.end (); ++notIt)
         {
            if (nameToIdx.containsKey (*notIt))
               isPresent = true;
         }

         if (isPresent) {
            validationError ("Presence of group '"
                             + key.getString ()
                             + "' excludes '" + notStr
                             + "'",
                             chIt->getProperty ("__sx_Value").getTag (),
                             dIt->getProperty ("__sx_Value").getTag ());
         }
      }

      SxString xorStr = chIt->hasProperty ("xor")?
                        chIt->getProperty ("xor").getString ()
                      : "";

      if (xorStr != "") {
         SxList<SxString> xorList = xorStr.tokenize (',');

         bool xorFound = false;
         for (auto xorIt = xorList.begin (); xorIt != xorList.end (); ++xorIt)
         {
            if (   nameToIdx.containsKey (*xorIt)
                && *xorIt != key.getString ())
               xorFound = true;
         }

         if (nItemsFound == 0 && xorFound == false) {
            validationError ("Missing one of these groups '"
                             + xorStr + "'",
                             chIt->getProperty ("__sx_Value").getTag (),
                             dIt->getProperty ("__sx_Value").getTag ());
         }
         if (nItemsFound > 0 && xorFound == true) {
            validationError ("Only one of these groups '"
                             + xorStr + "' are allowed",
                             chIt->getProperty ("__sx_Value").getTag (),
                             dIt->getProperty ("__sx_Value").getTag ());
         }
      } else {
         SxString opt = chIt->hasProperty ("optional")?
                        chIt->getProperty ("optional").getString ()
                      : "";

         if (nItemsFound == 0 && (opt == "" || opt == "false")) {
            validationError ("group item '"
                             + key.getString ()
                             + "' not found in '"
                             + dIt->getProperty ("__sx_Key").getString ()
                             + "'",
                             chIt->getProperty ("__sx_Value").getTag (),
                             dIt->getProperty ("__sx_Value").getTag ());
         }

         ssize_t n = chIt->hasProperty ("nItems")?
                     chIt->getProperty ("nItems").getInt ()
                   : max<ssize_t>();
         ssize_t maxItems = chIt->hasProperty ("maxItems")?
                            chIt->getProperty ("maxItems").getInt ()
                          : max<ssize_t>();
         ssize_t minItems = chIt->hasProperty ("minItems")?
                            chIt->getProperty ("minItems").getInt ()
                          : min<ssize_t>();

         if (n != max<ssize_t>()) {
            if (nItemsFound > n) {
               validationError ("too many groups with key '"
                                + key.getString ()
                                + "'",
                                chIt->getProperty ("__sx_Value").getTag (),
                                dIt->getProperty ("__sx_Value").getTag ());
            }
            if (nItemsFound < n) {
               validationError ("too few groups with key '"
                                + key.getString ()
                                + "'",
                                chIt->getProperty ("__sx_Value").getTag (),
                                dIt->getProperty ("__sx_Value").getTag ());
            }
         }

         if (nItemsFound > maxItems) {
            validationError ("too many groups with key '"
                             + key.getString ()
                             + "'",
                             chIt->getProperty ("__sx_Value").getTag (),
                             dIt->getProperty ("__sx_Value").getTag ());
         }

         if (nItemsFound < minItems) {
            validationError ("too few groups with key '"
                             + key.getString ()
                             + "'",
                             chIt->getProperty ("__sx_Value").getTag (),
                             dIt->getProperty ("__sx_Value").getTag ());
         }
      }

   }

   // any unexpected items
   for (ssize_t m = 0; m < dIt.getSizeOut (); ++m) {
      bool found = false;
      SxString key = dIt.out (m)->getProperty ("__sx_Key").getString ();
      for (ssize_t n = 0; n < sIt.getSizeOut (); ++n) {
         SxString k = sIt.out (n)->getProperty ("__sx_Key").getString ();
         if (key == k) {
            found = true;
            break;
         }
      }
      if (!found) {
         validationError ("unexpected key '"
                          + key + "' in group '"
                          + sIt->getProperty ("__sx_Key").getString ()
                          + "'",
                          sIt->getProperty ("__sx_Value").getTag (),
                          dIt.out (m)->getProperty ("__sx_Value").getTag ());
      }
   }

   auto elemsIt = dIt->getProperties ().getKeys ().begin ();
   for (; elemsIt.isValid (); ++elemsIt) {
      bool found = false;
      const SxString &key = *elemsIt;
      if (key.find ("__sx_") == 0) continue;
      for (ssize_t n = 0; n < sIt.getSizeOut (); ++n) {
         SxString k = sIt.out (n)->getProperty ("__sx_Key").getString ();
         if (key == k) {
            found = true;
            break;
         }
      }
      if (!found) {
         validationError ("unexpected key '"
                          + key + "' in group '"
                          + sIt->getProperty ("__sx_Key").getString ()
                          + "'",
                          sIt->getProperty ("__sx_Value").getTag (),
                          dIt->getProperty (key).getTag ());
      }
   }

   return true;
}

bool SxSchema::validate (const SxPtr<SxGraph<SxGProps> > &dataG_)
{
   SX_CHECK (dataG_.getPtr ());

   dataG = dataG_;

   auto dIt = dataG->begin ();
   auto sIt = schemaG->begin ();

   // skip dummy root
   ++sIt; ++dIt;
   Type t = (Type)sIt->getProperty ("type").getInt ();
   if (t == Type::Group)
      validateObj (sIt, dIt);
   else
      validateArray (sIt, dIt);

   return true;
}

SxGraph<SxGProps>::Iterator
SxSchema::get (const SxGraph<SxGProps>::Iterator &it,
               const SxString &key,
               const SxVariant &val)
{
   SX_CHECK (key != "");
   SX_CHECK (it.isValid ());

   for (ssize_t i = 0; i < it.getSizeOut (); ++i)
   {
      auto chIt = it.out (i);
      if (chIt->hasProperty (key)) {
         if (!val.isInitialized ())
            return chIt;
         else if (chIt->getProperty (key) == val)
            return chIt;
      }
   }
   return schemaG->end ();
}

void SxSchema::simplifyInt (SxGraph<SxGProps>::Iterator &sIt)
{
   Type t = (Type)SxVariant::getTypeId (sIt->getProperty ("type").getString ());

   sIt->setProperty ("type", (int)t);

   SX_CHECK (t == Type::Int);

   ssize_t minVal = min<ssize_t>();
   ssize_t maxVal = max<ssize_t>();

   if (sIt->hasProperty ("min")) {
      minVal = sIt->getProperty ("min").getInt ();
   }
   if (sIt->hasProperty ("max")) {
      maxVal = sIt->getProperty ("max").getInt ();
   }

   sIt->getProperty ("__sx_Value").setType (Type::Int);
   sIt->getProperty ("__sx_Value").set (minVal);
   sIt->getProperty ("__sx_Value").setLimits (minVal, maxVal);
   SX_CHECK (sIt->getProperty ("__sx_Value").getType() == Type::Int);
}

void SxSchema::simplifyDouble (SxGraph<SxGProps>::Iterator &sIt)
{
   Type t = (Type)SxVariant::getTypeId (sIt->getProperty ("type").getString ());

   sIt->setProperty ("type", (int)t);

   SX_CHECK (t == Type::Double);

   double minVal = min<double>();
   double maxVal = max<double>();

   if (sIt->hasProperty ("min")) {
      minVal = sIt->getProperty ("min").getDouble ();
   }
   if (sIt->hasProperty ("max")) {
      maxVal = sIt->getProperty ("max").getDouble ();
   }

   sIt->getProperty ("__sx_Value").setType (Type::Double);
   sIt->getProperty ("__sx_Value").set (minVal);
   sIt->getProperty ("__sx_Value").setLimits (minVal, maxVal);
   SX_CHECK (sIt->getProperty ("__sx_Value").getType() == Type::Double);
}

void SxSchema::simplifyString (SxGraph<SxGProps>::Iterator &sIt)
{
   Type t = (Type)SxVariant::getTypeId (sIt->getProperty ("type").getString ());

   sIt->setProperty ("type", (int)t);

   SX_CHECK (t == Type::String);

   ssize_t minVal = min<ssize_t>();
   ssize_t maxVal = max<ssize_t>();

   if (sIt->hasProperty ("min")) {
      minVal = sIt->getProperty ("min").getInt ();
   }
   if (sIt->hasProperty ("max")) {
      maxVal = sIt->getProperty ("max").getInt ();
   }

   sIt->getProperty ("__sx_Value").setType (Type::String);
   sIt->getProperty ("__sx_Value").set (SxString (""));
   sIt->getProperty ("__sx_Value").setLimits (minVal, maxVal);
   SX_CHECK (sIt->getProperty ("__sx_Value").getType() == Type::String);
}

void SxSchema::simplifyBool (SxGraph<SxGProps>::Iterator &sIt)
{
   Type t = (Type)SxVariant::getTypeId (sIt->getProperty ("type").getString ());
   sIt->setProperty ("type", (int)t);

   SX_CHECK (t == Type::Bool);
}


void SxSchema::simplifyArray (SxGraph<SxGProps>::Iterator &sIt)
{
   Type t = (Type)SxVariant::getTypeId (sIt->getProperty ("type").getString ());

   sIt->setProperty ("type", (int)t);

   SX_CHECK (t == Type::List);

   auto itemsTypeIt = get (sIt, "__sx_Key", "itemsType");
   SX_CHECK (itemsTypeIt.isValid ());

   Type itemType = (Type)SxVariant::getTypeId (itemsTypeIt->getProperty ("type").getString ());

   if (itemType == Type::Group) {
      simplifyObj (itemsTypeIt);
   } else if (itemType == Type::List) {
      itemsTypeIt->setProperty ("type", (int)itemType);
      auto subArrayIt = get (sIt, "__sx_Key", "subarray");
      simplifyArray (subArrayIt);
   } else if (itemType == Type::Int) {
      simplifyInt (itemsTypeIt);
   } else if (itemType == Type::Double) {
      simplifyDouble (itemsTypeIt);
   } else if (itemType == Type::Bool) {
      simplifyBool (itemsTypeIt);
   } else if (itemType == Type::String) {
      simplifyString (itemsTypeIt);
   }

}

void SxSchema::simplifyObj (SxGraph<SxGProps>::Iterator &sIt)
{
   Type t = (Type)SxVariant::getTypeId (sIt->getProperty ("type").getString ());

   sIt->setProperty ("type", (int)t);

   SX_CHECK (t == Type::Group);

   for (ssize_t i = 0; i < sIt.getSizeOut (); ++i) {
      auto chIt = sIt.out (i);
      Type elemType = (Type)SxVariant::getTypeId (chIt->getProperty ("type").getString ());

      if (elemType == Type::Group) {
         simplifyObj (chIt);
      } else if (elemType == Type::List) {
         simplifyArray (chIt);
      } else if (elemType == Type::Int) {
         simplifyInt (chIt);
      } else if (elemType == Type::Double) {
         simplifyDouble (chIt);
      } else if (elemType == Type::Bool) {
         simplifyBool (chIt);
      } else if (elemType == Type::String) {
         simplifyString (chIt);
      }
   }

}

void SxSchema::simplifySchemaAst ()
{
   SX_CHECK (schemaG.getPtr ());

   auto sIt = schemaG->begin ();
   // skip dummy root
   ++sIt;
   Type t = (Type)SxVariant::getTypeId (sIt->getProperty ("type").getString ());
   if (t == Type::Group)
      simplifyObj (sIt);
   else
      simplifyArray (sIt);
}

