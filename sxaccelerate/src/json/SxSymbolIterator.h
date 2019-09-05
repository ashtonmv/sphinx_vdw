#ifndef _SX_SYMBOL_ITERATOR_H_
#define _SX_SYMBOL_ITERATOR_H_

#include <SxParserAst.h>
#include <SxGProps.h>
#include <SxGraph.h>
#include <SxGQuery.h>

using namespace sx;
typedef typename SxGQExprBase::SelSet SelSet;
typedef typename SxVariantType::DataType  Type;

class SxSymbolIterator
{
   public:
      SxSymbolIterator () { }
      SxSymbolIterator (const SxPtr<SxGraph<SxGProps> > &gPtr_, ssize_t nodeId_)
      {
         SX_CHECK (gPtr_.getPtr ());
         SX_CHECK (nodeId_ >= 0, nodeId_);

         dataG   = gPtr_;
         value   = NULL;
         nodeId  = nodeId_;
         childId = -1;

         auto it = dataG->begin (nodeId);
         SX_CHECK (it.isValid ());

         name  = it->getProperty ("sxKey").toString ();
         value = &(it->getProperty ("sxValue"));
         type  = (Type)value->getType ();
      }

      SxSymbolIterator (const SxPtr<SxGraph<SxGProps> > &gPtr_, ssize_t nodeId_,
                        const SxString &key_)
      {
         SX_CHECK (gPtr_.getPtr ());
         SX_CHECK (nodeId_ >= 0, nodeId_);

         dataG   = gPtr_;
         value   = NULL;
         nodeId  = nodeId_;
         childId = -1;

         auto it = dataG->begin (nodeId);
         SX_CHECK (it.isValid ());
         SX_CHECK (it->hasProperty (key_));
         name  = key_;
         value = &(it->getProperty (name));
         type  = (Type)value->getType ();
      }

      SxSymbolIterator (const SxPtr<SxGraph<SxGProps> > &gPtr_, ssize_t nodeId_,
                        SxVariant *val_)
      {
         SX_CHECK (gPtr_.getPtr ());
         SX_CHECK (nodeId_ >= 0, nodeId_);

         dataG   = gPtr_;
         value   = NULL;
         nodeId  = nodeId_;
         childId = -1;

         auto it = dataG->begin (nodeId);
         SX_CHECK (it.isValid ());
         name  = "undefined";
         value = val_;
         type  = (Type)value->getType ();
      }

     ~SxSymbolIterator () { }

      void setChildId (ssize_t id_) {
         SX_CHECK (id_ >=0);
         SX_CHECK (id_ < dataG->begin (nodeId).in(0).getSizeOut ());
         childId = id_;
      }

      ssize_t getChildId () const {
         return childId;
      }

      SxSymbolIterator getElem (const SxString &key) const {
         auto it = dataG->begin (nodeId);
         SX_CHECK (it.isValid ());

         // element is of primitive type, use specific functions
         SX_CHECK (!it->hasProperty (key));

         for (ssize_t i = 0; i < it.getSizeOut (); ++i) {
            auto chIt = it.out (i);
            if (chIt->getProperty ("sxKey").toString () == key) {
               SxSymbolIterator s(dataG, chIt.getIdx ());
               s.setChildId (i);
               return s;
            }
         }
         return SxSymbolIterator ();
      }

      SxList<ssize_t> toIntList () const {
         SX_CHECK (type == Type::List);
         SX_CHECK (value != NULL);

         SxList<SxVariant>::Iterator it = value->begin ();
         SxList<ssize_t> lst;
         for (; it != value->end (); ++it) {
            lst.append (it->toInt ());
         }
         return lst;
      }

      SxList<SxString> toStringList () const {
         SX_CHECK (type == Type::List);
         SX_CHECK (value != NULL);

         SxList<SxVariant>::Iterator it = value->begin ();
         SxList<SxString> lst;
         for (; it != value->end (); ++it) {
            lst.append (it->toString ());
         }
         return lst;
      }

      SxList<double> toDoubleList () const {
         SX_CHECK (type == Type::List);
         SX_CHECK (value != NULL);

         SxList<SxVariant>::Iterator it = value->begin ();
         SxList<double> lst;
         for (; it != value->end (); ++it) {
            lst.append (it->toDouble ());
         }
         return lst;
      }

      SxSymbolIterator getNext () const {
         SX_CHECK (childId >=0);
         SX_CHECK (dataG->begin (nodeId).getSizeIn () > 0);
         auto parent = dataG->begin (nodeId).in (0);
         if (childId+1 < parent.getSizeOut ()) {
            SxSymbolIterator res (dataG, parent.out ((childId+1)).getIdx ());
            res.setChildId ((childId+1));
            return res;
         }
         return SxSymbolIterator ();
      }

      SxList<SxSymbolIterator> toList () const {
         SX_CHECK (type == Type::Group
                || type == Type::List, type);

         SxList<SxSymbolIterator> res;
         if (type == Type::Group) {
            auto nIt = dataG->begin (nodeId);
            auto kIt = nIt->getProperties ().getKeys ().begin ();
            for (;kIt.isValid (); ++kIt) {
               if ((*kIt)(0,1) != "sx") {
                  res.append (SxSymbolIterator (dataG, nodeId, *kIt));
               }
            }
         }

         SxGQuery q = (N("sxKey") == name) ^ (N("sxKey").any ());
         SelSet sels = q.matchAll (dataG, dataG->begin (nodeId).hops (0));

         for (auto it = sels->begin ();it != sels->end (); ++it) {
            auto sel = *it;
            res.append (SxSymbolIterator(dataG, (*sel)(1)));
         }

         if (type == Type::List) {
            auto elemIt = value->begin ();
            for (; elemIt.isValid (); ++elemIt) {
               res.append (SxSymbolIterator (dataG, nodeId, &(*elemIt)));
            }
         }

         return res;
      }

      SxString getName () const {
         return name;
      }

      Type getType () const {
         return type;
      }

      SxString toString () const {
         SX_CHECK (type == Type::String);
         SX_CHECK (value != NULL);
         return value->getString ();
      }

      ssize_t toInt () const {
         SX_CHECK (type == Type::Int);
         SX_CHECK (value != NULL);
         return value->getInt ();
      }

      double toDouble () const {
         SX_CHECK (type == Type::Double);
         SX_CHECK (value != NULL);
         return value->getDouble ();
      }

      bool toBool () const {
         SX_CHECK (type == Type::Bool);
         SX_CHECK (value != NULL);
         return value->getBool ();
      }

      bool isValid () const {
         if (type != Type::Undefined) return true;
         else                         return false;
      }

      void printObject (int lvl, int tabs = 3) const {
         SX_CHECK (type == Type::Group);
         SxList<SxSymbolIterator> elems = toList ();
         ssize_t c = 0;
         for (int i = 0; i < lvl*tabs; ++i)
            cout << " ";

         if (name != "")
            cout << name << ":";

         cout << "{";
         for (auto it = elems.begin (); it != elems.end (); ++it)
         {
            if (c != 0)
               cout << ",";
            cout << "\n";
            it->print (lvl+1);
            c++;
         }
         cout << "\n";
         for (int i = 0; i < lvl*tabs; ++i)
            cout << " ";
         cout << "}";
      }

      void printArray (int lvl, int tabs = 3) const {
         SX_CHECK (type == Type::List);
         SxList<SxSymbolIterator> elems = toList ();
         ssize_t c = 0;
         for (int i = 0; i < lvl*tabs; ++i)
            cout << " ";

         if (name != "")
            cout << name << ":";

         cout << "[";
         for (auto it = elems.begin (); it != elems.end (); ++it)
         {
            if (c != 0)
               cout << ",";
            cout << "\n";
            it->print (lvl+1);
            c++;
         }
         cout << "\n";
         for (int i = 0; i < lvl*tabs; ++i)
            cout << " ";
         cout << "]";

      }

      void printInt (int lvl, int tabs = 3) const {
         SX_CHECK (type == Type::Int);

         for (int i = 0; i < lvl*tabs; ++i)
            cout << " ";

         if (name != "")
            cout << name << ":";

         cout << toInt ();
      }

      void printDouble (int lvl, int tabs = 3) const {
         SX_CHECK (type == Type::Double);

         for (int i = 0; i < lvl*tabs; ++i)
            cout << " ";

         if (name != "")
            cout << name << ":";

         cout << toDouble ();
      }

      void printString (int lvl, int tabs = 3) const {
         SX_CHECK (type == Type::String);

         for (int i = 0; i < lvl*tabs; ++i)
            cout << " ";

         if (name != "")
            cout << name << ":";

         cout << toString ();
      }

      void printBool (int lvl, int tabs = 3) const {
         SX_CHECK (type == Type::Bool);

         for (int i = 0; i < lvl*tabs; ++i)
            cout << " ";

         if (name != "")
            cout << name << ":";
         if (toBool ())
            cout << "true";
         else
            cout << "false";
      }

      void print (int lvl) const {
         SX_CHECK (isValid ());
         if (type == Type::Group) {
            printObject (lvl);
         } else if (type == Type::List) {
            printArray (lvl);
         } else if (type == Type::Int) {
            printInt (lvl);
         } else if (type == Type::Double) {
            printDouble (lvl);
         } else if (type == Type::String) {
            printString (lvl);
         } else if (type == Type::Bool) {
            printBool (lvl);
         } else {
            SX_EXIT; // type unknown
         }
      }

   protected:
      SxString name; // key from json element
      Type type; // type of value
      SxPtr<SxGraph<SxGProps> > dataG;
      SxVariant *value;
      ssize_t nodeId;
      ssize_t childId;
};

#endif /* _SX_SYMBOL_ITERATOR_H_ */
