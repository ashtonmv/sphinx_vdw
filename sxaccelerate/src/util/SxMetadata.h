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

#ifndef _SX_METADATA_H_
#define _SX_METADATA_H_

#include <SxUtil.h>
#include <SxMap.h>
#include <SxString.h>
#include <typeinfo>

/*** \brief Auxiliary class for SxMetadata: the virtual base class
  */
class SxMetaBase {
   public:
      /// Virtual destructor
      virtual ~SxMetaBase () {/* empty */}
      /// Explicit type cast operator (crashes if type does not match)
      template <class T> inline operator T& ();
      /// Explicit type cast operator (crashes if type does not match)
      template <class T> inline operator const T& ()
      {
         return operator T& ();
      }
      /// Explicit type cast operator (crashes if type does not match)
      template <class T> inline operator const T& () const;
   protected:
      template <class T> friend class SxMetaMap;
      friend class SxMetadata;
      /// Get a copy (for replicating metadata)
      virtual SxMetaBase * getCopy () const = 0;
      const std::type_info *typeId;
      /// Check type compatibility
      bool typeMatch (const std::type_info &id) const { return *typeId == id; }
      /// Size of relevant type
      size_t typeSize;
};

/// Data containing template class
template<class T>
class SxMetaItem : public SxMetaBase
{
   public:
      /// The actual value
      T value;
      /// Constructor
      SxMetaItem (const T& in) : value(in) { /* empty */}
      /// Destructor
      virtual ~SxMetaItem () { /* empty */ }
   protected:
      /// Get a copy (for replicating metadata)
      virtual SxMetaBase* getCopy () const
      {
         return new SxMetaItem<T> (value);
      }
};

namespace {
template<class T> inline size_t mysizeof_ (T *) { return sizeof(T); }

/// template specialization for mysizeof: SxPtr
template<class T> inline size_t mysizeof_ (SxPtr<T> *) { return sizeof(T); }

/// template specialization for mysizeof: C pointers
template<class T> inline size_t mysizeof_(T **) { return sizeof(T); }

template <class T> inline size_t mysizeof () { return mysizeof_((T*)(NULL)); }
template<> inline size_t mysizeof<void> () { return 0; }

}

/*** \brief Auxiliary class for SxMetadata: mapping a key to a value

  @note The key type is the template parameter
  */
template<class Key>
class SxMetaMap : public SxMetaBase
{
   protected:
      /// The metadata container
      SxMap<Key, SxMetaBase*> map;
   public:
      /// Check for existence
      bool contains (const Key &key) const
      {
         return map.containsKey (key);
      }

      /// Type-aware get
      template<class T>
      T& get (const Key &key) const
      {
         SX_CHECK (map.containsKey (key));
         return map(key)->operator T& ();
      }

      /// Type-aware get pointer, return NULL if not available
      template<class T>
      const T* getPtr (const Key &key) const
      {
         if (map.containsKey (key))
            return &map(key)->operator T& ();
         return NULL;
      }

      /// Type-blind get (to be resolved by assignment / type cast)
      SxMetaBase& get (const Key &key) const {
         SX_CHECK (map.containsKey (key));
         return *(map(key));
      }

      /// Attach new data (key must not exist)
      template<class T>
      void attach (const Key &key, const T &value)
      {
         SX_CHECK (!map.containsKey (key));
         map(key) = new SxMetaItem<T> (value); // deallocated in destructor
         //std::cout << __PRETTY_FUNCTION__ << ": " << value << std::endl;
      }

      /// Update metadata if key exists, or attach new
      template<class T>
      void update (const Key &key, const T &value)
      {
         if (map.containsKey (key))  {
            SxMetaItem<T> *ptr = dynamic_cast<SxMetaItem<T>*>(map(key));
            SX_CHECK(ptr);
            ptr->value = value;
         } else {
            attach (key, value);
         }
      }

      /// Remove some metadata
      void remove (const Key &key)
      {
         if (map.containsKey (key))  {
            delete map(key);
            map.removeKey (key);
         }
      }

      /// Remove all metadata
      void removeAll ()
      {
         for (typename SxMap<Key, SxMetaBase*>::Iterator it = map.begin ();
              it != map.end ();
              it++)
         {
            delete it.getValue ();
         }
         map.removeAll ();
      }

      /// Destructor
      virtual ~SxMetaMap () { removeAll (); }

      SxMetaMap () {
         typeId = &typeid(Key);
         typeSize = mysizeof_((Key*)(NULL));
      }

   protected:
      friend class SxMetadata;
      /// Replicate data
      void replicate (const SxMetaMap<Key> &in)
      {
         removeAll ();
         for (typename SxMap<Key, SxMetaBase*>::ConstIterator it = in.map.begin ();
              it != in.map.end ();
              ++it)
         {
            map(it.getKey ()) = it.getValue ()->getCopy ();
         }
      }

      /// Get a copy of this metadata container
      virtual SxMetaBase* getCopy () const
      {
         SxMetaMap<Key> *newMap = new SxMetaMap<Key> ();
         newMap->replicate (*this);
         return newMap;
      }
};


/*** \brief Auxiliary class for SxMetadata: type is the key
  */
template<>
class SxMetaMap<void> : public SxMetaBase
{
   protected:
      /// The metadata container
      SxMap<size_t, SxList<SxMetaBase*> > map;

      /// Find the metadata
      template <class T>
      SxMetaItem<T>* find (SxList<SxMetaBase*>*& list) const
      {
         if (!map.containsKey (mysizeof<T> () )) {
            list = NULL;
            return NULL;
         }
         list = const_cast<SxList<SxMetaBase*>*> (&map(mysizeof<T> ()));
         for (SxList<SxMetaBase*>::Iterator it = list->begin ();
              it != list->end ();
              it++)
         {
            if (SxMetaItem<T>* ptr = dynamic_cast<SxMetaItem<T>*> (*it))  {
               return ptr;
            }
         }
         return NULL;
      }

      /// Wrapper for find
      template <class T> SxMetaItem<T>* find () const
      {
         SxList<SxMetaBase*> *dummy;
         return find<T> (dummy);
      }

   public:
      /// Check for metadata presence
      template<class T> bool contains () const
      {
         return find<T> ();
      }


      /// Get metadata of specified type (the type is the key!)
      template<class T> T& get () const
      {
         if (SxMetaItem<T>* ptr = find<T> ())  {
            return ptr->operator T& ();
         }
         std::cout << "Could not find metadata in " << __PRETTY_FUNCTION__ << std::endl;
         SX_EXIT;
      }

      /** \brief Get pointer to metadata of specified type (the type is the
                 key!), or NULL
      */
      template<class T> T* getPtr () const
      {
         if (SxMetaItem<T>* ptr = find<T> ())
            return &ptr->operator T& ();
         return NULL;
      }

      /// Attach new data (must not exist before)
      template<class T> void attach (const T &value)
      {
         SxList<SxMetaBase*> *list = NULL;
         if (find<T> (list)) { SX_EXIT; }
         if (!list) list = &map(mysizeof<T> ());
         *list << new SxMetaItem<T> (value);
         //std::cout << __PRETTY_FUNCTION__ << ": " << value << std::endl;
      }

      /// Update metadata if it exists, or attach new
      template<class T> void update (const T &value)
      {
         SxList<SxMetaBase*> *list;
         if (SxMetaItem<T>* ptr = find<T> (list))  {
            ptr->value = value;
         } else {
            if (!list) list = &map(mysizeof<T> ());
            *list << new SxMetaItem<T> (value);
         }
      }

      /// Remove metadata
      template<class T> void remove ()
      {
         if (!map.containsKey (mysizeof<T> () )) return;
         SxList<SxMetaBase*> &list = map(mysizeof<T> ());
         for (SxList<SxMetaBase*>::Iterator it = list.begin ();
              it != list.end ();
              it++)
         {
            if (SxMetaItem<T>* ptr = dynamic_cast<SxMetaItem<T>*> (*it))  {
               delete ptr;
               list.removeItem (&it);
               if (list.getSize () == 0) map.removeKey (mysizeof<T> ());
               return;
            }
         }
      }

      /// Remove all metadata
      void removeAll ()
      {
         for (SxMap<size_t, SxList<SxMetaBase*> >::Iterator listIt = map.begin ();
              listIt != map.end ();
              listIt++)
         {
            for (SxList<SxMetaBase*>::Iterator it = listIt.getValue ().begin ();
                 it != listIt.getValue ().end ();
                 ++it)
            {
               delete *it;
            }
         }
         map.removeAll ();
      }

      SxMetaMap () {
         typeId = &typeid(void);
         typeSize = 0;
      }

      /// Destructor
      virtual ~SxMetaMap () { removeAll (); }


   protected:
      friend class SxMetadata;
      /// Get copy of this metadata
      virtual SxMetaBase* getCopy () const
      {
         SxMetaMap<void> *newMap = new SxMetaMap<void> ();
         for (SxMap<size_t, SxList<SxMetaBase*> >::ConstIterator listIt = map.begin ();
              listIt != map.end ();
              listIt++)
         {
            for (SxList<SxMetaBase*>::ConstIterator it = listIt.getValue ().begin ();
                 it != listIt.getValue ().end ();
                 ++it)
            {
               newMap->map(listIt.getKey ()).append ((*it)->getCopy ());
            }
         }

         return newMap;
      }
};


/** \brief Flexible metadata

    \b SxClass = S/PHI/nX metadata class

    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_UTIL SxMetadata
{
   protected:
      /// List of possible key->value maps. Each map is specific
      /// for the key type
      mutable SxArray<SxMetaBase *> metaMaps;

      /// Auxiliary function: quickly find map for a key type from typeSize and type_info
      SxMetaBase* findMap (size_t typeSize, const std::type_info &typeId) const;

      template <class Key>
      SxMetaMap<Key>& getMap () const
      {
         static const std::type_info &typeId = typeid(Key);
         SxMetaBase *ptr = findMap (mysizeof<Key> (), typeId);
         if (ptr)  {
            SX_CHECK (dynamic_cast<SxMetaMap<Key>*> (ptr));
            return *static_cast<SxMetaMap<Key>*>(ptr);
         }
         SxMetaMap<Key> *newMap = new SxMetaMap<Key> ();
         metaMaps(metaMaps.getSize () - 1) = newMap;
         return *newMap;
      }

   public:
      /// Check for metadata presence
      template <class Key>
      bool contains(const Key &key) const {
         return getMap<Key> ().contains (key);
      }

      /// Get metadata with explicit type
      template <class Value, class Key>
      const Value& get(const Key &key) const {
         return getMap<Key> ().template get<Value> (key);
      }

      /// Get pointer to metadata with explicit type, or NULL if not available
      template <class Value, class Key>
      const Value* getPtr(const Key &key) const {
         return getMap<Key> ().template getPtr<Value> (key);
      }

      /// Get metadata with explicit type
      template <class Value, class Key>
      Value& get(const Key &key)  {
         return getMap<Key> ().template get<Value> (key);
      }

      /// Get metadata with implicit type (e.g. reference assignment)
      template <class Key>
      SxMetaBase &get (const Key &key)  {
         return getMap<Key> ().get (key);
      }

      /// Attach new metadata
      template <class Value, class Key>
      void attach (const Key &key, const Value &value) {
         getMap<Key> ().attach (key, value);
      }

      /// Update or attach new metadata
      template <class Value, class Key>
      void update (const Key &key, const Value &value) {
         getMap<Key> ().update (key, value);
      }

      /// Remove metadata
      template <class Key>
      void remove (const Key &key) {
         getMap<Key> ().remove (key);
      }


      // --- const char* => SxString
      /// Check for metadata presence
      bool contains(const char *key) const {
         return getMap<SxString> ().contains (key);
      }

      /// Get metadata with explicit type
      template<class Value>
      const Value& get(const char *key) const {
         return getMap<SxString> ().get<Value> (key);
      }

      /// Get pointer to metadata with explicit type, or NULL if not available
      template <class Value>
      const Value* getPtr(const char *key) const {
         return getMap<SxString> ().getPtr<Value> (key);
      }

      /// Get metadata with explicit type
      template<class Value>
      Value& get(const char *key)  {
         return getMap<SxString> ().get<Value> (key);
      }

      /// Get metadata with implicit type (e.g. reference assignment)
      SxMetaBase &get (const char *key)  {
         return getMap<SxString> ().get (key);
      }

      /// Attach new metadata
      template<class Value>
      void attach (const char *key, const Value &value) {
         getMap<SxString> ().attach (key, value);
      }

      /// Update or attach new metadata
      template<class Value>
      void update (const char *key, const Value &value) {
         getMap<SxString> ().update (key, value);
      }

      /// Remove metadata
      void remove (const char *key) {
         getMap<SxString> ().remove (key);
      }
      // --- Now the special case where type is used as key

      /// Check for presence of metadata (type is key)
      template <class Value> bool contains () const {
         return getMap<void> ().contains<Value> ();
      }

      /// Get metadata (type is key)
      template <class Value> const Value& get () const {
         return getMap<void> ().get<Value> ();
      }

      /// Get pointer to metadata with explicit type, or NULL if not available
      template <class Value>
      const Value* getPtr () const {
         return getMap<void> ().getPtr<Value> ();
      }

      /// Get metadata (type is key)
      template <class Value> Value& get()  {
         return getMap<void> ().get<Value> ();
      }

      /// Attach new metadata (type is key)
      template <class Value> void attach (const Value &value) {
         getMap<void> ().attach (value);
      }

      /// Update or attach new metadata (type is key)
      template <class Value> void update (const Value &value) {
         getMap<void> ().update (value);
      }

      /// Remove metadata (type is key)
      template <class Value> void remove () {
         getMap<void> ().remove<Value> ();
      }
      // -------------------------

      /// Remove all metadata
      void removeAll ();

      /// Destructor
      virtual ~SxMetadata () { removeAll (); }

      /// Constructor
      SxMetadata () {/* empty */}

      /// Copy assignment
      SxMetadata& operator= (const SxMetadata &in);

      /// Copy constructor
      SxMetadata (const SxMetadata &in) { (*this) = in; }
};

// --- implementation ---

template <class T>
SxMetaBase::operator T& ()
{
   //SxMetaItem<T> *ptr = dynamic_cast<SxMetaItem<T> *>(this);
   // Check for type mismatch
   //SX_CHECK (ptr);
   SX_CHECK (dynamic_cast<SxMetaItem<T> *>(this));
   SxMetaItem<T> *ptr = static_cast<SxMetaItem<T> *>(this);
   return ptr->value;
}

template <class T>
SxMetaBase::operator const T& () const
{
   //SxMetaItem<T> *ptr = dynamic_cast<SxMetaItem<T> *>(this);
   // Check for type mismatch
   //SX_CHECK (ptr);
   SX_CHECK (dynamic_cast<SxMetaItem<T> *>(this));
   SxMetaItem<T> *ptr = static_cast<SxMetaItem<T> *>(this);
   return ptr->value;
}

#endif /* _SX_METADATA_H_ */
