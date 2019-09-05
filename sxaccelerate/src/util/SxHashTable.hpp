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

template<class K, class V, class H>
SxHashTable<K,V,H>::SxHashTable ()
   : tableSize(0),
     tableUsed(0),
     tableBound(0)
{
   // empty
}

template<class K, class V, class H>
void SxHashTable<K,V,H>::remove (ssize_t idx_)
{
   this->elements[idx_] = &deleted;
   tableSize--;
}

template<class K, class V, class H>
void SxHashTable<K,V,H>::removeAll ()
{
   SxArray<typename SxList<K>::Node*>::resize (0);
   values.resize (0);
   tableSize = 0;
   tableUsed = 0;
   tableBound = 0;
}

template<class K, class V, class H>
ssize_t SxHashTable<K,V,H>::findPos (const K &elem_) const
{
   if (tableSize > 0)  {
      size_t h = H::hash (elem_);
      size_t m = static_cast<size_t>(this->getSize());
      size_t idx = h % m;

      if (this->elements[idx] &&
          this->elements[idx] != &deleted &&
          this->elements[idx]->elem == elem_)
      {
         return static_cast<ssize_t>(idx);
      }

      size_t step = h % (m - 1) + 1;
      size_t lastIdx = idx;

      idx = (idx + step) % m;
      while (idx != lastIdx) {
         if (this->elements[idx] == NULL)  {
            return -1;
         }  else if (this->elements[idx] != &deleted &&
                     this->elements[idx]->elem == elem_)
         {
            return static_cast<ssize_t>(idx);
         }
         idx = (idx + step) % m;
      }
   }

   return -1;
}

template<class K, class V, class H>
bool SxHashTable<K,V,H>::contains (const K &elem_) const
{
   return this->findPos (elem_) >= 0;
}

template<class K, class V, class H>
bool SxHashTable<K,V,H>::contains (const K &elem_, ssize_t *result_)
{
   SX_CHECK (result_);

   if (tableUsed >= tableBound)  {
      if (tableSize * 2 < this->getSize())  {
         resize (this->getSize() - 1);
      }  else  {
         resize (this->getSize() + 1);
      }
   }

   size_t h = H::hash (elem_);
   size_t m = static_cast<size_t>(this->getSize());
   size_t idx = h % m;
   size_t step = h % (m - 1) + 1;
   size_t lastIdx = idx;
   size_t replaceIdx = m;

   do {
      if (this->elements[idx] == NULL)  {
         // --- insert new element
         if (replaceIdx != m)  {
            idx = replaceIdx;
         }  else  {
            tableUsed++;
         }
         tableSize++;
         *result_ = static_cast<ssize_t>(idx);
         return false;
      }  else if (this->elements[idx] == &deleted)  {
         replaceIdx = idx;
      }  else if (this->elements[idx]->elem == elem_)  {
         *result_ = static_cast<ssize_t>(idx);
         return true;
      }
      idx = (idx + step) % m;
   } while (idx != lastIdx);

   // "Can't insert to HashTable."
   SX_CHECK (idx != lastIdx);
   return false;
}

template<class K, class V, class H>
void SxHashTable<K,V,H>::resize (ssize_t newSize_)
{
   // "HashTable tableSize is too big."
   SX_CHECK (newSize_ < 1610612741, newSize_);

   // --- http://planetmath.org/encyclopedia/GoodHashTablePrimes.html
   ssize_t prime[30] = {0, 3, 11, 23, 53, 97, 193, 389, 769, 1543, 3079, 6151,
                        12289, 24593, 49157, 98317, 196613, 393241, 786433,
                        1572869, 3145739, 6291469, 12582917, 25165843, 50331653,
                        100663319, 201326611, 402653189, 805306457, 1610612741};
   size_t idx = 0;
   while (prime[idx] <= newSize_)  {
      idx++;
   }
   newSize_ = prime[idx];
   ssize_t newBound = (ssize_t)((double)newSize_ * 0.77 + 0.5);
   ssize_t n = this->getSize ();
   ssize_t i;

   if (tableSize < newBound)  {
      SxArray<typename SxList<K>::Node*> newTable(newSize_);
      SxArray<typename SxList<V>::Node*> newValues(newSize_);
      if (newSize_ > 0)  {
         newTable.set (NULL);
         newValues.set (NULL);
      }

      for (i = 0; i < n; i++)  {
         if (this->elements[i] && this->elements[i] != &deleted)  {
            size_t h = H::hash (this->elements[i]->elem);
            size_t m = static_cast<size_t>(newSize_);
            idx = h % m;
            size_t step = h % (m - 1) + 1;
#ifndef NDEBUG
            size_t lastIdx = idx;
#endif

            while (newTable.elements[idx])  {
               idx = (idx + step) % m;
               SX_CHECK (idx != lastIdx);
            }
            newTable.elements[idx] = this->elements[i];
            newValues.elements[idx] = values.elements[i];
         }
      }

      tableUsed = tableSize;
      tableBound = newBound;
      SxArray<typename SxList<K>::Node*>::resize (newSize_, false);
      size_t len = sizeof(typename SxList<K>::Node*) * (size_t)this->getSize();
      ::memcpy (this->elements, newTable.elements, len);

      values.resize (newSize_, false);
      len = sizeof(typename SxList<V>::Node*) * (size_t)this->getSize();
      ::memcpy (values.elements, newValues.elements, len);
   }
}

template<class K, class V, class H>
void SxHashTable<K,V,H>::print () const
{
   ssize_t i;
   ssize_t n = this->getSize ();

   std::cout << "table[" << tableSize << "/" << n << "]: ";

   for (i = 0; i < n; i++)  {
      if (this->elements[i] == NULL)  {
         std::cout << i << ":NULL ";
      }  else if (this->elements[i] == &deleted)  {
         std::cout << i << ":DEL ";
      }  else  {
         std::cout << i << ":" << (void*)this->elements[i] << " ";
      }
   }
   std::cout << std::endl;
}

// ----------------------------------------------------------------------------

template<class K, class H>
SxHashTable<K,SxNull,H>::SxHashTable ()
   : tableSize(0),
     tableUsed(0),
     tableBound(0)
{
   // empty
}

template<class K, class H>
void SxHashTable<K,SxNull,H>::remove (ssize_t idx_)
{
   this->elements[idx_] = &deleted;
   tableSize--;
}

template<class K, class H>
void SxHashTable<K,SxNull,H>::removeAll ()
{
   SxArray<typename SxList<K>::Node*>::resize (0);
   tableSize = 0;
   tableUsed = 0;
   tableBound = 0;
}

template<class K, class H>
ssize_t SxHashTable<K,SxNull,H>::findPos (const K &elem_) const
{
   if (tableSize > 0)  {
      size_t h = H::hash (elem_);
      size_t m = static_cast<size_t>(this->getSize());
      size_t idx = h % m;
      size_t step = h % (m - 1) + 1;
      size_t lastIdx = idx;

      do {
         if (this->elements[idx] == NULL)  {
            return -1;
         }  else if (this->elements[idx] != &deleted &&
                     this->elements[idx]->elem == elem_)
         {
            return static_cast<ssize_t>(idx);
         }
         idx = (idx + step) % m;
      } while (idx != lastIdx);
   }
   return -1;
}

template<class K, class H>
bool SxHashTable<K,SxNull,H>::contains (const K &elem_) const
{
   return (this->findPos (elem_) >= 0) ? true : false;
}

template<class K, class H>
bool SxHashTable<K,SxNull,H>::contains (const K &elem_, ssize_t *result_)
{
   SX_CHECK (result_);

   if (tableUsed >= tableBound)  {
      if (tableSize * 2 < this->getSize())  {
         resize (this->getSize() - 1);
      }  else  {
         resize (this->getSize() + 1);
      }
   }
   SX_CHECK (this->getSize() > 0);

   size_t h = H::hash (elem_);
   size_t m = static_cast<size_t>(this->getSize());
   size_t idx = h % m;
   size_t step = h % (m - 1) + 1;
   size_t lastIdx = idx;
   size_t replaceIdx = m;

   do {
      if (this->elements[idx] == NULL)  {
         // --- insert new element
         if (replaceIdx != m)  {
            idx = replaceIdx;
         }  else  {
            tableUsed++;
         }
         tableSize++;
         *result_ = static_cast<ssize_t>(idx);
         return false;
      }  else if (this->elements[idx] == &deleted)  {
         replaceIdx = idx;
      }  else if (this->elements[idx]->elem == elem_)  {
         *result_ = static_cast<ssize_t>(idx);
         return true;
      }
      idx = (idx + step) % m;
   } while (idx != lastIdx);

   // "Can't insert to HashTable."
   SX_CHECK (idx != lastIdx);
   return false;
}

template<class K, class H>
void SxHashTable<K,SxNull,H>::resize (ssize_t newSize_)
{
   // "HashTable tableSize is too big."
   SX_CHECK (newSize_ < 1610612741, newSize_);

   // --- http://planetmath.org/encyclopedia/GoodHashTablePrimes.html
   ssize_t prime[30] = {0, 3, 11, 23, 53, 97, 193, 389, 769, 1543, 3079, 6151,
                        12289, 24593, 49157, 98317, 196613, 393241, 786433,
                        1572869, 3145739, 6291469, 12582917, 25165843, 50331653,
                        100663319, 201326611, 402653189, 805306457, 1610612741};
   size_t idx = 0;
   while (prime[idx] <= newSize_)  {
      idx++;
   }
   newSize_ = prime[idx];
   ssize_t newBound = (ssize_t)((double)newSize_ * 0.77 + 0.5);
   ssize_t n = this->getSize ();
   ssize_t i;

   if (tableSize < newBound)  {
      SxArray<typename SxList<K>::Node*> newTable(newSize_);
      if (newSize_ > 0)  {
         newTable.set (NULL);
      }

      for (i = 0; i < n; i++)  {
         if (this->elements[i] && this->elements[i] != &deleted)  {
            size_t h = H::hash (this->elements[i]->elem);
            size_t m = static_cast<size_t>(newSize_);
            idx = h % m;
            size_t step = h % (m - 1) + 1;
#ifndef NDEBUG
            size_t lastIdx = idx;
#endif

            while (newTable.elements[idx])  {
               idx = (idx + step) % m;
               SX_CHECK (idx != lastIdx);
            }
            newTable.elements[idx] = this->elements[i];
         }
      }

      tableUsed = tableSize;
      tableBound = newBound;
      SxArray<typename SxList<K>::Node*>::resize (newSize_, false);
      size_t len = sizeof(typename SxList<K>::Node*) * (size_t)this->getSize();
      ::memcpy (this->elements, newTable.elements, len);
   }
}

template<class K, class H>
void SxHashTable<K,SxNull,H>::print () const
{
   ssize_t i;
   ssize_t n = this->getSize ();

   std::cout << "table[" << tableSize << "/" << n << "]: ";

   for (i = 0; i < n; i++)  {
      if (this->elements[i] == NULL)  {
         std::cout << i << ":NULL ";
      }  else if (this->elements[i] == &deleted)  {
         std::cout << i << ":DEL ";
      }  else  {
         std::cout << i << ":" << (void*)this->elements[i] << " ";
      }
   }
   std::cout << std::endl;
}
