// Including the header-files
#include <SxBitArray.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
  
// Constructors
// -- Constructor  
SxBitArray::SxBitArray () : nBits (0)
{
   // empty
}

// -- Typecast
SxBitArray::SxBitArray (size_t nBits_) : nBits (nBits_)
{
   SX_CHECK (sizeof(unsigned char) == 1);
   size_t nBytes = (nBits >> 3) + (nBits & 7 ? 1 : 0);
   bytes.resize (static_cast<ssize_t>(nBytes));
   unsetAllBits ();     
}

SxBitArray::SxBitArray (const unsigned char *data_, size_t nBits_)
   : bytes (data_, static_cast<ssize_t>(nBits_)),
     nBits (nBits_)
{
   // empty
}

// -- Copy-Constructor     
SxBitArray::SxBitArray (const SxBitArray & rhs) : 
   bytes (rhs.bytes), nBits (rhs.nBits) 
{
   // empty
}

SxBitArray::SxBitArray (SxBitArray &&rhs) : 
   bytes (std::move(rhs.bytes)), nBits (rhs.nBits) 
{
   // empty
}

// -- Creates a bit array from int array
SxBitArray::SxBitArray (const SxArray<int> &in_)
{
   SX_CHECK (in_.getSize() >= 0, in_.getSize());

   ssize_t i;
   ssize_t n;
   
   n = in_.getSize ();
   resize (static_cast<size_t>(n));
   
   for (i = 0; i < n; ++i)  {
      if (in_(i) != 0)  {
         setBit (static_cast<size_t>(i));
      }
   }   
}

// -- Destructor      
SxBitArray::~SxBitArray ()
{
   // empty
}

void SxBitArray::resize (size_t nBits_, bool keep)
{
   nBits = nBits_;
   size_t nBytes = (nBits >> 3) + (nBits & 7 ? 1 : 0);
   bytes.resize (static_cast<ssize_t>(nBytes), keep);
   if (!keep)  {
      unsetAllBits ();
   }
}

SxBitArray & SxBitArray::set (bool b)
{
   SX_CHECK (nBits > 0);
   if (b)  {
      setAllBits ();
   }  else  {
      unsetAllBits ();
   }
   return *this;
}

SxBitArray & SxBitArray::set (size_t iBit, bool b)
{
   if (b)  {
      return setBit (iBit);
   } else  {
      return unsetBit (iBit);      
   }
}

SxBitArray & SxBitArray::setBit (size_t iBit)
{
   SX_CHECK (iBit < nBits);
   ssize_t iChr = iBit / (sizeof (unsigned char) * 8);
   size_t mod = iBit - static_cast<size_t>(iChr) * (sizeof (unsigned char) * 8);

   // using the or-operation to set the iBits bit 
   // example:
   //   10000110
   // | 00100000
   //   --------
   //   10100110
   bytes(iChr) = static_cast<unsigned char>(bytes(iChr) | (1<<mod));

   return *this;
}

SxBitArray & SxBitArray::unsetBit (size_t iBit)
{
   if (!(*this)(iBit)) return *this;
   
   SX_CHECK (iBit < nBits);
   ssize_t iChr = iBit / (sizeof (unsigned char) * 8);
   size_t mod = iBit - static_cast<size_t>(iChr) * (sizeof (unsigned char) * 8);

   // using the xor-operation to unset the iBits bit 
   // example:
   //   10100110
   // ^ 00100000
   //   --------
   //   10000110
   bytes(iChr) = static_cast<unsigned char>(bytes(iChr) ^ (1<<mod));

   return *this;
}

const size_t &SxBitArray::getSize () const
{
   return nBits;
}

// Operators
SxBitArray & SxBitArray::operator= (int rhs)
{
   // Performing a selftest
   if (rhs == 1)  {
      setAllBits ();    
   } else if (rhs == 0)  {
      unsetAllBits ();
   } else  {
      SX_EXIT;
   }

   return *this;
}

SxBitArray & SxBitArray::operator= (const SxBitArray &rhs)
{
   // Performing a selftest
   if (&rhs == this)  {
      return *this;
   }

   // Assigning members
   bytes = rhs.bytes;
   nBits = rhs.nBits;

   return *this;
}

SxBitArray &SxBitArray::operator= (SxBitArray &&rhs)
{
   // Performing a selftest
   if (&rhs == this)  { 
      return *this;
   }

   // Assigning members
   bytes = std::move (rhs.bytes);
   nBits = rhs.nBits;

   return *this;
}

bool SxBitArray::operator() (size_t iBit) const
{
   SX_CHECK (iBit < nBits, iBit, nBits);
   ssize_t iChr = iBit / (sizeof (unsigned char) * 8);
   size_t mod = iBit - static_cast<size_t>(iChr) * (sizeof (unsigned char) * 8);
   unsigned char testChr = bytes(iChr);

   return (testChr & (1 << mod)) == (1 << mod);
}

SxBitArray SxBitArray::operator& (const SxBitArray &rhs) const
{
   SX_CHECK (nBits == rhs.nBits);   
   SxBitArray ret (*this);
   for (ssize_t iByte = 0; iByte < bytes.getSize (); ++iByte)  {
      ret.bytes(iByte) = bytes(iByte)&(rhs.bytes(iByte));
   }

   return ret;
}

SxBitArray SxBitArray::operator| (const SxBitArray &rhs) const
{

   SX_CHECK (nBits == rhs.nBits);   
   SxBitArray ret (*this);
   for (ssize_t iByte = 0; iByte < bytes.getSize (); ++iByte)  {
      ret.bytes(iByte) = bytes(iByte)|(rhs.bytes(iByte));
   }

   return ret;
}
SxBitArray SxBitArray::operator^ (const SxBitArray &rhs) const
{
   
   SX_CHECK (nBits == rhs.nBits);   
   SxBitArray ret (*this);
   for (ssize_t iByte = 0; iByte < bytes.getSize (); ++iByte)  {
      ret.bytes(iByte) = bytes(iByte)^(rhs.bytes(iByte));
   }

   return ret;
}

SxBitArray SxBitArray::operator~ () const
{
   SxBitArray ret (*this);
   for (ssize_t iByte = 0; iByte <  bytes.getSize (); ++iByte)  {
      ret.bytes(iByte) = static_cast<unsigned char>(~ret.bytes(iByte));
   }

   return ret;
}

bool SxBitArray::operator== (const SxBitArray &rhs) const
{
   SX_CHECK (nBits == rhs.nBits);
   bool ret (true); 
   for (ssize_t iByte = 0;
        ret && iByte < bytes.getSize () - 1;
        ++iByte)  {

      // Comparing the bytes in which all bits are relevant
      ret = (bytes(iByte) == rhs.bytes(iByte));
      
   }

   // Comparing each bit of the last byte that belongs to the SxBitArray and
   // that not just exists to achieve parity
   for (size_t iBit = 0;
         ret &&
         iBit < nBits - (static_cast<size_t>(bytes.getSize ()) - 1) * (sizeof (unsigned char) * 8);
         ++iBit)  {
      ret = ((bytes(bytes.getSize () - 1) & (1 << iBit)) == 
            (rhs.bytes(bytes.getSize () - 1) & (1 << iBit)));

   }

   return ret;
}

SxBitArray SxBitArray::operator<< (int n) const
{
   SxBitArray ret (*this);
   for (ssize_t iByte = 0; iByte < ret.bytes.getSize (); ++iByte)  {
      ret.bytes(iByte) = static_cast<unsigned char>(ret.bytes(iByte) << n);
   }
   return ret;
}

SxBitArray SxBitArray::operator>> (int n) const
{
   SxBitArray ret (*this);
   for (ssize_t iByte = ret.bytes.getSize () - 1; iByte >= 0; --iByte)  {
      ret.bytes(iByte) = static_cast<unsigned char>(ret.bytes(iByte) >> n);
   }
   return ret;
}

SxBitArray & SxBitArray::operator<<= (int n)
{
   for (ssize_t iByte = 0; iByte < bytes.getSize (); ++iByte)  {
      bytes(iByte) = static_cast<unsigned char>(bytes(iByte) << n);
   }
   return *this;
}

SxBitArray & SxBitArray::operator>>= (int n)
{
   for (ssize_t iByte = 0; iByte < bytes.getSize (); ++iByte)  {
      bytes(iByte) = static_cast<unsigned char>(bytes(iByte) >> n);
   }
   return *this;
}

void SxBitArray::setAllBits ()
{
   memset(bytes.elements,'\xFF', static_cast<size_t>(bytes.getSize ()));
}

void SxBitArray::unsetAllBits ()
{
   memset(bytes.elements, '\0', static_cast<size_t>(bytes.getSize ()));
}

ostream & operator<< (ostream & os, const SxBitArray &rhs)
{
   for (size_t iBit = 0; iBit < (size_t) rhs.getSize (); ++iBit)  {
      os << rhs(iBit);
   }
   return os;
}
