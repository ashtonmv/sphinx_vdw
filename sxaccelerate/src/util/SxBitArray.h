// Headerguard
#ifndef _SX_BITS_H_
#define _SX_BITS_H_

// Including the header-files
// -- Qt-Headers
// -- SFHIngX- respectively PHInaX-Headers
#include <SxArray.h>
#include <SxUtil.h>

// -- C++-Standardheaders
// -- C-Standardheaders
#include <iostream>
using namespace std;


// SxBitArray-class
/**
  \author Thomas Uchdorf, t.uchdorf@mpie.de
  \brief
 
  \example
  \code
   SxBitArray bitArr1 (25);
   SxBitArray bitArr2 (12);
   SxBitArray bitArr3 (5);
   cout << "bitArr1:                   " << bitArr1 << endl;
   cout << "bitArr1.getSize ():        " << bitArr1.getSize () << endl;
   cout << "bitArr1.setBit (3):        " << bitArr1.setBit (3) << endl;
   cout << "bitArr1.setBit (3):        " << bitArr1.setBit (3) << endl;
   cout << "bitArr1:                   " << bitArr1 << endl;
   cout << "~bitArr1:                  " << ~bitArr1 << endl;
   cout << "(~bitArr1).unsetBit (5):   " <<  (~bitArr1).unsetBit (5) << endl;
   cout << "(~bitArr1).unsetBit (3):   " <<  (~bitArr1).unsetBit (3) << endl;
   cout << "bitArr1 = 0:               " <<  (bitArr1 = 0) << endl;
   cout << "bitArr2:                   " << bitArr2 << endl;
   cout << "bitArr2.getSize ():        " << bitArr2.getSize () << endl;
   cout << "bitArr2 = 1:               " << (bitArr2 = 1) << endl;
   cout << "bitArr3:                   " << bitArr3 << endl;
   cout << "bitArr3.setBit (4):        " << bitArr3.setBit (4) << endl;
   cout << "bitArr3.setBit (0):        " << bitArr3.setBit (0) << endl;
   cout << "bitArr3.setBit (1):        " << bitArr3.setBit (1) << endl;
   SxBitArray bitArr4 (bitArr3);
   cout << "bitArr4:                   " << bitArr4 << endl;
   bitArr3.resize (12, true); 
   cout << "bitArr3.resize (12, true): " << bitArr3 << endl; 
   bitArr4.resize (12);
   cout << "bitArr4.resize (12):       " << bitArr4 << endl;
   cout << "bitArr4.setBit (4):        " << bitArr4.setBit (4) << endl;
   cout << "bitArr4.setBit (0):        " << bitArr4.setBit (0) << endl;
   cout << "bitArr4.setBit (10):       " << bitArr4.setBit (10) << endl;
   cout << "bitArr4.setBit (7):        " << bitArr4.setBit (7) << endl;
   cout << "bitArr3 & bitArr4:         " << (bitArr3 & bitArr4) << endl;
   cout << "bitArr3 | bitArr4:         " << (bitArr3 | bitArr4) << endl;
   cout << "bitArr3 ^ bitArr4:         " << (bitArr3 ^ bitArr4) << endl;
   cout << "bitArr3:                   " << bitArr3 << endl;
   cout << "bitArr4 = (bitArr3 << 5):  " << (bitArr4 = (bitArr3 << 5)) << endl;
   cout << "(bitArr3 >> 2):            "<< (bitArr3 >> 2) << endl; 
   SxBitArray bitArr5 (12), bitArr6 (12);
   cout << "bitArr5:                   " << bitArr5 << endl;
   cout << "bitArr6:                   " << bitArr6 << endl;
   cout << "bitArr5 == bitArr6         " << (bitArr5 == bitArr6) << endl;
   cout << "bitArr5.setBit (9):        " << bitArr5.setBit (9) << endl;
   cout << "bitArr5 == bitArr6         " << (bitArr5 == bitArr6) << endl;
   cout << "bitArr6 = 1:               " << (bitArr6 = 1) << endl;
   for (int i = 0; i < 9; ++i)  {
      cout << "bitArr5.setBit (" << i << "):        " << bitArr5.setBit (i);
      cout << endl;
   }
   for (int i = 10; i < 12; ++i)  {
      cout << "bitArr5.setBit (" << i << "):       " << bitArr5.setBit (i); 
      cout << endl;
   } 
   cout << "bitArr5 == bitArr6         " << (bitArr5 == bitArr6) << endl;
   //output
   25
   12
   5
   bitArr1:                   0000000000000000000000000
   bitArr1.getSize ():        25
   bitArr1.setBit (3):        0001000000000000000000000
   bitArr1.setBit (3):        0001000000000000000000000
   bitArr1:                   0001000000000000000000000
   ~bitArr1:                  1110111111111111111111111
   (~bitArr1).unsetBit (5):   1110101111111111111111111
   (~bitArr1).unsetBit (3):   1110111111111111111111111
   bitArr1 = 0:               0000000000000000000000000
   bitArr2:                   000000000000
   bitArr2.getSize ():        12
   bitArr2 = 1:               111111111111
   bitArr3:                   00000
   bitArr3.setBit (4):        00001
   bitArr3.setBit (0):        10001
   bitArr3.setBit (1):        11001
   bitArr4:                   11001
   bitArr3.resize (12, true): 110010001101
   bitArr4.resize (12):       000000000000
   bitArr4.setBit (4):        000010000000
   bitArr4.setBit (0):        100010000000
   bitArr4.setBit (10):       100010000010
   bitArr4.setBit (7):        100010010010
   bitArr3 & bitArr4:         100010000000
   bitArr3 | bitArr4:         110010011111
   bitArr3 ^ bitArr4:         010000011111
   bitArr3:                   110010001101
   bitArr4 = (bitArr3 << 5):  000001100000
   (bitArr3 >> 2):            001000000110
   12
   12
   bitArr5:                   000000000000
   bitArr6:                   000000000000
   bitArr5 == bitArr6         1
   bitArr5.setBit (9):        000000000100
   bitArr5 == bitArr6         0
   bitArr6 = 1:               111111111111
   bitArr5.setBit (0):        100000000100
   bitArr5.setBit (1):        110000000100
   bitArr5.setBit (2):        111000000100
   bitArr5.setBit (3):        111100000100
   bitArr5.setBit (4):        111110000100
   bitArr5.setBit (5):        111111000100
   bitArr5.setBit (6):        111111100100
   bitArr5.setBit (7):        111111110100
   bitArr5.setBit (8):        111111111100
   bitArr5.setBit (10):       111111111110
   bitArr5.setBit (11):       111111111111
   bitArr5 == bitArr6         1

 */
class SX_EXPORT_UTIL SxBitArray 
{

   public:
      SxArray<unsigned char> bytes;

      // Constructors
      // -- Constructor
      SxBitArray ();

      SxBitArray (size_t nBits_);
      SxBitArray (const unsigned char *, size_t nBits_);

      // -- Copy-Constructor
      SxBitArray (const SxBitArray & rhs);
      SxBitArray (SxBitArray &&);

      // -- Creates a bit array from int array
      SxBitArray (const SxArray<int> &);

      // -- Destructor
      ~SxBitArray ();

      // Methods
      void resize (size_t nBits_, bool keep = false);

      SxBitArray & set (bool);

      SxBitArray & setBit (size_t iBit);

      SxBitArray & unsetBit (size_t iBit);

      SxBitArray & set (size_t iBit, bool b);

      const size_t & getSize () const;

      // Operators
      SxBitArray & operator= (const SxBitArray &rhs);
      SxBitArray & operator= (SxBitArray &&rhs);

      // bitArray = 1;  bitArray = 0
      SxBitArray & operator= (int rhs);

      bool operator() (size_t iBit) const;

        // binary and
      SxBitArray operator& (const SxBitArray &rhs) const;

      // binary or
      SxBitArray operator| (const SxBitArray &rhs) const;

      // binary exclusive or
      SxBitArray operator^ (const SxBitArray &rhs) const;

      // binary not
      SxBitArray operator~ () const;

      // comparison
      bool operator== (const SxBitArray &rhs) const;

      SxBitArray operator<< (int n) const;

      SxBitArray operator>> (int n) const;

      SxBitArray & operator<<= (int n);

      SxBitArray & operator>>= (int n);

   protected:

      // Methods
      void setAllBits ();

      void unsetAllBits ();

      // Members
      size_t nBits;
};

SX_EXPORT_UTIL ostream &operator<< (ostream & os, const SxBitArray &rhs);

#endif /* _SX_BITS_H_ */
