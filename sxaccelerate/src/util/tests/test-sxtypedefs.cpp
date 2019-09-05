
//  SxTypeDefs.h test

#include <SxUtil.h>
#include <SxString.h>
#include <SxTypeDefs.h>


#define SX_TEST(p)        \
           SX_CHECK(p);   \
           if (!(p))  {   \
              return -1;  \
           }


int main (int, char **)
{
#ifdef MSVC
   printf ("VisualC++\n");
#else
   printf ("stdint.h\n");
#endif
   
   SX_TEST (sizeof(int8_t)   == 1);
   SX_TEST (sizeof(uint8_t)  == 1);
   SX_TEST (sizeof(int16_t)  == 2);
   SX_TEST (sizeof(uint16_t) == 2);
   SX_TEST (sizeof(int32_t)  == 4);
   SX_TEST (sizeof(uint32_t) == 4);
   SX_TEST (sizeof(int64_t)  == 8);
   SX_TEST (sizeof(uint64_t) == 8);
   
   uint64_t u64 = (uint64_t)1 << 32;
   SX_TEST ((int)(u64 >> 32) == 1);
   
   return 0;
}
