#include <SxLog.h>
#include <iostream>

SxLog::SxLog ()
{
   // empty
}

SxLog::~SxLog ()
{
   // empty
}


SxLog &SxLog::getGlobalObj ()
{
   static SxLog sxLog;
   return sxLog;
}



// very simple but fast hash: http://www.cse.yorku.ca/~oz/hash.html
uint32_t SxLog::djb2 (const char *str)
{
   uint32_t hash = 5381, c = 0;
   while ( (c = static_cast<uint32_t>(*str++)) != 0 )  hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
   return hash;
}


bool SxLog::isEnabled (const uint32_t &compHash)
{
   if (compHash == 0)  return true;  // top-level apps
   bool res = false;
   SxLog::lock ();
   res = SxLog::getGlobalObj().compIDs.contains (compHash);
   std::cout << "isEnabled " << compHash << " = " << res << std::endl;
   SxLog::unlock ();
   return res;
}


void SxLog::enable (const char *soName)
{
   uint32_t hash = djb2 (soName);

   SxLog::getGlobalObj().compIDs << hash;
}


void SxLog::disable (const char *soName)
{
   uint32_t hash = djb2 (soName);

   SxLog::getGlobalObj().compIDs.removeElement (hash);
}


void SxLog::lock ()
{
#  ifdef SX_LOG
      SxUtil::getGlobalObj().lockLog ();
#  endif
}

void SxLog::unlock ()
{
#  ifdef SX_LOG
      SxUtil::getGlobalObj().unlockLog ();
#  endif
}
