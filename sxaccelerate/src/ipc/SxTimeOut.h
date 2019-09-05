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

#ifndef _SX_TIMEOUT_H_
#define _SX_TIMEOUT_H_

#include <SxTimer.h>
#include <SxThread.h>
#include <SxPtr.h>

class SxTimeOutException
{
   public:
      
};

/** \brief ...

    \b SxClass = S/PHI/nX ...

    try {
      SX_TIMEOUT (5, 
         {
            for (int i=0; i < 2; ++i)  {
               printf ("hi\n"); fflush (stdout); 
               SxThread::sleep (1);
            }
            printf ("END\n"); fflush (stdout);
         });
   } catch (SxTimeOutException e)  {
      printf ("exceeded!\n");
   }


    \author Sixten Boeck, boeck@sfhingx.de */
class SxTimeOut : public SxThread
{
   public:
      SxTimeOut (SxThread *, double sec);
      virtual ~SxTimeOut ();
      virtual void main ();

      bool hasExceeded () const;

   protected:
      SxThread *task;
      SxTimer  timer;
      double timeout;  // msec
      bool exceeded;
};

#define SX_UNIQUE_TIMEOUT_VAR(x)   sxtimeouthandle ## x
#define SX_UNIQUE_TIMEOUT_NAME(x)  sxtimeout ## x
#define SX_UNIQUE_TIMEOUT_CLASS(x) SxTimeOut ## x
#define SX_TIMEOUT(sec, expr)  \
   class SX_UNIQUE_TIMEOUT_CLASS(__LINE__) : public SxThread   \
   {\
      public: \
         virtual void main () \
         {\
            expr;  \
         }\
   } SX_UNIQUE_TIMEOUT_NAME(__LINE__);  \
   SxTimeOut SX_UNIQUE_TIMEOUT_VAR(__LINE__) (&SX_UNIQUE_TIMEOUT_NAME(__LINE__), sec); \
   SX_UNIQUE_TIMEOUT_VAR(__LINE__).start (); \
   SX_UNIQUE_TIMEOUT_VAR(__LINE__).wait ();  \
   if (SX_UNIQUE_TIMEOUT_VAR(__LINE__).hasExceeded ())  throw SxTimeOutException ();


#endif /* _SX_TIMEOUT_H_ */
