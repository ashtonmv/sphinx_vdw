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

#ifndef _SX_LOG_BUF_THREAD_H_
#define _SX_LOG_BUF_THREAD_H_

#include <SxIPC.h>
#include <SxLogBuf.h>
#include <SxThreadCond.h>
#include <SxSystemThread.h>

/** \brief Thread safe log bufer

    \b SxLogBufThread = S/PHI/nX Log Buffer Thread

    \author Sixten Boeck, boeck@gemmantics.com
    \author Vaclav Bubnik, bubnik@gemmantics.com */
class SX_EXPORT_IPC SxLogBufThread : public std::streambuf, 
                                     public SxSystemThread
{
   public:
      SxLogBufThread ();
      virtual ~SxLogBufThread ();

      void setFile (const SxString &, int mode = 0600); // #exception

      void stop ();
      virtual void main ();

   private:
      int overflow (int c);
        
   protected:
      FILE *fp;
      int maxWidth;

      SxArray<char> line;
      ssize_t lineSize;

      SxList<SxString> messages;
      SxString buffer;
      bool stopCmd;
      SxMutex mutex;
      SxThreadCond condition;

      SxString write (const SxString &, bool flush);
};

#endif /* _SX_LOG_BUF_THREAD_H_ */
