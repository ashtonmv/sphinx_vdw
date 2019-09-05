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

#ifndef _SX_LOG_BUF_H_
#define _SX_LOG_BUF_H_

#include <SxUtil.h>
#include <SxString.h>
#include <SxException.h>
#include <streambuf>
#include <sstream>
#include <ostream>
#include <string>
#include <stdarg.h>

/** \brief ...

    \b SxClass = S/PHI/nX ...

\example
\code
   SxLogBuf logBuf;
// logBuf.setWidth (60);
   logBuf.setFile ("abc.log");
   SxRedirect redirect (std::cout, &logBuf, true);
   redirect.enable ();
   cout << "abcdef ghijk lmnopqrs tuvwxyz ";
   cout << "ABCDEF GHIJK LMNOPQRS TUVWXYZ ";
   cout << "abcdef ghijk lmnopqrs tuvwxyz\n";
   SxTime::msleep (2500);
   cout << "012345 67890\n";
   cout << "abcdef ghijk lmnopqrs tuvwxyz ";
   cout << "ABCDEF GHIJK LMNOPQRS TUVWXYZ ";
   cout << "012345 67890 ";
   cout << endl;
\endcode

    \author Sixten Boeck */
class SX_EXPORT_UTIL SxLogBuf : public std::streambuf
{
   public:
      SxLogBuf ();
      SxLogBuf (int bufSize);

      virtual ~SxLogBuf ();

      void setWidth (int);
      void setPrefix (const SxString &);
      void setFile (const SxString &, int mode = 0600); // #exception

      virtual SxString write (const SxString &, bool flush=false);

   private:

      int overflow (int c);
      int sync ();

   protected:
      int maxWidth;
      SxString line;
      SxString prefix;
      FILE *fp;
      //bool inSxError;

      virtual SxString getPrefix ();
};

#endif /* _SX_LOG_BUF_H_ */
