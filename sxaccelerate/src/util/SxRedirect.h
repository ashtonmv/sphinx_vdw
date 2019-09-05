
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

#ifndef _SX_REDIRECT_H_
#define _SX_REDIRECT_H_

#include <SxUtil.h>
#include <SxString.h>
#include <ostream>
#include <streambuf>
#include <fstream>

/** \brief Redirect stream to a file

    \b SxClass = S/PHI/nX Tee Streamer

\code
   #include <fstream>
   ...
   std::ofstream log ("file.log");
   SxRedirect tee (std::cout, log, true);
   cout << "this goes to stdout and to file.log\n"
\endcode

or simpler

\code
   SxRedirect tee (std::cout, "file.log", true);
   cout << "this goes to stdout and to file.log\n"
\endcode

\code
   SxRedirect tee (std::cout, SxRedirect::DevNull");
   cout << "this goes nowhere. 'Identical to myProg > /dev/null'\n"
\endcode

    \author Sixten Boeck */
class SX_EXPORT_UTIL SxRedirect : public std::streambuf
{
   public:

      enum ZeroDevice { DevNull };  

      SxRedirect (const std::ostream &src, const char *filename, 
                  bool tee=false);
      SxRedirect (const std::ostream &src, std::ostream &dest, bool tee=false);
      SxRedirect (const std::ostream &src, std::ostream *dest, bool tee=false);
      SxRedirect (const std::ostream &src, std::streambuf *dest, bool tee=false);
      SxRedirect (const std::ostream &src, ZeroDevice);
      virtual ~SxRedirect ();

      void enable ();
      void disable ();

   protected:

      enum State { Disabled, Enabled };
      enum Mode { Quiet, Redirect, Tee };
      typedef std::char_traits<char> TCharTraits;

      std::ofstream   fStream;
      std::ostream   *srcStream;
      std::ostream   *destStream;
      std::streambuf *srcBuf;
      std::streambuf *destBuf;
      std::streambuf *origBuf;
      State state;
      Mode  mode;

      virtual int_type overflow (int_type c);
      virtual int sync ();
};

#endif /* _SX_REDIRECT_H_ */
