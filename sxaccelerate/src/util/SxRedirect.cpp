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

#include <SxRedirect.h>
#include <SxError.h>

SxRedirect::SxRedirect (const std::ostream &src, const char *filename, bool tee)
   : fStream (filename),
     srcStream (const_cast<std::ostream *>(&src)),
     destStream (&fStream),
     srcBuf(NULL),
     destBuf(NULL),
     origBuf(srcStream->rdbuf()),
     state(SxRedirect::Disabled),
     mode(SxRedirect::Redirect)
{
   if (tee)  mode = Tee;
   else      srcStream->rdbuf(0);
   enable ();
}

SxRedirect::SxRedirect (const std::ostream &src, std::ostream &dest, bool tee)
   : srcStream (const_cast<std::ostream *>(&src)),
     destStream (const_cast<std::ostream *>(&dest)),
     srcBuf(NULL),
     destBuf(NULL),
     origBuf(srcStream->rdbuf()),
     state(SxRedirect::Disabled),
     mode(SxRedirect::Redirect)
{
   if (tee)  mode = Tee;
   else      srcStream->rdbuf(0);
   enable ();
}

SxRedirect::SxRedirect (const std::ostream &src, std::ostream *dest, bool tee)
   : srcStream (const_cast<std::ostream *>(&src)),
     destStream (dest),
     srcBuf(NULL),
     destBuf(NULL),
     origBuf(srcStream->rdbuf()),
     state(SxRedirect::Disabled),
     mode(SxRedirect::Redirect)
{
   if (tee)  mode = Tee;
   else      srcStream->rdbuf(0);
   enable ();
}

SxRedirect::SxRedirect (const std::ostream &src, std::streambuf *dest, bool tee)
   : srcStream (const_cast<std::ostream *>(&src)),
     destStream (NULL),
     srcBuf(NULL),
     destBuf(dest),
     origBuf(srcStream->rdbuf()),
     state(SxRedirect::Disabled),
     mode(SxRedirect::Redirect)
{
   if (tee)  mode = Tee;
   else      srcStream->rdbuf(0);
   enable ();
}

SxRedirect::SxRedirect (const std::ostream &src, SxRedirect::ZeroDevice)
   : srcStream (const_cast<std::ostream *>(&src)),
     srcBuf(NULL),
     destBuf(NULL),
     origBuf(srcStream->rdbuf()),
     state(SxRedirect::Disabled),
     mode(SxRedirect::Quiet)
{
   enable ();
}


SxRedirect::~SxRedirect ()
{
   disable ();
}


void SxRedirect::enable ()
{
   if (state == Enabled)  return;
   state = Enabled;

   switch (mode)  {
      case Quiet    : SX_CHECK (srcStream);
                      srcStream->rdbuf(0);
                      break;
      case Redirect : SX_CHECK (srcStream);
                      srcBuf = srcStream->rdbuf();
                      if (destStream)  {
                         srcStream->rdbuf (destStream->rdbuf());
                      } else {
                         srcStream->rdbuf (destBuf);
                      }
                      break;
      case Tee      : SX_CHECK (srcStream);
                      SX_CHECK (destStream || destBuf);
                      srcBuf = srcStream->rdbuf();
                      if (destStream)  destBuf = destStream->rdbuf();
                      srcStream->rdbuf (this);
                      break;
   }
}

void SxRedirect::disable ()
{
   if (state == Disabled)  return;
   state = Disabled;

   srcStream->flush ();
   if (mode == Redirect || mode == Tee)  {
      if (destStream)  destStream->flush ();
   }

   srcStream->rdbuf (origBuf);
}


SxRedirect::int_type SxRedirect::overflow (int_type c)
{
   if (!TCharTraits::eq_int_type (c, TCharTraits::eof()))  {
//    c = static_cast<SxOut *>(buffer1)->sputc (c);
      c = srcBuf->sputc (static_cast<char>(c));
      if (!TCharTraits::eq_int_type (c, TCharTraits::eof()))
         c = static_cast<SxRedirect *>(destBuf)->overflow (c);
      return c;
   }  else  {
      return TCharTraits::not_eof (c);
   }
}

int SxRedirect::sync ()
{
// int rc = static_cast<SxOut *>(buffer1)->sync ();
   int rc = -1;
   if (srcBuf)  rc = srcBuf->pubsync ();
// if (rc != -1)  rc = static_cast<SxOut *>(buffer2)->sync ();
   if (rc != -1 && destBuf)  rc = destBuf->pubsync ();
   return rc;
}

