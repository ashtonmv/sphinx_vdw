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

#include <SxLogBuf.h>
#include <SxString.h>
#include <SxTime.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <SxConfig.h>
#include <unistd.h>

SxLogBuf::SxLogBuf ()
   : std::streambuf (),
     maxWidth (78),
     fp (stdout)
//     inSxError(false)
{
   // empty
}

SxLogBuf::SxLogBuf (int bufSize)
   : std::streambuf (),
     maxWidth (78),
     fp (stdout)
//     inSxError(false)
{
   if (bufSize)  {
      char *ptr = new char [static_cast<size_t>(bufSize)];
      setp (ptr, ptr+bufSize);
   }  else  {
     setp (0,0);
   } 
}

SxLogBuf::~SxLogBuf ()
{
   sync ();
   delete [] pbase();

   if (fp && fp != stdout) fclose (fp);
}

void SxLogBuf::setWidth (int w)
{
   maxWidth = w;
}

void SxLogBuf::setPrefix (const SxString &prefix_)
{
   prefix = prefix_;
}

void SxLogBuf::setFile (const SxString &name, int mode)
{
   if (fp && fp != stdout) fclose (fp);
   // --- set the access permissions (new or change the current)
   int fd = ::open (name.ascii(), O_WRONLY | O_CREAT, mode);
   if (fd < 0)  {
      SX_THROW ("Can't open log file '"+name+"': " + sxstrerror ());
   }
   ::close (fd);
   sxchmod (name.ascii(), mode);
   // --- open file
   if ( !(fp = fopen (name.ascii(), "a")) )  {
      SX_THROW ("Can't append to log file '"+name+"': " + sxstrerror ());
   }
   
}

SxString SxLogBuf::getPrefix ()
{
   if (prefix != "")  return prefix;

   return SxTime::strftime ("%m/%d/%y %H:%M:%S: ");
}

SxString SxLogBuf::write (const SxString &line_, bool flush)
{
   SxString linePrefix = getPrefix ();
   SxList<SxString> tokens = line_.tokenize ('\n');

   if (line_.tail(1) == "\n") flush = true;

   ssize_t i = 0;
   ssize_t n = tokens.getSize ();
   ssize_t max = flush ? n : n-1;
   SxList<SxString>::ConstIterator it;
   for (i=0, it = tokens.begin(); i < max; ++it, ++i)  {
      if (fp) fprintf (fp, "%s%s\n", linePrefix.ascii (), it->ascii ());
   }
   if (fp) fflush (fp);

   // add last line to buffer
   return (flush) ? "" : (*it);
}

int SxLogBuf::overflow (int c)
{
   sync ();
   if (c != EOF)  {
      if (pbase() == epptr())  {
         if (c == '\n')  {

            line = write (line, true);

         }  else  {
            line += char(c);

            if (line.getSize() >= maxWidth)  {
               line = write (line);
            }
         }
      }  else  {
         sputc (char(c));
      }
   }
   return 0;
}

int SxLogBuf::sync ()
{
   if (pbase() != pptr())  {
      int len = int(pptr() - pbase());
      line += SxString ( pbase(), len );
      line  = write (line);
      setp (pbase(), epptr());
   }
   return 0;
}

