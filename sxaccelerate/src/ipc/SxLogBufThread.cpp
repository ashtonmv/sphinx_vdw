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

#include <SxLogBufThread.h>

#include <SxString.h>
#include <SxTime.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <SxConfig.h>
#include <unistd.h>

SxLogBufThread::SxLogBufThread ()
   : fp(NULL),
     maxWidth(78),
     line(1000),
     lineSize(0),
     stopCmd(false)
{
   // empty
}

SxLogBufThread::~SxLogBufThread ()
{
   stop ();
   wait ();
   if (fp) fclose (fp);
}

void SxLogBufThread::setFile (const SxString &name, int mode)
{
   if (fp) fclose (fp);
   // --- set the access permissions (new or change the current)
   int fd = ::open (name.getElems (), O_WRONLY | O_CREAT, mode);
   if (fd < 0)  {
      SX_THROW ("Can't open log file '"+name+"' for writing: "
                + sxstrerror ());
   }
   ::close (fd);
   sxchmod (name.getElems (), mode);
   // --- open file
   if ( !(fp = fopen (name.getElems (), "a")) )  {
      SX_THROW ("Can't open log file '"+name+"' for writing: "
                + sxstrerror ());
   }
}

int SxLogBufThread::overflow (int c)
{
   if (c >= 0 && c < 256)  {
      SX_MUTEX (mutex)  {
         line(lineSize++) = static_cast<char>(c);
         if (c == 10 || lineSize >= line.getSize ())  {
            // --- write to file on newline or full line buffer
            messages << SxString(line.elements, lineSize);
            lineSize = 0;
            condition.wakeOne ();
         }
      }
   }
   return 0;
}

void SxLogBufThread::stop ()
{
   SX_MUTEX (mutex)  {
      stopCmd = true;
      condition.wakeAll ();
   }
}

void SxLogBufThread::main ()
{
   bool workingLoop = true;
   
   mutex.lock ();
   while (workingLoop)  {
      // --- wait for some messages
      while (!stopCmd && messages.getSize () < 1)  {
         condition.wait (&mutex);
      }
      
      // --- collect all messages
      buffer += SxString::join (messages);
      messages.removeAll ();
      
      // --- print mesages I/O
      if (buffer.contains ("\n"))  {
         mutex.unlock ();
         buffer = write (buffer, false);
         mutex.lock ();
      }
      
      // --- flush and exit thread main
      if (stopCmd)  {
         buffer += SxString(line.elements, lineSize);
         lineSize = 0;
         write (buffer, true);
         workingLoop = false;
      }
   }
   mutex.unlock ();
}

SxString SxLogBufThread::write (const SxString &line_, bool flush)
{
   SxString linePrefix = SxTime::strftime ("%m/%d/%y %H:%M:%S: ");
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

