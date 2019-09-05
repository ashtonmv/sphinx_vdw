#include <SxFileIO.h>
#include <SxFSAction.h>

#ifdef WIN32
#  include <windows.h>
#else
#  include <dirent.h>
#endif /* WIN32 */

#include <errno.h>


void SxFileIO::syncDir (const SxString &path_)
{
   SX_TRACE ();
#  ifndef WIN32
      DIR *dir = opendir (path_.getElems ());
      if (!dir)
         SX_THROW (SxFileIO::dirOpenErr (path_));

      // or int fd = open (path_.getElems (), O_RDONLY);
      int fd = dirfd (dir);
      if (fd < 0)
         SX_THROW (SxFileIO::dirFDErr (path_));

      if (fsync (fd) != 0)
         SX_THROW (SxFileIO::dirSyncErr (path_));

      if (closedir (dir) != 0)
         SX_THROW (SxFileIO::dirCloseErr (path_));
#  endif /* WIN32 */
}
// --------------------------------------------------------------------------

SxFileIO::SxFileIO ()
   : fp(NULL),
     fileMode(0),
     buffering(true),
     updateMode(false),
     prevOp (OpType::None),
     fileSize(0),
     fileOffset(0),
     path(""),
     lastBytesRead(0),
     lastBytesWritten(0)
{
   SX_TRACE ();
}

SxFileIO::~SxFileIO ()
{
   SX_TRACE ();
   close ();
}

void SxFileIO::enableMode (unsigned int mode_)
{
   SX_TRACE ();
   fileMode |= mode_;
}

void SxFileIO::disableMode (unsigned int mode_)
{
   SX_TRACE ();
   fileMode &= ~(mode_);
}

bool SxFileIO::hasMode (enum SxFileIO::Mode mode_) const
{
   SX_TRACE ();
   return (fileMode & static_cast<unsigned int>(mode_));
}

bool SxFileIO::isBuffering () const
{
   SX_TRACE ();
   return buffering;
}

void SxFileIO::setBuffering (bool isBuffering_)
{
   SX_TRACE ();
   buffering = isBuffering_;
}

void SxFileIO::setModeFlags ()
{
   SX_TRACE ();
   SxString m;
   if (hasMode (SxFileIO::ModeAppend)) {
      m += "a";
      if (hasMode (SxFileIO::ModeRead))
         m += "+";

      // Force binary mode
      //if (hasMode (SxFileIO::ModeBinary))
      m += "b";

      modeStr = m;
      return;

   } else if (   hasMode (SxFileIO::ModeRead)
             && !hasMode (SxFileIO::ModeCreate))
   {
      m += "r";
      if (hasMode (SxFileIO::ModeWrite))
         m += "+";

      // Force binary mode
      //if (hasMode (SxFileIO::ModeBinary))
      m += "b";

      modeStr = m;
      return;

   } else if (hasMode (SxFileIO::ModeWrite))
   {
      m += "w";
      if (hasMode (SxFileIO::ModeRead))
         m += "+";

      // Force binary mode
      //if (hasMode (SxFileIO::ModeBinary))
      m += "b";

      modeStr = m;
      return;
   }
}

void SxFileIO::setSize ()
{
   SX_TRACE ();
   int64_t lastOffset = static_cast<int64_t>(fileOffset);
   fSeek (0, SEEK_END);
   fileSize = fTell ();
   fSeek (lastOffset, SEEK_SET);
   fileOffset = static_cast<uint64_t>(lastOffset);
}

void SxFileIO::open (const SxString &path_,
                     const SxString &mode_,
                     int            permission_)
{
   SX_TRACE ();
   SX_DBG_MSG ("path='" << path_ <<
               "', mode='" << mode_ <<
               "', permission=" << SxString::sprintf("%04o", permission_));

   close ();
   path = path_;

   // --- get the filemode
   fileMode = 0;
   ssize_t i = 0;
   ssize_t n = mode_.getSize ();
   while (i < n)  {
      char c = mode_(i);
      if (c == 'r') {
         if (i + 1 < n && mode_(i + 1) == '+') {
            enableMode ((SxFileIO::ModeRead | SxFileIO::ModeWrite));
            updateMode = true;
            i++;
         } else {
            enableMode (SxFileIO::ModeRead);
         }
      } else if (c == 'w')  {
         if (i + 1 < n && mode_(i + 1) == '+') {
            enableMode ((SxFileIO::ModeRead  |
                         SxFileIO::ModeWrite |
                         SxFileIO::ModeCreate));
            updateMode = true;
            i++;
         } else {
            enableMode (SxFileIO::ModeWrite |
                        SxFileIO::ModeCreate);
         }
      } else if (c == 'a') {
         if (i + 1 < n && mode_(i + 1) == '+') {
            enableMode (SxFileIO::ModeRead   |
                        SxFileIO::ModeWrite  |
                        SxFileIO::ModeAppend |
                        SxFileIO::ModeCreate);
            updateMode = true;
            i++;
         } else {
            enableMode (SxFileIO::ModeWrite  |
                        SxFileIO::ModeAppend |
                        SxFileIO::ModeCreate);
         }
      } else if (c == 'b') {
         enableMode (SxFileIO::ModeBinary);
      } else {
         SX_EXIT; // invalid mode
      }
      i++;
   }
   setModeFlags ();

   // --- create file if it does not exists
   SxFileInfo f(path);
   path = f.getAbsPath ();
   if (!f.exists ()) {

      if (!hasMode (SxFileIO::ModeCreate))
         SX_THROW (SxFileIO::noSuchFileErr (path));

      // --- create new file in corresponding mode
#     ifdef WIN32
         if (path.isUnicode ())
            fp = _wfopen ((const wchar_t *)path.utf16 ().elements, L"w+b");
         else fp = fopen (path.getElems (), "w+b");
#     else
         fp = fopen (path.getElems (), "w+b");
#     endif

      if (!fp) {
         SX_THROW (SxFileIO::fileCreateErr (path));
      }
      if (fclose (fp) != 0) {
         fp = NULL;
         SX_THROW (SxFileIO::fileCloseErr (path));
      }
      fp = NULL;
      SxFSAction::chmod (permission_, path);
   }

   // --- open file in the correct mode
#  ifdef WIN32
      if (path.isUnicode ())
         fp = _wfopen ((const wchar_t *)path.utf16 ().elements,
                       (const wchar_t *)modeStr.utf16 ().elements);
      else fp = fopen (path.ascii (), modeStr.ascii ());
#  else
      fp = fopen (path.getElems (), modeStr.getElems ());
#  endif

   if (!fp) {
      SX_THROW (SxFileIO::fileOpenErr (path));
   }

   // set internal read/write caching to zero
   if (!buffering)  setvbuf (fp, NULL, _IONBF, 0);

   setSize ();
   tell (); // fetch/set current file pointer
   return;
}

SxArray<char> SxFileIO::readLines (int64_t nLines)
{
   SX_TRACE ();
   SX_CHECK (isOpen () && !hasMode (SxFileIO::ModeBinary));
   SX_CHECK (hasMode (SxFileIO::ModeRead));

   if (updateMode && (prevOp & OpType::WriteOp)) {
      // call SxFileIO.seek (0, SxFileIO::Seek::CUR) or Flush
      SX_CHECK ((prevOp & (OpType::FlushOp | OpType::SeekOp)) != 0);
   }

   if (nLines <= 0) {
      nLines = static_cast<int64_t>(fileSize); // max lines possible
   }

   if (fileSize <= 0) return SxArray<char> ();

   int64_t count = 0;
   int64_t pos = 0;

   SxArray<char> buffer (static_cast<ssize_t>(fileSize));
   int ch;
   while ((ch = fgetc (fp)) != EOF) {
      if ((char)ch == '\0') // ignore null characters
         continue;
      buffer (pos++) = (char)ch;
      if ((char)ch == '\n')  count++;
      if (count == nLines)  break;

   }
   if (ferror (fp)) {
      SX_THROW ("FileReadError", SxFileIO::fileReadErr (path));
   }

   prevOp = OpType::ReadOp;

   return SxArray<char> (buffer.elements, pos);
}

SxArray<char> SxFileIO::readTail (int64_t nLines)
{
   SX_TRACE ();
   SX_CHECK (isOpen () && !hasMode (SxFileIO::ModeBinary));
   SX_CHECK (hasMode (SxFileIO::ModeRead));

   if (updateMode && (prevOp & OpType::WriteOp)) {
      // call SxFileIO.seek (0, SxFileIO::Seek::CUR) or Flush
      SX_CHECK ((prevOp & (OpType::FlushOp | OpType::SeekOp)) != 0);
   }

   if (nLines <= 0) {
      nLines = static_cast<int64_t>(fileSize); // max lines possible
   }

   if (fileSize <= 0) return SxArray<char> ();

   int64_t count = 0;
   int64_t pos = 0;
   SxArray<char> buffer (static_cast<ssize_t>(fileSize));

   int ch;
   while ((ch = fgetc (fp)) != EOF) {
      if ((char)ch == '\0') // ignore null characters
         continue;
      if ((char)ch == '\n')  count++;
   }

   if (ferror (fp)) {
      SX_THROW ("FileReadError", SxFileIO::fileReadErr (path));
   }

   seek (0, Seek::BEG);
   int64_t tmpLines = 0;
   SxString res;
   ch = 0;
   while ((ch = fgetc (fp)) != EOF) {
      if ((char)ch == '\0')
         continue;
      if ((char)ch == '\n')  tmpLines++;
      if (tmpLines >= (count - nLines) && tmpLines != count ) {
         buffer (pos++) = (char)ch;
      }
   }
   if (ferror (fp)) {
      SX_THROW ("FileReadError", SxFileIO::fileReadErr (path));
   }

   prevOp = OpType::ReadOp;

   return SxArray<char> (buffer.elements, pos);
}

uint64_t SxFileIO::read (SxString *str, int64_t nBytes)
{
   SX_TRACE ();
   SX_CHECK (isOpen ());
   SX_CHECK (fileSize > 0);
   SX_CHECK (fileOffset < fileSize, fileOffset, fileSize);
   SX_DBG_MSG ("size='" << nBytes << "'");

   bool     sizeProvided = false;
   uint64_t nBytesLeft   = fileSize - fileOffset;

   // read rest of the file
   if (nBytes < 0)  nBytes = (int64_t)nBytesLeft;
   else             sizeProvided = true;

   uint64_t nBytesToRead = nBytesLeft < (uint64_t)nBytes ?
                           nBytesLeft : (uint64_t)nBytes;

   if (nBytes == 0 || nBytesLeft == 0) return 0;

   if (sizeProvided && (str->getNBytes () < (int64_t)nBytesToRead))
      SX_THROW ("ReadBufferError",
                SxFileIO::fileReadBufferErr (str->getNBytes ()));

   SxArray<char> buf(static_cast<ssize_t>(nBytesToRead));

   if (hasMode (SxFileIO::ModeRead)) {

      if (updateMode && (prevOp & OpType::WriteOp)) {
         // call SxFileIO.seek (0, SxFileIO::Seek::CUR) or Flush
         SX_CHECK ((prevOp & (OpType::FlushOp | OpType::SeekOp)) != 0);
      }

      size_t nBytesRead = fread (buf.elements, 1, nBytesToRead, fp);
      lastBytesRead = (nBytesRead > 0)? static_cast<uint64_t>(nBytesRead) : 0;
      if (nBytesRead != nBytesToRead)
         SX_THROW (SxFileIO::fileReadErr (path));
      (*str) = SxString(buf.elements, static_cast<ssize_t>(nBytesRead));
      fileOffset += lastBytesRead;

      prevOp = OpType::ReadOp;

      return (uint64_t)nBytesRead;

   } else { SX_EXIT; } // invalid mode
}

uint64_t SxFileIO::readBuffer (void *outPtr, uint64_t nBytes)
{
   SX_TRACE ();
   SX_CHECK (isOpen ());
   SX_CHECK (fileSize > 0); 
   SX_CHECK (fileOffset < fileSize, fileOffset, fileSize);
   SX_DBG_MSG ("size='" << nBytes << "'");

   uint64_t nBytesLeft = fileSize - fileOffset; 
   uint64_t nBytesToRead = nBytesLeft < nBytes ? nBytesLeft : nBytes;

   if (nBytes == 0 || nBytesLeft == 0) return 0;

   if (hasMode (SxFileIO::ModeRead)) {

      if (updateMode && (prevOp & OpType::WriteOp)) {
         // call SxFileIO.seek (0, SxFileIO::Seek::CUR) or Flush
         SX_CHECK ((prevOp & (OpType::FlushOp | OpType::SeekOp)) != 0);
      }

      size_t nBytesRead = fread (outPtr, 1, nBytesToRead, fp);
      lastBytesRead = (nBytesRead > 0)? static_cast<uint64_t>(nBytesRead) : 0;
      if (nBytesRead != nBytesToRead)
         SX_THROW (SxFileIO::fileReadErr (path));
      fileOffset += lastBytesRead;

      prevOp = OpType::ReadOp;

      return (uint64_t)nBytesRead;

   } else { SX_EXIT; } // invalid mode
}

uint64_t SxFileIO::write (const SxString &str, int64_t nBytes)
{
   SX_TRACE ();
   SX_DBG_MSG ("size='" << str.getSize () << "'");
   SX_CHECK (isOpen ());

   uint64_t nBytesToWrite;

   // nBytes == -1, write the whole str
   if (nBytes < 0)
      nBytesToWrite = static_cast<uint64_t>(str.getNBytes ());
   else
      nBytesToWrite = static_cast<uint64_t>(nBytes);

   if (str.getNBytes () < (int64_t)nBytesToWrite)
      SX_THROW ("WriteBufferError",
                SxFileIO::fileWriteBufferErr (str.getNBytes ()));
   if (hasMode(SxFileIO::ModeWrite) || hasMode(SxFileIO::ModeAppend)) {

      if (updateMode && (prevOp & OpType::ReadOp) ) {
         // call SxFileIO.seek (0, SxFileIO::Seek::CUR)
         SX_CHECK ((prevOp & OpType::SeekOp) != 0);
      }

      size_t nBytesWritten = fwrite (str.getElems (), 1, nBytesToWrite, fp);
      lastBytesWritten = (nBytesWritten > 0)? static_cast<uint64_t>(nBytesWritten) : 0;
      if (nBytesWritten != nBytesToWrite)
         SX_THROW ("FileWriteError", SxFileIO::fileWriteErr (path));

      if (hasMode(SxFileIO::ModeAppend)) {
         fileSize += lastBytesWritten;
      } else {
         fileOffset += lastBytesWritten;
         if (fileOffset > fileSize) {
            fileSize = fileOffset;
         }
      }

      prevOp = OpType::WriteOp;

      return (uint64_t)nBytesWritten;

   } else { SX_EXIT; } // invalid mode
}

uint64_t SxFileIO::writeBuffer (const void *inPtr, uint64_t nBytes)
{
   SX_TRACE ();
   SX_DBG_MSG ("size='" << nBytes << "'");
   SX_CHECK (isOpen ());

   if (hasMode(SxFileIO::ModeWrite) || hasMode(SxFileIO::ModeAppend)) {

      if (updateMode && (prevOp & OpType::ReadOp) ) {
         // call SxFileIO.seek (0, SxFileIO::Seek::CUR)
         SX_CHECK ((prevOp & OpType::SeekOp) != 0);
      }

      size_t nBytesWritten = fwrite (inPtr, 1, nBytes, fp);
      lastBytesWritten = (nBytesWritten > 0)? static_cast<uint64_t>(nBytesWritten) : 0;
      if (nBytesWritten != nBytes)
         SX_THROW ("FileWriteError", SxFileIO::fileWriteErr (path));

      if (hasMode(SxFileIO::ModeAppend)) {
         fileSize += lastBytesWritten;
      } else {
         fileOffset += lastBytesWritten;
         if (fileOffset > fileSize) {
            fileSize = fileOffset;
         }
      }

      prevOp = OpType::WriteOp;

      return (uint64_t)nBytesWritten;

   } else { SX_EXIT; } // invalid mode
}

void SxFileIO::cleanup ()
{
   SX_TRACE ();
   SX_CHECK (fp != NULL);
   fp = NULL;
   path = "";
   modeStr = "";
   fileMode = 0;
   fileSize = 0;
   fileOffset = 0;
   lastBytesRead = 0;
   lastBytesWritten = 0;
   buffering   = true;
   updateMode  = false;
   prevOp = OpType::None;

}

void SxFileIO::close ()
{
   SX_TRACE ();
   if (!isOpen()) return;
   SX_DBG_MSG ("path='" << path << "'");
   if (fclose (fp) != 0) {
      cleanup ();
      SX_THROW (SxFileIO::fileCloseErr (path));
   }
   cleanup ();
   //path = ""; keep for error messages
}

bool SxFileIO::isOpen () const
{
   SX_TRACE ();
   return (bool)(fp);
}

bool SxFileIO::isOpenForWriting () const
{
   SX_TRACE ();
   return (isOpen () && hasMode (SxFileIO::ModeWrite));
}

uint64_t SxFileIO::getSize () const
{
   SX_TRACE ();
   return fileSize;
}

uint64_t SxFileIO::getOffset () const
{
   SX_TRACE ();
   return fileOffset;
}

uint64_t SxFileIO::getLastRead () const
{
   SX_TRACE ();
   return lastBytesRead;
}

SxString SxFileIO::getFilePath () const
{
   SX_TRACE ();
   return path;
}

uint64_t SxFileIO::getLastWritten () const
{
   SX_TRACE ();
   return lastBytesWritten;
}

void SxFileIO::fFlush ()
{
   SX_TRACE ();
   if (fp) {
      if (fflush (fp) != 0)
         SX_THROW (SxFileIO::fileFlushErr (path));
   }
}

void SxFileIO::flush ()
{
   SX_TRACE ();
   SX_CHECK (isOpen ());

   fFlush ();
   prevOp |= static_cast<unsigned int>(OpType::FlushOp);
}

uint64_t SxFileIO::fTell ()
{
   SX_TRACE ();
   int64_t offset = 0;
#  if defined(WIN32)
      if ((offset = _ftelli64 (fp)) < 0) {
         SX_THROW (SxFileIO::fileTellErr (path));
      }
#  elif defined(MACOSX)
      if ((offset = ftello (fp)) < 0) {
         SX_THROW (SxFileIO::fileTellErr (path));
      }
#  elif defined(LINUX)
      if ((offset = ftello64 (fp)) < 0) {
         SX_THROW (SxFileIO::fileTellErr (path));
      }
#  else
      if ((offset = ftello (fp)) < 0) {
         SX_THROW (SxFileIO::fileTellErr (path));
      }
#  endif
   return (uint64_t)offset;
}

uint64_t SxFileIO::tell ()
{
   SX_TRACE ();
   SX_CHECK (isOpen ());
   fileOffset = fTell ();
   return fileOffset;
}

void SxFileIO::fSeek (int64_t offset, int origin)
{
   SX_TRACE ();
#  if defined(WIN32)
      if (_fseeki64 (fp, offset, origin) != 0)
         SX_THROW (SxFileIO::fileSeekErr (path));
#  elif defined(LINUX)
      if (fseeko64 (fp, offset, origin) != 0)
         SX_THROW (SxFileIO::fileSeekErr (path));
#  elif defined(MACOSX)
      if ((fileOffset = (uint64_t)fseeko (fp, offset, origin)) != 0)
         SX_THROW (SxFileIO::fileSeekErr (path));
#  else
#     error "fseek for LFS not yet ported."
#  endif
}

void SxFileIO::seek (int64_t offset, SxFileIO::Seek from)
{
   SX_TRACE ();
   SX_CHECK (isOpen ());
   if (offset > (int64_t)fileSize)
      SX_THROW ("Cannot set file '" + path + "' position to given offset");

   fSeek (offset, (int)from);
   fileOffset = fTell ();
   prevOp |= static_cast<unsigned int>(OpType::SeekOp);
}

// --------------------------------------------------------------------------
SxString SxFileIO::noSuchFileErr (const SxString &path)
{
   return "Cannot open file '" + path + "': No such file";
}

SxString SxFileIO::fileModeErr (const SxString &path,
                                const SxString &modes)
{
   return "Cannot open file '" + path + "': File modes '" + modes
          + "' are mutually exclusive.";
}

SxString SxFileIO::fileCreateErr (const SxString &path)
{
   return "Cannot create file '" + path + "': " + sxstrerror ();
}

SxString SxFileIO::fileOpenErr (const SxString &path)
{
   return "Cannot open file '" + path + "': " + sxstrerror ();
}

SxString SxFileIO::fileCloseErr (const SxString &path)
{
   return "Cannot close file '" + path + "': " + sxstrerror ();
}

SxString SxFileIO::fileReadErr (const SxString &path)
{
   return "Read error in file '" + path + "': " + sxstrerror ();
}

SxString SxFileIO::fileReadBufferErr (const SxString &path)
{
   return "Read error in file '" + path +
          "': Buffer size is smaller than requested nElems to read";
}

SxString SxFileIO::fileWriteErr (const SxString &path)
{
   return "Write error in file '" + path + "': " + sxstrerror ();
}

SxString SxFileIO::fileWriteBufferErr (const SxString &path)
{
   return "Write error in file '" + path +
          "': Buffer size is smaller than requested nElems to write";
}

SxString SxFileIO::fileSeekErr (const SxString &path)
{
   return "Cannot set file '" + path + "' position to given offset: "
          + sxstrerror ();
}

SxString SxFileIO::fileTellErr (const SxString &path)
{
   return "Cannot get file '" + path + "' position: " + sxstrerror ();
}

SxString SxFileIO::fileFlushErr (const SxString &path)
{
   return "Cannot flush file '" + path + "': " + sxstrerror ();
}

SxString SxFileIO::fileSyncErr (const SxString &path)
{
   return "Cannot sync file '" + path + "': " + sxstrerror ();
}

SxString SxFileIO::dirOpenErr (const SxString &path)
{
   return "Cannot open directory '" + path + "': " + sxstrerror ();
}

SxString SxFileIO::dirFDErr (const SxString &path)
{
   return "Cannot get descriptor from directory '" + path + "': "
          + sxstrerror ();
}

SxString SxFileIO::dirSyncErr (const SxString &path)
{
   return "Cannot sync directory '" + path + "': " + sxstrerror ();
}

SxString SxFileIO::dirCloseErr (const SxString &path)
{
   return "Cannot close directory '" + path + "': " + sxstrerror ();
}

// --- test functions

void SxFileIO::printModeStr () const
{
   sxprintf ("Mode string: %s\n", modeStr.getElems ());
}

// --- Streambuf functions

SxFileIO::int_type SxFileIO::overflow (int c)
{
   return sputc (static_cast<char>(c));
}

SxFileIO::int_type SxFileIO::sputc (char c)
{
   writeBuffer ((void *)(&c), 1);
   return c;
}

SxFileIO::int_type SxFileIO::uflow ()
{
   return sgetc ();
}

SxFileIO::int_type SxFileIO::underflow ()
{
   int_type c = sgetc ();
   seek (static_cast<int64_t>(--fileOffset), Seek::BEG);
   return c;
}

SxFileIO::int_type SxFileIO::sgetc ()
{
   if (fileOffset >= fileSize)
      return EOF;
   char c;
   read (&c);
   return c;
}
