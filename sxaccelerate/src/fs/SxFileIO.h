#ifndef _SX_FILE_IO_H_
#define _SX_FILE_IO_H_

#include <SxFS.h>
#include <SxString.h>
#include <SxError.h>
#include <SxException.h>
#include <streambuf>

#ifdef WIN32
#   include <windows.h>
#else
#   define _FILE_OFFSET_BITS 64
#   include <unistd.h>
#endif


/** \brief File I/O

   implemented: - mode 'f'
                - sync // only flushes for fileStream mode
                - syncDir // for linux


   \par Overview
   This class hides file implementation differences on Linux/Mac/Windows.
   Synchronous write on Windows requires different file handle and functions
   than on POSIX systems.

   \par Directory synchronization (syncDir)
   Synchronize directory to hard drive after creating or removing a file or
   a subdirectory.

   \par OS
   Different implementation for POSIX and Windows. Small details for MacOS.

   \par Open File Modes
\verbatim
       Open for             !Existing file   Stream position
       -------------------  ---------------  ---------------
   r   reading,                              beginning
   w   writing,             create/truncate, beginning
   a   writing,             created        , end
   r+  reading and writing,                , beginning
   w+  reading and writing, create/truncate, beginning
   a+  reading and writing, created        , end

   b   binary
\endverbatim
 */
class SX_EXPORT_FS SxFileIO : public std::streambuf
{
   public:
      enum Mode {
         ModeUndefined  = 0,
         ModeRead       = (1 << 0), // r, r+, w+, a+
         ModeWrite      = (1 << 1), // r+, w, w+, a, a+
         ModeAppend     = (1 << 2), // a, a+
         ModeCreate     = (1 << 3), // w, w+, a, a+
         ModeBinary     = (1 << 4), // b
      };

      enum class Seek
      {
         BEG = SEEK_SET,
         CUR = SEEK_CUR,
         END = SEEK_END
      };

      enum OpType {
         None = 0x00,
         ReadOp = 0x01,
         WriteOp = 0x02,
         SeekOp  = 0x04,
         FlushOp = 0x08
      };

      SxFileIO ();
     ~SxFileIO ();

      void open (const SxString &path,
                 const SxString &mode,
                 int            permission = 0600);

      void close ();

      bool isOpen () const;

      bool isOpenForWriting () const;

      bool hasMode (enum Mode) const;

      bool isBuffering () const;
      void setBuffering (bool isBuffering_ = true);

      uint64_t getSize () const;
      uint64_t getOffset () const;

      uint64_t getLastRead () const;
      uint64_t getLastWritten () const;

      SxString getFilePath () const;

      template<typename... Args>
      int printf (const char *format, Args... args);
      template<typename... Args>
      int scanf (const char *format, Args... args);

      template<typename T>
      uint64_t read (SxArray<T> *outBufPtr, int64_t nElems = -1);
      uint64_t read (SxString *str, int64_t nBytes = -1);
      uint64_t readBuffer (void *outPtr, uint64_t nBytes);
      template<typename T>
      uint64_t read (T *ptr);

      SxArray<char> readLines (int64_t nLines = 0);
      SxArray<char> readTail (int64_t nLines = 0);

      template<typename T>
      uint64_t write (const SxArray<T> &buffer, int64_t nElems = -1);
      uint64_t write (const SxString &str, int64_t nBytes = -1);
      uint64_t writeBuffer (const void *inPtr, uint64_t nBytes);
      template<typename T>
      uint64_t write (const T &obj);


      uint64_t tell ();
      void seek (int64_t offset, SxFileIO::Seek from = Seek::CUR);

      // --- flushes the application buffers
      void flush ();

      // syncDir in order to update directory after file creation
      static void syncDir (const SxString &path);

      // --- test functions
      void printModeStr () const;

      // static read/write functions
      static inline SxArray<char> readLines (const SxString &filename,
                                             int64_t nLines=0);
      static inline SxArray<char> readBinary (const SxString &filename,
                                              int64_t nBytes);
      static inline SxArray<char> readTail (const SxString &filename,
                                            int64_t nLines = 0);

      static inline void write (const SxString &str, const SxString &filename,
                                int mode);
      static inline void appendToFile (const SxString &str,
                                       const SxString &filename);


      operator istream& () {
         if (istreamPtr.getPtr () == NULL)
            istreamPtr = SxPtr<istream>::create (static_cast<streambuf*>(this));
         return *istreamPtr;
      }

      operator ostream& () {
         if (ostreamPtr.getPtr () == NULL)
            ostreamPtr = SxPtr<ostream>::create (static_cast<streambuf*>(this));
         return *ostreamPtr;
      }

      template<typename T>
      friend SxFileIO &operator<< (SxFileIO &f, const SxArray<T> &buffer) {
         try {
            f.write (buffer, buffer.getSize ());
         } catch (SxException ex) {
            SX_THROW (ex, "FileWriteError",
                      SxFileIO::fileReadErr (f.getFilePath ()));
         }
         return f;
      }

      template<typename T>
      friend SxFileIO &operator>> (SxFileIO &f, SxArray<T> &buffer) {
         try {
            f.read (&buffer, buffer.getSize ());
         } catch (SxException ex) {
            SX_THROW (ex, "FileReadError",
                      SxFileIO::fileWriteErr (f.getFilePath ()));
         }
         return f;
      }

   private:

      SxPtr<ostream> ostreamPtr;
      SxPtr<istream> istreamPtr;

   protected:


      // --- stream mode
      FILE *fp;

      void enableMode (unsigned int);
      void disableMode (unsigned int);
      void setModeFlags ();
      unsigned fileMode;
      SxString modeStr;    // cstdlib fopen ()

      bool buffering;
      bool updateMode;
      unsigned prevOp;


      void setSize ();
      uint64_t fileSize;
      uint64_t fileOffset;

      SxString path;

      uint64_t lastBytesRead;
      uint64_t lastBytesWritten;


      // --- file stream related
      uint64_t fTell ();
      void fSeek (int64_t offset, int origin);
      void fFlush ();

      void cleanup ();

      static SxString noSuchFileErr (const SxString &path);
      static SxString fileModeErr (const SxString &path,
                                   const SxString &modes);
      static SxString fileCreateErr (const SxString &path);
      static SxString fileOpenErr (const SxString &path);
      static SxString fileCloseErr (const SxString &path);
      static SxString fileReadErr (const SxString &path);
      static SxString fileReadBufferErr (const SxString &path);
      static SxString fileWriteErr (const SxString &path);
      static SxString fileWriteBufferErr (const SxString &path);
      static SxString fileSeekErr (const SxString &path);
      static SxString fileTellErr (const SxString &path);

      static SxString fileFlushErr (const SxString &path);
      static SxString fileSyncErr (const SxString &path);
      static SxString dirOpenErr (const SxString &path);
      static SxString dirFDErr (const SxString &path);
      static SxString dirSyncErr (const SxString &path);
      static SxString dirCloseErr (const SxString &path);


      // --- streambuf functions
      SxFileIO::int_type overflow (int_type c = EOF);
      SxFileIO::int_type sputc (char c);
      SxFileIO::int_type uflow ();
      SxFileIO::int_type underflow ();
      SxFileIO::int_type sgetc ();

};

SxArray<char> SxFileIO::readLines (const SxString &filename, int64_t nLines)
{
   //SX_TRACE ();

   SxFileIO f;
   f.open (filename, "r");
   return f.readLines (nLines);
}

SxArray<char> SxFileIO::readBinary (const SxString &filename,
                                    int64_t nBytesToRead)
{
   //SX_TRACE ();

   SxFileIO f;
   f.open (filename, "rb");
   uint64_t fSize = f.getSize ();
   if (nBytesToRead < 0) nBytesToRead = static_cast<int64_t>(fSize);
   if ((uint64_t)nBytesToRead > fSize) {
      f.close ();
      SX_THROW ("Cannot read " + SxString(nBytesToRead) + " bytes from file '"
                 + filename+"': File size is only "+SxString(fSize)+" bytes.");
   }
   SxArray<char> buffer (nBytesToRead);

   if (nBytesToRead == 0)  return SxArray<char>();

   uint64_t nBytesRead = f.read (&buffer, nBytesToRead);
   f.close ();
   return SxArray<char> (buffer.elements, static_cast<ssize_t>(nBytesRead));
}

SxArray<char> SxFileIO::readTail (const SxString &filename, int64_t nLines)
{
   //SX_TRACE ();

   SxFileIO f;
   f.open (filename, "r");
   return f.readTail (nLines);
}

void SxFileIO::write (const SxString &str, const SxString &filename,
                      int permissions)
{
   //SX_TRACE ();

   const int minPerm = 0;
   const int maxPerm = 7777;

   if (permissions < minPerm || permissions > maxPerm)  {
      SX_THROW ("Cannot open file '" + filename + "' for writing "
                "with permission mode 0" + SxString::sprintf("%04o", permissions)
                );
   }

   SxFileIO f;
   f.open (filename, "wb", permissions);

   int64_t nBytesToWrite = str.getNBytes ();

   SxArray<char> buffer (str.getElems (), nBytesToWrite);

   uint64_t nBytesWritten = f.write<char> (buffer, nBytesToWrite);
   if (nBytesWritten < static_cast<uint64_t>(nBytesToWrite)) {
      f.close ();
      SX_THROW ("FileWriteError", "Cannot write" + SxString(nBytesToWrite) +
                "to file: "+filename);
   }
   f.close ();
}

void SxFileIO::appendToFile (const SxString &str, const SxString &filename)
{
   SxFileIO f;
   f.open (filename, "a");

   int64_t nBytesToWrite = str.getNBytes ();

   SxArray<char> buffer (str.getElems (), nBytesToWrite);

   uint64_t nBytesWritten = f.write (buffer, nBytesToWrite);
   if (nBytesWritten < static_cast<uint64_t>(nBytesToWrite)) {
      f.close ();
      SX_THROW ("FileWriteError", "Cannot append" + SxString (nBytesToWrite) +
                "to file: " + filename);
   }
   f.close ();
}

template<typename... Args>
int SxFileIO::printf (const char *format, Args... args)
{
   //SX_TRACE ();
   SX_CHECK (isOpen () && !hasMode (SxFileIO::ModeBinary));

   if (hasMode (SxFileIO::ModeWrite) || hasMode(SxFileIO::ModeAppend)) {

      if (updateMode && (prevOp & OpType::ReadOp) ) {
         // call SxFileIO.seek (0, SxFileIO::Seek::CUR)
         SX_CHECK ((prevOp & OpType::SeekOp) != 0);
      }

      int res = ::fprintf (fp, format, args...);
      if (ferror (fp)) {
         SX_THROW ("FileWriteError", SxFileIO::fileWriteErr (path));
      }

      prevOp = OpType::WriteOp;

      if (hasMode(SxFileIO::ModeAppend)) {
         fileSize += (uint64_t)res;
      } else {
         fileOffset += (uint64_t)res;
         if (fileOffset > fileSize) {
            fileSize = fileOffset;
         }
      }
      return res;
   } else { SX_EXIT; } // invalid mode
   return 0;
}

template<typename... Args>
int SxFileIO::scanf (const char *format, Args... args)
{
   //SX_TRACE ();
   SX_CHECK (isOpen () && !hasMode (SxFileIO::ModeBinary));

   if (hasMode (SxFileIO::ModeRead)) {

      if (updateMode && (prevOp & OpType::WriteOp)) {
         // call SxFileIO.seek (0, SxFileIO::Seek::CUR) or Flush
         SX_CHECK ((prevOp & (OpType::FlushOp | OpType::SeekOp)) != 0);
      }

      int res = ::fscanf (fp, format, args...);
      if (ferror (fp)) {
         SX_THROW ("FileReadError", SxFileIO::fileReadErr (path));
      }
      fileOffset += (uint64_t)res;
      prevOp = OpType::ReadOp;

      return res;
   } else { SX_EXIT; } // invalid mode
}

template<typename T>
uint64_t SxFileIO::read (SxArray<T> *outBufPtr, int64_t nElems)
{
   //SX_TRACE ();
   SX_CHECK (isOpen ());
   SX_CHECK (fileSize > 0);
   SX_CHECK (fileOffset < fileSize, fileOffset, fileSize);
   //SX_DBG_MSG ("size='" << nElems * static_cast<int64_t>(sizeof (T)) << "'");

   constexpr uint64_t typeSize = sizeof (T);
   uint64_t nItemsLeft = (fileSize - fileOffset) / typeSize;

   // --- nElems == -1, buffer has zero size   : read rest of the file
   //     nElems == -1, buffer has finite size : read buffersize of file
   if ((nElems < 0) && (outBufPtr->getSize () == 0)) {
      nElems = (int64_t)nItemsLeft;
      outBufPtr->resize (nElems);
   } else if (nElems < 0)  {
      nElems = outBufPtr->getSize ();
   }

   // do not read more items than available
   uint64_t nItemsToRead = (nItemsLeft < static_cast<uint64_t>(nElems))
                         ? nItemsLeft : static_cast<uint64_t>(nElems);

   // Buffer length will always correspond to number of elements to read
   if ((uint64_t)outBufPtr->getSize () != nItemsToRead) outBufPtr->resize ((ssize_t)nItemsToRead);

   if (hasMode (SxFileIO::ModeRead)) {

      // nothing to read
      if (nItemsToRead == 0)  return 0;

      if (updateMode && (prevOp & OpType::WriteOp)) {
         // call SxFileIO.seek (0, SxFileIO::Seek::CUR) or Flush
         SX_CHECK ((prevOp & (OpType::FlushOp | OpType::SeekOp)) != 0);
      }

      size_t nItemsRead = fread (outBufPtr->elements, typeSize, nItemsToRead,
                                 fp);
      lastBytesRead = (nItemsRead > 0)
                    ? static_cast<uint64_t>(nItemsRead) * typeSize : 0;

      if (nItemsRead != nItemsToRead)
         SX_THROW (SxFileIO::fileReadErr (path));

      fileOffset += lastBytesRead;
      prevOp = OpType::ReadOp;

      return (uint64_t)nItemsRead;

   } else { SX_EXIT; } // invalid mode
}

template<typename T>
uint64_t SxFileIO::read (T *ptr)
{
   // return number of elements read (shall be one)
   return readBuffer (ptr, sizeof(T)) / sizeof(T);
}

template<typename T>
uint64_t SxFileIO::write (const SxArray<T> &buffer, int64_t nElems)
{
   //SX_TRACE ();
   //SX_DBG_MSG ("size='" << nElems * static_cast<int64_t>(sizeof (T)) << "'");
   SX_CHECK (isOpen ());

   // Trying to write empty array
   if (nElems == -1 && buffer.getSize () == 0)
      return 0;

   // nElems == -1, write the whole buffer
   if (nElems < 0)  nElems = buffer.getSize ();

   SX_CHECK (buffer.getSize () > 0);
   SX_CHECK (buffer.getSize () >= nElems, buffer.getSize (), nElems);

   constexpr uint64_t typeSize = sizeof(T);

   if (hasMode(SxFileIO::ModeWrite) || hasMode(SxFileIO::ModeAppend)) {

      if (updateMode && (prevOp & OpType::ReadOp) ) {
         // call SxFileIO.seek (0, SxFileIO::Seek::CUR)
         SX_CHECK ((prevOp & OpType::SeekOp) != 0);
      }

      size_t nItemsWritten = fwrite (buffer.elements, typeSize,
                                     static_cast<size_t>(nElems), fp);
      lastBytesWritten = (nItemsWritten > 0)
                       ? static_cast<uint64_t>(nItemsWritten) * typeSize : 0;

      if (hasMode(SxFileIO::ModeAppend)) {
         fileSize += lastBytesWritten;
      } else {
         fileOffset += lastBytesWritten;
         if (fileOffset > fileSize) {
            fileSize = fileOffset;
         }
      }
      if (static_cast<int64_t>(nItemsWritten) != nElems)
         SX_THROW ("FileWriteError", SxFileIO::fileWriteErr (path));

      prevOp = OpType::WriteOp;

      return (uint64_t)nItemsWritten;

   } else { SX_EXIT; } // invalid mode
}

template<typename T>
uint64_t SxFileIO::write (const T &obj)
{
   //return number of elements written (shall be one)
   return writeBuffer (&obj, sizeof(T)) / sizeof(T);
}

#endif /* _SX_FILE_IO_H_ */
