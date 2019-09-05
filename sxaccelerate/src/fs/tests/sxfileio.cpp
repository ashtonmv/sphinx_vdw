// --- Including header-files
#include <iostream>
#include <SxCLI.h>
#include <SxString.h>
#include <SxFile.h>
#include <SxFSAction.h>
#include <SxFileIO.h>
#include <SxTime.h>

#ifdef WIN32
#  include <io.h>
#  include <fcntl.h>
#  include <iosfwd>
#  include <string>
#  include <sstream>
#  include <iostream>
#endif

void SxFS_FileIO (const SxString &rootPath, const SxString &str);

int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   cli.finalize ();

   try  {

      // --- Langauges that will be a good test scenario.

      // English  - 1 byte  - [example]     [word]
      // Greek    - 2 bytes - [παράδειγμα]  [λέξη]
      // Russian  - 2 bytes - [пример]      [слово]
      // Chinease - 3 bytes - [例]          [字]
      // Thai     - 3 bytes - [ตัวอย่าง]      [คำ]

      // --------- Basic tests

      if (SxFileInfo::getTmpStr () == "")
         SX_THROW ("[TMPDIR] or [TMP] environment variable is was not found."
                   "Please verify that the environment variable is set.");
      // GET TMP
      SxDir tmpDir = SxDir(SxFileInfo::getTmpStr ());

      if (!tmpDir.exists ())
         SX_THROW ("[TMPDIR] or [TMP] directory not found. Please verify "
                   "that the directory exists");

      // GET HOME
      if (SxFileInfo::getHomeStr () == "")
         SX_THROW ("'HOME' environment variable is was not found."
                   "Please verify that the environment variable is set.");
      SxDir homeDir = SxDir (SxFileInfo::getHomeStr ());
      if (!homeDir.exists ())
         SX_THROW ("'HOME' directory not found. Please verify that the"
                   "directory exists");

      // Unicode folder paths
      SxList<SxString> unicodePaths;
      unicodePaths << "example"
                   << SxString::unicodeFromUtf8 (u8"παράδειγμα")
                   << SxString::unicodeFromUtf8 (u8"пример")
                   << SxString::unicodeFromUtf8 (u8"例")
                   << SxString::unicodeFromUtf8 (u8"ตัวอย่าง");

      if (unicodePaths.getSize () == 0)  {
        SX_THROW ("No language test cases specified.");
     }

     // Folder in which all tests will be performed
     tmpDir = tmpDir.getAbsPath () + SxString("/SxFileIOTest/");

      // Clean old tests if they exist
      if (SxFSAction::test_d (tmpDir.getAbsPath ())) {
         SxFSAction::rm_r (tmpDir);
      }

      // Check combinations of Unicode, with different encoding size
      for (const SxString &path1 : unicodePaths)  {
         for (const SxString &path2 : unicodePaths)  {

            // Skip repeating cases, except full ASCII test
            if (path1 == path2 && path1 != unicodePaths(0)) continue;

            // All the tests go here.
            SxString rootPath = tmpDir.getAbsPath () + "/" + path1 + "/"
                              + path2 + "/";
            // Just in case the rootPath is already unicode
            rootPath = SxString::unicodeFromUtf8 (rootPath.getElems ());

            SxFSAction::mkdir_p (rootPath, 0777);

            // SxFileIO
            SxFS_FileIO (rootPath, path2);

            // --- Remove the root test folder
            SxString cleanPath = tmpDir.getAbsPath () + "/" + path1 + "/";
            SxFSAction::rm_r (cleanPath);

            if (SxFSAction::test_d (rootPath))
               SX_THROW ("Failed to remove test folder '" + rootPath + "'");

            if (SxFSAction::test_d (cleanPath))
               SX_THROW ("Failed to remove root folder '" + cleanPath + "'");
         }
      }

   } catch (SxException e)  {
      cout << e.toString ().wrap ("ERROR: ") << endl;
      cerr << "Exception stack:" << endl;
      cerr << e.toString (SxException::DebugStack) << endl;
      SX_EXIT;
   }

   return 0;
}

class Vector3 {
public:
   Vector3 (int x_ = 0, int y_ = 0, int z_ = 0)
      : x(x_), y(y_), z(z_)
   {
      SX_TRACE ();
   }
  ~Vector3 (){};
   int x;
   int y;
   int z;
   bool operator== (const Vector3 &vec) const {
      return (this->x == vec.x && this->y == vec.y && this->z == vec.z);
   }
   bool operator!= (const Vector3 &vec) const {
      return !(this->x == vec.x && this->y == vec.y && this->z == vec.z);
   }
   SxString toString () const {
      return SxString::sprintf ("[%i,%i,%i]", x, y, z);
      //return "[" + x + ", " + y + ", " + in.z + "]";
   }
   friend std::ostream &operator<< (std::ostream &s, const Vector3 &in) {
      //s << "[" << in.x << ", " << in.y << ", " << in.z << "]";
      s << in.toString ().getElems ();
      return s;
   }
};


void SxFS_FileIO (const SxString &rootPath, const SxString &str)
{
   SX_TRACE ();
   try {

      SxString fileName = str + ".txt";
      SxString filePath = rootPath + fileName;
      SxString fileData = "abcde";
      SxFileIO file;
      uint64_t nBytesRead = 0;
      uint64_t nBytesWritten = 0;

      if (file.isOpen ())
         SX_THROW ("SxFileIO::open unexpected result");

      // File permissions are set only if file doesn't exist already
      file.open (filePath, "w", 0644);

#ifndef WIN32
         if (SxFileInfo(filePath).getMode () != 0644)
            SX_THROW ("SxFileIO::open/close unexpected file permissions");
#endif
      if (!file.isOpen ())
         SX_THROW ("SxFileIO::open failed");

      if (file.getFilePath () != filePath)
         SX_THROW ("SxFileIO::getFilePath unexpected result");

      if (file.getSize () != 0)
         SX_THROW ("SxFileIO::getSize unexpected result");
      file.close ();

      file.open (filePath, "w");
      if (!file.isOpenForWriting ())
         SX_THROW ("SxFileIO::isOpenForWriting unexpected result");

      // Write data to file
      nBytesWritten = file.write (fileData);

      if (nBytesWritten != (uint64_t)fileData.getSize ())
         SX_THROW ("SxFileIO::write (SxString) unexpected amount of "
                   "bytes written");
      if (file.getLastWritten () != (uint64_t)fileData.getSize ())
         SX_THROW ("SxFileIO::getLastWritten unexpected result");
      if (file.getLastRead () != 0)
         SX_THROW ("SxFileIO::getLastRead unexpected result");

      file.close ();
      if (file.isOpen ())
         SX_THROW ("SxFileIO::close failed");

      file.open (filePath, "r");
      SxString readData;
      // Read into emtpy buffer
      nBytesRead = file.read (&readData);
      if (nBytesRead != 0)
         SX_THROW ("SxFileIO::getSize expected '0' as the size of buffer passed "
                   "is 0, but got '" + SxString(nBytesRead) + "' instead");

      readData.resize (fileData.getSize ());

      nBytesRead = file.read (&readData);
      if (nBytesRead != (uint64_t)fileData.getSize ())
         SX_THROW ("SxFileIO::read (SxString) unexpected amount of byte read");

      if (file.getSize () != (uint64_t)fileData.getSize ())
         SX_THROW ("SxFileIO::getSize unexpected result");

      if (readData != fileData)
         SX_THROW ("SxFileIO::read unexpected result");

      if (file.getLastRead () != (uint64_t)readData.getSize ())
         SX_THROW ("SxFileIO::getLastRead unexpected result");
      if (file.getLastWritten () != 0)
         SX_THROW ("SxFileIO::getLastWritten unexpected result");

      file.close ();

      if (file.isOpen ())
         SX_THROW ("SxFileIO::close failed");

      // ---  Mode checking
      file.open (filePath, "rb");
      if (  !file.hasMode (SxFileIO::Mode::ModeRead)
         ||  file.hasMode (SxFileIO::Mode::ModeWrite)
         ||  file.hasMode (SxFileIO::Mode::ModeAppend)
         ||  file.hasMode (SxFileIO::Mode::ModeCreate)
         || !file.hasMode (SxFileIO::Mode::ModeBinary))
         SX_THROW ("SxFileIO::hasMode unexpected result");
      file.close ();

      file.open (filePath, "r+");
      if (  !file.hasMode (SxFileIO::Mode::ModeRead)
         || !file.hasMode (SxFileIO::Mode::ModeWrite)
         ||  file.hasMode (SxFileIO::Mode::ModeAppend)
         ||  file.hasMode (SxFileIO::Mode::ModeCreate)
         ||  file.hasMode (SxFileIO::Mode::ModeBinary))
         SX_THROW ("SxFileIO::hasMode unexpected result");
      file.close ();

      file.open (filePath, "w");
      if (   file.hasMode (SxFileIO::Mode::ModeRead)
         || !file.hasMode (SxFileIO::Mode::ModeWrite)
         ||  file.hasMode (SxFileIO::Mode::ModeAppend)
         || !file.hasMode (SxFileIO::Mode::ModeCreate)
         ||  file.hasMode (SxFileIO::Mode::ModeBinary))
         SX_THROW ("SxFileIO::hasMode unexpected result");
      file.close ();

      file.open (filePath, "w+");
      if (  !file.hasMode (SxFileIO::Mode::ModeRead)
         || !file.hasMode (SxFileIO::Mode::ModeWrite)
         ||  file.hasMode (SxFileIO::Mode::ModeAppend)
         || !file.hasMode (SxFileIO::Mode::ModeCreate)
         ||  file.hasMode (SxFileIO::Mode::ModeBinary))
         SX_THROW ("SxFileIO::hasMode unexpected result");
      file.close ();

      file.open (filePath, "a");
      if (   file.hasMode (SxFileIO::Mode::ModeRead)
         || !file.hasMode (SxFileIO::Mode::ModeWrite)
         || !file.hasMode (SxFileIO::Mode::ModeAppend)
         || !file.hasMode (SxFileIO::Mode::ModeCreate)
         ||  file.hasMode (SxFileIO::Mode::ModeBinary))
         SX_THROW ("SxFileIO::hasMode unexpected result");
      file.close ();

      file.open (filePath, "a+");
      if (  !file.hasMode (SxFileIO::Mode::ModeRead)
         || !file.hasMode (SxFileIO::Mode::ModeWrite)
         || !file.hasMode (SxFileIO::Mode::ModeAppend)
         || !file.hasMode (SxFileIO::Mode::ModeCreate)
         || file.hasMode (SxFileIO::Mode::ModeBinary))
         SX_THROW ("SxFileIO::hasMode unexpected result");
      file.close ();

      SxFileIO::write (fileData, filePath, 0644);
      if (SxString::read (filePath) != fileData)
         SX_THROW ("static SxFileIO::write (SxString) invalid file content");

      SxString firstLine = "12345\n";
      SxString secondLine = "67890\n";
      SxString thirdLine = "abcde\n";
      firstLine.write (filePath);
      SxFileIO::appendToFile (secondLine, filePath);
      SxFileIO::appendToFile (thirdLine, filePath);

      file.open (filePath, "r+");
      SxString readLine1 = file.readLines (1);
      file.close ();

      file.open (filePath, "r+");
      SxString readLine2 = file.readLines (2);
      file.close ();

      file.open (filePath, "r+");
      SxString readLine3 = file.readLines (3);
      file.close ();

      if (readLine1 != firstLine)
         SX_THROW ("SxFileIO::readLines unexpected result");

      if (readLine2 != (firstLine + secondLine))
         SX_THROW ("SxFileIO::readLines unexpected result");

      if (readLine3 != (firstLine + secondLine + thirdLine))
         SX_THROW ("SxFileIO::readLines unexpected result");

      if (readLine1 != SxFileIO::readLines (filePath, 1))
         SX_THROW ("static SxFileIO::readLines unexpected result");

      if (readLine2 != SxFileIO::readLines (filePath, 2))
         SX_THROW ("static SxFileIO::readLines unexpected result");

      if (readLine3 != SxFileIO::readLines (filePath, 3))
         SX_THROW ("static SxFileIO::readLines unexpected result");

      file.open (filePath, "r+");
      SxString readTail1 = file.readTail (1);
      file.close ();

      file.open (filePath, "r+");
      SxString readTail2 = file.readTail (2);
      file.close ();

      file.open (filePath, "r+");
      SxString readTail3 = file.readTail (3);
      file.close ();

      if ((readTail1 + '\n') != ('\n' + thirdLine))
         SX_THROW ("SxFileIO::readTail unexpected result");

      if ((readTail2 + '\n') != ('\n' + secondLine + thirdLine))
         SX_THROW ("SxFileIO::readTail unexpected result");

      if ((readTail3 + '\n') != (firstLine + secondLine + thirdLine))
         SX_THROW ("SxFileIO::readTail unexpected result");

      if (readTail1 != SxFileIO::readTail (filePath, 1))
         SX_THROW ("static SxFileIO::readLines unexpected result");

      if (readTail2 != SxFileIO::readTail (filePath, 2))
         SX_THROW ("static SxFileIO::readLines unexpected result");

      if (readTail3 != SxFileIO::readTail (filePath, 3))
         SX_THROW ("static SxFileIO::readLines unexpected result");

      SxString binaryFileStr = rootPath + "data.bin";
      SxFileIO binaryFile;
      int binaryFileSize = 32;
      SxArray<unsigned char> binaryData(binaryFileSize);
      for (int i = 0; i < binaryFileSize; ++i)  {
         binaryData(i) = (unsigned char)(rand() % 256);
      }

      binaryFile.open (binaryFileStr, "wb", 0644);
      nBytesWritten = binaryFile.writeBuffer ((void *)binaryData.elements, (uint64_t)binaryFileSize);
      binaryFile.close ();

      if (nBytesWritten != (uint64_t)binaryFileSize)
         SX_THROW ("SxFileIO::write (void *) unexpected amount of bytes written");

      SxArray<unsigned char> readBinaryData (binaryFileSize);
      binaryFile.open (binaryFileStr, "rb", 0644);
      nBytesRead = binaryFile.readBuffer ((void *)readBinaryData.elements, (uint64_t)binaryFileSize);
      binaryFile.close ();

      if (  nBytesRead != (uint64_t)binaryFileSize
         || nBytesRead != nBytesWritten)
         SX_THROW ("SxFileIO::read (void *) unexpected amount of bytes read");

      if (binaryData.getSize () != readBinaryData.getSize ())
         SX_THROW ("SxFileIO::write/read (void *) content size is not same");

      for (int i = 0; i < binaryFileSize; i++)  {
         if (binaryData(i) != readBinaryData(i))
            SX_THROW ("SxFileIO::write/read (void *) content is not same");
      }

      SxArray<char> staticReadData = SxFileIO::readBinary (binaryFileStr,
                                                  (int64_t)binaryFileSize);
      if (staticReadData.getSize () != binaryFileSize)
         SX_THROW ("static SxFileIO::write (void *) unexpected amount of bytes written");

      for (int i = 0; i < binaryFileSize; i++)  {
         if (binaryData(i) != (unsigned char)staticReadData(i))
            SX_THROW ("static SxFileIO::readBinary (void *) content is not same");
      }


      // --- Single obj

      SxString objFileStr = rootPath + "obj.dat";
      SxFileIO objFile;
      objFile.open (objFileStr, "w", 0644);
      Vector3 vecOrig (1, 2, 3);
      uint64_t nElemsWritten = objFile.write (vecOrig);
      objFile.close ();
      if (nElemsWritten != 1)
         SX_THROW ("SxFileIO::write (const T &obj) unexpected "
                   "amount of objects written");

      Vector3 vecRead;
      objFile.open (objFileStr, "r+", 0644);
      uint64_t nElemsRead = objFile.read (&vecRead);
      objFile.close ();

      if (nElemsWritten != nElemsRead)
         SX_THROW ("SxFileIO::read (const T &obj) unexpected "
                   "amount of objects read");
      if (vecOrig != vecRead)
         SX_THROW ("SxFileIO::write/read (const T &obj) content is not same");

      // --- Array of Object
      SxArray<Vector3> vecsOrig;
      vecsOrig.resize (3);
      vecsOrig(0) = Vector3 (1, 1, 1);
      vecsOrig(1) = Vector3 (2, 2, 2);
      vecsOrig(2) = Vector3 (3, 3, 3);

      // Remove old file
      SxFSAction::rm (objFileStr);

      objFile.open (objFileStr, "w", 0644);
      nElemsWritten = objFile.write (vecsOrig, 3);
      objFile.close ();
      if (nElemsWritten != (uint64_t)vecsOrig.getSize ())
        SX_THROW ("SxFileIO::write (const SxArray<T> &) unexpected "
                  "amount of elements written");

      SxArray<Vector3> readVecs;
      readVecs.resize ((ssize_t)nElemsWritten);
      objFile.open (objFileStr, "r+", 0644);
      nElemsRead = objFile.read (&readVecs, (int64_t)nElemsWritten);
      objFile.close ();
      if (nElemsWritten != nElemsRead)
         SX_THROW ("SxFileIO::read (const SxArray<T> &) unexpected "
                  "amount of elements read");
      if (!(vecsOrig == readVecs))
         SX_THROW ("SxFileIO::write/read (const SxArray<T> &) "
                   "content is not same");

      SxString operatorFileStr = rootPath + "objOperator.dat";
      SxFileIO operatorFile;
      operatorFile.open (operatorFileStr, "w");
      operatorFile << vecsOrig;
      operatorFile.close ();
      if (  SxFileInfo (operatorFileStr).getSize ()
         != (int64_t)((unsigned long)vecsOrig.getSize () * sizeof (Vector3)))
         SX_THROW ("SxFileIO &operator<< created file of unexpected size");

      SxArray<Vector3> readOperatorData;
      readOperatorData.resize (vecsOrig.getSize ());
      operatorFile.open (operatorFileStr, "r");
      operatorFile >> readOperatorData;
      operatorFile.close ();
      if (readOperatorData != vecsOrig)
         SX_THROW ("SxFileIO &operator>> read content different from original");

      operatorFile.open (operatorFileStr, "w");
      (ostream&)operatorFile << fileData << vecOrig;
      operatorFile.close ();

      SxString operatorReadData = SxString::read (operatorFileStr);
      if (operatorReadData != (fileData + vecOrig.toString ()))
         SX_THROW ("SxFileIO ostream &operator<< written content "
                   "is different from original");

      SxArray<char> operatorIStreamData;
      operatorIStreamData.resize (operatorReadData.getSize () + 1);
      operatorFile.open (operatorFileStr, "r");
      ((istream&)operatorFile) >> operatorIStreamData.elements;
      operatorFile.close ();
      if (SxString(operatorIStreamData) != operatorReadData)
         SX_THROW ("SxFileIO istream &operator>> read content "
                   "is different from original");

      SxString printfFileStr = rootPath + "printfFile.dat";
      SxString fmt = "%d%s";
      SxFileIO printfFile;
      printfFile.open (printfFileStr, "w");
      printfFile.printf (fmt.getElems (), fileData.getSize (), fileData.elements);
      printfFile.close ();
      SxString sxprintRes = SxString::sprintf (fmt.getElems (), fileData.getSize (), fileData.elements);
      if (sxprintRes != SxString::read (printfFileStr))
         SX_THROW ("SxFileIO printf write content is different from original");

      SxString scanfStr;
      scanfStr.resize (fileData.getSize () + 1);
      int scanfInt = 0;
      printfFile.open (printfFileStr, "r+");
      printfFile.scanf (fmt.getElems (), &scanfInt, scanfStr.elements);
      printfFile.close ();

      if (scanfInt != fileData.getSize ())
         SX_THROW ("SxFileIO scanf [int] content is different from original");

      // An additional '/0' may be appended
      if ((scanfStr.getSize () - 1) == fileData.getSize ())
         scanfStr.resize (fileData.getSize (), true);
      if (scanfStr != fileData)
         SX_THROW ("SxFileIO scanf [char *] content is different from original");


      // --- Seek
      file.open (filePath, "r+");
      if (file.getOffset () != 0)
         SX_THROW ("SxFileIO::getOffset unexpected result");
      file.seek (1, SxFileIO::Seek::CUR);
      if (file.getOffset () != 1)
         SX_THROW ("SxFileIO::getOffset unexpected result");
      file.seek (1, SxFileIO::Seek::CUR);
      if (file.getOffset () != 2)
         SX_THROW ("SxFileIO::getOffset unexpected result");
      file.seek (3, SxFileIO::Seek::BEG);
      if (file.getOffset () != 3)
         SX_THROW ("SxFileIO::getOffset unexpected result");
      file.seek (-2, SxFileIO::Seek::END);
      if (file.getOffset () != file.getSize () - 2)
         SX_THROW ("SxFileIO::getOffset unexpected result");

      file.close ();

      // tell and getOffSet is same?

      SxFileIO flushFile;
      flushFile.open (filePath, "w+");
      flushFile.write (fileData);
      if (SxFileInfo(filePath).getSize () != 0)
         SX_THROW ("SxFileIO::write already flushed data.");
      flushFile.flush ();
      if (SxFileInfo(filePath).getSize () != fileData.getSize ())
         SX_THROW ("SxFileIO::flush did not flush the content to file.");
      flushFile.close ();

      /*
      // syncDir in order to update directory after file creation
      static void syncDir (const SxString &path);
      */

   } catch (SxException e) {
      SX_THROW (e, "SxFS_FileIO Failed.");
   }
}
