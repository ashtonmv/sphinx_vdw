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

#include <SxFileParser.h>
#include <SxStack.h>
#include <stdarg.h>

SxFileParser::SxFileParser ()
   : fp(NULL),
     line(-1),
     fpos(-1),
     verbose (false)
{
   // empty
}

SxFileParser::SxFileParser (const SxString &fileName)
   : fp(NULL), fpos(-1), verbose (false)
{
   open (fileName);
}

SxFileParser::~SxFileParser ()
{
   if (fp) fclose (fp);
}

void SxFileParser::open (const SxString &fileName)
{
   SX_CHECK (!fp); // please close previous file first
   fp = fopen (fileName.ascii (), "r");
   if (!fp)  {
      cout << "Can't open file '" << fileName << "' for reading: "
           << sxstrerror() << endl;
      SX_QUIT;
   }
   name = fileName;
   line = 1;
   fpos = 0;
}

void SxFileParser::close ()
{
   SX_CHECK (fp);
   fclose (fp);
   fp = NULL;
   name = "";
   line = -1;
   fpos = -1;
}

void SxFileParser::updateLine ()
{
   SX_CHECK (fp);
   long currentPos = ftell (fp);
   if (currentPos <= fpos) return;
   // jump to last known line position
   fseek (fp, fpos, SEEK_SET);
   char buffer[1024];
   while (currentPos > fpos)  {
      size_t n = size_t (currentPos - fpos);
      if (n > 1024) n = 1024;
#ifndef NDEBUG
      size_t nn =
#endif
      fread (buffer, sizeof(char), n, fp);
      SX_CHECK (nn == n, nn, n);
      for (size_t i = 0; i < n; ++i)
         if (buffer[i] == '\n') line++;
      fpos += (long)n;
   }
   SX_CHECK (currentPos == ftell (fp), currentPos, ftell(fp));
   fpos = currentPos;
}

void SxFileParser::nextLine (int nLine)
{
   SX_CHECK (nLine > 0, nLine);
   updateLine ();
   while (!feof (fp))  {
     if (fgetc (fp) == '\n')  {
        line++;
        nLine--;
        if (nLine == 0) break;
     }
   }
   fpos = ftell (fp);
   noEof ();
}

void SxFileParser::where ()
{
   cout << "while reading ";
   if (parserTopic.getSize () > 0) cout << parserTopic << " from ";
   cout << name << " at line " << getLineNo ();
}

void SxFileParser::noEof ()
{
   SX_CHECK (fp);
   if (feof (fp))  {
      cout << "Unexpected end of file ";
      where ();
      cout << endl;
      SX_QUIT;
   }
}

double SxFileParser::getDouble ()
{
   noEof ();
   double x;
   int n = fscanf (fp, "%lf", &x);
   if (n != 1)  {
      noEof ();
      cout << "Failed to read real number ";
      where ();
      cout << "." << endl;
      cout << "Rest of line reads: " << getLine () << endl;
      SX_QUIT;
   }
   return x;
}

int SxFileParser::getInt ()
{
   noEof ();
   int x;
   int n = fscanf (fp, "%d", &x);
   if (n != 1)  {
      noEof ();
      cout << "Failed to read integer number ";
      where ();
      cout << "." << endl;
      cout << "Rest of line reads: " << getLine () << endl;
      SX_QUIT;
   }
   return x;
}

long SxFileParser::getLong ()
{
   noEof ();
   long x;
   int n = fscanf (fp, "%ld", &x);
   if (n != 1)  {
      noEof ();
      cout << "Failed to read integer number ";
      where ();
      cout << "." << endl;
      cout << "Rest of line reads: " << getLine () << endl;
      SX_QUIT;
   }
   return x;
}

SxString SxFileParser::getLine ()
{
   SX_CHECK (fp);
   noEof ();
   updateLine ();
   SxString res;
   char buffer[1024];

   // make sure that buffer is empty if fgets doesn't write to it
   buffer[0]='\0';
   while (fgets (buffer, 1024, fp))  {
      res += buffer;
      SX_CHECK (res.getSize () > 0);
      if (res(res.getSize () - 1) == '\n') {
         line++;
         fpos = ftell (fp);
         break;
      }
      if (feof (fp))  {
         cout << "Warning: last line of " << name << " has no newline" << endl;
         break;
      }
      // make sure that buffer is empty if fgets doesn't write to it
      buffer[0]='\0';
   }
   return res;
}

void SxFileParser::skipWhite ()
{
   SX_CHECK (fp);
   fscanf (fp, " ");
}

void SxFileParser::read (const char* what)
{
   noEof ();
   if (what[0] != ' ') skipWhite ();
   long pos = ftell (fp);
   for (ssize_t i = 0; what[i] != '\0'; i++)  {
      if (fgetc (fp) != what[i])  {
         fseek (fp, pos, SEEK_SET);
         cout << "Unexpected content ";
         where ();
         cout << endl;
         cout << "Expected   '" << what << "'" << endl;
         SxString currentLine = getLine ();
         ssize_t last = currentLine.getSize () - 1;
         if (currentLine(last) == '\n') currentLine(last) = 0;
         cout << "but found: '" << currentLine << "'" << endl;
         SX_QUIT;
      }
   }
}

void SxFileParser::readUntil (char what)
{
   SX_CHECK (fp);
   while (fgetc (fp) != what && !feof (fp))
   {
      // empty
   }
   noEof ();
}

bool SxFileParser::reads (const char* what)
{
   noEof ();
   if (what[0] != ' ') skipWhite ();
   long pos = ftell (fp);
   for (ssize_t i = 0; what[i] != '\0'; i++)  {
      if (fgetc (fp) != what[i])  {
         fseek (fp, pos, SEEK_SET);
         return false;
      }
   }
   return true;
}

SxArray<double> SxFileParser::getVector ()
{
   SX_CHECK (fp);
   noEof ();
   SxStack<double> res;
   double x;
   while (fscanf (fp, " %lf", &x) == 1)  {
      res << x;
   }
   return SxArray<double> (res);
}

SxArray<double> SxFileParser::getVector (ssize_t n)
{
   SX_CHECK (fp);
   SX_CHECK (n > 0, n);
   noEof ();
   SxArray<double> res(n);
   for (ssize_t i = 0; i < n; ++i)  {
      if (fscanf (fp, " %lf", &res(i)) != 1)  {
         where ();
         cout << endl;
         cout << "Tried to read vector with " << n << " elements," << endl;
         cout << "but read only " << i << " elements, and then found";
         SxString restOfLine =  getLine ();
         if (restOfLine.getSize () == 0 && feof(fp))
            cout << " nothing (end of file).";
         else
            cout << ": " << restOfLine;
         cout << endl;
         SX_QUIT;
      }
   }
   return res;
}

void SxFileParser::find (const char *what)
{
   SX_CHECK (what);
   SX_CHECK (what[0] != '\0');
   long pos = ftell (fp);
   while (!feof(fp))  {
      if (fgetc (fp) == what[0])  {
         if (what[1] == '\0') return;
         fseek (fp, -1, SEEK_CUR);
         if (reads (what)) return;
         fgetc (fp); // not found, move to next char
      }
   }
   cout << "Unexpected end of file while searching for '"
        << what << "'" << endl
        << "in " << name << endl;
   fseek (fp, pos, SEEK_SET);
   cout << "Started search in line " << getLineNo () << " at '"
        << getLine () << "'." << endl;
   SX_QUIT;
}
