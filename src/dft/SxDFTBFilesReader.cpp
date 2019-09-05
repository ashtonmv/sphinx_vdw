// ---------------------------------------------------------------------------
//
//      The ab-initio based multiscale library
//
//                  S / P H I / n X
//
//      Copyright:  Max-Planck-Institute for Iron Research
//                  40237 Duesseldorf, Germany
//
//      Contact:    https://sxlib.mpie.de
//      Authors:    see sphinx/AUTHORS
//      License:    see sphinx/LICENSE
//
// ---------------------------------------------------------------------------

#include <SxDFTBFilesReader.h>

double SxDFTBFilesReader::getNextReal ()
{
   SX_CHECK (fpIn);
   if (occurence > 0)  {
      occurence--;
      return val;
   }
   // read next word
   char buffer[1024];
   if (feof(fpIn))  {
      cout << "\n Unexpected end of file" << endl;
      SX_EXIT;
   }
   int nScan = fscanf (fpIn, "%1023s", buffer);
   if (nScan == 0)  {
      cout << "\n Unexpected end of file" << endl;
      SX_EXIT;
   }
   if (sscanf(buffer, "%d*%lf", &occurence, &val) == 2)  {
      if (occurence <= 0)  {
         cout << "\n Illegal format N*x. N=" << occurence << endl;
         SX_EXIT;
      }
      occurence--;
      return val;
   } else {
      char c;
      occurence = 0;
      if (sscanf(buffer, "%lf%c", &val, &c) != 1)  {
         cout << "\n Cannot convert to number: " << buffer << endl;
         SX_EXIT;
      }
      return val;
   }
}

