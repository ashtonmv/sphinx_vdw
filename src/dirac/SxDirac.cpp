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

#include <SxDirac.h>
#include <SxLoopMPI.h>

void writePlot (const SxString &fileName,
                const SxDiracVec<Double> &x,
                const SxDiracVec<Double> &y1,
                const SxDiracVec<Double> &y2,
                const SxDiracVec<Double> &y3)
{ 
   if (SxLoopMPI::me() != 0) return;
   FILE *out = fopen (fileName.ascii (), "w");
   if (!out)  {
      cout << "Failed to open " << fileName << " for writing" << endl;
      SX_EXIT;
   }
   SX_CHECK (y1.getSize () == x.getSize (), x.getSize (), y1.getSize ());
   if (y2.getSize () == 0)  {
      for (int i = 0; i < x.getSize (); ++i)
         fprintf (out, "%.12f\t%.16f\n", x(i), y1(i));
   } else if (y3.getSize () == 0) {
      SX_CHECK (y2.getSize () == x.getSize (), x.getSize (), y2.getSize ());
      for (int i = 0; i < x.getSize (); ++i)
         fprintf (out, "%.12f\t%.16f\t%.16f\n", x(i), y1(i), y2(i));
   } else {
      SX_CHECK (y2.getSize () == x.getSize (), x.getSize (), y2.getSize ());
      SX_CHECK (y3.getSize () == x.getSize (), x.getSize (), y3.getSize ());
      for (int i = 0; i < x.getSize (); ++i)
         fprintf (out, "%.12f\t%.16f\t%.16f\t%.16f\n", x(i), y1(i), y2(i), y3(i));
   }
   fclose (out);
}

void writePlot (const SxString &fileName,
                double dx,
                const SxDiracVec<Double> &y1,
                const SxDiracVec<Double> &y2)
{
   if (SxLoopMPI::me() != 0) return;
   FILE *out = fopen (fileName.ascii (), "w");
   if (!out)  {
      cout << "Failed to open " << fileName << " for writing" << endl;
      SX_EXIT;
   }
   if (y2.getSize () == 0)  {
      for (int i = 0; i < y1.getSize (); ++i)
         fprintf (out, "%f\t%.16f\n", i * dx, y1(i));
   } else {
      SX_CHECK(y2.getSize () == y1.getSize (), y1.getSize (), y2.getSize ());
      for (int i = 0; i < y1.getSize (); ++i)
         fprintf (out, "%f\t%.16f\t%.16f\n", i * dx, y1(i), y2(i));
   }
   fclose (out);
}

