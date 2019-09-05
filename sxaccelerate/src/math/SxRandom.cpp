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
// This is a random namber generator ("universal")
// Source: R. Wieczorkowski, R. Zielincki "Komputerowe generatory
// licyb losowych" WNT, Warszawa 1997, page 39-40
//----------------------------------------------------------------------------

#include <SxRandom.h>
#include <stdio.h>

double SxRandom::uu[97];
long int SxRandom::ip;
long int SxRandom::jp;
double SxRandom::cc;
double SxRandom::cd;
double SxRandom::cm;



void initRandom (int i, int j, int k, int l)
{
   int ii, jj, m, wi, wj, wk, wl;
   double s, t;

   SxRandom::ip = 97;
   SxRandom::jp = 33;
   SxRandom::cc = 362436.   / 16777216.;
   SxRandom::cd = 7654321.  / 16777216.;
   SxRandom::cm = 16777213. / 16777216.;

   wi = i;
   wj = j;
   wk = k;
   wl = l;
   for (ii = 0; ii < 97; ii++)
   {
      s = 0.;
      t = 0.5;
      for(jj = 1; jj <= 24; jj++)
      {
         m = (((wi * wj) % 179) * wk) % 179;
         wi = wj;
         wj = wk;
         wk = m;
         wl = (53 * wl + 1) % 169;
         if ((wl * m) % 64 >= 32)
            s += t;
         t *= 0.5;

      }
      SxRandom::uu[ii] = s;
   }
}

double SxRandom::get ()
{
   double tmp;

    tmp = uu[ip - 1] - uu[jp - 1];

   if(tmp < 0.)  tmp += 1.;
   uu[ip - 1] = tmp;

   ip--;
   if(ip == 0) ip = 97;

   jp--;
   if(jp == 0) jp = 97;

   cc -= cd;
   if(cc < 0.) cc += cm;

   tmp -= cc;
   if(tmp < 0.) tmp += 1;

   return tmp;
}

