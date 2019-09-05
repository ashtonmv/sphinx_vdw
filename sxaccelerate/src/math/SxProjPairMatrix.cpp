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

#include <SxProjPairMatrix.h>

void SxProjPairMatrix::append (const Vector& s, const Vector& t,
                               double prefactor,
                               double ass, double ast, double ats, double att)
{
   // --- add new projector pair
   sk << s;
   tk << t;
   akss << (prefactor * ass);
   akst << (prefactor * ast);
   akts << (prefactor * ats);
   aktt << (prefactor * att);
   if (sk.getSize () > nProj)  {
      // remove oldest projector pair
      sk.removeFirst ();
      tk.removeFirst ();
      akss.removeFirst ();
      akst.removeFirst ();
      akts.removeFirst ();
      aktt.removeFirst ();
   }
}

SxProjPairMatrix::Vector
SxProjPairMatrix::operator^(const Vector& in) const
{
   ssize_t n = sk.getSize ();
   Vector res = diag * in;
   double bs, bt;
   SxList<Vector>::ConstIterator 
      sIt = sk.begin (),
      tIt = tk.begin ();
   SxList<double>::ConstIterator
      ss = akss.begin (),
      st = akst.begin (),
      ts = akts.begin (),
      tt = aktt.begin (); 
   // --- loop over projector pairs
   for (int i = 0; i < n; ++i, ++sIt,++tIt,++ss,++st,++ts,++tt)
   {
      
      // --- right projection
      bs = dot(*sIt, in);
      bt = dot(*tIt, in);
      
      // --- left projection
      // res += s * (ass * bs + ast * bt);
      // res += t * (ats * bs + att * bt);
      res.plus_assign_ax (*ss * bs + *st * bt, *sIt);
      res.plus_assign_ax (*ts * bs + *tt * bt, *tIt);
   }
   return res;
}
