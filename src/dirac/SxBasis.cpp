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

#include <SxBasis.h>
#include <SxDirac.h>

bool SxBasis::registerBasis (const SxBasis &basis) const
{
   if (isRegistered (&basis)) return false;
   ssize_t n = registeredBases.getSize ();
   registeredBases.resize (n + 1, true);
   registeredBases(n) = &basis;
   basis.registerUnknown (*this);
   return true;
}

SxBasis::~SxBasis()  {
   // if this fails, you forgot to put deregisterAll () in your
   // destructor
   SX_CHECK (registeredBases.getSize () == 0);
}

void SxBasis::deregisterAll ()
{
   for (ssize_t n = registeredBases.getSize (); n > 0; --n)
      deregister (registeredBases(n-1));
   SX_CHECK (registeredBases.getSize () == 0);
}

void SxBasis::print () const
{
   sxprintf ("SxBasis: Type = %s\n", getType ().ascii ());
   sxprintf ("registrated peers: \n");
   for (int i = 0; i < registeredBases.getSize (); ++i)
      sxprintf ("   %s\n", registeredBases(i)->getType ().ascii ());
}

void SxBasis::deregisterBasis (const SxBasis *basis) const
{
   int i;
   // find basis
   for (i = 0; i < registeredBases.getSize (); ++i)
      if (registeredBases(i) == basis) break;
   if (i == registeredBases.getSize ()) return; // nothing to do
   for (++i ; i < registeredBases.getSize (); ++i)
      registeredBases(i-1) = registeredBases(i);
   registeredBases.resize (registeredBases.getSize () - 1, true);
   basis->deregister (this);
}

double SxBasis::scalarProduct (const SxDiracVec<Double> &x,
                               const SxDiracVec<Double> &y) const
{
   SX_CHECK (x.getSize () == getNElements (),
             x.getSize (), getNElements ());
   if (&x == &y) return x.normSqr ();
   return dot (x, y);
}

SxComplex16 SxBasis::scalarProduct (const SxDiracVec<Complex16> &x,
                                    const SxDiracVec<Complex16> &y) const
{
   SX_CHECK (x.getSize () == getNElements (),
             x.getSize (), getNElements ());
   if (&x == &y) return x.normSqr ();
   return dot (x, y);
}

void SxBasis::projectionFailed (const char *basis, const char *types) const
{
   sxprintf ("Projector from %s to %s not yet implemented. %s\n",
             getType ().ascii (), basis, types);
   SX_EXIT;
}

