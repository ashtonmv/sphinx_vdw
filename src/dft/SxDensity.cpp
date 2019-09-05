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

#include <SxDensity.h>
#include <SxBinIO.h>

SxDensity::SxDensity (const SxPtr<SxDensity> &inPtr)
: ptr(inPtr)
{
   // empty
}

void SxDensity::operator= (const SxDensity &in)
{
   // note: you must not assign a derived density to a general density
   // use an indirect density instead (e.g. via getCopy)
   SX_CHECK (in.ptr);
   ptr = in.ptr;
}

void SxDensity::operator += (const SxDensity &x)
{
   // indirection
   SX_CHECK (ptr);
   ptr->operator+= (x);
}

void SxDensity::operator -= (const SxDensity &x)
{
   // indirection
   SX_CHECK (ptr);
   ptr->operator-= (x);
}


void SxDensity::plus_assign_ax (double a, const SxDensity &x)
{
   // indirection
   SX_CHECK (ptr);
   return ptr->plus_assign_ax (a, x);
}

void SxDensity::plus_assign_aspin (double a, const SxDensity &x)
{
   // indirection
   SX_CHECK (ptr);
   return ptr->plus_assign_aspin (a, x);
}

double SxDensity::operator| (const SxDensity &x) const
{
   // indirection
   SX_CHECK (ptr);
   return ptr->operator| (x);
}

double SxDensity::normSqr () const
{
   // indirection
   if (ptr)
      return ptr->normSqr ();
   return (*this | *this);
}

SxDensity SxDensity::operator- (const SxDensity &x) const
{
   // indirection
   SX_CHECK (ptr);
   return ptr->operator- (x);
}

SxDensity SxDensity::getCopy () const
{
   // indirection
   SX_CHECK (ptr);
   return ptr->getCopy ();
}

SxDensity SxDensity::spin () const
{
   // indirection
   SX_CHECK (ptr);
   return ptr->spin ();
}

bool SxDensity::hasSpin () const
{
   // indirection
   SX_CHECK (ptr);
   return ptr->hasSpin ();
}

void SxDensity::renormalize ()
{
   // indirection
   SX_CHECK (ptr);
   ptr->renormalize ();
}

void SxDensity::readRho (const SxString &file)
{
   if (ptr)  {
      // indirection
      ptr->readRho (file);
   } else {
      SxBinIO io;
      try {
         io.open (file, SxBinIO::BINARY_READ_ONLY);
         readRho (io);
         io.close ();
      } catch (SxException e) {
         e.print ();
         SX_EXIT;
      }
   }
}

void SxDensity::writeRho (const SxString &file) const
{
   if (ptr)  {
      // indirection
      ptr->writeRho (file);
   } else {
      SxBinIO io;
      try {
         io.open (file, SxBinIO::BINARY_WRITE_ONLY);
         writeRho(io);
         io.setMode (SxBinIO::WRITE_DATA);
         writeRho(io);
         io.close ();
      } catch (SxException e) {
         e.print ();
         SX_EXIT;
      }
   }
}

void SxDensity::readRho (const SxBinIO &file)
{
   // indirection
   SX_CHECK(ptr);
   ptr->readRho (file);
}

void SxDensity::writeRho (SxBinIO &file) const
{
   // indirection
   SX_CHECK(ptr);
   ptr->writeRho (file);
}

void SxDensity::syncMPI ()
{
   SX_CHECK (ptr);
   ptr->syncMPI ();
}
