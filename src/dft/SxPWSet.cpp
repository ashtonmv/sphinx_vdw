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

#include <SxPWSet.h>

SxPWSet::SxPWSet () 
   : SxPsiSet (PW),
     gkBasisPtr() 
{ 
   /* empty */ 
}
      
SxPWSet::SxPWSet (SxPtr<SxGkBasis> inPtr)
   : SxPsiSet (PW),
     gkBasisPtr(inPtr) 
{ 
   /* empty */ 
}

SxGkBasis &SxPWSet::getGkBasis () const 
{
   SX_CHECK (gkBasisPtr);
   return *gkBasisPtr;
}

const SxGBasis::TPsi& SxPWSet::operator() (int, int) const
{
   // may not be implemented
   SX_EXIT;
   return *(SxGBasis::TPsi *)(NULL);
}

SxGBasis::TPsi& SxPWSet::operator() (int, int)
{
   // may not be implemented
   SX_EXIT;
   return *(SxGBasis::TPsi *)(NULL);
}

PsiG SxPWSet::getBlock (int, int, int, int)
{
   // may not be implemented
   SX_EXIT;
   return PsiG ();
}

SxPtr<SxPWSet> SxPWSet::getNew () const
{
   SX_EXIT;
   return SxPtr<SxPWSet> ();
}

#include <SxFermi.h>
void SxPWSet::writeWavesFile (const SxString          &filename,
                              const SxFermi           &fermi,
                              const SxAtomicStructure &structure,
                              bool                    /*writeParserBuffer*/) const
{
   SX_ALLOC_CACHE;
   const SxGkBasis &Gk = getGkBasis ();
   // --- write waves to disc
   double wavesSize;
   wavesSize  = getNStates ();
   wavesSize *= getNSpin ();
   wavesSize *= getNk ();
   wavesSize *= Gk.nGkMax;
   if (Gk.gBasis)
      wavesSize *= Gk.gBasis->getNComp();
   wavesSize *= sizeof(PrecCoeffG);
   
   SxBinIO::Mode ncFormat;
#ifdef USE_PARALLEL_NETCDF4
   ncFormat = SxBinIO::BINARY_WRITE_PARALLEL;    // parallel NetCDF4/HDF5 format
#else
   if (wavesSize > 1e9 /* = 1 GB */)
      ncFormat = SxBinIO::BINARY_WRITE_LARGE;    // large (64bit) format
   else
      ncFormat = SxBinIO::BINARY_WRITE_ONLY;     // classic format
#endif


   /* The version works as follows:
      Before decimal point: crucial change in format. Wave function files
      are not compatible at all. Changes in this digit require an
      implementation of the transfer in sxwavetransfer. Backward 
      compatibility code should be moved from read functions to sxwavetransfer.

      First digit: major change. Read functions may reject older formats if
      they differ in this digit. A workaround should be implemented in
      sxwavetransfer.

      2nd and 3rd digit: minor change. We use two digits here, which gives
      us space for for 100 minor changes per major change.

      Examples:
      1.000 SFHIngX 1.0 format
      1.1   SFHIngX 1.1 format
      1.101 first modification of SFHIngX 1.1 format
      1.102 second modification of SFHIngX 1.1 format
      1.234 34th change to 1.2 format, which is substantially different from
            SFHIngX 1.1 format
      2.0   Upgrade to the next generation binary format.

   */
   double version = 1.001;
   try  {
      SxBinIO io (filename, ncFormat);

      io.write ("formatVersion", version);
//    if (writeParserBuffer)
//       io.write ("input", SxParser_buffer);
      fermi.write (io);
      Gk.write (io);
      structure.write (io);
      write (io);

      io.setMode (SxBinIO::WRITE_DATA);
      io.write ("formatVersion", version);
//      if (writeParserBuffer)
//         io.write ("input", SxParser_buffer);
      fermi.write (io);
      Gk.write (io);
      structure.write (io);
      write (io);

      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}

void SxPWSet::write (SxBinIO &) const
{
   SX_EXIT;
}
