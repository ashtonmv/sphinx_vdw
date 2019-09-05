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

#include <SxCLI.h>
#include <SxPW.h>
#include <SxBinIO.h>
#include <SxFermi.h>
#include <SxPseudoPot.h>
#include <SxPerturbK.h>
#include <netcdf.h>


int main (int argc, char **argv)
{

   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   SxCLI cli (argc,argv);
   cli.authors = "C. Freysoldt";
   cli.preUsageMessage = "This add-ons computes the matrix elements "
                         "of the position operator. NO SPIN!";

   SxString wavesFile = cli.option ("-w|--waves","file","waves file")
                        .toString("waves.sxb");

   SxString outFile = cli.option ("-o|--output", "file", 
                                  "output file (ETSF/netCDF format)")
                     .toString ("rme-etsf.nc");

   int nL = cli.option("-l|--left","number","number of left states")
            .toInt(-1);
   cli.last ().defaultValue = "default: number of valence states";
   int nR = cli.option("-r|--right","number","number of right states")
            .toInt(-1);
   cli.last ().defaultValue = "default: number of conduction states";

   int offsetL = cli.option("--offset-left","number","left state offset")
                 .toInt(1) - 1;
   int offsetR = cli.option("--offset-right","number","right state offset")
                 .toInt(-1) - 1;
   cli.last ().defaultValue = "default: first conduction state";

   cli.newGroup ("nonlocal pseudopotential");
   SxString inFile = cli.option ("-i|--input", "file", "S/PHI/nX input file")
                     .toString ("input.sx");

   bool nonLocal 
      = cli.option ("--nonlocal","add non-local pseudopotential contributions")
        .toBool ();

   int blockSize
      = cli.option ("--blocksize","set blocksize. Smaller blocksizes may reduce "
                    "the memory consumption, larger blocksizes may speed up the "
                    " <phi|psi> part").toInt (64);


   cli.version ("0.9");
   cli.finalize ();

   SxPW waves(wavesFile, SxPW::ReadOnDemand);
   SxPtr<SxGkBasis> gkBasisPtr;

   SxAtomicStructure structure;

   try {
      SxBinIO io (wavesFile, SxBinIO::BINARY_READ_ONLY);

      if ((nL == -1 || nR == -1 || offsetR == -1))  {
         SxFermi fermi;
         fermi.read (io);
         if (! fermi.isSemiconductor () ) {
            cout << SxString(
                    "This is not a semiconductor. Can't detect occupied/"
                    "unoccupied states automatically. Please give -l, -r,"
                    " and --offset-right explicitly.").wrap ()
                 << endl;
            SX_QUIT;
         }

         if (nL == -1) nL = fermi.getNValenceBands (0,0);
         if (nR == -1) nR = fermi.getNStates(0) - fermi.getNValenceBands (0,0);
         if (offsetR == -1) offsetR = fermi.getNValenceBands (0,0);
      }

      gkBasisPtr = SxPtr<SxGkBasis>::create (io, false);
      waves.setGkBasisPtr (gkBasisPtr);
      structure.read(io);

      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   SxGkBasis &gkBasis = *gkBasisPtr;

   int nk = waves.getNk ();
   int iSpin = 0;

   SxDiracMat<TPrecCoeffG> npn(nk,3);

   SxPerturbK kp;
   if (nonLocal)  {
      SxParser parser;
      SxParser::Table table = parser.read(inFile);
      SxPseudoPot psPot(&*table);
      gkBasis.changeTau(structure);
      kp.set (psPot, gkBasis, structure);
      kp.blockSize = blockSize;
   }
   
   // --- create netcdf file
   int ncId, cartDirId, nLId, nRId, nkId, complexId, meId;
#define DOONCE(what) do { what ; } while (0) 
   
#define SX_CHECK_NC_CALL(call) DOONCE( \
   int res = call; \
   if (res != NC_NOERR) {\
      cout << "netCDF error occurred: " << nc_strerror(res) << endl; \
      SX_EXIT; \
   } )

#define NC_DIM_DEF(name,val,id) \
   SX_CHECK_NC_CALL(nc_def_dim(ncId,(name),(val),&(id)))
   
   // create output file
   SX_CHECK_NC_CALL( nc_create(outFile.ascii (), NC_WRITE, &ncId) );
   // define dimensions
   NC_DIM_DEF("real_or_complex",                2, complexId);
   NC_DIM_DEF("number_of_cartesian_directions", 3, cartDirId);
   NC_DIM_DEF("number_of_kpoints",             nk, nkId);
   NC_DIM_DEF("number_of_left_states",         nL, nLId);
   NC_DIM_DEF("number_of_right_states",        nR, nRId);
   
   // define variable
   int meDims[5];
   size_t meStart[5], meCount[5];
   for (size_t *ptr = meStart; ptr < meStart+5; ) *ptr++ = 0; 

   meDims[0] = nkId;      meCount[0] = 1; // k-point-wise write
   meDims[1] = nRId;      meCount[1] = nR;
   meDims[2] = nLId;      meCount[2] = nL;
   meDims[3] = cartDirId; meCount[3] = 3;
   meDims[4] = complexId; meCount[4] = 2;
   
   SX_CHECK_NC_CALL( 
      nc_def_var(ncId, "position_matrix_elements", NC_DOUBLE, 5, meDims, &meId)
   );

   // define attributes
   int offsetFileR = offsetR + 1,
       offsetFileL = offsetL + 1;
   SX_CHECK_NC_CALL( 
      nc_put_att_int(ncId, meId, "offset_right", NC_INT, 1, &offsetFileR)
   );
   SX_CHECK_NC_CALL( 
      nc_put_att_int(ncId, meId, "offset_left", NC_INT, 1, &offsetFileL)
   );
   SX_CHECK_NC_CALL( 
      nc_put_att_text(ncId, meId, "offset_hint", 19, "offset starts at 1")
   );

   // end definition mode
   SX_CHECK_NC_CALL( 
      nc_enddef(ncId);
   );

   
   SxDiracMat<TPrecCoeffG> kpMatElem;
   PsiG wavesL, wavesR;
   int ng;
   
   for (int ik = 0; ik < nk; ik++)  {
      (cout << "ik = " << ik << endl).flush ();
      
      // get left and right states
      wavesL = PsiG (); wavesR = PsiG (); // free reference
      ng = gkBasis(ik).ng;
      wavesL = waves(iSpin, ik)(SxIdx(offsetL*ng, (offsetL+nL)*ng - 1));
      wavesL.reshape(ng, nL);
      wavesR = waves(iSpin, ik)(SxIdx(offsetR*ng, (offsetR+nR)*ng - 1));
      wavesR.reshape(ng, nR);

      // compute matrix elements
      (cout << "Computing data..." << endl).flush ();
      kpMatElem = kp.getMatrixElements(wavesL, wavesR).transpose (); 
      
      // --- write matrix elements
      meStart[0] = ik;
      (cout << "Writing data..." << endl).flush ();
      SX_CHECK_NC_CALL (
         nc_put_vara_double(ncId, meId, meStart, meCount, 
                            (double *)kpMatElem.elements)
      );
            
   } // ik

   SX_CHECK_NC_CALL ( nc_close(ncId) );

   if (nonLocal)
      kp.printTimer ();

   return 0;
   
}

