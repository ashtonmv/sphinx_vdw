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

#ifndef _SX_WAVES_JOIN_H_
#define _SX_WAVES_JOIN_H_

#include <SxFermi.h>
#include <SxPW.h>
#include <SxExt.h>

/** \brief Join wave functions of different k-points

  \note The derivation from SxGkBasis allows to set up the G-bases
        on the fly.
 */
class SX_EXPORT_EXT SxWavesJoin : public SxGkBasis, protected SxPW,
                                  public SxThis<SxWavesJoin>
{
   public:
      /** Constructor
        @param kPoints   k-points to join
        @param fileNames files containing the data
        */
      SxWavesJoin (const SxKPoints &kPoints,
                   const SxList<SxString> &fileNames);

   protected:
      /// Waves files
      SxArray<SxBinIO> files;

      /// Map of k-points (which file and which index within file)
      class Map {
         public:
            /// File id
            int iFile;
            /// k-point index in file
            int ik;
            /// Empty constructor
            Map () : iFile(-1), ik (-1) {}
      };
      /// Map of k-points (which file and which index within file)
      SxArray<Map> kMap;
      /// Number of states
      SxVector<Int> nPerK;

      /// Fermi data
      SxFermi fermi;

      /// Get number of states
      virtual int getNStates (int ik = 0) const;

      /// Get a state
      virtual const SxGBasis::TPsi operator() (int i, int iSpin, int ik) const;


   public:
      /** \brief Do the actual joining to a new file */
      void write (const SxString &filename);
};

#endif /* _SX_WAVES_JOIN_H_ */
