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

#ifndef _SX_DEFECTALIGN_UTIL_H_
#define _SX_DEFECTALIGN_UTIL_H_

#include <SxDFT.h>
#include <SxAtomicStructure.h>
#include <SxTypes.h>
#include <SxCLI.h>

// --- Utility functions for sxdefectalign tools

namespace SxDefectAlignUtil {

enum FileType {
   None,
   sxb,
   VASP_LOCPOT,
   Socorro,
   QuantumEspresso
};


SxMeshR getPot SX_EXPORT_DFT 
   (SxCell &cell, SxMesh3D &mesh, SxString &fileName,
    FileType fileType, SxAtomicStructure *structure);

SxVector<Double> readLine SX_EXPORT_DFT
                          (const SxCell &potCell,
                           const SxMesh3D &potMesh,
                           const SxMeshR &potData,
                           int idir,
                           ssize_t n,
                           const SxCell &cell,
                           const SxString &file);

SxVector<Double> average SX_EXPORT_DFT (const SxVector<Double> &x, double w);

int getMappedAtom SX_EXPORT_DFT (const Coord &pos,
                                 SxAtomicStructure &structure,
                                 const int iSpecies);

FileType getFileType SX_EXPORT_DFT (SxCLI &cli);

}

#endif /* _SX_DEFECTALIGN_UTIL_H_ */
