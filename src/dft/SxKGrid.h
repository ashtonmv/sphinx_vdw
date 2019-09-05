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

#ifndef _SX_K_GRID_H_
#define _SX_K_GRID_H_

#include <SxKPoints.h>
#include <SxFFT3d.h>
#include <SxDFT.h>

/** \brief Regular k-grid

    \b SxClass = S/PHI/nX k-point grid 

    This is a service class for Monkhorst-Pack grids where access
    to the unsymmetrized grid is needed. The unsymmetrized grid is
    generated from the symmetrized one.

    \author C. Freysoldt freyso@fhi-berlin.mpg.de */
class SX_EXPORT_DFT SxKGrid : public SxKPoints
{
   public:

      /// Empty constructor
      SxKGrid ();

      /** \brief Constructor from SxKpoints
        \param in a SxKPoints object or a derived class, e.g. SxGKBasis
        \param realSpaceCell the cell containing the real space symmetries,
               i.e. structure.cell
        */
      SxKGrid (const SxKPoints &in, const SxCell &realSpaceCell);

      /** \brief Constructor from symbol table
          \param str structure
          \param table the symbol table (top level)
        */
      SxKGrid (const SxAtomicStructure &str,
               const SxSymbolTable *table);
      
      /** \brief Constructor from sxb file
        \param io open sxb file in BINARY_READ_ONLY mode
        \param realSpaceCell the cell containing the real space symmetries,
               i.e. structure.cell
        */
      SxKGrid (const SxBinIO &io, const SxCell &realSpaceCell);

      /// Destructor
      ~SxKGrid () {}

      /** \brief Derive grid from symmetrized points
          \param cell The real space cell with the symmetries
        */
      void setup (const SxCell &cell);
      
      /// The grid cell
      SxCell subCell;

      /// An R<->k FFT
      SxFFT3d fftRk;

      /// Request FFT setup
      const SxFFT3d& setupFFT (const SxFFT::Directions dir);

      /** \brief The k -> symmetrized k map. Size: fftRk.getSize ()

          This stores the mapping fftIdx <-> ik in the following way:
          - symK(i) =   ik+1  if S ^ getVec(i) + k0 = + getK(ik)
          - symK(i) = -(ik+1) if S ^ getVec(i) + k0 = - getK(ik)
          \example
          \code
int meshSize = kgrid.fftRk.getSize ();
SxDiracVec<Complex16> fftMesh(meshSize), values(nk);
for (int i = 0; i < meshSize; ++i)
   if ((ik = kGrid.symK(i)) < 0) 
      fftMesh = values(ik-1);
   else
      fftMesh = values(-ik -1).conj ();
\endcode
      */
      SxVector<Int> symK;

      void setMesh (SxDiracVec<Complex16> *mesh, 
                    const SxDiracVec<Complex16> &values);

      const SxKGrid& operator= (const SxKGrid &)
      { SX_EXIT; /* to be implemented */ return *this; }
      const SxKGrid& operator= (const SxKPoints &)
      { SX_EXIT; /* to be implemented */ return *this; }

};

#endif /* _SX_K_GRID_H_ */
