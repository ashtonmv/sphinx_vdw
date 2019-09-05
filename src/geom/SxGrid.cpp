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

#include<SxGrid.h>
#include<SxAtomicStructure.h>

SxGrid::SxGrid (const SxAtomicStructure &str,
                const SxVector3<Int> &gridMeshIn)
  : regCell(str.cell.getRegularCell ()),
    gridMesh(gridMeshIn),
    gridCell(regCell.basis(0) / gridMesh(0),
             regCell.basis(1) / gridMesh(1),
             regCell.basis(2) / gridMesh(2))
{
   setup (str);
}


SxGrid::SxGrid (const SxAtomicStructure &str,
                int atomsPerGridCell)
  : regCell(str.cell.getRegularCell ()),
    gridMesh(SxGrid::suggestMesh(str.cell, 
                                 str.getNAtoms () / atomsPerGridCell)),
    gridCell(regCell.basis(0) / gridMesh(0),
             regCell.basis(1) / gridMesh(1),
             regCell.basis(2) / gridMesh(2))
{
   setup(str);
}

void SxGrid::setup (const SxAtomicStructure &str)
{
   int nGrid = gridMesh.product ();
   SX_CHECK (nGrid > 0);
   nAtPerGridCell.resize (nGrid);
   nAtPerGridCell.set (0);
   
   int nAtoms = str.getNAtoms ();
   SxArray<int> gridId(nAtoms);
   mappedStr.resize (nAtoms);
   
   // --- map atoms into grid cells (sort by atom)
   Coord pos;
   int iGrid;
   SxVector3<Int> gridVec;
   for (int ia = 0; ia < nAtoms; ++ia) {
      pos = regCell.getMapped (str.constRef (ia), SxCell::Positive);
      mappedStr(ia) = pos;
      gridVec = floor(gridCell.carToRel(pos));
      // --- handle boundary numerical noise
      for (int dim = 0; dim < 3; ++dim) {
         if (gridVec(dim) == -1           ) gridVec(dim) = 0;
         if (gridVec(dim) == gridMesh(dim)) gridVec(dim) = gridMesh(dim)-1;
      }
      gridId(ia) = iGrid = (int)gridMesh.getMeshIdx(gridVec,SxMesh3D::Positive);
      nAtPerGridCell(iGrid)++;
   }
   
   // ---  setup offset
   offset.resize (nGrid);
   int sum = 0;
   SxVector<Int>::Iterator offsetIt = offset.begin (),
                           nAtIt    = nAtPerGridCell.begin ();
   for (iGrid = 0; iGrid < nGrid; ++iGrid)  {
      *offsetIt++ = sum;
      sum += *nAtIt++;
   }
   
   // --- setup idx (same as gridId, but sorted by cells)
   idx.resize(nAtoms, false, -1);
   SxVector<Int> iAtPerCell;
   iAtPerCell.copy (offset);
   for (int ia = 0; ia < nAtoms; ++ia)
      idx(iAtPerCell(gridId(ia))++) = ia;
}

SxVector3<Int> 
SxGrid::suggestMesh(const SxCell &cell, int reduction)
{
   SX_CHECK (cell.volume > 0.);
   SX_CHECK (reduction  >= 0);
   SxVector3<Int> mesh(1,1,1);

   SxCell regularCell = cell.getRegularCell ();
   double l1[3], l[3];
   for (int dim = 0; dim < 3; ++dim) 
      l1[dim] = l[dim] = regularCell.basis(dim).norm ();

   // --- basic idea: increase mesh for the longest grid cell vector
   //                 until the mesh is large enough
   int longest;
   while (mesh.product () < reduction)  {
      longest = 0;
      // find longest basis vector
      for (int dim = 1; dim < 3; ++dim)
         if (l[dim] > l[longest]) longest = dim;
      mesh(longest)++;
      l[longest] = l1[longest] / mesh(longest);
   }
   
   return mesh;
   
}

