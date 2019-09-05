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


#include <SxHDF5.h>

#ifdef USE_HDF5

void SxHDF5::open (const SxString &filename, Mode mode)
{
   if (mode == BINARY_READ_ONLY)  {
            file_id = H5Fopen(filename.ascii (), H5F_ACC_RDONLY, H5P_DEFAULT);
   }

   if (mode == BINARY_CREATE)  {
            file_id = H5Fcreate(filename.ascii (), 
                  H5F_ACC_TRUNC, 
                  H5P_DEFAULT, 
                  H5P_DEFAULT);
   }

   if (mode == BINARY_READ_WRITE)  {
            file_id = H5Fopen(filename.ascii (), H5F_ACC_RDWR, H5P_DEFAULT);
   }

   hid_t id = H5Gopen(file_id,"/",H5P_DEFAULT);
   group_ids.append(id);
}

void SxHDF5::close ()
{
   leaveAllGroups ();
   H5Gclose (group_ids.last ());
   H5Fclose (file_id);
}

SxComplex<double> SxHDF5::getWaveCoeff(const SxString &name, 
      int iSpin, 
      int iState, 
      int iCoeff)
{
   SxComplex<double> result;
   hsize_t nPoints = 2;
   SxArray<hsize_t> points(4 * nPoints);
   int memRank = 1;
   SxArray<hsize_t> memDim (memRank);
   memDim(0) = 2;
   
   //Specify points
   points(0) = 0;       points(4) = 1;
   points(1) = iCoeff;  points(5) = iCoeff;
   points(2) = iState;  points(6) = iState;
   points(3) = iSpin;   points(7) = iSpin;

   hid_t dataset = H5Dopen(group_ids.last (),name.ascii(),H5P_DEFAULT);
   hid_t space = H5Dget_space(dataset);
   H5Sselect_elements(space, H5S_SELECT_SET, nPoints, points.elements);
   hid_t type = H5Dget_type(dataset);
   hid_t memspace = H5Screate_simple(memRank,memDim.elements,memDim.elements);
   H5Dread(dataset,type,memspace,space,H5P_DEFAULT,&result);
   H5Tclose(type);
   H5Sclose(space);
   H5Sclose(memspace);
   H5Dclose(dataset);

   return result;
}

SxDiracVec<Complex16> SxHDF5::getWaveState(const SxString &name, 
      int iSpin, 
      int iState)
{
   //Open dataset 
   hid_t dataset = H5Dopen(group_ids.last (),name.ascii(),H5P_DEFAULT);
   hid_t dataspace = H5Dget_space(dataset);
   
   // get dataset shape and type
   int rank = H5Sget_simple_extent_ndims(dataspace);
   SxArray<hsize_t> nDims(rank);
   H5Sget_simple_extent_dims(dataspace,nDims.elements,NULL);
   hid_t type = H5Dget_type(dataset);

   // current order in nDims: [iSpin][iState][iCoeff][iComplex]
   hsize_t nCoeffs = nDims(2);

   // specify hyper rectangle
   SxArray<hsize_t> offset = SxList<hsize_t> () << iSpin << iState << 0 << 0;
   SxArray<hsize_t> count = SxList<hsize_t> () << 1 << 1 << nCoeffs << 2;

   // select hyperslab
   H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, 
                       offset.elements, NULL, 
                       count.elements, NULL);
   
   // define memspace
   int memRank = 2;
   SxArray<hsize_t> memDim = SxList<hsize_t> () << nCoeffs << 2;
   hid_t memspace = H5Screate_simple(memRank,memDim.elements,memDim.elements);

   // read data
   SxDiracVec<Complex16> result (nCoeffs);
   H5Dread(dataset,type,memspace,dataspace,H5P_DEFAULT,result.elements);
   
   // close ids
   H5Tclose(type);
   H5Sclose(dataspace);
   H5Sclose(memspace);
   H5Dclose(dataset);

   return result;
}

SxVector3<Double> SxHDF5::getCoord(const SxString &name, int iAtom)
{
   //Open dataset 
   hid_t dataset = H5Dopen(group_ids.last (),name.ascii(),H5P_DEFAULT);
   hid_t dataspace = H5Dget_space(dataset);

   // get dataset shape and type
   int rank = H5Sget_simple_extent_ndims(dataspace);
   SxArray<hsize_t> nDims(rank);
   H5Sget_simple_extent_dims(dataspace,nDims.elements,NULL);
   hid_t type = H5Dget_type(dataset);

   // current order in nDims: [iAtom][xyz]
   hsize_t nXYZ = nDims(1);

   // specify hyper rectangle
   SxArray<hsize_t> offset = SxList<hsize_t> () << iAtom << 0;
   SxArray<hsize_t> count = SxList<hsize_t> () << 1 << nXYZ;

   // select hyperslab
   H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, 
                       offset.elements, NULL, 
                       count.elements, NULL);
   
   // define memspace
   int memRank = 1;
   SxArray<hsize_t> memDim = SxList<hsize_t> () << nXYZ;
   hid_t memspace = H5Screate_simple(memRank,memDim.elements,memDim.elements);

   // read data
   SxVector3<Double> result;
   H5Dread(dataset,type,memspace,dataspace,H5P_DEFAULT,result.v);
   
   // close ids
   H5Tclose(type);
   H5Sclose(dataspace);
   H5Sclose(memspace);
   H5Dclose(dataset);

   return result;
}



void SxHDF5::enterGroupByIdx(int idx)
{
   SxString name = getGroupNameByIdx (idx);
   enterGroup(name);
}

SxString SxHDF5::getGroupNameByIdx (int idx)
{
   size_t length = -1;
   length = H5Gget_objname_by_idx(group_ids.last (),idx,NULL,length) + 1;

   SxString result;
   result.resize(length - 1);
   
   H5Gget_objname_by_idx(group_ids.last (),idx,result.elements,length);

   return result;
}

void SxHDF5::writeVector(const SxString &name, const SxVector<Double> &vector)
{
   int rank = 1;
   SxArray<hsize_t> spaceDim (rank);
   spaceDim(0) = vector.getSize ();
   hid_t space = H5Screate_simple(rank, spaceDim.elements, spaceDim.elements);
   hid_t type = H5T_NATIVE_DOUBLE;
   hid_t data = H5Dcreate(group_ids.last (), name.ascii(), type, space,
         H5P_DEFAULT,
         H5P_DEFAULT,
         H5P_DEFAULT);
   hid_t memSpace = H5S_ALL;

   H5Dwrite(data,type,memSpace,space,
         H5P_DEFAULT,
         vector.elements);

   H5Dclose(data);
   H5Sclose(space);
}

void SxHDF5::writeMatrix(const SxString &name, const SxMatrix<Double> &matrix)
{
   int rank = 2;
   SxArray<hsize_t> spaceDim (rank);
   spaceDim(0) = matrix.nRows ();
   spaceDim(1) = matrix.nCols ();
   hid_t space = H5Screate_simple(rank, spaceDim.elements, spaceDim.elements);
   hid_t type = H5T_NATIVE_DOUBLE;
   hid_t data = H5Dcreate(group_ids.last (), name.ascii(), type, space,
         H5P_DEFAULT,
         H5P_DEFAULT,
         H5P_DEFAULT);
   hid_t memSpace = H5S_ALL;

   H5Dwrite(data, type, memSpace, space,
         H5P_DEFAULT,
         matrix.elements);

   H5Dclose(data);
   H5Sclose(space);
}
#endif
