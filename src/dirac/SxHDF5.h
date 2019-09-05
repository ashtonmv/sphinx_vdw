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
#ifndef _SX_HDF5_H_
#define _SX_HDF5_H_

#include <SxDirac.h>
#include <SxMatrix3.h>
#include <SxMatrix.h>

#ifdef USE_HDF5
#include <H5Cpp.h>

#ifndef H5_NO_NAMESPACE
   using namespace H5;
#endif

class SX_EXPORT_DIRAC SxHDF5
{
   public:

      enum Mode { /// Read hdf5 binary file
                  BINARY_READ_ONLY,
                  /// Write hdf5 binary file (classic format)
                  BINARY_CREATE,
                  /// Open as a scratch file (allow reading and writing)
                  BINARY_READ_WRITE};

      /// Constructor
      SxHDF5 () { /*empty*/ };
      /// Constructor using input file and atomic strusture
      SxHDF5 (const SxString &filename, Mode mode=BINARY_READ_ONLY)
      {
         open (filename,mode);
      };
      /// Destructor
      ~SxHDF5 () { close (); };

      void open (const SxString &filename, Mode mode=BINARY_READ_ONLY);
      void close ();

      template <class T> T getData(const SxString &name);
      template <class T> SxVector3<T> getVector3(const SxString &name);
      template <class T> SxMatrix3<T> getMatrix3(const SxString &name);
      template <class T> SxVector<T> getVector(const SxString &name, 
            hsize_t dim, 
            hsize_t offset = 0);
      void writeVector(const SxString &name,const SxVector<Double> &vector);
      void writeMatrix(const SxString &name,const SxMatrix<Double> &matrix);
      template <class T> SxMatrix<T> getMatrix(const SxString &name,
            int nRows, 
            int nCols);
      template <class T> SxArray<T> getArray(const SxString &name,
                                             int length);
      template <class T> SxDiracVec<T> getVBlock(const SxString &name, 
            const SxArray<hsize_t> &lower,
            const SxArray<hsize_t> &upper,
            const hsize_t dim);
      SxComplex16 getWaveCoeff(const SxString &name, int iSpin, 
            int iState, int iCoeff);
      SxDiracVec<Complex16> getWaveState(const SxString &name, 
            int iSpin, int iState);
      SxVector3<Double> getCoord(const SxString &coordName, int iAtom);


      /// create group mkdir name
      void createGroup(const SxString &name)
      {
         hid_t id = H5Gcreate(group_ids.last(),
               name.ascii(),
               H5P_DEFAULT,
               H5P_DEFAULT,
               H5P_DEFAULT);
         H5Gclose(id);
      };
      /// enter group: cd name
      void enterGroup(const SxString &name)
      {
         hid_t id = H5Gopen(group_ids.last(),name.ascii(),H5P_DEFAULT);
         group_ids.append(id);
      };
      /// enter group by idx
      void enterGroupByIdx(int idx);
      /// get group name by idx
      SxString getGroupNameByIdx (int idx);
      /// leave group: cd .. 
      void leaveGroup()
      {
         hid_t id = group_ids.last ();
         H5Gclose(id);
         group_ids.removeLast ();
      };
      /// leave all groups: cd /
      void leaveAllGroups()
      {
         SxList<hid_t>::Iterator it;
         for (it = group_ids.end(); it != NULL; it--)
            H5Gclose (*it);

         hid_t id = H5Gopen(file_id,"/",H5P_DEFAULT);
         group_ids.append(id);
      };

   protected:
      hid_t file_id;
      SxList<hid_t> group_ids;
};


template <class T>
T SxHDF5::getData(const SxString &name)
{
   T result;
   hid_t dataset_id = H5Dopen(group_ids.last(),name.ascii(),H5P_DEFAULT);
   hid_t type = H5Dget_type(dataset_id);
   H5Dread(dataset_id,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,&result);
   H5Tclose(type);
   H5Dclose(dataset_id);

   return result;
}

template <class T>
SxVector3<T> SxHDF5::getVector3(const SxString &name)
{
   SxVector3<T> result;
   hid_t dataset_id = H5Dopen(group_ids.last (),name.ascii(),H5P_DEFAULT);
   hid_t type = H5Dget_type(dataset_id);
   H5Dread(dataset_id,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,result.v);
   H5Tclose(type);
   H5Dclose(dataset_id);

   return result;
}

template <class T>
SxMatrix3<T> SxHDF5::getMatrix3(const SxString &name)
{
   SxMatrix3<T> result;
   hid_t dataset_id = H5Dopen(group_ids.last (),name.ascii(),H5P_DEFAULT);
   hid_t type = H5Dget_type(dataset_id);
   H5Dread(dataset_id,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,result.m);
   H5Tclose(type);
   H5Dclose(dataset_id);

   return result;
}

template<class T>
SxVector<T> SxHDF5::getVector(const SxString &name, hsize_t dim, 
      hsize_t offset)
{
   SxVector<T> result(dim);
   SxArray<hsize_t> points(dim);
   for (hsize_t i = 0; i < dim; i++) points(i) = offset + i;
   hid_t dataset = H5Dopen(group_ids.last (),name.ascii(),H5P_DEFAULT);
   hid_t type = H5Dget_type(dataset);
   hid_t space = H5Dget_space(dataset);
   hid_t memspace = H5Screate_simple(1,&dim,&dim);
   H5Sselect_elements(space, H5S_SELECT_SET, dim, points.elements);

   H5Dread(dataset,type,memspace,space,H5P_DEFAULT,result.elements);
   H5Tclose(type);
   H5Sclose(space);
   H5Sclose(memspace);
   H5Dclose(dataset);

   return result;
}

template <class T>
SxMatrix<T> SxHDF5::getMatrix(const SxString &name, int nRows, int nCols)
{
   // Internal Transpose in SxMatix 
   SxMatrix<T> result(nCols,nRows);
   hid_t dataset_id = H5Dopen(group_ids.last (),name.ascii(),H5P_DEFAULT);
   hid_t type = H5Dget_type(dataset_id);
   H5Dread(dataset_id,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,result.elements);
   H5Tclose(type);
   H5Dclose(dataset_id);

   return result.transpose ();
}

template <class T>
SxArray<T> SxHDF5::getArray(const SxString &name, int length)
{
   SxArray<T> result(length);
   hid_t dataset_id = H5Dopen(group_ids.last (),name.ascii(),H5P_DEFAULT);
   hid_t type = H5Dget_type(dataset_id);
   H5Dread(dataset_id,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,result.elements);
   H5Tclose(type);
   H5Dclose(dataset_id);

   return result;
}

template <class T>
SxDiracVec<T> SxHDF5::getVBlock(const SxString &name, 
      const SxArray<hsize_t> &lower,
      const SxArray<hsize_t> &upper,
      const hsize_t dim)
{
   SX_CHECK (lower.getSize() == upper.getSize());
   bool isVector = true;
   for (int i = 0; i < lower.getSize (); i++)  {
      hsize_t dims = upper(i) - lower(i);
      if (dims <= 0) isVector = false;
      if (dims > 1) {
         if (isVector) isVector = false;
         else  {
            cout << "Fetching block has not vector dimension!" << endl;
            cout << "Lower boundary: [";
            for (int j = 0; j < lower.getSize (); j++)  {
               if (j > 0) cout << ",";
               cout << lower(i);
            }
            cout << "]" << endl;
            cout << "Upper boundary: [";
            for (int j = 0; j < upper.getSize (); j++)  {
               if (j > 0) cout << ",";
               cout << upper(i) << ",";
            }
            cout << "]" << endl;
            SX_EXIT;
         }
      }
   }

   //Open dataset 
   hid_t dataset = H5Dopen(group_ids.last (),name.ascii(),H5P_DEFAULT);
   hid_t dataspace = H5Dget_space(dataset);

   // get dataset type
   hid_t type = H5Dget_type(dataset);

   // select hyperslab
   H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, 
                       lower.elements, NULL, 
                       upper.elements, NULL);
   // define memspace
   int memRank = 1;
   SxArray<hsize_t> memDim = SxList<hsize_t> () << dim;
   hid_t memspace = H5Screate_simple(memRank,memDim.elements,memDim.elements);

   // read data
   SxDiracVec<T> result (dim);
   H5Dread(dataset,type,memspace,dataspace,H5P_DEFAULT,result.elements);
   
   // close ids
   H5Tclose(type);
   H5Sclose(dataspace);
   H5Sclose(memspace);
   H5Dclose(dataset);

   return result;
}
#endif /* USE_HDF5 */
#endif /* _SX_HDF5_H_ */
