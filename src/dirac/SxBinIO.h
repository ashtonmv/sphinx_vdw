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

#ifndef _SX_BIN_IO_H_
#define _SX_BIN_IO_H_

#include <stdio.h>
#include <SxString.h>
#include <SxException.h>
#include <SxDirac.h>
//#include <SxTypes.h>
#include <SxMatrix.h>
#include <SxMatrix3.h>
#ifdef _USE_H5
#  include <hdf5.h>
#endif /* _USE_H5 */
#include <SxFFT3d.h>
#include <SxDirac.h>


/**
  @ingroup Communication
  */
class SX_EXPORT_DIRAC SxBinIO
{
   public:
      SxString filename;

      int ncId;
      FILE *fp;
      bool isOpen;
      mutable int err;
      /** The mode determines which kind of file is created, as well
          as the allowed operations.
        */
      enum Mode { /// Read netcdf binary file
                  BINARY_READ_ONLY,
                  /// Write netcdf binary file (classic format)
                  BINARY_WRITE_ONLY,
                  /// Write netcdf binary file (large file format, 64bit)
                  BINARY_WRITE_LARGE,
                  /// Write parallel NetCDF4 file using MPI-IO
                  BINARY_WRITE_PARALLEL,
                  ///
                  ASCII_READ_ONLY,  ASCII_WRITE_ONLY,
                  SX_WRITE_ONLY, SX_APPEND_ONLY,
                  /// Open as a scratch file (allow reading and writing)
                  SCRATCH_READ_WRITE,
                  /// Uninitialized
                  UNKNOWN };
      enum NcMode { WRITE_HEADER=0, WRITE_DATA=1 };
      enum Mode mode;
      mutable int  ncMode;   // 0 - write header, 1 - write data


      SxBinIO ();
      SxBinIO (const SxBinIO &);
      SxBinIO (const SxString &filename, Mode mode=ASCII_READ_ONLY);
      ~SxBinIO ();
      void open (const SxString &filename, Mode mode=ASCII_READ_ONLY);
           
      void close ();

      void ncError (int) const;
      void ncError () const;

      SxBinIO &operator= (const SxBinIO &);

      void writeXYPlot (const SxVector<TReal8> &yVec);
      void writeXYPlot (const SxVector<Complex16> &yVec);
      void writeXYPlot (const SxVector<TReal8> &xVec,
                        const SxVector<TReal8> &yVec);
      void writeXYPlot (const SxDiracVec<TReal8> &xVec,
                        const SxDiracVec<TReal8> &yVec)
      {
         writeXYPlot (toVector(xVec), toVector(yVec));
      };
      void writeXYPlot (const SxDiracVec<TReal8> &xVec,
                        const SxDiracVec<Complex16> &yVec)
      {
         writeXYPlot (toVector(xVec), toVector(yVec));
      };
      void writeXYPlot (const SxVector<TReal8>    &xVec,
                        const SxVector<Complex16> &yVec);
      SxMatrix<Double> readXYPlot ();

      /**
        \brief XY-plot with specified precision
        \param xVec x-values
        \param yVec y-values
        \param xDecimals number of decimals for x values
        \param yDecimals number of decimals for y values
                         If it is not given, or negative, xDecimals is used

        \note The number of x-values and y-values must be the same.
        */
      void writeXYPlot (const SxVector<TReal8> &xVec,
                        const SxVector<TReal8> &yVec,
                        int   xDecimals,
                        int   yDecimals = -1);

      void writeNXYPlot (const SxMatrix<Double> &yMat);


      /** \author Abdullah Alsharif
        */
      void writeXYZPlot (const SxMatrix<TPrecCoeffG> &);

      /** \brief Creates a unique file name for temprary files

          For creation of temprary files you have to ensure to use a
          filename, that is unique. This function returns such a
          filename.
          \param tmpDir  location of temprary directory
        */
      static SxString createTempName (const SxString &tmpDir="/tmp");

      /** \brief Remove file from the file system

          If the provided file exists and with proper file/folder
          permissions this function deletes the file.
          If the file name refers to a symbolic link the link will be
          removed.
          \depecate Use SxFile::remove
       */
      static void deleteFile (const SxString &file);
      /** \brief remove a file

          \deprecate Use SxFileInfo::exists */
      static bool fileExists (const SxString &file);

      int  getSize (const SxString &name) const;
      void setMode (int) const;
      void addDimension (const SxString &dimName, int size) const;
      int  getDimension (const SxString &dimName) const;
      bool contains (const SxString &varName) const;
      bool containsDim (const SxString &dimName) const;

      /// \name Scalar values
      //@{
      /// Write a string
      void write (const SxString &varName, const SxString &val) const;
      /// Write a double to a netCDF file
      void write (const SxString &varName, double val) const;
      /// Read a double from a netCDF file
      void read (const SxString &varName, double *val) const;
      //@}

      void read  (const SxString &varName, SxString *val) const;
      // --- SxList<int>
      void read  (const SxString &varName, SxList<int> *val, int nElem) const;

      /// \name SxVector3
      //@{
      /// write SxVector3<Int>
      void write (const SxString &varName,
                  const SxVector3<Int> &val,
                  const SxString &dimName) const;
      /// read SxVector3<Int>
      void read  (const SxString &varName,
                  SxVector3<Int> *val) const;
      /// write SxVector3<Double>
      void write (const SxString &varName,
                  const SxVector3<Double> &val,
                  const SxString &dimName) const;
      /// read SxVector3<Double>
      void read  (const SxString &varName,
                  SxVector3<Double> *val) const;
      //@}

      /// \name SxMatrix3
      //@{
      /// Write SxMatrix3 to netCDF
      void write (const SxString &varName,
                  const SxMatrix3<Double> &val,
                  const SxString &dimName) const;
      /// Read SxMatrix3 from netCDF
      void read  (const SxString &varName,
                  SxMatrix3<Double> *val) const;
      /** @param dim1Name dimension of array, i.e. val.getSize ()
          @param dim2Name dimension of matrix, i.e. 3
          @par Example:
        \code
   SxArray<SxMatrix3<Double> > array;

   ...

   io.addDimension ("arraySize", array.getSize ());
   io.addDimension ("xyz", 3);
   io.write ("array", array, "arraySize", "xyz");
        \endcode
          @brief Writes an array of 3x3 matrices
        */
      void write (const SxString &varName,
                  const SxArray<SxMatrix3<Double> > &val,
                  const SxString &dim1Name,
                  const SxString &dim2Name) const;
      /** The array must have the correct size.
          @par Example:
        \code
   int size = io.getDimension ("arraySize");
   SxArray<SxMatrix3<Double> > array;
   array.resize (size);
   io.read ("array", &array);
        \endcode
          @brief Reads in an array of 3x3 matrices
        */
      void read  (const SxString &varName,
                  SxArray<SxMatrix3<Double> > *val) const;
      //@}
      // --- SxVector
      // --- changes for NetCDF4 parallel IO, please check ------------------
      void write (const SxString &varName,
                  const SxVector<Int> &val,
                  const SxString &dimName,
                  int offset=0, int nelem=0, int localOffset=0) const;
      // --------------------------------------------------------------------
      void write (const SxString &varName,
                  const SxDiracVec<Int> &val,
                  const SxString &dimName, int offset=0) const
      { write (varName, toVector(val), dimName, offset); }
      void read  (const SxString &varName,
                        SxVector<Int> *val,
                        int nElem, int offset=0) const;
      void read  (const SxString &varName,
                        SxDiracVec<Int> *val,
                        int nElem, int offset=0) const;
      // --- changes for NetCDF4 parallel IO, please check ------------------
      void write (const SxString &varName,
                  const SxVector<Double> &val,
                  const SxString &dimName,
                  int offset=0, int nelem=0, int localOffset=0) const;
      // --------------------------------------------------------------------
      void write (const SxString &varName,
                  const SxDiracVec<Double> &val,
                  const SxString &dimName, int offset=0) const
      { write (varName, toVector(val), dimName, offset); }
      void read  (const SxString &varName,
                        SxVector<Double> *val,
                        int nElem, int offset=0) const;
      void read  (const SxString &varName,
                        SxDiracVec<Double> *val,
                        int nElem, int offset=0) const;
      void write (const SxString &varName,
                  const SxVector<Complex16> &val,
                  const SxString &dimName, int offset=0) const;
      // --- changes for NetCDF4 parallel IO, please check ------------------
      void write (const SxString &varName,
                  const SxDiracVec<Complex16> &val,
                  const SxString &dimName, int offset=0) const
     { write (varName, toVector(val), dimName, offset); }
      // --------------------------------------------------------------------
      void read  (const SxString &varName,
                        SxVector<Complex16> *val,
                        int nElem, int offset=0) const;
      void read  (const SxString &varName,
                        SxDiracVec<Complex16> *val,
                        int nElem, int offset=0) const;
      /// \name SxMatrix
      //@{
      // --- changes for NetCDF4 parallel IO, please check ------------------
      void write (const SxString &varName,
                  const SxMatrix<Double> &val,
                  const SxString &dim1Name, const SxString &dim2Name,
                  int rowOffset=0, int rowNelem=0, int localOffset=0) const;
      // --------------------------------------------------------------------
      /**
        @param val pointer to a nElem1 x nElem2 matrix.
                   It must have the correct size.
        @param nElemRow number of rows
        @param nElemCol number of columns
        @param offsetRow
        @param offsetCol offsets in a larger matrix
        @brief Reads matrix from netCDF file
        */
      void read (const SxString &varName,
                 SxMatrix<Double> *val,
                 int nElemRow, int nElemCol,
                 int offsetRow = 0, int offsetCol = 0) const;
      /**
        @param varName    netCDF variable name
        @param val        matrix to be written. The matrix may be a submatrix
                          of a larger netCDF matrix
        @param dimRowName row    (1st index) name of the netCDF matrix
        @param dimColName column (2nd index) name of the netCDF matrix
        @param offsetRow
        @param offsetCol  offsets in a larger matrix
        @brief Writes matrix to netCDF file.
        */
      void write (const SxString &varName,
                  const SxMatrix<Int> &val,
                  const SxString &dimRowName, const SxString &dimColName,
                  int offsetRow = 0, int offsetCol = 0) const;
      /**
        @param val pointer to a nElem1 x nElem2 matrix.
                   It must have the correct size.
        @param nElemRow number of rows
        @param nElemCol number of columns
        @param offsetRow
        @param offsetCol offsets in a larger matrix
        @brief Reads matrix from netCDF file
        */
      void read (const SxString &varName,
                 SxMatrix<Int> *val,
                 int nElemRow, int nElemCol,
                 int offsetRow = 0, int offsetCol = 0) const;
      //@}

      /// \name SxDiracMat
      //@{
      /**
        @param varName    netCDF variable name
        @param val        matrix to be written. The matrix may be a submatrix
                          of a larger netCDF matrix
        @param dim1Name   row    (1st index) name of the netCDF matrix
        @param dim2Name   column (2nd index) name of the netCDF matrix
        @param offsetRow
        @param offsetCol offsets in a larger matrix
        @brief Writes matrix to netCDF file.
        */
      void write (const SxString &varName,
                  const SxDiracMat<Double> &val,
                  const SxString &dim1Name, const SxString &dim2Name,
                  int offsetRow = 0, int offsetCol = 0) const;

      /**
        @param val pointer to a nElem1 x nElem2 matrix.
                   It must have the correct size.
        @param nElemRow number of rows
        @param nElemCol number of columns
        @param offsetRow
        @param offsetCol offsets in a larger matrix
        @brief Reads matrix from netCDF file
        */
      void read (const SxString &varName,
                 SxDiracMat<Double> *val,
                 int nElemRow, int nElemCol,
                 int offsetRow = 0, int offsetCol = 0) const;

      /**
        @param varName    netCDF variable name
        @param val        matrix to be written. The matrix may be a submatrix
                          of a larger netCDF matrix
        @param dim1Name   row    (1st index) name of the netCDF matrix
        @param dim2Name   column (2nd index) name of the netCDF matrix
        @param offsetRow
        @param offsetCol offsets in a larger matrix
        @brief Writes matrix to netCDF file.
        */
      void write (const SxString &varName,
                  const SxDiracMat<Complex16> &val,
                  const SxString &dim1Name, const SxString &dim2Name,
                  int offsetRow = 0, int offsetCol = 0) const;

      /**
        @param val pointer to a nElem1 x nElem2 matrix.
                   It must have the correct size.
        @param nElemRow number of rows
        @param nElemCol number of columns
        @param offsetRow
        @param offsetCol offsets in a larger matrix
        @brief Reads matrix from netCDF file
        */
      void read (const SxString &varName,
                 SxDiracMat<Complex16> *val,
                 int nElemRow, int nElemCol,
                 int offsetRow = 0, int offsetCol = 0) const;
      //@}

      /// \name meshes
      //@{
      /// Write 3d mesh
      void writeMesh (const SxVector<Double>  &mesh,
                      const SxMatrix3<Double> &cell,
                      const SxVector3<Int>    &dim,
#                     ifdef USE_VECLIB
                         bool isRowMajor = false
#                     elif USE_ESSL
                         bool isRowMajor = true
#                     else /* FFTW */
                         bool isRowMajor = true
#                     endif
      );
      /// Write 3d mesh
      void writeMesh (const SxDiracVec<Double>  &mesh,
                      const SxMatrix3<Double> &cell,
                      const SxVector3<Int>    &dim,
#                     ifdef USE_VECLIB
                         bool isRowMajor = false
#                     elif USE_ESSL
                         bool isRowMajor = true
#                     else /* FFTW */
                         bool isRowMajor = true
#                     endif
      )
      {
         writeMesh (toVector(mesh), cell, dim, isRowMajor);
      }
      /// Write 3d meshes
      void writeMesh (const SxArray<SxVector<Double> > &meshes,
                      const SxMatrix3<Double>          &cell,
                      const SxVector3<Int>             &dim,
                      bool  isRowMajor = FFT_ROW_MAJOR) const;
      /// Write 3d meshes
      void writeMesh (const SxArray<SxDiracVec<Double> > &meshes,
                      const SxMatrix3<Double>          &cell,
                      const SxVector3<Int>             &dim,
                      bool  isRowMajor = FFT_ROW_MAJOR) const;
      /// Read 3d meshes
      SxArray<SxDiracVec<Double> >  readMesh (SxMatrix3<Double> *cellPtr=NULL,
                                            SxVector3<Int>    *dimPtr =NULL,
                                            bool  isRowMajor = FFT_ROW_MAJOR)
                                            const;
      //@}


      /** \brief Add a double matrix variable
        \note This is meant for cases where the data is not available
              as a matrix, but will be written with the writeRow routine
        */
      void addDoubleVar (const SxString &varName,
                         const SxString &dim1Name,
                         const SxString &dim2Name);

      /** \brief Write a row in a double matrix variable
        @param varName name of the variable
        @param val     data to be written
        @param ir      row index
        @param offset  optional column offset

        \note This is meant for cases where the data is not available
              as a matrix.
        \note Note that netCDF handles rows more efficiently than columns
              since we use the C interface.
        */
      void writeRow (const SxString &varName,
                     const SxVector<Double> &val,
                     int ir,
                     int offset = 0);

      /** \brief Read a row from a double matrix variable
        @param varName name of the variable
        @param val     data to be read
        @param ir      row index
        @param offset  optional column offset

        \note This is meant for cases where the data is not stored
              as a matrix.
        \note Note that netCDF handles rows more efficiently than columns
              since we use the C interface.
        */
      void readRow (const SxString &varName,
                    SxVector<Double> *val,
                    int ir,
                    int offset = 0) const;

      /*\brief  I/O routines for fhi98 structure format, needed for
                compatiblity purposes

        \deprecate
       */
      SxList<SxList<SxVector3<Double> > >  loadStructureFHI98 () const;

      /** \deprecate */
      SxList<SxString>  loadSpeciesFHI98 () const;

      /** \deprecate */
      SxMatrix3<Double> loadCellFHI98 () const;

};
#endif // _SX_BIN_IO_H_

