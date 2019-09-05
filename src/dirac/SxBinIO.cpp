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

#include <SxConfig.h>
#include <SxBinIO.h>
#include <stdio.h>
#include <stdlib.h>
#include <SxError.h>
#include <SxVector.h>
#include <unistd.h>

#ifdef USE_PARALLEL_NETCDF4
#include <SxLoopMPI.h>
#include <mpi.h>
#include <netcdf_par.h>
#endif

#include <netcdf.h>
#undef  NC_MAX_VARS
#define NC_MAX_VARS 10
#undef  MAX_NC_VARS
#define MAX_NC_VARS 10

/** Check if mode is one of the netcdf modes */
#define IS_NETCDF_MODE(mode)  (   mode == BINARY_READ_ONLY \
                               || mode == BINARY_WRITE_ONLY \
                               || mode == BINARY_WRITE_LARGE \
                               || mode == BINARY_WRITE_PARALLEL)
/** Check if mode is one of the netcdf write modes */
#define IS_NETCDF_WRITE(mode) (   mode == BINARY_WRITE_ONLY \
                               || mode == BINARY_WRITE_LARGE \
                               || mode == BINARY_WRITE_PARALLEL)

#ifdef CYGWIN
   extern "C" int mkstemp (char *);
#endif 

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
SxBinIO::SxBinIO ()
   : ncId(-1), fp(NULL), isOpen(false), err(0), mode(UNKNOWN)
{
   ncMode = WRITE_HEADER;
}


SxBinIO::SxBinIO (const SxString &filename_, Mode mode_)
   : ncId(-1), fp(NULL), isOpen(false), err(0), mode(UNKNOWN)
{
   open (filename_, mode_);
}


SxBinIO::SxBinIO (const SxBinIO &in)
   : ncId(-1), fp(NULL), isOpen(false), err(0), mode(UNKNOWN)
{
   // Using this routine may be a mistake, in most cases due to
   // foo (SxBinIO io) instead of foo(SxBinIO &io)
   // in general, it is a bad idea to have two filehandles for the
   // same file
   cout << "Warning: file handle copy constructor SxBinIO(const SxBinIO&) called.\n";
   open (in.filename, in.mode);
}

SxBinIO::~SxBinIO ()
{
   close ();
}


SxBinIO &SxBinIO::operator= (const SxBinIO &in)
{
   if (this == &in)  return *this;

   close ();
   // Using this routine may be a mistake, in most cases due to
   // io = SxBinIO (...)  instead of io.open (...)
   // in general, it is a bad idea to have two filehandles for the
   // same file
   cout << "Warning: file handle copied via SxBinIO::operator=" << endl;
   open (in.filename, in.mode);

   return *this;
}

void SxBinIO::writeXYPlot (const SxVector<TReal8> &yVec)
{
   ssize_t n = yVec.getSize();
   SxVector<TReal8> xVec (n);
   SxVector<TReal8>::Iterator x = xVec.begin();
   for (ssize_t i=0; i<n; i++)  {
      *x = (Real8)i;
      x++;
   }
   writeXYPlot (xVec, yVec);
}


void SxBinIO::writeXYPlot (const SxVector<TReal8> &xVec,
                        const SxVector<TReal8> &yVec)
{
   SX_CHECK (xVec.getSize() == yVec.getSize(),
             xVec.getSize(),   yVec.getSize());
   if (isOpen)  {
      SxVector<TReal8>::Iterator x, y;
      for (x = xVec.begin(), y=yVec.begin(); x != xVec.end(); x++, y++)  {
         fprintf (fp, "%15.12f\t%15.12f\n", *x, *y);
      }
   }  else  {
      sxprintf ("no filename given.\n");
      SX_EXIT;
   }
}

SxMatrix<Double> SxBinIO::readXYPlot ()
{
   if (isOpen)  {
      SxList<double> values;
      double x,y;
      char c;
      while (!feof(fp)) {
         while (fscanf(fp," #%c",&c) == 1) {
            while (fgetc(fp) != '\n' && !feof(fp));
         }
         int nRead = fscanf(fp, "%lf %lf",&x,&y);
            if (nRead == 2) {
               values << x << y;
            } else break;
      }
      ssize_t dim = values.getSize();
      if (dim == 0 ) {
         sxprintf ("noData read.\n");
         SX_EXIT;
      }
      if (dim % 2 == 1) {
         sxprintf ("incomplete XYPlot.\n");
         SX_EXIT;
      }

      SxMatrix<Double> result (dim/2,2,values);
      return result;

   }  else  {
      sxprintf ("no filename given.\n");
      SX_EXIT;
   }
}

void SxBinIO::writeXYPlot (const SxVector<TReal8> &xVec,
                        const SxVector<TReal8> &yVec,
                        int   xDecimals,
                        int   yDecimals)
{
   SX_CHECK (xVec.getSize() == yVec.getSize(),
             xVec.getSize(),   yVec.getSize());
   SX_CHECK (xDecimals >= 0 && xDecimals <= 30, xDecimals);
   if (yDecimals < 0) yDecimals = xDecimals;
   SX_CHECK (yDecimals >= 0 && yDecimals <= 30, yDecimals);

   if (isOpen)  {
      SxVector<TReal8>::Iterator x, y;
      for (x = xVec.begin(), y=yVec.begin(); x != xVec.end(); ++x, ++y)  {
         fprintf (fp, "%.*f\t%.*f\n", xDecimals, *x, yDecimals, *y);
      }
   }  else  {
      throw SxException(("File '" + filename + "' not open!").ascii (),
                        __FILE__, __LINE__ );
   }
}


void SxBinIO::writeXYPlot (const SxVector<Complex16> &yVec)
{
   ssize_t n = yVec.getSize();
   SxVector<TReal8> xVec (n);
   SxVector<TReal8>::Iterator x = xVec.begin();
   for (ssize_t i=0; i<n; i++)  {
      *x = (Real8)i;
      x++;
   }
   writeXYPlot (xVec, yVec);
}




void SxBinIO::writeXYPlot (const SxVector<TReal8>    &xVec,
                        const SxVector<Complex16> &yVec)
{
   SX_CHECK (xVec.getSize() == yVec.getSize(),
             xVec.getSize(),   yVec.getSize());
   if (isOpen)  {
      fprintf (fp, "# x | y.re | y.im | y^2\n");
      SxVector<TReal8>::Iterator    x;
      SxVector<Complex16>::Iterator y;
      for (x = xVec.begin(), y=yVec.begin(); x != xVec.end(); x++, y++)  {
         fprintf (fp, "%15.12f\t%15.12f\t%15.12f\t%15.12f\n",
                  *x, (*y).re, (*y).im, (*y).absSqr ());
      }
   }  else  {
      sxprintf ("no filename given.\n");
      SX_EXIT;
   }
}

void SxBinIO::writeNXYPlot (const SxMatrix<Double> &yMat)
{
   if (isOpen)  {
      ssize_t nRows = yMat.nRows();
      ssize_t nCols = yMat.nCols();
      for (ssize_t iRow = 0; iRow < nRows; iRow++)  {
         fprintf (fp, "%i ", (int)iRow + 1);
         for (ssize_t iCol = 0; iCol < nCols; iCol++)  {
            fprintf (fp, "%.6f ", yMat(iRow, iCol));
         }
         fprintf (fp, "\n");
      }
   } else {
      sxprintf ("no filename given\n");
      SX_EXIT;
   }

}


//----------------------------------------------------------------------------
// Surface plot utility
//----------------------------------------------------------------------------
void SxBinIO::writeXYZPlot (const SxMatrix<TPrecCoeffG> &M)
{
   if (isOpen)  {
      ssize_t iG, iGp, nG = M.nRows();

      for (iG = 0; iG < nG; iG++)  {
         for (iGp= 0; iGp< nG; iGp++)   {
            fprintf (fp, "%i  %i  %g  %g \n",
                          (int)iG, (int)iGp, M(iG,iGp).re, M(iG,iGp).im );
         }
         fprintf (fp, "\n");
      }

   } else {
      sxprintf ("no filename given\n");
      SX_EXIT;
   }
}



SxString SxBinIO::createTempName (const SxString &tmpDir)
{
   SxString res = tmpDir + "/tmpsxXXXXXX";  // template
#  ifdef HAVE_MKSTEMP
      char buffer[10240];
      memcpy (buffer, res.ascii(), res.getSize()+1);
      int fd = mkstemp (buffer);
      if (fd < 0)  {
         sxprintf ("Can not create temporary filename %s.\n", res.ascii());
         SX_EXIT;
      }
      res = buffer;
#  else
      cout << "Missing mkstemp. Needs workaround." << endl;
      SX_EXIT;  // not yet implemented
#  endif /* WIN32 */
// --- on non-BSD systems use:
// SxString res;
// res.resize (L_tmpnam);   or tempnam (uses TMPDIR evironment variable)
// tmpnam (res.ascii());

   // --- check, whether temporary filename is valid
   FILE *fpTmp = fopen (res.ascii(), "w");
   if ( !fpTmp )  {
      sxprintf ("Can not write to temporary directory %s.\n", tmpDir.ascii());
      SX_EXIT;
   }
   fclose (fpTmp);
   unlink (res.ascii());

   return res;
}


void SxBinIO::deleteFile (const SxString &file)
{
   unlink (file.ascii());
}

bool SxBinIO::fileExists (const SxString &file)
{
   // --- try to open file
   FILE *fp = fopen (file.ascii(), "r");
   if (!fp)  return false;

   fclose (fp);
   return true;
}


void SxBinIO::open (const SxString &filename_, Mode mode_)
{
   if (isOpen)  close ();

   bool verbose = false;

   filename = filename_;
   mode     = mode_;

   switch (mode)  {
      case ASCII_READ_ONLY       : fp = fopen (filename.ascii(), "r");
                                   isOpen = (fp != NULL);
                                   break;
      case ASCII_WRITE_ONLY      : fp = fopen (filename.ascii(), "w");
                                   isOpen = (fp != NULL);
                                   break;
      case SX_WRITE_ONLY         : fp = fopen (filename.ascii(), "w");
                                   isOpen = (fp != NULL);
                                   break;
      case SX_APPEND_ONLY        : fp = fopen (filename.ascii(), "a");
                                   isOpen = (fp != NULL);
                                   break;
      case BINARY_READ_ONLY  :     err = nc_open (filename.ascii(),
                                              NC_NOWRITE,
                                              &ncId);

                               if (err != NC_NOERR)  ncError (err);
                               isOpen = true;
                               ncMode = WRITE_HEADER;
                               break;

#ifndef USE_NETCDF4
      // use conventional NetCDF3 output
      case BINARY_WRITE_ONLY : err = nc_create (filename.ascii(),
                                                NC_CLOBBER,
                                                &ncId);
                               if (err != NC_NOERR)  ncError (err);
                               isOpen = true;
                               ncMode = WRITE_HEADER;
                               if (verbose)
                                  cout << "### SxBinIO::open() : BINARY_WRITE_ONLY" << endl;
                               break;
      case BINARY_WRITE_LARGE:
#  ifdef NC_64BIT_OFFSET
                               err = nc_create (filename.ascii(),
                                                NC_CLOBBER | NC_64BIT_OFFSET,
                                                &ncId);
#  else
    /* This is to support netcdf versions < 3.6 for some transition
       time. If you have problems with this part, update your netcdf
       library to 3.6 (see numlib/packages/ )
     */
                               err = nc_create (filename.ascii(),
                                                NC_CLOBBER,
                                                &ncId);
                               cout << endl   
 << "Warning: large netcdf files are not supported by netcdf library." << endl
 << "         The written file may be corrupted. Please update the"    << endl
 << "         library and recompile the S/PHI/ngX package if you need" << endl
 << "         large file (>2GB) support."                              << endl;
#  endif
                               if (err != NC_NOERR)  ncError (err);
                               isOpen = true;
                               ncMode = WRITE_HEADER;
                               if (verbose)
                                  cout << "### SxBinIO::open() : BINARY_WRITE_LARGE" << endl;
                               break;
#else
      // use NetCDF4 output
      case BINARY_WRITE_ONLY :
      case BINARY_WRITE_LARGE:
                                 err = nc_create (filename.ascii(),
                                       NC_NETCDF4,
                                       &ncId);
                                 if (err != NC_NOERR)  ncError (err);
                                 isOpen = true;
                                 ncMode = WRITE_HEADER;
                                 if (verbose)
                                    cout << "### SxBinIO::open() : NetCDF4-serial" << endl;
                                 break;
#endif

      case BINARY_WRITE_PARALLEL:
#ifdef USE_PARALLEL_NETCDF4
         /* NetCDF4/HDF5 and MPIIO */
                               err = nc_create_par (filename.ascii(),
                                        NC_NETCDF4 | NC_MPIIO,
                                        SxLoopMPI::MpiCommWorld(),
                                        SxLoopMPI::MpiInfoNull(),
                                        &ncId);
                               if (err != NC_NOERR)  ncError (err);
                               isOpen = true;
                               ncMode = WRITE_HEADER;
#else
                               cout << "Error: Sx was built without support for parallel NetCDF4 IO." << endl;
                               SX_EXIT;
#endif
                               if (verbose)
                                  cout << "### SxBinIO::open() : NetCDF4-parallel" << endl;
                               break;
                               
                               
      case SCRATCH_READ_WRITE: fp = fopen (filename.ascii(), "w+");
                               isOpen = (fp != NULL);
                               break;
      default                : cout << "Unknown mode in SxBinIO::open." << endl;
                               SX_EXIT;
                               break;
   }
}

void SxBinIO::close ()
{
   if (isOpen)  {
      if (IS_NETCDF_MODE(mode))  {
         err = nc_close (ncId);
         if (err != NC_NOERR)  ncError (err);
      }  else  {
         fclose (fp);
      }
   }
   fp     = NULL;
   ncId   = -1;
   isOpen = false;
   mode   = UNKNOWN;
}


void SxBinIO::ncError (int errNo) const
{
   if (mode == BINARY_READ_ONLY)
      throw SxException (( SxString("I/O Error: ") + nc_strerror(errNo)
                        + " occured during reading S/PHI/nX binary file "
                        + filename
                        ).ascii(),
                        __FILE__, __LINE__);
   else
      throw SxException (( SxString("I/O Error: ") + nc_strerror(errNo)
                        + " occured during writing S/PHI/nX binary file "
                        + filename
                        ).ascii(),
                        __FILE__, __LINE__);
}

void SxBinIO::ncError () const
{
   if (mode == BINARY_READ_ONLY)  {
      throw SxException (( SxString ("I/O Error: ") +  "Corrupt input file! "
                          "Use 'ncdump' to inspect the file:\n"
                        + filename).ascii(),
                        __FILE__, __LINE__);
   } else {
      SX_EXIT;
   }
}


void SxBinIO::addDimension (const SxString &dimName, int size) const
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);

   int dimId=0;

   // does dimension exist already and is dimension of correct size?
   if (nc_inq_dimid (ncId, dimName.ascii(), &dimId) == NC_NOERR)  {
      size_t n = 0;
      err = nc_inq_dimlen  (ncId, dimId, &n);
      if (err != NC_NOERR)  ncError (err);
      if ((size_t)size != n)  {
         throw SxException (( SxString("I/O Error: ")
                           + "Overwriting of dimensions "
                           + dimName + " is not allowed.").ascii(),
                           __FILE__, __LINE__);
      }

      // dimension exists already and has the correct size
      return;
   }
   err = nc_def_dim (ncId, dimName.ascii(), size, &dimId);
   if (err != NC_NOERR)  ncError (err);
}


int SxBinIO::getDimension (const SxString &dimName) const
{
   int dimId = 0;
   err = nc_inq_dimid (ncId, dimName.ascii(), &dimId);
   if (err != NC_NOERR)  ncError (err);

   size_t size = 0;
   err = nc_inq_dimlen  (ncId, dimId, &size);
   if (err != NC_NOERR)  ncError (err);

   return (int)size;
}


bool SxBinIO::contains (const SxString &varName) const
{
   int varId = 0;
   // variable?
   err = nc_inq_varid (ncId, varName.ascii(), &varId);
   if (err == NC_NOERR)  return true;
   // global attribute?
   err = nc_inq_attid (ncId, NC_GLOBAL, varName.ascii (), &varId);
   if (err != NC_NOERR)  return false;
   else                  return true;
}


bool SxBinIO::containsDim (const SxString &dimName) const
{
   int dimId = 0;
   err = nc_inq_dimid (ncId, dimName.ascii(), &dimId);
   if (err != NC_NOERR)  return false;
   else                  return true;
}


void SxBinIO::setMode (int  m) const
{
   SX_CHECK (ncMode == WRITE_HEADER, ncMode);
   ncMode = m;
   if (m == WRITE_DATA)  {
      err = nc_enddef (ncId);
      if (err != NC_NOERR)  ncError (err);
   }
}


int SxBinIO::getSize (const SxString &name) const
{
   size_t size;

   // --- retrieve dimension
   int dimId=0;
   err = nc_inq_dimid (ncId, name.ascii(), &dimId);
   if (err != NC_NOERR)  ncError (err);

   err = nc_inq_dimlen (ncId, dimId, &size);
   if (err != NC_NOERR)  ncError (err);

   return int(size);
}


void SxBinIO::read (const SxString &varName, SxString *val) const
{
   SX_CHECK     (val);
   SX_CHECK (mode == BINARY_READ_ONLY, mode);

   char buffer[10240];
   size_t len = 0;
   err = nc_inq_attlen (ncId, NC_GLOBAL, varName.ascii (), &len);
   char *buf = (len < 10239) ? buffer : new char[len+1];

   err = nc_get_att_text (ncId, NC_GLOBAL, varName.ascii(), buf);
   if (err != NC_NOERR)    ncError (err);
   buf[len] = 0; // trailing 0

   *val = SxString(buffer);
   if (buf != buffer) delete [] buf;
}






void SxBinIO::read (const SxString &varName, SxList<int> *val, int nElem) const
{
   SX_CHECK     (val);
   SX_CHECK (mode == BINARY_READ_ONLY, mode);

   int id = 0;
   size_t start[] = {0};
   size_t count[] = {(size_t)nElem};

   err = nc_inq_varid (ncId, varName.ascii(), &id);
   if (err != NC_NOERR)  ncError (err);

   SxVector<Int> buffer (nElem);

   err = nc_get_vara_int (ncId, id, start, count, buffer.elements);
   if (err != NC_NOERR)    ncError (err);
   if (!buffer.isValid())  ncError ();

   val->removeAll ();
   for (ssize_t i=0; i < nElem; i++)  *val << buffer(i);
}



void SxBinIO::write (const SxString &varName, const SxString &val) const
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);

   if (ncMode == WRITE_HEADER)  {
      err = nc_put_att_text (ncId, NC_GLOBAL, varName.ascii(),
                             val.getSize(), val.ascii());
      if (err != NC_NOERR)  ncError (err);

   }
}


void SxBinIO::write (const SxString &varName, double val) const
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);

   if (ncMode == WRITE_HEADER)  {
      err = nc_put_att_double (ncId, NC_GLOBAL, varName.ascii(),
                               NC_DOUBLE, 1, &val);
      if (err != NC_NOERR)  ncError (err);
   }
}

void SxBinIO::read (const SxString &varName, double *val) const
{
   SX_CHECK (val);
   SX_CHECK (mode == BINARY_READ_ONLY, mode);

   size_t len;

   err = nc_inq_attlen (ncId, NC_GLOBAL, varName.ascii (), &len);
   if (err != NC_NOERR)  ncError (err);
   if (len != 1)
      throw SxException (( "Type missfit for " + varName
                         + ": Wanted number, found vector (" + double(len)
                         + " elements).").ascii());

   err = nc_get_att_double (ncId, NC_GLOBAL, varName.ascii (), val);
   if (err != NC_NOERR)  ncError (err);
}

void SxBinIO::read (const SxString &varName, SxVector3<Int> *val) const
{
   SX_CHECK     (val);
   SX_CHECK (mode == BINARY_READ_ONLY, mode);

   int id = 0;
   size_t start[] = {0};
   size_t count[] = {3};

   err = nc_inq_varid (ncId, varName.ascii(), &id);
   if (err != NC_NOERR)  ncError (err);

   err = nc_get_vara_int (ncId, id, start, count, val->v);
   if (err != NC_NOERR)  ncError (err);

}


void SxBinIO::write (const SxString &varName,
                  const SxVector3<Int> &val,
                  const SxString &dimName) const
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);

   int id = 0;

   // --- retrieve dimension
   int dimId=0;
   err = nc_inq_dimid (ncId, dimName.ascii(), &dimId);
   if (err != NC_NOERR)  ncError (err);

   if (ncMode == WRITE_HEADER)  {

      int dims[] = {dimId};

      err = nc_def_var (ncId, varName.ascii(), NC_INT, 1, dims, &id);
      if (err != NC_NOERR)  ncError (err);

   }  else  {

      err = nc_inq_varid (ncId, varName.ascii(), &id);
      if (err != NC_NOERR)  ncError (err);

      size_t start[] = {0};
      size_t count[] = {3};

      err = nc_put_vara_int (ncId, id, start, count,
                            (const int *)val.v);
      if (err != NC_NOERR)  ncError (err);

   }
}



void SxBinIO::write (const SxString &varName,
                  const SxVector3<Double> &val,
                  const SxString &dimName) const
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);

   int id = 0;

   // --- retrieve dimension
   int dimId=0;
   err = nc_inq_dimid (ncId, dimName.ascii(), &dimId);
   if (err != NC_NOERR)  ncError (err);

   if (ncMode == WRITE_HEADER)  {

      int dims[] = {dimId};

      err = nc_def_var (ncId, varName.ascii(), NC_DOUBLE, 1, dims, &id);
      if (err != NC_NOERR)  ncError (err);

   }  else  {

      err = nc_inq_varid (ncId, varName.ascii(), &id);
      if (err != NC_NOERR)  ncError (err);

      size_t start[] = {0};
      size_t count[] = {3};

      err = nc_put_vara_double (ncId, id, start, count,
                               (const double *)val.v);
      if (err != NC_NOERR)  ncError (err);

   }
}


void SxBinIO::read (const SxString &varName, SxVector3<Double> *val) const
{
   SX_CHECK     (val);
   SX_CHECK (mode == BINARY_READ_ONLY, mode);

   int id = 0;
   size_t start[] = {0};
   size_t count[] = {3};

   err = nc_inq_varid (ncId, varName.ascii (), &id);
   if (err != NC_NOERR)  ncError (err);

   err = nc_get_vara_double (ncId, id, start, count, val->v);
   if (err != NC_NOERR)  ncError (err);
}


void SxBinIO::write (const SxString &varName,
                  const SxMatrix3<Double> &val,
                  const SxString &dimName) const
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);

   int id = 0;

   // --- retrieve dimension
   int dimId=0;
   err = nc_inq_dimid (ncId, dimName.ascii(), &dimId);
   if (err != NC_NOERR)  ncError (err);

   if (ncMode == WRITE_HEADER)  {

      int dims[] = {dimId, dimId};

      err = nc_def_var (ncId, varName.ascii(), NC_DOUBLE, 2, dims, &id);
      if (err != NC_NOERR)  ncError (err);

   }  else  {

      err = nc_inq_varid (ncId, varName.ascii(), &id);
      if (err != NC_NOERR)  ncError (err);

      size_t start[] = {0, 0};
      size_t count[] = {3, 3};

      err = nc_put_vara_double (ncId, id, start, count,
                               (const double *)val.m);
      if (err != NC_NOERR)  ncError (err);

   }
}


void SxBinIO::read (const SxString &varName, SxMatrix3<Double> *val) const
{
   SX_CHECK     (val);
   SX_CHECK (mode == BINARY_READ_ONLY, mode);

   int id = 0;
   size_t start[] = {0, 0};
   size_t count[] = {3, 3};

   err = nc_inq_varid (ncId, varName.ascii(), &id);
   if (err != NC_NOERR)  ncError (err);

   err = nc_get_vara_double (ncId, id, start, count, (double *)&val->m);
   if (err != NC_NOERR)  ncError (err);
}

void SxBinIO::write (const SxString &varName,
                  const SxArray<SxMatrix3<Double> > &val,
                  const SxString &dim1Name,
                  const SxString &dim2Name) const
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);

   int id = 0;

   // --- retrieve dimension
   int dim1Id=0, dim2Id = 0;
   err = nc_inq_dimid (ncId, dim1Name.ascii (), &dim1Id);
   if (err != NC_NOERR)  ncError (err);
   err = nc_inq_dimid (ncId, dim2Name.ascii (), &dim2Id);
   if (err != NC_NOERR)  ncError (err);

   if (ncMode == WRITE_HEADER)  {

      int dims[] = {dim1Id, dim2Id, dim2Id};

      err = nc_def_var (ncId, varName.ascii (), NC_DOUBLE, 3, dims, &id);
      if (err != NC_NOERR)  ncError (err);

   }  else  {

      err = nc_inq_varid (ncId, varName.ascii (), &id);
      if (err != NC_NOERR)  ncError (err);

      size_t start[] = {0, 0, 0};
      size_t count[] = {1, 3, 3};

      // --- write one SxMatrix3 at a time
      for (ssize_t i = 0; i < val.getSize (); i++)  {
         start[0] = (size_t)i;
         err = nc_put_vara_double (ncId, id, start, count,
                                   (const double *)val(i).m);
         if (err != NC_NOERR)  ncError (err);
      }

   }
}


void SxBinIO::read (const SxString &varName,
                 SxArray<SxMatrix3<Double> > *val) const
{
   SX_CHECK     (val);
   SX_CHECK (mode == BINARY_READ_ONLY, mode);

   int id = 0;
   size_t start[] = {0, 0, 0};
   size_t count[] = {1, 3, 3};

   err = nc_inq_varid (ncId, varName.ascii (), &id);
   if (err != NC_NOERR)  ncError (err);

   for (ssize_t i = 0; i < val->getSize (); i++)  {
      start[0] = (size_t)i;
      err = nc_get_vara_double (ncId, id, start, count, (double *)&(*val)(i).m);
      if (err != NC_NOERR)  ncError (err);
   }
}


// TODO   Check support for NetCDF4 carefully.
void SxBinIO::write (const SxString &varName,
                  const SxVector<Int> &val,
                  const SxString &dimName,
                  int offset, int nelem, int localOffset) const
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);
   SX_CHECK ((offset >= 0), offset);           // negative offsets are used
   SX_CHECK ((localOffset >= 0), localOffset); // during debugging

   int id = 0;

   // --- retrieve dimension
   int dimId=0;
   err = nc_inq_dimid (ncId, dimName.ascii(), &dimId);
   if (err != NC_NOERR)  ncError (err);

   if (ncMode == WRITE_HEADER)  {

      int dims[] = {dimId};

      err = nc_def_var (ncId, varName.ascii(), NC_INT, 1, dims, &id);
      if (err != NC_NOERR && !offset)  ncError (err);

   }  else  {

      err = nc_inq_varid (ncId, varName.ascii(), &id);
      if (err != NC_NOERR)  ncError (err);

      size_t start[] = {(size_t)offset};
      size_t count[] = {size_t(nelem == 0 ? val.getSize() : nelem)};

      int * rowPtr = val.elements;
      for (int i = 0; i < localOffset; i++) rowPtr++;

      err = nc_put_vara_int (ncId, id, start, count, rowPtr);
      if (err != NC_NOERR)  ncError (err);

   }
}


void SxBinIO::read (const SxString &varName,
                 SxVector<Int> *val, int nElem, int offset) const
{
   SX_CHECK     (val);
   SX_CHECK (mode == BINARY_READ_ONLY, mode);

   int id = 0;
   size_t start[] = {(size_t)offset};
   size_t count[] = {(size_t)nElem};

   err = nc_inq_varid (ncId, varName.ascii(), &id);
   if (err != NC_NOERR)  ncError (err);

   // val->resize (nElem);
   SX_CHECK (val->getSize () >= nElem, val->getSize (), nElem);

   err = nc_get_vara_int (ncId, id, start, count, val->elements);
   if (err != NC_NOERR)  ncError (err);
   if (!val->isValid())  ncError ();
}

void SxBinIO::read (const SxString &varName,
                 SxDiracVec<Int> *val, int nElem, int offset) const
{
   SX_CHECK     (val);
   SX_CHECK (mode == BINARY_READ_ONLY, mode);

   int id = 0;
   size_t start[] = {(size_t)offset};
   size_t count[] = {(size_t)nElem};

   err = nc_inq_varid (ncId, varName.ascii(), &id);
   if (err != NC_NOERR)  ncError (err);

   // val->resize (nElem);
   SX_CHECK (val->getSize () >= nElem, val->getSize (), nElem);

   err = nc_get_vara_int (ncId, id, start, count, val->elements);
   if (err != NC_NOERR)  ncError (err);
   if (!val->isValid())  ncError ();
}




// TODO   check support for NetCDF4
void SxBinIO::write (const SxString &varName,
                  const SxVector<Double> &val,
                  const SxString &dimName,
                  int offset, int nelem, int localOffset) const
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);
   SX_CHECK ((offset >= 0), offset);
   SX_CHECK ((localOffset >= 0), localOffset);

   int id = 0;

   // --- retrieve dimension
   int dimId=0;
   err = nc_inq_dimid (ncId, dimName.ascii(), &dimId);
   if (err != NC_NOERR)  ncError (err);

   if (ncMode == WRITE_HEADER)
   {
      int dims[] = {dimId};

      err = nc_def_var (ncId, varName.ascii(), NC_DOUBLE, 1, dims, &id);
      if (err != NC_NOERR && !offset)  ncError (err);
   }
   else
   {
      err = nc_inq_varid (ncId, varName.ascii(), &id);
      if (err != NC_NOERR)  ncError (err);

      size_t start[] = {(size_t)offset};
      size_t count[] = {size_t(nelem == 0 ? val.getSize() : nelem)};

      double * rowPtr = val.elements;
      for (int i = 0; i < localOffset; i++) rowPtr++;

      err = nc_put_vara_double (ncId, id, start, count, rowPtr);
      if (err != NC_NOERR)  ncError (err);
   }
}


void SxBinIO::read (const SxString &varName,
                 SxVector<Double> *val, int nElem, int offset) const
{
   SX_CHECK     (val);
   SX_CHECK (mode == BINARY_READ_ONLY, mode);

   int id = 0;
   size_t start[] = {(size_t)offset};
   size_t count[] = {(size_t)nElem};

   err = nc_inq_varid (ncId, varName.ascii(), &id);
   if (err != NC_NOERR)  ncError (err);

   // val->resize (nElem);
   SX_CHECK (val->getSize () >= nElem, val->getSize (), nElem);

   err = nc_get_vara_double (ncId, id, start, count, val->elements);
   if (err != NC_NOERR)  ncError (err);
   if (!val->isValid())  ncError ();
}


void SxBinIO::read (const SxString &varName,
                 SxDiracVec<Double> *val, int nElem, int offset) const
{
   SX_CHECK     (val);
   SX_CHECK (mode == BINARY_READ_ONLY, mode);

   int id = 0;
   size_t start[] = {(size_t)offset};
   size_t count[] = {(size_t)nElem};

   err = nc_inq_varid (ncId, varName.ascii(), &id);
   if (err != NC_NOERR)  ncError (err);

   // allow a partial read
   SX_CHECK (val->getSize () >= nElem, val->getSize (), nElem);
   // val->resize (nElem);

   err = nc_get_vara_double (ncId, id, start, count, val->elements);
   if (err != NC_NOERR)  ncError (err);
   if (!val->isValid())  ncError ();
}


// TODO   Check support for NetCDF4 parallel IO
void SxBinIO::write (const SxString &varName,
                  const SxVector<Complex16> &val,
                  const SxString &dimName, int offset) const
{


   if (IS_NETCDF_WRITE(mode)) {
      int idRe=0, idIm=0;

      // --- retrieve dimension
      int dimId=0;
      err = nc_inq_dimid (ncId, dimName.ascii(), &dimId);
      if (err != NC_NOERR)  ncError (err);

      if (ncMode == WRITE_HEADER)  {

         int dims[] = {dimId};

         err = nc_def_var (ncId, (varName+".re").ascii(),
               NC_DOUBLE, 1, dims, &idRe);
         if (err != NC_NOERR && !offset)  ncError (err);
         err = nc_def_var (ncId, (varName+".im").ascii(),
               NC_DOUBLE, 1, dims, &idIm);
         if (err != NC_NOERR && !offset)  ncError (err);

      }  else  {

         err = nc_inq_varid (ncId, (varName+".re").ascii(), &idRe);
         if (err != NC_NOERR)  ncError (err);
         err = nc_inq_varid (ncId, (varName+".im").ascii(), &idIm);
         if (err != NC_NOERR)  ncError (err);

         size_t start[] = {(size_t)offset};
         size_t count[] = {(size_t)val.getSize()};

         err = nc_put_vara_double (ncId, idRe, start, count,
               val.real().elements);
         if (err != NC_NOERR)  ncError (err);
         err = nc_put_vara_double (ncId, idIm, start, count,
                                val.imag().elements);
         if (err != NC_NOERR)  ncError (err);

      }
   }

   if (( mode == SX_WRITE_ONLY ) || ( mode == SX_APPEND_ONLY)) {
      ssize_t j;
      if (mode == SX_WRITE_ONLY)
         fprintf (fp, "format matrix;\n");
      if (mode == SX_APPEND_ONLY)
          fprintf (fp, "\n");


      fprintf (fp,"%s = [\n", varName.ascii ());
         SxString line;
         line = SxString("");
         for (j = 0; j < val.getSize (); j++) {
            line += SxString (val(j).re, SxString("(%8.7f, "));
            line += SxString (val(j).im,SxString("%8.7f)"));
            if (j < val.getSize() - 1) line += SxString(",");
            line += SxString (val(j).im,SxString("\n"));
         }
      fprintf ( fp,"%s", line.ascii ());
      fprintf (fp, "];\n");
   }

}


void SxBinIO::read (const SxString &varName,
                 SxVector<Complex16> *val, int nElem, int offset) const
{
   SX_CHECK      (val);
   SX_CHECK (val->getSize() == nElem, val->getSize(), nElem);
   SX_CHECK  (mode == BINARY_READ_ONLY, mode);

   int idRe=0, idIm=0;
   size_t start[] = {(size_t)offset};
   size_t count[] = {(size_t)nElem};

   err = nc_inq_varid (ncId, (varName+".re").ascii(), &idRe);
   if (err != NC_NOERR)  ncError (err);
   err = nc_inq_varid (ncId, (varName+".im").ascii(), &idIm);
   if (err != NC_NOERR)  ncError (err);

// val->resize (nElem);

   SxVector<Double> vecRe(nElem), vecIm(nElem);

   err = nc_get_vara_double (ncId, idRe, start, count, vecRe.elements);
   if (err != NC_NOERR)  ncError (err);
   err = nc_get_vara_double (ncId, idIm, start, count, vecIm.elements);
   if (err != NC_NOERR)  ncError (err);

   *val <<= SxVector<Complex16> (vecRe, vecIm);
   if (!val->isValid())  ncError();
}

void SxBinIO::read (const SxString &varName,
                 SxDiracVec<Complex16> *val, int nElem, int offset) const
{
   SX_CHECK      (val);
   SX_CHECK (val->getSize() == nElem, val->getSize(), nElem);
   SX_CHECK  (mode == BINARY_READ_ONLY, mode);

   int idRe=0, idIm=0;
   size_t start[] = {(size_t)offset};
   size_t count[] = {(size_t)nElem};

   err = nc_inq_varid (ncId, (varName+".re").ascii(), &idRe);
   if (err != NC_NOERR)  ncError (err);
   err = nc_inq_varid (ncId, (varName+".im").ascii(), &idIm);
   if (err != NC_NOERR)  ncError (err);

// val->resize (nElem);

   SxDiracVec<Double> vecRe(nElem), vecIm(nElem);

   err = nc_get_vara_double (ncId, idRe, start, count, vecRe.elements);
   if (err != NC_NOERR)  ncError (err);
   err = nc_get_vara_double (ncId, idIm, start, count, vecIm.elements);
   if (err != NC_NOERR)  ncError (err);

   *val <<= SxDiracVec<Complex16> (vecRe, vecIm);
   if (!val->isValid())  ncError();
}

// --- SxDiracMat
void SxBinIO::write (const SxString &varName,
                  const SxDiracMat<Double> &val,
                  const SxString &dim1Name, const SxString &dim2Name,
                  int offsetRow, int offsetCol) const
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);

   int id = 0;

   // --- retrieve dimension
   int dim1Id=0, dim2Id=0;
   err = nc_inq_dimid (ncId, dim1Name.ascii(), &dim1Id);
   if (err != NC_NOERR)  ncError (err);
   err = nc_inq_dimid (ncId, dim2Name.ascii(), &dim2Id);
   if (err != NC_NOERR)  ncError (err);

   if (ncMode == WRITE_HEADER)  {
      int dims[] = {dim1Id, dim2Id};

      err = nc_def_var (ncId, varName.ascii(), NC_DOUBLE, 2, dims, &id);
      if (err != NC_NOERR)  ncError (err);

   }  else  {

      err = nc_inq_varid (ncId, varName.ascii(), &id);
      if (err != NC_NOERR)  ncError (err);

      size_t size1, size2;
      err = nc_inq_dimlen (ncId, dim1Id, &size1);
      if (err != NC_NOERR)  ncError (err);
      SX_CHECK (size_t(offsetRow + val.nRows()) <= size1,
                offsetRow, val.nRows(), size1);
      err = nc_inq_dimlen (ncId, dim2Id, &size2);
      if (err != NC_NOERR)  ncError (err);
      SX_CHECK (size_t(offsetCol + val.nCols()) <= size2,
                offsetCol, val.nCols(), size2);

      size_t start[] = {(size_t)offsetRow, (size_t)offsetCol};
      size_t count[] = {(size_t)val.nRows(), 1};

      // write matrix columnwise, because SxMatrix is stored columnwise
      for (ssize_t col = 0; col < val.nCols (); col++)  {
         err = nc_put_vara_double (ncId, id, start, count,
               val.colRef(col).elements);
         if (err != NC_NOERR)  ncError (err);
         start[1]++; // next column
      }
   }
}

void SxBinIO::read (const SxString &varName,
                 SxDiracMat<Double> *val,
                 int nElemRow, int nElemCol,
                 int offsetRow, int offsetCol) const
{
   SX_CHECK      (val);
   SX_CHECK (val->nRows() == nElemRow, val->nRows(), nElemRow);
   SX_CHECK (val->nCols() == nElemCol, val->nCols(), nElemCol);
   SX_CHECK  (mode == BINARY_READ_ONLY, mode);


   int id = 0;
   size_t start[] = {(size_t)offsetRow, (size_t)offsetCol};
   size_t count[] = {(size_t)nElemRow, 1}; // read 1 column at a time

   err = nc_inq_varid (ncId, varName.ascii (), &id);
   if (err != NC_NOERR)  ncError (err);

   // read matrix columnwise, because SxMatrix is stored columnwise
   for (ssize_t col = 0; col < nElemCol; col++)  {
      err = nc_get_vara_double (ncId, id, start, count,
                                val->colRef(col).elements);
      if (err != NC_NOERR)  ncError (err);
      start[1]++; // next column
   }
}

void SxBinIO::write (const SxString &varName,
                  const SxDiracMat<Complex16> &val,
                  const SxString &dim1Name, const SxString &dim2Name,
                  int offsetRow, int offsetCol) const
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);

   int idRe=0, idIm=0;

   // --- retrieve dimensions
   int dimRowId=0, dimColId=0;
   err = nc_inq_dimid (ncId, dim1Name.ascii(), &dimRowId);
   if (err != NC_NOERR)  ncError (err);
   err = nc_inq_dimid (ncId, dim2Name.ascii(), &dimColId);
   if (err != NC_NOERR)  ncError (err);

   if (ncMode == WRITE_HEADER)  {
      int dims[] = {dimRowId, dimColId};

      err = nc_def_var (ncId, (varName+".re").ascii(),
            NC_DOUBLE, 2, dims, &idRe);
      if (err != NC_NOERR)  ncError (err);
      err = nc_def_var (ncId, (varName+".im").ascii(),
            NC_DOUBLE, 2, dims, &idIm);
      if (err != NC_NOERR)  ncError (err);

   }  else  {

      err = nc_inq_varid (ncId, (varName+".re").ascii(), &idRe);
      if (err != NC_NOERR)  ncError (err);
      err = nc_inq_varid (ncId, (varName+".im").ascii(), &idIm);
      if (err != NC_NOERR)  ncError (err);

      size_t sizeRow, sizeCol;
      err = nc_inq_dimlen (ncId, dimRowId, &sizeRow);
      if (err != NC_NOERR)  ncError (err);
      SX_CHECK (size_t(offsetRow + val.nRows()) <= sizeRow,
                offsetRow, val.nRows(), sizeRow);
      err = nc_inq_dimlen (ncId, dimColId, &sizeCol);
      if (err != NC_NOERR)  ncError (err);
      SX_CHECK (size_t(offsetCol + val.nCols()) <= sizeCol,
                offsetCol, val.nCols(), sizeCol);

      size_t start[] = {(size_t)offsetRow, (size_t)offsetCol};
      size_t count[] = {(size_t)val.nRows(), 1};

      // --- write matrix columnwise, because matrices are stored columnwise
      for (ssize_t col = 0; col < val.nCols(); col++)  {
         err = nc_put_vara_double (ncId, idRe, start, count,
               val.colRef(col).real().elements);
         if (err != NC_NOERR)  ncError (err);
         err = nc_put_vara_double (ncId, idIm, start, count,
               val.colRef(col).imag().elements);
         if (err != NC_NOERR)  ncError (err);
         start[1]++;  // next column
      }
   }
}

void SxBinIO::read (const SxString &varName,
                 SxDiracMat<Complex16> *val,
                 int nElemRow, int nElemCol,
                 int offsetRow, int offsetCol) const
{
   SX_CHECK      (val);
   SX_CHECK (val->nRows() == nElemRow, val->nRows(), nElemRow);
   SX_CHECK (val->nCols() == nElemCol, val->nCols(), nElemCol);

   SX_CHECK  (mode == BINARY_READ_ONLY, mode);

   int idRe=0, idIm=0;
   size_t start[] = {(size_t)offsetRow, (size_t)offsetCol};
   size_t count[] = {(size_t)nElemRow, 1};  // read 1 column at a time

   err = nc_inq_varid (ncId, (varName+".re").ascii(), &idRe);
   if (err != NC_NOERR)  ncError (err);
   err = nc_inq_varid (ncId, (varName+".im").ascii(), &idIm);
   if (err != NC_NOERR)  ncError (err);

   SxDiracVec<Double>    colRe(nElemRow), colIm(nElemRow);

   // --- read matrices columnwise, because they are stored columnwise
   for (ssize_t col = 0; col < nElemCol; col++)  {
      err = nc_get_vara_double (ncId, idRe, start, count, colRe.elements);
      if (err != NC_NOERR)  ncError (err);
      err = nc_get_vara_double (ncId, idIm, start, count, colIm.elements);
      if (err != NC_NOERR)  ncError (err);

      (*val).colRef(col) <<= SxDiracVec<Complex16> (colRe, colIm);
      start[1]++;  // next column
   }

   if (!val->isValid())  ncError();
}


// --- SxMatrix
void SxBinIO::write (const SxString &varName,
                  const SxMatrix<Int> &val,
                  const SxString &dimRowName, const SxString &dimColName,
                  int offsetRow, int offsetCol) const
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);

   int id = 0;

   // --- retrieve dimension
   int dimRowId=0, dimColId=0;
   err = nc_inq_dimid (ncId, dimRowName.ascii (), &dimRowId);
   if (err != NC_NOERR)  ncError (err);
   err = nc_inq_dimid (ncId, dimColName.ascii (), &dimColId);
   if (err != NC_NOERR)  ncError (err);

   if (ncMode == WRITE_HEADER)  {

      int dims[] = {dimRowId, dimColId};

      err = nc_def_var (ncId, varName.ascii (), NC_INT, 2, dims, &id);
      if (err != NC_NOERR)  ncError (err);

   }  else  {

      err = nc_inq_varid (ncId, varName.ascii (), &id);
      if (err != NC_NOERR)  ncError (err);

      size_t sizeRow, sizeCol;
      err = nc_inq_dimlen (ncId, dimRowId, &sizeRow);
      if (err != NC_NOERR)  ncError (err);
      SX_CHECK (size_t(offsetRow + val.nRows ()) <= sizeRow,
                offsetRow, val.nRows (), sizeRow);
      err = nc_inq_dimlen (ncId, dimColId, &sizeCol);
      if (err != NC_NOERR)  ncError (err);
      SX_CHECK (size_t(offsetCol + val.nCols ()) <= sizeCol,
                offsetCol, val.nCols (), sizeCol);

      size_t start[] = {(size_t)offsetRow, (size_t)offsetCol};
      size_t count[] = {(size_t)val.nRows(), 1}; // 1 column at a time

      // write matrix columnwise, because SxMatrix is stored columnwise
      for (ssize_t col = 0; col < val.nCols (); col++)  {
         err = nc_put_vara_int (ncId, id, start, count,
                                val.colRef(col).elements);
         if (err != NC_NOERR)  ncError (err);
         start[1]++; // next column
      }

   }
}


void SxBinIO::read (const SxString &varName,
                 SxMatrix<Int> *val,
                 int nElemRow, int nElemCol,
                 int offsetRow, int offsetCol) const
{
   SX_CHECK      (val);
   SX_CHECK (val->nRows() == nElemRow, val->nRows(), nElemRow);
   SX_CHECK (val->nCols() == nElemCol, val->nCols(), nElemCol);
   SX_CHECK  (mode == BINARY_READ_ONLY, mode);


   int id = 0;
   size_t start[] = {(size_t)offsetRow, (size_t)offsetCol};
   size_t count[] = {(size_t)nElemRow, 1}; // read 1 column at a time

   err = nc_inq_varid (ncId, varName.ascii (), &id);
   if (err != NC_NOERR)  ncError (err);

   // read matrix columnwise, because SxMatrix is stored columnwise
   for (ssize_t col = 0; col < nElemCol; col++)  {
      err = nc_get_vara_int (ncId, id, start, count, val->colRef(col).elements);
      if (err != NC_NOERR)  ncError (err);
      start[1]++; // next column
   }
}


void SxBinIO::read (const SxString &varName,
                 SxMatrix<Double> *val,
                 int nElemRow, int nElemCol,
                 int offsetRow, int offsetCol) const
{
   SX_CHECK      (val);
   SX_CHECK (val->nRows() == nElemRow, val->nRows(), nElemRow);
   SX_CHECK (val->nCols() == nElemCol, val->nCols(), nElemCol);

   if ( mode == BINARY_READ_ONLY ) {

      int id = 0;
      size_t start[] = {(size_t)offsetRow, (size_t)offsetCol};
      size_t count[] = {(size_t)nElemRow, 1}; // read 1 column at a time

      err = nc_inq_varid (ncId, varName.ascii (), &id);
      if (err != NC_NOERR)  ncError (err);

      // read matrix columnwise, because SxMatrix is stored columnwise
      for (ssize_t col = 0; col < nElemCol; col++)  {
         err = nc_get_vara_double (ncId, id, start, count,
                                   val->colRef(col).elements);
         if (err != NC_NOERR)  ncError (err);
         start[1]++; // next column
      }
   }

   if ( mode == ASCII_READ_ONLY ) {
      const int BUFLEN = 1024000;
      ssize_t i, j;
      char buffer[BUFLEN];
      char *bs = NULL;
      SxList<SxString> list;
      SxString tk;

      for (i = 0; i < nElemRow; i++) {
         bs = fgets (buffer, BUFLEN, fp);
         tk = SxString(bs);
         list = tk.tokenize (' ');
         if (list.last () == "\n") list.removeLast ();
         if (list.getSize () != nElemCol)  {
            cout << "Unexpected matrix size in " << filename << endl;
            cout << "In row " << (i+1) << " of " << nElemRow << endl;
            cout << "Expected  " << nElemCol << " colums, but found "
                 << list.getSize () << endl;
            SX_QUIT;
         }
         for (j = 0; j < nElemCol; j++)  {
            try {
               (*val)(i, j)=list(j).toDouble ();
            } catch (SxException e)  {
               e.print ();
               SX_QUIT;
            }
         }
      }
   }
}


// TODO   Check support for parallel IO
void SxBinIO::write (const SxString &varName,
                  const SxMatrix<Double> &val,
                  const SxString &dim1Name, const SxString &dim2Name,
                  int rowOffset, int rowNelem, int localOffset) const
{
   SX_CHECK ((rowOffset >= 0), rowOffset);
   SX_CHECK ((localOffset >= 0), localOffset);

   int id = 0;


   if ( IS_NETCDF_WRITE(mode) ) {
      // --- retrieve dimension
      int dim1Id=0, dim2Id=0;
      err = nc_inq_dimid (ncId, dim1Name.ascii(), &dim1Id);
      if (err != NC_NOERR)  ncError (err);
      err = nc_inq_dimid (ncId, dim2Name.ascii(), &dim2Id);
      if (err != NC_NOERR)  ncError (err);

      if (ncMode == WRITE_HEADER)  {

         int dims[] = {dim1Id, dim2Id};

         err = nc_def_var (ncId, varName.ascii(), NC_DOUBLE, 2, dims, &id);
         if (err != NC_NOERR)  ncError (err);

      }  else  {

         err = nc_inq_varid (ncId, varName.ascii(), &id);
         if (err != NC_NOERR)  ncError (err);

         size_t size1, size2;
         err = nc_inq_dimlen (ncId, dim1Id, &size1);
         if (err != NC_NOERR)  ncError (err);
         err = nc_inq_dimlen (ncId, dim2Id, &size2);
         if (err != NC_NOERR)  ncError (err);

         // rowOffset and rowNelem are default parameters initialized to zero
         size_t start[] = {(size_t)rowOffset, 0};
         size_t count[] = {size_t(rowNelem == 0 ? val.nRows() : rowNelem), 1};

         // write matrix columnwise, because SxMatrix is stored columnwise
         for (ssize_t col = 0; col < val.nCols (); col++)
         {
            double * rowPtr = val.colRef(col).elements;
            rowPtr += localOffset;

            err = nc_put_vara_double (ncId, id, start, count, rowPtr );
            if (err != NC_NOERR)  ncError (err);

            start[1]++; // next column
         }
      }
   }

   if ( mode == ASCII_WRITE_ONLY ) {
       ssize_t i, j;
       SxString line;
       for (i = 0; i < val.nRows (); i++) {
          line = "";
          for (j = 0; j < val.nCols (); j++)  {
             line += SxString ((double)val(i, j), "%13.12e");
             line += " ";
          }
          fprintf (fp, "%s\n", line.ascii ());
       }
   }

   if (( mode == SX_WRITE_ONLY ) || ( mode == SX_APPEND_ONLY)) {
      ssize_t i, j;
      if (mode == SX_WRITE_ONLY)
         fprintf (fp, "format matrix;\n%s = [\n", varName.ascii ());
      if (mode == SX_APPEND_ONLY) {
         cout << val << endl;
         cout << val.nRows () << endl;
         cout << val.nCols () << endl;
      }
      SxString line;

      for (i = 0; i < val.nCols (); i++) {
         line = SxString("[");
         for (j = 0; j < val.nCols (); j++) {
            line += SxString (val(i, j), SxString("%13.12e"));
            if (j < (val.nCols () - 1))
               line += SxString (", ");
         }
         line += SxString("]");
         if (i < (val.nCols () - 1))
            line += SxString(",");
         fprintf (fp, "%s\n", line.ascii ());
      }
      fprintf (fp, "];\n");

    }

}


void SxBinIO::writeMesh (const SxVector<Double>  &mesh,
                      const SxMatrix3<Double> &cell,
                      const SxVector3<Int>    &dim,
                      bool  isRowMajor)
{
   SxArray<SxVector<Double> > tmpArray (1);
   tmpArray(0) = mesh;
   writeMesh (tmpArray, cell, dim, isRowMajor);
}

void SxBinIO::writeMesh (const SxArray<SxVector<Double> > &meshes,
                      const SxMatrix3<Double>          &cell,
                      const SxVector3<Int>             &dim,
                      bool  isRowMajor) const
{
   int nMeshes = int(meshes.getSize());
   int nElem = int(meshes(0).getSize ());

   // --- write header, create dimensions
   addDimension ("nMeshes",  nMeshes);
   addDimension ("xyz",      3);
   addDimension ("meshSize", nElem);

   // --- write data
   write ("dim", SxVector<Int>(SxList<int> () << dim(0) << dim(1) << dim(2)),
          "xyz");
   write ("cell", cell, "xyz");

   SxVector<Double> out;

   SxString varName;

   for (int i=0; i < nMeshes; i++)  {
      varName = "mesh";
      if (nMeshes > 1)  varName += SxString("-") + i;

      if (isRowMajor)  out = meshes(i);
      else             out = meshes(i).toRowMajor (dim);

      write (varName, out, "meshSize");
   }
}

void SxBinIO::writeMesh (const SxArray<SxDiracVec<Double> > &meshes,
                      const SxMatrix3<Double>          &cell,
                      const SxVector3<Int>             &dim,
                      bool  isRowMajor) const
{
   int nMeshes = int(meshes.getSize());
   int nElem = int(meshes(0).getSize ());

   // --- write header, create dimensions
   addDimension ("nMeshes",  nMeshes);
   addDimension ("xyz",      3);
   addDimension ("meshSize", nElem);

   // --- write data
   write ("dim", SxVector<Int>(SxList<int> () << dim(0) << dim(1) << dim(2)),
          "xyz");
   write ("cell", cell, "xyz");

   SxDiracVec<Double> out;

   SxString varName;

   for (int i=0; i < nMeshes; i++)  {
      varName = "mesh";
      if (nMeshes > 1)  varName += SxString("-") + i;

      if (isRowMajor)  out = meshes(i);
      else             out = meshes(i).toRowMajor (dim);

      write (varName, out, "meshSize");
   }
}


SxArray<SxDiracVec<Double> > SxBinIO::readMesh (SxMatrix3<Double> *cellPtr,
                                             SxVector3<Int>    *dimPtr,
                                             bool isRowMajor) const
{
   int nMeshes = getDimension ("nMeshes");
   int nElem   = getDimension ("meshSize");
   SxVector3<Int> dim;  read("dim", &dim);

   SxArray<SxDiracVec<Double> >  meshes (nMeshes);
   SxString file;
   for (int i=0; i < nMeshes; i++)  {
      file = "mesh";
      if (nMeshes > 1)  file += SxString("-") + i;

      meshes(i).resize (nElem);
      read (file, &meshes(i), nElem);
      if (!meshes(i).isValid())  ncError();

      if (!isRowMajor) meshes(i) = meshes(i).toColMajor (dim);
   }

   if (cellPtr)  read ("cell", cellPtr);
   if (dimPtr)   read ("dim",  dimPtr);

   return meshes;
}


SxList<SxList<SxVector3<Double> > >
SxBinIO::loadStructureFHI98 () const
{
   const int BUFLEN = 1024000;
   int i, is, ia, ic;
   char buffer[BUFLEN];
   char *bs = NULL;
   int nSpecies, nAtoms;
   SxString tk;
   SxList<SxList<SxVector3<Double> > >  returnValue;
   SxList<SxString> list;

   for (i = 0; i < 3; i++)  bs = fgets (buffer, 1024, fp);

   bs = fgets (buffer, 1024, fp);
   tk = SxString(bs);
   nSpecies = tk.toInt ();

   returnValue.resize (nSpecies);

   for (is=0; is < nSpecies; is++) {
      bs = fgets (buffer, 1024, fp);
      tk = SxString(bs);
      nAtoms = tk.toInt ();
      returnValue(is).resize (nAtoms);

      bs = fgets (buffer, 1024, fp);

      for (ia=0; ia < nAtoms; ia++)  {
         bs = fgets(buffer, 1024, fp);
         tk = SxString(bs);
         list = tk.tokenize (' ');
         for (ic = 0; ic < 3; ic++) {
            returnValue(is)(ia)(ic)=list(ic).toDouble ();
         }
      }
   }
   return returnValue;
}


SxList<SxString> SxBinIO::loadSpeciesFHI98 () const
{
   cout << "FHI98 support should be replaced by standard formats\n"
              "Use add-on converters instead" << endl;
   SX_EXIT;
   const int BUFLEN = 1024000;
   int i, is, ia;
   char buffer[BUFLEN];
   char *bs = NULL;
   int nSpecies, nAtoms;
   SxString tk;
   SxList<SxString>  returnValue;


   for (i = 0; i < 3; i++)  bs = fgets (buffer, 1024, fp);

   bs = fgets (buffer, 1024, fp);
   tk = SxString(bs);
   nSpecies = tk.toInt ();

   for (is=0; is < nSpecies; is++) {
      bs = fgets (buffer, 1024, fp);
      tk = SxString(bs);

      nAtoms = tk.toInt ();

      bs = fgets (buffer, 1024, fp);
      SxString help (bs);
      SxList<SxString> helpList = help.tokenize('\n');
      returnValue.append (helpList(0));

      for (ia=0; ia < nAtoms; ia++)  {
         bs = fgets(buffer, 1024, fp);
      }
   }
   return returnValue;
}

SxMatrix3<Double> SxBinIO::loadCellFHI98 () const
{
   cout << "FHI98 support should be replaced by standard formats\n"
              "Use add-on converters instead" << endl;
   SX_EXIT;
   const int BUFLEN = 1024000;
   int i, ic;
   char buffer[BUFLEN];
   char *bs = NULL;
   SxString tk;
   SxMatrix3<Double>   returnValue;


   for (i = 0; i < 3; i++)  {
         bs = fgets(buffer, 1024, fp);
         tk = SxString(bs);
         SxList<SxString> list;
         list = tk.tokenize (' ');
         for (ic = 0; ic < 3; ic++) {
            try {
               returnValue(i, ic)=list(ic).toDouble ();
            } catch (SxException e)  {
               e.print ();
               SX_QUIT;
            }
         }
   }
   return returnValue;
}

void SxBinIO::addDoubleVar (const SxString &varName,
                         const SxString &dim1Name,
                         const SxString &dim2Name)
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);
   SX_CHECK (ncMode == WRITE_HEADER);

   // --- retrieve dimension
   int dim1Id=0, dim2Id = 0;
   err = nc_inq_dimid (ncId, dim1Name.ascii (), &dim1Id);
   if (err != NC_NOERR)  ncError (err);
   err = nc_inq_dimid (ncId, dim2Name.ascii (), &dim2Id);
   if (err != NC_NOERR)  ncError (err);

   int id = 0;
   int dims[] = {dim1Id, dim2Id};

   err = nc_def_var (ncId, varName.ascii (), NC_DOUBLE, 2, dims, &id);
   if (err != NC_NOERR)  ncError (err);

}

void SxBinIO::writeRow (const SxString &varName,
                     const SxVector<Double> &val,
                     int ir,
                     int offset)
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);
   SX_CHECK (ncMode == WRITE_DATA);

   int id = 0;
   err = nc_inq_varid (ncId, varName.ascii (), &id);
   if (err != NC_NOERR)  ncError (err);

   size_t start[] = {(size_t)ir, (size_t)offset};
   size_t count[] = {1, (size_t)val.getSize ()};

   // --- write column
   err = nc_put_vara_double (ncId, id, start, count, val.elements);
   if (err != NC_NOERR)  ncError (err);
}

void SxBinIO::readRow (const SxString &varName,
                    SxVector<Double> *val,
                    int ir,
                    int offset) const
{
   SX_CHECK  (mode == BINARY_READ_ONLY, mode);
   SX_CHECK (val);

   int id = 0;
   err = nc_inq_varid (ncId, varName.ascii (), &id);
   if (err != NC_NOERR)  ncError (err);

   size_t start[] = {(size_t)ir, (size_t)offset};
   size_t count[] = {1, (size_t)val->getSize ()};

   // --- write column
   err = nc_get_vara_double (ncId, id, start, count, val->elements);
   if (err != NC_NOERR)  ncError (err);
}
