# ----------------------------------------------------------------------------
#
#           The general purpose cross platform C/C++ framework
#
#                       S x A c c e l e r a t e
#
#           Home:       https://www.sxlib.de
#           License:    Apache 2
#           Authors:    see src/AUTHORS
#
# ----------------------------------------------------------------------------

NAME = "math_benchmarks"
MODE = "withBin"
NAME_HASH = 0

SX_DEPLIBS  = $(SXMATH) $(SXUTIL) $(SXNONSTD)
#if COND_MPI
SX_DEPLIBS += $(SXMPI)
#endif
NUM_DEPLIBS =  $(FFT) $(BLAS)

# ----------------------------------------------------------------------------

include $(top_srcdir)/system/make.am.filelist
include $(top_srcdir)/system/make.am.apps

AM_CXXFLAGS += @CXX_PEDANTIC@

$(top_builddir)/math/blas/libsxblas.la:
	$(MAKE) -C ../blas libsxblas.la
