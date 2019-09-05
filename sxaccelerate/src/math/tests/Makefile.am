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

NAME = "math_tests"
MODE = "withBin"
NAME_HASH = 0

SX_DEPLIBS  = $(SXMATH) $(SXUTIL) $(SXNONSTD) $(SXMPI) $(SXIO)
NUM_DEPLIBS = $(BLAS) $(FFT)

# ----------------------------------------------------------------------------

include $(top_srcdir)/system/make.am.filelist
include $(top_srcdir)/system/make.am.apps


$(top_builddir)/math/libsxmath.la:
	$(MAKE) -C .. libsxmath.la

# enforce CPU-specific compilation, needed e.g. for SxGemmm, SxFFT2d1d etc.
AM_CXXFLAGS += @CXX_M_ARCH@
