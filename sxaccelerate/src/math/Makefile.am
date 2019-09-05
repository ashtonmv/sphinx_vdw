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

NAME      = "sxmath"
MODE      = "libOnly"
VERSION   = "1:0:0"
NAME_HASH = 483996378  # sxloghash


SX_DEPLIBS  = $(SXNONSTD) $(SXUTIL)
# comment out conditional to get stable Make-inc.am
#if COND_MPI
SX_DEPLIBS += $(SXMPI)
#endif
SX_DEPLIBS += $(SXIO) $(SXMPI)
NUM_DEPLIBS = $(FFT) $(BLAS)

SUBDIRS= . benchmarks

# pulled to main src-level Makefile.am
#if COND_BUILD_EXAMPLES
#   SUBDIRS += examples
#endif
#if COND_BUILD_TESTS
#   SUBDIRS += tests
#endif

filelist: subdir_filelist
subdir_filelist:
	for d in $(SUBDIRS); do \
           test $$d != "." && $(MAKE) -C $(top_builddir)/math/$$d filelist; \
        done

# ----------------------------------------------------------------------------

include $(top_srcdir)/system/make.am.filelist

AM_CXXFLAGS += @CXX_PEDANTIC@
# enforce CPU-specific compilation, needed e.g. for SxGemmm, SxFFT2d1d etc.
AM_CXXFLAGS += @CXX_M_ARCH@
