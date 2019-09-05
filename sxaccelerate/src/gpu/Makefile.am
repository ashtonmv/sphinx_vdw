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

NAME      = "sxgpu"
MODE      = "libOnly"
VERSION   = "1:0:0"
NAME_HASH = 274961564  # sxloghash


SX_DEPLIBS = $(SXUTIL) $(SXMATH)

NUM_DEPLIBS = $(PTHREADS)

#SUBDIRS =
#if COND_BUILD_EXAMPLES
#   SUBDIRS += examples
#endif


# hack to make libtool work with cuda, any suggestions for improvement? ------

CCLD = @CXXLD@
LDFLAGS += $(CUDA_LDFLAGS) $(CUDA_LIBS)

.cu.o:
	$(NVCC) $(CUDA_CFLAGS) $(AM_CXXFLAGS) -o $@ -c $<

.cu.lo:
	$(top_srcdir)/system/cudalt.py $@ $(NVCC) $(CUDA_CFLAGS) $(AM_CXXFLAGS) -c $<

# ----------------------------------------------------------------------------

include $(top_srcdir)/system/make.am.filelist
EXTRA_DIST=test
