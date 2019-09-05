# ---------------------------------------------------------------------------
#
#      The ab-initio based multiscale library
#
#                  S / P H I / n X
#
#      Copyright:  Max-Planck-Institute for Iron Research
#                  40237 Duesseldorf, Germany
#
#      Contact:    https://sxlib.mpie.de
#      Authors:    see sphinx/AUTHORS
#      License:    see sphinx/LICENSE
#
# ---------------------------------------------------------------------------

NAME = "sxdft"
MODE = "libOnly"
VERSION = "1:0:0"

SX_DEPLIBS  = $(SXGEOM) $(SXDIRAC) $(SXIO) $(SXFS) $(SXIPC)
SX_DEPLIBS += $(SXMATH) $(SXUTIL) $(SXNONSTD) $(SXMPI)
NUM_DEPLIBS = $(FFT) $(BLAS)

#SUBDIRS =
#if COND_BUILD_EXAMPLES
#   SUBDIRS += examples
#endif


# ----------------------------------------------------------------------------

include $(top_srcdir)/system/make.am.filelist

AM_CPPFLAGS += @CXX_PEDANTIC@
