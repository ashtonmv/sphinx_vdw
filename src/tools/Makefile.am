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

NAME = "sxtools"
MODE = "withBin"

SX_DEPLIBS  = $(SXEXT) $(SXSTRUCT) $(SXDFT) $(SXGEOM) $(SXEXX) $(SXCLASSIC)
SX_DEPLIBS += $(SXDIRAC) $(SXIO) $(SXIPC) $(SXFS) $(SXMATH)
SX_DEPLIBS += $(SXUTIL) $(SXNONSTD)
SX_DEPLIBS += $(SXMPI)
NUM_DEPLIBS = $(FFT) $(BLAS) $(REGEX) $(NETCDFT) $(SHORTCUT)

# ----------------------------------------------------------------------------

include $(top_srcdir)/system/make.am.filelist
include $(top_srcdir)/system/make.am.apps

AM_CPPFLAGS += @CXX_PEDANTIC@

AM_CXXFLAGS += -Wno-long-long

sxposmatel_LDADD += -lnetcdf
