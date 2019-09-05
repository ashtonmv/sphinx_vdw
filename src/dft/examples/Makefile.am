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

NAME = "dft_examples"
MODE = "withBin"

SX_DEPLIBS  = $(SXDFT) $(SXGEOM) $(SXDIRAC) $(SXIO) $(SXFS) $(SXIPC)
SX_DEPLIBS += $(SXMATH) $(SXUTIL) $(SXNONSTD)
NUM_DEPLIBS = $(FFT) $(BLAS)

include $(top_srcdir)/system/make.am.filelist
include $(top_srcdir)/system/make.am.apps

$(top_builddir)/dft/libsxdft.la:
	$(MAKE) -C .. libsxdft.la
EXTRA_DIST=std
