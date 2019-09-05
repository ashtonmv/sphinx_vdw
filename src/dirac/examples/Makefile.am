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

NAME = "dirac_examples"
MODE = "withBin"

SX_DEPLIBS  = $(SXDIRAC) $(SXIO) $(SXIPC) $(SXMATH) $(SXFS) 
SX_DEPLIBS += $(SXUTIL) $(SXNONSTD)
NUM_DEPLIBS = $(NETCDF)

# ----------------------------------------------------------------------------

include $(top_srcdir)/system/make.am.filelist
include $(top_srcdir)/system/make.am.apps


$(top_builddir)/dirac/libsxdirac.la:
	$(MAKE) -C .. libsxdirac.la
