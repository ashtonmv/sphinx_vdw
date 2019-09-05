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

NAME = "sxext"
MODE = "withApps"
VERSION = "1:0:0"

SX_DEPLIBS  = $(SXSTRUCT) $(SXDFT) $(SXGEOM) $(SXCLASSIC) $(SXEXX)
SX_DEPLIBS += $(SXDIRAC) $(SXIO) $(SXIPC) $(SXFS)
SX_DEPLIBS += $(SXMATH) $(SXUTIL) $(SXNONSTD)
SX_DEPLIBS += $(SXMPI)
NUM_DEPLIBS = $(FFT) $(BLAS) $(NETCDF) $(REGEX) $(SHORTCUT)

# ----------------------------------------------------------------------------
ACLOCAL_AMFLAGS = -I system/m4

include $(top_srcdir)/system/make.am.filelist
include $(top_srcdir)/system/make.am.apps

AM_CPPFLAGS += @CXX_PEDANTIC@

scripts_DATA=sxget.rc
scripts_SCRIPTS+=sxget
EXTRA_DIST=sxget.rc hamiltonian.vim sphinx.vim sxcomplete.pl
