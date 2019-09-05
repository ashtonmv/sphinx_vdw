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

NAME = "struct_examples"
MODE = "withBin"

SX_DEPLIBS  = $(SXSTRUCT) $(SXCLASSIC) $(SXDFT) $(SXGEOM) $(SXDIRAC)
SX_DEPLIBS += $(SXIO) $(SXIPC) $(SXFS) $(SXMATH)
SX_DEPLIBS += $(SXUTIL) $(SXNONSTD)
NUM_DEPLIBS =

# ----------------------------------------------------------------------------

include $(top_srcdir)/system/make.am.filelist
include $(top_srcdir)/system/make.am.apps
