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

NAME      = "sxlibexec"
MODE      = "withLibExec"
VERSION   = "1:0:0"
NAME_HASH = 0

SX_DEPLIBS  = $(SXSTRUCT) $(SXCLASSIC) $(SXDFT) $(SXGEOM)
SX_DEPLIBS += $(SXDIRAC) $(SXIO) $(SXIPC) $(SXFS)
SX_DEPLIBS += $(SXUTIL) $(SXNONSTD)
NUM_DEPLIBS = $(REGEX)

# ----------------------------------------------------------------------------

include $(top_srcdir)/system/make.am.filelist
include $(top_srcdir)/system/make.am.apps

AM_CXXFLAGS += @CXX_PEDANTIC@
