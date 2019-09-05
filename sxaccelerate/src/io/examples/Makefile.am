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

NAME = "io_examples"
MODE = "withBin"
NAME_HASH = 0

SX_DEPLIBS  = $(SXIO) $(SXMATH) $(SXFS) $(SXIPC) $(SXUTIL) $(SXNONSTD)
NUM_DEPLIBS = 

# ----------------------------------------------------------------------------
SUBDIRS=plugins

include $(top_srcdir)/system/make.am.filelist
include $(top_srcdir)/system/make.am.apps

$(top_builddir)/io/libsxio.la:
	$(MAKE) -C .. libsxio.la
EXTRA_DIST=createTestResults.sh parser.std parser.sx
