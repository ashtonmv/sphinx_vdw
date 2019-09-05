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

NAME = "fs_examples"
MODE = "withBin"
NAME_HASH = 0

SX_DEPLIBS  = $(SXFS) $(SXUTIL) $(SXNONSTD)
NUM_DEPLIBS = $(REGEX) $(SHORTCUT)


# ----------------------------------------------------------------------------

include $(top_srcdir)/system/make.am.filelist
include $(top_srcdir)/system/make.am.apps


AM_CXXFLAGS += @CXX_PEDANTIC@

$(top_builddir)/fs/libsxfs.la:
	$(MAKE) -C .. libsxfs.la
