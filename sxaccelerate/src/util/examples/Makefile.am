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

NAME = "util_examples"
MODE = "withBin"
NAME_HASH = 0

SX_DEPLIBS  = $(SXUTIL) $(SXNONSTD)
NUM_DEPLIBS = $(REGEX)

# ----------------------------------------------------------------------------

include $(top_srcdir)/system/make.am.filelist
include $(top_srcdir)/system/make.am.apps

AM_CXXFLAGS += @CXX_PEDANTIC@

$(top_builddir)/util/libsxutil.la:
	$(MAKE) -C .. libsxutil.la
