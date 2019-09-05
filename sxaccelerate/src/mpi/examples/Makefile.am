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

NAME = "mpi_examples"
MODE = "withBin"
NAME_HASH = 0

SX_DEPLIBS  = $(SXIPC) $(SXIO) $(SXMPI) $(SXUTIL) $(SXNONSTD)
NUM_DEPLIBS = $(REGEX)

# ----------------------------------------------------------------------------

include $(top_srcdir)/system/make.am.filelist
include $(top_srcdir)/system/make.am.apps


$(top_builddir)/mpi/libsxmpi.la:
	$(MAKE) -C .. libsxmpi.la
