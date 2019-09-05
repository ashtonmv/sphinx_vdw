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

NAME      = "sxfs"
MODE      = "libOnly"
VERSION   = "1:0:0"
NAME_HASH = 2090740521  # sxloghash


SX_DEPLIBS  = $(SXNONSTD) $(SXUTIL)
NUM_DEPLIBS = $(REGEX) $(SHORTCUT)

# pulled to main src-level Makefile.am
#SUBDIRS =
#if COND_BUILD_EXAMPLES
#   SUBDIRS += examples
#endif
#if COND_BUILD_TESTS
#   SUBDIRS += tests
#endif

# ----------------------------------------------------------------------------

ACLOCAL_AMFLAGS = -I system/m4

include $(top_srcdir)/system/make.am.filelist

AM_CXXFLAGS += @CXX_PEDANTIC@
