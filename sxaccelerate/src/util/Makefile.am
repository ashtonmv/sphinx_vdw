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

NAME      = "sxutil"
MODE      = "libOnly"
VERSION   = "1:0:0"
NAME_HASH = 484304206  # sxloghash

SX_DEPLIBS  = $(SXNONSTD)
NUM_DEPLIBS = $(REGEX)

SUBDIRS =
#if COND_BUILD_EXAMPLES
#   SUBDIRS += examples
#endif
if COND_BUILD_TESTS
   SUBDIRS += tests
endif

# ----------------------------------------------------------------------------

ACLOCAL_AMFLAGS = -I system/m4

include $(top_srcdir)/system/make.am.filelist
nodist_include_HEADERS = SxUnixConfig.h

AM_CXXFLAGS += @CXX_PEDANTIC@

# --- iOS: NSLog support in SxLog
if COND_BUILD_IOS
libsxutil_la_SOURCES += SxLog-ios.mm
endif

