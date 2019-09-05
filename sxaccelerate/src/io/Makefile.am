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

NAME      = "sxio"
MODE      = "libOnly"
VERSION   = "1:0:0"
NAME_HASH = 2090740616  # sxloghash



SX_DEPLIBS  = $(SXIPC) $(SXFS) $(SXUTIL) $(SXNONSTD)
NUM_DEPLIBS = $(REGEX)

#SUBDIRS =
#if COND_BUILD_EXAMPLES
#   SUBDIRS += examples
#endif

# ----------------------------------------------------------------------------

include $(top_srcdir)/system/make.am.filelist

AM_CXXFLAGS += @CXX_PEDANTIC@

# --- LEX and YACC support
AUX_SRC       = SxParser.lex.cpp
AUX_SRC      += SxParser.tab.cpp SxParser.tab.hpp
AM_YFLAGS     = -pSxParser_ -d
BUILT_SOURCES = $(AUX_SRC)
EXTRA_DIST    = SxParser.lpp SxParser.ypp
CLEANFILES   += $(AUX_SRC) SxParser.output


%.lex.lo: AM_CXXFLAGS += @LEX_CXXFLAGS@ -Wno-unused
%.tab.lo: AM_CXXFLAGS += @LEX_CXXFLAGS@

%.lex.cpp: %.lpp %.tab.cpp
	$(LEX) -PSxParser_ $(LEXFLAGS) -o $@ $<

%.tab.cpp %.tab.hpp: %.ypp
	$(YACC) -pSxParser_ -d $< -o $@
