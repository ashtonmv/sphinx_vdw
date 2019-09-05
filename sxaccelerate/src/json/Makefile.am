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

NAME      = "sxjson"
MODE      = "libOnly"
VERSION   = "1:0:0"
NAME_HASH = 483908010 # sxloghash



SX_DEPLIBS  = $(SXPARSERKIT) $(SXGRAPH) $(SXFS) $(SXUTIL) $(SXNONSTD)
NUM_DEPLIBS = $(REGEX) $(PTHREADS)

SUBDIRS =
if COND_BUILD_EXAMPLES
   SUBDIRS += examples
endif
#if COND_BUILD_TESTS
#   SUBDIRS += tests
#endif


# ----------------------------------------------------------------------------

include $(top_srcdir)/system/make.am.filelist

AM_CXXFLAGS += @CXX_PEDANTIC@

# --- LEX and YACC support
AUX_SRC       = SxJSONParser.lex.cpp
AUX_SRC      += SxJSONParser.tab.cpp SxJSONParser.tab.hpp SxJSONParser.grammar
BUILT_SOURCES = $(AUX_SRC)
EXTRA_DIST  = SxJSONParser.lpp SxJSONParser.ypp
CLEANFILES += $(AUX_SRC) SxJSONParser.output

%.lex.lo: AM_CXXFLAGS += @LEX_CXXFLAGS@ -Wno-unused
%.tab.lo: AM_CXXFLAGS += @LEX_CXXFLAGS@

%.lex.cpp: %.lpp %.tab.cpp %.grammar
	$(LEX) -o $@ $<

%.tab.cpp %.tab.hpp: %.ypp
	$(YACC) --report=itemset $< -o $@

%.grammar: %.tab.cpp
	@head -1 $(@:.grammar=.output) 2>/dev/null | grep '^Grammar' > $@ \
      || ( echo "Conflicts in $(@:.grammar=.ypp). Inspect $(@:.grammar=.output)."; rm $@; exit 1)
