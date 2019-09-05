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

NAME      = "sxparserkit_demo"
MODE      = "withLibExec"
VERSION   = "1:0:0"
NAME_HASH = 0 # sxloghash


SX_DEPLIBS  = $(SXPARSERKIT) $(SXGRAPH) $(SXFS) $(SXUTIL) $(SXNONSTD)
NUM_DEPLIBS = $(REGEX)

SUBDIRS =
SUBDIRS += 01_hello 02_stateful_scanner 03_unicode 04_include_files 05_parser_stack
SUBDIRS += 06_schema

# ----------------------------------------------------------------------------

include $(top_srcdir)/system/make.am.filelist

AM_CXXFLAGS += @CXX_PEDANTIC@

AM_CPPFLAGS += -I$(srcdir)/01_hello
AM_CPPFLAGS += -I$(srcdir)/02_stateful_scanner
AM_CPPFLAGS += -I$(srcdir)/03_unicode
AM_CPPFLAGS += -I$(srcdir)/04_include_files
AM_CPPFLAGS += -I$(srcdir)/05_parser_stack
AM_CPPFLAGS += -I$(srcdir)/06_schema
SX_DEPLIBS  += $(top_builddir)/parserkit/examples/01_hello/libsxdemo1.la
SX_DEPLIBS  += $(top_builddir)/parserkit/examples/02_stateful_scanner/libsxdemo2.la
SX_DEPLIBS  += $(top_builddir)/parserkit/examples/03_unicode/libsxdemo3.la
SX_DEPLIBS  += $(top_builddir)/parserkit/examples/04_include_files/libsxdemo4.la
SX_DEPLIBS  += $(top_builddir)/parserkit/examples/05_parser_stack/libsxdemo5.la
SX_DEPLIBS  += $(top_builddir)/parserkit/examples/06_schema/libsxdemo6.la
