NAME      = "sxtestbed"
MODE      = "withLibExec"
VERSION   = "1:0:0"
NAME_HASH = 0

libexecdir  = $(prefix)/libexec/tests/bin

SXDOC_PATH  = $(abs_top_srcdir)/share/manual/sxdoc
SXDOC_SRC  := $(wildcard $(SXDOC_PATH)/*.sxdoc)
SXDOC_MAIN := $(SXDOC_PATH)/Manual.sxdoc

SXGENTEST       = $(SXQDIR)/bin/sxgentests
SXGENTEST_OPTS  = --sxdoc $(SXDOC_MAIN) 
SXGENTEST_OPTS += --path $(SXDOC_PATH)
SXGENTEST_OPTS += --builddir .
SXGENTEST_OPTS += --destdir ..

BUILT_SOURCES = sxgentests-inc.h sxgentests-map.h

AUX_SRC = $(SXDOC_SRC)

SX_DEPLIBS = $(SXIO) $(SXFS) $(SXIPC) $(SXUTIL) $(SXNONSTD)

EXT_LIBS =


# ----------------------------------------------------------------------------

include $(top_srcdir)/system/make.am.filelist


#AM_CXXFLAGS += @CXX_PEDANTIC@
AM_CXXFLAGS += -I$(SXQDIR)/include

CLEANFILES += sxgentests-inc.h sxgentests-map.h core*

sxgentests-inc.h sxgentests-map.h: $(SXDOC_SRC)
	$(SXGENTEST) $(SXGENTEST_OPTS)
