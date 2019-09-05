SUBDIRS = bin

SXDOC_PATH  = $(abs_top_srcdir)/share/manual/sxdoc
SXDOC_SRC  := $(wildcard $(SXDOC_PATH)/*.sxdoc)
SXDOC_MAIN := $(SXDOC_PATH)/API.sxdoc

SXTEST = $(SXQDIR)/bin/sxtest
SXTEST_OPTS  = --sxpath $(abs_top_srcdir)/share:$(SXQDIR)/share
SXTEST_OPTS += --testdir $(abs_builddir)
SXTEST_ENV   = PATH=$(PATH):./bin:${abs_top_builddir}/apps/libexec

BUILT_SOURCES = sxtest.success 

AUX_SRC = $(SXDOC_SRC)

SX_DEPLIBS = $(SXIO) $(SXFS) $(SXIPC) $(SXUTIL) $(SXNONSTD)

EXT_LIBS =

filelist:
	for d in $(SUBDIRS); do \
      $(MAKE) -C $(top_builddir)/$$d $@; \
   done

msfilelist:
	for d in $(SUBDIRS); do \
      $(MAKE) -C $(top_builddir)/$$d $@; \
   done



#AM_CXXFLAGS += @CXX_PEDANTIC@
#AM_CXXFLAGS += -I$(SXQDIR)/include

CLEANFILES = sxtest.success *.sx core*

bin/sxtestfn:
	$(MAKE) -C bin

sxtest.success: bin/sxtestfn
	$(SXTEST_ENV) $(SXTEST) $(SXTEST_OPTS) --run && touch $@

clean-local:
	@rm -rf *-*-*-* cases sxtest.success

install-exec-hook:
	mkdir -p    $(libexecdir)/tests
	cp -r cases $(libexecdir)/tests
	cp *.sx     $(libexecdir)/tests
