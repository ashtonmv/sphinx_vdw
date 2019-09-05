SUBDIRS =
if COND_BUILD_APPS
   SUBDIRS += libexec
endif
if COND_BUILD_SXDOC
   SUBDIRS += tests
endif

filelist:
	for d in $(SUBDIRS); do \
      $(MAKE) -C $(top_builddir)/$$d $@; \
   done

msfilelist:
	for d in $(SUBDIRS); do \
      $(MAKE) -C $(top_builddir)/$$d $@; \
   done

