# ---------------------------------------------------------------------------
#
#      The ab-initio based multiscale library
#
#                  S / P H I / n X
#
#      Copyright:  Max-Planck-Institute for Iron Research
#                  40237 Duesseldorf, Germany
#
#      Contact:    https://sxlib.mpie.de
#      Authors:    see sphinx/AUTHORS
#      License:    see sphinx/LICENSE
#
# ---------------------------------------------------------------------------

# Top level Automake input file. Invoke with
#    automake -a
#    autoreconf -i
#    configure ...

SUBDIRS = 
SUBDIRS += dirac geom dft classic struct exx
SUBDIRS += add-ons
if COND_BUILD_TOOLS
   SUBDIRS += tools
endif
SUBDIRS += share

if COND_BUILD_EXAMPLES
   SUBDIRS += dft/examples
   SUBDIRS += dirac/examples
   SUBDIRS += geom/examples
   SUBDIRS += struct/examples
endif

if COND_PLAYGROUND
   SUBDIRS += playground
endif
dist-hook:
	rm -f $(distdir)/playground/*.h
	rm -f $(distdir)/playground/*.cpp
	rm -f $(distdir)/playground/Makefile.in
	echo "#empty" > $(distdir)/playground/Make-inc.am

if COND_BUILD_DOC
   SUBDIRS += doc
endif

AM_CPPFLAGS = -I$(srcdir)/util

ACLOCAL_AMFLAGS = -I system/m4 -I ../sxaccelerate/src/system/m4

CLEANFILES = MAKE.log

filelist:
	for d in $(SUBDIRS); do \
      $(MAKE) -C $(top_builddir)/$$d $@; \
   done

msfilelist:
	for d in $(SUBDIRS); do \
      $(MAKE) -C $(top_builddir)/$$d $@; \
   done

# copy LICENSE for "real" install (not top-level in-place)
install-exec-local:
	-(test -f ${prefix}/LICENSE \
     && diff $(srcdir)/../LICENSE ${prefix}/LICENSE) \
   || cp $(srcdir)/../LICENSE ${prefix}/LICENSE

EXTRA_DIST=setup system api
