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

# Top level Automake input file. Invoke with
#    automake -a
#    autoreconf -i
#    configure ...



SUBDIRS = 
SUBDIRS += nonstd util
SUBDIRS += graph
SUBDIRS += fs
SUBDIRS += ipc
SUBDIRS += parserkit json
SUBDIRS += io
SUBDIRS += mpi
if COND_BUILD_SXMATH
   SUBDIRS += math
endif
SUBDIRS += apps system
if COND_GPU
SUBDIRS += gpu
endif
SUBDIRS += share

# build tests after all the library parts have been built to satisfy dependencies
if COND_BUILD_TESTS
   SUBDIRS += fs/tests
   SUBDIRS += ipc/tests
   SUBDIRS += math/tests
   SUBDIRS += mpi/tests
   SUBDIRS += util/tests
#   SUBDIRS += tests
endif

# build examples after all the library parts have been built to satisfy dependencies
if COND_BUILD_EXAMPLES
   SUBDIRS += fs/examples
   SUBDIRS += io/examples
   SUBDIRS += ipc/examples
   SUBDIRS += math/examples
   SUBDIRS += mpi/examples
   SUBDIRS += util/examples
endif

AM_CXXFLAGS = -I$(srcdir)/util

ACLOCAL_AMFLAGS = -I system/m4

CLEANFILES = MAKE.log

doc:
	doxygen system/Doxyfile

filelist:
	for d in $(SUBDIRS); do \
      $(MAKE) -C $(top_builddir)/$$d $@; \
   done

msfilelist:
	for d in $(SUBDIRS); do \
      $(MAKE) -C $(top_builddir)/$$d $@; \
   done

lib:
	$(MAKE) install

# create sxconfig.dat and LICENSE for "real" install (not in-place)
install-exec-local:
	-diff ${prefix}/system/$(SXBINTARGET)/sxconfig.dat system/$(SXBINTARGET)/sxconfig.dat \
   || ( mkdir -p ${prefix}/system/$(SXBINTARGET)  \
        && grep -v ac_cv_with_sxaccelerate system/$(SXBINTARGET)/sxconfig.dat > ${prefix}/system/$(SXBINTARGET)/sxconfig.dat \
        && echo "sxaccel_build=false" >>${prefix}/system/$(SXBINTARGET)/sxconfig.dat \
      )
	-diff $(srcdir)/../../LICENSE ${prefix}/LICENSE \
   || cp $(srcdir)/../../LICENSE ${prefix}/LICENSE

EXTRA_DIST=setup system/cudalt.py system/WinBison.props system/WinBison.rules system/WinBison.targets system/WinBison.xml system/WinFlex.props system/WinFlex.rules system/WinFlex.targets system/WinFlex.xml share
