NAME="sxplugins"
MODE="withPlugins"


SXLIBS  = $(top_builddir)/io/libsxio.la
SXLIBS += $(top_builddir)/util/libsxutil.la
SXLIBS += @LIBS@



pkglib_LTLIBRARIES = sxplugina.la sxpluginb.la

sxplugina_la_SOURCES = SxPluginA.cpp
sxplugina_la_LDFLAGS = -module -no-undefined -avoid-version
sxplugina_la_LIBADD  = $(SXLIBS)

sxpluginb_la_SOURCES = SxPluginB.cpp
sxpluginb_la_LDFLAGS = -module -no-undefined -avoid-version
sxpluginb_la_LIBADD  = $(SXLIBS)

AM_CPPFLAGS = -I../..

include $(top_srcdir)/system/make.am.apps
EXTRA_DIST=SxDemoPluginBase.h SxPluginA.h SxPluginB.h
