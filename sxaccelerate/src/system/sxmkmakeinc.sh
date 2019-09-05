#!/bin/sh
#
# call this script with "(g)make filelist"
# or rerun in the srcdir (store data in $srcdir/sources)

getCanonical()
{
   stage1=`echo $1      | tr "[:upper:]" "[:lower:]"`
   stage2=`echo $stage1 | sed 's/\./_/g'`
   stage3=`echo $stage2 | sed 's/\-/_/g'`
   echo $stage3
}

# --- check if arguments come from commandline
if [ $# -ne 6 ]; then
   # --- rerun (parse ./sources for calling parameters)
   srcdir=.
   if test -e sources; then
     libname=`sed -ne's/# libname=//p' sources`
     out=`sed -ne's/# out=//p' sources`
     mode=`sed -ne's/# mode=//p' sources`
     depDirs=`sed -ne's/# depDirs=//p' sources`
     numlibs=`sed -ne's/# numlibs=//p' sources`
   else
     echo "Need file 'sources' or arguments."
     exit 1
   fi
else
   srcdir=$1     # .
   libname=$2    # pxsxext
   out=$3        # Make-inc.am
   mode=$4       # libOnly, withApps or withBin
   depDirs=$5    # sx-depdirs
   numlibs=$6    # numlibs

   cd $srcdir

   # --- update $srcdir/sources
   echo "# libname=$libname" > sources
   echo "# out=$out"         >> sources
   echo "# mode=$mode"       >> sources
   echo "# depDirs=$depDirs" >> sources
   echo "# numlibs=$numlibs" >> sources
   # --- get *.cpp but exclude automatic {.yy,.tab}.cpp files
   # setting LC_ALL garantuees a uniform sorting order on all machines
   export LC_ALL=C
   ls -1 | grep -e '.\.cpp$' -e '.\.cu$' | grep -v -e '\.yy\.cpp' -e '\.tab\.cpp' | sort -f >> sources
fi

sources=`sed -e'/^#/d' sources`
headers=`ls -1 *.h *.hpp 2>/dev/null | grep -v -e .tab.hpp -e SxUnixConfig.h 2>/dev/null | sed -ne'/./s/$/ \\\\/p' | sed -e'2,$s/^/                  /;$s/ *\\\\$//'`

rm $out   2>/dev/null
rm $out.1 2>/dev/null
rm $out.2 2>/dev/null

# --- identify MOC sources
moc_sources=""
for f in $sources; do
   hdr=`echo $f | sed 's/\.cpp/.h/'`
   if [ -f $hdr ]; then
      if [ -n "`grep Q_OBJECT $hdr`" ]; then
         fMoc=`echo $f | sed 's/\.cpp$/.moc.cpp/'`
         moc_sources="$moc_sources $fMoc"
      fi
   fi
done

programs=""
bundles=""
ltlibs=""
scripts=""
for f in $sources; do

   # (1) File list for subdirectory add-ons/
   # --- define wrapper
   if [ $mode = "withApps" ]; then
      if [ -z "`grep main $f`" ]; then
         echo "$f doesn't contain main()! Skipping..."
         sources=`echo $sources | sed "s/$f//g"`
      else
         x_source=`echo $f`
         x_target=`echo $f | sed 's/\.cpp$//g' | tr "[:upper:]" "[:lower:]"`
         am_x_source=` getCanonical \`echo $x_target | sed 's/\.cpp//g'\` `
         programs="$programs${x_target} "
         # --- define library sources
         echo "${am_x_source}_SOURCES   = \$(AUX_SRC) $x_source" >> $out.1
         if [ -n "$moc_sources" ] ; then
            echo "${am_x_source}_SOURCES  += $moc_sources"          >> $out.1
         fi
         if [ -n "`grep SX_STANDALONE $f | grep 'ifn*def'`" ]; then
            bundles="$bundles${x_target}.app "
            echo "${am_x_source}_CPPFLAGS  = \$(AM_CPPFLAGS) -DSX_STANDALONE" >> $out.1
#           echo "${am_x_source}_CPPFLAGS += -DSXEXECNAME=\"\\\"${am_x_source}\\\"\"" >> $out.1
#           echo "${am_x_source}_LDADD    = lib$libname.la \$(SXLIBS) \$(SXNUMLIBS) \$(XTRALIBS)" >> $out.1
         else
            echo "$f is not an add-on class! Skipping for lib$libname and treating as simple program..."
            sources=`echo $sources | sed "s/$f//g"`
            echo "${am_x_source}_CPPFLAGS  = \$(AM_CPPFLAGS)" >> $out.1
         fi
         echo "${am_x_source}_CPPFLAGS += -DSXLIBTOOL=\"\\\"@SXLIBTOOL@\\\"\"" >> $out.1
         echo "${am_x_source}_LDADD   = lib$libname.la \$(SXIMPORTLIBS)" >> $out.1
         echo "${am_x_source}_LDADD  += $numlibs" >> $out.1
#         echo "${am_x_source}_LDADD  += lib$libname.la" >> $out.1
         echo "${am_x_source}_LDADD  += \$(SXIMPORTOBJS)" >> $out.1
         echo                                                             >> $out.1
      fi
   fi

   # (2) File list for subdirectories tests/ and examples/
   if [ $mode = "withBin" -o $mode = "withSBin" -o $mode = "withGui" -o \
        $mode = "withCGI" -o $mode = "withFCGI" -o $mode = "withLibExec" ]; then
     x_source=`echo $f`
     x_target=`echo $f | sed 's/\.cpp$//g' | tr "[:upper:]" "[:lower:]"`
     if test $mode = "withCGI"; then
        x_target="${x_target}.cgi"
     fi
     if test $mode = "withFCGI"; then
        x_target="${x_target}.fcgi"
     fi
     am_x_source=` getCanonical \`echo $x_target | sed 's/\.cpp//g'\` `
     am_x_target=` getCanonical \`echo $x_target\` `
     programs="$programs${x_target} "
     bundles="$bundles${x_target}.app "
     echo "${am_x_target}_SOURCES   = \$(AUX_SRC)"  >> $out.1
     echo "${am_x_target}_SOURCES  += $x_source$moc_sources"  >> $out.1
     echo "${am_x_source}_CPPFLAGS  = \$(AM_CPPFLAGS) -DSXLIBTOOL=\"\\\"@SXLIBTOOL@\\\"\"" >> $out.1
#     echo "${am_x_source}_CPPFLAGS += -DSXEXECNAME=\"\\\"${am_x_source}\\\"\"" >> $out.1
     echo "${am_x_source}_LDADD     =  \$(SXIMPORTLIBS)"   >> $out.1
     echo "${am_x_source}_LDADD    +=  $numlibs"   >> $out.1
     echo "${am_x_source}_LDADD    +=  \$(SXIMPORTOBJS)"   >> $out.1
     echo                                            >> $out.1
   fi

   # (3) File list for plugin folder
   if [ $mode = "withPlugins" -o $mode = "withCondPlugins" ]; then
     x_source=`echo $f`
     moc_source=""
     hdr=`echo $f | sed 's/\.cpp/.h/'`
     if [ -f $hdr ]; then
        if [ -n "`grep Q_OBJECT $hdr`" ]; then
           moc_source=`echo $f | sed 's/\.cpp$/.moc.cpp/'`
        fi
     fi

     x_target=`echo $f | sed 's/\.cpp$//g' | tr "[:upper:]" "[:lower:]"`
     am_x_source=` getCanonical \`echo $x_target | sed 's/\.cpp//g'\` `
     am_x_target=` getCanonical \`echo $x_target\` `
     ltlibs="$ltlibs${am_x_source}.la "
     echo "${am_x_target}_la_SOURCES  = \$(AUX_SRC)"    >> $out.1
     echo "${am_x_target}_la_SOURCES += $x_source $moc_source"    >> $out.1
     echo "${am_x_target}_la_LIBADD   =  \$(SXIMPORTLIBS)"   >> $out.1
     echo "${am_x_target}_la_LDFLAGS  = -module -avoid-version" >> $out.1
     echo "${am_x_target}_la_LDFLAGS += -version-info \$(VERSION)" >> $out.1
     echo "if COND_BUILD_WIN32"                           >> $out.1
     echo "   ${am_x_target}_la_LDFLAGS += -no-undefined"    >> $out.1
     echo "endif"                                         >> $out.1
     echo                                            >> $out.1
   fi
done

sources_1by1=`echo "$sources" | sed -ne'/./s/$/ \\\\/p' | sed -e'2,$s/^/   /;$s/ *\\\\$//'`

# --- specify include paths
for d in $depDirs; do
   echo $d | grep '^@' > /dev/null
   if test $? -eq 0; then
      echo "AM_CPPFLAGS += -I$d" >> $out.2
   else
      char=`echo $d | awk '{print substr($0,0,1)}'`
      if test x"$char" = x"/"; then
         echo "AM_CPPFLAGS += -I$d" >> $out.2
      else
         echo "AM_CPPFLAGS += -I\$(top_builddir)/$d -I\$(top_srcdir)/$d" >> $out.2
      fi
   fi
done
echo "AM_CPPFLAGS += -I\$(abs_builddir) -I\$(srcdir)" >> $out.2
echo                                                  >> $out.2



echo "# This file has been generated automatically." >> $out
echo "# Do not edit it. Instead invoke"              >> $out
echo "#    (g)make filelist"                         >> $out
echo "# From the source folder"                      >> $out
echo                                                 >> $out
#echo "# Command line:"                               >> $out
#echo "#    $0 $*"                                    >> $out

if [ $mode = "empty" ]; then
echo                                                 >> $out
test $out.2 && cat $out.2                            >> $out
rm -f $out.1 $out.2
exit 0
fi

echo                                                 >> $out
echo "include_HEADERS = $headers"                    >> $out
echo                                                 >> $out

if [ $mode = "withApps" ]; then
echo "if COND_BUILD_ADDONS"                          >> $out
echo "   bin_PROGRAMS     = $programs"               >> $out
echo "endif"                                         >> $out
fi
if [ $mode = "withBin" -o $mode = "withGui" ]; then
echo "bin_PROGRAMS     = $programs"                  >> $out
if [ $mode = "withGui" ]; then
echo                                                 >> $out
echo "if BUILD_DARWIN"                               >> $out
echo "   scripts_SCRIPTS += $bundles"                >> $out
echo "endif"                                         >> $out
fi
echo                                                 >> $out
cat $out.1                                           >> $out
cat $out.2                                           >> $out
rm -f $out.1 $out.2
exit 0
fi
if [ $mode = "withSBin" ]; then
echo "sbin_PROGRAMS     = $programs"                  >> $out
echo                                                  >> $out
cat $out.1                                            >> $out
cat $out.2                                            >> $out
rm -f $out.1 $out.2
exit 0
fi
if [ $mode = "withLibExec" ]; then
echo "libexec_PROGRAMS  = $programs"                  >> $out
echo                                                  >> $out
cat $out.1                                            >> $out
cat $out.2                                            >> $out
rm -f $out.1 $out.2
exit 0
fi
if [ $mode = "withCGI" ]; then
echo "cgibin_PROGRAMS      = $programs"               >> $out
echo                                                  >> $out
cat $out.1                                            >> $out
cat $out.2                                            >> $out
rm -f $out.1 $out.2
exit 0
fi
if [ $mode = "withFCGI" ]; then
echo "fcgibin_PROGRAMS     = $programs"               >> $out
echo                                                  >> $out
cat $out.1                                            >> $out
cat $out.2                                            >> $out
rm -f $out.1 $out.2
exit 0
fi
if [ $mode = "withPlugins" ]; then
echo "pkglib_LTLIBRARIES = $ltlibs"                  >> $out
fi
if [ $mode = "withPlugins" -o $mode = "withCondPlugins" ]; then
echo                                                 >> $out
cat $out.1                                           >> $out
cat $out.2                                           >> $out
rm -f $out.1 $out.2
exit 0
fi
echo "lib_LTLIBRARIES       = lib$libname.la"        >> $out
echo "lib${libname}_la_SOURCES  = \$(AUX_SRC)"       >> $out
echo "lib${libname}_la_SOURCES += ${sources_1by1}$moc_sources" >> $out
echo                                                 >> $out
#echo "lib${libname}_la_LDFLAGS  = \$(SXLIB_LDFLAGS)"  >> $out
echo "lib${libname}_la_LIBADD   = \$(SXIMPORTLIBS)"   >> $out
echo "lib${libname}_la_LIBADD  += $numlibs"           >> $out
echo "lib${libname}_la_LIBADD  += \$(SXIMPORTOBJS)"   >> $out
echo "lib${libname}_la_LDFLAGS  = -version-info \$(VERSION)" >> $out
echo "if COND_BUILD_WIN32"                           >> $out
echo "   lib${libname}_la_LDFLAGS += -no-undefined"   >> $out
echo "endif"                                         >> $out
for m in $moc_sources; do
   source=`echo $m | sed 's/\.moc\.cpp/.cpp/g'`
   am_moc_source=`echo $m | sed 's/\.moc\.cpp//g'`
   echo >> $out
   echo "${am_moc_source}_MOC = ${source}" >> $out
done
echo                                                 >> $out
cat $out.2                                           >> $out
rm -f $out.2



if [ $mode = "withApps" ]; then
cat $out.1                                           >> $out
rm -f $out.1 $out.2
fi


