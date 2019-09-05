#!/bin/sh
#
# Repair make environment
#
# 2008/06/15   Sixten Boeck, boeck@mpie.de

if test -f configure.ac; then
   echo "Cleaning up Makefiles..."
   find . -name "Makefile" -exec rm {} \;
   find . -name "Makefile.in" -exec rm {} \;

   echo "Creating template Make-inc.am files..."
   system/sxmkmakeinc.sh -r

   echo "Running autoreconf..."
   ./setup

   echo "Reconfigure make environment..."
   if test -f config.status; then
      ./config.status
   else
      echo "Please invoke:"
      echo "   1. configure"
      echo "   2. make filelist"
      exit 0
   fi

   echo "Repairing dependencies..."
   cwd=`pwd`
   cmd="$cwd/$0"
   dirs=`find . -type d \( ! -regex '.*/\..*' \)`
   for d in $dirs; do
      if test -f $d/Makefile.am; then
         pushd $d
         $cwd/config.status Makefile depfiles
         popd
      fi
   done
   exit

   echo "Generating Make-inc.am files..."
   make filelist

   echo "Reconfigure make environment..."
   if test -f config.status; then
      ./config.status
   else
      echo "Please invoke:"
      echo "   configure"
      exit 0
   fi

else
   echo "ERROR: Please invoke $0 from the src/ directory."
   exit 1
fi
