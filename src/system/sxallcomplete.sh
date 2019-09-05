#!/bin/sh

SX_HOME=`echo $0 | sed -e's#[^/]*$#..#'`
ADDONS=`cd $SX_HOME/add-ons && echo *.x | sed -e's/\.x//g' | tr 'A-Z' 'a-z'`

for addon in $ADDONS; do
  echo "Analysing $addon ..." >&2
  sxcomplete $addon
done
