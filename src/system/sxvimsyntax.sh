#!/bin/sh
# author: C. Freysoldt, 2004/12/01 and later
# uses 'fold' for wrapping. Is that portable?

# syntax file (this should not be changed)
syntaxfile=${HOME}/.vim/after/syntax/cpp.vim
# the macros are searched in these headers (must not be empty)
sxbase=../sxaccelerate/src
sxmacroheaders="
  ${sxbase}/math/SxComplex.h
  ${sxbase}/util/SxAutoLoop.h
  ${sxbase}/util/SxError.h
  ${sxbase}/util/SxMemConsumer.h
  ${sxbase}/util/SxOperator.h
  ${sxbase}/util/SxTimer.h
  ${sxbase}/io/SxSimpleParser.h
  ${sxbase}/ipc/SxMPIMgr.h
  ${sxbase}/math/SxLoopMPI.h
  dirac/SxProjMatrix.h
  dirac/SxBasis.h
"

# determine home directory
SX_HOME=$1
if [ -z $SX_HOME ] ; then
  SX_HOME=.
fi

sxprecision_h=${SX_HOME}/${sxbase}/math/SxPrecision.h

# create vim syntax directory if necessary
test -e $HOME/.vim/after/syntax/ || mkdir -p $HOME/.vim/after/syntax/

# test for SX_HOME and complain if necessary
if test -r ${SX_HOME}/add-ons/SPHInX.cpp ; then :; else
cat >&2 <<USAGE_END
--------------------------------------------------------------------------------
This script creates a vim syntax file for the SPHInX source code:
$syntaxfile
SX_HOME_DIR is missing or wrong. Usage

$0 [<SX_HOME_DIR>]

The default is the current directory.
--------------------------------------------------------------------------------

USAGE_END
exit 1
fi

# Create file and write header
echo "Creating $syntaxfile" >&2

cat > $syntaxfile <<END_HEADER
" This file has been created by '$0'
" DO NOT EDIT! It will be overwritten next time.
"read possible extra-stuff (if you have some: put it there)
ru syntax/cpp_after.vim

"Skip rest of file if not SPHInX file or b:sphinx set
if strlen(matchstr(bufname ("%"), "SPHInX.cpp"))
   let b:sphinx = "SPHInX"
endif
if strlen(matchstr(bufname ("%"), "Sx.*")) == 0 && ! exists ("b:sphinx")
finish
endif
END_HEADER

if test -e $syntaxfile ; then :; else
  echo "Can't create $syntaxfile"
  exit 1
fi

# run through header files in all directories
for dir in */ ${sxbase}/*/ ${sxbase}/*/*/
do

# check for Sx*.h files
echo ${SX_HOME}/${dir}*.h | grep -q 'Sx' || continue

# get all classes class Sx...
echo "Getting SPHInX classes from ${SX_HOME}/${dir}*.h" >&2
echo '" From '${dir} >> $syntaxfile
sed -ne's/^ *class *\(Sx[^ ]*\) *.*/\1/p;s/^ *class *SX_EXPORT_[^ ]* *\(Sx[^ ]*\) *.*/\1/p' ${SX_HOME}/${dir}*.h \
| grep -v ';' \
| grep -v '<' \
| xargs | grep . | fold -s -w55 \
| sed -e's/^/syn keyword sxClass /' >> $syntaxfile

# Get typedefs Sx... type;
echo "Getting SPHInX typedefs from ${SX_HOME}/${dir}*.h" >&2
sed -ne's/^ *typedef Sx.* \([^ ]*\);.*/\1/p' ${SX_HOME}/${dir}*.h \
| xargs | grep . | fold -w55 -s \
| sed -e's/^/syn keyword sxTypedef /' >> $syntaxfile

#Get enum values
echo "Getting SPHInX enum values from ${SX_HOME}/${dir}*.h" >&2
timerCpp=`grep -l REGISTER_TIMERS ${SX_HOME}/${dir}*.cpp 2> /dev/null`
for cppFile in $timerCpp ; do
   echo " ... and timer enum values from $cppFile" >&2
done
sed -ne'
/^ *enum .*{/{
   :xx;
   /}/{
      # finalize: remove newlines, print, jump to end
      s/\n//g;
      p;
      b;
   };
   # get next line and skip comments
   N;
   s#//.*##;
   /\*\//{s#/\*.*\*/# #};
   b xx
};' ${SX_HOME}/${dir}*.h $timerCpp \
| sed -e 's/.*{//;s/}.*//;s/=[^,]*,/ /g;s/,/ /g;s/=.*$//' | xargs | fold -s -w 55 \
| sed -ne'/./s/^/syn keyword sxEnum /p' >> $syntaxfile
done

# Get precision typedefs
echo "Getting typedefs from $sxprecision_h" >&2
echo "\" From $sxprecision_h" >> $syntaxfile
sed -ne's/^ *typedef .* \([^ ]*\);.*/\1/p' $sxprecision_h \
  | xargs | fold -w55 -s \
  | sed -e's/^/syn keyword sxPrecType /' >> $syntaxfile

# preparser macros 
echo "Getting macros from ${SX_HOME}/{$sxmacroheaders}" >&2
echo 'From '$sxmacroheaders | fold -w70 -s | sed -e's/^/" /' >> $syntaxfile
headers_fullpath=`for header in ${sxmacroheaders} ; do echo "${SX_HOME}/$header" ; done`
grep '^ *# *define' $headers_fullpath \
  | grep -v 'SX_.*_H_' | sed -e's/^.*define *//;s/[ (].*$//;' | sort | uniq \
  | xargs | fold -w55 -s \
  | sed -e's/^/syn keyword sxMacro /' >> $syntaxfile

# things that aren't detected by automatic commands:
# - vector class names hidden in macro
# - VALIDATE_VECTOR in SxVec.h, but too many unwanted macros there
# - ssize_t isn't standard
cat <<END_STD_PART >> $syntaxfile

syn keyword sxClass SxVector SxDiracVec SxMatrix SxDiracMat
syn keyword sxClass SxDiracSymMat SxSymMatrix

syn keyword sxMacro VALIDATE_VECTOR SX_NO_MPI

syn keyword Type ssize_t
END_STD_PART
  
# highlighting
cat <<END_HIGHLIGHT >> $syntaxfile
hi link sxClass Type
hi link sxTypedef Type
hi link sxPrecType Type
hi link sxMacro Statement
hi link sxEnum SpecialChar
END_HIGHLIGHT

cat >&2 <<VIM_MSG_END
--------------------------------------------------------------------------------
  Your SPHInX source syntax file '$syntaxfile' 
  has been updated.
--------------------------------------------------------------------------------
VIM_MSG_END

