AC_DEFUN([SXQ_TEST_LATEX], [
   rm -rf .tmps_latex
   mkdir .tmps_latex
   cd .tmps_latex
   $2="no"; export $2;
   ifelse($#, 2, [
      cat > testconf.tex <<EOF
      $1
EOF
      ],$#,3,[
      echo "\\documentclass{$3}" > testconf.tex
      cat >> testconf.tex <<EOF
      $1
EOF
      ],$#,4,[
      echo "\\documentclass{$3}" >  testconf.tex
      echo "\\usepackage{$4}"    >> testconf.tex
      cat >> testconf.tex <<EOF
      $1
EOF
   ])
   cat testconf.tex | $ac_cv_path_LATEX 2>&1 >/dev/null && $2=yes; export $2;
   cd ..
   rm -rf .tmps_latex
])



AC_DEFUN([SXQ_LATEX_CLASS],[
   AC_CACHE_CHECK([for LaTeX class $1],[ac_cv_latex_class_]translit($1,[-],[_]),[
      SXQ_TEST_LATEX([
         \begin{document}
         \end{document}
      ],[ac_cv_latex_class_]translit($1,[-],[_]),$1)
   ])
   hasClass=$[ac_cv_latex_class_]translit($1,[-],[_])
#  test x"$hasClass" = x"yes" || AC_MSG_ERROR([LaTeX class $1 not installed.])
])

AC_DEFUN([SXQ_LATEX_PKG],[
   if test "$[ac_cv_latex_class_]translit($2,[-],[_])" = "" ;
   then
      SXQ_LATEX_CLASS($2,classparams)
      export classparams;
   else
      classparams=$[ac_cv_latex_class_]translit($2,[-],[_]) ;
      export classparams;
   fi;
   if test $classparams = "no" ;
   then
       AC_MSG_ERROR([Unable to find $2 class])
   fi
   AC_CACHE_CHECK([for $1 in class $2],[ac_cv_latex_]translit($1,[-],[_])[_]translit($2,[-],[_]),[
   SXQ_TEST_LATEX([
      \documentclass{$2}
      \usepackage{$1}
      \begin{document}
      \end{document}
   ],[ac_cv_latex_]translit($1,[-],[_])[_]translit($2,[-],[_]))
   ])
   hasPkg=$[ac_cv_latex_]translit($1,[-],[_])[_]translit($2,[-],[_])
   test x"$hasPkg" = x"yes" || AC_MSG_ERROR([LaTeX package $1 not found in class $2.])
])



AC_DEFUN([SXQ_CHECK_DOC], [

   SX_ARG_ENABLE( [sxqdoc], [USE_SXQDOC], [no],
                  [build SxQDOC documentation])

   useSxQDoc="no"
   SXQDIR=""

   if test x"$ac_cv_enable_sxqdoc" = x"yes"; then

      if test -n "$SXQ"; then
         SXQPATH="$SXQ/bin:$SXQ/libexec:$SXQ/apps:$PATH"
      else
         SXQPATH="${PATH}"
         if test -f /opt/sxq/bin/sxpublisher; then
            SXQPATH="$SXQPATH:/opt/sxq/bin:/opt/sxq/libexec"
         fi
      fi

#     SxQ:
      AC_PATH_PROG([SXPUBLISHER], [sxpublisher], [], [$SXQPATH])
      AC_PATH_PROG([SXCLIDOC],    [sxclidoc],    [], [$SXQPATH])
      AC_PATH_PROG([SXSTDDOC],    [sxstddoc],    [], [$SXQPATH])
      AC_PATH_PROG([SXSTDVIM],    [sxstdvim],    [], [$SXQPATH])

      test -z "$ac_cv_path_SXPUBLISHER" \
         && AC_MSG_ERROR([SxQDoc command sxpublisher not found.])
      test -z "$ac_cv_path_SXCLIDOC" \
         && AC_MSG_ERROR([SxQDoc command sxclidoc not found.])
      test -z "$ac_cv_path_SXSTDDOC" \
         && AC_MSG_ERROR([SxQDoc command sxstddoc not found.])
      test -z "$ac_cv_path_SXSTDVIM" \
         && AC_MSG_ERROR([SxQDoc command sxstdvim not found.])

      sxqBinDir=$(dirname $ac_cv_path_SXPUBLISHER)
      if test -d "$sxqBinDir/../share/manual"; then
         SXQDIR=$(cd $sxqBinDir/..; pwd)
      elif test -d "$sxqBinDir/../../share/manual"; then
         SXQDIR=$(cd $sxqBinDir/../..; pwd)
      else
         AC_MSG_ERROR([Could not find SxQ share folder.])
      fi

#     LaTeX:
      AC_PATH_PROG([LATEX],       [latex],       [])

      test -z "$ac_cv_path_LATEX" \
         && AC_MSG_ERROR([LaTeX interpreter not found.])

      SXQ_LATEX_CLASS(book)
      SXQ_LATEX_PKG(fontenc,book)
      SXQ_LATEX_PKG(lmodern,book)
      SXQ_LATEX_PKG(xcolor,book)
      SXQ_LATEX_PKG(calc,book)
      SXQ_LATEX_PKG(enumitem,book)
      SXQ_LATEX_PKG(hyperref,book)
      SXQ_LATEX_PKG(graphicx,book)
      SXQ_LATEX_PKG(geometry,book)
      SXQ_LATEX_PKG(framed,book)
      SXQ_LATEX_PKG(lineno,book)

#     PS, DVI
      AC_PATH_PROG([DVIPS],  [dvips], [])
      test -z "$ac_cv_path_DVIPS" && AC_MSG_ERROR([dvips found.])

      AC_PATH_PROG([PS2PDF], [ps2pdf], [])
      test -z "$ac_cv_path_PS2PDF" && AC_MSG_ERROR([ps2pdf found.])

      useSxQDoc="yes"

   fi

   AM_CONDITIONAL([COND_BUILD_SXDOC], [test x"$useSxQDoc" = x"yes"])


#  define variables in header
   AC_DEFINE_UNQUOTED(SXQDIR, "$SXQDIR", [SxQ directory])

#  substitute variables in makefile
   AC_SUBST(SXQDIR)
])
