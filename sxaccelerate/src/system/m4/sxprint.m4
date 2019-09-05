# --- Define macro for print statements
#     1: title of package
AC_DEFUN([SX_TITLE], [
   cat <<SX_TITLE_END
+----------------------------------------------------------------------------
| Package: [$1]
+----------------------------------------------------------------------------
SX_TITLE_END
])
