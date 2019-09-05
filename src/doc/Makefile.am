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

doc:
	@DOXYGEN@ doxygen.cnf

fastdoc:   
	@DOXYGEN@ doxygen-fast.cnf

wikidoc:
	@DOXYGEN@ doxygen-wiki.cnf

upload:
	ssh adam1.mpie.de rm -rf doxygen/sphinx/html
	scp -r wikidoc adam1.mpie.de:doxygen/sphinx/html
EXTRA_DIST=doxygen.cnf doxygen-fast.cnf doxygen-wiki.cnf _htaccess.ext
