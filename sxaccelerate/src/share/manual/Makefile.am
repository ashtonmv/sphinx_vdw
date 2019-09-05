SUBDIRS =

MANSECTION  = 1
AUTHOR      = "SxAccelerate Consortium"

SXDOC_SRC   = $(wildcard $(abs_srcdir)/sxdoc/*.sxdoc)
SXDOC_MAIN  = $(abs_srcdir)/sxdoc/Manual.sxdoc

SXPUBDBG    = --no-cleanup

CLEANFILES  = manual_*.* *.aux *.pdf *.log *.out
CLEANFILES += manual.$(MANSETION).gz

manualdir   = $(prefix)/share/manual
manual_DATA = manual_a4.pdf manual_a5.pdf manual.$(MANSECTION).gz

LATEX_ENV  = SXQTMPLDIR=$(abs_top_srcdir)/share/manual

SXQDOC_ENV = 

LATEX       = $(LATEX_ENV) latex
LATEX_ARGS  = -interaction=errorstopmode

manual_a4.pdf: $(SXDOC_SRC)
	$(SXQDOC_ENV) $(SXPUBLISHER) $(SXPUBDBG) --sxdoc $(SXDOC_MAIN) --author \"$(AUTHOR)\" --a4 --sxpath $(abs_srcdir)/sxdoc --output manual

manual_a5.pdf: $(SXDOC_SRC)
	$(SXQDOC_ENV) $(SXPUBLISHER) $(SXPUBDBG) --sxdoc $(SXDOC_MAIN) --author \"$(AUTHOR)\" --a5 --sxpath $(abs_srcdir)/sxdoc --output manual

manual.$(MANSECTION).gz: $(SXDOC_SRC)
	$(SXQDOC_ENV) $(SXPUBLISHER) --sxdoc $(SXDOC_MAIN) --author \"$(AUTHOR)\" --man --section $(MANSECTION) --sxpath $(abs_srcdir)/sxdoc --output manual
