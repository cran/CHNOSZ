# Makefile for building vignettes in CHNOSZ package,
# used by R CMD build
#
# significant tips from Writing R Extensions:
# - 'make' is called after Sweaving the files
# - the first target is for creation of PDF files
# - then if Makefile has a 'clean' target, 'make clean' is run
#
# this file would be unnecessary except that
# 1) CHNOSZ-formation, which includes graphics in postscript format,
#    doesn't build easily using pdflatex, and
# 2) R CMD build, sans Makefile, builds PDFs using texidvi(pdf=TRUE), 
#    which in turn uses pdflatex
# so, we have to provide rules to build CHNOSZ-formation using dvipdf
#
# because we have this Makefile we also have to provide rules
# for the other vignettes; for those we use pdflatex

# here we do not use
# all: pdf clean
# because 'clean' target might be invoked
# prematurely during parallel make
all: pdf

pdf: CHNOSZ-formation.pdf CHNOSZ-intro.pdf CHNOSZ-temperature.pdf

%.tex:	%.Rnw
	echo "Sweave(\"$<\", debug=TRUE)" | R --slave

CHNOSZ-formation.pdf: CHNOSZ-formation.tex
	latex CHNOSZ-formation
	bibtex CHNOSZ-formation
	latex CHNOSZ-formation
	# run bibtex twice? yes, yes!
	bibtex CHNOSZ-formation
	latex CHNOSZ-formation
	dvipdf CHNOSZ-formation.dvi

%.pdf:	%.tex
	pdflatex $*
	bibtex $*
	pdflatex $*
	bibtex $*
	pdflatex $*

clean:
	# remove files created during LaTeX runs
	rm -f *.png *.dvi *.aux *.toc *.blg *.log *.out *.snm *.bbl *.nav
	# remove figures created by Sweaving CHNOSZ-formation
	rm -f Gprotein*.ps
	rm -f Zprotein*.ps
	rm -f Gresidue*.ps
	rm -f O2*.ps
	rm -f Eh*.ps
	rm -f T-pH*.ps
	rm -f 80*.ps
	rm -f Rplots.pdf
	# remove figures created by Sweaving the other vignettes
	rm -f CHNOSZ-intro-*.eps
	rm -f CHNOSZ-intro-*.pdf
	rm -f CHNOSZ-temperature-*.eps
	rm -f CHNOSZ-temperature-*.pdf
	# don't remove figures that are included by the vignettes,
	# otherwise users with only the built package (e.g. CRAN) 
	# won't be able to build the vignettes themselves
	#rm -f chnosz_new.pdf  # used by CHNOSZ-intro
	#rm -f rxn_*.ps        # used by CHNOSZ-formation
	#rm -f 1TUI_*.ps       # used by CHNOSZ-formation
