REVAL=TRUE

all: pdf clean

pdf: CHNOSZ-formation.pdf CHNOSZ-intro.pdf CHNOSZ-temperature.pdf

%.tex:	%.Rnw
	echo "Sweave(\"$<\", debug=TRUE, eval=$(REVAL))" | R --slave

CHNOSZ-formation.pdf: CHNOSZ-formation.tex
	latex CHNOSZ-formation
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
	rm -f *.png *.dvi *.aux *.toc *.blg *.log *.out *.snm *.bbl *.nav
	rm -f Gprotein*.ps
	rm -f Zprotein*.ps
	rm -f Gresidue*.ps
	rm -f O2*.ps
	rm -f Eh*.ps
	rm -f T-pH*.ps
	rm -f 80*.ps
	rm -f Rplots.pdf
	rm -f CHNOSZ-intro-*.eps
	rm -f CHNOSZ-intro-*.pdf
	rm -f CHNOSZ-temperature-*.eps
	rm -f CHNOSZ-temperature-*.pdf
