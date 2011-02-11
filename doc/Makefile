doc = regevscape-manual
all: $(doc).pdf 

$(doc).pdf: $(doc).tex refs.bib plos2009.bst
	echo "s" | pdflatex $<
	echo "s" | pdflatex $<
	bibtex $(doc)
	bibtex $(doc)
	echo "s" | pdflatex $<
	pdflatex $<

clean:
	-rm $(doc).out $(doc).aux $(doc).bbl $(doc).blg $(doc).log $(doc).pdf
