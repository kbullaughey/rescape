doc = regevscape-manual
all: $(doc).pdf 

$(doc).pdf: $(doc).tex refs.bib plos2009.bst
	echo "s" | pdflatex $<
	echo "s" | pdflatex $<
	bibtex $(doc)
	bibtex $(doc)
	echo "s" | pdflatex $<
	pdflatex $<
