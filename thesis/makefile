main.pdf: main.tex imports/titlepage.tex clean
	pdflatex --shell-escape main.tex
	bibtex main
	pdflatex --shell-escape main.tex
	pdflatex --shell-escape main.tex
	cp main.pdf ../thesis.pdf
	$(MAKE) pseudo_clean

pseudo_clean:
	rm -f *.aux *.blg *.log *.bbl *.toc *.out *.bcf *.xml
	rm -f imports/*.aux
	rm -rf _minted-main

clean: pseudo_clean
	rm -f main.pdf
