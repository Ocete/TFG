main.pdf: main.tex imports/titlepage.tex clean
	pdflatex --shell-escape main.tex
	bibtex main
	pdflatex --shell-escape main.tex
	pdflatex --shell-escape main.tex
	cp main.pdf ../tfg.pdf

clean:
	-rm main.pdf *.aux *.blg *.log *.bbl *.toc *.out *.bcf *.xml
	rm -rf _minted-main

pdf_clean:
	-rm *.aux *.blg *.log *.bbl *.toc *.out *.bcf *.xml
	rm -rf _minted-main
