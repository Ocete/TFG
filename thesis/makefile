main.pdf: main.tex imports/titlepage.tex
	pdflatex --shell-escape main.tex
	bibtex main
	pdflatex --shell-escape main.tex
	pdflatex --shell-escape main.tex

clean:
	rm main.pdf *.aux *.blg *.log *.bbl *.toc *.out
	rm -rf _minted-main
