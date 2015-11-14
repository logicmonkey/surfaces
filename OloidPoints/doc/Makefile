../OloidPoints.pdf: *.png *.lyx
	lyx -e "pdflatex" OloidPoints.lyx
	lyx -e "pdf2" OloidPoints.lyx
	mv OloidPoints.pdf ../OloidPoints.pdf

clean:
	rm -f OloidPoints.tex ./OloidPoints.pdf ../OloidPoints.pdf
