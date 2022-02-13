pdflatex -shell-escape main.tex && bibtex main && pdflatex -shell-escape main.tex && zathura ./main.pdf
