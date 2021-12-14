gcc 	../src/EllipticCurves/memory.c \
	../src/EllipticCurves/models.c \
	../src/EllipticCurves/pretty_print.c \
	../src/EllipticCurves/arithmetic.c \
	../src/EllipticCurves/auxiliary.c \
	main.c \
	-lgmp -lflint -o main && ./main
