gcc	../confs/config.c \
	../src/EllipticCurves/memory.c \
	../src/EllipticCurves/models.c \
	../src/EllipticCurves/pretty_print.c \
	../src/EllipticCurves/arithmetic.c \
	../src/EllipticCurves/auxiliary.c \
	../src/Polynomials/binary_trees.c \
	main.c \
	-lgmp -lflint -o main && ./main
