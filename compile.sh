gcc 	EllipticCurves/memory.c \
	EllipticCurves/models.c \
	EllipticCurves/pretty_print.c \
	EllipticCurves/arithmetic.c \
	main.c \
	-lgmp -lflint -o main && ./main
