gcc	../confs/config.c \
	../src/EllipticCurves/memory.c \
	../src/EllipticCurves/models.c \
	../src/EllipticCurves/pretty_print.c \
	../src/EllipticCurves/arithmetic.c \
	../src/EllipticCurves/auxiliary.c \
	../src/Polynomials/binary_trees.c \
	../src/Polynomials/roots.c \
	../src/Isogeny/radical.c \
	../src/Isogeny/velu.c \
	test_rad_5.c \
	-g -lgmp -lflint -o test_rad_5 && ./test_rad_5
