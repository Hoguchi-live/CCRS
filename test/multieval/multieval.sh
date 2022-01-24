gcc \
	../../src/EllipticCurves/memory.c \
	../../src/EllipticCurves/models.c \
	../../src/EllipticCurves/pretty_print.c \
	../../src/EllipticCurves/arithmetic.c \
	../../src/EllipticCurves/auxiliary.c \
	../../src/Polynomials/binary_trees.c \
	../../src/Polynomials/roots.c \
	../../src/Isogeny/radical.c \
	../../src/Isogeny/velu.c \
	../../src/Isogeny/walk.c \
	multieval.c \
	-g -lgmp -lflint -o multieval && ./multieval
