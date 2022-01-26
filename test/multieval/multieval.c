#include <stdio.h>
#include <stdlib.h>

#include "../../src/EllipticCurves/models.h"
#include "../../src/EllipticCurves/memory.h"
#include "../../src/EllipticCurves/arithmetic.h"
#include "../../src/EllipticCurves/pretty_print.h"

#include "../../src/Polynomials/binary_trees.h"
#include "../../src/Polynomials/multieval.h"
#include "../../src/Polynomials/roots.h"

#include "../../src/Isogeny/radical.h"
#include "../../src/Isogeny/walk.h"

//#include "../../src/Exchange/setup.h"

#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fq.h>
#include <flint/fq_poly.h>
#include <flint/fq_poly_factor.h>

int main() {

	fmpz_t p;
	fq_ctx_t F;
	fq_poly_t P;
	fq_t tmp;
	fq_t *roots;
	fq_t *res;

	uint len = 6;
	fmpz_init_set_ui(p, 7);
	fq_ctx_init(F, p, 1, "g");
	fq_init(tmp, F);
	fq_poly_init(P, F);

	fq_set_ui(tmp, 1, F);
	fq_poly_set_coeff(P, 0, tmp, F);
	fq_poly_set_coeff(P, 1, tmp, F);
	fq_poly_set_coeff(P, 2, tmp, F);

	res = malloc(sizeof(fq_t)*len);
	roots = malloc(sizeof(fq_t)*len);

	for(int i=0; i < len; i++) {
		fq_init(roots[i], F);
		fq_set_ui(roots[i], i, F);

		fq_init(res[i], F);
	}

	fq_poly_multieval(res, roots, P, len, &F);

	//// Check results
	for(int i=0; i < len; i++) {
		printf("Evaluating ");
		fq_poly_print_pretty(P, "X", F);
		printf(" at ");
		fq_print_pretty(roots[i], F);
		printf(" yields ");
		fq_print_pretty(res[i], F);
		printf("\n");
	}
}
