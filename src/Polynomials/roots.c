/// @file binary_trees.h

#include "roots.h"

/**
  Extract square root of op by factoring polynomial X^2-op.
  Returns 1 if successful and 0 otherwise.
*/
int fq_sqr_from_polyfact(fq_t rop, fq_t op, const fq_ctx_t F) {

	fq_t one, tmp, lead;
	fq_poly_t pol;
	fq_poly_factor_t fac;

	fq_init(one, F);
	fq_init(tmp, F);
	fq_init(lead, F);
	fq_poly_init(pol, F);
	fq_poly_factor_init(fac, F);

	fq_neg(tmp, op, F);

	fq_set_ui(one, 1, F);
	fq_poly_set_coeff(pol, 2, one, F);
	fq_poly_set_coeff(pol, 0, tmp, F);

	fq_poly_factor(fac, lead, pol, F);

	if(fac->num < 2) return 0;

	fq_poly_get_coeff(rop, fac->poly, 0, F);

	fq_poly_factor_clear(fac, F);
	fq_poly_clear(pol, F);
	fq_clear(lead, F);
	fq_clear(tmp, F);
	fq_clear(one, F);

	return 1;
}
