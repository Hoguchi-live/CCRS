#include "arithmetic.h"

void SW_j_invariant(fq_t *output, SW_curve *E) {

	fq_ctx_t F;					// base field
	//fmpz_t p;					// base field characteristic
	fmpz_t fmpz_3, fmpz_2;				// exponents as fmpz_t integers
	fq_t delta, j_invariant;	 		// discriminant, j-invariant
	fq_t tmp1, tmp2, tmp3;				// temporary variables
	fq_t fq_minus_16, fq_4, fq_27, fq_minus_1728;	// constants as fq_t finite field elements

	memcpy(&F, E->F, sizeof(fq_ctx_t));

	// init constants
	fmpz_init(fmpz_2);
	fmpz_init(fmpz_3);

	fq_init(fq_minus_16, F);
	fq_init(fq_4, F);
	fq_init(fq_27, F);
	fq_init(fq_minus_1728, F);

	fq_init(tmp1, F);
	fq_init(tmp2, F);
	fq_init(tmp3, F);

	fq_init(delta, F);
	fq_init(j_invariant, F);

	// set constants
	fmpz_set_si(fmpz_2, 2);
	fmpz_set_si(fmpz_3, 3);
	fq_set_si(fq_minus_16, -16, *(E->F));
	fq_set_si(fq_4, 4, F);
	fq_set_si(fq_27, 27, F);
	fq_set_si(fq_minus_1728, -1728, F);

	//fmpz_init_set(p, fq_ctx_prime(F));

	// compute delta
	fq_pow(tmp1, E->a,  fmpz_3, F);
	fq_mul(tmp1, tmp1, fq_4, F);

	fq_pow(tmp2, E->b,  fmpz_2, F);
	fq_mul(tmp2, tmp2, fq_27, F);

	fq_add(delta, tmp1, tmp2, F);
	fq_mul(delta, delta, fq_minus_16, F);

	// compute j
	fq_mul(tmp1, tmp1, fq_4, F);
	fq_mul(tmp1, tmp1, fq_4, F);
	fq_mul(j_invariant, fq_minus_1728, tmp1, F);
	fq_inv(tmp3, delta, F);
	fq_mul(j_invariant, j_invariant, tmp3, F);

	memcpy(output, &j_invariant, sizeof(fq_t));

	// clear memory
	fmpz_clear(fmpz_2);
	fmpz_clear(fmpz_3);

	fq_clear(fq_minus_16, F);
	fq_clear(fq_4, F);
	fq_clear(fq_27, F);
	fq_clear(fq_minus_1728, F);

	fq_clear(tmp1, F);
	fq_clear(tmp2, F);
	fq_clear(tmp3, F);

	fq_clear(delta, F);
	fq_clear(j_invariant, F);
}
