// @file radical.c
#include "radical.h"

/**
  Extract n-th root of op using the following trick:
  	we always have
  	op ^ (p + 1) = op ^ 2
	and assuming
	p + 1 = 2*l*k
	for some natural number k, we get
	(op ^ e)^l = +/- op
	for e = (p + 1) / (2 * l).
  We get rid of the +/- sign by testing via fast exponentiation.
  A NAS condition for this to work is p + 1 = 0 % 2 * l.
  This is about 1000 times faster than bruteforcing X^l - op.
**/
void fq_nth_root_trick(fq_t rop, fq_t op, fmpz_t l, const fq_ctx_t F) {

	fmpz_t tmp, e, p;
	fq_t sgn_check, alpha;

	fq_init(alpha, F);
	fq_init(sgn_check, F);
	fmpz_init(p);
	fmpz_init(e);
	fmpz_init(tmp);
	fmpz_set(p, fq_ctx_prime(F));

	//// Compute e
	fmpz_mul_ui(tmp, l, 2);
	fmpz_add_ui(e, p, 1);
	fmpz_fdiv_q(e, e, tmp);

	//// Compute alpha = op ^ e
	fq_pow(alpha, op, e, F);

	//// Check for sign
	fq_pow(sgn_check, alpha, l, F);
	if(!fq_equal(op, sgn_check, F)) fq_neg(alpha, alpha, F);

	//// Copy buffer
	fq_set(rop, alpha, F);

	fq_clear(alpha, F);
	fq_clear(sgn_check, F);
	fmpz_clear(p);
	fmpz_clear(e);
	fmpz_clear(tmp);
}

void fq_nth_root_trick_ui(fq_t rop, fq_t op, slong l, const fq_ctx_t F) {

	fmpz_t ll;
	fmpz_init_set_ui(ll, l);
	fq_nth_root_trick(rop, op, ll, F);
	fmpz_clear(ll);
}

void radical_isogeny_5(TN_curve_t *rop, TN_curve_t *op, fmpz_t k) {

	fmpz_t l;
	fq_t b, alpha, tmp1, tmp2, tmp3, num, den, res;
	fq_t alpha_pow[4];

	const fq_ctx_t *F = op->F;
	fq_init(b, *F);
	fq_init(alpha, *F);
	fq_init(tmp1, *F);
	fq_init(tmp2, *F);
	fq_init(tmp3, *F);
	fq_init(res, *F);
	fq_init(num, *F);
	fq_init(den, *F);
	fmpz_init_set_ui(l, 5);

	// Init b = op->b
	fq_set(b, op->b, *F);

	// Main loop that goes through k isogeny steps
	for(int step=0; fmpz_cmp_ui(k, step) > 0; step++) {

		//// Extract root of rho = b
		fq_nth_root_trick_ui(alpha, b, 5, *F);

		//// Store alpha ^ i for i = 1 to 4
		for(int i=0; i< 4; i++){
			fq_init(alpha_pow[i], *F);
			if(i == 0) fq_set(alpha_pow[i], alpha, *F);
			else fq_mul(alpha_pow[i], alpha_pow[i-1], alpha, *F);
		}

		//// Precompute 2*alpha, 3*alpha, 4*alpha^2
		fq_mul_ui(tmp1, alpha, 2, *F);
		fq_add(tmp2, tmp1, alpha, *F);
		fq_mul_ui(tmp3, alpha_pow[1], 4, *F);

		// Compute base shared by numerator and denominator: alpha^4 + 4alpha^2 + 1
		fq_add(num, alpha_pow[3], tmp3, *F);
		fq_add_ui(num, num, 1, *F);
		fq_set(den, num, *F);

		//// Finish num/den
		fq_add(num, num, tmp1, *F);
		fq_mul(tmp3, tmp2, alpha_pow[1], *F); // 3alpha * alpha ^ 2
		fq_add(num, num, tmp3, *F);

		fq_sub(den, den, tmp2, *F);
		fq_mul(tmp1, tmp1, alpha_pow[1], *F); // 2alpha * alpha ^ 2
		fq_sub(den, den, tmp1, *F);

		//// Finally compute alpha * num / den
		fq_mul(res, alpha, num, *F);
		fq_inv(tmp1, den, *F);
		fq_mul(res, res, tmp1, *F);

		//// Copy buffer
		fq_set(b, res, *F);
	}
	//// Set curve (here b = c)
	TN_curve_set(rop, b, b, l, F);

	fmpz_clear(l);
	fq_clear(alpha, *F);
	fq_clear(tmp1, *F);
	fq_clear(tmp2, *F);
	fq_clear(tmp3, *F);
	fq_clear(res, *F);
	fq_clear(num, *F);
	fq_clear(den, *F);
	for(int i=0; i< 4; i++) fq_clear(alpha_pow[i], *F);
}

