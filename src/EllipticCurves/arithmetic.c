#include "arithmetic.h"

// Elliptic curves

void SW_j_invariant(fq_t *output, SW_curve *E) {

	fq_t tmp1, tmp2, tmp3;				// temporary registers
	fq_t delta, j_invariant;	 		// discriminant, j-invariant

	fq_init(tmp1, *(E->F));
	fq_init(tmp2, *(E->F));
	fq_init(tmp3, *(E->F));

	fq_init(delta, *(E->F));
	fq_init(j_invariant, *(E->F));

	// compute delta
	fq_pow_ui(tmp1, E->a, 3, *(E->F));
	fq_mul_si(tmp1, tmp1, 4, *(E->F));

	fq_pow_ui(tmp2, E->b,  2, *(E->F));
	fq_mul_si(tmp2, tmp2, 27, *(E->F));

	fq_add(delta, tmp1, tmp2, *(E->F));
	fq_mul_si(delta, delta, -16, *(E->F));

	// compute j
	fq_mul_si(tmp1, tmp1, 4, *(E->F));
	fq_mul_si(tmp1, tmp1, 4, *(E->F));
	fq_mul_si(j_invariant, tmp1, -1728, *(E->F));
	fq_inv(tmp3, delta, *(E->F));
	fq_mul(j_invariant, j_invariant, tmp3, *(E->F));

	// set output
	fq_set(*output, j_invariant, *(E->F));

	// clear memory
	fq_clear(tmp1, *(E->F));
	fq_clear(tmp2, *(E->F));
	fq_clear(tmp3, *(E->F));
	fq_clear(delta, *(E->F));
	fq_clear(j_invariant, *(E->F));
}


void MG_j_invariant(fq_t *output, MG_curve *E) {

	fq_t tmp1, tmp2;				// temporary registers
	fq_t j_invariant;	 		// discriminant, j-invariant

	fq_init(tmp1, *(E->F));
	fq_init(tmp2, *(E->F));
	fq_init(j_invariant, *(E->F));

	// numerator
	fq_pow_ui(tmp1, E->A, 2, *(E->F));
	fq_add_si(tmp1, tmp1, -3, *(E->F));
	fq_pow_ui(tmp1, tmp1, 3, *(E->F));
	fq_mul_ui(tmp1, tmp1, 256, *(E->F));

	// denominator
	fq_pow_ui(tmp2, E->A, 2, *(E->F));
	fq_sub_ui(tmp2, tmp2, 4, *(E->F));

	fq_div(j_invariant, tmp1, tmp2, *(E->F));

	// set output
	fq_set(*output, j_invariant, *(E->F));

	// clear memory
	fq_clear(tmp1, *(E->F));
	fq_clear(tmp2, *(E->F));
	fq_clear(j_invariant, *(E->F));
}

// Points on elliptic curves

void SW_point_valid(bool *output, SW_point *P) {

	const fq_ctx_t *F;
	fq_t res, tmp1, tmp2, tmp3, tmp4;

	F = P->E->F;
	fq_init(res, *F);
	fq_init(tmp1, *F);
	fq_init(tmp2, *F);
	fq_init(tmp3, *F);
	fq_init(tmp4, *F);

	// Y^2 * Z
	//fq_pow_ui(tmp1, P->y, 2, *F);
	fq_mul(tmp1, tmp1, P->z, *F);
	// X^3
	fq_pow_ui(tmp2, P->x, 3, *F);
	// aX*Z^2
	fq_pow_ui(tmp3, P->z, 2, *F);
	fq_mul(tmp3, tmp3, P->x, *F);
	fq_mul(tmp3, tmp3, P->E->a, *F);
	// bZ^3
	fq_pow_ui(tmp4, P->z, 3, *F);
	fq_mul(tmp4, tmp4, P->E->b, *F);

	fq_add(res, tmp2, tmp3, *F);
	fq_add(res, res, tmp4, *F);
	fq_sub(res, res, tmp1, *F);

	*output = fq_is_zero(res, *F);

	// clear
	fq_clear(res, *F);
	fq_clear(tmp1, *F);
	fq_clear(tmp2, *F);
	fq_clear(tmp3, *F);
	fq_clear(tmp4, *F);
}

// Montgomery curve arithmetic

void MG_xADD(MG_point *output, MG_point P, MG_point Q, MG_point D) {

	const fq_ctx_t *F;
	F = P->E->F;

	fq_t v0, v1, v2, v3;
	fq_init(v0, *F);
	fq_init(v1, *F);
	fq_init(v2, *F);
	fq_init(v3, *F);

	fq_add(v0, P->x, P->z, *F);
	fq_sub(v1, Q->x, Q->z, *F);
	fq_mul(v1, v1, v0, *F);
	fq_sub(v0, P->x, P->z, *F);
	fq_add(v2, Q->x, Q->z, *F);
	fq_mul(v2, v2, v0, *F);
	fq_add(v3, v1, v2, *F);
	fq_sqr(v3, v3, *F);
	fq_sub(v1, v1, v2, *F);
	fq_sqr(v1, v1, *F);
	fq_mul(output->x, D->z, v3, *F);
	fq_mul(output->z, D->x, v1, *F);

	// clear memory
	fq_clear(v0, *F);
	fq_clear(v1, *F);
	fq_clear(v2, *F);
	fq_clear(v3, *F);
}

void MG_xDBL(MG_point *output, MG_point P) {

	const fq_ctx_t *F;
	F = P->E->F;
	

	fq_t v1, v2, v3;
	fq_init(v1, *F);
	fq_init(v2, *F);
	fq_init(v3, *F);

	fq_add(v1, P->x, P->z, *F);
	fq_sqr(v1, v1, *F);
	fq_sub(v2, P->x, P->z, *F);
	fq_sqr(v2, v2, *F);
	fq_mul(output->x, v1, v2, *F);
	fq_sub(v1, v1, v2, *F);
	fq_add_ui(v3, *(P->E->A), 2, *F);
	fq_div_ui(v3, v3, 4, *F);
	fq_mul(v3, v3, v1, *F);
	fq_add(v3, v3, v2, *F);
	fq_mul(output->z, v1, v3, *F)

	// clear memory
	fq_clear(v1, *F);
	fq_clear(v2, *F);
	fq_clear(v3, *F);
}