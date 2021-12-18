/// @file arithmetic.c
#include "arithmetic.h"

/**
  Sets output to 1 if point P belongs to the underlying curve and 0 otherwise.
TODO: check if point P is infinity and otherwise work in affine.
*/
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
	fq_sub_ui(tmp1, tmp1, 3, *(E->F));
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

/*************************************
  Points on elliptic curves
*************************************/
/**
  Sets output to 1 if point P=[0, 1, 0] and 0 otherwise.
  Warning: this does not check whether P is a valid projective point or not.
*/
void SW_point_isinfinity(bool *output, SW_point *P) {

	*output = fq_is_zero(P->x, *(P->E->F)) && fq_is_zero(P->z, *(P->E->F));
}

/**
  Sets output to 1 if point P=(1, 0) and 0 otherwise.
  Warning: this does not check whether P is a valid projective point or not.
*/
void MG_point_isinfinity(bool *output, MG_point *P) {

	*output = fq_is_zero(P->Z, *(P->E->F));
}

/**
  Returns -1 if P is not a projective point.
  If it is, sets output to 1 if point P belongs to the underlying curve and 0 otherwise.
  Finally returns 0.
*/
int SW_point_isvalid(bool *output, SW_point *P) {

	// Check if P is projective
	if(fq_is_zero(P->x, *(P->E->F)) && fq_is_zero(P->y, *(P->E->F)) && fq_is_zero(P->z, *(P->E->F))) return -1;

	const fq_ctx_t *F;
	fq_t res, tmp1, tmp2, tmp3, tmp4;

	F = P->E->F;
	fq_init(res, *F);
	fq_init(tmp1, *F);
	fq_init(tmp2, *F);
	fq_init(tmp3, *F);
	fq_init(tmp4, *F);

	// Y^2 * Z
	fq_pow_ui(tmp1, P->y, 2, *F);
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

	return 0;
}

/**
  Returns -1 if P is not a projective point.
  If it is, sets output to 1 if point P belongs to the underlying curve and 0 otherwise.
  Finally returns 0.
*/
int MG_point_isvalid(bool *output, MG_point *P) {

	// Check if P is projective
	if(fq_is_zero(P->X, *(P->E->F)) && fq_is_zero(P->Z, *(P->E->F))) return -1;
	return 0;
}

void SW_point_rand(SW_point *P) {

	fq_t x, tmp1, tmp2;
	const fq_ctx_t *F;

	F = P->E->F;

	fq_init(x, *F);
	fq_init(tmp1, *F);
	fq_init(tmp2, *F);

	flint_rand_t state;
	flint_randinit(state);

	flint_randclear(state);
}

void MG_point_rand(MG_point *P) {

}

/******************************
  Montgomery Arithmetics
******************************/

void MG_xADD(MG_point *output, MG_point P, MG_point Q, MG_point D) {

	const fq_ctx_t *F;
	F = (P.E)->F;

	fq_t v0, v1, v2, v3;
	fq_init(v0, *F);
	fq_init(v1, *F);
	fq_init(v2, *F);
	fq_init(v3, *F);

	fq_add(v0, P.X, P.Z, *F);
	fq_sub(v1, Q.X, Q.Z, *F);
	fq_mul(v1, v1, v0, *F);
	fq_sub(v0, P.X, P.Z, *F);
	fq_add(v2, Q.X, Q.Z, *F);
	fq_mul(v2, v2, v0, *F);
	fq_add(v3, v1, v2, *F);
	fq_sqr(v3, v3, *F);
	fq_sub(v1, v1, v2, *F);
	fq_sqr(v1, v1, *F);
	fq_mul(output->X, D.Z, v3, *F);
	fq_mul(output->Z, D.X, v1, *F);

	// clear memory
	fq_clear(v0, *F);
	fq_clear(v1, *F);
	fq_clear(v2, *F);
	fq_clear(v3, *F);
}

void MG_xDBL(MG_point *output, MG_point P) {

	const fq_ctx_t *F;
	F = (P.E)->F;

	fq_t v1, v2, v3;
	fq_init(v1, *F);
	fq_init(v2, *F);
	fq_init(v3, *F);

	fq_add(v1, P.X, P.Z, *F);
	fq_sqr(v1, v1, *F);
	fq_sub(v2, P.X, P.Z, *F);
	fq_sqr(v2, v2, *F);
	fq_mul(output->X, v1, v2, *F);
	fq_sub(v1, v1, v2, *F);
	fq_add_ui(v3, (P.E)->A, 2, *F);
	fq_div_ui(v3, v3, 4, *F);
	fq_mul(v3, v3, v1, *F);
	fq_add(v3, v3, v2, *F);
	fq_mul(output->Z, v1, v3, *F);

	// clear memorY
	fq_clear(v1, *F);
	fq_clear(v2, *F);
	fq_clear(v3, *F);
}

void MG_ladder_rec(MG_point *X0, MG_point *X1, fmpz_t k, MG_point P, const fq_ctx_t *F) {

	// base case
	if (fmpz_is_one(k)) {
		fq_set(X0->X, P.X, *F);
		fq_set(X0->Z, P.Z, *F);
		MG_xDBL(X1, P);
	}

	fmpz_t rem;
	fmpz_init(rem);
	fmpz_t two;
	fmpz_init_set_ui(two, 2);
	fmpz_tdiv_qr(k, rem, k, two);	// the value of k is modified here
	fmpz_clear(two);

	MG_ladder_rec(X0, X1, k, P, F); // recursive call

	if (fmpz_is_zero(rem)) {
		// copy *x0 into tmp
		MG_point *tmp;
		MG_point_init(tmp, P.E);
		fq_set(tmp->X, X0->X, *F);
		fq_set(tmp->Z, X0->Z, *F);

		MG_xDBL(X0, *tmp);
		MG_xADD(X1, *tmp, *X1, P);

		MG_point_clear(tmp);
	}

	else {
		MG_xADD(X0, *X0, *X1, P);
		MG_xDBL(X1, *X1);
	}

	// clear memory
	fmpz_clear(rem);
}

void MG_ladder(MG_point *X0, fmpz_t k, MG_point P) {

	const fq_ctx_t *F;
	F = (P.E)->F;

	MG_point *X1;
	MG_point_init(X1, P.E);

	MG_ladder_rec(X0, X1, k, P, F);

	MG_point_clear(X1);
}

void MG_ladder_iter(MG_point *X0, MG_point *X1, fmpz_t k, MG_point P, fq_ctx_t *F) {

	fq_set(X0->X, P.X, *F);
	fq_set(X0->Z, P.Z, *F);
	MG_xDBL(X1, P);

	ulong l;
	l = fmpz_sizeinbase(k, 2);

	for (ulong i = l-2; l>=0; l--) {
		if (fmpz_tstbit(k, i)) {
			MG_xADD(X0, *X0, *X1, P);
			MG_xDBL(X1, *X1);
		}
		else {
			MG_xADD(X1, *X0, *X1, P);
			MG_xDBL(X0, *X0);
		}
	}
}

