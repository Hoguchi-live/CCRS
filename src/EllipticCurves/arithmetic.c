// @file arithmetic.c
#include "arithmetic.h"

/**
  Sets output to 1 if point P belongs to the underlying curve and 0 otherwise.
TODO: check if point P is infinity and otherwise work in affine.
*/
void SW_j_invariant(fq_t *output, SW_curve_t *E) {

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


void MG_j_invariant(fq_t *output, MG_curve_t *E) {

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
void SW_point_isinfinity(bool *output, SW_point_t *P) {

	*output = fq_is_zero(P->x, *(P->E->F)) && fq_is_zero(P->z, *(P->E->F));
}

/**
  Sets output to 1 if point P=(1, 0) and 0 otherwise.
  Warning: this does not check whether P is a valid projective point or not.
*/
void MG_point_isinfinity(bool *output, MG_point_t *P) {

	*output = fq_is_zero(P->Z, *(P->E->F));
}

/**
  Returns -1 if P is not a projective point.
  If it is, sets output to 1 if point P belongs to the underlying curve and 0 otherwise.
  Finally returns 0.
*/
int SW_point_isvalid(bool *output, SW_point_t *P) {

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
TODO: try by square root to determine if point is in predefined extension.
  Returns -1 if P is not a projective point.
  If it is, sets output to 1 if point P belongs to the underlying curve and 0 otherwise.
  Finally returns 0.
*/
int MG_point_isvalid(bool *output, MG_point_t *P) {

	// Check if P is projective
	if(fq_is_zero(P->X, *(P->E->F)) && fq_is_zero(P->Z, *(P->E->F))) return -1;
	return 0;
}


void MG_point_isinfty(bool *output, MG_point_t *P) {

	*output = fq_is_zero(P->Z, *(P->E->F));
}

bool MG_point_isinfty_(MG_point_t *P) {

	return fq_is_zero(P->Z, *(P->E->F));
}

/**
  Normalizes point coordinate to (X/Z, 1) or (1, 0) if P is infinity.
*/
void MG_point_normalize(MG_point_t *P) {

	bool isinfty;
	MG_point_isinfty(&isinfty, P);
	if(!isinfty) {
		fq_div(P->X, P->X, P->Z, *(P->E->F));
		fq_one(P->Z, *(P->E->F));
	}
}

void SW_point_rand_ninfty(SW_point_t *P) {

	fq_t x, y, tmp1, tmp2;
	flint_rand_t state;

	const fq_ctx_t *F = P->E->F;

	fq_init(x, *F);
	fq_init(y, *F);
	fq_init(tmp1, *F);
	fq_init(tmp2, *F);
	flint_randinit(state);

	// Main loop
	int ret = 0;
	while(ret == 0) {
		// Find random x in base field
		fq_randtest(x, state, *F);

		// Compute T := X^3 + aX + b
		fq_pow_ui(tmp1, x, 3, *F);
		fq_mul(tmp2, P->E->a, x, *F);
		fq_add(tmp2, tmp2, P->E->b, *F);
		fq_add(tmp1, tmp1, tmp2, *F);

		// Extract root if exists, otherwise fail with 0.
		ret = fq_sqr_from_polyfact(y, tmp1, *F);
	}

	// Create corresponding SW_point_t, is not infinity
	fq_set(P->x, x, *F);
	fq_set(P->y, y, *F);
	fq_set_ui(P->z, 1, *F);

	flint_randclear(state);
	fq_clear(tmp2, *F);
	fq_clear(tmp1, *F);
	fq_clear(y, *F);
	fq_clear(x, *F);
}

void MG_point_rand_ninfty(MG_point_t *P) {

	fq_t X, Y, tmp1, tmp2;
	flint_rand_t state;

	const fq_ctx_t *F = P->E->F;

	fq_init(X, *F);
	fq_init(Y, *F);
	fq_init(tmp1, *F);
	fq_init(tmp2, *F);
	flint_randinit(state);

	// Main loop
	int ret = 0;
	while(ret == 0) {
		// Find random x in base field
		fq_randtest(X, state, *F);

		// Compute T := B^-1 * x * (x^2 + Ax + 1)
		fq_pow_ui(tmp1, X, 2, *F);
		fq_mul(tmp2, P->E->A, X, *F);
		fq_add_ui(tmp2, tmp2, 1, *F);
		fq_add(tmp1, tmp1, tmp2, *F);
		fq_mul(tmp1, tmp1, X, *F);

		fq_inv(tmp2, P->E->B, *F);
		fq_mul(tmp1, tmp1, tmp2, *F);

		// Extract root if exists, otherwise fail with 0. TODO: JUST CHECK IF IT EXISTS?
		ret = fq_sqr_from_polyfact(tmp2, tmp1, *F);
	}
	// Create corresponding SW_point_t, is not infinity
	fq_set(P->X, X, *F);
	fq_set_ui(P->Z, 1, *F);

	flint_randclear(state);
	fq_clear(tmp2, *F);
	fq_clear(tmp1, *F);
	fq_clear(Y, *F);
	fq_clear(X, *F);
}

/******************************
  Montgomery Arithmetics
******************************/

/**
  Bruteforce recovery of a possible y-coordinate of P into rop.
  Returns 1 if successful, 0 otherwise.
*/
int MG_rec_y(fq_t rop, MG_point_t *P){

	const fq_ctx_t *F = P->E->F;

	// normalize P and check if equals O
	MG_point_normalize(P);
	if(MG_point_isinfty_(P)){
		fq_set_ui(rop, 1, *F);
		return 1;
	}

	fq_t tmp1, tmp2;
	fq_init(tmp1, *F);
	fq_init(tmp2, *F);

	// Get B^-1 * x(x^2 + Ax + 1)
	fq_pow_ui(tmp1, P->X, 2, *F);
	fq_mul(tmp2, P->E->A, P->X, *F);
	fq_add_ui(tmp2, tmp2, 1, *F);
	fq_add(tmp1, tmp1, tmp2, *F);
	fq_mul(tmp1, tmp1, P->X, *F);

	fq_inv(tmp2, P->E->B, *F);
	fq_mul(tmp1, tmp1, tmp2, *F);

	int ret = fq_sqr_from_polyfact(rop, tmp1, *F);

	fq_clear(tmp2, *F);
	fq_clear(tmp1, *F);

	return ret;
}

void MG_xADD(MG_point_t *output, MG_point_t P, MG_point_t Q, MG_point_t D) {

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

void MG_xDBL(MG_point_t *output, MG_point_t P) {

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

	// clear memory
	fq_clear(v1, *F);
	fq_clear(v2, *F);
	fq_clear(v3, *F);
}

void MG_ladder_rec(MG_point_t *X0, MG_point_t *X1, fmpz_t k, MG_point_t P, const fq_ctx_t *F) {

	// base case
	if (fmpz_is_one(k)) {
		fq_set(X0->X, P.X, *F);
		fq_set(X0->Z, P.Z, *F);
		MG_xDBL(X1, P);
	}

	fmpz_t rem;
	fmpz_t two;

	fmpz_init(rem);
	fmpz_init_set_ui(two, 2);
	fmpz_tdiv_qr(k, rem, k, two);	// the value of k is modified here

	fmpz_clear(two);

	MG_ladder_rec(X0, X1, k, P, F); // recursive call

	if (fmpz_is_zero(rem)) {
		// copy *x0 into tmp
		MG_point_t *tmp;
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

void MG_ladder(MG_point_t *X0, fmpz_t k, MG_point_t P) {

	const fq_ctx_t *F;
	F = (P.E)->F;

	MG_point_t X1;
	MG_point_init(&X1, P.E);

	MG_ladder_rec(X0, &X1, k, P, F);

	MG_point_clear(&X1);
}


void MG_ladder_iter_(MG_point_t *rop, fmpz_t k, MG_point_t *op) {
	// Check if k <0
	//TODO

	MG_curve_t *E = op->E;
	const fq_ctx_t *F = E->F;

	// Check if k = 0
	if(fmpz_is_zero(k)){
		fq_one(rop->X, *F);
		fq_zero(rop->Z, *F);
		return;
	}

	// Check if P = O
	bool isinfty;
	MG_point_isinfty(&isinfty, op);
	if(isinfty) return;

	// Buffers
	MG_point_t X0, X1;

	MG_point_init(&X0, E);
	MG_point_init(&X1, E);

	fq_set(X0.X, op->X, *F);
	fq_set(X0.Z, op->Z, *F);
	MG_xDBL(&X1, *op);

	int l;
	l = fmpz_sizeinbase(k, 2);

	for (int i = l-2; i>=0; i--) {
		if (fmpz_tstbit(k, i)) {
			MG_xADD(&X0, X0, X1, *op);
			MG_xDBL(&X1, X1);
		}
		else {
			MG_xADD(&X1, X0, X1, *op);
			MG_xDBL(&X0, X0);
		}
	}

	fq_set(rop->X, X0.X, *F);
	fq_set(rop->Z, X0.Z, *F);

	MG_point_clear(&X0);
	MG_point_clear(&X1);
}

void MG_ladder_iter(MG_point_t *X0, MG_point_t *X1, fmpz_t k, MG_point_t P, fq_ctx_t *F) {
	//TODO NEED P AS POINTER

	// Check if k <0
	//TODO

	// Check if k = 0
	if(fmpz_is_zero(k)){
		fq_one(X0->X, *(P.E->F));
		fq_zero(X0->Z, *(P.E->F));
		return;
	}

	// Check if P = O
	bool isinfty;
	MG_point_isinfty(&isinfty, &P);
	if(isinfty) return;

	fq_set(X0->X, P.X, *F);
	fq_set(X0->Z, P.Z, *F);
	MG_xDBL(X1, P);

	int l;
	l = fmpz_sizeinbase(k, 2);

	for (int i = l-2; i>=0; i--) {
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

/**
   Sets rop to the frobenius' trace for the CRS base curve.
   rop must be initialized.
   One should call PARI to turn this into a real function.
*/
void MG_curve_trace(fmpz_t rop) {

	char trace[] = "-147189550172528104900422131912266898599387555512924231762107728432541952979290";
	fmpz_set_str(rop, trace, 10);
}

/**
   Sets rop to the cardinal of E(F_p) where p is the caracteristic of the base field.
   rop must be initialized.
*/
void MG_curve_card_base(fmpz_t rop, MG_curve_t *E) {

	fmpz_t trace;
	fmpz_init(trace);
	MG_curve_trace(trace);

	fmpz_set(rop, fq_ctx_prime(*(E->F)));
	fmpz_add_ui(rop, rop, 1);
	fmpz_sub(rop, rop, trace);

	fmpz_clear(trace);
}

/**
   Sets rop to the cardinal of E(F_p^r) where p is the caracteristic of the base field.
   rop must be initialized.
   r must be positive and greater than 1.
   ref: https://perso.univ-rennes1.fr/christophe.ritzenthaler/cours/elliptic-curve-course.pdf page 11.
	TODO: Work with q instead of p, problem with quad twist?
*/
void MG_curve_card_ext(fmpz_t rop, MG_curve_t *E, fmpz_t r) {

	fmpz_t t, p, c, s, s_, tmp;

	fmpz_init(t);	// trace
	fmpz_init(p);	// caracteristic
	fmpz_init(c);	// counter
	fmpz_init(s);	// s_n
	fmpz_init(s_);	// s_{n-1}
	fmpz_init(tmp);

	fmpz_set(p, fq_ctx_prime(*(E->F)));
	MG_curve_trace(t);
	fmpz_set(c, r);

	// Initialization
	fmpz_set_ui(s_, 2);
	fmpz_set(s, t);

	while(!fmpz_is_one(c)) {

		// compute new s_n
		fmpz_set(tmp, s);
		fmpz_mul(s, s, t);
		fmpz_mul(s_, s_, p);
		fmpz_sub(s, s, s_);

		// set s_{n-1} to previous s_n
		fmpz_set(s_, tmp);

		fmpz_sub_ui(c, c, 1);
	}
	// Compute cardinal as q^r + 1 - s_
	fmpz_pow_ui(rop, p, fmpz_get_ui(r));
	fmpz_add_ui(rop, rop, 1);
	fmpz_sub(rop, rop, s);

	fmpz_clear(t);
	fmpz_clear(p);
	fmpz_clear(c);
	fmpz_clear(s);
	fmpz_clear(s_);
	fmpz_clear(tmp);
}

/**
   Sets P to a random l-torsion point on the underlying curve and returns 1.
   The point P will be strictly in E(F_q^r).

   not implemented yet:
   If r % 2 == 0, the x-coordinate of P will be in F_q^r//2 (x-only arithmetics).
   Variable card holds the cardinal of E(F_q) and can be computed using MG_curve_card.
   TODO: card will be hardcoded and held in a struct.
   Returns 0 in case of failure (no such point on E).
*/
int MG_curve_rand_torsion(MG_point_t *P, fmpz_t l, fmpz_t r, fmpz_t card) {

	fmpz_t val, l_pow, cofactor, e;
	MG_point_t Q, R;
	bool isinfty = 1;

	fmpz_init(val);
	fmpz_init(cofactor);
	fmpz_init(e);
	MG_point_init(&Q, P->E);
	MG_point_init(&R, P->E);

	MG_point_set_infty(&Q);

	fmpz_val_q(val, cofactor, card, l);
	if(fmpz_is_zero(val)) return 0;

	while(isinfty) {

		MG_point_rand_ninfty(&R);
		MG_ladder_iter_(&Q, cofactor, &R);
		MG_point_isinfty(&isinfty, &Q);
	};

	// Extract l-torsion point from possibly l^val-torsion point.
	// Here R acts as a temporary variable for l*Q
	MG_ladder_iter_(&R, l, &Q);
	MG_point_isinfty(&isinfty, &R);
	fmpz_set_ui(e, 1);

	// While l*Q != O do Q := l*Q
	while(!isinfty && 0 >= fmpz_cmp(e, val)) {
		MG_point_set_(&Q, &R);
		MG_ladder_iter_(&R, l, &Q);
		MG_point_isinfty(&isinfty, &R);

		fmpz_add_ui(e, e, 1);
	}

	MG_point_set_(P, &Q);

	MG_point_clear(&R);
	MG_point_clear(&Q);
	fmpz_clear(e);
	fmpz_clear(cofactor);
	fmpz_clear(val);

}

