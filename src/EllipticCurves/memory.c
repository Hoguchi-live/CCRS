/// @file memory.c
#include "memory.h"

/*****************************************
   Short Weierstrass memory management
*****************************************/

/**
  Initializes E for use, with context F, and setting its coefficients to zero.
  A corresponding call to SW_curve_clear() must be made after finishing with the SW_curve to free the memory used by the curve.
*/
void SW_curve_init(SW_curve *E, const fq_ctx_t *F) {

	fq_init(E->a, *F);
	fq_init(E->b, *F);
}

/**
 Sets E to elliptic curve over F in Weierstrass form with coefficients a and b.
 Curve parameters are given as elements of F.
*/
void SW_curve_set(SW_curve *E, const fq_ctx_t *F, const fq_t a, const fq_t b) {

	E->F = F;
	fq_set(E->a, a, *F);
	fq_set(E->b, b, *F);
}

/**
 See SW_curve_set().
 Curve coefficients are given as signed integers.
*/
void SW_curve_set_si(SW_curve *E, const fq_ctx_t *F, const slong a, const slong b) {

	fq_t aa, bb;

	fq_init(aa, *F);
	fq_init(bb, *F);

	fq_set_si(aa, a, *F);
	fq_set_si(bb, b, *F);

	SW_curve_set(E, F, aa, bb);
}

/**
 See SW_curve_set().
 Curve coefficients are given as unsigned integers.
*/
void SW_curve_set_su(SW_curve *E, const fq_ctx_t *F, const ulong a, const ulong b) {

	fq_t aa, bb;

	fq_init(aa, *F);
	fq_init(bb, *F);

	fq_set_ui(aa, a, *F);
	fq_set_ui(bb, b, *F);

	SW_curve_set(E, F, aa, bb);

	fq_clear(aa, *F);
	fq_clear(bb, *F);
}


/**
 See SW_curve_set().
 Parameters are given as null-terminated strings, in base b. The base b can vary between 2 and 62, inclusive. Returns 0 if the string contain valid inputs and −1 otherwise.
*/
int SW_curve_set_str(SW_curve *E, const fq_ctx_t *F, const char *str_a, const char *str_b, int b) {

	fmpz_t fmpz_a, fmpz_b;
	fq_t fq_a, fq_b;
	int ret;

	fmpz_init(fmpz_a);
	fmpz_init(fmpz_b);
	fq_init(fq_a, *F);
	fq_init(fq_b, *F);

	ret = fmpz_set_str(fmpz_a, str_a, b);
	if(ret) return -1;
	ret = fmpz_set_str(fmpz_b, str_b, b);
	if(ret) return -1;

	fq_set_fmpz(fq_a, fmpz_a, *F);
	fq_set_fmpz(fq_b, fmpz_b, *F);

	SW_curve_set(E, F, fq_a, fq_b);

	fmpz_clear(fmpz_a);
	fmpz_clear(fmpz_b);
	fq_clear(fq_a, *F);
	fq_clear(fq_b, *F);

	return 0;
}

/**
  Clears the given curve, releasing any memory used. It must be reinitialised in order to be used again.
*/
void SW_curve_clear(SW_curve *E) {

	fq_clear(E->a, *(E->F));
	fq_clear(E->b, *(E->F));
}

/*********************************************
   Short Weierstrass points memory management
**********************************************/

/**
  Initializes P for use, with context F, and setting its coefficients to zero.
  A corresponding call to SW_point_clear() must be made after finishing with the SW_point to free the memory used by the curve.
TODO: swap F for E in parameters. A point is member of E not of F.
*/
void SW_point_init(SW_point *P, SW_curve *E) {

	fq_init(P->x, *(E->F));
	fq_init(P->y, *(E->F));
	fq_init(P->z, *(E->F));

	P->E = E;
}

/**
 Sets P to point of elliptic curve E with coordinates x, y, z.
 Point parameters are given as elements of F.
TODO: Check if P and E's fields are correct
*/
void SW_point_set(SW_point *P, const fq_t x, const fq_t y, const fq_t z, SW_curve *E) {

	fq_set(P->x, x, *(E->F));
	fq_set(P->y, y, *(E->F));
	fq_set(P->z, z, *(E->F));

	P->E = E;
}

/**
 See SW_point_set().
 Point coordinates are given as signed integers.
*/
void SW_point_set_si(SW_point *P, const slong x, const slong y, const slong z, SW_curve *E) {

	fq_t xx, yy, zz;

	fq_init(xx, *(E->F));
	fq_init(yy, *(E->F));
	fq_init(zz, *(E->F));

	fq_set_si(xx, x, *(E->F));
	fq_set_si(yy, y, *(E->F));
	fq_set_si(zz, z, *(E->F));

	SW_point_set(P, xx, yy, zz, E);

	fq_clear(xx, *(E->F));
	fq_clear(yy, *(E->F));
	fq_clear(zz, *(E->F));
}

/**
 See SW_point_set().
 Point coordinates are given as unsigned integers.
*/
void SW_point_set_ui(SW_point *P, const ulong x, const ulong y, const ulong z, SW_curve *E) {

	fq_t xx, yy, zz;

	fq_init(xx, *(E->F));
	fq_init(yy, *(E->F));
	fq_init(zz, *(E->F));

	fq_set_ui(xx, x, *(E->F));
	fq_set_ui(yy, y, *(E->F));
	fq_set_ui(zz, z, *(E->F));

	SW_point_set(P, xx, yy, zz, E);

	fq_clear(xx, *(E->F));
	fq_clear(yy, *(E->F));
	fq_clear(zz, *(E->F));
}

/**
  Clears the given point, releasing any memory used. It must be reinitialised in order to be used again.
*/
void SW_point_clear(SW_point *P) {

	fq_clear(P->x, *(P->E->F));
	fq_clear(P->y, *(P->E->F));
	fq_clear(P->z, *(P->E->F));
}


/**************************************
   Montgomery curves memory management
**************************************/

/**
  Initializes E for use, with context F, and setting its coefficients to zero.
  A corresponding call to MG_curve_clear() must be made after finishing with the MG_curve to free the memory used by the curve.
*/
void MG_curve_init(MG_curve *E, const fq_ctx_t *F) {

	fq_init(E->A, *F);
	fq_init(E->B, *F);
}

/**
 Sets E to elliptic curve over F in Montgomery form with coefficients A and B.
 Curve parameters are given as elements of F.
*/
void MG_curve_set(MG_curve *E, const fq_ctx_t *F, const fq_t A, const fq_t B) {

	E->F = F;
	fq_set(E->A, A, *F);
	fq_set(E->B, B, *F);
}

/**
 See MG_curve_set().
 Curve coefficients are given as signed integers.
*/
void MG_curve_set_si(MG_curve *E, const fq_ctx_t *F, const slong A, const slong B) {

	fq_t AA, BB;

	fq_init(AA, *F);
	fq_init(BB, *F);

	fq_set_si(AA, A, *F);
	fq_set_si(BB, B, *F);

	MG_curve_set(E, F, AA, BB);

	fq_clear(AA, *F);
	fq_clear(BB, *F);
}

/**
 See MG_curve_set().
 Curve coefficients are given as unsigned integers.
*/
void MG_curve_set_ui(MG_curve *E, const fq_ctx_t *F, const ulong A, const ulong B) {

	fq_t AA, BB;

	fq_init(AA, *F);
	fq_init(BB, *F);

	fq_set_ui(AA, A, *F);
	fq_set_ui(BB, B, *F);

	MG_curve_set(E, F, AA, BB);

	fq_clear(AA, *F);
	fq_clear(BB, *F);
}

/**
 See MG_curve_set().
 Parameters are given as null-terminated strings, in base b. The base b can vary between 2 and 62, inclusive. Returns 0 if the string contain valid inputs and −1 otherwise.
*/
int MG_curve_set_str(MG_curve *E, const fq_ctx_t *F, const char *str_A, const char *str_B, int b) {

	fmpz_t fmpz_A, fmpz_B;
	fq_t fq_A, fq_B;
	int ret;

	fmpz_init(fmpz_A);
	fmpz_init(fmpz_B);
	fq_init(fq_A, *F);
	fq_init(fq_B, *F);

	ret = fmpz_set_str(fmpz_A, str_A, b);
	if(ret) return -1;
	ret = fmpz_set_str(fmpz_B, str_B, b);
	if(ret) return -1;

	fq_set_fmpz(fq_A, fmpz_A, *F);
	fq_set_fmpz(fq_B, fmpz_B, *F);

	MG_curve_set(E, F, fq_A, fq_B);

	fmpz_clear(fmpz_A);
	fmpz_clear(fmpz_B);
	fq_clear(fq_A, *F);
	fq_clear(fq_B, *F);

	return 0;
}

/**
  Clears the given curve, releasing any memory used. It must be reinitialised in order to be used again.
*/
void MG_curve_clear(MG_curve *E) {

	fq_clear(E->A, *(E->F));
	fq_clear(E->B, *(E->F));
}

/*********************************************
   Montgomery curve points memory management
**********************************************/

/**
  Initializes P for use, with context F, and setting its coefficients to zero.
  A corresponding call to SW_point_clear() must be made after finishing with the SW_point to free the memorZ used bZ the curve.
TODO: swap F for E in parameters. A point is member of E not of F.
*/
void MG_point_init(MG_point *P, MG_curve *E) {

	fq_init(P->X, *(E->F));
	fq_init(P->Z, *(E->F));

	P->E = E;
}

/**
 Sets P to point of elliptic curve E with coordinates X, Z, z.
 Point parameters are given as elements of F.
TODO: Check if P and E's fields are correct
*/
void MG_point_set(MG_point *P, const fq_t X, const fq_t Z, MG_curve *E) {

	fq_set(P->X, X, *(E->F));
	fq_set(P->Z, Z, *(E->F));

	P->E = E;
}

/**
 See MG_point_set().
 Point coordinates are given as signed integers.
*/
void MG_point_set_si(MG_point *P, const slong X, const slong Z, MG_curve *E) {

	fq_t XX, ZZ, zz;

	fq_init(XX, *(E->F));
	fq_init(ZZ, *(E->F));

	fq_set_si(XX, X, *(E->F));
	fq_set_si(ZZ, Z, *(E->F));

	MG_point_set(P, XX, ZZ, E);

	fq_clear(XX, *(E->F));
	fq_clear(ZZ, *(E->F));
}

/**
 See MG_point_set().
 Point coordinates are given as unsigned integers.
*/
void MG_point_set_ui(MG_point *P, const ulong X, const ulong Z, MG_curve *E) {

	fq_t XX, ZZ;

	fq_init(XX, *(E->F));
	fq_init(ZZ, *(E->F));

	fq_set_ui(XX, X, *(E->F));
	fq_set_ui(ZZ, Z, *(E->F));

	MG_point_set(P, XX, ZZ, E);

	fq_clear(XX, *(E->F));
	fq_clear(ZZ, *(E->F));
}

/**
  Clears the given point, releasing anZ memorZ used. It must be reinitialised in order to be used again.
*/
void MG_point_clear(MG_point *P) {

	fq_clear(P->X, *(P->E->F));
	fq_clear(P->Z, *(P->E->F));
}

