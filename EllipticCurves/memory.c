#include <stdio.h>
#include <stdlib.h>

#include "memory.h"
#include "models.h"

/*
   Short Weierstrass elliptic curve
*/

void SW_curve_init(SW_curve *E, const fq_ctx_t *F) {

	fq_init(E->a, *F);
	fq_init(E->b, *F);
}

void SW_curve_set(SW_curve *E, const fq_ctx_t *F, const fq_t a, const fq_t b) {

	E->F = F;
	fq_set(E->a, a, *F);
	fq_set(E->b, b, *F);
}

void SW_curve_set_si(SW_curve *E, const fq_ctx_t *F, const ulong a, const ulong b) {

	fq_t aa, bb;

	fq_init(aa, *F);
	fq_init(bb, *F);

	fq_set_si(aa, a, *F);
	fq_set_si(bb, b, *F);

	SW_curve_set(E, F, aa, bb);
}

void SW_curve_clear(SW_curve *E) {

	fq_clear(E->a, *(E->F));
	fq_clear(E->b, *(E->F));
}

// Points on elliptic curves
void SW_point_init(SW_point *P, const fq_ctx_t *F) {

	fq_init(P->x, *F);
	fq_init(P->y, *F);
	fq_init(P->z, *F);

	P->E = NULL;
}

void SW_point_set(SW_point *P, const fq_t x, const fq_t y, const fq_t z, SW_curve *E) {

	fq_set(P->x, x, *(E->F));
	fq_set(P->y, y, *(E->F));
	fq_set(P->z, z, *(E->F));

	P->E = E;
}

void SW_point_set_si(SW_point *P, const ulong x, const ulong y, const ulong z, SW_curve *E) {

	fq_t xx, yy, zz;

	fq_init(xx, *(E->F));
	fq_init(yy, *(E->F));
	fq_init(zz, *(E->F));

	fq_set_si(xx, x, *(E->F));
	fq_set_si(yy, y, *(E->F));
	fq_set_si(zz, z, *(E->F));

	SW_point_set(P, xx, yy, zz, E);
}

void SW_point_clear(SW_point *P) {

	fq_clear(P->x, *(P->E->F));
	fq_clear(P->y, *(P->E->F));
	fq_clear(P->z, *(P->E->F));
}


/*
   Montgomery elliptic curve
*/

void MG_curve_init(MG_curve *E, const fq_ctx_t *F) {

	fq_init(E->A, *F);
	fq_init(E->B, *F);
}

void MG_curve_set(MG_curve *E, const fq_ctx_t *F, const fq_t A, const fq_t B) {

	E->F = F;
	fq_set(E->A, A, *F);
	fq_set(E->B, B, *F);
}

void MG_curve_set_si(MG_curve *E, const fq_ctx_t *F, const ulong A, const ulong B) {

	fq_t AA, BB;

	fq_init(AA, *F);
	fq_init(BB, *F);

	fq_set_si(AA, A, *F);
	fq_set_si(BB, B, *F);

	MG_curve_set(E, F, AA, BB);
}

void MG_curve_clear(MG_curve *E) {

	fq_clear(E->A, *(E->F));
	fq_clear(E->B, *(E->F));
}
