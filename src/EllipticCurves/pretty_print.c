#include "pretty_print.h"

/**
  Prints a representation of E to stdout.
*/
void SW_curve_print(SW_curve *E) {

	fmpz_t p;
	long d;

	fmpz_init_set(p, fq_ctx_prime(*(E->F)));
	d = fq_ctx_degree(*(E->F));

	printf("Elliptic curve Y^3 = X^3 + ");
	fq_print_pretty(E->a, *(E->F));
	printf("X + ");
	fq_print_pretty(E->b, *(E->F));
	printf(" in short Weierstrass form over F_p^d with p = ");
	fmpz_print(p);
	printf(" and d = %ld\n", d);

	fmpz_clear(p);
}

/**
  Prints a representation of P to stdout.
*/
void SW_point_print(SW_point *P) {

	printf("Point [");
	fq_print_pretty(P->x, *(P->E->F));
	printf(", ");
	fq_print_pretty(P->y, *(P->E->F));
	printf(", ");
	fq_print_pretty(P->z, *(P->E->F));
	printf("] on ");
	SW_curve_print(P->E);
}

/**
  Prints a representation of E to stdout.
*/
void MG_curve_print(MG_curve *E) {

	fmpz_t p;
	long d;

	fmpz_init_set(p, fq_ctx_prime(*(E->F)));
	d = fq_ctx_degree(*(E->F));

	printf("Elliptic curve ");
	fq_print_pretty(E->B, *(E->F));
	printf("Y^3 = X^3 + ");
	fq_print_pretty(E->A, *(E->F));
	printf("X^2 + X");
	printf(" in Montgomery form over F_p^d with p = ");
	fmpz_print(p);
	printf(" and d = %ld\n", d);

	fmpz_clear(p);
}

