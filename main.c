#include <stdio.h>
#include <stdlib.h>

#include "EllipticCurves/models.h"
#include "EllipticCurves/memory.h"
#include "EllipticCurves/arithmetic.h"
#include "EllipticCurves/pretty_print.h"

#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fq.h>

int main() {

	// Elliptic curve parameters
	char p_str[] = "0x13";
	slong d = 2;

	fmpz_t p_fmpz;
	fmpz_init(p_fmpz);

	fmpz_set_str(p_fmpz, p_str, 0);

	// Finite field context
	fq_ctx_t F;
	char *Fgen = "g";
	fq_ctx_init(F, p_fmpz, d, Fgen);

	// Base curve
	SW_curve E;
	SW_curve_init(&E, &F);
	SW_curve_set_si(&E, &F, 1, 7);

	// Print base curve
	SW_curve_print(&E);

	// Compute j-invariant
	fq_t j;
	fq_init(j, F);
	SW_j_invariant(&j, &E);
	printf("j-invariant of E: ");
	fq_print_pretty(j, F);
	printf("\n");

	// Test points
	SW_point P;
	SW_point_init(&P, &F);
	SW_point_set_si(&P, 0, 1, 0, &E);

	// clear
	SW_point_clear(&P);
	SW_curve_clear(&E);
	fq_ctx_clear(F);

	return 0;
}

