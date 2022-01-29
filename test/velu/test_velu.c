#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../../src/EllipticCurves/models.h"
#include "../../src/EllipticCurves/memory.h"
#include "../../src/EllipticCurves/arithmetic.h"
#include "../../src/EllipticCurves/pretty_print.h"

#include "../../src/Polynomials/binary_trees.h"
#include "../../src/Polynomials/roots.h"

#include "../../src/Isogeny/radical.h"
#include "../../src/Isogeny/walk.h"

#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fq.h>
#include <flint/fq_poly.h>
#include <flint/fq_poly_factor.h>

#define BASE_p "12037340738208845034383383978222801137092029451270197923071397735408251586669938291587857560356890516069961904754171956588530344066457839297755929645858769"
#define BASE_A	"10861338504649280383859950140772947007703646408372831934324660566888732797778932142488253565145603672591944602210571423767689240032829444439469242521864171"
#define BASE_B	"1"
#define BASE_t	"-147189550172528104900422131912266898599387555512924231762107728432541952979290"

int main() {

	// Elliptic curve parameters
	char base_p_str[] = BASE_p;
	slong d = 1;
	fmpz_t base_p;
	fmpz_init(base_p);
	fmpz_set_str(base_p, base_p_str, 0);
	// Finite field context
	fq_ctx_t F;
	char *Fgen = "g";
	fq_ctx_init(F, base_p, d, Fgen);
	/**********************************
		Base Curve
	**********************************/
	MG_curve_t E;
	MG_curve_init(&E, &F);
	MG_curve_set_str(&E, &F, BASE_A, BASE_B, 10);

	/**********************************
		Walk on the 5-graph
	**********************************/
	int ec;
	fmpz_t k, l;
	fq_t j_inv;
	MG_curve_t E_tmp1;

	fmpz_init_set_ui(l, 17);
	fmpz_init_set_ui(k, 1);
	fq_init(j_inv, F);
	MG_curve_init(&E_tmp1, &F);

	/********************************
	  	CLOCK START
	********************************/

	clock_t start = clock(), diff;
	//// Walk
	ec = walk_velu(&E_tmp1, &E, l, k);

	/********************************
		CLOCK STOP
	********************************/
	diff = clock() - start;
	int msec = diff * 1000 / CLOCKS_PER_SEC;
	printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
	printf("walk_velu returned: %d\n", ec);

	////// Target curve
	printf("target curve: ");
	MG_curve_print(&E_tmp1);
	printf("\n");

	//// J-invariant
	MG_j_invariant(&j_inv, &E_tmp1);
	printf("target curve j-invariant: ");
	fq_print_pretty(j_inv, F);
	printf("\n\n");
	/**********************************
		Clear Memory
	**********************************/
	MG_curve_clear(&E_tmp1);
	fmpz_clear(l);
	fmpz_clear(k);

	fq_ctx_clear(F);
	fmpz_clear(base_p);
	MG_curve_clear(&E);

}

