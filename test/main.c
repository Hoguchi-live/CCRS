#include <stdio.h>
#include <stdlib.h>

#include "../src/EllipticCurves/models.h"
#include "../src/EllipticCurves/memory.h"
#include "../src/EllipticCurves/arithmetic.h"
#include "../src/EllipticCurves/pretty_print.h"

#include "../src/Polynomials/binary_trees.h"
#include "../src/Polynomials/roots.h"

#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fq.h>
#include <flint/fq_poly.h>
#include <flint/fq_poly_factor.h>


#define BASE_p "12037340738208845034383383978222801137092029451270197923071397735408251586669938291587857560356890516069961904754171956588530344066457839297755929645858769"
#define BASE_q "12037340738208845034383383978222801137092029451270197923071397735408251586669938291587857560356890516069961904754171956588530344066457839297755929645858769"

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
	MG_curve E;
	MG_curve_init(&E, &F);
	MG_curve_set_str(&E, &F, BASE_A, BASE_B, 10);
	//MG_curve_print(&E);

	/**********************************
		Arithmetics Tests
	**********************************/
	// Compute j-invariant
	//fq_t j;
	//fq_init(j, F);
	//MG_j_invariant(&j, &E);
	//printf("j-invariant of E: ");
	//fq_print_pretty(j, F);
	//printf("\n");

	// Test points
	//SW_point P;
	//SW_point_init(&P, &E);
	//SW_point_set_si(&P, 0, 1, 0, &E);
	//SW_point_print(&P);

	//bool *is_valid;
	//is_valid = malloc(sizeof(bool));
	//SW_point_valid(is_valid, &P);
	//printf("Point P is valid? %d\n", *is_valid);
	//free(is_valid);

	/**********************************
		Weierstrass Points
	**********************************/
	//SW_point P;
	//SW_point_init(&P, &E);
	//SW_point_print(&P);
	//SW_point_clear(&P);

	/**********************************
		Montgomery Points
	**********************************/
	fmpz_t k;
	MG_point P, X0, X1;

	fmpz_init(k);
	MG_point_init(&P, &E);
	MG_point_init(&X0, &E);
	MG_point_init(&X1, &E);

	fmpz_set_ui(k, 1);
	MG_point_rand_ninfty(&P);

	MG_ladder_iter(&X0, &X1, k, P, &F);

	fmpz_clear(k);
	MG_point_clear(&P);
	MG_point_clear(&X1);
	MG_point_clear(&X0);

	/**********************************
		     SQRT
	**********************************/

	/**********************************
		Clear Memory
	**********************************/
	fq_ctx_clear(F);
	fmpz_clear(base_p);
	MG_curve_clear(&E);

	return 0;
}
