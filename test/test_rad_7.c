#include <stdio.h>
#include <stdlib.h>

#include "../src/EllipticCurves/models.h"
#include "../src/EllipticCurves/memory.h"
#include "../src/EllipticCurves/arithmetic.h"
#include "../src/EllipticCurves/pretty_print.h"

#include "../src/Polynomials/binary_trees.h"
#include "../src/Polynomials/roots.h"

#include "../src/Isogeny/radical.h"

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
	MG_curve_t E;
	MG_curve_init(&E, &F);
	MG_curve_set_str(&E, &F, BASE_A, BASE_B, 10);
	/**********************************
		     Cardinals
	**********************************/
	// Card over extensions
	fmpz_t card,card_quad, r, r_quad;

	fmpz_init(card);
	fmpz_init(card_quad);
	fmpz_init(r);
	fmpz_init(r_quad);

	fmpz_set_ui(r, 1);
	fmpz_set_ui(r_quad, 2);

	MG_curve_card_ext(card, &E, r);
	MG_curve_card_ext(card_quad, &E, r_quad);

	/**********************************
		     Torsion points
	**********************************/
	//Random torsion point over Fp
	MG_point_t P;
	fmpz_t l;

	MG_point_init(&P, &E);
	fmpz_init(l);

	fmpz_set_ui(l, 7);

	int ret = MG_curve_rand_torsion(&P, l, r, card);
	printf("\nTorsion returned %d\n", ret);
	MG_point_normalize(&P);
	printf("\n");
	/**********************************
		MG -> TN conversion
	**********************************/
	fq_t j_inv;
	TN_curve_t E_TN;

	fq_init(j_inv, F);
	TN_curve_init(&E_TN, &F);

	MG_get_TN(&E_TN, &E, &P, l);
	TN_j_invariant(&j_inv, &E_TN);
	printf("\nMG--->TN gave j_invariant: \n");
	fq_print_pretty(j_inv, F);
	/**********************************
	 	Isogeny forward
	**********************************/
	fmpz_t k;
	TN_curve_t E_for;

	fmpz_init_set_ui(k, 1);
	TN_curve_init(&E_for, &F);

	radical_isogeny_7(&E_for, &E_TN, k);
	TN_j_invariant(&j_inv, &E_for);
	printf("\n\nTN--radical-->TN: \n");
	TN_curve_print(&E_for);
	printf("\n\nTN--radical-->TN gave j_invariant: \n");
	fq_print_pretty(j_inv, F);
	printf("\n\n");
	/**********************************
		Conversion to MG
	**********************************/
	MG_curve_t E_for_MG;
	MG_curve_init(&E_for_MG, &F);
	ret = TN_get_MG(&E_for_MG, &E_for);
	printf("\nSuccessfully converted TN->MG: %d\n", ret);
	MG_curve_print(&E_for_MG);

	MG_j_invariant(&j_inv, &E_for_MG);
	printf("\nMG target curve gave j_invariant: \n");
	fq_print_pretty(j_inv, F);
	printf("\n\n");
	///**********************************
	//	MG to TN again
	//**********************************/
	//Random torsion point over Fp^2
	MG_point_t Q, Q_test;

	MG_point_init(&Q, &E_for_MG);
	MG_point_init(&Q_test, &E_for_MG);

	int ret_quad = MG_curve_rand_torsion_(&Q, l, r, card_quad);
	printf("\nTorsion returned %d\n", ret_quad);
	MG_point_normalize(&Q);

	MG_ladder_iter_(&Q_test, l, &Q);
	MG_point_normalize(&Q_test);
	MG_point_print(&Q_test);
	/**********************************
		MG -> TN conversion
	**********************************/
	fq_t j_inv_quad;
	TN_curve_t E_for_TN_tors;

	fq_init(j_inv_quad, F);
	TN_curve_init(&E_for_TN_tors, &F);

	MG_get_TN(&E_for_TN_tors, &E_for_MG, &Q, l);
	TN_j_invariant(&j_inv_quad, &E_for_TN_tors);
	printf("\nMG--->TN gave j_invariant: \n");
	fq_print_pretty(j_inv_quad, F);
	/**********************************
	 	Isogeny backward
	**********************************/
	TN_curve_t E_back;

	TN_curve_init(&E_back, &F);

	radical_isogeny_7(&E_back, &E_for_TN_tors, k);
	TN_j_invariant(&j_inv_quad, &E_back);
	printf("\n\nTN--radical-->TN: \n");
	TN_curve_print(&E_back);
	printf("\n\nTN--radical-->TN gave j_invariant: \n");
	fq_print_pretty(j_inv_quad, F);
	printf("\n\n");

	/**********************************
		Clear Memory
	**********************************/
	fmpz_clear(r);
	fmpz_clear(l);
	fmpz_clear(card);
	MG_point_clear(&P);
	fmpz_clear(k);
	TN_curve_clear(&E_for);

	fq_ctx_clear(F);
	fmpz_clear(base_p);
	MG_curve_clear(&E);

	return 0;
}

