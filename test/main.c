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
	MG_curve_t E;
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
	//SW_point_t P;
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
	//SW_point_t P;
	//SW_point_init(&P, &E);
	//SW_point_print(&P);
	//SW_point_clear(&P);

	/**********************************
		Montgomery Points
	**********************************/
	//fmpz_t k;
	//MG_point_t P, Q, R, res;

	//fmpz_init(k);
	//MG_point_init(&P, &E);
	//MG_point_init(&Q, &E);
	//MG_point_init(&R, &E);
	//MG_point_init(&res, &E);

	//fmpz_set_ui(k, 2);
	//MG_point_rand_ninfty(&P);

	// Rand Sage test
	//fmpz_t AA;
	//fmpz_set_str(AA, "5271924811707382594219474826841227070541619870303527507321729363247235083901432543590160620034604195244770376879479244051028554117448522106375394422417497", 10);
	//fq_set_fmpz(P.X, AA, F);
	//fq_set_ui(P.Z, 1, F);
	//fmpz_clear(AA);

	/***************************
		Test xADD
	***************************/
	// P
	//fmpz_t tmp;
	//fmpz_set_str(tmp, "1706202762055133895294293812437047635215141946710319269825765792863929768061371366694612669405467945104077853682860729195547221027289546642680421098590919", 10);
	//fq_set_fmpz(P.X, tmp, F);
	//fq_set_ui(P.Z, 1, F);
	//fmpz_clear(tmp);
	//// Q = infty
	//MG_point_set_ui(&Q, 1, 0, &E);
	//// Q = P
	//MG_point_set(&R, P.X, P.Z, &E);

	//MG_xADD(&res, P, Q, R);
	//MG_point_normalize(&res);
	//MG_point_print(&res);

	/***************************
		Test xDBL
	***************************/
	//MG_xDBL(&res, P);
	//MG_point_normalize(&res);
	//MG_point_print(&res);

	//fmpz_clear(k);
	//MG_point_clear(&P);
	//MG_point_clear(&Q);
	//MG_point_clear(&R);
	//MG_point_clear(&res);

	/**********************************
		     Torsion
	**********************************/

	// Card over extension
	fmpz_t card, r;
	fmpz_init(card);
	fmpz_init(r);

	fmpz_set_ui(r, 1);

	MG_curve_card_ext(card, &E, r);

	//Random torsion point
	MG_point_t P, Q;
	bool isinfty;
	fmpz_t l;

	MG_point_init(&P, &E);
	MG_point_init(&Q, &E);
	fmpz_init(l);

	fmpz_set_ui(l, 7);

	int ret = MG_curve_rand_torsion(&P, l, r, card);
	printf("\nTorsion returned %d\n", ret);
	MG_ladder_iter_(&Q, l, &P);
	MG_point_normalize(&P);
	MG_point_normalize(&Q);
	printf("Test torsion: \n");
	MG_point_print(&P);
	printf("\n");
	printf("\n");
	/**********************************
		MG -> TN
	**********************************/
	fq_t j_inv;
	TN_curve_t E_TN;
	MG_curve_t E_recover;

	fq_init(j_inv, F);
	TN_curve_init(&E_TN, &F);
	MG_curve_init(&E_recover, &F);

	// Print original j-invariant
	MG_j_invariant(&j_inv, &E);

	// MG to TN
	MG_get_TN(&E_TN, &E, &P, l);
	TN_curve_print(&E_TN);
	TN_j_invariant(&j_inv, &E_TN);
	printf("\nMG--->TN gave j_invariant: \n");
	fq_print_pretty(j_inv, F);

	//// TN to MG
	//int ret = TN_get_MG(&E_recover, &E_TN);
	//printf("\nSuccessfully converted TN->MG: %d\n", ret);
	//MG_curve_print(&E_recover);

	//// Print recovered j-invariant
	//MG_j_invariant(&j_inv, &E_recover);
	//printf("\n");
	//fq_print_pretty(j_inv, F);

	// clear curves
	MG_curve_clear(&E_recover);
	TN_curve_clear(&E_TN);
	fq_clear(j_inv, F);

	// clear torsion
	fmpz_clear(r);
	fmpz_clear(l);
	fmpz_clear(card);
	MG_point_clear(&P);

	/**********************************
		Clear Memory
	**********************************/
	fq_ctx_clear(F);
	fmpz_clear(base_p);
	MG_curve_clear(&E);

	return 0;
}

