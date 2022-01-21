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
//#define BASE_A	"6994408703168739098746713699616203766621720654593622191587805269977974214071997593327742781251142902555902895849326720402577188961785807254661257646688533"
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
	fmpz_t card, r;

	fmpz_init(card);
	fmpz_init(r);

	fmpz_set_ui(r, 1);

	MG_curve_card_ext(card, &E, r);

	/**********************************
		     Torsion
	**********************************/

	fmpz_t l;
	MG_point_t Q, Q_test;

	fmpz_init_set_ui(l, 5);
	MG_point_init(&Q, &E);
	MG_point_init(&Q_test, &E);

	int ret = MG_curve_rand_torsion(&Q, l, r, card);
	printf("Res: %d\n", ret);

	return 0;
}

