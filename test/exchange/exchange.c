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

#include "../../src/Exchange/setup.h"
#include "../../src/Exchange/keygen.h"
#include "../../src/Exchange/dh.h"
#include "../../src/Exchange/info.h"

#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fq.h>
#include <flint/fq_poly.h>
#include <flint/fq_poly_factor.h>

int main() {

	cfg_t *cfg = cfg_init_set();

	const fq_ctx_t *F = cfg->F;
	MG_curve_t E = *(cfg->E);

	fq_t j_inv;
	fmpz_t k, l;
	MG_curve_t E_tmp1, E_tmp2;

	fq_init(j_inv, *F);
	fmpz_init_set_ui(l, 17);
	fmpz_init_set_si(k, 1);
	MG_curve_init(&E_tmp1, F);
	MG_curve_init(&E_tmp2, F);

	//// Config
	cfg_print(cfg);

	//// Keygen
	key__t *key = keygen_(cfg);
	key_print(key);

	//// Apply key to base curve
	int ret = apply_key(&E_tmp1, cfg->E, key, cfg);



	//MG_j_invariant(&j_inv, &E_tmp2);
	//printf("Found jinv: \n");
	//fq_print_pretty(j_inv, *F);

	cfg_clear(cfg);
	MG_curve_clear(&E_tmp1);
	MG_curve_clear(&E_tmp2);
}
