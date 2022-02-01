// @file setup.c
#include "dh.h"

/**
  Walk in the isogeny graph starting from op following the path given by the key.
  The curve rop is set to the endpoint of this walk.
  rop must be initialized.
  The l-primes are ordered according to their working extension.
  That way we just change the extension degree when needed.
  Returns 1 if successful and 0 if an error occured during a walk.
*/
int apply_key(MG_curve_t *rop, MG_curve_t *op, key__t *key, cfg_t *cfg) {

	int ec = 1;
	int r = 1;
	lprime_t *lp;
	fmpz_t *steps;
	MG_curve_t tmp1, tmp2;

	const fq_ctx_t *F = cfg->F;

	MG_curve_init(&tmp1, F);
	MG_curve_init(&tmp2, F);

	//// Init tmp1 at base curve op
	MG_curve_set_(&tmp1, op);

	for(int i = 0; i < key->nb_primes; i++) {

		lp = (key->lprimes) + i;
		steps = (key->steps) + i;

		if( lp->type == 1 ) ec = walk_rad(&tmp2, &tmp1, lp->l, *steps);
		else ec = walk_velu(&tmp2, &tmp1, lp->l, *steps);
		#ifdef VERBOSE
		print_verbose_walk_rad(lp->type, lp->l, *steps, ec);
		#endif
	}

	MG_curve_set_(rop, &tmp2);

	MG_curve_clear(&tmp1);
	MG_curve_clear(&tmp2);

	return ec;
}
