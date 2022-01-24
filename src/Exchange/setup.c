// @file setup.c
#include "setup.h"

/******************************
  	   l-prime
******************************/
void lprime_init(lprime_t *op){

	fmpz_init(op->l);
}

lprime_t *lprime_init_(){

	lprime_t *rop = malloc(sizeof(lprime_t));
	return rop;
}

void lprime_set(lprime_t *op, fmpz_t l, uint type, uint lbound, uint hbound, uint lorder, uint horder, uint r, uint bkw){

	fmpz_set(op->l, l);

	op->type = type;
	op->lbound = lbound;
	op->hbound = hbound;
	op->lorder = lorder;
	op->horder = horder;
	op->r = r;
	op->bkw = bkw;
}

void lprime_clear(lprime_t *op) {

	fmpz_clear(op->l);
}

/******************************
  	   config
******************************/
/**
  Returns a pointer to a cfg_t config context.
*/
cfg_t cgf_init_set() {

	//// BASE CURVE PARAMETERS
	char base_p_str[] = BASE_p;

	fmpz_t base_p;
	fmpz_init(base_p);
	fmpz_set_str(base_p, base_p_str, 0);

	fq_ctx_t F;
	char *Fgen = "g";
	fq_ctx_init(F, base_p, 1, Fgen);

	//// GLOBAL PROTOCOL PARAMETERS
	uint l_PRIMES_int[NB_PRIMES] = {3, 5, 7};
	uint l_PRIMES_LBOUNDS[NB_PRIMES] = {1000, 1000, 1000};
	uint l_PRIMES_HBOUNDS[NB_PRIMES] = {1000, 1000, 1000};
	uint l_PRIMES_LORDER[NB_PRIMES] = {};
	uint l_PRIMES_HORDER[NB_PRIMES] = {};
	uint l_PRIMES_R[NB_PRIMES] = {1, 1, 1};
	uint l_PRIMES_BKW[NB_PRIMES] = {1, 1, 1};

	//// Alloc config struct
	cfg_t *cfg = malloc(sizeof(cfg_t));

	//// Alloc lprimes array
	lprime_t **lprimes;
	lprimes = malloc(sizeof(lprime_t *) * NB_PRIMES);

	cfg->F = F;

	fq_init(cfg->A, *F);
	fq_init(cfg->B, *F);

	fq_set(cfg->A, A, *F);
	fq_set(cfg->B, B, *F);
	cfg->nb_primes = nb_primes;

	//// Create and set l-primes according to global variables
	lprime_t *lp;
	fmpz_t l_fmpz;
	uint type, lbound, hbound, lorder, horder, r, bkw;

	fmpz_init(l_fmpz);
	for(int i=0; i < nb_primes; i++) {

		uint l = l_PRIMES_int[i];
		lp = lprime_init_();
		fmpz_set_ui(l_fmpz, l);

		switch(l) {
			case 3: case 5: case 7:
				type = 1;	//radical isogeny
				lbound = l_PRIMES_LBOUNDS[i];
				hbound = l_PRIMES_HBOUNDS[i];
				lorder = l_PRIMES_LORDER[i];
				horder = l_PRIMES_HORDER[i];
				r = l_PRIMES_R[i];
				bkw = l_PRIMES_BKW[i];

				break;
			default:
				break;
		}
		lprime_set(lp, l_fmpz, type, lbound, hbound, lorder, horder, r, bkw);

		cfg->lprimes[i] = lp;

	}

	fmpz_clear(l_fmpz);
	return cfg;
}

void cfg_clear(cfg_t *op) {

	free(op->lprimes);

	fq_clear(op->A, *(op->F));
	fq_clear(op->B, *(op->F));
	fq_ctx_clear(*(op->F));
}

