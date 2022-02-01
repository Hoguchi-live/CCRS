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
	lprime_init(rop);
	return rop;
}

void lprime_set(lprime_t *op, fmpz_t l, uint type, uint lbound, uint hbound, uint r, uint bkw){

	fmpz_set(op->l, l);

	op->type = type;
	op->lbound = lbound;
	op->hbound = hbound;
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
  Everything is hardcoded here.
*/
cfg_t* cfg_init_set() {

	//// Alloc config struct
	cfg_t *cfg = malloc(sizeof(cfg_t));

	//// BASE FIELD
	fq_ctx_t *F;
	fmpz_t base_p;
	char *Fgen = "x";
	char base_p_str[] = BASE_p;

	fmpz_init(base_p);
	fmpz_set_str(base_p, base_p_str, 0);

	F = malloc(sizeof(fq_ctx_t));
	fq_ctx_init(*F, base_p, 1, Fgen);

	cfg->F = F;

	//// BASE CURVE PARAMETERS
	MG_curve_t *E;

	E = malloc(sizeof(MG_curve_t));

	MG_curve_init(E, F);
	MG_curve_set_str(E, F, BASE_A, BASE_B, 10);

	cfg->E = E;

	//// GLOBAL PROTOCOL PARAMETERS
	uint l_PRIMES_int[NB_PRIMES] = {3, 5, 7, 11, 13, 17, 103, 523, 821, 947, 1723,    19, 661,     1013, 1181,     31, 61, 1321,
					29, 71, 547,
					881,
					37, 1693};
	uint l_PRIMES_LBOUNDS[NB_PRIMES] = {1000, 1000, 1000, 100, 100, 100, 100, 100, 100, 100, 100, 10, 10, 10, 10,    10, 10, 10,
					5, 5, 5,
					5,
					5, 5};
	uint l_PRIMES_HBOUNDS[NB_PRIMES] = {1000, 1000, 1000, 100, 100, 100, 100, 100, 100, 100, 100, 10, 10, 10, 10,    10, 10, 10,
					5, 5, 5,
					5,
					5, 5};
	uint l_PRIMES_R[NB_PRIMES] =   {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,    3, 3,    4, 4,   5, 5, 5,
					7, 7, 7,
					8,
					9, 9};
	uint l_PRIMES_BKW[NB_PRIMES] = {1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1,    1, 0,    0, 0,   1, 1, 0,
					1, 1, 0,
					0,
					1, 0};

	//// Alloc lprimes array
	cfg->lprimes = (lprime_t *)malloc(sizeof(lprime_t) * NB_PRIMES);
	cfg->nb_primes = NB_PRIMES;

	//// Create and set l-primes accordingly
	fmpz_t l_fmpz;
	uint type, lbound, hbound, r, bkw;

	fmpz_init(l_fmpz);

	for(int i=0; i < NB_PRIMES; i++) {

		uint l = l_PRIMES_int[i];
		fmpz_set_ui(l_fmpz, l);

		switch(l) {
			case 3: case 5: case 7:
				type = 1;	//radical isogeny TODO: SHOULD USE ENUM
				lbound = l_PRIMES_LBOUNDS[i];
				hbound = l_PRIMES_HBOUNDS[i];
				r = l_PRIMES_R[i];
				bkw = l_PRIMES_BKW[i];

				break;
			default:
				type = 2;
				lbound = l_PRIMES_LBOUNDS[i];
				hbound = l_PRIMES_HBOUNDS[i];
				r = l_PRIMES_R[i];
				bkw = l_PRIMES_BKW[i];
				break;
		}
		lprime_init(&(cfg->lprimes)[i]);
		lprime_set(&(cfg->lprimes)[i], l_fmpz, type, lbound, hbound, r, bkw);
	}

	//// BASE FIELD EXTENSIONS
	//// Alloc fields array
	cfg->fields = (fq_ctx_t *)malloc(sizeof(fq_ctx_t) * MAX_EXTENSION_DEGREE);
	char gen[] = "x";

	//// Initialize extensions
	for(int i=1; i < MAX_EXTENSION_DEGREE + 1; i++) {

		fq_ctx_init( (cfg->fields)[i-1], base_p, i , gen);
	}

	//// RANDOM SEED FOR KEY GENERATION
	cfg->seed = 0;

	fmpz_clear(l_fmpz);
	fmpz_clear(base_p);
	return cfg;
}

void cfg_print(cfg_t *cfg) {

	printf("*** Config structure ***\n Finite field Fp^d \np = ");
	fmpz_print(fq_ctx_prime(*(cfg->F)));
	printf("\nd = %d", fq_ctx_degree(*(cfg->F)));
	printf("\n Base elliptic curve BY^2 = X^3 + AX^2 + X\nA = ");
	fq_print_pretty(cfg->E->A, *(cfg->F));
	printf("\nB = ");
	fq_print_pretty(cfg->E->B, *(cfg->F));
	printf("\n l-primes \nnb_primes = %d\n", cfg->nb_primes);
	printf("************************\n");
}

void cfg_clear(cfg_t *op) {

	//// Free l-primes
	for(int i = 0; i < NB_PRIMES; i++) lprime_clear( &(op->lprimes)[i] );
	free(op->lprimes);

	//// Free fields
	for(int i = 0; i < MAX_EXTENSION_DEGREE; i++) fq_ctx_clear( (op->fields)[i] );

	//// Free structures
	MG_curve_clear(op->E);
	fq_ctx_clear(*((fq_ctx_t*)(op->F)));
}
