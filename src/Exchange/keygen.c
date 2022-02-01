// @file arithmetic.c
#include "keygen.h"

/**
  Initialize key for use with keygen and cfg.
*/
void key_init(key__t *key, cfg_t *cfg) {

	key->nb_primes = cfg->nb_primes;

	key->lprimes = malloc(sizeof(lprime_t) * key->nb_primes);
	memcpy(key->lprimes, cfg->lprimes, sizeof(lprime_t) * key->nb_primes);

	key->steps = malloc(sizeof(fmpz_t) * key->nb_primes);
	for(int i = 0; i < key->nb_primes; i++) fmpz_init((key->steps)[i]);
}

key__t *key_init_(cfg_t *cfg) {

	key__t *key = malloc(sizeof(key__t));
	key_init(key, cfg);
	return key;
}


void keygen(key__t *key, cfg_t *cfg) {

	fmpz_t mod;
	fmpz_t steps;
	flint_rand_t state;

	fmpz_init(mod);
	fmpz_init(steps);
	srand(cfg->seed);
	flint_randinit(state); // Should depend on cfg->seed

	for(int i = 0; i < key->nb_primes; i++) {

		// direction (0 backward, 1 forward)
		if(rand() % 2) {
			fmpz_set_ui(mod, (key->lprimes)[i].hbound + 1);
			fmpz_randtest_mod(steps, state, mod);

			fmpz_set((key->steps)[i], steps);
		}
		else {
			fmpz_set_ui(mod, (key->lprimes)[i].hbound + 1);
			fmpz_randtest_mod(steps, state, mod);
			fmpz_neg(steps, steps);

			fmpz_set((key->steps)[i], steps);
		}
	}

	fmpz_clear(mod);
	fmpz_clear(steps);
	flint_randclear(state);
}

key__t *keygen_(cfg_t *cfg) {
	key__t *key = key_init_(cfg);
	keygen(key, cfg);
	return key;
}

void key_clear(key__t *key) {

	free(key->lprimes);
	free(key->steps);
}

void key_print(key__t *key){

	printf("*** Key structure ***\n");
	for(int i = 0; i < key->nb_primes; i++) {
		fmpz_print((key->lprimes)[i].l);
		printf(" ---> ");
		fmpz_print((key->steps)[i]);
		printf("\n");
	}
	printf("*********************\n");
}
