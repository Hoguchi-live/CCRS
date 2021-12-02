#include <stdio.h>
#include <stdlib.h>

#include "memory.h"
#include "models.h"

void SW_init(SW_curve *E, const fq_ctx_t F) {

	fq_init(E->a, F);
	fq_init(E->b, F);
}

void SW_set(SW_curve *E, const fq_ctx_t *F, const fmpz_t a, const fmpz_t b) {

	E->F = F;
	fq_set_fmpz(E->a, a, *F);
	fq_set_fmpz(E->b, b, *F);
}

void SW_free(SW_curve *E) {

	fq_clear(E->a, *(E->F));
	fq_clear(E->b, *(E->F));
}
