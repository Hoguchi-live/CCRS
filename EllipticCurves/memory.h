#ifndef _MEMORY_H_
#define _MEMORY_H_

#include <string.h>

#include "models.h"

#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fq.h>

void SW_init(SW_curve *, const fq_ctx_t);
void SW_set(SW_curve *, const fq_ctx_t *, const fmpz_t, const fmpz_t);
void SW_clear(SW_curve *);

#endif
