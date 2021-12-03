#ifndef _MEMORY_H_
#define _MEMORY_H_

#include <string.h>

#include "models.h"

#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fq.h>

void SW_curve_init(SW_curve *, const fq_ctx_t *);
void SW_curve_set(SW_curve *, const fq_ctx_t *, const fq_t, const fq_t);
void SW_curve_set_si(SW_curve *, const fq_ctx_t *, const ulong, const ulong);
void SW_curve_clear(SW_curve *);

void SW_point_init(SW_point *, const fq_ctx_t *);
void SW_point_set(SW_point *, const fq_t, const fq_t, const fq_t, SW_curve *);
void SW_point_set_si(SW_point *, const ulong, const ulong, const ulong, SW_curve *);
void SW_point_clear(SW_point *);

#endif

