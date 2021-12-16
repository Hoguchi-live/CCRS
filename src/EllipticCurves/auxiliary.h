/// @file memory.h
#ifndef _AUXILIARY_H_
#define _AUXILIARY_H_

#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fq.h>

void fq_add_ui(fq_t, fq_t, ulong, const fq_ctx_t);
void fq_add_si(fq_t, fq_t, slong, const fq_ctx_t);
void fq_sub_ui(fq_t, fq_t, ulong, const fq_ctx_t);
void fq_sub_si(fq_t, fq_t, slong, const fq_ctx_t);

void fq_div_si(fq_t, fq_t, slong, const fq_ctx_t);
void fq_div_ui(fq_t, fq_t, ulong, const fq_ctx_t);
#endif

