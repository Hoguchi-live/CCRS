/// @file memory.h
#ifndef _MEMORY_H_
#define _MEMORY_H_

#include <string.h>

#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fq.h>

#include "models.h"


/*********************************************
   Short Weierstrass curves memory management
*********************************************/
void SW_curve_init(SW_curve *, const fq_ctx_t *);
void SW_curve_set(SW_curve *, const fq_ctx_t *, const fq_t, const fq_t);
void SW_curve_set_si(SW_curve *, const fq_ctx_t *, const slong, const slong);
void SW_curve_set_ui(SW_curve *, const fq_ctx_t *, const ulong, const ulong);
void SW_curve_clear(SW_curve *);


/*********************************************
   Short Weierstrass points memory management
*********************************************/
void SW_point_init(SW_point *, const fq_ctx_t *);
void SW_point_set(SW_point *, const fq_t, const fq_t, const fq_t, SW_curve *);
void SW_point_set_si(SW_point *, const ulong, const ulong, const ulong, SW_curve *);
void SW_point_set_ui(SW_point *, const slong, const slong, const slong, SW_curve *);
void SW_point_clear(SW_point *);


/**************************************
   Montgomery curves memory management
**************************************/
void MG_curve_init(MG_curve *, const fq_ctx_t *);
void MG_curve_set(MG_curve *, const fq_ctx_t *, const fq_t, const fq_t);
void MG_curve_set_si(MG_curve *, const fq_ctx_t *, const ulong, const ulong);
void MG_curve_clear(MG_curve *);

#endif

