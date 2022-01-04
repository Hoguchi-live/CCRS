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
void SW_curve_init(SW_curve_t *, const fq_ctx_t *);
void SW_curve_set(SW_curve_t *, const fq_ctx_t *, const fq_t, const fq_t);
void SW_curve_set_si(SW_curve_t *, const fq_ctx_t *, const slong, const slong);
void SW_curve_set_ui(SW_curve_t *, const fq_ctx_t *, const ulong, const ulong);
int SW_curve_set_str(SW_curve_t *, const fq_ctx_t *, const char *, const char *, int);
void SW_curve_clear(SW_curve_t *);


/*********************************************
   Short Weierstrass points memory management
*********************************************/
void SW_point_init(SW_point_t *, SW_curve_t *);
void SW_point_set(SW_point_t *, const fq_t, const fq_t, const fq_t, SW_curve_t *);
void SW_point_set_si(SW_point_t *, const slong, const slong, const slong, SW_curve_t *);
void SW_point_set_ui(SW_point_t *, const ulong, const ulong, const ulong, SW_curve_t *);
void SW_point_clear(SW_point_t *);


/**************************************
   Montgomery curves memory management
**************************************/
void MG_curve_init(MG_curve_t *, const fq_ctx_t *);
void MG_curve_set(MG_curve_t *, const fq_ctx_t *, const fq_t, const fq_t);
void MG_curve_set_si(MG_curve_t *, const fq_ctx_t *, const slong, const slong);
int MG_curve_set_str(MG_curve_t *, const fq_ctx_t *, const char *, const char *, int);
void MG_curve_clear(MG_curve_t *);


/*********************************************
   Montgomery points memory management
*********************************************/
void MG_point_init(MG_point_t *, MG_curve_t *);
void MG_point_set(MG_point_t *, const fq_t, const fq_t, MG_curve_t *);
void MG_point_set_si(MG_point_t *, const slong, const slong, MG_curve_t *);
void MG_point_set_ui(MG_point_t *, const ulong, const ulong, MG_curve_t *);
void MG_point_clear(MG_point_t *);

#endif

