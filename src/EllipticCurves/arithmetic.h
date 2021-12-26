#ifndef _ARITHMETIC_H_
#define _ARITHMETIC_H_

#include "string.h"
#include "stdbool.h"

#include "memory.h"
#include "models.h"
#include "auxiliary.h"

#include "../Polynomials/roots.h"

// Elliptic curves
void SW_j_invariant(fq_t *, SW_curve *);
void MG_j_invariant(fq_t *, MG_curve *);

// Points on elliptic curves
void SW_point_isinfinity(bool *, SW_point *);
void MG_point_isinfinity(bool *, MG_point *);
int SW_point_isvalid(bool *, SW_point *);
int MG_point_isvalid(bool *, MG_point *);
void MG_point_isinfty(bool *, MG_point *);
void MG_point_normalize(MG_point *);
void SW_point_rand_ninfty(SW_point *);
void MG_point_rand_ninfty(MG_point *);

// Montgomery curve arithmetic
void MG_xADD(MG_point *, MG_point, MG_point, MG_point);
void MG_xDBL(MG_point *, MG_point);

// Montgomery ladder
void MG_ladder_rec(MG_point *, MG_point *, fmpz_t, MG_point, const fq_ctx_t *);
void MG_ladder(MG_point *x0, fmpz_t k, MG_point P);
void MG_ladder_iter(MG_point *, MG_point *, fmpz_t, MG_point, fq_ctx_t *);
void MG_ladder_iter_(MG_point *, fmpz_t, MG_point *);

#endif
