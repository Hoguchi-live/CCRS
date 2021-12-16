#ifndef _ARITHMETIC_H_
#define _ARITHMETIC_H_

#include "string.h"
#include "stdbool.h"

#include "models.h"
#include "auxiliary.h"

// Elliptic curves
void SW_j_invariant(fq_t *, SW_curve *);
void MG_j_invariant(fq_t *, MG_curve *);

// Points on elliptic curves
void SW_point_valid(bool *, SW_point *);

// Montgomery curve arithmetic
void MG_xADD(MG_point *, MG_point, MG_point, MG_point);
void MG_xDBL(MG_point *, MG_point);

#endif

