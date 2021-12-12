#ifndef _ARITHMETIC_H_
#define _ARITHMETIC_H_

#include "string.h"
#include "stdbool.h"

#include "models.h"


// Elliptic curves
void SW_j_invariant(fq_t *, SW_curve *);
void MG_j_invariant(fq_t *, MG_curve *);

// Points on elliptic curves
void SW_point_valid(bool *, SW_point *);

#endif
