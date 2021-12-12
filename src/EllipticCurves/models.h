#ifndef _MODELS_H_
#define _MODELS_H_

#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fq.h>

typedef struct SW_curve{

	const fq_ctx_t *F;	// base field
	fq_t a, b;		// curve parameters
} SW_curve;

typedef struct SW_point {

	fq_t x, y, z;	// coordinates
	SW_curve *E;	// base curve
} SW_point;

typedef struct MG_curve{

	const fq_ctx_t *F;
	fq_t A, B;
} MG_curve;


#endif
