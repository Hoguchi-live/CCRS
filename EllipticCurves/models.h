#ifndef _MODELS_H_
#define _MODELS_H_

#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fq.h>

typedef struct SW_curve{

	const fq_ctx_t *F;	// base field

	fq_t a, b;		// curve parameters
} SW_curve;
#endif
