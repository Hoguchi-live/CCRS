/// @file auxiliary.c
#include "auxiliary.h"

/**
  Sets rop to the sum of op and x, where x is an ulong considered as an element of F.
*/
void fq_add_ui(fq_t rop, fq_t op, ulong x, const fq_ctx_t F){

	fq_t xx;
	fq_init(xx, F);
	fq_set_ui(xx, x, F);

	fq_add(rop, op, xx, F);

	fq_clear(xx, F);
}


/**
  Sets rop to the sum of op and x, where x is a slong considered as an element of F.
*/
void fq_add_si(fq_t rop, fq_t op, slong x, const fq_ctx_t F){

	fq_t xx;
	fq_init(xx, F);
	fq_set_si(xx, x, F);

	fq_add(rop, op, xx, F);

	fq_clear(xx, F);
}

/**
  Sets rop to the difference of op and x, where x is an ulong considered as an element of F.
*/
void fq_sub_ui(fq_t rop, fq_t op, ulong x, const fq_ctx_t F){

	fq_t xx;
	fq_init(xx, F);
	fq_set_ui(xx, x, F);

	fq_sub(rop, op, xx, F);

	fq_clear(xx, F);
}


/**
  Sets rop to the sum of op and x, where x is a slong considered as an element of F.
*/
void fq_sub_si(fq_t rop, fq_t op, slong x, const fq_ctx_t F){

	fq_t xx;
	fq_init(xx, F);
	fq_set_si(xx, x, F);

	fq_sub(rop, op, xx, F);

	fq_clear(xx, F);
}

