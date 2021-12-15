/// @file binary_trees.h
#ifndef _BINARY_TREES_H_
#define _BINARY_TREES_H_

#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fq.h>
#include <flint/fq_poly.h>


/*********************************
  Structures
*********************************/
/**
  Representation of binary links with fq_poly as nodes.
*/
typedef struct fq_poly_blink_t {
	const fq_ctx_t *F;
	fq_poly_t data;
	//fq_poly_blink_t *parent;
	struct fq_poly_blink_t *left;
	struct fq_poly_blink_t *right;
} fq_poly_blink_t;

/**
  Representation of binary tree with fq_poly_blink_t nodes.
*/
typedef struct fq_poly_btree_t {
	const fq_ctx_t *F;
	fq_poly_blink_t *head;
} fq_poly_btree_t;


/*********************************
  Functions
*********************************/
void fq_poly_blink_init(fq_poly_blink_t *, const fq_ctx_t *);
void fq_poly_blink_set(fq_poly_blink_t *, fq_poly_t);
int fq_poly_blink_set_right(fq_poly_blink_t *, fq_poly_blink_t *);
int fq_poly_blink_set_left(fq_poly_blink_t *, fq_poly_blink_t *);
void fq_poly_blink_clear(fq_poly_blink_t *);

void fq_poly_btree_init(fq_poly_btree_t *, const fq_ctx_t *);
void fq_poly_btree_set(fq_poly_btree_t *, fq_poly_blink_t *);
void fq_poly_btree_clear(fq_poly_btree_t *);

#endif
