/// @file multieval.c
#include "multieval.h"

/**
  Build a remainder tree from a fq_t roots array.
  This is supposed to be wrapped by the remainderTree function.
  The offset is used to slice the array (represent the position of element 0).
  This function is not responsible for the handling of the memory associated with the bcells.
*/
void remainderCell(fq_poly_bcell_t *rop, fq_t **roots, uint len, uint offset, fq_ctx_t *F) {

	uint new_offset;
	fq_t tmp;
	fq_poly_t p;

	fq_init(tmp, *F);
	fq_poly_init(p, *F);

	if(len == 0) {

		//// Set rop to a unique cell containing the 1 polynomial
		fq_set_ui(tmp, 1, *F);
		fq_poly_set_fq(p, tmp, *F);
		fq_poly_bcell_set(rop, p);
	}
	else if(len == 1) {

		//// Set rop to a unique cell containing the X-root polynomial
		fq_neg(tmp, (*roots)[offset], *F);
		fq_poly_set_coeff(p, 0, tmp, *F);

		fq_set_ui(tmp, 1, *F);
		fq_poly_set_coeff(p, 1, tmp, *F);

		fq_poly_bcell_set(rop, p);
	}
	else {

		//// Recursively construct the tree
		fq_poly_bcell_t *left, *right;

		left = malloc(sizeof(fq_poly_bcell_t));
		right = malloc(sizeof(fq_poly_bcell_t));

		fq_poly_bcell_init(left, F);
		fq_poly_bcell_init(right, F);

		remainderCell(left, roots, len/2, offset, F);
		remainderCell(right, roots, len - len/2, offset + len/2, F);

		//// Compute product of child polynomials
		fq_poly_mul(p, left->data, right->data, *F);

		fq_poly_bcell_set_(rop, left, right, p);
	}

	fq_clear(tmp, *F);
	fq_poly_clear(p, *F);
}

/**
  See remainderCell.
  T should be initialized.
*/
void remainderTree(fq_poly_btree_t *T, fq_t **roots, uint len, fq_ctx_t *F) {

	remainderCell(T->head, roots, len, 0, F);
}

void fq_poly_multieval(fq_t ** rop, fq_t ** op, fq_poly_t P, uint len, fq_ctx_t *F) {

	fq_poly_btree_t T;

	fq_poly_btree_init(&T, F);

	remainderTree(&T, op, len, F);

	fq_poly_btree_clear(&T);

}

