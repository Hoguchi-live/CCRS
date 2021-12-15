/// @file binary_trees.c
#include "binary_trees.h"

/*******************************
  blink
*******************************/
/**
  Initializes b for use, with context F, and setting its childs to NULL.
  A corresponding call to fq_poly_blink_clear() must be made after finishing with the fq_poly_blink_t to free the memory used by the link.
*/
void fq_poly_blink_init(fq_poly_blink_t *b, const fq_ctx_t *F) {

	b->F = F;
	fq_poly_init(b->data, *F);
	//b->parent = NULL;
	b->right = NULL;
	b->left = NULL;
}

/**
 Sets b to blink with data p.
*/
void fq_poly_blink_set(fq_poly_blink_t *b, fq_poly_t p) {

	fq_poly_set(b->data, p, *(b->F));
}

/**
  Sets right child of b1 to b2.
  Returns -1 if b1 already has a right child.
*/
int fq_poly_blink_set_right(fq_poly_blink_t *b1, fq_poly_blink_t *b2) {

	if(b1->right  != NULL) return -1;
	b1->right = b2;

	return 0;
}

/**
  Sets left child of b1 to b2.
  Returns -1 if b1 already has a left child.
*/
int fq_poly_blink_set_left(fq_poly_blink_t *b1, fq_poly_blink_t *b2) {

	if(b1->left  != NULL) return -1;
	b1->left = b2;

	return 0;
}

/**
  Recursively clears the given blink, releasing any memory used. It must be reinitialised in order to be used again.
*/
void fq_poly_blink_clear(fq_poly_blink_t *b) {

	printf("Freeing blink!\n");
	if(b == NULL) return;

	fq_poly_clear(b->data, *(b->F));

	if(b->right != NULL) fq_poly_blink_clear(b->right);
	if(b->left != NULL) fq_poly_blink_clear(b->left);
}

/*******************************
  btree
*******************************/
/**
  Initializes t for use, with context F.
  A corresponding call to fq_poly_btree_clear() must be made after finishing with the fq_poly_btree_t to free the memory used by the tree.
*/
void fq_poly_btree_init(fq_poly_btree_t *t, const fq_ctx_t *F) {

	t->F = F;
	t->head = NULL;
}

/**
 Sets t to tree with head blink b.
TODO: check ctx
*/
void fq_poly_btree_set(fq_poly_btree_t *t, fq_poly_blink_t *b) {

	t->head = b;
}

/**
  Recursively clears the given tree and its blinks, releasing any memory used. It must be reinitialised in order to be used again.
*/
void fq_poly_btree_clear(fq_poly_btree_t *t) {

	fq_poly_blink_clear(t->head);
}
