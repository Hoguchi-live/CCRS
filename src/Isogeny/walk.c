// @file walk.c
#include "walk.h"

/**
  Take k steps in the l-isogeny graph using radical isogeny.
TODO: harcode cardinals, update MG_curve_set_ to MG_curve_set. Update with card_ext.
	MG_get_TN should return an int error code.
	radical_isogeny should return an int error code.
**/
int walk_rad(MG_curve_t *rop, MG_curve_t *op, fmpz_t l, fmpz_t k) {

	int ec = 1;

	// Nothing to do
	if(fmpz_equal_ui(k, 0)) {
		MG_curve_set_(rop, op);
		return ec;
	}

	//// Init variables
	fmpz_t k_local;
	MG_point_t P;
	TN_curve_t E_TN_tmp1, E_TN_tmp2;
	fmpz_t card, r;

	fmpz_init_set(k_local, k);
	MG_point_init(&P, op);
	TN_curve_init(&E_TN_tmp1, l, op->F);
	TN_curve_init(&E_TN_tmp2, l, op->F);
	fmpz_init(r);

	//// Direction of the walk
	if(fmpz_cmp_ui(k, 0) >= 0) {
		// case k>0
		fmpz_set_ui(r, 1);
		MG_curve_card_ext(card, op, r);
		ec = MG_curve_rand_torsion(&P, l, card);
		printf("walk::MG->TN result: %d\n", ec);
	}
	else {
		// case k<0
		fmpz_set_ui(r, 2);
		fmpz_neg(k_local, k_local);
		MG_curve_card_ext(card, op, r);
		ec = MG_curve_rand_torsion_(&P, l, card);
		printf("walk::MG->TN result: %d\n", ec);
	}

	//// Transform op in Tate-normal form
	MG_get_TN(&E_TN_tmp1, op, &P, l);

	/// TEST
	fq_t j;
	fq_init(j, *(op->F));
	printf("walk::op in MG j-inv: ");
	MG_j_invariant(&j, op);
	fq_print_pretty(j, *(op->F));
	printf("\nwalk::op in TN j-inv: ");
	TN_j_invariant(&j, &E_TN_tmp1);
	fq_print_pretty(j, *(op->F));
	printf("\n");
	/// TEST

	//// Walk
	if(fmpz_equal_ui(l, 3)) radical_isogeny_3(&E_TN_tmp2, &E_TN_tmp1, k_local);
	else if(fmpz_equal_ui(l, 5)) radical_isogeny_5(&E_TN_tmp2, &E_TN_tmp1, k_local);
	else if(fmpz_equal_ui(l, 7)) radical_isogeny_7(&E_TN_tmp2, &E_TN_tmp1, k_local);
	else return 0;

	/// TEST
	printf("walk::rop in TN j-inv: ");
	TN_j_invariant(&j, &E_TN_tmp2);
	fq_print_pretty(j, *(op->F));
	printf("\n");
	/// TEST

	//// Transform result back into Mongomery form
	ec = TN_get_MG(rop, &E_TN_tmp2);
	printf("walk::TN get MG success: %d\n", ec);

	/// TEST
	printf("walk::rop in MG j-inv: ");
	MG_j_invariant(&j, rop);
	fq_print_pretty(j, *(op->F));
	printf("\n\n");
	/// TEST

	//// Clear
	fmpz_clear(k_local);
	MG_point_clear(&P);
	TN_curve_clear(&E_TN_tmp1);
	TN_curve_clear(&E_TN_tmp2);
	fmpz_clear(r);

	return ec;
}

/**
  Take k steps in the l-isogeny graph using the sqrt-velu algorithm.
TODO: harcode cardinals, update MG_curve_set_ to MG_curve_set. Update with card_ext.
**/
int walk_velu(MG_curve_t *rop, MG_curve_t *op, uint l, uint r, uint k) {

	int ec = 1;

	// Nothing to do
	if(k == 0) {
		MG_curve_set_(rop, op);
		return ec;
	}

	//// Init variables
	uint b, bprime, lenK;

	//// Compute constants
	_init_lengths(&b, &bprime, &lenK, l);

	//// TEST
	printf("Found b = %d, b' = %d\n", b, bprime);
	//// TEST

	//// Compute I, J and K


	return ec;
}

