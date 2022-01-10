#include <stdio.h>
#include <stdlib.h>

#include "models.h"

/**
  Sets rop to the Tate-normal form of op relative to point P
  rop must be initialized.
TODO: Check torsion automatically.
*/
void MG_get_TN(TN_curve_t *rop, MG_curve_t op, MG_point_t P, fmpz_t l){

	// Case l == 3

	// Case l >= 4
	fq_t b2, b3, b4, c1, c2;

	// Translation step


}

