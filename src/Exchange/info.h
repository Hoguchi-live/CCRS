#ifndef _info_H_
#define _info_H_

#include <stdio.h>
#include <stdlib.h>

#include "../../src/EllipticCurves/models.h"
#include "../../src/Exchange/keygen.h"
#include "../../src/Exchange/setup.h"

#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fq.h>

void print_verbose_walk_rad(uint, fmpz_t, fmpz_t, int);

#endif
