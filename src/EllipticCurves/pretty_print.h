#ifndef _PRETTY_PRINT_H_
#define _PRETTY_PRINT_H_

#include <string.h>

#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fq.h>

#include "models.h"

void SW_curve_print(SW_curve *);
void SW_point_print(SW_point *);
void MG_curve_print(MG_curve *);
void MG_point_print(MG_point *);

void SW_curve_print_full(SW_curve *);
void SW_point_print_full(SW_point *);
void MG_curve_print_full(MG_curve *);
void MG_point_print_full(MG_point *);

//compact

#endif
