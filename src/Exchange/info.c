// @file info.c
#include "info.h"

void print_verbose_walk_rad(uint type, fmpz_t l, fmpz_t k, int ec) {

	printf("VERBOSE::apply_key:Taking ");
	fmpz_print(k);
	printf(" step(s) in ");
	fmpz_print(l);
	printf("-isogeny graph using ");

	if(type == 1) printf("radical isogeny algorithm ");
	else printf("sqrt-Velu algorithm ");

	printf("resulted in error code %d ", ec);

	if(ec) printf("(success)\n");
	else printf("(failure)\n");

}
