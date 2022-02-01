// @file info.c
#include "info.h"

void print_verbose_walk_rad(uint type, fmpz_t l, fmpz_t k, int ec, int msec) {

	printf("VERBOSE::apply_key:Taking ");
	fmpz_print(k);
	printf(" step(s) in the ");
	fmpz_print(l);
	printf("-isogeny graph using ");

	if(type == 1) printf("[radical isogeny] algorithm ");
	else printf("[sqrt-Velu] algorithm ");

	printf("resulted in error code %d ", ec);

	if(ec) printf("(success) ");
	else printf("(failure) ");

	if(msec > 1000) printf("in <%ds> \n", msec/1000);
	else printf("in <%dms> \n", msec);

}

void print_verbose_walk_total_time(int total_time){


	printf("VERBOSE::apply_key:Total user-time: %ds\n", total_time/1000);
}
