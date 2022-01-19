#include "velu.h"

void KPS(MG_point_t P, int l) {
	//TODO: comment
	//TODO: I,J,K as pointers, define, initialize and compute their length b,bprime and lenK in another function
	
	MG_curve_t *E = P.E;

	//computing b and bprime
	int b, bprime;
	bprime = l-1;
	b = 1;
	while ((b+1)*(b+1) <= bprime) {
		b++;
	}
	b = b/2;
	bprime = bprime/(4*b);

	//computing J = {(2j+1)*P for j = 1, ..., b-1}
	//if l>17 then bprime >= b >= 2 therefore J has at least two elements
	MG_point_t J[b];
	MG_point_t P2; //P2 = 2*P
	MG_point_init(&P2, E);
	MG_xDBL(&P2, P);
	J[0] = P;
	MG_point_init(&J[1], E);
	MG_xADD(&J[1], P, P2, P); //J[1] = 3*P
	for (int j=2; j<b; j++) {
		MG_point_init(&J[j], E);
		MG_xADD(&J[j], J[j-1], P2, J[j-2]);
	}

	//computing I = {2b(2i+1)*P for i = 1, ..., bprime-1}
	MG_point_t I[bprime];
	MG_point_t P4, P4b; 
	MG_point_init(&P4, E);
	MG_xDBL(&P4, P2); //P4 = 4*P
	
	MG_point_init(&I[0], E); //I[0] = 2b*P
	if (b%2) {
		MG_xADD(&I[0], J[b/2], J[b-(b/2)], P2);
	}
	else {
		MG_xADD(&I[0], J[b/2], J[b-(b/2)], P4);
	}

	MG_point_init(&P4b, E);
	MG_xDBL(&P4b, I[0]); // P4b = 4b*P
	MG_point_init(&I[1], E);
	MG_xADD(&I[1], P4b, I[0], I[0]); // I[1] = 6b*P = 4b*P + 2b*P

	for (int i=2; i<bprime; i++) {
		MG_point_init(&I[i], E);
		MG_xADD(&I[i], I[i-1], P4b, I[i-2]);
	}

	//computing K = {i*P for i = 4*b*bprime+1, ..., l-4, l-2}
	int lenK = (l-1-4*b*bprime)/2;
	MG_point_t K[lenK];
	if (lenK>0) {
		MG_point_init(&K[lenK-1], E);
		K[lenK-1] = P2; // (l-2)*P = -2*P
		if (lenK>1) {
			MG_point_init(&K[lenK-2], E);
			K[lenK-2] = P4; // (l-4)*P = -4*P
		}
	}
	for (int i = lenK-3; i>0; i--) {
		MG_point_init(&K[i], E);
		MG_xADD(&K[i], K[i+1], P2, K[i+2]);
	}

	// TODO: Memory management
}

void xISOG(fq_t *A2, MG_point_t P, int l, MG_point_t I[], MG_point_t J[], MG_point_t K[], int b, int bprime) {

	const fq_ctx_t *F;
	F = (P.E)->F;

	fq_poly_t h, E0, E1, R0, R1, M0, M1;
	fq_poly_init(h, *F); //TODO: initialize with length brpime = #I
	fq_poly_one(h, *F);

	fq_poly_t f, c;
	fq_poly_init(f, *F);
	fq_poly_init(c, *F);
	for (int i=0; i<bprime; i++) {
		fq_poly_gen(f, *F);
		fq_poly_set_fq(c, I[i].X, *F);
		fq_poly_sub(f, f, c, *F);
		fq_poly_mul(h, h, f, *F);
	}
	fq_poly_clear(f, *F);
	fq_poly_clear(c, *F);
}

void _F0(fq_poly_t *rop, MG_point_t P, fq_ctx_t ctx) {
	fq_poly_zero(*rop, ctx);
	fq_t tmp;
	fq_init(tmp, ctx);
	fq_one(tmp, ctx);
	fq_poly_set_coeff(*rop, 2, tmp, ctx); //coeff X^2 = 1
	fq_mul_si(tmp, P.X, -2, ctx);
	fq_poly_set_coeff(*rop, 1, tmp, ctx); //coeff X^1 = -2*x
	fq_sqr(tmp, P.X, ctx);
	fq_poly_set_coeff(*rop, 0, tmp, ctx); //coeff const = x^2
	fq_clear(tmp, ctx);
}

void _F1(fq_poly_t *rop, MG_point_t P, fq_ctx_t ctx) {
}

void _F2(fq_poly_t *rop, MG_point_t P, fq_ctx_t ctx) {
}