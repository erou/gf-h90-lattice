#ifndef MINPOLY_H
#define MINPOLY_H

#include <flint/fq_nmod.h>
#include <flint/nmod_poly_mat.h>

void minpoly(nmod_poly_t result,
	     const fq_nmod_t f, const fq_nmod_ctx_t ctx);

#endif
