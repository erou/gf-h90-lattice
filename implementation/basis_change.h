#ifndef BASIS_CHANGE_H
#define BASIS_CHANGE_H

#include <flint/fq_nmod.h>

/*
 * Given a ∈ ctx_from, and given the image g of the canonical
 * generator of ctx_to in ctx_from, compute the image of a in ctx_to.
 *
 * ctx_to can be smaller than ctx_from, but not bigger. In this case,
 * deg(ctx_from)/deg(ctx_to) must not be divisible by the
 * characteristic, othwerise 0 is returned.
 */
void change_basis_inverse(fq_nmod_t res,
			  const fq_nmod_t a, const fq_nmod_t g,
			  const fq_nmod_ctx_t ctx_from, const fq_nmod_ctx_t ctx_to);

/*
 * Given a ∈ ctx_from, and given the image g of the canonical
 * generator of ctx_from in ctx_to, compute the image of a in ctx_to.
 *
 * I.e., compute a(g) mod ctx_to.
 *
 * ctx_from is unused and a NULL pointer can be passed in its place.
 */
void change_basis_direct(fq_nmod_t res,
			 const fq_nmod_t a, const fq_nmod_ctx_t ctx_from_unused,
			 const fq_nmod_t g, const fq_nmod_ctx_t ctx_to);

#endif
