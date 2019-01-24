#ifndef TENSOR_H
#define TENSOR_H

#include <flint/fq_nmod_poly.h>

/*
 * A struct representing an algebra L⊗R, with the convention that
 * deg(R) = ord(p mod deg(L)).
 *
 * Note that if R is trivial (i.e. if deg(L) divides p-1), the
 * defining polynomial must be a degree 1 polynomial X-ζ, with ζ a
 * deg(L)-th root of unity.
 */
struct tensor {
  fq_nmod_ctx_t L;
  fq_nmod_ctx_t R;
};

static inline slong tensor_degree(struct tensor A) {
  return fq_nmod_ctx_degree(A.L);
}
static inline slong tensor_level(struct tensor A) {
  return fq_nmod_ctx_degree(A.R);
}
static inline mp_limb_t tensor_p(struct tensor A) {
  return A.L->mod.n;
}
static inline mp_limb_t tensor_pinv(struct tensor A) {
  return A.L->mod.ninv;
}

/*
 * Rewrites a, seen as an element of (k[x]/from)[y]/to to an element
 * of (k[y]/to)[x]/from.
 *
 * Requires res to be a zero polynomial.
 */
void _transpose(fq_nmod_poly_t res, const fq_nmod_poly_t a,
		const fq_nmod_ctx_t ctx_from, const fq_nmod_ctx_t ctx_to);

/*
 * Computes a*b
 */
void tensor_mul(fq_nmod_poly_t res,
		const fq_nmod_poly_t a, const fq_nmod_poly_t b,
		struct tensor A);

#endif
