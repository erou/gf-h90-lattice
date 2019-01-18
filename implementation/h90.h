#ifndef H90_H
#define H90_H

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

static inline slong tensor_order(struct tensor A) {
  return fq_nmod_ctx_degree(A.L);
}
static inline slong tensor_ext_degree(struct tensor A) {
  return fq_nmod_ctx_degree(A.R);
}
static inline mp_limb_t tensor_p(struct tensor A) {
  return A.L->mod.n;
}
static inline mp_limb_t tensor_pinv(struct tensor A) {
  return A.L->mod.ninv;
}

/*
 * Compute the polynomial Σ a_i X^i of degree deg(A.R) with
 *
 *   a_{k-1} = a
 *   a_i = a_{i+1}^p + A.R[i+1] a
 *
 * Can use to reconstruct a H90 solution from its (k-1)-th coefficient.
 */
void lift_h90(fq_nmod_poly_t res,
	      const fq_nmod_t a, const struct tensor A);

/*
 * Find a solution res to the hilbert 90 equation
 * 
 *   (σ⊗1)(res) = ζ · res
 *
 * in the algebra L⊗R, where ζ is the canonical root of R
 */
void solve_h90(fq_nmod_poly_t res, flint_rand_t state,
	       const struct tensor A);
/*
 * Tests if x∈A is a soultion of Hilbert 90 of order ord(A).
 */
int is_h90(const fq_nmod_poly_t x, const struct tensor A);

#endif
