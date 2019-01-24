#include <flint/fq_nmod_poly.h>
#include <flint/fq_nmod_poly_factor.h>

/*
 * Returns an nth root of a.
 */
void nth_root(fq_nmod_t res, const fq_nmod_t a, ulong n, const fq_nmod_ctx_t ctx) {
  fq_nmod_poly_t P;
  flint_rand_t state;
  flint_randinit(state);

  fq_nmod_poly_init(P, ctx);

  fq_nmod_poly_one(P, ctx);
  fq_nmod_poly_neg(P, P, ctx);
  fq_nmod_poly_shift_left(P, P, n, ctx);
  fq_nmod_poly_set_coeff(P, 0, a, ctx);
  fq_nmod_poly_neg(P, P, ctx);
  fq_nmod_poly_factor_equal_deg_prob(P, state, P, 1, ctx);
  fq_nmod_poly_get_coeff(res, P, 0, ctx);
  fq_nmod_neg(res, res, ctx);
    
  fq_nmod_poly_clear(P, ctx);
  flint_randclear(state);
}
