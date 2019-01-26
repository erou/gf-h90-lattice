#include <flint/fmpz.h>
#include "basis_change.h"
#include "minpoly.h"

/**
 * Computes the dual representation of the given polynomial {@code a} modulo {@code ctx}
 */
void monomial_to_dual(mp_limb_t *dual, const nmod_poly_t a, const fq_nmod_ctx_t ctx) {
  nmod_poly_t temp;
  nmod_poly_init(temp, ctx->modulus->mod.n);
  
  slong m = nmod_poly_degree(ctx->modulus);
  
  nmod_poly_derivative(temp, ctx->modulus);
  nmod_poly_mulmod(temp, temp, a, ctx->modulus);
  nmod_poly_reverse(temp, temp, m);
  nmod_poly_mullow(temp, temp, ctx->inv, m);
  
  for (slong i = 0; i < m; i++)
    dual[i] = nmod_poly_get_coeff_ui(temp, i);
  
  nmod_poly_clear(temp);
}

/**
 * Computes the representation of the given dual {@code dual}
 * in monomial basis modulo {@code modulus} 
 */
void dual_to_monomial(nmod_poly_t result, const mp_limb_t *dual, const nmod_poly_t modulus) {
  nmod_poly_t temp1;
  nmod_poly_t temp2;
  nmod_poly_init(temp1, modulus->mod.n);
  nmod_poly_init(temp2, modulus->mod.n);

  slong m = nmod_poly_degree(modulus);
  for (slong i = 0; i < m; i++)
    nmod_poly_set_coeff_ui(temp1, i, dual[i]);
  
  nmod_poly_reverse(temp2, modulus, m + 1);
  nmod_poly_mullow(temp1, temp1, temp2, m);
  
  nmod_poly_reverse(temp2, temp1, m);
  
  // compute 1 / modulus'
  // This could be precomputed
  nmod_poly_derivative(temp1, modulus);
  nmod_poly_invmod(temp1, temp1, modulus);
  
  nmod_poly_mulmod(result, temp1, temp2, modulus);
  
  nmod_poly_clear(temp1);
  nmod_poly_clear(temp2);
}

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
			  const fq_nmod_ctx_t ctx_from, const fq_nmod_ctx_t ctx_to) {
  slong degree = fq_nmod_ctx_degree(ctx_from);
  slong deg_low = fq_nmod_ctx_degree(ctx_to);
  fmpz_t deg_ratio;
  
  mp_limb_t *dual = _nmod_vec_init(degree); //new mp_limb_t[degree];
  for (slong i = 0; i < degree; i++)
    dual[i] = 0;

  // compute the dual basis image of a
  monomial_to_dual(dual, a, ctx_from);
  
  // compute the power projection <dual, ctx_to>
  project_powers(dual, dual, degree, g, ctx_from->modulus, ctx_from->inv);

  // compute result such that result(f) = g
  dual_to_monomial(res, dual, ctx_to->modulus);

  // correct scalar for proper embeddings
  fmpz_init_set_ui(deg_ratio, degree / deg_low);
  fmpz_invmod(deg_ratio, deg_ratio, &(ctx_to->p));
  fq_nmod_mul_fmpz(res, res, deg_ratio, ctx_to);
  
  _nmod_vec_clear(dual);
}

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
			 const fq_nmod_t g, const fq_nmod_ctx_t ctx_to) {
  nmod_poly_compose_mod_brent_kung_preinv(res, a, g, ctx_to->modulus, ctx_to->inv);
}
