#include "h90.h"
#include "AE.h"

/*
 * Compute the polynomial Σ a_i X^i of degree deg(A.R) with
 *
 *   a_{k-1} = a
 *   a_i = a_{i+1}^p + A.R[i+1] a
 *
 * Can use to reconstruct a H90 solution from its (k-1)-th coefficient.
 */
void lift_h90(fq_nmod_poly_t res,
	      const fq_nmod_t a, const struct tensor A) {
  slong k = tensor_ext_degree(A);

  if (fq_nmod_is_zero(a, A.L)) {
    fq_nmod_poly_zero(res, A.L);
    return;
  }
  
  fq_nmod_t temp;
  fq_nmod_init(temp, A.L);
  
  fq_nmod_poly_set_coeff(res, k-1, a, A.L);
  for (slong i = k-2; i >= 0; i--) {
    fq_nmod_frobenius(temp, res->coeffs + i + 1, 1, A.L);
    fq_nmod_mul_ui(res->coeffs + i, a, nmod_poly_get_coeff_ui(A.R->modulus, i+1), A.L);
    fq_nmod_add(res->coeffs + i, res->coeffs + i, temp, A.L);
  }

  fq_nmod_clear(temp, A.L);
}


/*
 * Evaluate
 * 
 *   Σ σ^i(a) ⊗ ζ^(-i)
 *
 * where a ∈ L, σ is the Frobenius of L, and ζ is the canonical root
 * of R.
 *
 * phi_div_R is the polynomial (X^n - 1) / R, where n is the order of
 * ζ and the degree of L.
 */
void _eval_cycle(fq_nmod_poly_t res, const fq_nmod_t a,
		 const struct tensor A, const nmod_poly_t phi_div_R) {
  nmod_poly_t b;
  
  nmod_poly_init(b, tensor_p(A));
  // eval phi_div_R(σ)(a) mod L
  compose(b, phi_div_R, a, A.L->modulus, A.L->inv);
  lift_h90(res, b, A);

  nmod_poly_clear(b);
}

/*
 * Find a solution res to the hilbert 90 equation
 * 
 *   (σ⊗1)(res) = ζ · res
 *
 * in the algebra L⊗R, where ζ is the canonical root of R
 */
void solve_h90(fq_nmod_poly_t res, flint_rand_t state,
	       const struct tensor A) {
  fq_nmod_t a;
  nmod_poly_t g;
  slong n = tensor_order(A);

  // Set g to (X^n - 1) / R
  // Can't this be improved?
  nmod_poly_init2_preinv(g, tensor_p(A), tensor_pinv(A), n+1);
  nmod_poly_set_coeff_ui(g, n, 1);
  nmod_poly_set_coeff_ui(g, 0, -1);
  nmod_poly_div(g, g, A.R->modulus);

  fq_nmod_init(a, A.L);

  do {
    fq_nmod_randtest(a, state, A.L);
    _eval_cycle(res, a, A, g);
  } while(fq_nmod_poly_is_zero(res, A.L));

  fq_nmod_clear(a, A.L);
  nmod_poly_clear(g);
}

/*
 * Tests if x∈A is a soultion of Hilbert 90 of order ord(A).
 */
int is_h90(const fq_nmod_poly_t x, const struct tensor A) {
  int ret = 1;
  slong k = tensor_ext_degree(A);
  fq_nmod_t temp1, temp2;
  fq_nmod_init(temp1, A.L);
  fq_nmod_init(temp2, A.L);
  
  for (slong i = 0; i < k; i++) {
    fq_nmod_mul_ui(temp1, x->coeffs + k - 1,
		   nmod_poly_get_coeff_ui(A.R->modulus, i), A.L);
    if (i > 0) {
      fq_nmod_sub(temp1, temp1, x->coeffs + i - 1, A.L);
    }
    fq_nmod_frobenius(temp2, x->coeffs + i, 1, A.L);
    fq_nmod_add(temp1, temp1, temp2, A.L);
    if (!fq_nmod_is_zero(temp1, A.L)) {
      ret = 0;
      break;
    }
  }
  
  fq_nmod_clear(temp2, A.L);
  fq_nmod_clear(temp1, A.L);
  return ret;
}
