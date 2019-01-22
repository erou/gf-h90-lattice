#include <math.h>

#include "minpoly.h"

/**
 * Computes the minimal polynomial of the sequence {@code sequence}.
 * 
 * @param result	the minimal polynomial of degree <= degree
 * @param sequence	a sequence of length >= {@code 2 * degree}
 */
void minpoly_seq(nmod_poly_t result, const mp_limb_t *sequence, slong degree) {
  slong length = 2 * degree;

  slong slength = length;
  while(slength > 0 && sequence[slength-1] == 0)
    slength--;

  mp_limb_t *xpow;
  xpow = (mp_limb_t *) malloc((length+1)*sizeof(mp_limb_t));
  for (slong i = 0; i < length; i++)
    xpow[i] = 0;
  xpow[length] = 1;

  mp_limb_t *A;
  A = (mp_limb_t *) malloc((length+1)*sizeof(mp_limb_t));
  slong lenA;
  mp_limb_t *B;
  B = (mp_limb_t *) malloc(slength*sizeof(mp_limb_t));
  slong lenB;

  mp_limb_t *Minv[4];
  for (slong i = 0; i < 4; i++)
    Minv[i] = (mp_limb_t *) malloc((length+1)*sizeof(mp_limb_t));
  slong lenMinv[4];

  _nmod_poly_hgcd(Minv, lenMinv, A, &lenA, B, &lenB, xpow, length+1, sequence, slength, result->mod);

  for (slong i = 0; i < lenMinv[0]; i++)
    nmod_poly_set_coeff_ui(result, lenMinv[0]-1-i, Minv[0][i]);
  nmod_poly_make_monic(result, result);

  free(xpow);
  free(A);
  free(B);
  for (slong i = 0; i < 4; i++)
    free(Minv[i]);
}

/**
 * Computes the inner product of {@code a} and the coefficient vector of {@code f} 
 */
mp_limb_t inner_product(const mp_limb_t *a, const nmod_poly_t f) {
  slong len = nmod_poly_degree(f)+1;
  slong nlimbs = _nmod_vec_dot_bound_limbs(len, f->mod);
  return _nmod_vec_dot(a, f->coeffs, len, f->mod, nlimbs);
}

/*
 * Code for simple extension.
 */

/**
 * Computes the transposed product of {@code a} and {@code b}.
 * It is assumed that deg(a) <= deg(b) + m.
 * 
 * @param result	the transposed product of degree <= m
 */
void transposed_mul(nmod_poly_t result, const nmod_poly_t a, const nmod_poly_t b,
		    slong m) {

  slong n = nmod_poly_degree(b);
  if (n == -1) {
    nmod_poly_zero(result);
    return;
  }

  nmod_poly_t temp;
  nmod_poly_init(temp, b->mod.n);

  nmod_poly_reverse(temp, b, n + 1);
  nmod_poly_mullow(temp, a, temp, m + n + 1);
  nmod_poly_shift_right(temp, temp, n);
  nmod_poly_set(result, temp);

  nmod_poly_clear(temp);
}

/**
 * Computes the transposed remainder of {@code b} and {@code mod}.
 * 
 * @param a			the remainder of degree < n where n = deg(mod)
 * @param mod			the modulus of degree n
 * @param mod_inv_rev		1/rev(mod, n) mod x^{n - 1} where n = deg(mod)
 * @param result	the transposed remainder of degree <= m
 */
void transposed_rem(nmod_poly_t result, const nmod_poly_t a,
		    const nmod_poly_t mod,
		    const nmod_poly_t mod_inv_rev,
		    slong m) {

  nmod_poly_t temp;
  nmod_poly_init(temp, mod->mod.n);

  slong n = nmod_poly_degree(mod);

  transposed_mul(temp, a, mod, m - n);
  nmod_poly_mullow(temp, temp, mod_inv_rev, m - n + 1);
  nmod_poly_neg(temp, temp);
  nmod_poly_shift_left(temp, temp, n);
  nmod_poly_add(temp, a, temp);
  nmod_poly_set(result, temp);

  nmod_poly_clear(temp);
}

/**
 * Computes the transposed modular product of {@code a} and {@code b} modulo {@code mod}: b°a.
 * 
 * @param a			a vector of length deg(mod)
 * @param b			a polynomial of degree < deg(mod)
 * @param mod			the modulus
 * @param mod_inv_rev		1/rev(mod, n) mod x^{n - 1} where n = deg(mod)
 */
void transposed_mulmod(mp_limb_t *result, const mp_limb_t *a, const nmod_poly_t b, const nmod_poly_t mod,
		       const nmod_poly_t mod_rev_inv) {

  slong n = nmod_poly_degree(mod);
  slong m = nmod_poly_degree(b);

  nmod_poly_t temp;
  nmod_poly_init2(temp, mod->mod.n, n);

  for (slong i = 0; i < n; i++)
    nmod_poly_set_coeff_ui(temp, i, a[i]);

  //transposed_mulmod(temp, temp, b, mod, mod_rev_inv);
  if (m == -1) {
    nmod_poly_zero(temp);
  } else {

    nmod_poly_t temp2;
    nmod_poly_init(temp2, mod->mod.n);

    transposed_rem(temp2, temp, mod, mod_rev_inv, m + n - 1);
    transposed_mul(temp, temp2, b, n - 1);

    nmod_poly_clear(temp2);
  }
  /**/
	
  for (slong i = 0; i < n; i++)
    result[i] = nmod_poly_get_coeff_ui(temp, i);

  nmod_poly_clear(temp);
}

/**
 * Computes the power projection:
 * \[\langle v, 1 \rangle, \langle v, h \rangle, \dots, \langle v, h^{l - 1} \rangle\]
 * where $\langle ., . \rangle$ is the inner product, considering $h^i$ as the vector
 * of its coefficients.
 * 
 * @param modulus_inv_rev	1 / rev(m + 1, modulus) mod x^{m - 1} where m = deg(modulus)
 */
void project_powers(mp_limb_t *result, const mp_limb_t *a, slong l, const nmod_poly_t h,
		    const nmod_poly_t modulus, const nmod_poly_t modulus_inv_rev) {
  slong k = n_sqrt(l);
  slong m = (slong) ceil((double) l / (double) k);
  slong degree = nmod_poly_degree(modulus);

  // a temp for a to do computations without changing a
  mp_limb_t *temp_a = _nmod_vec_init(degree);
  for (slong i = 0; i < degree; i++)
    temp_a[i] = a[i];

  nmod_poly_t *h_powers = malloc((k+1) * sizeof(nmod_poly_t)); //new nmod_poly_t[k + 1];
  // initials h_powers
  for (slong i = 0; i <= k; i++)
    nmod_poly_init(h_powers[i], h->mod.n);

  // compute 1, h, h^2, ... h^k
  nmod_poly_set_coeff_ui(h_powers[0], 0, 1);
  nmod_poly_set(h_powers[1], h);
  for (slong i = 2; i <= k; i++)
    nmod_poly_mulmod_preinv(h_powers[i], h_powers[i - 1], h, modulus, modulus_inv_rev);
    
  slong base = 0;
  for (slong i = 0; i < m; i++) {

    for (slong j = 0; (j < k) && (base + j < l); j++)
      result[base + j] = inner_product(temp_a, h_powers[j]);

    // todo: preconditioning on h_powers[k]
    transposed_mulmod(temp_a, temp_a, h_powers[k], modulus, modulus_inv_rev);
    base += k;
  }

  _nmod_vec_clear(temp_a);
  for (slong i = 0; i <= k; i++)
    nmod_poly_clear(h_powers[i]);
  free(h_powers);
}

/**
 * Computes the minimal polynomial of {@code f} modulo {@code ctx}.
 * 
 * @param result	the minimal polynomial of {@code f}
 */

void minpoly(nmod_poly_t result,
	     const fq_nmod_t f, const fq_nmod_ctx_t ctx,
	     flint_rand_t state) {
  const nmod_poly_struct * modulus = ctx->modulus;
  const nmod_poly_struct * modulus_inv_rev = ctx->inv;
  nmod_poly_t g;
  nmod_poly_t temp_g;
  nmod_poly_t tau;

  slong degree = nmod_poly_degree(modulus);
  if (fq_nmod_is_zero(f, ctx)) {
    nmod_poly_zero(result);
    return;
  }

  nmod_poly_init(g, f->mod.n);
  nmod_poly_init(temp_g, f->mod.n);
  nmod_poly_init(tau, f->mod.n);

  nmod_poly_one(g);
  nmod_poly_one(tau);

  mp_limb_t *sequence = _nmod_vec_init(2 * degree);

  slong l = 0;
  while (1) {
    _nmod_vec_randtest(sequence, state, 2 * degree, f->mod);

    transposed_mulmod(sequence, sequence, tau, modulus, modulus_inv_rev);
    l = degree - nmod_poly_degree(g);
    project_powers(sequence, sequence, l * 2, f, modulus, modulus_inv_rev);
    minpoly_seq(temp_g, sequence, l);

    nmod_poly_mul(g, g, temp_g);
    if (nmod_poly_degree(g) == degree)
      break;

    nmod_poly_compose_mod(temp_g, temp_g, f, modulus);
    nmod_poly_mulmod(tau, tau, temp_g, modulus);
    if (nmod_poly_is_zero(tau))
      break;
  }

  nmod_poly_set(result, g);

  nmod_poly_clear(g);
  nmod_poly_clear(temp_g);
  nmod_poly_clear(tau);
  _nmod_vec_clear(sequence);
}
