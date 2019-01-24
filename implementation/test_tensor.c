#include <time.h>
#include "tensor.h"

int main() {
  struct tensor A;
  flint_rand_t state;
  fq_nmod_poly_t a, b, c, d;
  nmod_poly_t r;
  fq_nmod_poly_t mod, inv;
  fmpz_t coeff;
  int primes[3] = {3, 11, 31};
  time_t diff;
  
  flint_printf("Testing tensor\n");
  flint_randinit(state);
  for (int i = 100 ; i < 150 ; i++) {
    fmpz_init(coeff);
    fmpz_set_si(coeff, primes[i % 3]);
    fq_nmod_ctx_init(A.L, coeff, i, "x");
    nmod_poly_init(r, tensor_p(A));
    nmod_poly_randtest_monic(r, state, tensor_degree(A) + 1);
    fq_nmod_ctx_init_modulus(A.R, r, "z");

    fq_nmod_poly_init(a, A.L);
    fq_nmod_poly_init(b, A.L);
    fq_nmod_poly_init(c, A.L);
    fq_nmod_poly_init(d, A.L);
    
    fq_nmod_poly_randtest(a, state, tensor_level(A), A.L);
    fq_nmod_poly_randtest(b, state, tensor_level(A), A.L);
    diff = clock();
    tensor_mul(c, a, b, A);
    diff += -clock();

    fq_nmod_poly_init(mod, A.L);
    for (slong i = 0; i < A.R->modulus->length; i++) {
      fmpz_set_si(coeff, A.R->modulus->coeffs[i]);
      fq_nmod_poly_set_coeff_fmpz(mod, i, coeff, A.L);
    }
    fq_nmod_poly_init(inv, A.L);
    for (slong i = 0; i < A.R->inv->length; i++) {
      fmpz_set_si(coeff, A.R->inv->coeffs[i]);
      fq_nmod_poly_set_coeff_fmpz(inv, i, coeff, A.L);
    }
    diff -= clock();
    fq_nmod_poly_mulmod_preinv(d, a, b, mod, inv, A.L);
    diff += clock();
    
    if (!fq_nmod_poly_equal(c, d, A.L)) {
      flint_printf("Polynomials differ: ");
      fq_nmod_poly_print_pretty(c, "Z", A.L);
      flint_printf("   !=   ");
      fq_nmod_poly_print_pretty(d, "Z", A.L);
      flint_printf("\n");
    } else {
      flint_printf(".");
    }
    if (diff < -1000) {
      flint_printf("Warning: naive mul is faster by %g clocks", (double)diff);
    }

    nmod_poly_clear(r);
    fmpz_clear(coeff);
    fq_nmod_poly_clear(mod, A.L);
    fq_nmod_poly_clear(inv, A.L);
    fq_nmod_poly_clear(a, A.L);
    fq_nmod_poly_clear(b, A.L);
    fq_nmod_poly_clear(c, A.L);
    fq_nmod_poly_clear(d, A.L);
    fq_nmod_ctx_clear(A.L);
    fq_nmod_ctx_clear(A.R);
  }

  flint_randclear(state);
  flint_printf("done\n");
}
