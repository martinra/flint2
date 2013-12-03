/* rref_mod.c ---  */

/* Copyright (C) 2013 Martin Raum */

/* Author: Martin Raum */

/* This program is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU General Public License */
/* as published by the Free Software Foundation; either version 3 */
/* of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program. If not, see <http://www.gnu.org/licenses/>. */


#include "nfz_mat.h"

slong
nfz_mat_rref_mod(nfz_mat_t B, fmpz_t den, const nfz_mat_t A, const nfz_ctx_t ctx, int (*next_prime)(int))
{
  int ret;

  nf_nmod_ctx_t nmod_ctx;

  nmod_poly_t nmod_modulus;
  nmod_poly_factor_t nmod_modulus_factored;

  nmod_mat_t * fp_rrefs;
  fq_nmod_mat_t * fq_rrefs;

  ret = 0;

  while (True)
    {
      // todo: on 64 bit, use 27 instead
      p = next_prime(1 << 10);
     

      // todo: move this into nmod ctx
      nmod_poly_init(nmod_modulus, p);
      fmpz_poly_get_nmod_poly(nmod_modulus, ctx->modulus);

      nmod_poly_factor_init(nmod_modulus_factored);
      nmod_poly_factor(nmod_modulus_factored, nmod_modulus);
      if (!nmod_poly_factor_squarefree(nmod_modulus))
	continue;

      nf_nmod_ctx_init_by_nfz_ctx(nmod_ctx, ctx);

      nf_nmod_mat_init(A_nmod, nmod_ctx);
      nfz_mat_get_nmod_mat(A_nmod, A, ctx, nmod_ctx);

      nfz_nmod_decompose_to_simple_parts(&fp_rrefs, &fq_rrefs , A_nmod, nmod_ctx);

      nfz_nmod_reconstruct_from_simple_parts(A_nmod, fp_rrefs, fq_rrefs, nmod_ctx);

      // todo: next collect the reductions, compute pivots, and lift matrices
	  
      // todo: clear all variables
    }

  return ret;
}









