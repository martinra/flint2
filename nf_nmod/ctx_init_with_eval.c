/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

nnnn=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Martin Raum
 
******************************************************************************/

#include <string.h>
#include "nmod_poly.h"
#include "nmod_poly_factor.h"
#include "nf_nmod.h"
#include "nf_nmod_mat.h"

void
nf_nmod_ctx_init_with_eval(nf_nmod_ctx_t ctx, const nmod_poly_t modulus, const nmod_mat_t evl_mat, const nmod_mat_t intrpl_mat, const char *var)
{
  /* note that we keep the reference to evl_mat and intrpl_mat.  So
   * the caller may not free them. */
  nmod_t mod = modulus->mod;

  /* set the variable name */
  ctx->var = flint_malloc(strlen(var) + 1);
  strcpy(ctx->var, var);

  /* set the modulus */
  nmod_poly_init(ctx->modulus, mod.n);
  nmod_poly_set(ctx->modulus, modulus);
  ctx->deg = nmod_poly_degree(modulus);
  ctx->separable = nmod_poly_is_squarefree(modulus);

  slong deg = ctx->deg;
  ulong big_size = 2 * deg - 1;

  /* set the evaluation and interpolation matrices */
  ctx->evl_mat = evl_mat;
  ctx->intrpl_mat = intrpl_mat;

  /* compute decomposition into copies of Fp and Fq */
  if (ctx->separable)
    {
      int i, j;
      mp_limb_t c;
      nmod_poly_factor_t modulus_factored;
      nmod_mat_t decomp_wdw;


      ctx->nfp = 0;
      ctx->nfq = 0;
      ctx->fp_moduli = (mp_limb_t *)flint_malloc(sizeof(mp_limb_t) * ctx->deg);
      ctx->fq_moduli = (nmod_poly_t *)flint_malloc(sizeof(nmod_poly_t *) * ctx->deg);

      nmod_mat_init(ctx->decomp_mat, deg, deg, mod.n);
      nmod_mat_init(ctx->recons_mat, deg, deg, mod.n);

      nmod_poly_factor_init(modulus_factored);


      nmod_poly_factor(modulus_factored, modulus);

      /* sort the factors according to whether they correspond to Fp or Fq */
      for (i = 0; i < modulus_factored->num; ++i)
	{
	  if (nmod_poly_degree(modulus_factored->p + i) == 1)
	    {
	      ctx->fp_moduli[ctx->nfp] = nmod_poly_get_coeff_ui(modulus_factored->p + i, 0);
	      ++ctx->nfp;
	    }
	  else
	    {
	      nmod_poly_init(ctx->fq_moduli[ctx->nfq]);
	      nmod_poly_set(ctx->fq_moduli[ctx->nfq], modulus_factored->p + i);
	      ++ctx->nfq;
	    }
	}

      /* set the decomposition matrix */
      for (i = 0; i < ctx->nfp; ++i)
	{
	  c = 1;
	  for (j = 0; j < ctx->deg; ++j)
	    {
	      nmod_mat_entry(ctx->decomp_mat, i, j) = c;
	  
	      if (j + 1 < ctx->deg)
		c = nmod_mul(c, ctx->fp_moduli[i], modulus->mod);
	    }
	}

      for (i = 0; i < ctx->nfq; ++i)
	{
	  /* todo: if evl_mat is defined, it might be faster to use
	     evaluation in order to compute reductions */
	  /* todo: implement nmod_mat_window_init and check (with
	     flint-devel) whether implemnetation is correct */
	  nmod_mat_window_init(decomp_wdw, ctx->decomp_mat, 0, ctx->nfp + i, ctx->deg, nmod_poly_degree(ctx->fq_moduli[i]));
	  _nf_nmod_reduction_mat(decomp_wdw, ctx->fq_moduli[i], ctx->deg);
	  nmod_mat_window_clear(decomp_wdw);
	}

      /* set the reconstruction matrix */
      nmod_mat_inv(ctx->recons_mat, ctx->decomp_mat);

      /* cleanup */
      nmod_poly_factor_clear(modulus_factored);
    }
}
