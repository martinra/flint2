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

#include "nf_nmod.h"

void
nf_nmod_ctx_init_with_eval(nf_nmod_ctx_t ctx, const nmod_poly_t modulus, const nmod_mat_t ev_mat, const nmod_mat_t int_mat, const char *var)
{
  /* set the variable name */
  ctx->var = flint_malloc(strlen(var) + 1);
  strcpy(ctx->var, var);

  /* set the modulus */
  nmod_poly_init_set(ctx->modulus, modulus);
  ctx->deg = nmod_poly_degree(modulus);
  ctx->separable = nmod_poly_is_squarefree(modulus);

  /* set the evaluation and interpolation matrices */
  if (ev_mat)
    nmod_mat_init_set(ctx->ev_mat, ev_mat);
  else
    ctx->ev_mat = NULL;

  if (int_mat)
    nmod_mat_init_set(ctx->int_mat, int_mat);
  else
    ctx->int_mat = NULL;

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

      nmod_mat_init(ctx->decomp_mat, deg, deg);
      nmod_mat_init(ctx->reconst_mat, deg, deg);

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
	      nmod_poly_init_set(ctx->fq_moduli[ctx->nfq], modulus_factored->p + i);
	      ++ctx->nfq;
	    }
	}

      /* set the decomposition matrix */
      for (i = 0; i < ctx->nfp; ++i)
	{
	  // todo: we probably have to implement nmod_mat_set_entry_ui
	  c = 1;
	  for (j = 0; j < ctx->deg; ++j)
	    {
	      nmod_mat_set_entry_ui(ctx->decomp_mat, i, c);
	  
	      if (j + 1 < ctx->deg)
		c = nmod_mul(c, ctx->fp_moduli[i], modulus->mod);
	    }
	}

      for (i = 0; i < ctx->nfq; ++i)
	{
	  /* todo: if ev_mat is defined, it might be faster to use
	     evaluation in order to compute reductions */
	  /* todo: implement nmod_mat_window_init and check (with
	     flint-devel) whether implemnetation is correct */
	  nmod_mat_window_init(decomp_wdw, ctx->decomp_mat, 0, ctx->nfp + i, ctx->deg, nmod_poly_degree(ctx->fq_moduli + i));
	  _nf_nmod_reduction_mat(decomp_wdw, ctx->fq_moduli + i, ctx->deg);
	  nmod_mat_window_clear(decomp_wdw);
	}

      /* set the reconstruction matrix */
      nmod_mod_mat_inv(ctx->reconst_mat, ctx->decomp_mat);


      /* cleanup */
      nmod_poly_factor_clean(modulus_factored);
    }
}
