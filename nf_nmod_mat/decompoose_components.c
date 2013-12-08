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

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Martin Raum
 
******************************************************************************/

#include "nf_nmod_mat.h"

void
nf_nmod_mat_decompose_components(nmod_mat_t * fp_comps, fq_nmod_mat_t * fq_comps, fq_ctx_t * fq_ctxs, const nf_nmod_mat_t A, const nf_nmod_ctx_t ctx)
{
  int i, j, k, l, n;

  slong fq_deg;

  nmod_mat_t coeff_mat;
  nmod_mat_t fq_comp_coeff_mat;

  mp_limb_t * decomp_row;

  decomp_row = ctx->decomp_mat->rows[0];
  for (i = 0; i < ctx->nfp; ++i)
    {
      nmod_mat_zero(fp_comps[i]);
      for (n = 0; n < ctx->deg; ++n)
	{
	  nf_nmod_mat_init_coeff_mat(coeff_mat, A, n, ctx);
	  /* todo: implement nmod_mat_add_scalar_mul */
	  nmod_mat_add_scalar_mul(fp_comps[i], fp_comps[i], decomp_row[n], coeff_mat);
	  nf_nmod_mat_clear_coeff_mat(coeff_mat, ctx);
	}

      ++decomp_row;
    }

  nmod_mat_init(fq_comp_coeff_mat, A->r, A->c, ctx->modulus->mod.n);
  for (i = 0; i < ctx->nfq; ++i)
    {
      fq_deg = nmod_poly_degree(ctx->fq_moduli[i]);
      for (j = 0; j < fq_deg; ++j)
	{
	  nmod_mat_zero(fq_comp_coeff_mat);
	  for (n = 0; n < ctx->deg; ++n)
	    {
	      nf_nmod_mat_init_coeff_mat(coeff_mat, A, n, ctx);
	      nmod_mat_add_scalar_mul(fq_comp_coeff_mat, fq_comp_coeff_mat, decomp_row[n], coeff_mat);
	      nf_nmod_mat_clear_coeff_mat(coeff_mat, ctx);
	    }
	  for (k = 0; k < A->r; ++k)
	    for (l = 0; l < A->c; ++l)
	      nmod_poly_set_coeff_ui(fq_nmod_mat_entry(fq_comps[i], k, l), n, nmod_mat_entry(coeff_mat, k, l));

	  ++decomp_row;
	}
    }
  nmod_mat_clear(fq_comp_coeff_mat);
}
