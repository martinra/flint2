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
nf_nmod_mat_reconstruct_components(nf_nmod_mat_t A, const nmod_mat_t * fp_comps, const fq_nmod_mat_t * fq_comps, const fq_ctx_t fq_ctx, const nf_nmod_ctx_t ctx)
{
  int i, j, k, l, n;
  nmod_mat_t coeff_mat;
  nmod_mat_t fq_comp_coeff_mat;

  nf_nmod_mat_zero(B);

  for (i = 0; i < ctx->nfp; ++i)
    {
      for (n = 0; n < ctx->deg; ++n)
	{
	  nf_nmod_mat_init_coeff_mat(coeff_mat, A, n, ctx);
	  nmod_mat_add_scalar_mul(coeff_mat, coeff_mat, nmod_mat_entry(ctx->reconst_mat, i, n), fp_comps[i]);
	  nf_nmod_mat_clear_coeff_mat(coeff_mat);
	}
    }

  reconst_mat_row = ctx->nfp;
  nmod_mat_init(fq_comp_coeff_mat, A->r, A->c, ctx->mod.n);
  for (i = 0; i < ctx->nfp; ++i)
    for (j = 0; j < fq_ctx[i]->deg; ++j)
      {
	nmod_mat_init(fq_comp_coeff_mat, A->r, A->c, ctx->mod.n);
	for (k = 0; k < A->r; ++k)
	  for (l = 0; l < A->c; ++l)
	    nmod_mat_set_entry_ui(fq_comp_coeff_mat, k, l, fmpz_poly_get_coeff_ui(fq_nmod_entry(fq_comps[i], k, l), j));

	for (n = 0; n < ctx->deg; ++n)
	  {
	    nf_nmod_mat_init_coeff_mat(coeff_mat, A, n, ctx);
	    nmod_mat_add_scalar_mul(coeff_mat, coeff_mat, nmod_mat_entry(ctx->reconst_mat,  reconst_mat_row, n), fq_comp_coeff_mat);
	    nf_nmod_mat_clear_coeff_mat(coeff_mat);
	  }
	++reconst_mat_row;
      }
  nmod_mat_clear(fq_comp_coeff_mat);
}
