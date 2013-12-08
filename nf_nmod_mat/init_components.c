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
_nf_nmod_mat_init_components(nmod_mat_t * fp_comps, fq_nmod_mat_t * fq_comps, fq_ctx_t * fq_ctxs, const nf_nmod_mat_t A, const nf_nmod_ctx_t ctx)
{
  int i, j;
  slong fq_deg;
  fmpz_poly_t fq_modulus;


  fp_comps = (nmod_mat_t *) flint_malloc(sizeof(nmod_mat_t) * ctx->nfp);
  fq_comps = (fq_nmod_mat_t *) flint_malloc(sizeof(fq_nmod_mat_t) * ctx->nfq);
  fq_ctxs = (fq_ctx_t *) flint_malloc(sizeof(fq_ctx_t) * ctx->nfq);

  for (i = 0; i < ctx->nfp; ++i)
    nmod_mat_init(fp_comps[i], A->r, A->c, ctx->modulus->mod.n);

  fmpz_poly_init(fq_modulus);
  for (i = 0; i < ctx->nfq; ++i)
    {
      fq_deg = nmod_poly_degree(ctx->fq_moduli[i]);
      fmpz_poly_zero(fq_modulus);
      for (j = 0; j < fq_deg; ++i)
	fmpz_poly_set_coeff_ui(fq_modulus, j, nmod_poly_get_coeff_ui(ctx->fq_moduli[i], j));

      fq_ctx_init_modulus(fq_ctxs[i], fq_modulus, ctx->var);
      fq_nmod_mat_init(fq_comps[i], A->r, A->c, fq_ctxs[i]);
    }
  fmpz_poly_clear(fq_modulus);
}
