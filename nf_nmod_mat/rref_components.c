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

slong
nf_nmod_rref_components(nf_nmod_mat_t B, const nf_nmod_mat_t A, const nf_nmod_ctx_t ctx)
{
  slong rank, rank_comp;

  nmod_mat_t * fp_rrefs;
  fq_nmod_mat_t * fq_rrefs;
  fq_ctx_nmod_t * fq_ctxs;

  fp_rrefs = (nmod_mat_t *) flint_malloc(sizeof(nmod_mat_t) * ctx->nfq);
  fq_rrefs = (fq_nmod_mat_t *) flint_malloc(sizeof(fq_nmod_mat_t) * ctx->nfq);
  fq_ctxs = (fq_ctx_nmod_t *) flint_malloc(sizeof(fq_ctx_t) * ctx->nfq);

  // note: decompose into elements corresponding to simple
  // components of the algebra Fp[X] / p(x)
  nf_nmod_mat_init_components(fp_rrefs, fq_rrefs, fq_ctxs, A->r, A->c, ctx);
  nf_nmod_mat_decompose_comp(fp_rrefs, fq_rrefs, fq_ctxs, A, ctx);

  for (int i = 0; i < ctx->nfq; ++i)
    {
      rank_comp = nmod_mat_rref(fp_rrefs[i]);
      if (i == 0)
	rank = rank_comp;
      else if (rank != rank_comp)
	rank = -1;
    }

  for (int i = 0; i < nmb_fqs; ++i)
    {
      rank_comp = fq_nmod_mat_rref(fq_rrefs[i], fq_ctxs[i]);
      if (rank != rank_comp)
	rank = -1;
    }

  nf_nmod_reconstruct_comp(B, fp_rrefs, fq_rrefs, ctx);
  nf_nmod_mat_clear_components(fp_rrefs, fq_rrefs, fq_ctxs, ctx);

  return rank;
}
