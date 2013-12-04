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
  bool break_modn_cycle;

  nf_nmod_ctx_t nmod_ctx;

  nmod_poly_t nmod_modulus;
  nmod_poly_factor_t nmod_modulus_factored;

  nmod_mat_t * fp_rrefs;
  fq_nmod_mat_t * fq_rrefs;

  size_t * nfz_rank_profile;
  size_t * nmod_rank_profile;
  size_t * nmod_component_rank_profile;
  nfz_rank_profile = flint_malloc(sizeof(size_t) * nfz_mat_nrows(A));
  nmod_rank_profile = flint_malloc(sizeof(size_t) * nfz_mat_nrows(A));
  nmod_component_rank_profile = flint_malloc(sizeof(size_t) * nfz_mat_nrows(A));
  nfz_rank_profile[0] = -1;
  nmod_rank_profile[0] = -1;

  while (True)
    {
      break_modn_cycle = false;
      p = next_prime(1 << 27);

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

      // todo: prepare context for fq_nmod_mat

      nmb_fps = nf_nmod_ctx_number_fp_components();
      nmb_fqs = nf_nmod_ctx_number_fq_components();

      fp_rrefs = flint_malloc(sizeof(nmod_mat_t) * nmb_fps);
      fq_rrefs = flint_malloc(sizeof(fq_nmod_mat_t) * nmb_fqs);
      fq_ctxs = flint_malloc(sizeof(fq_ctx_t) * nmb_fqs);

      nf_nmod_decompose_to_simple_parts(fp_rrefs, fq_rrefs, fq_ctxs, A_nmod, nmod_ctx);

      for (int i = 0; i < nmb_fps; ++i)
	{
	  nmod_mat_rref(fp_rrefs[i]);
	  nmod_component_rank_profile = nmod_mat_rank_profile(fp_rrefs[i]);

	  if (i == 0)
	    _flint_cp_rank_profile(nmod_rank_profile, nmod_component_rank_profile);

	  int cmp = _flint_cmp_rank_profile(nmod_component_rank_profile, nfz_rank_profile);
	  if (cmp > 0)
	    {
	      _flint_cp_rank_profile(nfz_rank_profile, nmod_component_rank_profile);
	      if (i == 0)
		{
		  // todo: reinitialize lifts to nfz
		}
	    }
	  if (cmp < 0 || (i != 0 && cmp > 0))
	    {
	      // todo: cleanup
		 
	      break_modn_cycle = True;
	      break;
	    }
	}
      if (break_modn_cycle)
	break;

      for (int i = 0; i < nmb_fqs; ++i)
	{
	  fq_nmod_mat_rref(fq_rrefs[i], fq_ctxs[i]);
	  nmod_component_rank_profile = fq_nmod_mat_rank_profile(fp_rrefs[i]);

	  int cmp = _flint_cmp_rank_profile(nmod_component_rank_profile, nfz_rank_profile);
	  if (cmp > 0)
	    _flint_cp_rank_profile(nfz_rank_profile, nmod_component_rank_profile);
	  if (cmp != 0)
	    {
	      // todo: cleanup
		 
	      break_modn_cycle = True;
	      break;
	    }
	}
      if (break_modn_cycle)
	break;

      nf_nmod_reconstruct_from_simple_parts(A_nmod, fp_rrefs, fq_rrefs, nmod_ctx);

      // todo: when implementing this, we should think about the
      // (exceptional) case that mp_limt_t is not ulong.  Then we need
      // to translate in some way
      nfz_mat_CRT_nmod(B, B, p_product, A_nmod, p, 0);
      p_product = p_product * p;
      if (nfq_mat_reconstruct_nfz(B_nfq, B, p_product) == 0)
	continue;
      
      // todo: check hight bounds and break loop if solution was found

      // todo: clear all variables
    }
}









