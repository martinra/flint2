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

  ulong p;
  fmpz_t p_prod;

  nf_nmod_ctx_t nmod_ctx;
  nf_nmod_mat_t A_nmod;

  nmod_mat_t * fp_rrefs;
  fq_nmod_mat_t * fq_rrefs;
  fq_nmod_ctx_t * fq_ctxs;

  size_t * nfz_rk_prof;
  size_t * nmod_rk_prof;
  size_t * nmod_comp_rk_prof;

  fmpz_t A_height_bd;
  fmpz_t B_height_bd;

  nfq_mat_t B_nfq;


  fmpz_init(p_prod);

  nfz_rk_prof = flint_malloc(sizeof(size_t) * nfz_mat_nrows(A));
  nmod_rk_prof = flint_malloc(sizeof(size_t) * nfz_mat_nrows(A));
  nmod_comp_rk_prof = flint_malloc(sizeof(size_t) * nfz_mat_nrows(A));

  fmpz_init(A_height_bd);
  fmpz_init(B_height_bd);


  // note: this make the rank profile smaller than any other rank profile
  // todo: check with later implementation
  nfz_rk_prof[0] = -1;
  nmod_rk_prof[0] = -1;

  A_height_bd = nfz_mat_height_bd(A);

  while (True)
    {
      break_modn_cycle = false;
      p = next_prime(1 << 27);

      // todo: move this into nmod ctx
      /* nmod_poly_init(nmod_modulus, p); */
      /* fmpz_poly_get_nmod_poly(nmod_modulus, ctx->modulus); */

      /* nmod_poly_factor_init(nmod_modulus_factored); */
      /* nmod_poly_factor(nmod_modulus_factored, nmod_modulus); */
      /* if (!nmod_poly_factor_squarefree(nmod_modulus)) */
      /* 	continue; */

      nf_nmod_ctx_init_by_nfz_ctx(nmod_ctx, ctx);
      if (!nf_nmod_ctx_is_separable(nmod_ctx))
	{
	  nf_nmod_ctx_clear(nmod_ctx);
	  continue;
	}

      nf_nmod_mat_init(A_nmod, nmod_ctx);
      nfz_mat_get_nmod_mat(A_nmod, A, ctx, nmod_ctx);

      // note: number of p and q components
      nmb_fps = nf_nmod_ctx_npcomp();
      nmb_fqs = nf_nmod_ctx_nqcomp();

      fp_rrefs = flint_malloc(sizeof(nmod_mat_t) * nmb_fps);
      fq_rrefs = flint_malloc(sizeof(fq_nmod_mat_t) * nmb_fqs);
      fq_ctxs = flint_malloc(sizeof(fq_ctx_t) * nmb_fqs);

      // note: decompose into elements corresponding to simple
      // components of the algebra Fp[X] / p(x)
      nf_nmod_mat_decompose_comp(fp_rrefs, fq_rrefs, fq_ctxs, A_nmod, nmod_ctx);
      nf_nmod_mat_clear(A_nmod, nmod_ctx);


      // todo: rename _flint_cp/cmp_rk_prof to some matrix function;
      // probably fmpz_mat is the best place for this
      for (int i = 0; i < nmb_fps; ++i)
	{
	  nmod_mat_rref(fp_rrefs[i]);
	  nmod_comp_rk_prof = nmod_mat_rk_prof(fp_rrefs[i]);

	  if (i == 0)
	    _flint_cp_rk_prof(nmod_rk_prof, nmod_comp_rk_prof);

	  int cmp = _flint_cmp_rk_prof(nmod_comp_rk_prof, nfz_rk_prof);
	  if (cmp > 0)
	    {
	      _flint_cp_rk_prof(nfz_rk_prof, nmod_comp_rk_prof);
	      if (i == 0)
		{
		  // todo: reinitialize lifts to nfz
		}
	    }
	  if (cmp < 0 || (i != 0 && cmp > 0))
	    {
	      nf_nmod_ctx_clear(nmod_ctx);

	      for (int j = 0; j < nmb_fps; ++j)
		nmod_mat_clear(fp_rrefs[j]);
	      flint_free(fp_rrefs);

	      for (int j = 0; j < nmb_fqs; ++j)
		{
		  fq_nmod_mat_clear(fq_rrefs[j], fq_ctxs[j]);
		  fq_nmod_ctx_clear(fq_ctxs[j]);
		}
	      flint_free(fq_rrefs);
	      flint_free(fq_ctxs);
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
	  nmod_comp_rk_prof = fq_nmod_mat_rk_prof(fp_rrefs[i]);

	  int cmp = _flint_cmp_rk_prof(nmod_comp_rk_prof, nfz_rk_prof);
	  if (cmp > 0)
	    _flint_cp_rk_prof(nfz_rk_prof, nmod_comp_rk_prof);
	  if (cmp != 0)
	    {
	      nf_nmod_ctx_clear(nmod_ctx);

	      for (int j = 0; j < nmb_fps; ++j)
		nmod_mat_clear(fp_rrefs[j]);
	      flint_free(fp_rrefs);

	      for (int j = 0; j < nmb_fqs; ++j)
		{
		  fq_nmod_mat_clear(fq_rrefs[j], fq_ctxs[j]);
		  fq_nmod_ctx_clear(fq_ctxs[j]);
		}
	      flint_free(fq_rrefs);
	      flint_free(fq_ctxs);
	      // todo: cleanup
		 
	      break_modn_cycle = True;
	      break;
	    }
	}
      if (break_modn_cycle)
	break;

      nf_nmod_mat_init(A_nmod, nmod_ctx);
      nf_nmod_reconstruct_comp(A_nmod, fp_rrefs, fq_rrefs, nmod_ctx);

      // todo: when implementing this, we should think about the
      // (exceptional) case that mp_limt_t is not ulong.  Then we need
      // to translate in some way
      nfz_mat_CRT_nmod(B, B, p_prod, A_nmod, p, 0);
      fmpz_mul_ui(p_prod, p_prod, p);

      nf_nmod_mat_clear(A_nmod, nmod_ctx);
      nf_nmod_ctx_clear(nmod_ctx);

      nfq_mat_init(B_nfq);
      if (nfq_mat_reconstruct_nfz(B_nfq, B, p_prod) == 0)
	{
	  nfq_mat_clear(B_nfq);

	  continue;
	}
      
      nfq_mat_get_nfz_mat_matwise(B, den, B_nfq);
      B_height_bd = B_den * nfz_mat_height_bd(B);

      if (p_prod >= nfz_mat_ncols(A) * A_height_bd * B_height_bd)
	{
	  slong rk = 0;	  
	  for (int i = 0; i < c; ++i)
	    {
	      if (nfz_rk_prof[i] != -1)
		++rk;
	      else
		break;
	    }

	  nfq_mat_clear(B_nfq);

	  fmpz_clear(p_prod);

	  flint_free(nfz_rk_prof);
	  flint_free(nmod_rk_prof);
	  flint_free(nmod_comp_rk_prof);

	  fmpz_clear(A_height_bd);
	  fmpz_clear(B_height_bd);

	  return rk;
	}
    }
}
