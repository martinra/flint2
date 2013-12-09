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

#include "rank_profile.h"
#include "nfz_mat.h"
#include "nf_nmod_mat.h"

slong
_nfz_mat_rref_mod_prime_generator(nfz_mat_t B, fmpz_t den, const nfz_mat_t A, const nfz_ctx_t ctx, int (*next_prime)(const mp_limb_t))
{
  int cmp;

  ulong p;
  fmpz_t p_prod;

  slong rank, nmod_rank;

  nf_nmod_ctx_t ctx_nmod;
  nf_nmod_mat_t A_nmod;

  rank_profile_t nfz_rk_prof;
  rank_profile_t nmod_rk_prof;

  fmpz_t A_height_bd;
  fmpz_t B_height_bd;
  fmpz_t height_bd;

  nfq_mat_t B_nfq;


  fmpz_init(p_prod);

  rank_profile_init(nfz_rk_prof, A->r);
  rank_profile_init(nmod_rk_prof, A->r);

  fmpz_init(A_height_bd);
  fmpz_init(B_height_bd);
  fmpz_init(height_bd);


  rank_profile_entry(nfz_rk_prof, 0) = -1;
  rank_profile_entry(nmod_rk_prof, 0) = -1;

  /* todo: improve this.  the height bound here is (coeff bound reduction matrix)
     (coeff bound A) (coeff bound den * B) (A->c) */
  nfz_mat_coeff_bound(A_height_bd, A, ctx);
  _nfz_reduction_coeff_bound(height_bd, ctx->modulus, 2 * ctx->deg - 1);
  fmpz_mul(A_height_bd, A_height_bd, height_bd);
  fmpz_mul_ui(A_height_bd, A_height_bd, A->c);

  p = next_prime(1 << 27);
  while (true)
    {
      // todo: rename _flint_cp/cmp_rk_prof to some matrix function;
      // probably fmpz_mat is the best place for this
      p = next_prime(p);

      nfz_ctx_get_nmod_ctx(ctx_nmod, ctx, p);
      if (!nf_nmod_ctx_is_separable(ctx_nmod))
	{
	  nf_nmod_ctx_clear(ctx_nmod);
	  continue;
	}

      nf_nmod_mat_init(A_nmod, A->r, A->c, ctx_nmod);
      nfz_mat_get_nmod_mat(A_nmod, A, ctx_nmod, ctx);

      nf_nmod_mat_rref_components(A_nmod, A_nmod, ctx_nmod);
      nmod_rank = nf_nmod_mat_rank_profile(nmod_rk_prof, A_nmod, ctx_nmod);
      if (nmod_rank == -1)
	{
	  nf_nmod_ctx_clear(ctx_nmod);
	  nf_nmod_mat_clear(A_nmod, ctx_nmod);
	  continue;
	}
      cmp = rank_profile_cmp(nmod_rk_prof, nfz_rk_prof);
      if (cmp > 0)
	{
	  rank = nmod_rank;
	  rank_profile_set(nfz_rk_prof, nmod_rk_prof);

	  fmpz_set_ui(p_prod, 1);
	}
      else if (cmp < 0)
	{
	  nf_nmod_ctx_clear(ctx_nmod);
	  nf_nmod_mat_clear(A_nmod, ctx_nmod);
	}


      nfz_mat_CRT_nmod(B, B, p_prod, A_nmod, 0, ctx_nmod, ctx);
      fmpz_mul_ui(p_prod, p_prod, p);

      nf_nmod_mat_clear(A_nmod, ctx_nmod);
      nf_nmod_ctx_clear(ctx_nmod);

      nfq_mat_init(B_nfq);
      if (nfq_mat_reconstruct_nfz(B_nfq, B, p_prod) == 0)
	{
	  nfq_mat_clear(B_nfq);

	  continue;
	}
      
      nfq_mat_get_nfz_mat_matwise(B, den, B_nfq);
      nfz_mat_coeff_bound(B_height_bd, B, ctx);
      fmpz_mul(B_height_bd, B_height_bd, den);

      fmpz_mul(height_bd, A_height_bd, B_height_bd);
      if (fmpz_cmp(p_prod, height_bd) >= 0)
	{
	  slong rk = 0;	  
	  for (int i = 0; i < A->r; ++i)
	    {
	      if (rank_profile_entry(nfz_rk_prof, i) != -1)
		++rk;
	      else
		break;
	    }

	  nfq_mat_clear(B_nfq);

	  fmpz_clear(p_prod);

	  rank_profile_clear(nfz_rk_prof);
	  rank_profile_clear(nmod_rk_prof);

	  fmpz_clear(A_height_bd);
	  fmpz_clear(B_height_bd);
	  fmpz_clear(height_bd);

	  return rk;
	}
    }
}
