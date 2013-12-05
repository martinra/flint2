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

#include "nfz_mat.h"

slong
_nfz_mat_rref_mod_prime_generator(nfz_mat_t B, fmpz_t den, const nfz_mat_t A, const nfz_ctx_t ctx, int (*next_prime)(const int))
{
  ulong p;
  fmpz_t p_prod;

  nf_ctx_nmod_t ctx_nmod;
  nf_nmod_mat_t A_nmod;

  size_t * nfz_rk_prof;
  size_t * nmod_rk_prof;

  fmpz_t A_height_bd;
  fmpz_t B_height_bd;

  nfq_mat_t B_nfq;


  fmpz_init(p_prod);

  nfz_rk_prof = flint_malloc(sizeof(size_t) * nfz_mat_nrows(A));
  nmod_rk_prof = flint_malloc(sizeof(size_t) * nfz_mat_nrows(A));

  fmpz_init(A_height_bd);
  fmpz_init(B_height_bd);


  // note: this make the rank profile smaller than any other rank profile
  // todo: check with later implementation
  nfz_rk_prof[0] = -1;
  nmod_rk_prof[0] = -1;

  A_height_bd = nfz_mat_height_bd(A);

  p = next_prime(1 << 27);
  while (True)
    {
      // todo: rename _flint_cp/cmp_rk_prof to some matrix function;
      // probably fmpz_mat is the best place for this
      p = next_prime(p);

      nf_ctx_nmod_init_by_nfz_ctx(ctx_nmod, ctx);
      if (!nf_ctx_nmod_is_separable(ctx_nmod))
	{
	  nf_ctx_nmod_clear(ctx_nmod);
	  continue;
	}

      nf_nmod_mat_init(A_nmod, ctx_nmod);
      nfz_mat_get_nmod_mat(A_nmod, A, ctx, ctx_nmod);

      nf_nmod_mat_rref_components(A_nmod, A_nmod, ctx_nmod);
      if (!nf_nmod_mat_rank_profile(nmod_rk_prof, A_nmod, ctx_nmod))
	{
	  nf_ctx_nmod_clear(ctx_nmod);
	  nf_nmod_mat_clear(A_nmod, ctx_nmod);
	  continue;
	}
      cmp = _flint_cmp_rk_prof(nmod_rk_prof, nfz_rk_prof);
      if (cmp > 0)
	{
	  _flint_cp_rk_prof(nfz_rk_prof, nmod_comp_rk_prof);

	  fmpz_set_ui(p_prod, 1);
	}
      else if (cmp < 0)
	{
	  nf_ctx_nmod_clear(ctx_nmod);
	  nf_nmod_mat_clear(A_nmod, ctx_nmod);
	}


      nfz_mat_CRT_nmod(B, B, p_prod, A_nmod, p, 0);
      fmpz_mul_ui(p_prod, p_prod, p);

      nf_nmod_mat_clear(A_nmod, ctx_nmod);
      nf_ctx_nmod_clear(ctx_nmod);

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
