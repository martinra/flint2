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

#include "nmod_mat.h"
#include "perm.h"
#include "nf_nmod.h"
#include "nf_nmod_mat.h"

void
nf_nmod_ctx_init(nf_nmod_ctx_t ctx, const nmod_poly_t modulus, const char *var)
{
  /* we compute the evaluation and interpolation matrix and then call
   * _nf_nmod_ctx_init_with_eval */

  /* note: this is the same code as int nfz/ctx_init.c converted to
     nmod matrices */

  slong deg = nmod_poly_degree(modulus);
  ulong big_size = 2 * deg - 1;

  nmod_t mod = modulus->mod;

  if (big_size > mod.n)
    {
      /* in this case we cannot garentee that the evaluation and
	 interpolation will work */
    _nf_nmod_ctx_init_with_eval(ctx, modulus, NULL, NULL, var);
    return;
    }

  slong i, j;

  mp_limb_t c, s;

  nmod_mat_t evl_mat;
  nmod_mat_t intrpl_mat;
  nmod_mat_t evl_big_mat;
  nmod_mat_t intrpl_big_mat;
  nmod_mat_t intrpl_big_red_mat;
  nmod_mat_t red_mat;

  nf_nmod_mat_rank_profile_t rk_prof;
  slong * perm;

  nmod_mat_init(evl_mat, deg, deg, mod.n);
  nmod_mat_init(intrpl_mat, deg, deg, mod.n);

  nmod_mat_init(evl_big_mat, big_size, big_size, mod.n);
  nmod_mat_init(intrpl_big_mat, big_size, big_size, mod.n);
  nmod_mat_init(intrpl_big_red_mat, deg, big_size, mod.n);
  nmod_mat_init(red_mat, big_size, deg, mod.n);

  nf_nmod_mat_init_rank_profile(rk_prof, deg);

  // todo: use symmetry of these values to speed this up (hardly worth
  // it, if deg is not too large)
  for (i = 0; i < big_size; ++i)
    {
      c = 1;
      NMOD_RED(s, i - deg + 1, mod);
      for (j = 0; j < big_size; ++j)
	{
	  *(evl_mat->rows[i] + j) = c;
	  if (j + 1 < big_size)
	    c = nmod_mul(c, s, mod);
	}
    }
  nmod_mat_inv(intrpl_big_mat, evl_big_mat);

  _nf_nmod_reduction_mat(red_mat, modulus, big_size);
  /* The result is not the reduction matrix, but is stored in the same
   * variable to save space */
  nmod_mat_mul(red_mat, intrpl_big_mat, red_mat);
  nmod_mat_transpose(intrpl_big_red_mat, red_mat);

  /* c is used as a temporary variable */
  perm = _perm_init(intrpl_big_red_mat->r);
  nmod_mat_lu(perm, intrpl_big_red_mat, 0);
  _perm_clear(perm);
  /* todo: implement nmod_mat_rank_profile */
  nmod_mat_rank_profile(rk_prof, intrpl_big_red_mat, ctx);

  /* set the evaluation and interpolation matrices */
  for (i = 0; i < deg; ++i)
    for (j = 0; j < deg; ++j)
      {
	*(evl_mat->rows[i] + j) = nmod_mat_entry(evl_big_mat, i, rk_prof[j]);
	*(intrpl_mat->rows[i] + j) = nmod_mat_entry(intrpl_mat, rk_prof[i], j);
      }

  /* set the variable name */
  nmod_mat_clear(evl_big_mat);
  nmod_mat_clear(intrpl_big_mat);
  nmod_mat_clear(intrpl_big_red_mat);
  nmod_mat_clear(red_mat);

  flint_free(rk_prof);


  _nf_nmod_ctx_init_with_eval(ctx, modulus, evl_mat, intrpl_mat, var);

  nmod_mat_clear(evl_mat);
  nmod_mat_clear(intrpl_mat);
}
