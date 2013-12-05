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

#include "nf_nmod.h"
#include "nmod_mat.h"
#include "perm.h"

void
nf_nmod_ctx_init(nf_nmod_ctx_t ctx, const nmod_poly_t modulus, const char *var)
{
  /* we compute the evaluation and interpolation matrix and then call
   * _nf_nmod_ctx_init_with_eval */

  /* note: this is the same code as int nfz/ctx_init.c converted to
     nmod matrices */

  slong deg = nmod_poly_degree(modulus);
  slong big_size = 2 * deg - 1;

  nmod_t mod = modulus->mod;

  if (big_size > mod->n)
    {
      /* in this case we cannot garentee that the evaluation and
	 interpolation will work */
      ctx->ev_mat = NULL;
      ctx->int_mat = NULL;

    _nf_nmod_ctx_init_with_eval(ctx, modulus, ev_mat, int_mat, *var);
    return;
    }

  slong i, j;

  mp_limb_t c, s;

  nmod_mat_t ev_mat;
  nmod_mat_t int_mat;
  nmod_mat_t ev_big_mat;
  nmod_mat_t int_big_mat;
  nmod_mat_t int_big_red_mat;
  nmod_mat_t red_mat;

  slong * rk_prof, perm;

  nmod_mat_init(ev_mat, deg, deg);
  nmod_mat_init(int_mat, deg, deg);

  nmod_mat_init(ev_big_mat, big_size, big_size);
  nmod_mat_init(int_big_mat, big_size, big_size);
  nmod_mat_init(int_big_red_mat, deg, big_size);
  nmod_mat_init(red_mat, big_size, deg);

  rk_prof = flint_malloc(sizeof(slong) * deg);

  // todo: use symmetry of these values to speed this up (hardly worth
  // it, if deg is not too large)
  for (i = 0; i < big_size; ++i)
    {
      c = 1;
      NMOD_RED(s, i - deg + 1, mod)
      for (j = 0; j < big_size; ++j)
	{
	  *(ev_mat->rows[i] + j) = c;
	  if (j + 1 < big_size)
	    c = nmod_mul(c, s);
	}
    }
  nmod_mat_inv(int_big_mat, int_den, ev_big_mat);

  _nf_nmod_reduction_mat(red_mat, modulus, big_size);
  /* The result is not the reduction matrix, but is stored in the same
   * variable to save space */
  nmod_mat_mul(red_mat, int_big_mat, red_mat);
  nmod_mat_transpose(int_big_red_mat, red_mat);

  /* c is used as a temporary variable */
  perm = _perm_init(nmod_mat_nrows(A));
  nmod_mat_lu(perm, int_big_red_mat, 0);
  _perm_clear(P);
  nmod_mat_rank_profile(rk_prof, int_big_red_mat);

  /* set the evaluation and interpolation matrices */
  for (i = 0; i < deg; ++i)
    for (j = 0; j < deg; ++j)
      {
	*(ev_mat->rows[i] + j) = nmod_mat_entry(ev_mat_big, i, rk_prof[j]);
	*(int_mat->rows[i] + j) = nmod_mat(int_mat, rk_prof[i], j);
      }

  /* set the variable name */
  nmod_mat_clear(ev_big_mat);
  nmod_mat_clear(int_big_mat);
  nmod_mat_clear(int_big_red_mat);
  nmod_mat_clear(red_mat);
  
  flint_free(rk_prof);


  _nf_nmod_ctx_init_with_eval(ctx, modulus, ev_mat, int_mat, *var);

  nmod_mat_clear(ev_mat);
  nmod_mat_clear(int_mat);
}
