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

#include "nfz.h"

void
nfz_ctx_init(nfz_ctx_t ctx, const fmpz_poly_t modulus, const char *var)
{
  // todo: test irreducibility? Probably we don't need that, but it
  // could be nice to know that the polynomial is separable.

  slong i, j;
  slong deg = fmpz_poly_degree(modulus);
  slong big_size = 2 * deg - 1;

  fmpz_t c;

  fmpz_mat_t ev_big_mat;
  fmpz_mat_t int_big_mat;
  fmpz_mat_t int_big_red_mat;
  fmpz_mat_t red_mat;

  slong * rk_prof;

  fmpz_mat_init(ctx->ev_mat, deg, deg);
  fmpz_mat_init(ctx_>int_mat, deg, deg);
  fmpz_init(ctx->int_den);

  fmpz_init(c);

  fmpz_mat_init(ev_big_mat, big_size, big_size);
  fmpz_mat_init(int_big_mat, big_size, big_size);
  fmpz_mat_init(int_big_red_mat, deg, big_size);
  fmpz_mat_inti(red_mat, big_size, deg);

  rk_prof = flint_malloc(sizeof(slong) * deg);

  // todo: use symmetry of these values to speed this up (hardly worth
  // it, if deg is not too large)
  for (i = 0; i < big_size; ++i)
    {
      fmpz_set_ui(c, 1);
      for (j = 0; j < big_size; ++j)
	{
	  fmpz_set(fmpz_mat_entry(ev_mat, i, j), c);

	  if (j + 1 < big_size)
	    fmpz_mul_si(c, c, i - deg + 1);
	}
    }
  fmpz_mat_inv(int_big_mat, int_den, ev_big_mat);

  _nfz_reduction_mat(red_mat, modulus, big_size);
  /* The result is not the reduction matrix, but is stored in the same
   * variable to save space */
  fmpz_mat_mul(red_mat, int_big_mat, red_mat);
  fmpz_mat_transpose(int_big_red_mat, red_mat);

  /* c is used as a temporary variable */
  fmpz_mat_fflu(int_big_red_mat, c, NULL, int_big_red_mat, 0);
  fmpz_mat_rank_profile(rk_prof, int_big_red_mat);

  /* set the evaluation and interpolation matrices */
  for (i = 0; i < deg; ++i)
    for (j = 0; j < deg; ++j)
      {
	fmpz_set(fmpz_mat_entry(ev_mat_final, i, j), fmpz_mat_entry(ev_mat_big, i, rk_prof[j]));
	fmpz_set(fmpz_mat_entry(int_mat_final, i, j), fmpz_mat(int_mat, rk_prof[i], j));
      }

  /* set the variable name */
  ctx->var = flint_malloc(strlen(var) + 1);
  strcpy(ctx->var, var);

  /* set the modulus */
  fmpz_poly_init_set(ctx->modulus, modulus);
  ctx->deg = deg;

  fmpz_clear(c);
  fmpz_mat_clear(ev_big_mat);
  fmpz_mat_clear(int_big_mat);
  fmpz_mat_clear(int_big_red_mat);
  fmpz_mat_clear(red_mat);
  
  flint_free(rk_prof);
}

