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

    Copyright (C) 2013-2014 Martin Raum
 
******************************************************************************/

#include <string.h>
#include "fmpz_poly.h"
#include "nfz.h"

void
nfz_ctx_init(nfz_ctx_t ctx, const fmpz_poly_t modulus, const char *var)
{
  long i, j;

  slong deg = fmpz_poly_degree(modulus);
  slong evl_size = 2 * deg - 1;

  fmpz_t c;

  fmpz_mat_init(ctx->evl_mat, evl_size, evl_size);
  fmpz_mat_init(ctx->intrpl_mat, evl_size, evl_size);
  fmpz_init(ctx->intrpl_den);
  
  // todo: use symmetry of these values to speed this up (hardly worth
  // it, if deg is not too large)
  for (i = 0; i < evl_size; ++i) {
      fmpz_set_ui(c, 1);
      for (j = 0; j < evl_size; ++j)
	{
	  fmpz_set(fmpz_mat_entry(ctx->evl_mat, i, j), c);

	  if (j + 1 < evl_size)
	    fmpz_mul_si(c, c, i - deg + 1);
	}
    }
  fmpz_mat_inv(ctx->intrpl_mat, ctx->intrpl_den, ctx->evl_mat);


  /* set the variable name */
  ctx->var = flint_malloc(strlen(var) + 1);
  strcpy(ctx->var, var);

  /* set the modulus */
  fmpz_poly_init(ctx->modulus);
  fmpz_poly_set(ctx->modulus, modulus);

  ctx->deg = deg;

  fmpz_clear(c);
}
