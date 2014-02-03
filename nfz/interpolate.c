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

#include "fmpz_vec.h"

void
_nfz_interpolate(nfz_t f, const fmpz * evl, const nfz_ctx_t ctx)
{
  if (f->alloc < ctx->deg)
    fmpz_poly_realloc(f, ctx->deg);

  _fmpz_poly_set_length(f, 0);
  for (long i = ctx->deg - 1; i >= 0; --i) {
    fmpz_zero(f->coeffs + i);

    fmpz * intrpl_row = ctx->intrpl_mat->rows[i];
    for (long j = 0; j < ctx->deg; ++j)
      fmpz_addmul(f->coeffs + i, evl + j, intrpl_row + j);

    if (!fmpz_is_zero(f->coeffs + i))
      _fmpz_poly_set_length(f, i + 1);
  }
}
