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

#include "fmpz_vec.h"

#include "nfz.h"

void
nfz_mul(nfz_t f, const nfz_t g, const nfz_t h, const nfz_ctx_t ctx)
{
  fmpz * g_evl = _fmpz_vec_init(ctx->deg);
  fmpz * h_evl = _fmpz_vec_init(ctx->deg);

  for (long i = 0; i < ctx->deg; ++i) {
    fmpz_zero(g_evl + i);
    fmpz_zero(h_evl + i);

    fmpz * evl_row = ctx->evl_mat->rows + i;
    for (long j = 0; j < ctx->deg; ++j) {
      fmpz_addmul(g_evl + i, g->coeffs + j, evl_row + j);
      fmpz_addmul(h_evl + i, h->coeffs + j, evl_row + j);
    }

    fmpz_mul(g_evl + i, g_evl + i, h_evl + i);
  }

  if (f->alloc < ctx->deg)
    fmpz_poly_realloc(f, ctx->deg);

  fmpz_poly_set_length(f, 0);
  for (long i = ctx->deg - 1; i >= 0; --i) {
    fmpz_zero(g->coeffs + i);

    fmpz * intrpl_row = ctx->intrpl_mat + i;
    for (long j = 0; j < ctx->deg; ++j)
      fmpz_addmul(f->coeffs + i, g_evl + j, intrpl_row + j);

    if (!fmpz_is_zero(f->coeffs + i))
      fmpz_poly_set_length(f, i + 1);
  }
}
