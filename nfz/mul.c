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
  fmpz * g_evl = _fmpz_vec_init(ctx->evl_mat->r);
  fmpz * h_evl = _fmpz_vec_init(ctx->evl_mat->r);

  _nfz_eval(g_evl, g->coeffs, fmpz_poly_length(g), ctx);
  _nfz_eval(h_evl, h->coeffs, fmpz_poly_length(h), ctx);

  for (long i = 0; i < ctx->deg; ++i)
    fmpz_mul(g_evl + i, g_evl + i, h_evl + i);

  long f_length;
  if (f->alloc < ctx->deg)
    fmpz_poly_realloc(f, ctx->deg);
  _nfz_interpolate(f->coeffs, &f_length, g_evl, ctx);
  _fmpz_poly_set_length(f, f_length);

  _fmpz_vec_clear(g_evl, ctx->evl_mat->r);
  _fmpz_vec_clear(h_evl, ctx->evl_mat->r);
}
