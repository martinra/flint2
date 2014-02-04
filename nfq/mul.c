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

#include "nfq.h"

void
nfq_mul(nfq_t f, const nfq_t g, const nfq_t h, const nfz_ctx_t ctx)
{
  fmpz * g_evl = _fmpz_vec_init(ctx->deg);
  fmpz * h_evl = _fmpz_vec_init(ctx->deg);

  _nfz_eval(g_evl, g->coeffs, g->length, ctx);
  _nfz_eval(h_evl, h->coeffs, h->length, ctx);

  for (long i = 0; i < ctx->deg; ++i)
    fmpz_mul(g_evl + i, g_evl + i, h_evl + i);

  long f_length;
  if (f->alloc < ctx->deg)
    fmpq_poly_realloc(f, ctx->deg);
  _nfz_interpolate(f->coeffs, &f_length, g_evl, ctx);
  _fmpq_poly_set_length(f, f_length);

  fmpz_t den;
  fmpz_init(den);
  fmpz_mul(den, g->den, h->den);
  fmpz_one(f->den);
  fmpq_poly_scalar_div_fmpz(f, f, den);
  fmpz_clear(den);

  _fmpz_vec_clear(g_evl, ctx->deg);
  _fmpz_vec_clear(h_evl, ctx->deg);
}
