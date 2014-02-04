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
nfq_div(nfq_t f, const nfq_t g, const nfq_t h, const nfz_ctx_t ctx)
{
  fmpq * f_evl = _fmpq_vec_init(ctx->deg);
  fmpz * g_evl = _fmpz_vec_init(ctx->deg);
  fmpz * h_evl = _fmpz_vec_init(ctx->deg);

  fmpq_t den_quot;
  fmpq_init(den_quot);
  fmpq_set_fmpz_frac(den_quot, h->den, g->den);

  _nfz_eval(g_evl, g->coeffs, g->length, ctx);
  _nfz_eval(h_evl, h->coeffs, h->length, ctx);

  fmpz_t den;
  fmpz_init(den);
  fmpz_one(den);

  for (long i = 0; i < ctx->deg; ++i) {
    fmpq_set_fmpz_frac(f_evl + i, g_evl + i, h_evl + i);
    fmpq_mul(f_evl + i, f_evl + i, den_quot);
    fmpz_lcm(den, den, fmpq_denref(f_evl + i));
  }

  for (long i = 0; i < ctx->deg; ++i)
    fmpz_mul(g_evl + i, fmpq_numref(f_evl + i), den);

  long f_length;
  if (f->alloc < ctx->deg)
    fmpq_poly_realloc(f, ctx->deg);
  _nfz_interpolate(f->coeffs, &f_length, g_evl, ctx);
  fmpz_set(f->den, den);
  _fmpq_poly_set_length(f, f_length);

  _fmpq_vec_clear(f_evl, ctx->deg);
  _fmpz_vec_clear(g_evl, ctx->deg);
  _fmpz_vec_clear(h_evl, ctx->deg);

  fmpq_clear(den_quot);
  fmpz_clear(den);
}
